"""Join a sample's depth calls that one event left in abutting pieces.

Since v1.12.0, `svengers call` joins same-CN depth segments split by a short
interruption (`segment.hmm.bridge_segments`: a single outlier bin, or the
2 kb CN0 of a common deletion on the remaining allele inside a terminal
deletion). Calls made before (the production cohorts were called earlier)
still carry the pieces: a 22q13 terminal deletion comes out as 43.6-49.32 Mb
and 49.322 Mb-qter. ``bridge_calls`` applies the same rule to the raw calls,
from the evidence store alone, before the filter scores them:

- pieces are depth-only calls (SOURCE=DEPTH, NBINS >= 1, not MaxN) of one
  contig, type and integer CN, taken in start order; a call supported by a
  junction is never joined (its breakpoints are evidence that the event ends
  there), nor a mostly masked one (it would carry the event into a repeat: a
  de novo 16p11.2 deletion lost its PASS when joined to a piece of the BP4
  segmental duplication);
- the next piece joins the run when the gap is at most
  ``segment.bridge_max_gap_bins`` bins and ``segment.bridge_max_gap_frac`` of
  the shorter side, and the mean copy ratios differ by at most
  ``segment.bridge_max_ratio_diff`` (as at calling time); joins chain;
- the joined call spans the gap, keeps the first piece's ID and CIPOS and the
  last piece's CIEND, gets NBINS / MASKFRAC from the store and GT, GQ, CN, CNQ,
  RD, BAFN, BAFDEV from the store genotyper, BRIDGED = number of pieces, and
  its own call-time tags (NoEvidence, MaxN, LowPosterior, MinSize; LowPosterior
  when pieces so tagged make up most of its bins); calls inside the gap are
  kept unchanged.

On calls made with bridging the pass is a no-op.
"""

from __future__ import annotations

import logging
import math
from pathlib import Path

import numpy as np
import pysam

from ..cohort.genotyper import StoreGenotyper
from ..config import Config
from ..extract.store import EvidenceStore
from ..variants import SOURCE_DEPTH, CNVCall
from ..vcfout import span, spec_coords
from .regenotype import write_genotype

logger = logging.getLogger(__name__)

BRIDGED_INFO = ("BRIDGED", "1", "Integer",
                "Number of abutting depth calls of one event joined into this one at the filter "
                "stage (calls made before call-time bridging)")
REGENOTYPED_KEYS = ("svengers_regenotyped", "vedette_regenotyped")
CALL_TAGS = ("NoEvidence", "MaxN", "LowPosterior", "MinSize")   # merge_candidates.segment_to_call


def _piece(rec, spec: bool):
    """(start, end, cn, ratio) of a joinable depth call, else None."""
    if rec.info.get("SVTYPE") not in ("DEL", "DUP") or rec.info.get("SOURCE") != SOURCE_DEPTH:
        return None
    if int(rec.info.get("NBINS") or 0) < 1 or "MaxN" in rec.filter.keys():
        return None   # a mostly masked piece would carry the event into the repeat
    fmt = rec.samples[0]
    cn, rd = fmt.get("CN"), fmt.get("RD")
    if cn is None or cn < 0 or rd is None or not math.isfinite(rd):
        return None
    start, end = span(rec, spec)
    return start, end, int(cn), float(rd)


def join_runs(pieces: list[tuple[int, int, int, float]], bin_size: int, max_gap_bins: int,
              max_gap_frac: float, max_ratio_diff: float) -> list[list[int]]:
    """Runs (lists of indices into ``pieces``, start-ordered, >= 2) of
    same-CN pieces to join. ``pieces`` are (start, end, cn, ratio) of one
    contig and type."""
    runs: list[list[int]] = []
    by_cn: dict[int, list[int]] = {}
    for i, p in enumerate(pieces):
        by_cn.setdefault(p[2], []).append(i)
    for idx in by_cn.values():
        idx.sort(key=lambda i: (pieces[i][0], pieces[i][1]))
        run = [idx[0]]
        run_bins = (pieces[idx[0]][1] - pieces[idx[0]][0]) / bin_size
        run_ratio = pieces[idx[0]][3]
        for i in idx[1:] + [None]:
            if i is not None:
                s, e, _cn, r = pieces[i]
                gap = (s - pieces[run[-1]][1]) / bin_size
                n = (e - s) / bin_size
                if (0 <= gap <= max_gap_bins and gap <= max_gap_frac * min(run_bins, n)
                        and abs(run_ratio - r) <= max_ratio_diff):
                    run_ratio = (run_ratio * run_bins + r * n) / (run_bins + n)
                    run_bins += n
                    run.append(i)
                    continue
            if len(run) > 1:
                runs.append(run)
            if i is not None:
                run = [i]
                run_bins = (pieces[i][1] - pieces[i][0]) / bin_size
                run_ratio = pieces[i][3]
    return runs


def bridge_calls(vcf_in: str | Path, store_path: str | Path, vcf_out: str | Path,
                 config: Config | None = None) -> dict[str, int]:
    """Write ``vcf_out``: ``vcf_in`` with abutting depth pieces joined.
    Returns counts: calls (in), joined (runs), absorbed (pieces removed)."""
    config = config or Config()
    store = EvidenceStore(store_path)
    counts = {"calls": 0, "joined": 0, "absorbed": 0}
    with pysam.VariantFile(str(vcf_in)) as vin:
        spec = spec_coords(vin.header)
        recs = list(vin)
        hdr = vin.header.copy()
        _write(recs, hdr, spec, store, config, vcf_out, counts)
    pysam.tabix_index(str(vcf_out), preset="vcf", force=True)
    if counts["joined"]:
        logger.info("bridge: %d runs of abutting depth calls joined (%d pieces absorbed)",
                    counts["joined"], counts["absorbed"])
    return counts


def _write(recs, hdr, spec, store, config, vcf_out, counts) -> None:
    seg = config.segment
    bs = store.bin_size
    counts["calls"] = len(recs)
    rewritten = any(r.key in REGENOTYPED_KEYS for r in hdr.records)
    groups: dict[tuple[str, str], list[int]] = {}
    pieces: dict[int, tuple[int, int, int, float]] = {}
    if seg.bridge_max_gap_bins > 0:
        for i, rec in enumerate(recs):
            p = _piece(rec, spec)
            if p is not None:
                pieces[i] = p
                groups.setdefault((rec.contig, rec.info["SVTYPE"]), []).append(i)
    joined: dict[int, list[int]] = {}     # first record index -> run (record indices)
    for ids in groups.values():
        for run in join_runs([pieces[i] for i in ids], bs, seg.bridge_max_gap_bins,
                             seg.bridge_max_gap_frac, seg.bridge_max_ratio_diff):
            members = sorted((ids[k] for k in run), key=lambda i: pieces[i][0])
            joined[members[0]] = members
    absorbed = {i for m in joined.values() for i in m[1:]}
    if joined:
        name, number, vtype, desc = BRIDGED_INFO
        if name not in hdr.info:
            hdr.info.add(name, number, vtype, desc)
        sg = StoreGenotyper(store, config, rewritten_baf=rewritten)
    with pysam.VariantFile(str(vcf_out), "wz", header=hdr) as vout:
        for i, rec in enumerate(recs):
            if i in absorbed:
                continue
            rec.translate(hdr)
            if i in joined:
                members = joined[i]
                start = pieces[members[0]][0]
                end = max(pieces[m][1] for m in members)
                last = recs[members[-1]]
                svtype = rec.info["SVTYPE"]
                rec.stop = end if spec else end + 1
                rec.info["SVLEN"] = (end - start) * (1 if svtype == "DUP" else -1)
                if "CIEND" in last.info:
                    rec.info["CIEND"] = last.info["CIEND"]
                _raw, _corr, mask = store.depth(rec.contig, start, end)
                rec.info["NBINS"] = int(len(mask))
                rec.info["MASKFRAC"] = round(float(np.asarray(mask, bool).mean()), 4) if len(mask) else 0.0
                rec.info["BRIDGED"] = len(members)
                # call-time tags of the joined segment (`segment_to_call`), not the first piece's;
                # LowPosterior when pieces so tagged make up most of the joined bins (the
                # posterior itself is not in the VCF)
                low_post = sum(pieces[m][1] - pieces[m][0] for m in members
                               if "LowPosterior" in recs[m].filter.keys())
                tags = [f for f in rec.filter.keys() if f not in CALL_TAGS]
                if len(mask) < seg.min_bins:
                    tags.append("NoEvidence")
                if rec.info["MASKFRAC"] > seg.max_masked_frac:
                    tags.append("MaxN")
                if low_post > 0.5 * sum(pieces[m][1] - pieces[m][0] for m in members):
                    tags.append("LowPosterior")
                if end - start < config.breakpoint.min_size:
                    tags.append("MinSize")
                rec.filter.clear()
                for t in tags or ["PASS"]:
                    rec.filter.add(t)
                call = CNVCall(contig=rec.contig, start=start, end=end, svtype=svtype,
                               source=SOURCE_DEPTH)
                data = sg._contig_data(rec.contig)
                if data is not None and data.ploidy > 0:
                    sg.genotyper.genotype(call, data)
                    write_genotype(rec.samples[0], call)
                counts["joined"] += 1
                counts["absorbed"] += len(members) - 1
            vout.write(rec)
