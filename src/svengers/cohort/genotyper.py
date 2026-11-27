"""`svengers genotype`: force-genotype cohort sites in one sample, evidence-store only.

Reads the sites VCF from `svengers merge` and this sample's evidence store —
never the CRAM — and emits one genotype row per site. Sample-parallel by
construction: run one job per sample, then assemble with `svengers cohort-filter`.
"""

from __future__ import annotations

import csv
import gzip
import logging
from pathlib import Path

import numpy as np
import pysam

from ..config import Config
from ..extract.store import PE_DEL, PE_DUP, SR_LEFT_CLIP, SR_RIGHT_CLIP, EvidenceStore
from ..genotype.likelihood import ContigData, Genotyper
from ..segment.hmm import fit_sigma0
from ..variants import SOURCE_BOTH, SOURCE_DEPTH, CNVCall
from ..vcfout import span, spec_coords

logger = logging.getLogger(__name__)

GENOTYPE_COLUMNS = [
    "site_id", "contig", "start", "end", "svtype",
    "sample", "gt", "gq", "cn", "cnq", "rd", "pe", "sr", "cig", "baf_n", "baf_dev",
    "own",  # 1: this sample's own passing call nominated the site (sites.members.tsv.gz)
]

_SR_TOL = 100  # bp tolerance around site breakpoints for SR support


class StoreGenotyper:
    def __init__(self, store: EvidenceStore, config: Config,
                 rewritten_baf: bool | None = None):
        self.store = store
        self.config = config
        self.original_baf = not (config.genotype.rewritten_baf if rewritten_baf is None
                                 else rewritten_baf)
        self.genotyper = Genotyper(config.genotype, config.segment.max_cn)
        self.mean_depth = float(store.meta["mean_depth"])
        self.coverage = float(store.meta.get("coverage", self.mean_depth))
        self.ploidy: dict[str, int] = {
            k: int(v) for k, v in dict(store.meta.get("ploidy", {})).items()
        }
        insert_mu = float(store.meta.get("insert_mu", 0.0))
        insert_sigma = float(store.meta.get("insert_sigma", 0.0))
        self.pe_window = int(insert_mu + 3 * insert_sigma) if insert_mu else 1000
        self._data: dict[str, ContigData] = {}
        self._has_junctions: dict[str, bool] = {}
        self._junc: dict[str, dict] = {}

    def _contig_data(self, contig: str) -> ContigData | None:
        if contig not in self._data:
            if not self.store.has_contig(contig):
                self._data[contig] = None
            else:
                _raw, corr, mask = self.store.depth(contig)
                ratio = corr / self.mean_depth
                length = len(ratio) * self.store.bin_size
                baf = self.store.baf(contig, 0, length, original=self.original_baf)
                n_het = int(((baf["vaf"] > 0.15) & (baf["vaf"] < 0.85)).sum())
                callable_bp = int((~mask).sum()) * self.store.bin_size
                self._data[contig] = ContigData(
                    ratio=ratio,
                    mask=mask,
                    baf_pos=baf["pos"].astype(np.int64),
                    baf_vaf=baf["vaf"],
                    bin_size=self.store.bin_size,
                    ploidy=self.ploidy.get(contig, 2),
                    mean_depth=self.coverage,
                    sigma0=fit_sigma0(ratio, mask, self.ploidy.get(contig, 2)),
                    het_rate=n_het / callable_bp if callable_bp else None,
                    extra_sd=self.store.pon_sd(contig),
                )
                self._has_junctions[contig] = self.store.has_junction_evidence(contig)
        return self._data[contig]

    def _junction_arrays(self, contig: str) -> dict:
        """Whole-contig PE/SR/CIG arrays, loaded once: per-site ranged zarr
        reads decompress the same chunks over and over (profiled at ~9 sites/s
        on 238k cohort sites — days per cohort); in-memory searchsorted slices
        are ~100x faster for a few hundred MB per contig at most."""
        if contig not in self._junc:
            length = 2**62  # whole contig
            pe = self.store.pe(contig, 0, length)
            sr = self.store.sr(contig, 0, length)
            cig = self.store.cig(contig, 0, length)
            self._junc[contig] = {
                "pe_pos1": pe["pos1"].astype(np.int64),
                "pe_pos2": pe["pos2"].astype(np.int64),
                "pe_orient": pe["orient"],
                "sr_pos": sr["pos"].astype(np.int64),
                "sr_side": sr["side"],
                "cig_pos": cig["pos"].astype(np.int64),
                "cig_svlen": cig["svlen"].astype(np.int64),
            }
        return self._junc[contig]

    def _junction_support(
        self, contig: str, start: int, end: int, svtype: str
    ) -> tuple[int, int, int]:
        """(pe, sr, cig) support for the adjacency start->end in this sample."""
        j = self._junction_arrays(contig)
        w = self.pe_window
        orient = PE_DEL if svtype == "DEL" else PE_DUP
        n_pe = 0
        lo = int(np.searchsorted(j["pe_pos1"], max(start - w, 0), side="left"))
        hi = int(np.searchsorted(j["pe_pos1"], start + w, side="right"))
        if hi > lo:
            pos2 = j["pe_pos2"][lo:hi]
            sel = (
                (j["pe_orient"][lo:hi] == orient)
                & (pos2 >= end - w) & (pos2 <= end + w)
            )
            n_pe = int(sel.sum())
        n_sr = 0
        lo = int(np.searchsorted(j["sr_pos"], max(start - _SR_TOL, 0), side="left"))
        hi = int(np.searchsorted(j["sr_pos"], start + _SR_TOL, side="right"))
        if hi > lo:
            n_sr += int((j["sr_side"][lo:hi] == SR_RIGHT_CLIP).sum())
        lo = int(np.searchsorted(j["sr_pos"], max(end - _SR_TOL, 0), side="left"))
        hi = int(np.searchsorted(j["sr_pos"], end + _SR_TOL, side="right"))
        if hi > lo:
            n_sr += int((j["sr_side"][lo:hi] == SR_LEFT_CLIP).sum())
        n_cig = 0
        if svtype == "DEL":
            lo = int(np.searchsorted(j["cig_pos"], max(start - _SR_TOL, 0), side="left"))
            hi = int(np.searchsorted(j["cig_pos"], start + _SR_TOL, side="right"))
            if hi > lo:
                size = max(end - start, 1)
                sel = np.abs(j["cig_svlen"][lo:hi] - size) <= max(10, 0.2 * size)
                n_cig = int(sel.sum())
        return n_pe, n_sr, n_cig

    def nonref(
        self, contig: str, start: int, end: int, svtype: str, source: str,
        pe: int, sr: int, cig: int,
    ) -> tuple[tuple[int, int], int] | None:
        """(GT, GQ) of the most likely non-reference genotype of a call with
        this junction support, or None where the sample has no copy (chrY in
        a female) or the store lacks the contig."""
        data = self._contig_data(contig)
        if data is None or data.ploidy == 0:
            return None
        call = CNVCall(contig=contig, start=start, end=end, svtype=svtype, source=source,
                       pe_support=pe, sr_support=sr, cig_support=cig)
        self.genotyper.genotype(call, data)
        return call.gt_alt, call.gq_alt

    def genotype_site(
        self, site_id: str, contig: str, start: int, end: int, svtype: str,
        own: bool = False,
    ) -> dict | None:
        """Force-genotype one site. ``own``: the sample's own passing call
        nominated the site, so the site is known to be present in it and the
        genotype is the best non-reference one when the prior alone made it 0/0."""
        data = self._contig_data(contig)
        if data is None:
            return None
        if data.ploidy == 0:
            return None  # chrY in a female sample
        has_junc = self._has_junctions.get(contig, False)
        pe, sr, cig = (
            self._junction_support(contig, start, end, svtype) if has_junc else (0, 0, 0)
        )
        call = CNVCall(
            contig=contig, start=start, end=end, svtype=svtype,
            source=SOURCE_BOTH if has_junc else SOURCE_DEPTH,
            pe_support=pe, sr_support=sr, cig_support=cig,
        )
        self.genotyper.genotype(call, data)
        if own and call.gt == (0, 0) and call.gt_alt != (-1, -1):
            call.gt, call.gq = call.gt_alt, call.gq_alt
        gt = "./." if call.gt == (-1, -1) else f"{call.gt[0]}/{call.gt[1]}"
        return {
            "site_id": site_id, "contig": contig, "start": start, "end": end,
            "svtype": svtype, "sample": self.store.sample,
            "gt": gt, "gq": call.gq, "cn": call.cn,
            "cnq": round(call.cn_posterior, 4) if call.cn >= 0 else "",
            "rd": round(call.mean_ratio, 4) if not np.isnan(call.mean_ratio) else "",
            "pe": pe, "sr": sr, "cig": cig, "baf_n": call.baf_n,
            "baf_dev": round(call.baf_dev, 4) if not np.isnan(call.baf_dev) else "",
            "own": int(own),
        }


def read_own_sites(members: str | Path, sample: str) -> set[str]:
    """Site IDs nominated by ``sample``'s own passing calls, from the
    ``site_id<TAB>sample`` table `svengers merge` writes beside the sites."""
    tail = f"\t{sample}\n"
    with gzip.open(members, "rt") if str(members).endswith(".gz") else open(members) as fh:
        return {line.split("\t", 1)[0] for line in fh if line.endswith(tail)}


def run_genotype(
    sites_vcf: str | Path,
    store_path: str | Path,
    out_tsv: str | Path,
    config: Config,
    members: str | Path | None = None,
) -> int:
    store = EvidenceStore(store_path)
    gtr = StoreGenotyper(store, config)
    own = read_own_sites(members, store.sample) if members else set()
    n = 0
    with (
        pysam.VariantFile(str(sites_vcf)) as sites,
        open(out_tsv, "w", newline="") as fh,
    ):
        writer = csv.DictWriter(fh, fieldnames=GENOTYPE_COLUMNS, delimiter="\t")
        writer.writeheader()
        spec = spec_coords(sites.header)
        for rec in sites:
            row = gtr.genotype_site(rec.id, rec.contig, *span(rec, spec), rec.info["SVTYPE"],
                                    own=rec.id in own)
            if row is not None:
                writer.writerow(row)
                n += 1
    logger.info("genotype: %s -> %d site genotypes (%d nominated by its own calls)",
                store.sample, n, len(own))
    return n
