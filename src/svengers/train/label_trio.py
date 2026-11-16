"""Transmission labels for a truth-less sample from its parents' genotypes.

For a technically atypical child (wavy coverage) in a trio whose PARENTS are
clean, cohort re-genotyping gives near-ground-truth labels without GIAB:

- TP: a call matching a cohort site carried by >=1 parent at GQ >= min_gq —
  the event is transmitted, so the child's call is real. Parents must be
  non-wavy (the caller enforces nothing here; select trios accordingly),
  otherwise batch-shared artifacts would be labeled TP.
- FP: a call whose overlapping cohort sites are all genotyped 0/0 at
  GQ >= min_gq in BOTH parents — an apparent de novo, which in a wavy genome
  is an artifact with overwhelming probability (the genuine de novo CNV rate
  is <0.1/genome; wavy children show hundreds to thousands).
- ambiguous (dropped): no overlapping cohort site, or parents uncertain.
- common (dropped, with ``known_sites``): the call's site in a reference
  cohort is flagged as an artefact (CohortArtifact / HighCarrierFrac /
  CNDiscordant) or carried by >= ``max_ref_carrier_frac`` of it. A recurrent
  artefact is present in a parent too and would read as "transmitted": on
  4 wavy children, 99 % of the small junction DUPs labeled TP sat on such
  sites (TP rate 91 % there, 5-22 % elsewhere, 19 % against GIAB truth), and
  a model trained on them passed 12-17k calls per genome instead of ~5k.

Matching mirrors the merge compatibility: same svtype, reciprocal overlap
>= 0.5 against the parent genotype tables from ``svengers genotype``.
"""

from __future__ import annotations

import csv
import logging
from bisect import bisect_left
from pathlib import Path

logger = logging.getLogger(__name__)

_ARTEFACT_FLAGS = ("CohortArtifact", "HighCarrierFrac", "CNDiscordant")


def _load_parent_sites(gt_tsv: str | Path, min_gq: int) -> dict[tuple[str, str], list]:
    """{(contig, svtype): sorted [(start, end, carrier, confident_ref)]}."""
    out: dict[tuple[str, str], list] = {}
    with open(gt_tsv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            gt, gq = r["gt"], int(r["gq"] or 0)
            if gt == "./.":
                continue
            carrier = "1" in gt and gq >= min_gq
            confident_ref = gt == "0/0" and gq >= min_gq
            if not carrier and not confident_ref:
                continue
            out.setdefault((r["contig"], r["svtype"]), []).append(
                (int(r["start"]), int(r["end"]), carrier, confident_ref)
            )
    for rows in out.values():
        rows.sort()
    return out


def _overlapping(sites: list, start: int, end: int, min_ro: float = 0.5):
    """Sites reciprocally overlapping [start, end) by >= min_ro."""
    size = max(end - start, 1)
    lo = bisect_left(sites, (start - 10 * size - 1_000_000,))
    hits = []
    for s, e, carrier, confident_ref in sites[lo:]:
        if s >= end:
            break
        inter = min(e, end) - max(s, start)
        if inter > 0 and inter / max(size, e - s) >= min_ro:
            hits.append((carrier, confident_ref))
    return hits


def _load_own_calls(vcf_path: str | Path) -> dict[tuple[str, str], list]:
    """{(contig, svtype): sorted [(start, end, True, True)]}: every call of the
    sample's own VCF, whatever its FILTER (the call exists)."""
    import pysam

    from ..vcfout import span, spec_coords

    out: dict[tuple[str, str], list] = {}
    with pysam.VariantFile(str(vcf_path)) as vf:
        spec = spec_coords(vf.header)
        for rec in vf:
            svtype = rec.info.get("SVTYPE")
            if svtype in ("DEL", "DUP"):
                out.setdefault((rec.contig, svtype), []).append((*span(rec, spec), True, True))
    for rows in out.values():
        rows.sort()
    return out


def label_by_transmission(
    features_tsv: str | Path,
    parent_gt_tsvs: list[str | Path],
    out_tsv: str | Path,
    min_gq: int = 20,
    min_ro: float = 0.5,
    known_sites: str | Path | None = None,
    max_ref_carrier_frac: float = 0.5,
    parent_calls: list[str | Path] | None = None,
) -> dict[str, int]:
    """Label a child's feature rows TP/FP by parental transmission.

    With ``parent_calls`` (each parent's own VCF, same order as
    ``parent_gt_tsvs``), re-genotyping alone no longer decides: TP needs a
    parent that is a carrier AND has its own call there, FP needs both
    parents confidently 0/0 AND neither with a call. Carriage by forced
    re-genotyping without an own call was false 33/33 in the IGV review of
    PMS (ROADMAP 1septies), so it would label artefacts as transmitted.
    """
    parents = [_load_parent_sites(p, min_gq) for p in parent_gt_tsvs]
    own = [_load_own_calls(p) for p in parent_calls] if parent_calls else None
    if own is not None and len(own) != len(parents):
        raise ValueError("parent_calls must list one VCF per parent genotype table")
    ref = None
    if known_sites:
        from ..cohort.metrics import _KnownSites

        ref = _KnownSites(known_sites, min_ro=min_ro)
    counts = {"tp": 0, "fp": 0, "ambiguous": 0, "common": 0}
    with open(features_tsv) as fin, open(out_tsv, "w", newline="") as fout:
        reader = csv.DictReader(fin, delimiter="\t")
        writer = csv.DictWriter(
            fout, fieldnames=[*reader.fieldnames, "label"], delimiter="\t"
        )
        writer.writeheader()
        for r in reader:
            contig, svtype = r["contig"], r["svtype"]
            start, end = int(r["start"]), int(r["end"])
            carrier = False
            confident_ref_both = True
            any_overlap = False
            for k, psites in enumerate(parents):
                hits = _overlapping(
                    psites.get((contig, svtype), []), start, end, min_ro
                )
                called = own is not None and bool(
                    _overlapping(own[k].get((contig, svtype), []), start, end, min_ro)
                )
                if hits:
                    any_overlap = True
                if any(h[0] for h in hits) and (own is None or called):
                    carrier = True
                if not any(h[1] for h in hits) or called:
                    confident_ref_both = False
            if (carrier or (any_overlap and confident_ref_both)) and ref is not None:
                hit = ref.lookup(contig, start, end, svtype)
                if hit is not None and (
                    hit["carrier_frac"] >= max_ref_carrier_frac
                    or any(f in hit["filter"] for f in _ARTEFACT_FLAGS)
                ):
                    counts["common"] += 1
                    continue
            if carrier:
                r["label"] = 1
                counts["tp"] += 1
            elif any_overlap and confident_ref_both:
                r["label"] = 0
                counts["fp"] += 1
            else:
                counts["ambiguous"] += 1
                continue
            writer.writerow(r)
    logger.info(
        "label-trio %s: %d TP, %d FP, %d ambiguous, %d common/artefact (dropped)",
        features_tsv, counts["tp"], counts["fp"], counts["ambiguous"], counts["common"],
    )
    return counts
