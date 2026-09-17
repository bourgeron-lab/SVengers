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

Matching mirrors the merge compatibility: same svtype, reciprocal overlap
>= 0.5 against the parent genotype tables from ``vedette genotype``.
"""

from __future__ import annotations

import csv
import logging
from bisect import bisect_left
from pathlib import Path

logger = logging.getLogger(__name__)


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


def label_by_transmission(
    features_tsv: str | Path,
    parent_gt_tsvs: list[str | Path],
    out_tsv: str | Path,
    min_gq: int = 20,
    min_ro: float = 0.5,
) -> dict[str, int]:
    """Label a child's feature rows TP/FP by parental transmission."""
    parents = [_load_parent_sites(p, min_gq) for p in parent_gt_tsvs]
    counts = {"tp": 0, "fp": 0, "ambiguous": 0}
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
            for psites in parents:
                hits = _overlapping(
                    psites.get((contig, svtype), []), start, end, min_ro
                )
                if hits:
                    any_overlap = True
                if any(h[0] for h in hits):
                    carrier = True
                if not any(h[1] for h in hits):
                    confident_ref_both = False
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
        "label-trio %s: %d TP, %d FP, %d ambiguous (dropped)",
        features_tsv, counts["tp"], counts["fp"], counts["ambiguous"],
    )
    return counts
