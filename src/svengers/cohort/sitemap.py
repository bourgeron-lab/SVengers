"""Correspondence of cohort sites between two runs (issue #9).

A cohort site id (``svengers_site_N``) is positional and renumbered at every
merge, and a site's bounds move when the merge sees other breakpoints. A
curator's validation keyed on a site must follow the event to the new run.
``map_sites`` pairs the sites of two cohort VCFs, same type, in three passes:

1. same bounds: ``same``;
2. the best match of each side by reciprocal overlap (>= ``min_ro``, 0.5 as
   the merge's own compatibility rule) when it is mutual: ``moved``;
3. a site left over whose best match on the other side is already paired:
   an old site cut in several (``split``) or several old sites joined in one
   (``merged``); the pairs of that group take the label;

and ``lost`` / ``new`` for the sites left without a match. Mutual matching
first matters: in a real cohort many sites overlap neighbours of the same
type by more than half (nested and stacked clusters), and pairing every
overlap read two identical PMS VCFs as 92,004 "complex" sites.

Keys are the family export's ``site_key`` (``chr:start-end:del|dup``,
1-based inclusive), read from either coordinate convention.
"""

from __future__ import annotations

import bisect
import csv
from collections import defaultdict
from pathlib import Path

import pysam

from ..vcfout import span as vcf_span
from ..vcfout import spec_coords

MAP_COLUMNS = ["old_site_id", "old_site_key", "new_site_id", "new_site_key", "relation",
               "overlap"]


def read_sites(vcf_path: str | Path) -> dict[tuple[str, str], list[tuple[int, int, str, str]]]:
    """{(contig, svtype): sorted [(start, end, site_id, site_key)]}, 0-based half-open."""
    out: dict[tuple[str, str], list] = defaultdict(list)
    with pysam.VariantFile(str(vcf_path)) as vf:
        spec = spec_coords(vf.header)
        for rec in vf.fetch() if vf.index is not None else vf:
            svtype = rec.info.get("SVTYPE")
            if svtype not in ("DEL", "DUP"):
                continue
            start, end = vcf_span(rec, spec)
            out[(rec.contig, svtype)].append(
                (start, end, rec.id or "", f"{rec.contig}:{start + 1}-{end}:{svtype.lower()}"))
    for v in out.values():
        v.sort()
    return out


def _ro(a0: int, a1: int, b0: int, b1: int) -> float:
    inter = min(a1, b1) - max(a0, b0)
    return inter / max(a1 - a0, b1 - b0, 1) if inter > 0 else 0.0


def map_sites(old_vcf: str | Path, new_vcf: str | Path,
              min_ro: float = 0.5) -> tuple[list[dict], dict[str, int]]:
    """(rows, counts by relation) of the correspondence old -> new."""
    old, new = read_sites(old_vcf), read_sites(new_vcf)
    pairs: list[tuple[tuple, tuple, str, float]] = []
    lost, born = [], []
    for key in sorted(set(old) | set(new)):
        olds, news = old.get(key, []), new.get(key, [])
        # 1. same bounds
        by_bounds: dict[tuple[int, int], list] = defaultdict(list)
        for n in news:
            by_bounds[(n[0], n[1])].append(n)
        free_old, used_new = [], set()
        for o in olds:
            twins = by_bounds.get((o[0], o[1]))
            if twins:
                n = twins.pop()
                used_new.add(n)
                pairs.append((o, n, "same", 1.0))
            else:
                free_old.append(o)
        free_new = [n for n in news if n not in used_new]
        # 2. mutual best matches among the rest
        starts = [n[0] for n in free_new]
        cand: dict[tuple, list[tuple[float, tuple]]] = defaultdict(list)
        for o in free_old:
            # a match of reciprocal overlap >= min_ro starts within one size
            lo = bisect.bisect_left(starts, o[0] - (o[1] - o[0]))
            for n in free_new[lo:]:
                if n[0] >= o[1]:
                    break
                ro = _ro(o[0], o[1], n[0], n[1])
                if ro >= min_ro:
                    cand[o].append((ro, n))
        best_new = {o: max(c)[1] for o, c in cand.items()}
        best_old: dict[tuple, tuple[float, tuple]] = {}
        for o, c in cand.items():
            for ro, n in c:
                if n not in best_old or (ro, o) > best_old[n]:
                    best_old[n] = (ro, o)
        paired_old, paired_new = {}, {}
        for o, n in best_new.items():
            if best_old[n][1] == o:
                ro = best_old[n][0]
                paired_old[o], paired_new[n] = n, o
                pairs.append((o, n, "moved", ro))
        # 3. leftovers whose best match is already paired: split / merged
        split, merged = set(), set()
        for o, n in best_new.items():
            if o not in paired_old and n in paired_new:
                merged.add(n)
                pairs.append((o, n, "merged", max(c for c in cand[o])[0]))
        for n, (ro, o) in best_old.items():
            if n not in paired_new and o in paired_old:
                split.add(o)
                pairs.append((o, n, "split", ro))
        matched_old = {p[0] for p in pairs}
        matched_new = {p[1] for p in pairs}
        lost += [o for o in olds if o not in matched_old]
        born += [n for n in news if n not in matched_new]
        # the primary pair of a split or merged group takes its label
        for k, (o, n, rel, ro) in enumerate(pairs):
            if rel == "moved" and (o in split or n in merged):
                pairs[k] = (o, n, "complex" if o in split and n in merged
                            else "split" if o in split else "merged", ro)
    rows, counts = [], defaultdict(int)
    for o, n, rel, ro in pairs:
        counts[rel] += 1
        rows.append({"old_site_id": o[2], "old_site_key": o[3], "new_site_id": n[2],
                     "new_site_key": n[3], "relation": rel, "overlap": f"{ro:.3f}"})
    for side, sites, rel in (("old", lost, "lost"), ("new", born, "new")):
        for s_ in sites:
            counts[rel] += 1
            rows.append({f"{side}_site_id": s_[2], f"{side}_site_key": s_[3], "relation": rel})
    return rows, dict(counts)


def write_map(rows: list[dict], out_tsv: str | Path) -> None:
    def order(r):
        k = r.get("old_site_key") or r.get("new_site_key")
        contig, rng, _ = k.split(":")
        return contig, int(rng.split("-")[0])
    with open(out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=MAP_COLUMNS, delimiter="\t", lineterminator="\n",
                           restval="")
        w.writeheader()
        w.writerows(sorted(rows, key=order))
