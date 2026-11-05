"""ClinGen dosage sensitivity map for the impact score.

Reads the two ClinGen curation lists (GRCh38), as published on
ftp://ftp.clinicalgenome.org:

- ``ClinGen_region_curation_list_GRCh38.tsv``: recurrent and other regions
  (15q11.2q13 PWS/AS, 16p11.2, 22q11.2, 7q11.23, Xq28 MECP2...), with a
  haploinsufficiency (HI) and a triplosensitivity (TS) score;
- ``ClinGen_gene_curation_list_GRCh38.tsv``: the same scores per gene.

Scores: 3 sufficient evidence, 2 some, 1 little, 0 none, 30 autosomal
recessive, 40 dosage sensitivity unlikely. "Not yet evaluated", "-1" and
locations still "tbd" are read as missing.

Why regions matter: GeneTrek has no triplosensitivity, and a recurrent
syndrome region is curated as a region, not through one gene. In the v5 IGV
review, a de novo 15q11.2q13 duplication at copy number 4 in an affected
child (ISCA-37478, TS 3) was on its family page with an impact of +0.15
("TS non évaluée"), below hundreds of rows.
"""

from __future__ import annotations

import bisect
import csv
from dataclasses import dataclass
from pathlib import Path

# a recurrent region is bounded by segmental duplications, where depth calls
# stop short: covering 90 % of it counts as complete overlap (ClinGen 2A)
COMPLETE_FRAC = 0.90
# partial overlap worth a reader's attention (the v5 review's 15q site
# covered 35 % of its region, the child's own call 21 %)
PARTIAL_FRAC = 0.20


@dataclass(frozen=True)
class Region:
    isca: str
    name: str
    contig: str
    start: int  # 0-based half-open
    end: int
    hi: int | None
    ts: int | None

    def short_name(self) -> str:
        """The region name without its "(includes ...)" tail, at most 60 characters."""
        name = self.name.split(" (includes")[0]
        return name if len(name) <= 60 else name[:59] + "…"


def _score(v: str) -> int | None:
    try:
        x = int(v)
    except (TypeError, ValueError):
        return None
    return x if x >= 0 else None


def _read(path: str | Path, first_col: str) -> list[dict]:
    """Rows of a ClinGen list: comment lines, then a header starting with
    ``#<first_col>``."""
    with open(path) as fh:
        lines = [ln for ln in fh if not ln.startswith("#") or ln.startswith(f"#{first_col}")]
    if not lines or not lines[0].startswith(f"#{first_col}"):
        raise ValueError(f"{path}: no '#{first_col}' header line (ClinGen curation list expected)")
    lines[0] = lines[0][1:]
    return list(csv.DictReader(lines, delimiter="\t"))


def _location(loc: str) -> tuple[str, int, int] | None:
    """'chr15:23465365-28134728' (1-based inclusive) -> ('chr15', 23465364, 28134728)."""
    try:
        contig, rng = loc.strip().split(":")
        a, b = (int(x.replace(",", "")) for x in rng.split("-"))
    except ValueError:
        return None
    return contig, a - 1, b


def _norm(contig: str) -> str:
    return contig.removeprefix("chr")


class ClinGen:
    """Region and gene dosage scores, with overlap lookups."""

    def __init__(self, regions: list[Region], genes: dict[str, tuple[int | None, int | None]]):
        self.genes = genes
        self._by_contig: dict[str, list[Region]] = {}
        for r in sorted(regions, key=lambda r: (r.contig, r.start)):
            self._by_contig.setdefault(_norm(r.contig), []).append(r)
        self._starts = {c: [r.start for r in rs] for c, rs in self._by_contig.items()}
        # longest region per contig bounds the left search
        self._maxlen = {c: max(r.end - r.start for r in rs) for c, rs in self._by_contig.items()}

    @classmethod
    def load(cls, genes_tsv: str | Path | None = None,
             regions_tsv: str | Path | None = None) -> ClinGen:
        regions = []
        if regions_tsv:
            for row in _read(regions_tsv, "ISCA ID"):
                loc = _location(row.get("Genomic Location", ""))
                if loc is None:
                    continue
                regions.append(Region(
                    row["ISCA ID"], row.get("ISCA Region Name", ""), *loc,
                    _score(row.get("Haploinsufficiency Score")),
                    _score(row.get("Triplosensitivity Score")),
                ))
        genes = {}
        if genes_tsv:
            for row in _read(genes_tsv, "Gene Symbol"):
                genes[row["Gene Symbol"]] = (_score(row.get("Haploinsufficiency Score")),
                                             _score(row.get("Triplosensitivity Score")))
        return cls(regions, genes)

    @property
    def n_regions(self) -> int:
        return sum(len(rs) for rs in self._by_contig.values())

    def gene(self, symbol: str) -> tuple[int | None, int | None] | None:
        """(HI, TS) of a gene symbol, None if ClinGen has no curation for it."""
        return self.genes.get(symbol)

    def regions(self, contig: str, start: int, end: int) -> list[tuple[Region, float, bool]]:
        """[(region, fraction of the region covered, event inside the region)]
        of the regions overlapping [start, end)."""
        c = _norm(contig)
        rs = self._by_contig.get(c)
        if not rs:
            return []
        lo = bisect.bisect_left(self._starts[c], start - self._maxlen[c])
        hi = bisect.bisect_left(self._starts[c], end)
        out = []
        for r in rs[lo:hi]:
            ov = min(end, r.end) - max(start, r.start)
            if ov > 0:
                out.append((r, ov / (r.end - r.start), start >= r.start and end <= r.end))
        return out


def region_points(clingen: ClinGen | None, contig: str, start: int, end: int,
                  dup: bool) -> tuple[float, str, str]:
    """(points, rule, benign rule) of the ClinGen regions a CNV overlaps, with
    the HI score for a loss and the TS score for a gain (Riggs et al. 2020,
    section 2):

    - established region (score 3) covered >= COMPLETE_FRAC: +1.00 (2A);
    - established region covered >= PARTIAL_FRAC: +0.45. ClinGen scores a
      partial overlap 0 pending review (2B); the triage score keeps it
      visible, since depth calls stop inside the flanking repeats;
    - benign rule (else empty): the CNV lies inside a region curated "dosage sensitivity
      unlikely" (score 40) for its type (2F for losses, 2C-2D for gains).
    """
    if clingen is None:
        return 0.0, "", ""
    best, rule, benign = 0.0, "", ""
    for r, frac, inside in clingen.regions(contig, start, end):
        score = r.ts if dup else r.hi
        if score == 40 and inside and not benign:
            benign = (f"dans la région ClinGen {r.short_name()} "
                      f"({r.isca}, {'TS' if dup else 'HI'} 40)")
        if score != 3 or frac < PARTIAL_FRAC:
            continue
        pts = 1.00 if frac >= COMPLETE_FRAC else 0.45
        if pts > best:
            kind = "complète" if pts == 1.00 else f"partielle {frac:.0%}".replace("%", " %")
            best = pts
            rule = (f"{'DUP' if dup else 'DEL'} {kind} de la région ClinGen "
                    f"{r.short_name()} ({r.isca}, {'TS' if dup else 'HI'} 3)")
    return best, rule, benign
