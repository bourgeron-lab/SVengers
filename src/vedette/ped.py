"""Pedigree parsing and family structure.

Handles standard 6-column PED/TSV files, with or without a header line
(sniffed). Column semantics: family_id, sample_id, father_id, mother_id,
sex (1=male, 2=female, 0=unknown), phenotype (kept verbatim). Missing
parents are encoded as "0", "-9" or "".

Mendelian-consistency and de novo scoring live in the cohort stage; this
module only builds the structure they operate on.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path

logger = logging.getLogger(__name__)

_MISSING = {"0", "-9", "", "NA", "."}


@dataclass(frozen=True)
class Individual:
    family_id: str
    sample_id: str
    father_id: str | None
    mother_id: str | None
    sex: int          # 1 male, 2 female, 0 unknown
    phenotype: str


@dataclass
class Pedigree:
    individuals: dict[str, Individual] = field(default_factory=dict)

    @property
    def families(self) -> dict[str, list[Individual]]:
        fams: dict[str, list[Individual]] = {}
        for ind in self.individuals.values():
            fams.setdefault(ind.family_id, []).append(ind)
        return fams

    def trios(self) -> list[tuple[Individual, Individual, Individual]]:
        """(child, father, mother) triples where both parents are in the pedigree."""
        out = []
        for ind in self.individuals.values():
            if (
                ind.father_id
                and ind.mother_id
                and ind.father_id in self.individuals
                and ind.mother_id in self.individuals
            ):
                out.append(
                    (ind, self.individuals[ind.father_id], self.individuals[ind.mother_id])
                )
        return out

    def founders(self) -> list[Individual]:
        """Individuals with no parent present in the pedigree (for HWE/AF)."""
        return [
            ind
            for ind in self.individuals.values()
            if (not ind.father_id or ind.father_id not in self.individuals)
            and (not ind.mother_id or ind.mother_id not in self.individuals)
        ]

    def sex_of(self, sample_id: str) -> int:
        ind = self.individuals.get(sample_id)
        return ind.sex if ind else 0


def _is_header(fields: list[str]) -> bool:
    lowered = [f.lower() for f in fields]
    return any(tok in lowered for tok in ("fid", "sample_id", "iid", "fatherbarcode"))


def parse_ped(path: str | Path) -> Pedigree:
    """Parse a 6+ column pedigree TSV, sniffing an optional header line."""
    ped = Pedigree()
    with open(path) as fh:
        lines = [ln.rstrip("\n") for ln in fh if ln.strip()]
    if not lines:
        raise ValueError(f"Empty pedigree file: {path}")

    first = lines[0].split("\t")
    if len(first) < 6:
        raise ValueError(
            f"{path}: expected >=6 tab-separated columns, got {len(first)} in first line"
        )
    body = lines[1:] if _is_header(first) else lines

    for i, line in enumerate(body, start=2 if _is_header(first) else 1):
        fields = line.split("\t")
        if len(fields) < 6:
            raise ValueError(f"{path}:{i}: expected >=6 columns, got {len(fields)}")
        fam, sid, father, mother, sex, pheno = (f.strip() for f in fields[:6])
        if sid in ped.individuals:
            raise ValueError(f"{path}:{i}: duplicate sample_id {sid!r}")
        try:
            sex_code = int(sex)
        except ValueError:
            sex_code = 0
        if sex_code not in (0, 1, 2):
            sex_code = 0
        ped.individuals[sid] = Individual(
            family_id=fam,
            sample_id=sid,
            father_id=None if father in _MISSING else father,
            mother_id=None if mother in _MISSING else mother,
            sex=sex_code,
            phenotype=pheno,
        )

    # referenced parents absent from the file: warn once, keep the reference
    missing = {
        pid
        for ind in ped.individuals.values()
        for pid in (ind.father_id, ind.mother_id)
        if pid and pid not in ped.individuals
    }
    if missing:
        logger.warning(
            "%s: %d referenced parent(s) not present as individuals (e.g. %s)",
            path, len(missing), sorted(missing)[0],
        )
    logger.info(
        "Pedigree %s: %d individuals, %d families, %d trios",
        Path(path).name, len(ped.individuals), len(ped.families), len(ped.trios()),
    )
    return ped
