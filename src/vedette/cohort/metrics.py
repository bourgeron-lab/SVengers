"""`vedette cohort-filter`: assemble per-sample genotype tables into a
multi-sample VCF with cohort metrics (AC/AF/NS, carrier-vs-caller discordance,
founder HWE, artifact flag).
"""

from __future__ import annotations

import csv
import logging
import math
from pathlib import Path

import pysam

from ..ped import Pedigree
from .. import __version__

logger = logging.getLogger(__name__)


def hwe_exact(het: int, hom_alt: int, hom_ref: int) -> float:
    """Wigginton exact Hardy-Weinberg p-value (SNP-HWE)."""
    n = het + hom_alt + hom_ref
    if n == 0:
        return 1.0
    rare = 2 * min(hom_alt, hom_ref) + het
    probs = [0.0] * (rare + 1)
    mid = rare * (2 * n - rare) // (2 * n)
    if (mid % 2) != (rare % 2):
        mid += 1
    probs[mid] = 1.0
    total = 1.0
    curr_het, curr_hom_r = mid, (rare - mid) // 2
    curr_hom_c = n - curr_het - curr_hom_r
    h, hr, hc = curr_het, curr_hom_r, curr_hom_c
    while h > 1:
        probs[h - 2] = probs[h] * h * (h - 1) / (4.0 * (hr + 1) * (hc + 1))
        total += probs[h - 2]
        h, hr, hc = h - 2, hr + 1, hc + 1
    h, hr, hc = curr_het, curr_hom_r, curr_hom_c
    while h <= rare - 2:
        probs[h + 2] = probs[h] * 4.0 * hr * hc / ((h + 2) * (h + 1))
        total += probs[h + 2]
        h, hr, hc = h + 2, hr - 1, hc - 1
    obs = probs[het] if het <= rare else 0.0
    p = sum(pr for pr in probs if pr <= obs * 1.0000001) / total
    return min(p, 1.0)


def _read_genotype_table(path: str | Path) -> dict[str, dict]:
    """{site_id: row} for one sample."""
    rows: dict[str, dict] = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            rows[row["site_id"]] = row
    return rows


def run_cohort_filter(
    sites_vcf: str | Path,
    genotype_tsvs: list[str | Path],
    out_vcf: str | Path,
    pedigree: Pedigree | None = None,
    artifact_af: float = 0.02,
    artifact_gtr_discordance: float = 0.5,
) -> dict[str, int]:
    tables = {}
    for path in genotype_tsvs:
        table = _read_genotype_table(path)
        if not table:
            logger.warning("Empty genotype table: %s", path)
            continue
        sample = next(iter(table.values()))["sample"]
        tables[sample] = table
    samples = sorted(tables)
    if not samples:
        raise ValueError("No genotype tables to assemble")
    founders = (
        {ind.sample_id for ind in pedigree.founders()} if pedigree else set(samples)
    )

    counts = {"sites": 0, "artifact": 0}
    with pysam.VariantFile(str(sites_vcf)) as sites:
        header = sites.header.copy()
        header.add_line(f"##source=vedette-cohort-{__version__}")
        for line in (
            ("AC", "A", "Integer", "Alternate allele count"),
            ("AN", "1", "Integer", "Total genotyped alleles"),
            ("AF", "A", "Float", "Alternate allele frequency"),
            ("NS", "1", "Integer", "Samples with a genotype"),
            ("NCARRIER", "1", "Integer", "Samples carrying the variant after re-genotyping"),
            ("GTR_DISC", "1", "Float",
             "Discordance between nominating-sample count and re-genotyped carriers"),
            ("HWE_P", "1", "Float", "Founder Hardy-Weinberg exact p-value"),
        ):
            if line[0] not in header.info:
                header.info.add(*line)
        if "CohortArtifact" not in header.filters:
            header.filters.add(
                "CohortArtifact", None, None,
                "High-frequency or genotype-discordant recurrent locus",
            )
        for line in (
            ("GT", "1", "String", "Genotype"),
            ("GQ", "1", "Integer", "Genotype quality"),
            ("CN", "1", "Integer", "Integer copy number"),
            ("RD", "1", "Float", "Mean copy ratio"),
            ("PE", "1", "Integer", "Discordant pairs"),
            ("SR", "1", "Integer", "Split reads"),
        ):
            if line[0] not in header.formats:
                header.formats.add(*line)
        for s in samples:
            header.add_sample(s)

        out_path = str(out_vcf)
        mode = "wz" if out_path.endswith(".gz") else "w"
        with pysam.VariantFile(out_path, mode, header=header) as out:
            for rec in sites:
                counts["sites"] += 1
                new = out.new_record(
                    contig=rec.contig, start=rec.start, stop=rec.stop,
                    alleles=rec.alleles, id=rec.id,
                )
                for k, v in rec.info.items():
                    new.info[k] = v

                ac = an = ns = ncar = 0
                het = hom_alt = hom_ref = 0
                for s in samples:
                    row = tables[s].get(rec.id)
                    fmt = new.samples[s]
                    if row is None or row["gt"] == "./.":
                        fmt["GT"] = (None, None)
                        continue
                    alleles = tuple(int(a) for a in row["gt"].split("/"))
                    fmt["GT"] = alleles
                    fmt["GQ"] = int(row["gq"])
                    if row["cn"] not in ("", "-1"):
                        fmt["CN"] = int(row["cn"])
                    if row["rd"]:
                        fmt["RD"] = float(row["rd"])
                    fmt["PE"] = int(row["pe"])
                    fmt["SR"] = int(row["sr"])
                    ns += 1
                    an += 2
                    n_alt = sum(alleles)
                    ac += n_alt
                    if n_alt:
                        ncar += 1
                    if s in founders:
                        if n_alt == 1:
                            het += 1
                        elif n_alt == 2:
                            hom_alt += 1
                        else:
                            hom_ref += 1

                new.info["AC"] = (ac,)
                new.info["AN"] = an
                new.info["AF"] = (round(ac / an, 6) if an else 0.0,)
                new.info["NS"] = ns
                new.info["NCARRIER"] = ncar
                nominators = int(rec.info.get("NSAMPLES", 0))
                disc = (
                    abs(nominators - ncar) / max(nominators, 1) if nominators else 0.0
                )
                new.info["GTR_DISC"] = round(disc, 4)
                hwe_p = hwe_exact(het, hom_alt, hom_ref)
                new.info["HWE_P"] = round(hwe_p, 6)

                af = ac / an if an else 0.0
                is_artifact = (
                    (af > artifact_af and hwe_p < 1e-6)
                    or (nominators >= 3 and disc > artifact_gtr_discordance)
                )
                if is_artifact:
                    new.filter.add("CohortArtifact")
                    counts["artifact"] += 1
                else:
                    new.filter.add("PASS")
                out.write(new)
    if mode == "wz":
        pysam.tabix_index(out_path, preset="vcf", force=True)
    logger.info(
        "cohort-filter: %d sites, %d flagged CohortArtifact, %d samples",
        counts["sites"], counts["artifact"], len(samples),
    )
    return counts
