"""`vedette cohort-filter`: assemble per-sample genotype tables into a
multi-sample VCF with cohort metrics (AC/AF/NS, carrier-vs-caller discordance,
founder HWE, artifact flag).
"""

from __future__ import annotations

import csv
import logging
from pathlib import Path

import pysam

from .. import __version__
from ..ped import Pedigree

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


def _load_table_arrays(path: str | Path, site_idx: dict[str, int], n_sites: int):
    """One sample's genotype table as compact numpy columns (site-indexed).

    Dict-of-dict loading costs ~1-2 kB/row in Python objects — at 2676
    samples x 80k sites that is >200 GB. Arrays cost 12 bytes/site.
    """
    import numpy as np

    nalt = np.full(n_sites, -1, dtype=np.int8)   # -1 = no-call
    gq = np.zeros(n_sites, dtype=np.int16)
    cn = np.full(n_sites, -1, dtype=np.int8)
    rd = np.full(n_sites, np.nan, dtype=np.float32)
    pe = np.zeros(n_sites, dtype=np.uint16)
    sr = np.zeros(n_sites, dtype=np.uint16)
    sample = ""
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            sample = row["sample"]
            i = site_idx.get(row["site_id"])
            if i is None or row["gt"] == "./.":
                continue
            nalt[i] = sum(int(a) for a in row["gt"].split("/"))
            gq[i] = min(int(row["gq"]), 32767)
            if row["cn"] not in ("", "-1"):
                cn[i] = min(int(row["cn"]), 127)
            if row["rd"]:
                rd[i] = float(row["rd"])
            pe[i] = min(int(row["pe"]), 65535)
            sr[i] = min(int(row["sr"]), 65535)
    return sample, {"nalt": nalt, "gq": gq, "cn": cn, "rd": rd, "pe": pe, "sr": sr}


def run_cohort_filter(
    sites_vcf: str | Path,
    genotype_tsvs: list[str | Path],
    out_vcf: str | Path,
    pedigree: Pedigree | None = None,
    artifact_af: float = 0.02,
    artifact_gtr_discordance: float = 0.5,
) -> dict[str, int]:
    import numpy as np

    # site order defines the array index
    site_ids: list[str] = []
    with pysam.VariantFile(str(sites_vcf)) as sv:
        for rec in sv:
            site_ids.append(rec.id)
    site_idx = {sid: i for i, sid in enumerate(site_ids)}
    n_sites = len(site_ids)

    tables: dict[str, dict] = {}
    for path in genotype_tsvs:
        sample, arrays = _load_table_arrays(path, site_idx, n_sites)
        if not sample:
            logger.warning("Empty genotype table: %s", path)
            continue
        tables[sample] = arrays
    samples = sorted(tables)
    if not samples:
        raise ValueError("No genotype tables to assemble")
    logger.info(
        "cohort-filter: %d samples x %d sites (%.1f MB in arrays)",
        len(samples), n_sites, len(samples) * n_sites * 12 / 1e6,
    )
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

        is_founder = np.array([s in founders for s in samples])
        gt_of_nalt = {0: (0, 0), 1: (0, 1), 2: (1, 1)}
        out_path = str(out_vcf)
        mode = "wz" if out_path.endswith(".gz") else "w"
        with pysam.VariantFile(out_path, mode, header=header) as out:
            for i, rec in enumerate(sites):
                counts["sites"] += 1
                new = out.new_record(
                    contig=rec.contig, start=rec.start, stop=rec.stop,
                    alleles=rec.alleles, id=rec.id,
                )
                for k, v in rec.info.items():
                    new.info[k] = v

                col_nalt = np.array([tables[s]["nalt"][i] for s in samples])
                called = col_nalt >= 0
                for s_i, s in enumerate(samples):
                    fmt = new.samples[s]
                    t = tables[s]
                    if not called[s_i]:
                        fmt["GT"] = (None, None)
                        continue
                    fmt["GT"] = gt_of_nalt[int(col_nalt[s_i])]
                    fmt["GQ"] = int(t["gq"][i])
                    if t["cn"][i] >= 0:
                        fmt["CN"] = int(t["cn"][i])
                    if not np.isnan(t["rd"][i]):
                        fmt["RD"] = float(t["rd"][i])
                    fmt["PE"] = int(t["pe"][i])
                    fmt["SR"] = int(t["sr"][i])

                ns = int(called.sum())
                an = 2 * ns
                ac = int(col_nalt[called].sum())
                ncar = int((col_nalt > 0).sum())
                f = is_founder & called
                het = int((col_nalt[f] == 1).sum())
                hom_alt = int((col_nalt[f] == 2).sum())
                hom_ref = int((col_nalt[f] == 0).sum())

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
