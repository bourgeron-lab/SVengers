"""`vedette cohort-filter`: assemble per-sample genotype tables into a
multi-sample VCF with cohort metrics (AC/AF/NS, carrier-vs-caller discordance,
founder HWE, artifact flag).

Post-processing flags added on IGV feedback (family FAM-04, 2026-09-24:
10/13 PASS sites false, all carried by 8-9 members out of 9 and all
CohortArtifact at 46-100 % carriers in the 2670-sample cohort):

- ``HighCarrierFrac``: carried by >= ``carrier_frac_max`` of the genotyped
  samples — reference/mappability artifacts that HWE cannot see (everyone 1/1
  is in equilibrium) and that small family runs have no frequency to catch;
- ``CNDiscordant``: the re-genotyped copy number of most carriers contradicts
  the SV type (DUP with CN <= 1, DEL with CN >= 3);
- ``CommonRef``: with ``--known-sites`` (a reference cohort VCF, e.g. the
  2670-sample run), the matching site (same type, RO >= 0.5) is carried by
  >= ``known_common_frac`` of that cohort. REF_AF / REF_CARRIER_FRAC /
  REF_FILTER are annotated whenever a match exists.
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


class _KnownSites:
    """Reference-cohort lookup: best reciprocal-overlap match of the same type."""

    def __init__(self, path: str | Path, min_ro: float = 0.5):
        # INFO and FILTER only: decoding thousands of samples per record made
        # each lookup cost seconds on a 2670-sample reference
        self.vcf = pysam.VariantFile(str(path), drop_samples=True)
        self.min_ro = min_ro
        self._contigs = set(self.vcf.header.contigs)

    def lookup(self, contig: str, start: int, stop: int, svtype: str):
        if contig not in self._contigs:
            return None
        size = max(stop - start, 1)
        best, best_ro = None, 0.0
        try:
            it = self.vcf.fetch(contig, max(0, start - size), stop + size)
        except ValueError:
            return None
        for rec in it:
            if rec.info.get("SVTYPE") != svtype:
                continue
            ov = min(rec.stop, stop) - max(rec.start, start)
            if ov <= 0:
                continue
            ro = ov / max(rec.stop - rec.start, size)
            if ro >= self.min_ro and ro > best_ro:
                best, best_ro = rec, ro
        if best is None:
            return None
        af = best.info.get("AF", (0.0,))
        af = af[0] if isinstance(af, tuple) else af
        ns = int(best.info.get("NS", 0))
        ncar = int(best.info.get("NCARRIER", 0))
        return {
            "af": float(af),
            "carrier_frac": ncar / ns if ns else 0.0,
            # ';' separates INFO fields: joining on it would truncate REF_FILTER
            # to its first value and leak the others as undeclared INFO keys
            "filter": "|".join(best.filter.keys()) or "PASS",
        }


def run_cohort_filter(
    sites_vcf: str | Path,
    genotype_tsvs: list[str | Path],
    out_vcf: str | Path,
    pedigree: Pedigree | None = None,
    artifact_af: float = 0.02,
    artifact_gtr_discordance: float = 0.5,
    carrier_frac_max: float = 0.85,
    cn_discordant_frac: float = 0.5,
    known_sites: str | Path | None = None,
    known_common_frac: float = 0.5,
) -> dict[str, int]:
    import numpy as np

    known = _KnownSites(known_sites) if known_sites else None

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

    counts = {"sites": 0, "artifact": 0, "high_carrier": 0, "cn_discordant": 0,
              "common_ref": 0, "pass": 0}
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
             "Absolute discordance between nominating-sample count and re-genotyped carriers "
             "(only the 'fewer carriers than nominations' direction is flagged)"),
            ("HWE_P", "1", "Float", "Founder Hardy-Weinberg exact p-value"),
            ("CARRIER_FRAC", "1", "Float", "Fraction of genotyped samples carrying the variant"),
            ("CN_DISC", "1", "Float",
             "Fraction of carriers whose copy number contradicts SVTYPE (DUP with CN<=1, DEL with CN>=3)"),
            ("REF_AF", "1", "Float", "Allele frequency of the matching site in the reference cohort (--known-sites)"),
            ("REF_CARRIER_FRAC", "1", "Float", "Carrier fraction of the matching site in the reference cohort"),
            ("REF_FILTER", "1", "String", "FILTER of the matching site in the reference cohort, '|'-separated"),
        ):
            if line[0] not in header.info:
                header.info.add(*line)
        for name, desc in (
            ("CohortArtifact", "High-frequency or genotype-discordant recurrent locus"),
            ("HighCarrierFrac",
             f"Carried by >= {carrier_frac_max:g} of the genotyped samples (reference/mappability artifact)"),
            ("CNDiscordant",
             f"Copy number of >= {cn_discordant_frac:g} of the carriers contradicts the SV type"),
            ("CommonRef",
             f"Matching site carried by >= {known_common_frac:g} of the reference cohort (--known-sites)"),
        ):
            if name not in header.filters:
                header.filters.add(name, None, None, desc)
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
                col_cn = np.array([tables[s]["cn"][i] for s in samples])
                called = col_nalt >= 0
                svtype = rec.info.get("SVTYPE", "")
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
                # Directional: only "re-genotyping reproduces FEWER carriers
                # than the callers nominated" is suspicious. The other
                # direction is re-genotyping being more sensitive than
                # per-sample calling, which is the whole point of genotyping
                # from the stores — on the PMS cohort it filtered the causal
                # SHANK3 deletion (6 nominators, 64 carriers) and 36328 sites
                # in that direction vs 425 in the other. The artifact case it
                # used to cover ("the genotyper calls everyone") is now
                # HighCarrierFrac's job.
                disc_down = (nominators - ncar) / nominators if nominators else 0.0
                new.info["GTR_DISC"] = round(disc, 4)
                hwe_p = hwe_exact(het, hom_alt, hom_ref)
                new.info["HWE_P"] = round(hwe_p, 6)

                af = ac / an if an else 0.0
                tags = []
                if (af > artifact_af and hwe_p < 1e-6) or (
                    nominators >= 3 and disc_down > artifact_gtr_discordance
                ):
                    tags.append("CohortArtifact")

                carrier_frac = ncar / ns if ns else 0.0
                new.info["CARRIER_FRAC"] = round(carrier_frac, 4)
                if ns >= 4 and carrier_frac >= carrier_frac_max:
                    tags.append("HighCarrierFrac")

                carrier_cn = col_cn[(col_nalt > 0) & (col_cn >= 0)]
                if len(carrier_cn) >= 2:
                    contra = (carrier_cn <= 1) if svtype == "DUP" else (carrier_cn >= 3)
                    cn_disc = float(contra.mean())
                    new.info["CN_DISC"] = round(cn_disc, 4)
                    if cn_disc >= cn_discordant_frac:
                        tags.append("CNDiscordant")

                if known is not None:
                    hit = known.lookup(rec.contig, rec.start, rec.stop, svtype)
                    if hit is not None:
                        new.info["REF_AF"] = round(hit["af"], 6)
                        new.info["REF_CARRIER_FRAC"] = round(hit["carrier_frac"], 4)
                        new.info["REF_FILTER"] = hit["filter"]
                        if hit["carrier_frac"] >= known_common_frac:
                            tags.append("CommonRef")

                for tag in tags or ["PASS"]:
                    new.filter.add(tag)
                counts["artifact"] += "CohortArtifact" in tags
                counts["high_carrier"] += "HighCarrierFrac" in tags
                counts["cn_discordant"] += "CNDiscordant" in tags
                counts["common_ref"] += "CommonRef" in tags
                counts["pass"] += not tags
                out.write(new)
    if mode == "wz":
        pysam.tabix_index(out_path, preset="vcf", force=True)
    logger.info(
        "cohort-filter: %d sites, %d PASS, %d CohortArtifact, %d HighCarrierFrac, "
        "%d CNDiscordant, %d CommonRef, %d samples",
        counts["sites"], counts["pass"], counts["artifact"], counts["high_carrier"],
        counts["cn_discordant"], counts["common_ref"], len(samples),
    )
    return counts
