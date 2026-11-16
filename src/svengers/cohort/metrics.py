"""`svengers cohort-filter`: assemble per-sample genotype tables into a
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
- ``CohortArtifact`` from the reference (``known_artifacts``, INFO
  REF_ARTIFACT): a site common here (AF > ``artifact_af``) whose matching
  reference site is CohortArtifact. A small cohort has neither the power nor
  the stability of the founders' HWE test: PMS (206 founders) flags 0.46 % of
  its HWE failures differently whenever one founder is left out, against
  0.04 % on dataset_paper (1,688), and of PMS's 54,492 common sites matched
  in dataset_paper, dataset_paper flags 10,767 that PMS does not (76 the
  other way). Rare sites keep their own verdict.

On chrX / chrY a man carries one copy (his genotypes are 0/0 or 1/1, his
reference CN 1): HWE is tested on the female founders only (none on chrY),
HighCarrierFrac on each sex's carrier fraction, and CNDiscordant against
each carrier's own ploidy. Each sample's chrX ploidy is read from the copy
number of its reference genotypes on chrX — the stores' ploidy, which a
pedigree sex can contradict. Pooled, the two sexes hid chrX artifacts:
a junction artifact carried by 96 % of the women and half of the men sat
at 0.6-0.75 of the cohort, below ``carrier_frac_max``.

- ``ChrYNoDepth``: on chrY outside the PARs, the copy number of fewer than
  half of the carriers supports the SV type (DUP above, DEL below the
  carrier's ploidy). dataset_paper: 453 of 525 chrY PASS sites were small
  junction-only calls in the Y repeats, carried by up to two thirds of the
  men at their baseline copy number, and a son disagreed with his father in
  40-54 % of the father/son pairs carrying them (20 % on the 70 deletions
  the copy number supports). Sites shorter than one depth bin
  (``CHRY_MIN_SIZE``) are flagged too, their carriers' CN being set by the
  junctions: on the 72 sites left, father/son discordance 34 % (DEL) and
  86 % (DUP) below 1 kb, 10 % and 24 % above.

Two trio flags (with a pedigree, autosomes only; a trio counts at a site when
the three are genotyped with GQ >= ``trio_min_gq``):

- ``LowTransmission``: among the trios with one heterozygous parent and the
  other 0/0 (INFO TR_OPP), the child carries in < ``low_transmission`` of
  them (INFO TR_RATE), with >= ``trio_min_opp`` such trios. A true variant
  is transmitted half of the time; parents 'carrying' a noise genotype
  transmit nothing.
- ``DeNovoExcess``: >= ``denovo_excess_n`` carrier children have both
  parents 0/0 and a copy ratio that does not support the SV type (DEL with
  RD > ``denovo_rd_del``, DUP with RD < ``denovo_rd_dup``; INFO DN_NODEPTH),
  and they are >= ``denovo_excess_frac`` of the carrier children (DN_FRAC
  counts every de novo child). A site where children keep carrying what
  neither parent has, without the depth of the SV, is a genotyping
  artifact of the children's genomes: on dataset_paper, the de novo
  children of such DEL sites had a median RD of 0.91, 74 % of them in wavy
  genomes, their parents 0.97. The depth condition keeps the true recurrent
  de novo: the 22q13.3 terminal deletion (3 unrelated children, RD 0.48-0.50)
  and a 74 kb deletion at 11q13.3 (3 children, RD 0.50-0.51) were flagged
  without it, and in a Phelan-McDermid cohort most probands carry a de novo
  SHANK3 deletion.

dataset_paper (665 trios), model v14b, PASS sites >= 1 kb: the DUP >= 1 kb
the v12b model did not have were transmitted 35 % of the time at 1-10 kb and
41 % at >= 10 kb, vs 46 % for the DEL and 44-45 % for the shared DUP; 27 %
of the new 1-10 kb DUP sites had a site rate < 0.2 (8-10 % of the DEL), and
their carrier children had both parents 0/0 twice as often as the DEL's.
"""

from __future__ import annotations

import csv
import logging
from pathlib import Path

import pysam

from .. import __version__
from ..ped import Pedigree
from ..vcfout import span, spec_coords
from .mendel import site_mode

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
        self._spec = spec_coords(self.vcf.header)

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
            r_start, r_end = span(rec, self._spec)
            ov = min(r_end, stop) - max(r_start, start)
            if ov <= 0:
                continue
            ro = ov / max(r_end - r_start, size)
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


# set by run_cohort_filter before its pool forks: shared copy-on-write
_SHARED: dict = {}


SEX_CONTIGS = ("chrX", "X", "chrY", "Y")
CHRY_MIN_SIZE = 1000   # one bin of the depth tracks (fast_bin_size)


def hemizygous_x(tables: dict, samples: list[str], x_sites):
    """Per sample (in ``samples`` order): one chrX copy, from the median copy
    number of its 0/0 genotypes on the chrX sites (``x_sites``: indices)."""
    import numpy as np

    out = np.zeros(len(samples), dtype=bool)
    if len(x_sites) == 0:
        return out
    for k, s in enumerate(samples):
        t = tables[s]
        cn = t["cn"][x_sites]
        ref = cn[(t["nalt"][x_sites] == 0) & (cn >= 0)]
        out[k] = len(ref) >= 10 and float(np.median(ref)) < 1.5
    return out


def _load_one(path: str | Path):
    return _load_table_arrays(path, _SHARED["site_idx"], _SHARED["n_sites"])


TRIO_DEFAULTS = {
    "trio_min_gq": 20, "trio_min_opp": 10, "low_transmission": 0.2,
    "denovo_excess_n": 3, "denovo_excess_frac": 0.5,
    "denovo_rd_del": 0.75, "denovo_rd_dup": 1.25,
    # FORMAT/FT=WeakSupport (weak_genotype): no discordant pair and a copy
    # ratio that barely moves; the trio tools read such a genotype as missing
    "weak_min_size": 1000, "weak_large_size": 10_000,
    "weak_del_rd": 0.75, "weak_del_rd_large": 0.62, "weak_dup_rd": 1.25,
}


def weak_genotype(svtype: str, length: int, pe, rd, t: dict | None = None):
    """True where a non-reference genotype has no support of its own: no
    discordant pair, and a copy ratio that does not show the SV type (DEL
    above weak_del_rd, weak_del_rd_large from weak_large_size; DUP below
    weak_dup_rd), from weak_min_size (below one depth bin the ratio says
    nothing). Vectorised over samples (pe, rd arrays; NaN rd is not weak).

    Calibrated on dataset_paper (665 trios, PASS sites, GQ >= 20): de novo
    genotypes of children are split-read-only calls on a copy ratio near 1
    (DEL 1-10 kb: PE 0, RD 0.74-0.78, GQ ~48; inherited ones PE ~10, RD 0.6,
    GQ 99), concentrated on some flowcells (three children of one flowcell
    carry 1,600-3,200 each). The rule removes 36-72 % of the de novo DEL
    1-10 kb (35 % in wave >= 0.20) for 2-8 % of the inherited genotypes, and
    halves the de novo DEL >= 10 kb for 4-6 %; DUP >= 1 kb de novo -49 to -80 %."""
    import numpy as np

    t = {**TRIO_DEFAULTS, **(t or {})}
    pe = np.asarray(pe)
    rd = np.asarray(rd, dtype=float)
    if length < t["weak_min_size"] or svtype not in ("DEL", "DUP"):
        return np.zeros(pe.shape, dtype=bool)
    if svtype == "DEL":
        thr = t["weak_del_rd_large"] if length >= t["weak_large_size"] else t["weak_del_rd"]
        off = rd > thr
    else:
        off = rd < t["weak_dup_rd"]
    return (pe == 0) & off & ~np.isnan(rd)


def _cohort_header(sites_vcf, samples, carrier_frac_max, cn_discordant_frac,
                   known_common_frac, trio: dict | None = None) -> pysam.VariantHeader:
    t = {**TRIO_DEFAULTS, **(trio or {})}
    with pysam.VariantFile(str(sites_vcf)) as sites:
        header = sites.header.copy()
    header.add_line(f"##source=svengers-cohort-{__version__}")
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
        ("REF_ARTIFACT", "0", "Flag", "CohortArtifact taken from the matching reference site (common site, --known-sites)"),
        ("TR_OPP", "1", "Integer",
         f"Trios (GQ >= {t['trio_min_gq']}) with one heterozygous parent and the other 0/0 (autosomes)"),
        ("TR_RATE", "1", "Float", "Fraction of the TR_OPP trios whose child carries the variant"),
        ("DN_TRIO", "1", "Integer",
         f"Carrier children (trio GQ >= {t['trio_min_gq']}) whose parents are both 0/0 (autosomes)"),
        ("DN_FRAC", "1", "Float", "DN_TRIO over the carrier children of the trios"),
        ("DN_NODEPTH", "1", "Integer",
         f"DN_TRIO children whose RD does not support the SV type (DEL RD > {t['denovo_rd_del']:g}, "
         f"DUP RD < {t['denovo_rd_dup']:g})"),
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
        ("ChrYNoDepth",
         "chrY outside the PARs: shorter than one depth bin, or the copy number of fewer "
         "than half of the carriers supports the SV type (junction-only call in the Y repeats)"),
        ("LowTransmission",
         f"Transmitted to < {t['low_transmission']:g} of the children of a heterozygous x 0/0 "
         f"couple, over >= {t['trio_min_opp']} such trios (TR_OPP, TR_RATE)"),
        ("DeNovoExcess",
         f">= {t['denovo_excess_n']} carrier children with both parents 0/0 and no depth support, >= "
         f"{t['denovo_excess_frac']:g} of the carrier children (DN_NODEPTH)"),
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
        ("FT", "1", "String",
         f"WeakSupport: non-reference genotype >= {t['weak_min_size']} bp with no discordant "
         f"pair and a copy ratio that does not show the SV type (DEL RD > {t['weak_del_rd']:g}, "
         f"> {t['weak_del_rd_large']:g} from {t['weak_large_size']} bp; DUP RD < "
         f"{t['weak_dup_rd']:g}); read as missing by the trio tools"),
    ):
        if line[0] not in header.formats:
            header.formats.add(*line)
    for s in samples:
        header.add_sample(s)
    return header


def _write_chunk(task: tuple[int, int, str]) -> dict[str, int]:
    """Sites [i0, i1) of the sites VCF as cohort records, to ``part`` (BGZF
    or plain text, no header); returns the FILTER counts."""
    import numpy as np
    from pysam.libcbgzf import BGZFile

    i0, i1, part = task
    sh = _SHARED
    tables, samples, is_founder = sh["tables"], sh["samples"], sh["is_founder"]
    male = sh["male"]
    artifact_af = sh["artifact_af"]
    artifact_gtr_discordance = sh["artifact_gtr_discordance"]
    carrier_frac_max = sh["carrier_frac_max"]
    cn_discordant_frac = sh["cn_discordant_frac"]
    known_common_frac = sh["known_common_frac"]
    known_artifacts = sh.get("known_artifacts", True)
    trio = {**TRIO_DEFAULTS, **sh.get("trio", {})}
    tc, tf, tm = sh.get("trio_idx", (np.zeros(0, int),) * 3)
    header = _cohort_header(
        sh["sites_vcf"], samples, carrier_frac_max, cn_discordant_frac,
        known_common_frac, trio,
    )
    known = _KnownSites(sh["known_sites"]) if sh["known_sites"] else None
    gt_of_nalt = {0: (0, 0), 1: (0, 1), 2: (1, 1)}
    counts = {"sites": 0, "artifact": 0, "ref_artifact": 0, "high_carrier": 0, "cn_discordant": 0,
              "common_ref": 0, "chry_no_depth": 0, "low_transmission": 0,
              "denovo_excess": 0, "pass": 0}
    fo = BGZFile(part, "wb") if part.endswith(".gz") else open(part, "wb")
    with pysam.VariantFile(str(sh["sites_vcf"])) as sites, fo:
        spec = spec_coords(sites.header)
        for i, rec in enumerate(sites):
            if i < i0:
                continue
            if i >= i1:
                break
            counts["sites"] += 1
            start, end = span(rec, spec)
            new = header.new_record(
                contig=rec.contig, start=rec.start, stop=rec.stop,
                alleles=rec.alleles, id=rec.id,
            )
            for k, v in rec.info.items():
                new.info[k] = v

            col_nalt = np.array([tables[s]["nalt"][i] for s in samples])
            col_cn = np.array([tables[s]["cn"][i] for s in samples])
            col_gq = np.array([tables[s]["gq"][i] for s in samples])
            col_rd = np.array([tables[s]["rd"][i] for s in samples])
            called = col_nalt >= 0
            svtype = rec.info.get("SVTYPE", "")
            sex_contig = rec.contig in SEX_CONTIGS
            ploidy = np.where(male, 1, 2) if sex_contig else np.full(len(samples), 2)
            col_pe = np.array([tables[s]["pe"][i] for s in samples])
            weak = (col_nalt > 0) & weak_genotype(svtype, end - start, col_pe, col_rd, trio)
            # FT on every sample of a record that has one (a string FORMAT left
            # unset next to set ones reads back as padding bytes)
            any_weak = bool(weak.any())
            for s_i, s in enumerate(samples):
                fmt = new.samples[s]
                t = tables[s]
                if any_weak:
                    fmt["FT"] = "WeakSupport" if weak[s_i] else ("PASS" if called[s_i] else ".")
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
            if sex_contig:
                f &= ~male    # women only: a man is never heterozygous
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
            hwe_p = 1.0 if rec.contig in ("chrY", "Y") else hwe_exact(het, hom_alt, hom_ref)
            new.info["HWE_P"] = round(hwe_p, 6)

            af = ac / an if an else 0.0
            tags = []
            if (af > artifact_af and hwe_p < 1e-6) or (
                nominators >= 3 and disc_down > artifact_gtr_discordance
            ):
                tags.append("CohortArtifact")

            carrier_frac = ncar / ns if ns else 0.0
            new.info["CARRIER_FRAC"] = round(carrier_frac, 4)
            fracs = [carrier_frac] if ns >= 4 else []
            if sex_contig:   # each sex on its own
                fracs = []
                for grp in (male, ~male):
                    n_grp = int((called & grp).sum())
                    if n_grp >= 4:
                        fracs.append(int(((col_nalt > 0) & grp).sum()) / n_grp)
            if fracs and max(fracs) >= carrier_frac_max:
                tags.append("HighCarrierFrac")

            sel = (col_nalt > 0) & (col_cn >= 0)
            carrier_cn, carrier_pl = col_cn[sel], ploidy[sel]
            if len(carrier_cn) >= 2:
                contra = (
                    (carrier_cn < carrier_pl) if svtype == "DUP" else (carrier_cn > carrier_pl)
                )
                cn_disc = float(contra.mean())
                new.info["CN_DISC"] = round(cn_disc, 4)
                if cn_disc >= cn_discordant_frac:
                    tags.append("CNDiscordant")
            if len(carrier_cn) and site_mode(rec.contig, start, end) == "y":
                support = (
                    (carrier_cn > carrier_pl) if svtype == "DUP" else (carrier_cn < carrier_pl)
                )
                # below one depth bin the carrier CN comes from the junctions
                if support.mean() < 0.5 or end - start < CHRY_MIN_SIZE:
                    tags.append("ChrYNoDepth")

            if known is not None:
                hit = known.lookup(rec.contig, start, end, svtype)
                if hit is not None:
                    new.info["REF_AF"] = round(hit["af"], 6)
                    new.info["REF_CARRIER_FRAC"] = round(hit["carrier_frac"], 4)
                    new.info["REF_FILTER"] = hit["filter"]
                    if hit["carrier_frac"] >= known_common_frac:
                        tags.append("CommonRef")
                    if (known_artifacts and af > artifact_af and "CohortArtifact" not in tags
                            and "CohortArtifact" in hit["filter"].split("|")):
                        tags.append("CohortArtifact")
                        new.info["REF_ARTIFACT"] = True
                        counts["ref_artifact"] += 1

            if len(tc) and not sex_contig:
                gc, gf, gm = col_nalt[tc], col_nalt[tf], col_nalt[tm]
                ok = ((gc >= 0) & (gf >= 0) & (gm >= 0)
                      & (np.minimum(np.minimum(col_gq[tc], col_gq[tf]), col_gq[tm])
                         >= trio["trio_min_gq"]))
                opp = ok & (((gf == 1) & (gm == 0)) | ((gf == 0) & (gm == 1)))
                n_opp, n_tr = int(opp.sum()), int((opp & (gc > 0)).sum())
                car = ok & (gc > 0)
                dn = car & (gf == 0) & (gm == 0)
                n_car, n_dn = int(car.sum()), int(dn.sum())
                rd_c = col_rd[tc]   # NaN: no depth, no support
                depth = (rd_c <= trio["denovo_rd_del"]) if svtype == "DEL" else (
                    rd_c >= trio["denovo_rd_dup"])
                n_nodepth = int((dn & ~depth).sum())
                new.info["TR_OPP"] = n_opp
                new.info["DN_TRIO"] = n_dn
                new.info["DN_NODEPTH"] = n_nodepth
                if n_opp:
                    new.info["TR_RATE"] = round(n_tr / n_opp, 4)
                    if (n_opp >= trio["trio_min_opp"]
                            and n_tr / n_opp < trio["low_transmission"]):
                        tags.append("LowTransmission")
                if n_car:
                    new.info["DN_FRAC"] = round(n_dn / n_car, 4)
                    if (n_nodepth >= trio["denovo_excess_n"]
                            and n_nodepth / n_car >= trio["denovo_excess_frac"]):
                        tags.append("DeNovoExcess")

            for tag in tags or ["PASS"]:
                new.filter.add(tag)
            counts["artifact"] += "CohortArtifact" in tags
            counts["high_carrier"] += "HighCarrierFrac" in tags
            counts["cn_discordant"] += "CNDiscordant" in tags
            counts["common_ref"] += "CommonRef" in tags
            counts["chry_no_depth"] += "ChrYNoDepth" in tags
            counts["low_transmission"] += "LowTransmission" in tags
            counts["denovo_excess"] += "DeNovoExcess" in tags
            counts["pass"] += not tags
            # str(record) is the line VariantFile.write would emit
            fo.write(str(new).encode())
    if known is not None:
        known.vcf.close()
    return counts


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
    known_artifacts: bool = True,
    trio: dict | None = None,
    jobs: int = 1,
) -> dict[str, int]:
    """``jobs`` > 1 loads the genotype tables and writes the records in
    parallel (5 h on one core for 2670 samples): each worker writes a
    contiguous block of sites, the blocks are joined in site order, so the
    records are identical to the serial run. ``trio`` overrides
    ``TRIO_DEFAULTS`` (LowTransmission / DeNovoExcess thresholds); the trio
    flags need a ``pedigree``."""
    import shutil
    import tempfile

    import numpy as np

    from ..parallel import concat_bgzf, concat_text, pmap

    # site order defines the array index
    site_ids: list[str] = []
    x_sites: list[int] = []
    with pysam.VariantFile(str(sites_vcf)) as sv:
        for rec in sv:
            if rec.contig in ("chrX", "X"):
                x_sites.append(len(site_ids))
            site_ids.append(rec.id)
    site_idx = {sid: i for i, sid in enumerate(site_ids)}
    n_sites = len(site_ids)

    _SHARED.clear()
    _SHARED.update(site_idx=site_idx, n_sites=n_sites)
    tables: dict[str, dict] = {}
    try:
        for path, (sample, arrays) in zip(
            genotype_tsvs, pmap(_load_one, genotype_tsvs, jobs), strict=True
        ):
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
        male = hemizygous_x(tables, samples, np.array(x_sites, dtype=np.int64))
        pos = {s: k for k, s in enumerate(samples)}
        trios = [
            (pos[c.sample_id], pos[f.sample_id], pos[m.sample_id])
            for c, f, m in (pedigree.trios() if pedigree else [])
            if all(x.sample_id in pos for x in (c, f, m))
        ]
        trio_idx = tuple(np.array(col, dtype=np.int64) for col in zip(*trios, strict=True)) \
            if trios else (np.zeros(0, dtype=np.int64),) * 3
        logger.info("cohort-filter: %d complete trios (LowTransmission, DeNovoExcess)",
                    len(trios))
        logger.info("cohort-filter: %d samples with one chrX", int(male.sum()))
        _SHARED.update(
            tables=tables, samples=samples, male=male,
            is_founder=np.array([s in founders for s in samples]),
            sites_vcf=str(sites_vcf),
            known_sites=str(known_sites) if known_sites else None,
            artifact_af=artifact_af,
            artifact_gtr_discordance=artifact_gtr_discordance,
            carrier_frac_max=carrier_frac_max,
            cn_discordant_frac=cn_discordant_frac,
            known_common_frac=known_common_frac,
            known_artifacts=known_artifacts,
            trio={**TRIO_DEFAULTS, **(trio or {})}, trio_idx=trio_idx,
        )

        out_path = str(out_vcf)
        gz = out_path.endswith(".gz")
        mode = "wz" if gz else "w"
        tmp = tempfile.mkdtemp(prefix=".cohort_parts_", dir=Path(out_path).parent)
        try:
            n_chunks = 1 if jobs <= 1 else min(max(n_sites, 1), jobs * 4)
            bounds = np.linspace(0, n_sites, n_chunks + 1).astype(int)
            ext = ".vcf.gz" if gz else ".vcf"
            tasks = [
                (int(bounds[k]), int(bounds[k + 1]), f"{tmp}/part{k:05d}{ext}")
                for k in range(n_chunks)
            ]
            counts = {"sites": 0, "artifact": 0, "ref_artifact": 0, "high_carrier": 0,
                      "cn_discordant": 0, "common_ref": 0, "chry_no_depth": 0,
                      "low_transmission": 0, "denovo_excess": 0, "pass": 0}
            for c in pmap(_write_chunk, tasks, jobs):
                for k, v in c.items():
                    counts[k] += v
            head = f"{tmp}/header{ext}"
            header = _cohort_header(
                sites_vcf, samples, carrier_frac_max, cn_discordant_frac,
                known_common_frac, trio,
            )
            with pysam.VariantFile(head, mode, header=header):
                pass
            parts = [head] + [t[2] for t in tasks]
            (concat_bgzf if gz else concat_text)(parts, out_path)
        finally:
            shutil.rmtree(tmp, ignore_errors=True)
    finally:
        _SHARED.clear()
    if gz:
        pysam.tabix_index(out_path, preset="vcf", force=True)
    logger.info(
        "cohort-filter: %d sites, %d PASS, %d CohortArtifact (%d from the reference), "
        "%d HighCarrierFrac, %d CNDiscordant, %d CommonRef, %d ChrYNoDepth, "
        "%d LowTransmission, %d DeNovoExcess, %d samples",
        counts["sites"], counts["pass"], counts["artifact"], counts["ref_artifact"],
        counts["high_carrier"],
        counts["cn_discordant"], counts["common_ref"], counts["chry_no_depth"],
        counts["low_transmission"], counts["denovo_excess"],
        len(samples),
    )
    return counts
