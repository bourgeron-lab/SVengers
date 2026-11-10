"""Sex-chromosome karyotype of one genome, from its evidence store.

The mosaic scan covers the autosomes only, and the caller infers a sex to set
the X/Y baseline (male X = 1 copy) — after which an extra X is an ordinary
gain, a missing one invisible. A dataset_paper genome with chrX at CN ~2.8
(XXX) surfaced only as fragmented PASS duplications. This pass reads, per
genome:

- ``x_cn`` / ``y_cn``: twice the median raw depth over the informative bins of
  chrX / chrY, relative to the autosomal median. The bins come from the
  panel of normals' sex-stratified profiles: X bins at ~1.0 in women and ~0.5
  in men (the pseudo-autosomal regions, at 1.0 in both, drop out), Y bins at
  ~0.5 in men and ~0 in women (the X-transposed and ampliconic repeats drop
  out). Raw depth, not the corrected track: the corrected one is divided by
  the panel of the INFERRED sex, which is exactly what is tested here.
- ``x_het``: heterozygous-site rate on the chrX bins over the autosomal rate
  (one X: ~0, sequencing errors only; two or more: ~0.7, the X's lower
  diversity);
- ``x_baf``: median |BAF - 0.5| of those heterozygous sites (XX ~0.07 at 30x,
  XXX ~0.17: VAF 1/3 and 2/3);
- ``gc_contrast``: median copy number of the GC-rich autosomes (16, 17, 19,
  22) minus that of the AT-rich ones (4, 5, 13, 18). Raw depth carries the
  genome's GC bias, and chrX, AT-rich, reads like chr4: a dataset_paper man
  (wave 0.50, contrast 1.23) had chrX at 0.84 and chrY at 0.65, a PMS woman
  (contrast 1.06) chrX at 1.70 — two "mosaics" that were GC alone. Medians of
  four, so one trisomic chromosome (13, 18) does not move it. Of 2989
  genomes, 14 above 0.5; of the 27 the karyotype pass flagged, these two
  above 0.3 and the rest below 0.28.

``classify`` turns them into XX, XY, X0, XXY, XXX, XYY, or an
"intermediate" level (a mosaic such as 45,X/46,XX, or a partial loss of Y),
and flags a pedigree sex that disagrees.
"""

from __future__ import annotations

import numpy as np

AUTOSOMES = tuple(f"chr{i}" for i in range(1, 23))
AUTOSOME_STRIDE = 7          # every 7th bin: the autosomal median, ~0.4 M bins
MIN_BINS = 500
# integer copy-number tolerance; between two integers the level is reported
# as intermediate (mosaic or partial loss)
CN_TOL = 0.25
X_HET_ONE = 0.15             # x_het below this: a single X
X_BAF_THREE = 0.12           # x_baf above this with x_cn ~3: three X
GC_RICH = ("chr16", "chr17", "chr19", "chr22")
AT_RICH = ("chr4", "chr5", "chr13", "chr18")
GC_CONTRAST_MAX = 0.5        # above this an abnormal X/Y level is not read


def informative_bins(pon_root, contig: str) -> np.ndarray | None:
    """Bin indices of chrX / chrY whose panel profile tells the sexes apart."""
    try:
        g = pon_root[contig]
        med_m = np.asarray(g["median_male"][:], dtype=np.float64)
        med_f = np.asarray(g["median_female"][:], dtype=np.float64)
        mad_m = np.asarray(g["mad_male"][:], dtype=np.float64)
    except KeyError:
        return None
    if contig == "chrX":
        ok = (np.abs(med_f - 1.0) < 0.15) & (np.abs(med_m - 0.5) < 0.1) & (mad_m < 0.15)
    else:
        ok = (np.abs(med_m - 0.5) < 0.15) & (med_f < 0.05) & (mad_m < 0.2)
    return np.nonzero(ok)[0]


def _raw(store, contig: str) -> tuple[np.ndarray, np.ndarray]:
    raw, _corr, mask = store.depth(contig)
    return np.asarray(raw, dtype=np.float64), np.asarray(mask, dtype=bool)


def _het(store, contig: str, bins: np.ndarray | None = None) -> tuple[int, int, np.ndarray]:
    """(sites, heterozygous sites, |BAF - 0.5| of those) on a contig, or on
    the given bins of it."""
    b = store.baf(contig, 0, 1 << 30)
    pos = np.asarray(b["pos"], dtype=np.int64)
    vaf = np.asarray(b["vaf"], dtype=np.float64)
    if bins is not None:
        keep = np.isin(pos // store.bin_size, bins)
        pos, vaf = pos[keep], vaf[keep]
    het = (vaf > 0.1) & (vaf < 0.9)
    return len(vaf), int(het.sum()), np.abs(vaf[het] - 0.5)


def sex_chromosome_stats(store, pon_root) -> dict | None:
    """x_cn, y_cn, x_het, x_baf of one store (None without chrX or bins)."""
    if not store.has_contig("chrX"):
        return None
    auto = []
    per_chrom = {}
    n_sites = n_het = 0
    for c in AUTOSOMES:
        if not store.has_contig(c):
            continue
        raw, mask = _raw(store, c)
        v = raw[::AUTOSOME_STRIDE][~mask[::AUTOSOME_STRIDE]]
        auto.append(v[v > 0])
        if len(auto[-1]):
            per_chrom[c] = float(np.median(auto[-1]))
        if c in ("chr2", "chr7", "chr12"):
            s, h, _ = _het(store, c)
            n_sites += s
            n_het += h
    auto = np.concatenate(auto) if auto else np.array([])
    if len(auto) < MIN_BINS:
        return None
    base = float(np.median(auto))
    gc = [2.0 * per_chrom[c] / base for c in GC_RICH if c in per_chrom]
    at = [2.0 * per_chrom[c] / base for c in AT_RICH if c in per_chrom]
    out = {"sample": store.sample, "x_cn": None, "y_cn": None, "x_het": None, "x_baf": None,
           "x_bins": 0, "y_bins": 0,
           "gc_contrast": round(float(np.median(gc) - np.median(at)), 3) if gc and at else None}
    for c in ("chrX", "chrY"):
        bins = informative_bins(pon_root, c)
        if bins is None or not store.has_contig(c):
            continue
        raw, _mask = _raw(store, c)
        bins = bins[bins < len(raw)]
        if len(bins) < MIN_BINS:
            continue
        key = "x" if c == "chrX" else "y"
        out[f"{key}_cn"] = round(2.0 * float(np.median(raw[bins])) / base, 3)
        out[f"{key}_bins"] = int(len(bins))
        if c == "chrX" and n_sites:
            s, h, dev = _het(store, c, bins)
            if s:
                out["x_het"] = round((h / s) / (n_het / n_sites), 3)
                out["x_baf"] = round(float(np.median(dev)), 4) if len(dev) >= 50 else None
    return out


def _near(v: float | None, k: int) -> bool:
    return v is not None and abs(v - k) <= CN_TOL


def classify(st: dict, ped_sex: str | None = None) -> dict:
    """{'karyotype', 'note'} from sex_chromosome_stats; ped_sex 'M'/'F'/None."""
    x, y, het, baf = st.get("x_cn"), st.get("y_cn"), st.get("x_het"), st.get("x_baf")
    if x is None:
        return {"karyotype": "?", "note": "chrX non lisible"}
    nx = next((k for k in (1, 2, 3, 4) if _near(x, k)), None)
    ny = 0 if (y is None or y < CN_TOL) else next((k for k in (1, 2) if _near(y, k)), None)
    one_x = het is not None and het < X_HET_ONE
    notes = []
    if nx is None:
        kar = "intermédiaire"
        notes.append(f"chrX à {x:.2f} copies (mosaïque ?)")
    elif ny is None:
        kar = "intermédiaire"
        notes.append(f"chrY à {y:.2f} copie (perte partielle du Y ?)")
    else:
        kar = "X" * nx + "Y" * ny if (nx, ny) != (1, 0) else "X0"
        if nx >= 2 and one_x:
            notes.append(f"{nx} copies de chrX sans hétérozygotie (isodisomie ?)")
        if nx == 1 and het is not None and not one_x:
            # a single X cannot be heterozygous: on dataset_paper, 4 men at
            # x_het 0.22-0.57 with BAF ~0.25-0.30, three of them with
            # consecutive barcodes — DNA from another person in the sample
            notes.append(f"un seul chrX mais hétérozygote ({het:.2f}) : contamination de l'ADN ?")
        if nx == 3 and baf is not None and baf < X_BAF_THREE:
            notes.append(f"BAF de chrX {baf:.3f}, attendue ~0,17 pour XXX")
    gcc = st.get("gc_contrast")
    if gcc is not None and abs(gcc) > GC_CONTRAST_MAX and kar not in ("XX", "XY"):
        # an XX/XY reading survives the bias; anything else may be the bias
        return {"karyotype": "?", "note": f"biais GC du génome (contraste {gcc:.2f} copie entre "
                                          "autosomes riches et pauvres en GC) : niveaux X/Y non "
                                          "interprétables"}
    if ped_sex in ("M", "F") and kar not in ("?", "intermédiaire"):
        has_y = "Y" in kar
        if (ped_sex == "M") != has_y:
            notes.append(f"sexe du pedigree {'masculin' if ped_sex == 'M' else 'féminin'}")
    return {"karyotype": kar, "note": " ; ".join(notes)}


def is_anomaly(kar: dict) -> bool:
    return kar["karyotype"] not in ("XX", "XY") or bool(kar["note"])


def aneuploid_sex_chromosomes(karyotype: str) -> set[str]:
    """Sex chromosomes ("X", "Y") whose copy number the karyotype changes
    against the sex it implies: X for XXY, XXX or X0, Y for XYY. The whole
    chromosome then reads as a de novo gain or loss against the parents."""
    if karyotype in ("XX", "XY", "?", "intermédiaire"):
        return set()
    nx, ny = karyotype.count("X"), karyotype.count("Y")
    out = set()
    if nx != (1 if ny else 2):
        out.add("X")
    if ny > 1:
        out.add("Y")
    return out


def sex_ploidy(karyotype: str) -> dict[str, int] | None:
    """Baseline copy number of chrX / chrY for a karyotype read by ``classify``
    (XXY -> X 2, Y 1); None when the reading is not an integer karyotype
    ("?", "intermédiaire"): the caller then keeps its own sex inference."""
    if not karyotype or set(karyotype) - {"X", "Y", "0"} or "X" not in karyotype:
        return None
    return {"chrX": karyotype.count("X"), "chrY": karyotype.count("Y")}


class ArraySource:
    """The store interface ``sex_chromosome_stats`` reads (``sample``,
    ``bin_size``, ``has_contig``, ``depth``, ``baf``), over the in-memory
    arrays of ``vedette call`` — the karyotype then sets the X/Y baseline
    before segmentation instead of being read from the finished store."""

    def __init__(self, sample: str, bin_size: int, raw: dict, mask: dict, baf: dict):
        self.sample, self.bin_size = sample, bin_size
        self._raw, self._mask, self._baf = raw, mask, baf

    def has_contig(self, contig: str) -> bool:
        return contig in self._raw

    def depth(self, contig: str):
        return self._raw[contig], self._raw[contig], self._mask[contig]

    def baf(self, contig: str, start: int, end: int) -> dict:
        pos, vaf = (self._baf.get(contig) or (np.array([], np.uint32), np.array([], np.float32)))[:2]
        keep = (pos >= start) & (pos < end)
        return {"pos": pos[keep], "vaf": vaf[keep]}


# ---------------------------------------------------------------- autosomes
# Whole-chromosome aneuploidy of an autosome (trisomy 21, a trisomic chr12):
# with a diploid baseline the extra copy is called as fragmented duplications
# (dataset_paper: 15 on a trisomy 21, 117 on a chr12 at 2.92) and a deletion
# of one of the three copies sits at CN 2 and is missed (that chr12: 93 PASS
# DEL against a cohort median of 149). The call then takes the chromosome's
# own copy number as its baseline; the aneuploidy itself stays reported by
# the mosaic scan (CHROMOSOMIQUE rows), which reads the depth, not the ploidy.
ANEU_WINDOW_BP = 1_000_000
ANEU_MIN_WINDOWS = 20
ANEU_WINDOW_TOL = 0.3        # a window "at" k copies: |2 x median ratio - k| <= this
ANEU_MIN_FRAC = 0.85         # of the informative windows, at the same level
ANEU_GAP_CN = 0.25           # a window below this is an uncovered gap, not a copy number
ANEU_BAF_THREE = 0.12        # median |BAF - 0.5| above this: VAF 1/3, 2/3 (diploid ~0.06)
ANEU_HET_ONE = 0.15          # het rate / genome rate below this: a single copy


def autosomal_aneuploidies(
    ratio: dict[str, np.ndarray], mask: dict[str, np.ndarray], baf: dict, bin_size: int,
    autosomes: list[str],
) -> dict[str, int]:
    """{autosome: copy number} for whole autosomes at 1 or 3 copies, confirmed
    by the BAF; ``ratio`` is the GC-corrected depth over the autosomal mean."""
    w = max(1, ANEU_WINDOW_BP // bin_size)

    def hets(c):
        pos, vaf = baf.get(c, (np.array([]), np.array([])))[:2]
        vaf = np.asarray(vaf, dtype=np.float64)
        h = vaf[(vaf > 0.1) & (vaf < 0.9)]
        return len(vaf), h

    n_all = n_het = 0
    for c in autosomes:
        s, h = hets(c)
        n_all, n_het = n_all + s, n_het + len(h)
    out = {}
    for c in autosomes:
        r, m = np.asarray(ratio[c], dtype=np.float64), np.asarray(mask[c], dtype=bool)
        cns = []
        for i in range(0, len(r) - w + 1, w):
            v = r[i:i + w][~m[i:i + w]]
            if len(v) >= w // 2:
                cn = 2.0 * float(np.median(v))
                # the reference mask leaves some pericentromeric / acrocentric
                # gaps unmasked (no reads at all): no copy-number information.
                # On chr21 (~37 windows) four such windows flanking the
                # centromere put a trisomy just below the agreement threshold
                if cn >= ANEU_GAP_CN:
                    cns.append(cn)
        if len(cns) < ANEU_MIN_WINDOWS:
            continue
        cns = np.array(cns)
        k = int(round(float(np.median(cns))))
        if k not in (1, 3) or np.mean(np.abs(cns - k) <= ANEU_WINDOW_TOL) < ANEU_MIN_FRAC:
            continue
        s, h = hets(c)
        if not s or not n_all or not n_het:
            continue
        if k == 3 and len(h) >= 200 and float(np.median(np.abs(h - 0.5))) > ANEU_BAF_THREE:
            out[c] = 3
        elif k == 1 and (len(h) / s) / (n_het / n_all) < ANEU_HET_ONE:
            out[c] = 1
    return out
