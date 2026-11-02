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
  XXX ~0.17: VAF 1/3 and 2/3).

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
    n_sites = n_het = 0
    for c in AUTOSOMES:
        if not store.has_contig(c):
            continue
        raw, mask = _raw(store, c)
        v = raw[::AUTOSOME_STRIDE][~mask[::AUTOSOME_STRIDE]]
        auto.append(v[v > 0])
        if c in ("chr2", "chr7", "chr12"):
            s, h, _ = _het(store, c)
            n_sites += s
            n_het += h
    auto = np.concatenate(auto) if auto else np.array([])
    if len(auto) < MIN_BINS:
        return None
    base = float(np.median(auto))
    out = {"sample": store.sample, "x_cn": None, "y_cn": None, "x_het": None, "x_baf": None,
           "x_bins": 0, "y_bins": 0}
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
