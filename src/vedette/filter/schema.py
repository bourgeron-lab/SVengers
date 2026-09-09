"""Feature schema for the calibrated filter model.

The ORDER of FEATURES is part of the model contract: a trained model ships a
manifest carrying ``schema_hash``, and scoring refuses a mismatch. Add new
features only by appending and bumping the model version.
"""

from __future__ import annotations

import hashlib

FEATURES: list[str] = [
    # call-level
    "size_log10",
    "is_dup",
    "precise",
    "src_both",
    "src_junction",
    "maskfrac",
    # depth
    "dhffc",             # call median depth / flank median depth
    "fc_chrom",          # call median depth / chromosome median depth
    "fc_residual",       # |mean copy ratio - expected under called CN|
    "depth_cv",          # std/mean of in-call bins
    "frac_consistent",   # fraction of bins within 0.25 of expected ratio
    "n_bins_informative",
    "flank_asym",        # |left flank - right flank| / mean flank depth
    # breakpoint
    "pe",
    "sr",
    "cig",               # reads carrying the event as a CIGAR deletion
    "ins_z",             # insert-size-shift z-score (0 unless SOURCE=INSERT)
    "ins_pairs",         # spanning pairs behind the shift estimate
    "evidence_ratio",    # (pe + sr + cig) / sequence coverage
    "homlen",            # junction micro-homology (-1 = not measured)
    # BAF
    "baf_n",
    "baf_dev",           # mean |VAF - 0.5| of in-call clear hets (-1 if none)
    # context
    "segdup_ovl",
    "lcr_ovl",
    "exclude_ovl",
    "popcov_stdev",
    # caller internals
    "gq",
    "cn",
    "cnq",
]


def schema_hash() -> str:
    return hashlib.sha256("|".join(FEATURES).encode()).hexdigest()[:16]
