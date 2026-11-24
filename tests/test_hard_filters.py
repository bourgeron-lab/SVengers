from svengers.config import load_config
from svengers.filter.hard_filters import hard_filter_tags


def _feats(**kw):
    base = {"exclude_ovl": 0.0, "maskfrac": 0.0, "segdup_ovl": 0.0, "pe": 0.0, "sr": 0.0,
            "cig": 0.0, "ins_pairs": 0.0, "n_bins_informative": 5.0, "is_dup": 0.0}
    return {**base, **kw}


def test_pon_masked_deletion_with_junction_evidence_is_not_maxn():
    """Bins masked by the PoN (population-variable) must not hard-filter a junction-supported
    deletion below 10 kb; depth-only calls, larger calls, duplications and N/exclude overlaps
    still get MaxN."""
    f = load_config(None).filter
    masked = _feats(maskfrac=1.0, n_bins_informative=0.0)
    assert "MaxN" not in hard_filter_tags(dict(masked, pe=4, sr=6), 3_000, f)
    assert "MaxN" not in hard_filter_tags(dict(masked, cig=2), 9_999, f)
    assert "MaxN" in hard_filter_tags(masked, 3_000, f)                        # no junction evidence
    assert "MaxN" in hard_filter_tags(dict(masked, pe=5, sr=2), 2_150_000, f)  # pericentromeric artefact
    assert "MaxN" in hard_filter_tags(dict(masked, pe=4, sr=6), 10_000, f)
    assert "MaxN" in hard_filter_tags(_feats(exclude_ovl=0.5, pe=4, sr=6), 3_000, f)
    assert "MaxN" in hard_filter_tags(dict(masked, pe=4, sr=6, is_dup=1.0), 600, f)  # DUP: no exemption
    assert "MaxN" not in hard_filter_tags(_feats(maskfrac=0.7), 3_000, f)


def test_large_informative_deletion_is_not_lowmap():
    """A deletion between segmental-duplication blocks (1q21.1) with >= 500
    unmasked depth bins keeps clear of LowMap; a DUP or a short run does not."""
    f = load_config(None).filter
    sd = dict(segdup_ovl=0.8, maskfrac=0.67)
    assert "LowMap" not in hard_filter_tags(_feats(n_bins_informative=1950.0, **sd), 5_900_000, f)
    assert "LowMap" in hard_filter_tags(_feats(n_bins_informative=200.0, **sd), 600_000, f)
    assert "LowMap" in hard_filter_tags(_feats(n_bins_informative=1950.0, is_dup=1.0, **sd),
                                        5_900_000, f)
