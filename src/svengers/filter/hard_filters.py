"""Tier-1 hard filters: transparent, rule-based, never silently dropping calls."""

from __future__ import annotations

from ..config import FilterConfig


def hard_filter_tags(features: dict[str, float], size: int, config: FilterConfig) -> list[str]:
    tags: list[str] = []
    if size < config.min_size:
        tags.append("MinSize")
    # PoN-masked bins: a junction-supported deletion < 10 kb is exempt (FilterConfig)
    junction = features["pe"] + features["sr"] + features["cig"] > 0
    masked = features["maskfrac"] > config.max_maskfrac and not (
        junction and features["is_dup"] == 0 and size < config.maxn_junction_exempt_max_size
    )
    if features["exclude_ovl"] > config.max_excluded_overlap or masked:
        tags.append("MaxN")
    if (
        features["segdup_ovl"] > config.lowmap_segdup_frac
        and features["sr"] == 0
        and features["cig"] == 0
        and not (features["is_dup"] == 0
                 and features["n_bins_informative"] >= config.informative_core_bins)
    ):
        tags.append("LowMap")
    if (
        features["pe"] == 0
        and features["sr"] == 0
        and features["cig"] == 0
        and features["ins_pairs"] == 0
        and features["n_bins_informative"] < config.min_informative_bins
    ):
        tags.append("NoEvidence")
    return tags
