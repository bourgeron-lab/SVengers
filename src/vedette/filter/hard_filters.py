"""Tier-1 hard filters: transparent, rule-based, never silently dropping calls."""

from __future__ import annotations

from ..config import FilterConfig


def hard_filter_tags(features: dict[str, float], size: int, config: FilterConfig) -> list[str]:
    tags: list[str] = []
    if size < config.min_size:
        tags.append("MinSize")
    if features["exclude_ovl"] > config.max_excluded_overlap or features["maskfrac"] > 0.7:
        tags.append("MaxN")
    if features["segdup_ovl"] > config.lowmap_segdup_frac and features["sr"] == 0:
        tags.append("LowMap")
    if (
        features["pe"] == 0
        and features["sr"] == 0
        and features["n_bins_informative"] < config.min_informative_bins
    ):
        tags.append("NoEvidence")
    return tags
