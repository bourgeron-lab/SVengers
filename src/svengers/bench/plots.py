"""Publication figures from the bench tables (matplotlib, no styling deps)."""

from __future__ import annotations

import logging
from pathlib import Path

logger = logging.getLogger(__name__)

_SIZE_ORDER = ["50-300bp", "300bp-1kb", "1-10kb", "10-100kb", "100kb-1Mb", ">1Mb"]


def plot_size_stratified(
    rows: list[dict], out_png: str | Path, metric: str = "f1", svtype: str | None = None
) -> None:
    """Grouped bars of a metric per size class, one group per caller.

    ``rows`` is the tidy output of bench.compare (must contain a 'caller' key)
    or bench.truth (a single implicit caller).
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    svtypes = sorted({r["svtype"] for r in rows}) if svtype is None else [svtype]
    fig, axes = plt.subplots(
        1, len(svtypes), figsize=(6.5 * len(svtypes), 4.2), squeeze=False
    )
    for ax, st in zip(axes[0], svtypes, strict=True):
        sub = [
            r for r in rows
            if r["svtype"] == st and r["stratum_kind"] == "size" and r[metric] != ""
        ]
        callers = sorted({r.get("caller", "svengers") for r in sub})
        width = 0.8 / max(len(callers), 1)
        for k, caller in enumerate(callers):
            vals, xs = [], []
            for i, size_class in enumerate(_SIZE_ORDER):
                match = [
                    r for r in sub
                    if r.get("caller", "svengers") == caller and r["stratum"] == size_class
                ]
                if match:
                    xs.append(i + k * width)
                    vals.append(float(match[0][metric]))
            ax.bar(xs, vals, width=width, label=caller)
        ax.set_xticks(range(len(_SIZE_ORDER)))
        ax.set_xticklabels(_SIZE_ORDER, rotation=30, ha="right")
        ax.set_ylim(0, 1)
        ax.set_ylabel(metric)
        ax.set_title(f"{st} — {metric} by size class")
        ax.legend(fontsize=8)
        ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    logger.info("figure -> %s", out_png)


def plot_trio_sweep(rows: list[dict], out_png: str | Path) -> None:
    """MVR and transmission rate vs quality threshold (twin axes)."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    th = [r["threshold"] for r in rows]
    mvr = [float(r["mvr"]) if r["mvr"] != "" else None for r in rows]
    trans = [float(r["transmission"]) if r["transmission"] != "" else None for r in rows]
    field = rows[0]["field"] if rows else "GQ"

    fig, ax1 = plt.subplots(figsize=(6, 4))
    ax1.plot(th, mvr, "o-", color="tab:red", label="Mendelian violation rate")
    ax1.set_xlabel(f"min {field}")
    ax1.set_ylabel("MVR", color="tab:red")
    ax1.tick_params(axis="y", labelcolor="tab:red")
    ax2 = ax1.twinx()
    ax2.plot(th, trans, "s-", color="tab:blue", label="Transmission rate")
    ax2.axhline(0.5, color="tab:blue", linestyle=":", alpha=0.5)
    ax2.set_ylabel("Transmission rate", color="tab:blue")
    ax2.tick_params(axis="y", labelcolor="tab:blue")
    ax1.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    logger.info("figure -> %s", out_png)
