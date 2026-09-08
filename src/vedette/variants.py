"""The internal CNV call model shared by all stages."""

from __future__ import annotations

from dataclasses import dataclass, field

SOURCE_DEPTH = "DEPTH"
SOURCE_JUNCTION = "JUNCTION"
SOURCE_BOTH = "BOTH"


@dataclass
class CNVCall:
    contig: str
    start: int              # 0-based inclusive
    end: int                # 0-based exclusive
    svtype: str             # "DEL" | "DUP"
    source: str             # DEPTH | JUNCTION | BOTH

    # depth evidence
    cn: int = -1            # MAP integer copy number (-1 = not evaluated)
    cn_posterior: float = float("nan")
    mean_ratio: float = float("nan")
    n_bins: int = 0
    masked_frac: float = 0.0

    # breakpoint evidence
    pe_support: int = 0
    sr_support: int = 0
    cipos: tuple[int, int] = (0, 0)
    ciend: tuple[int, int] = (0, 0)
    precise: bool = False

    # BAF evidence
    baf_n: int = 0          # informative het sites inside the call
    baf_dev: float = float("nan")  # mean |VAF - 0.5| of those sites

    # genotype
    gt: tuple[int, int] = (-1, -1)   # (-1,-1) = ./.
    gq: int = 0
    sq: float = float("nan")         # calibrated quality from the ML filter

    filters: list[str] = field(default_factory=list)
    info_extra: dict = field(default_factory=dict)

    @property
    def size(self) -> int:
        return self.end - self.start

    @property
    def svlen(self) -> int:
        return self.size if self.svtype == "DUP" else -self.size
