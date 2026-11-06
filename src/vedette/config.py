"""Configuration for Vedette.

All tunable parameters live in one frozen dataclass tree, overridable from a TOML
file (``vedette --config vedette.toml``) and/or CLI options. Defaults are chosen
for 30-40x Illumina PCR-free WGS on GRCh38.
"""

from __future__ import annotations

import dataclasses
import tomllib
from dataclasses import dataclass, field
from pathlib import Path


@dataclass(frozen=True)
class ExtractConfig:
    """Evidence extraction from the CRAM."""

    bin_size: int = 100              # self-contained depth bin size (bp)
    fast_bin_size: int = 1000        # bin size of precomputed depth bedgraphs
    min_mapq_depth: int = 1          # MAPQ threshold for the raw depth counter
    min_mapq_evidence: int = 10      # MAPQ threshold for PE/SR evidence
    min_clip_len: int = 15           # minimum soft-clip length to record
    min_clip_baseq: float = 10.0     # minimum mean base quality of the clip
    min_cigar_sv: int = 30           # CIGAR deletion (D) ops >= this are junction evidence
    cigar_merge_gap: int = 30        # D ops closer than this within one read are summed
    insert_bin_size: int = 100       # bin size of the local insert-size track
    insert_fit_pairs: int = 1_000_000  # proper pairs used to fit the insert-size model
    insert_z_discordant: float = 4.0   # |insert - mu| > z*sigma => discordant
    max_discordant_insert: int = 10_000_000  # ignore pairs beyond this span (chromosome arms)
    baf_min_dp: int = 8              # minimum depth at a SNV site to use its BAF


@dataclass(frozen=True)
class SegmentConfig:
    """CN-state HMM segmentation."""

    max_cn: int = 5                  # states are CN 0..max_cn (last = ">=max_cn")
    transition: float = 1e-8         # per-bin probability of leaving the current CN state
    min_bins: int = 3                # minimum bins for a depth-only segment to be reported
    min_cn_posterior: float = 0.90   # posterior required for a depth-only (no junction) call
    max_masked_frac: float = 0.7     # segments with more masked bins than this are dropped
    dispersion_prior: float = 0.15   # initial NB overdispersion (refit per sample)
    baf_weight: float = 1.0          # weight of the BAF emission channel (0 disables)
    # bridging of a large segment broken by a short interruption: the Gaussian
    # depth emission has light tails, so one 1 kb bin at ratio ~1.3 inside a
    # het deletion (or the 2 kb CN0 of a common deletion on the other allele)
    # costs more under CN1 than two state switches and splits the event. Two
    # same-CN segments are joined when the interruption is at most
    # bridge_max_gap_bins long AND at most bridge_max_gap_frac of the shorter
    # side, and their mean ratios differ by at most bridge_max_ratio_diff.
    # The interrupting segments are kept as their own (nested) calls.
    bridge_max_gap_bins: int = 3     # 0 disables bridging
    bridge_max_gap_frac: float = 0.05
    bridge_max_ratio_diff: float = 0.15


@dataclass(frozen=True)
class BreakpointConfig:
    """PE/SR clustering and junction resolution."""

    min_pe: int = 3                  # discordant pairs to seed a PE cluster
    min_sr: int = 2                  # clipped reads at the same position to seed an SR breakpoint
    min_cigar_support: int = 3       # reads sharing a CIGAR deletion to form a junction
    # everted-pair clusters with zero split-read/realignment support are ~99.8%
    # false (HG002: 37689 FP / 74 TP): junction-only DUPs require SR evidence
    dup_junction_requires_sr: bool = True
    # insert-size-shift detector (SurVIndel2-style) for sub-insert deletions
    insert_shift_z: float = 3.5      # per-bin z-score threshold on the local mean insert
    insert_shift_min_bins: int = 3   # consecutive shifted bins required
    insert_shift_min_pairs: int = 10 # spanning pairs required across the run
    sr_pos_tolerance: int = 5        # bp tolerance when piling up clip positions
    max_junction_size: int = 20_000_000   # ignore candidate junctions larger than this
    min_size: int = 50               # minimum SV size reported
    junction_only_max_size: int = 10_000  # junction without depth support kept below this size
    # hidden-split-read recovery (clip realignment)
    hidden_splits: bool = True
    resolve_window: int = 2_000      # realignment window around each clip pileup (bp)
    resolve_min_clip: int = 20       # clips shorter than this are unmappable
    resolve_min_identity: float = 0.85
    resolve_max_clips: int = 3       # longest clips realigned per pileup
    resolve_min_support: int = 2     # concordant realigned clips to call a junction
    # per-contig budget of fit alignments: wavy samples fragment the HMM and
    # flood the recovery with pileups (5x median runtime); healthy 40x genomes
    # stay far below this, pathological ones keep their strongest pileups
    resolve_max_alignments: int = 150_000


@dataclass(frozen=True)
class GenotypeConfig:
    """Unified genotype likelihood."""

    mismap_rate: float = 0.005       # P(supporting read | no variant)
    # junction-read counts are overdispersed (mapping efficiency varies per
    # locus): negative-binomial variance = mu + alpha*mu^2. Poisson (alpha=0)
    # over-trusts the count and miscalls het vs hom on small events — the
    # dominant Mendelian-violation mode measured on the CANDY_ghent cohort.
    junction_overdispersion: float = 0.75
    baf_overdispersion: float = 50.0  # beta-binomial concentration for BAF terms
    het_prior: float = 1e-4          # prior P(0/1) at a random locus
    hom_prior: float = 1e-6          # prior P(1/1)
    min_gq_report: int = 0           # report all, filtering happens downstream
    # genotype on the BAF rewritten by `vedette baf` (no QUAL threshold)
    # instead of the store's original one (GHFC track, QUAL >= 30). Off: the
    # likelihood counts every het inside a deletion as a mismap, and in
    # paralogous loci the rewritten BAF has ~20 per kb, so a common deletion
    # genotyped 0/1 in 98 % of dataset_paper fell to 0/0 in half of them,
    # passed the cohort artifact filters, and made de novo calls (chr12:
    # 3,025 -> 3,653 with the rewritten BAF, 2,840 with the original one).
    rewritten_baf: bool = False


@dataclass(frozen=True)
class FilterConfig:
    """Tier-1 hard filters (Tier-2 ML model has its own versioned manifest)."""

    min_size: int = 50
    max_excluded_overlap: float = 0.30   # fraction of call overlapping N-gaps/exclude mask
    lowmap_segdup_frac: float = 0.70     # segdup overlap above which SR support is required
    min_informative_bins: int = 3
    low_sq_threshold: float = 5.0        # phred (P(true) ~0.68): operating point of models/v5n
    # Large deletions sit outside the model's training range: only 131 of
    # 132555 labelled DEL examples exceed 1 Mb (0.1%), and those all have
    # degenerate features (dhffc=0, fc_chrom=0, no CN) because they are
    # total-dropout regions. A real multi-megabase heterozygous deletion looks
    # nothing like them, and worse, once the HMM splits it into adjacent
    # segments each segment's flanks fall inside the neighbouring segment, so
    # dhffc ~ 1 and the model sees no event. On PMS that rejected half of all
    # 1-10 Mb deletions, including the causal terminal 22q13 ones. fc_chrom
    # (ratio to the chromosome median) stays informative there, so a call this
    # large whose chromosome-level depth agrees with a copy loss keeps its
    # PASS regardless of SQ. Retains 131/131 true and 0 false on every
    # labelled call we have.
    # 200 kb, not 1 Mb: on the labelled calls this rescues 5 rejected
    # deletions, all 5 true; at 100 kb it starts taking false ones (13 true /
    # 10 false), and at 1 Mb it misses the 229 kb first segment of the PMS
    # terminal deletions.
    large_del_rescue_size: int = 200_000
    large_del_rescue_fc_chrom: float = 0.75
    large_del_rescue_max_mask: float = 0.5
    # Large duplications fare worse: GIAB's depth-only DUPs are artefacts, so
    # the model scores every depth-only gain near 0 — a trisomy 21, a 15q11-q13
    # duplication, a 7q11.23 duplication and a trisomy 13 were all LowSQ with
    # SQ 0.01-0.09 (2026-10-15). Coverage waves mimic gains, so depth cannot
    # vouch for them; the allelic balance can: at CN 3 the heterozygous sites
    # sit at 1/3 and 2/3. DUP pieces inside the scan-confirmed constitutional
    # gains had mean |VAF - 0.5| 0.116-0.153 (91 % >= 0.12) against 0.053-0.09
    # elsewhere; segmental duplications fake that imbalance (paralogous
    # variants), so they are excluded. Rescues one DUP per genome in median
    # (p90 2), the extremes being the real gains.
    large_dup_rescue_size: int = 100_000
    large_dup_rescue_min_baf_n: int = 30
    large_dup_rescue_min_baf_dev: float = 0.12
    large_dup_rescue_max_segdup: float = 0.5
    large_dup_rescue_max_mask: float = 0.5


@dataclass(frozen=True)
class ReferenceConfig:
    """Reference genome and companion tracks."""

    fasta: str = ""
    exclude_bed: str = ""            # e.g. exclude.cnvnator_100bp.GRCh38.20170403.bed
    segdup_bed: str = ""
    lcr_bed: str = ""
    pop_cov_stdev_bw: str = ""       # population coverage stdev bigWig (optional)
    # contigs used for calling; everything else only feeds normalization
    call_contigs: tuple[str, ...] = tuple(
        [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    )


@dataclass(frozen=True)
class Config:
    reference: ReferenceConfig = field(default_factory=ReferenceConfig)
    extract: ExtractConfig = field(default_factory=ExtractConfig)
    segment: SegmentConfig = field(default_factory=SegmentConfig)
    breakpoint: BreakpointConfig = field(default_factory=BreakpointConfig)
    genotype: GenotypeConfig = field(default_factory=GenotypeConfig)
    filter: FilterConfig = field(default_factory=FilterConfig)


_SECTIONS = {
    "reference": ReferenceConfig,
    "extract": ExtractConfig,
    "segment": SegmentConfig,
    "breakpoint": BreakpointConfig,
    "genotype": GenotypeConfig,
    "filter": FilterConfig,
}


def load_config(path: str | Path | None = None, **overrides: dict) -> Config:
    """Build a Config from defaults, an optional TOML file, then keyword overrides.

    ``overrides`` maps section name -> {field: value}.
    """
    data: dict = {}
    if path is not None:
        with open(path, "rb") as fh:
            data = tomllib.load(fh)
    for section, values in overrides.items():
        if values:
            data.setdefault(section, {}).update(values)

    kwargs = {}
    for section, cls in _SECTIONS.items():
        section_data = data.get(section, {})
        valid = {f.name for f in dataclasses.fields(cls)}
        unknown = set(section_data) - valid
        if unknown:
            raise ValueError(f"Unknown option(s) in [{section}]: {sorted(unknown)}")
        if "call_contigs" in section_data:
            section_data["call_contigs"] = tuple(section_data["call_contigs"])
        kwargs[section] = cls(**section_data)
    return Config(**kwargs)
