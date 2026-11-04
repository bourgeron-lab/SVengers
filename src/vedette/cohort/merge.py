"""`vedette merge`: cross-sample clustering of filtered calls into cohort sites.

Anti-AF-inflation design (the SURVIVOR lesson): no unbounded transitive
chaining. A call joins a cluster only if it is compatible with the cluster's
running CENTROID (median start/end of members), so cluster diameter stays
bounded and distinct nearby events are not fused. DEL and DUP cluster
independently; precise (split-read-resolved) and imprecise calls use
different tolerances.

Output: a sites-only VCF (no genotypes) consumed by `vedette genotype`.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from functools import partial
from pathlib import Path

import numpy as np
import pysam

from .. import __version__
from ..parallel import pmap
from ..reference import Contig
from ..vcfout import add_coords_line, clamp_end, span, spec_coords, sv_record

logger = logging.getLogger(__name__)


@dataclass
class MemberCall:
    sample: str
    start: int
    end: int
    precise: bool
    gq: int
    passing: bool = True   # PASS (or kept-filter) in the single-sample VCF
    sq: float = 0.0


@dataclass
class SiteCluster:
    svtype: str
    members: list[MemberCall] = field(default_factory=list)
    # cached centroid: recomputing the median on EVERY start/end access made
    # cluster_calls quadratic-ish (py-spy'd at 6h+ on the 2670-sample cohort —
    # each sweep step probes every active cluster's centroid)
    _centroid: tuple[int, int] | None = None

    def add(self, member: MemberCall) -> None:
        self.members.append(member)
        self._centroid = None

    def _compute(self) -> tuple[int, int]:
        if self._centroid is None:
            self._centroid = (
                int(np.median([m.start for m in self.members])),
                int(np.median([m.end for m in self.members])),
            )
        return self._centroid

    @property
    def start(self) -> int:
        return self._compute()[0]

    @property
    def end(self) -> int:
        return self._compute()[1]

    @property
    def spread(self) -> int:
        starts = [m.start for m in self.members]
        ends = [m.end for m in self.members]
        return int((max(starts) - min(starts)) + (max(ends) - min(ends)))


def _compatible(call: MemberCall, cluster: SiteCluster) -> bool:
    c_start, c_end = cluster.start, cluster.end
    size = max(call.end - call.start, c_end - c_start, 1)
    inter = min(call.end, c_end) - max(call.start, c_start)
    if inter <= 0:
        return False
    ro = inter / size  # reciprocal overlap vs the larger event
    if call.precise:
        bp_tol = max(150, min(int(0.1 * size), 1000))
        return (
            ro >= 0.8
            and abs(call.start - c_start) <= bp_tol
            and abs(call.end - c_end) <= bp_tol
        )
    bp_tol = max(2000, int(0.2 * size))
    return (
        ro >= 0.7
        and abs(call.start - c_start) <= bp_tol
        and abs(call.end - c_end) <= bp_tol
    )


def cluster_calls(
    calls: list[tuple[MemberCall, str]],  # (call, svtype), sorted by start
) -> list[SiteCluster]:
    """Greedy sweep clustering with centroid-compatibility (per contig)."""
    open_clusters: dict[str, list[SiteCluster]] = {"DEL": [], "DUP": []}
    done: list[SiteCluster] = []
    horizon = 1_000_000  # clusters ending this far before the sweep line are closed

    for call, svtype in calls:
        active = open_clusters[svtype]
        still_open = []
        for cl in active:
            if cl.end + horizon < call.start:
                done.append(cl)
            else:
                still_open.append(cl)
        open_clusters[svtype] = active = still_open

        best, best_ro = None, 0.0
        for cl in active:
            if _compatible(call, cl):
                inter = min(call.end, cl.end) - max(call.start, cl.start)
                ro = inter / max(call.end - call.start, cl.end - cl.start)
                if ro > best_ro:
                    best, best_ro = cl, ro
        if best is not None:
            best.add(call)
        else:
            active.append(SiteCluster(svtype=svtype, members=[call]))

    for svtype in ("DEL", "DUP"):
        done.extend(open_clusters[svtype])
    return done


def read_sample_calls(
    vcf_path: str | Path,
    min_gq: int = 0,
    pass_only: bool = True,
    keep_filters: tuple[str, ...] = ("LowPosterior",),
    rescue_min_sq: float | None = 1.0,
) -> dict[tuple[str, str], list[MemberCall]]:
    """{(contig, svtype): calls} from one filtered single-sample VCF.

    A LowSQ-only call with SQ >= rescue_min_sq is loaded as a NON-passing
    member: it cannot nominate a site alone, but recurrent sub-threshold
    calls across samples can (see run_merge) — a familial CNV in a hard
    region scores SQ ~2 in every carrier and would otherwise vanish from
    the cohort (real case: 30 kb DEL over the SHANK3 start, 5 carriers).
    """
    out: dict[tuple[str, str], list[MemberCall]] = {}
    with pysam.VariantFile(str(vcf_path)) as vcf:
        sample = list(vcf.header.samples)[0]
        spec = spec_coords(vcf.header)
        for rec in vcf:
            svtype = rec.info.get("SVTYPE")
            if svtype not in ("DEL", "DUP"):
                continue
            fmt = rec.samples[sample]
            sq = fmt.get("SQ")
            sq = float(sq) if sq is not None else 0.0
            extra = set(rec.filter.keys()) - {"PASS"} - set(keep_filters)
            passing = not extra
            if pass_only and not passing:
                # only the calibrated-score filter is rescuable; hard-filter
                # failures (NoEvidence, MaxN, ...) stay out
                if (
                    rescue_min_sq is None
                    or extra != {"LowSQ"}
                    or sq < rescue_min_sq
                ):
                    continue
            fmt = rec.samples[sample]
            gt = fmt.get("GT")
            if gt is None or not any(a == 1 for a in gt if a is not None):
                continue  # non-carrier records don't nominate sites
            gq = int(fmt.get("GQ") or 0)
            if gq < min_gq:
                continue
            start, end = span(rec, spec)
            out.setdefault((rec.contig, svtype), []).append(
                MemberCall(
                    sample=sample,
                    start=start,
                    end=end,
                    precise="IMPRECISE" not in rec.info,
                    gq=gq,
                    passing=passing or not pass_only,
                    sq=sq,
                )
            )
    return out


_BY_KEY: dict[tuple[str, str], list[MemberCall]] | None = None


def _contig_sites(work: tuple) -> list[tuple[int, int, str, list[tuple]]]:
    """Cluster one contig's calls: (start, end, svtype, INFO items) per site,
    in output order."""
    contig, length, rescue_min_samples, rescue_total_sq = work
    by_key = _BY_KEY
    merged: list[tuple[int, SiteCluster]] = []
    for svtype in ("DEL", "DUP"):
        calls = by_key.get((contig, svtype), [])
        if not calls:
            continue
        calls.sort(key=lambda m: m.start)
        for cl in cluster_calls([(m, svtype) for m in calls]):
            merged.append((cl.start, cl))
    merged.sort(key=lambda t: t[0])
    sites = []
    for _pos, cl in merged:
        rescued = not any(m.passing for m in cl.members)
        # best SQ per carrier: rescue criterion, and SQMAX/SQMEAN
        # carried to the cohort VCF (the per-sample SQ is otherwise
        # lost past the merge)
        best_sq: dict[str, float] = {}
        for m in cl.members:
            best_sq[m.sample] = max(best_sq.get(m.sample, 0.0), m.sq)
        if rescued and (
            len(best_sq) < rescue_min_samples
            or sum(best_sq.values()) < rescue_total_sq
        ):
            continue
        # member calls from VCFs written before vcfout clamped END
        # overhang the contig end by up to one depth bin
        start = cl.start
        end = clamp_end(start, cl.end, length)
        info = [
            ("NMEMBERS", len(cl.members)),
            ("NSAMPLES", len({m.sample for m in cl.members})),
            ("SPREAD", cl.spread),
            ("PRECISE_FRAC", round(
                sum(1 for m in cl.members if m.precise) / len(cl.members), 3
            )),
            ("SQMAX", round(max(best_sq.values()), 2)),
            ("SQMEAN", round(sum(best_sq.values()) / len(best_sq), 2)),
        ]
        if rescued:
            info.append(("RESCUE", "RECURRENT"))
        sites.append((start, end, cl.svtype, info))
    return sites


def run_merge(
    vcf_paths: list[str | Path],
    contigs: list[Contig],
    out_vcf: str | Path,
    min_gq: int = 0,
    pass_only: bool = True,
    rescue_min_sq: float | None = 1.0,
    rescue_min_samples: int = 2,
    rescue_total_sq: float = 8.0,
    jobs: int = 1,
) -> int:
    """Merge per-sample VCFs into a sites-only VCF; returns the site count.

    A cluster with no passing member is still emitted (INFO/RESCUE=RECURRENT)
    when >= rescue_min_samples distinct carriers support it and their summed
    best SQ reaches rescue_total_sq — treating per-sample calibrated scores
    as accumulating evidence for recurrent sub-threshold sites.

    ``jobs`` > 1 reads the VCFs and clusters the contigs in parallel (4 h on
    one core for 2670 samples); input order is kept at both steps, so the
    output is identical to the serial run.
    """
    global _BY_KEY
    by_key: dict[tuple[str, str], list[MemberCall]] = {}
    read = partial(
        read_sample_calls, min_gq=min_gq, pass_only=pass_only,
        rescue_min_sq=rescue_min_sq,
    )
    for calls_of in pmap(read, vcf_paths, jobs):
        for key, calls in calls_of.items():
            by_key.setdefault(key, []).extend(calls)

    header = pysam.VariantHeader()
    header.add_line(f"##source=vedette-merge-{__version__}")
    add_coords_line(header)
    for c in contigs:
        header.contigs.add(c.name, length=c.length)
    header.add_line('##ALT=<ID=DEL,Description="Deletion">')
    header.add_line('##ALT=<ID=DUP,Description="Duplication">')
    header.info.add("SVTYPE", "1", "String", "Type of structural variant")
    header.info.add("END", "1", "Integer", "End position")
    header.info.add("SVLEN", "1", "Integer", "Signed length")
    header.info.add("NMEMBERS", "1", "Integer", "Number of clustered per-sample calls")
    header.info.add("NSAMPLES", "1", "Integer", "Number of distinct carrier samples")
    header.info.add("SPREAD", "1", "Integer", "Breakpoint spread among members (bp)")
    header.info.add("PRECISE_FRAC", "1", "Float", "Fraction of members with precise breakpoints")
    header.info.add("SQMAX", "1", "Float", "Highest per-sample SQ among the clustered calls")
    header.info.add("SQMEAN", "1", "Float", "Mean of the per-sample best SQ among carriers")
    header.info.add(
        "RESCUE", "1", "String",
        "Site discovery mode when no member passed the single-sample filter "
        "(RECURRENT: rescued by cross-sample recurrence of sub-threshold calls)",
    )

    contig_order = {c.name: i for i, c in enumerate(contigs)}
    contig_length = {c.name: c.length for c in contigs}
    work = [
        (contig, contig_length.get(contig), rescue_min_samples, rescue_total_sq)
        for contig in sorted(
            {k[0] for k in by_key}, key=lambda c: contig_order.get(c, 1_000_000)
        )
    ]
    # workers read the calls from this module global, shared copy-on-write
    # by fork instead of pickled to each of them
    # (all sites computed before the output opens: no htslib handle may be
    # open in the parent when the pool forks)
    _BY_KEY = by_key
    try:
        per_contig = list(pmap(_contig_sites, work, jobs))
    finally:
        _BY_KEY = None
    n_sites = 0
    out_path = str(out_vcf)
    mode = "wz" if out_path.endswith(".gz") else "w"
    with pysam.VariantFile(out_path, mode, header=header) as out:
        for (contig, *_), sites in zip(work, per_contig, strict=True):
            for start, end, svtype, info in sites:
                n_sites += 1
                rec = sv_record(out, contig, start, end, svtype, f"vedette_site_{n_sites}")
                for k, v in info:
                    rec.info[k] = v
                out.write(rec)
    if mode == "wz":
        pysam.tabix_index(out_path, preset="vcf", force=True)
    logger.info("merge: %d input VCFs -> %d sites", len(vcf_paths), n_sites)
    return n_sites
