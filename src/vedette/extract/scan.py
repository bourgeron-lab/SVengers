"""Single-pass CRAM evidence extraction (the self-contained path).

One coordinate-sorted sweep per contig (contigs run in parallel processes)
collecting, per record and with minimal Python work:

- binned depth: read-start counts per bin (MAPQ >= min_mapq_depth, primary,
  non-duplicate) — scale-free, the GC/ratio normalization downstream is
  invariant to counting reads vs coverage;
- discordant pairs: orientation class DEL/DUP/INV/INTER (recorded once, from
  the leftmost read of the pair);
- split reads: soft-clip positions >= min_clip_len with decent clip base
  quality;
- CIGAR deletions: D operations >= min_cigar_sv bp in primary alignments —
  the dominant (often only) evidence for 50-300 bp deletions, which bwa
  aligns THROUGH rather than clipping (raw-candidate recall for that size
  class was 29% before this channel existed).

The insert-size model is fit beforehand on a sample of proper pairs.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from multiprocessing import get_context

import numpy as np
import pysam

from ..config import ExtractConfig
from .store import PE_DEL, PE_DUP, PE_INTER, PE_INV, SR_LEFT_CLIP, SR_RIGHT_CLIP, InsertModel

logger = logging.getLogger(__name__)

# SAM flags
_PAIRED = 0x1
_UNMAPPED = 0x4
_MATE_UNMAPPED = 0x8
_REVERSE = 0x10
_MATE_REVERSE = 0x20
_SECONDARY = 0x100
_QCFAIL = 0x200
_DUP = 0x400
_SUPPLEMENTARY = 0x800

_SKIP_ALWAYS = _UNMAPPED | _SECONDARY | _QCFAIL | _DUP


@dataclass
class ContigEvidence:
    contig: str
    depth: np.ndarray                     # float32 read-start counts per bin
    pe: dict[str, np.ndarray] = field(default_factory=dict)
    sr: dict[str, np.ndarray] = field(default_factory=dict)
    sr_seq: list[bytes] = field(default_factory=list)  # clip sequences, aligned with sr rows
    cig: dict[str, np.ndarray] = field(default_factory=dict)  # large CIGAR deletions
    ins: dict[str, np.ndarray] = field(default_factory=dict)  # binned insert-size track
    n_reads: int = 0


def fit_insert_model(
    cram: str, fasta: str, contigs: list[str], config: ExtractConfig
) -> InsertModel:
    """Robust insert-size model (median/MAD) from proper pairs."""
    sizes: list[int] = []
    read_lengths: list[int] = []
    target = config.insert_fit_pairs
    # mode "r" auto-detects BAM vs CRAM; reference_filename is ignored for BAM
    with pysam.AlignmentFile(cram, "r", reference_filename=fasta) as af:
        for contig in contigs:
            for rec in af.fetch(contig):
                flag = rec.flag
                if flag & (_SKIP_ALWAYS | _SUPPLEMENTARY):
                    continue
                if not (flag & _PAIRED) or not rec.is_proper_pair:
                    continue
                tlen = rec.template_length
                if 0 < tlen < 100_000:
                    sizes.append(tlen)
                    if len(read_lengths) < 100_000:
                        read_lengths.append(rec.query_length)
                    if len(sizes) >= target:
                        break
            if len(sizes) >= target:
                break
    if len(sizes) < 10_000:
        raise ValueError(f"{cram}: only {len(sizes)} proper pairs for the insert model")
    arr = np.asarray(sizes)
    mu = float(np.median(arr))
    sigma = float(np.median(np.abs(arr - mu)) * 1.4826)
    read_length = float(np.median(read_lengths)) if read_lengths else 150.0
    logger.info(
        "Insert model: mu=%.1f sigma=%.1f readlen=%.0f (n=%d)", mu, sigma, read_length, len(arr)
    )
    return InsertModel(mu=mu, sigma=sigma, n_pairs=len(arr), read_length=read_length)


def _mean_clip_qual(qualities, start: int, end: int) -> float:
    if qualities is None:
        return 40.0  # no base qualities stored: do not penalize
    q = qualities[start:end]
    return float(sum(q) / max(len(q), 1))


def scan_contig(
    cram: str,
    fasta: str,
    contig: str,
    contig_length: int,
    insert: InsertModel,
    config: ExtractConfig,
) -> ContigEvidence:
    bin_size = config.bin_size
    n_bins = -(-contig_length // bin_size)
    depth = np.zeros(n_bins, dtype=np.float32)

    max_insert = insert.discordant_threshold(config.insert_z_discordant)
    min_mapq_ev = config.min_mapq_evidence
    min_mapq_depth = config.min_mapq_depth
    min_clip = config.min_clip_len
    min_clip_q = config.min_clip_baseq
    max_span = config.max_discordant_insert

    pe_pos1: list[int] = []
    pe_pos2: list[int] = []
    pe_mate_tid: list[int] = []
    pe_orient: list[int] = []
    pe_mapq: list[int] = []
    sr_pos: list[int] = []
    sr_side: list[int] = []
    sr_mapq: list[int] = []
    sr_cliplen: list[int] = []
    sr_seq: list[bytes] = []
    cig_pos: list[int] = []
    cig_len: list[int] = []
    min_cigar = config.min_cigar_sv
    merge_gap = config.cigar_merge_gap

    # local insert-size track: mean template length of spanning proper pairs
    # per small bin — the collective shift detects sub-insert deletions that
    # leave no clean junction (VNTR contractions)
    ins_bin = config.insert_bin_size
    n_ins_bins = -(-contig_length // ins_bin)
    ins_sum = np.zeros(n_ins_bins, dtype=np.float64)
    ins_n = np.zeros(n_ins_bins, dtype=np.uint32)
    buf_mid: list[int] = []
    buf_tlen: list[int] = []

    def _flush_insert_buffer() -> None:
        if buf_mid:
            bins = np.asarray(buf_mid, np.int64) // ins_bin
            np.add.at(ins_sum, bins, np.asarray(buf_tlen, np.float64))
            np.add.at(ins_n, bins, 1)
            buf_mid.clear()
            buf_tlen.clear()

    n_reads = 0
    with pysam.AlignmentFile(cram, "r", reference_filename=fasta) as af:
        tid = af.get_tid(contig)
        for rec in af.fetch(contig):
            flag = rec.flag
            if flag & _SKIP_ALWAYS:
                continue
            mapq = rec.mapping_quality
            is_supp = bool(flag & _SUPPLEMENTARY)
            pos = rec.reference_start
            n_reads += 1

            # ---- depth (primary only)
            if not is_supp and mapq >= min_mapq_depth:
                depth[pos // bin_size] += 1.0

            if mapq < min_mapq_ev:
                continue

            # ---- split reads: soft clips on primary alignments
            if not is_supp:
                cigar = rec.cigartuples
                if cigar:
                    op0, len0 = cigar[0]
                    if op0 == 4 and len0 >= min_clip:  # leading soft clip
                        if _mean_clip_qual(rec.query_qualities, 0, len0) >= min_clip_q:
                            sr_pos.append(pos)
                            sr_side.append(SR_LEFT_CLIP)
                            sr_mapq.append(mapq)
                            sr_cliplen.append(min(len0, 65535))
                            seq = rec.query_sequence
                            sr_seq.append(seq[:len0].encode() if seq else b"")
                    # large CIGAR deletions (op 2 = D): walk ref-consuming ops.
                    # Nearby D ops within one read are summed — aligners often
                    # fragment one VNTR-contraction deletion into several gaps.
                    if len(cigar) > 2:
                        ref_cursor = pos
                        run_start, run_len, run_end = -1, 0, -1
                        for op, oplen in cigar:
                            if op == 2:  # D
                                if run_start >= 0 and ref_cursor - run_end <= merge_gap:
                                    run_len += oplen
                                else:
                                    if run_len >= min_cigar:
                                        cig_pos.append(run_start)
                                        cig_len.append(run_len)
                                    run_start, run_len = ref_cursor, oplen
                                ref_cursor += oplen
                                run_end = ref_cursor
                            elif op in (0, 3, 7, 8):  # M, N, =, X consume reference
                                ref_cursor += oplen
                        if run_len >= min_cigar:
                            cig_pos.append(run_start)
                            cig_len.append(run_len)
                    opn, lenn = cigar[-1]
                    if opn == 4 and lenn >= min_clip:  # trailing soft clip
                        qlen = rec.query_length
                        if (
                            _mean_clip_qual(rec.query_qualities, qlen - lenn, qlen)
                            >= min_clip_q
                        ):
                            sr_pos.append(rec.reference_end)
                            sr_side.append(SR_RIGHT_CLIP)
                            sr_mapq.append(mapq)
                            sr_cliplen.append(min(lenn, 65535))
                            seq = rec.query_sequence
                            sr_seq.append(seq[qlen - lenn :].encode() if seq else b"")

            # ---- local insert-size track (proper pairs, leftmost read)
            if not is_supp and flag & _PAIRED and rec.is_proper_pair:
                tlen = rec.template_length
                if 0 < tlen < 100_000:
                    buf_mid.append(pos + tlen // 2)
                    buf_tlen.append(tlen)
                    if len(buf_mid) >= 2_000_000:
                        _flush_insert_buffer()

            # ---- discordant pairs (record once, from the leftmost read)
            if is_supp or not (flag & _PAIRED) or flag & _MATE_UNMAPPED:
                continue
            mate_tid = rec.next_reference_id
            mate_pos = rec.next_reference_start
            if mate_tid != tid:
                if tid < mate_tid:  # record the pair once, from the lower tid
                    pe_pos1.append(pos)
                    pe_pos2.append(mate_pos)
                    pe_mate_tid.append(mate_tid)
                    pe_orient.append(PE_INTER)
                    pe_mapq.append(mapq)
                continue
            if pos > mate_pos or (pos == mate_pos and flag & _REVERSE):
                continue  # rightmost read of the pair: already recorded
            span = mate_pos - pos
            if span > max_span:
                continue
            self_rev = bool(flag & _REVERSE)
            mate_rev = bool(flag & _MATE_REVERSE)
            if self_rev == mate_rev:
                orient = PE_INV
            elif self_rev and not mate_rev:
                orient = PE_DUP  # everted: leftmost read is reverse (tandem dup)
            elif span > max_insert:
                orient = PE_DEL  # normal orientation, stretched insert
            else:
                continue  # concordant
            pe_pos1.append(pos)
            pe_pos2.append(mate_pos)
            pe_mate_tid.append(mate_tid)
            pe_orient.append(orient)
            pe_mapq.append(mapq)

    return ContigEvidence(
        contig=contig,
        depth=depth,
        pe={
            "pos1": np.asarray(pe_pos1, np.uint32),
            "pos2": np.asarray(pe_pos2, np.uint32),
            "mate_tid": np.asarray(pe_mate_tid, np.int32),
            "orient": np.asarray(pe_orient, np.uint8),
            "mapq": np.asarray(pe_mapq, np.uint8),
        },
        sr={
            "pos": np.asarray(sr_pos, np.uint32),
            "side": np.asarray(sr_side, np.uint8),
            "mapq": np.asarray(sr_mapq, np.uint8),
            "cliplen": np.asarray(sr_cliplen, np.uint16),
        },
        sr_seq=sr_seq,
        cig={
            "pos": np.asarray(cig_pos, np.uint32),
            "svlen": np.asarray(cig_len, np.uint32),
        },
        ins=(_flush_insert_buffer() or {"sum": ins_sum.astype(np.float32), "n": ins_n}),
        n_reads=n_reads,
    )


def _scan_worker(args) -> ContigEvidence:
    cram, fasta, contig, length, insert, config = args
    return scan_contig(cram, fasta, contig, length, insert, config)


def scan_cram(
    cram: str,
    fasta: str,
    contigs: dict[str, int],
    config: ExtractConfig,
    threads: int = 8,
) -> tuple[InsertModel, dict[str, ContigEvidence]]:
    """Fit the insert model, then scan all contigs (in parallel)."""
    fit_contigs = [c for c in contigs if c not in ("chrY", "chrM")][:4]
    insert = fit_insert_model(cram, fasta, fit_contigs, config)

    jobs = [
        (cram, fasta, contig, length, insert, config) for contig, length in contigs.items()
    ]
    results: dict[str, ContigEvidence] = {}
    if threads <= 1 or len(jobs) == 1:
        for job in jobs:
            ev = _scan_worker(job)
            results[ev.contig] = ev
            logger.info("scan %s: %d reads", ev.contig, ev.n_reads)
    else:
        # ProcessPoolExecutor, not multiprocessing.Pool: when a worker is
        # OOM-killed, Pool.imap_unordered blocks forever on the result queue
        # (observed: 6h at 0% CPU on a deep genome), while the executor raises
        # BrokenProcessPool so the task fails fast and SLURM retries with more
        # memory. spawn context: htslib handles must not cross fork.
        from concurrent.futures import ProcessPoolExecutor, as_completed

        ctx = get_context("spawn")
        with ProcessPoolExecutor(
            max_workers=min(threads, len(jobs)), mp_context=ctx
        ) as pool:
            futures = [pool.submit(_scan_worker, job) for job in jobs]
            for fut in as_completed(futures):
                ev = fut.result()
                results[ev.contig] = ev
                logger.info("scan %s: %d reads", ev.contig, ev.n_reads)
    return insert, results
