"""Fair DUP evaluation: isolate the tandem-duplication subset of an
INS-represented truth set.

Assembly-based truth sets (GIAB v5.0q) encode duplications as sequence-explicit
insertions. Most truth insertions are novel sequence (mobile elements, etc.)
that a DEL/DUP caller cannot represent — evaluating <DUP> calls against ALL
insertions understates recall by construction. An insertion at position p with
sequence S is a TANDEM DUPLICATION when S duplicates the adjacent reference:
S ~ ref[p : p+|S|] (or the upstream flank), possibly with small indel jitter.

Classification: exact-slice identity first (cheap), then a fit alignment of S
against ref[p-|S| : p+2|S|] for the jittered cases (bounded size).
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pysam

from ..breakpoints.resolve import _fit_align, _seq_to_arr

logger = logging.getLogger(__name__)


def _slice_identity(s: np.ndarray, ref: np.ndarray) -> float:
    n = min(len(s), len(ref))
    if n == 0:
        return 0.0
    return float((s[:n] == ref[:n]).mean())


def is_tandem_ins(
    seq: str,
    contig: str,
    pos: int,
    fasta: pysam.FastaFile,
    min_identity: float = 0.80,
    max_align_len: int = 5000,
) -> bool:
    s = _seq_to_arr(seq.encode())
    length = len(s)
    contig_len = fasta.get_reference_length(contig)

    down = _seq_to_arr(
        fasta.fetch(contig, pos, min(pos + length, contig_len)).upper().encode()
    )
    up_start = max(pos - length, 0)
    up = _seq_to_arr(fasta.fetch(contig, up_start, pos).upper().encode())
    # cheap path: direct slice identity against either flank
    if _slice_identity(s, down) >= min_identity:
        return True
    if len(up) and _slice_identity(s[-len(up):], up) >= min_identity:
        return True
    # multi-copy tandem: the first repeat unit duplicates downstream sequence
    for unit in (length // 2, length // 3):
        if unit >= 30 and _slice_identity(s[:unit], down[:unit]) >= min_identity:
            return True

    if length > max_align_len:
        return False
    # jittered path: fit-align S in a +/- |S| window around the insertion point
    win_start = max(pos - length, 0)
    window = _seq_to_arr(
        fasta.fetch(contig, win_start, min(pos + 2 * length, contig_len)).upper().encode()
    )
    if len(window) < 30:
        return False
    score, _rs, _re = _fit_align(s, window)
    return score >= 2 * min_identity * length


def make_tandem_truth(
    truth_vcf: str | Path,
    fasta_path: str | Path,
    out_vcf: str | Path,
    min_size: int = 50,
    min_identity: float = 0.80,
) -> dict[str, int]:
    """Write the tandem-INS (+ explicit DUP) subset of a truth VCF.

    Deletions are dropped (this subset is the DUP-space denominator).
    """
    counts = {"ins_total": 0, "tandem": 0, "other_ins": 0}
    with (
        pysam.VariantFile(str(truth_vcf)) as vin,
        pysam.FastaFile(str(fasta_path)) as fasta,
        pysam.VariantFile(str(out_vcf), "wz", header=vin.header) as vout,
    ):
        for rec in vin:
            if not rec.alts:
                continue
            ref_len, alt_len = len(rec.ref), len(rec.alts[0])
            svtype = rec.info.get("SVTYPE", "")
            if svtype.startswith("DUP"):
                vout.write(rec)
                counts["tandem"] += 1
                continue
            if alt_len - ref_len < min_size or "<" in rec.alts[0]:
                continue  # not a sequence-explicit insertion of >= min_size
            counts["ins_total"] += 1
            inserted = rec.alts[0][ref_len:] if rec.alts[0].startswith(rec.ref) \
                else rec.alts[0][1:]
            if is_tandem_ins(inserted, rec.contig, rec.stop, fasta, min_identity):
                vout.write(rec)
                counts["tandem"] += 1
            else:
                counts["other_ins"] += 1
    pysam.tabix_index(str(out_vcf), preset="vcf", force=True)
    logger.info(
        "tandem truth: %(tandem)d tandem / %(ins_total)d INS>=min_size "
        "(%(other_ins)d non-tandem dropped)", counts,
    )
    return counts
