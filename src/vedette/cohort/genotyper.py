"""`vedette genotype`: force-genotype cohort sites in one sample, evidence-store only.

Reads the sites VCF from `vedette merge` and this sample's evidence store —
never the CRAM — and emits one genotype row per site. Sample-parallel by
construction: run one job per sample, then assemble with `vedette cohort-filter`.
"""

from __future__ import annotations

import csv
import logging
from pathlib import Path

import numpy as np
import pysam

from ..config import Config
from ..extract.store import PE_DEL, PE_DUP, SR_LEFT_CLIP, SR_RIGHT_CLIP, EvidenceStore
from ..genotype.likelihood import ContigData, Genotyper
from ..segment.hmm import fit_sigma0
from ..variants import SOURCE_BOTH, SOURCE_DEPTH, CNVCall

logger = logging.getLogger(__name__)

GENOTYPE_COLUMNS = [
    "site_id", "contig", "start", "end", "svtype",
    "sample", "gt", "gq", "cn", "cnq", "rd", "pe", "sr", "cig", "baf_n", "baf_dev",
]

_SR_TOL = 100  # bp tolerance around site breakpoints for SR support


class StoreGenotyper:
    def __init__(self, store: EvidenceStore, config: Config):
        self.store = store
        self.config = config
        self.genotyper = Genotyper(config.genotype, config.segment.max_cn)
        self.mean_depth = float(store.meta["mean_depth"])
        self.coverage = float(store.meta.get("coverage", self.mean_depth))
        self.ploidy: dict[str, int] = {
            k: int(v) for k, v in dict(store.meta.get("ploidy", {})).items()
        }
        insert_mu = float(store.meta.get("insert_mu", 0.0))
        insert_sigma = float(store.meta.get("insert_sigma", 0.0))
        self.pe_window = int(insert_mu + 3 * insert_sigma) if insert_mu else 1000
        self._data: dict[str, ContigData] = {}
        self._has_junctions: dict[str, bool] = {}
        self._junc: dict[str, dict] = {}

    def _contig_data(self, contig: str) -> ContigData | None:
        if contig not in self._data:
            if not self.store.has_contig(contig):
                self._data[contig] = None
            else:
                _raw, corr, mask = self.store.depth(contig)
                ratio = corr / self.mean_depth
                length = len(ratio) * self.store.bin_size
                baf = self.store.baf(contig, 0, length)
                self._data[contig] = ContigData(
                    ratio=ratio,
                    mask=mask,
                    baf_pos=baf["pos"].astype(np.int64),
                    baf_vaf=baf["vaf"],
                    bin_size=self.store.bin_size,
                    ploidy=self.ploidy.get(contig, 2),
                    mean_depth=self.coverage,
                    sigma0=fit_sigma0(ratio, mask),
                )
                self._has_junctions[contig] = self.store.has_junction_evidence(contig)
        return self._data[contig]

    def _junction_arrays(self, contig: str) -> dict:
        """Whole-contig PE/SR/CIG arrays, loaded once: per-site ranged zarr
        reads decompress the same chunks over and over (profiled at ~9 sites/s
        on 238k cohort sites — days per cohort); in-memory searchsorted slices
        are ~100x faster for a few hundred MB per contig at most."""
        if contig not in self._junc:
            length = 2**62  # whole contig
            pe = self.store.pe(contig, 0, length)
            sr = self.store.sr(contig, 0, length)
            cig = self.store.cig(contig, 0, length)
            self._junc[contig] = {
                "pe_pos1": pe["pos1"].astype(np.int64),
                "pe_pos2": pe["pos2"].astype(np.int64),
                "pe_orient": pe["orient"],
                "sr_pos": sr["pos"].astype(np.int64),
                "sr_side": sr["side"],
                "cig_pos": cig["pos"].astype(np.int64),
                "cig_svlen": cig["svlen"].astype(np.int64),
            }
        return self._junc[contig]

    def _junction_support(
        self, contig: str, start: int, end: int, svtype: str
    ) -> tuple[int, int, int]:
        """(pe, sr, cig) support for the adjacency start->end in this sample."""
        j = self._junction_arrays(contig)
        w = self.pe_window
        orient = PE_DEL if svtype == "DEL" else PE_DUP
        n_pe = 0
        lo = int(np.searchsorted(j["pe_pos1"], max(start - w, 0), side="left"))
        hi = int(np.searchsorted(j["pe_pos1"], start + w, side="right"))
        if hi > lo:
            pos2 = j["pe_pos2"][lo:hi]
            sel = (
                (j["pe_orient"][lo:hi] == orient)
                & (pos2 >= end - w) & (pos2 <= end + w)
            )
            n_pe = int(sel.sum())
        n_sr = 0
        lo = int(np.searchsorted(j["sr_pos"], max(start - _SR_TOL, 0), side="left"))
        hi = int(np.searchsorted(j["sr_pos"], start + _SR_TOL, side="right"))
        if hi > lo:
            n_sr += int((j["sr_side"][lo:hi] == SR_RIGHT_CLIP).sum())
        lo = int(np.searchsorted(j["sr_pos"], max(end - _SR_TOL, 0), side="left"))
        hi = int(np.searchsorted(j["sr_pos"], end + _SR_TOL, side="right"))
        if hi > lo:
            n_sr += int((j["sr_side"][lo:hi] == SR_LEFT_CLIP).sum())
        n_cig = 0
        if svtype == "DEL":
            lo = int(np.searchsorted(j["cig_pos"], max(start - _SR_TOL, 0), side="left"))
            hi = int(np.searchsorted(j["cig_pos"], start + _SR_TOL, side="right"))
            if hi > lo:
                size = max(end - start, 1)
                sel = np.abs(j["cig_svlen"][lo:hi] - size) <= max(10, 0.2 * size)
                n_cig = int(sel.sum())
        return n_pe, n_sr, n_cig

    def genotype_site(
        self, site_id: str, contig: str, start: int, end: int, svtype: str
    ) -> dict | None:
        data = self._contig_data(contig)
        if data is None:
            return None
        if data.ploidy == 0:
            return None  # chrY in a female sample
        has_junc = self._has_junctions.get(contig, False)
        pe, sr, cig = (
            self._junction_support(contig, start, end, svtype) if has_junc else (0, 0, 0)
        )
        call = CNVCall(
            contig=contig, start=start, end=end, svtype=svtype,
            source=SOURCE_BOTH if has_junc else SOURCE_DEPTH,
            pe_support=pe, sr_support=sr, cig_support=cig,
        )
        self.genotyper.genotype(call, data)
        gt = "./." if call.gt == (-1, -1) else f"{call.gt[0]}/{call.gt[1]}"
        return {
            "site_id": site_id, "contig": contig, "start": start, "end": end,
            "svtype": svtype, "sample": self.store.sample,
            "gt": gt, "gq": call.gq, "cn": call.cn,
            "cnq": round(call.cn_posterior, 4) if call.cn >= 0 else "",
            "rd": round(call.mean_ratio, 4) if not np.isnan(call.mean_ratio) else "",
            "pe": pe, "sr": sr, "cig": cig, "baf_n": call.baf_n,
            "baf_dev": round(call.baf_dev, 4) if not np.isnan(call.baf_dev) else "",
        }


def run_genotype(
    sites_vcf: str | Path,
    store_path: str | Path,
    out_tsv: str | Path,
    config: Config,
) -> int:
    store = EvidenceStore(store_path)
    gtr = StoreGenotyper(store, config)
    n = 0
    with (
        pysam.VariantFile(str(sites_vcf)) as sites,
        open(out_tsv, "w", newline="") as fh,
    ):
        writer = csv.DictWriter(fh, fieldnames=GENOTYPE_COLUMNS, delimiter="\t")
        writer.writeheader()
        for rec in sites:
            row = gtr.genotype_site(
                rec.id, rec.contig, rec.start, rec.stop, rec.info["SVTYPE"]
            )
            if row is not None:
                writer.writerow(row)
                n += 1
    logger.info("genotype: %s -> %d site genotypes", store.sample, n)
    return n
