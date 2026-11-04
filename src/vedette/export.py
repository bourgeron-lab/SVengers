"""Machine-readable family export, written next to each events page.

``families/<FID>/``:

- ``<FID>.vedette.tsv``: one row per row of ``events.html`` (cohort sites
  collapsed into events, sub-threshold, mosaic and chromosomal rows), with
  typed, coded columns (format: docs/EXPORT.md);
- ``<FID>.vedette.qc.tsv``: one row per member, the per-sample QC behind the
  row alerts;
- ``<FID>.vedette.json``: schema and caller versions, member order.

Coordinates: ``start``/``end`` are the first and last affected bases,
1-based inclusive (as IGV and UCSC show them); ``vcf_pos`` = ``start - 1`` is
the VCF POS (padding base). A missing value is an empty field.
"""

from __future__ import annotations

import csv
import json
import re
import time
from collections import Counter
from pathlib import Path

from . import __version__

# bump when a column is renamed, removed or changes meaning (adding one: no bump)
SCHEMA_VERSION = 1

ROW_COLUMNS = [
    "row_id", "kind", "site_ids", "contig", "start", "end", "vcf_pos", "svtype", "size",
    "filter", "rescue", "n_sites", "sqmax", "cohort_carrier_frac", "ref_carrier_frac",
    "lcr_ovl", "segdup_ovl", "support", "own_call", "mosaic_cn",
    "impact_score", "impact_level", "impact_reason", "impact",
    "exonic_genes", "exonic_gene_ids", "genic_genes", "genic_gene_ids",
    "carriers", "qc_alerts",
    # display / compatibility with the pre-v1.4.0 export: prefer the columns above
    "genes", "inheritance",
]
MEMBER_PREFIXES = ("GT", "CN", "GQ", "INH", "DNQ")
QC_COLUMNS = [
    "sample", "wave_score", "pass_sites_carried", "carrier_excess", "denovo_candidates",
    "karyotype", "flags", "flags_text",
]

IMPACT_LEVELS = {"élevé": "high", "modéré": "moderate", "faible": "low", "aucun gène codant": "none"}
_IMPACT_RE = re.compile(
    r"^\s*([+-]\d+(?:\.\d+)?) — (" + "|".join(IMPACT_LEVELS) + r")(?: · (.*))?$"
)
_OWN = {"PASS": "pass", "filtré": "filtered", "": "none", "LowSQ": "lowsq"}
_SUPPORT = {"fort": "strong", "faible": "weak"}
_GT = {"./.", "0/0", "0/1", "1/1"}


def impact_fields(label: str) -> tuple[str, str, str]:
    """(score, level, reason) of an impact label written by
    ``_GeneIndex.impact`` ("+0.30 — faible · DEL exonique SAMD11 (...)")."""
    m = _IMPACT_RE.match(label or "")
    if not m:
        return "", "", ""
    return f"{float(m.group(1)):.2f}", IMPACT_LEVELS[m.group(2)], m.group(3) or ""


def row_kind(e: dict) -> str:
    if e.get("sub"):
        return "subthreshold"
    return {"mosaïque": "mosaic", "chromosomique": "chromosomal"}.get(e.get("own") or "", "site")


def _inheritance(labels: list[tuple[str, str]]) -> tuple[str, str]:
    """(code, dnq) of one child over the sites of a row: the label most of
    them agree on; the best DNQ of its de novo sites."""
    if not labels:
        return "", ""
    top = Counter(inh for inh, _ in labels).most_common(1)[0][0]
    dnq = [int(d) for inh, d in labels if inh == "denovo" and d]
    return top, (str(max(dnq)) if dnq else "")


def family_rows(members: list[str], events: list[dict], gene_lists_fn=None) -> list[dict]:
    rows = []
    for e in events:
        kind = row_kind(e)
        score, level, reason = impact_fields(e.get("impact") or "")
        exonic, genic = ([], [])
        if gene_lists_fn is not None:
            exonic, genic = gene_lists_fn(e["contig"], e["pos"] - 1, e["end"])
        row = {
            "row_id": e.get("id") or "",
            "kind": kind,
            "site_ids": ";".join(e.get("_ids") or [e["id"]]) if kind == "site" and e.get("id") else "",
            "contig": e["contig"], "start": e["pos"], "end": e["end"], "vcf_pos": e["pos"] - 1,
            "svtype": e["svtype"], "size": e["end"] - e["pos"] + 1,
            "filter": e.get("filter") or "", "rescue": e.get("rescue") or "",
            "n_sites": e.get("n_sites", 1) if kind == "site" else "",
            "sqmax": e.get("sqmax") or "", "cohort_carrier_frac": e.get("cfrac") or "",
            "ref_carrier_frac": e.get("refcf") or "",
            "lcr_ovl": e.get("lcr") or "", "segdup_ovl": e.get("segdup") or "",
            "support": _SUPPORT.get(e.get("support") or "", ""),
            "own_call": "" if kind in ("mosaic", "chromosomal") else _OWN.get(e.get("own") or "", ""),
            "mosaic_cn": f"{e['mosaic_cn']:.2f}" if e.get("mosaic_cn") is not None else "",
            "impact_score": score, "impact_level": level, "impact_reason": reason,
            "impact": e.get("impact") or "",
            "exonic_genes": ";".join(n for n, _ in exonic),
            "exonic_gene_ids": ";".join(i for _, i in exonic),
            "genic_genes": ";".join(n for n, _ in genic),
            "genic_gene_ids": ";".join(i for _, i in genic),
            "carriers": ";".join(e.get("carriers") or []),
            "qc_alerts": e.get("qc") or "",
            "genes": e.get("genes") or "",
        }
        inh_by = e.get("_inh_by") or {}
        inherit = []
        for m, cell in zip(members, e.get("gt") or [("", "", "")] * len(members), strict=False):
            gt, gq, cn = cell
            integer = gt in _GT
            row[f"GT_{m}"] = gt if integer else ""
            row[f"GQ_{m}"] = gq if integer and gq else ""
            row[f"CN_{m}"] = cn if integer and cn else ""
            inh, dnq = _inheritance(inh_by.get(m, []))
            row[f"INH_{m}"], row[f"DNQ_{m}"] = inh, dnq
            if inh:
                inherit.append(f"{m}:{inh}")
        row["inheritance"] = "; ".join(inherit)
        rows.append(row)
    return rows


def write_family_export(fam_dir: Path, fid: str, members: list[str], events: list[dict],
                        gene_lists_fn=None, qc_rows: list[dict] | None = None,
                        meta: dict | None = None) -> None:
    columns = ROW_COLUMNS + [f"{p}_{m}" for m in members for p in MEMBER_PREFIXES]
    rows = family_rows(members, events, gene_lists_fn)
    rows.sort(key=lambda r: (r["contig"], r["start"], r["end"]))
    _write_tsv(fam_dir / f"{fid}.vedette.tsv", columns, rows)
    _write_tsv(fam_dir / f"{fid}.vedette.qc.tsv", QC_COLUMNS, qc_rows or [])
    info = {
        "schema_version": SCHEMA_VERSION,
        "vedette_version": __version__,
        "family_id": fid,
        "members": members,
        "rows": len(rows),
        "coordinates": "start/end: first/last affected base, 1-based inclusive; "
                       "vcf_pos = start - 1 (VCF POS, padding base)",
        "missing_value": "",
        "generated": time.strftime("%Y-%m-%dT%H:%M:%S"),
        **(meta or {}),
    }
    tmp = fam_dir / f".{fid}.vedette.json.partial"
    tmp.write_text(json.dumps(info, ensure_ascii=False, indent=1) + "\n")
    tmp.replace(fam_dir / f"{fid}.vedette.json")


def qc_row(sample: str, info: dict, flags_text: list[str]) -> dict:
    """One member row of ``<FID>.vedette.qc.tsv``; ``info``: wave_score,
    pass_sites_carried, carrier_excess, denovo_candidates, karyotype, flags
    (codes)."""
    def num(v, fmt):
        return "" if v is None else format(v, fmt)
    return {
        "sample": sample,
        "wave_score": num(info.get("wave_score"), ".3f"),
        "pass_sites_carried": num(info.get("pass_sites_carried"), "d"),
        "carrier_excess": num(info.get("carrier_excess"), ".2f"),
        "denovo_candidates": num(info.get("denovo_candidates"), "d"),
        "karyotype": info.get("karyotype") or "",
        "flags": ";".join(info.get("flags") or []),
        "flags_text": " · ".join(flags_text),
    }


def _write_tsv(path: Path, columns: list[str], rows: list[dict]) -> None:
    tmp = path.with_name(f".{path.name}.partial")
    with open(tmp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", lineterminator="\n",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    tmp.replace(path)
