"""Self-contained HTML cohort report + per-family workspaces.

``vedette report`` reads a cohort run directory (the Nextflow wrapper layout:
cohort/cohort.vcf.gz, cohort/mendel_*.tsv, cohort/trio_sweep.tsv, genotypes/)
and writes:

- ``report.html`` — one static page, figures inlined as base64 PNG: run
  summary, size/AF spectra, per-sample table, MVR/transmission sweep,
  per-family sections with de novo candidates and violations, rare sites;
- ``families/<FID>/`` — one workspace per family: the family-subset VCF
  (carriers of at least one member) and its Mendel rows.
"""

from __future__ import annotations

import base64
import csv
import html
import io
import logging
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pysam

from . import __version__
from .ped import Pedigree

logger = logging.getLogger(__name__)

_CSS = """
body{font-family:system-ui,sans-serif;margin:2em auto;max-width:1100px;color:#1a1a2e;background:#fafafa}
h1{border-bottom:3px solid #2a6f97}h2{color:#2a6f97;margin-top:2em}h3{color:#444}
table{border-collapse:collapse;margin:1em 0;font-size:.9em;background:#fff}
th,td{border:1px solid #ccc;padding:4px 10px;text-align:right}
th{background:#2a6f97;color:#fff;text-align:center}
td:first-child,th:first-child{text-align:left}
tr:nth-child(even){background:#f0f4f8}
.kpi{display:inline-block;background:#fff;border:1px solid #ccc;border-radius:8px;
     padding:10px 22px;margin:6px;text-align:center}
.kpi b{display:block;font-size:1.7em;color:#2a6f97}
img{max-width:100%;background:#fff;border:1px solid #ddd;padding:4px}
.warn{color:#b23}
details{margin:.6em 0}summary{cursor:pointer;font-weight:600}
"""


def _fig_to_b64(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=110, bbox_inches="tight")
    import matplotlib.pyplot as plt

    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


def _img(fig) -> str:
    return f'<img src="data:image/png;base64,{_fig_to_b64(fig)}"/>'


def _table(headers: list[str], rows: list[list], limit: int | None = None) -> str:
    body_rows = rows[:limit] if limit else rows
    out = ["<table><tr>" + "".join(f"<th>{html.escape(str(h))}</th>" for h in headers) + "</tr>"]
    for row in body_rows:
        out.append("<tr>" + "".join(f"<td>{html.escape(str(v))}</td>" for v in row) + "</tr>")
    out.append("</table>")
    if limit and len(rows) > limit:
        out.append(f"<p><i>… {len(rows) - limit} lignes supplémentaires non affichées</i></p>")
    return "\n".join(out)


def _read_cohort_vcf(path: Path):
    """One dict per site + per-sample GT map."""
    sites, gt_by_site = [], {}
    with pysam.VariantFile(str(path)) as vcf:
        samples = list(vcf.header.samples)
        for rec in vcf:
            af = rec.info.get("AF", (0.0,))
            af = af[0] if isinstance(af, tuple) else af
            gts = {}
            for s in samples:
                gt = rec.samples[s].get("GT")
                if gt is not None and None not in gt:
                    gts[s] = (sum(gt), int(rec.samples[s].get("GQ") or 0))
            sites.append({
                "id": rec.id, "contig": rec.contig, "start": rec.start, "end": rec.stop,
                "svtype": rec.info.get("SVTYPE", "?"), "size": rec.stop - rec.start,
                "af": float(af), "ns": int(rec.info.get("NS", 0)),
                "ncarrier": int(rec.info.get("NCARRIER", 0)),
                "hwe": float(rec.info.get("HWE_P", 1.0)),
                "filter": ";".join(rec.filter.keys()) or "PASS",
            })
            gt_by_site[rec.id] = gts
    return samples, sites, gt_by_site


def _spectra_figures(sites: list[dict]):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    passing = [s for s in sites if s["filter"] == "PASS"]
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.6))
    for svtype, color in (("DEL", "#c1121f"), ("DUP", "#2a6f97")):
        sizes = [s["size"] for s in passing if s["svtype"] == svtype]
        if sizes:
            axes[0].hist(np.log10(sizes), bins=40, alpha=0.65, label=svtype, color=color)
    axes[0].set_xlabel("log10(taille pb)")
    axes[0].set_ylabel("sites PASS")
    axes[0].legend()
    axes[0].set_title("Spectre de tailles")
    afs = [s["af"] for s in passing if s["af"] > 0]
    axes[1].hist(afs, bins=40, color="#588157")
    axes[1].set_xlabel("AF cohorte")
    axes[1].set_title("Spectre de fréquences alléliques")
    fig.tight_layout()
    return fig


def _sweep_figure(sweep_rows: list[dict]):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    th = [float(r["threshold"]) for r in sweep_rows]
    mvr = [float(r["mvr"]) if r["mvr"] else np.nan for r in sweep_rows]
    dn = [int(r["n_denovo"]) for r in sweep_rows]
    fig, ax1 = plt.subplots(figsize=(6.5, 3.6))
    ax1.plot(th, mvr, "o-", color="#c1121f")
    ax1.set_xlabel("seuil GQ")
    ax1.set_ylabel("taux de violations Mendel", color="#c1121f")
    ax2 = ax1.twinx()
    ax2.plot(th, dn, "s--", color="#2a6f97")
    ax2.set_ylabel("candidats de novo", color="#2a6f97")
    ax2.set_yscale("log")
    ax1.grid(alpha=0.3)
    fig.tight_layout()
    return fig


def generate_report(
    run_dir: str | Path,
    pedigree: Pedigree | None,
    out_html: str | Path | None = None,
    cohort_name: str = "",
    make_family_workspaces: bool = True,
    min_denovo_gq: int = 20,
) -> Path:
    run_dir = Path(run_dir)
    cohort_dir = run_dir / "cohort"
    vcf_path = cohort_dir / "cohort.vcf.gz"
    if not vcf_path.exists():
        raise FileNotFoundError(f"{vcf_path} absent — lancez le pipeline cohorte d'abord")
    out_html = Path(out_html) if out_html else run_dir / "report.html"
    cohort_name = cohort_name or run_dir.name

    samples, sites, gt_by_site = _read_cohort_vcf(vcf_path)
    passing = [s for s in sites if s["filter"] == "PASS"]
    n_artifact = sum(1 for s in sites if "CohortArtifact" in s["filter"])
    # "rare" by carrier count (an AF cutoff is meaningless in small cohorts:
    # one carrier out of 44 samples is already AF 1.1%)
    rare = sorted(
        (s for s in passing if 0 < s["ncarrier"] <= 2),
        key=lambda s: (s["ncarrier"], -s["size"]),
    )

    sweep_rows, mendel_rows, summary_rows = [], [], []
    if (cohort_dir / "trio_sweep.tsv").exists():
        sweep_rows = list(csv.DictReader(open(cohort_dir / "trio_sweep.tsv"), delimiter="\t"))
    if (cohort_dir / "mendel_report.tsv").exists():
        mendel_rows = list(csv.DictReader(open(cohort_dir / "mendel_report.tsv"), delimiter="\t"))
    if (cohort_dir / "mendel_summary.tsv").exists():
        summary_rows = list(csv.DictReader(open(cohort_dir / "mendel_summary.tsv"), delimiter="\t"))

    # ---- per-sample stats
    sample_stats = {
        s: {"genotyped": 0, "carrier": 0, "gq": []} for s in samples
    }
    for gts in gt_by_site.values():
        for s, (nalt, gq) in gts.items():
            st = sample_stats[s]
            st["genotyped"] += 1
            st["gq"].append(gq)
            if nalt:
                st["carrier"] += 1

    # ---- family workspaces
    family_sections = []
    if pedigree is not None:
        fam_of = {i.sample_id: i.family_id for i in pedigree.individuals.values()}
        families = defaultdict(list)
        for s in samples:
            families[fam_of.get(s, "sans_famille")].append(s)

        for fid in sorted(families):
            members = families[fid]
            fam_dir = run_dir / "families" / fid
            fam_viol = [r for r in mendel_rows if r["family_id"] == fid and r["violation"] == "1"]
            denovo = [
                r for r in mendel_rows
                if r["family_id"] == fid and r["inheritance"] == "denovo"
                and min(int(r["child_gq"]), int(r["father_gq"]), int(r["mother_gq"])) >= min_denovo_gq
            ]
            if make_family_workspaces:
                fam_dir.mkdir(parents=True, exist_ok=True)
                # family-subset VCF: sites where >=1 member carries
                with pysam.VariantFile(str(vcf_path)) as vin:
                    vin.subset_samples(members)
                    with pysam.VariantFile(
                        str(fam_dir / f"{fid}.vcf.gz"), "wz", header=vin.header
                    ) as vout:
                        for rec in vin:
                            if any(
                                rec.samples[m].get("GT") and 1 in (rec.samples[m]["GT"] or ())
                                for m in members
                            ):
                                vout.write(rec)
                pysam.tabix_index(str(fam_dir / f"{fid}.vcf.gz"), preset="vcf", force=True)
                if fam_viol or denovo:
                    with open(fam_dir / "mendel.tsv", "w", newline="") as fh:
                        writer = csv.DictWriter(
                            fh, fieldnames=list(mendel_rows[0]), delimiter="\t"
                        )
                        writer.writeheader()
                        writer.writerows(fam_viol + denovo)

            dn_rows = [
                [r["site_id"], f'{r["contig"]}:{int(r["start"]) + 1}-{r["end"]}',
                 r["svtype"], int(r["end"]) - int(r["start"]),
                 r["child"], r["child_gq"], r["dnq"]]
                for r in denovo
            ]
            family_sections.append(f"""
<details><summary>Famille {html.escape(fid)} — {len(members)} membres,
{len(denovo)} de novo (GQ≥{min_denovo_gq}), {len(fam_viol)} violations</summary>
<p>Membres : {", ".join(map(html.escape, members))} —
workspace : <code>families/{html.escape(fid)}/</code></p>
{_table(["site", "locus", "type", "taille", "enfant", "GQ", "DNQ"], dn_rows, limit=20)
 if dn_rows else "<p>Aucun candidat de novo au seuil.</p>"}
</details>""")

    # ---- HTML assembly
    svc = Counter(s["svtype"] for s in passing)
    sample_rows = sorted(
        (
            [s, st["genotyped"], st["carrier"],
             round(float(np.median(st["gq"])), 1) if st["gq"] else 0]
            for s, st in sample_stats.items()
        ),
        key=lambda r: r[0],
    )
    rare_rows = [
        [s["id"], f'{s["contig"]}:{s["start"] + 1}-{s["end"]}', s["svtype"], s["size"],
         f'{s["af"]:.4f}', s["ncarrier"]]
        for s in rare
    ]
    sweep_table = _table(
        ["seuil GQ", "génotypés", "violations", "MVR", "transmission", "de novo"],
        [[r["threshold"], r["n_genotyped"], r["n_violation"], r["mvr"],
          r["transmission"], r["n_denovo"]] for r in sweep_rows],
    ) if sweep_rows else "<p class=warn>Pas de trios / sweep indisponible.</p>"

    page = f"""<!DOCTYPE html><html lang="fr"><head><meta charset="utf-8">
<title>Vedette — {html.escape(cohort_name)}</title><style>{_CSS}</style></head><body>
<h1>Rapport de cohorte Vedette — {html.escape(cohort_name)}</h1>
<p>Généré le {time.strftime("%Y-%m-%d %H:%M")} par vedette {__version__} —
répertoire : <code>{html.escape(str(run_dir))}</code></p>

<div>
<span class=kpi><b>{len(samples)}</b>échantillons</span>
<span class=kpi><b>{len(sites)}</b>sites</span>
<span class=kpi><b>{len(passing)}</b>PASS</span>
<span class=kpi><b>{svc.get("DEL", 0)}</b>DEL</span>
<span class=kpi><b>{svc.get("DUP", 0)}</b>DUP</span>
<span class=kpi><b>{n_artifact}</b>artefacts cohorte</span>
<span class=kpi><b>{len(rare)}</b>sites rares (≤2 porteurs)</span>
</div>

<h2>Spectres (sites PASS)</h2>
{_img(_spectra_figures(sites))}

<h2>Qualité — trios</h2>
{_img(_sweep_figure(sweep_rows)) if sweep_rows else ""}
{sweep_table}
{_table(["famille", "enfant", "génotypés", "porteur", "de novo", "violations", "MVR", "transmission"],
        [[r["family_id"], r["child"], r["n_genotyped"], r["n_child_carrier"],
          r["n_denovo"], r["n_violation"], r["violation_rate"], r["transmission_rate"]]
         for r in summary_rows]) if summary_rows else ""}

<h2>Familles</h2>
{"".join(family_sections) if family_sections else "<p>Pedigree non fourni.</p>"}

<h2>Sites rares (≤ 2 porteurs)</h2>
{_table(["site", "locus", "type", "taille", "AF", "porteurs"], rare_rows, limit=50)}

<h2>Échantillons</h2>
{_table(["échantillon", "génotypés", "porteur", "GQ médian"], sample_rows)}

</body></html>"""

    out_html.write_text(page)
    logger.info(
        "rapport: %s (%d familles, %d sites rares)",
        out_html, len(family_sections), len(rare),
    )
    return out_html
