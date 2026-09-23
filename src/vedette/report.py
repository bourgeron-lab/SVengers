"""Self-contained HTML cohort report + per-family workspaces.

``vedette report`` reads a cohort run directory (the Nextflow wrapper layout:
cohort/cohort.vcf.gz, cohort/mendel_*.tsv, cohort/trio_sweep.tsv, genotypes/)
and writes:

- ``report.html`` — one static page, figures inlined as base64 PNG: run
  summary, size/AF spectra, per-sample table, MVR/transmission sweep,
  per-family sections with de novo candidates and violations, rare sites;
- ``families/<FID>/`` — one workspace per family: the family-subset VCF
  (carriers of at least one member) and its Mendel rows; with
  ``family_events``, an interactive ``events.html`` listing every DEL/DUP
  with chromosome/type selectors, genes, NDD impact and per-member GT:GQ:SQ.
"""

from __future__ import annotations

import base64
import csv
import html
import io
import logging
import re
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pysam

from . import __version__
from .ped import Pedigree
from .qc import WAVE_WARN

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


_EVENTS_CSS = """
body{font-family:system-ui,sans-serif;font-size:13px;margin:16px;color:#1a1a2e;background:#fafafa}
h1{font-size:1.3rem}
.bar{position:sticky;top:0;background:#fafafa;padding:8px 0;border-bottom:1px solid #dde2de;
     display:flex;gap:1.2em;align-items:center;flex-wrap:wrap;z-index:2}
.bar select{font-size:13px;padding:2px 6px}
table{border-collapse:collapse;width:100%;font-variant-numeric:tabular-nums;background:#fff}
th,td{border-bottom:1px solid #dde2de;padding:3px 8px;text-align:right;white-space:nowrap}
th{position:sticky;top:44px;background:#eef1ee;font-size:11px;text-transform:uppercase;cursor:pointer;user-select:none}
th.sorted::after{content:" \\25b4"} th.sorted.desc::after{content:" \\25be"}
td:first-child,th:first-child{text-align:left}
tr.del td.t{color:#b3362b;font-weight:600} tr.dup td.t{color:#2a6f97;font-weight:600}
tr.nopass{opacity:.55} .g{text-align:left;max-width:26ch;overflow:hidden;text-overflow:ellipsis}
.c{text-align:left;font-family:monospace;font-size:11px} .muted{color:#5c6763}
"""

_EVENTS_JS = """
const rows = Array.from(document.querySelectorAll("tbody tr"));
const chrSel = document.getElementById("chr"), typeSel = document.getElementById("type"),
      passCb = document.getElementById("pass"), count = document.getElementById("count");
function apply(){
  const c = chrSel.value, t = typeSel.value, p = passCb.checked;
  let n = 0;
  for (const r of rows){
    const ok = (c === "all" || r.dataset.chr === c)
            && (t === "all" || r.dataset.type === t)
            && (!p || r.dataset.pass === "1");
    r.hidden = !ok;
    if (ok) n++;
  }
  count.textContent = n.toLocaleString("fr-FR") + " / " + rows.length.toLocaleString("fr-FR");
}
chrSel.onchange = typeSel.onchange = passCb.onchange = apply;
apply();

// column sorting: click a header to sort, click again to reverse
const tbody = document.querySelector("tbody");
const NUM = /^-?[0-9]+([.][0-9]+)?$/;
function cellKey(r, i){
  const t = r.children[i].textContent.trim().replace(/[\\s,\\u202f]/g, "");
  return NUM.test(t) ? parseFloat(t) : r.children[i].textContent.trim().toLowerCase();
}
document.querySelectorAll("th").forEach((th, i) => {
  th.onclick = () => {
    const desc = th.classList.contains("sorted") && !th.classList.contains("desc");
    document.querySelectorAll("th").forEach(h => h.classList.remove("sorted", "desc"));
    th.classList.add("sorted");
    if (desc) th.classList.add("desc");
    const sorted = rows.slice().sort((a, b) => {
      const ka = cellKey(a, i), kb = cellKey(b, i);
      const c = (typeof ka === "number" && typeof kb === "number")
        ? ka - kb : String(ka).localeCompare(String(kb), "fr");
      return desc ? -c : c;
    });
    tbody.append(...sorted);
  };
});
"""

_INHERIT_FR = {
    "denovo": "de novo",
    "denovo_lowconf": "de novo (GQ faible)",
    "maternal": "maternel",
    "paternal": "paternel",
    "biparental": "biparental",
    "none": "aucun parent",
}


def _event_row(rec, members: list[str], genes_fn, impact_fn, inherit: str = "") -> dict:
    """Event dict for `_write_events_page` from a family-subset VCF record."""
    carriers, gts = [], []
    for m in members:
        fmt = rec.samples[m]
        gt = fmt.get("GT")
        gt_s = "./." if gt is None or gt[0] is None else "/".join(str(a) for a in gt)
        gq = fmt.get("GQ")
        sq = fmt.get("SQ")
        gts.append((gt_s, gq if gq is not None else "", f"{float(sq):.1f}" if sq is not None else ""))
        if gt_s not in ("0/0", "./."):
            carriers.append(m)
    return {
        "contig": rec.contig,
        "pos": rec.pos,
        "end": rec.stop,
        "svtype": rec.info.get("SVTYPE", "?"),
        "size": rec.stop - rec.pos,
        "filter": ";".join(rec.filter.keys()) or "PASS",
        "rescue": rec.info.get("RESCUE", "") or "",
        "genes": genes_fn(rec.contig, rec.start, rec.stop),
        "impact": impact_fn(rec.contig, rec.start, rec.stop)[1],
        "inherit": inherit,
        "carriers": carriers,
        "gt": gts,
    }


def _write_events_page(
    out_file: Path, title: str, members: list[str], events: list[dict]
) -> None:
    """One interactive HTML page listing every DEL/DUP of a family workspace:
    chromosome and type selectors plus a PASS-only toggle filter the table
    client-side (no server, opens from disk). Genes/impact and per-member
    GT:GQ:SQ on each row.

    ``events`` rows: contig, pos (1-based), end, svtype, size, filter, rescue,
    genes, impact, carriers (list), gt (list of (GT, GQ, SQ) per member).
    """
    by_chr: dict[str, int] = Counter(e["contig"] for e in events)
    contigs = sorted(by_chr, key=lambda c: (len(c), c))
    n_del = sum(1 for e in events if e["svtype"] == "DEL")
    n_dup = sum(1 for e in events if e["svtype"] == "DUP")

    chr_opts = '<option value="all">tous</option>' + "".join(
        f'<option value="{html.escape(c)}">{html.escape(c)} ({by_chr[c]})</option>'
        for c in contigs
    )
    head = (
        "<tr><th>chr</th><th>pos</th><th>end</th><th>type</th><th>taille</th><th>FILTER</th>"
        "<th>rescue</th><th>gènes</th><th>impact</th><th>héritage</th><th>porteurs</th>"
        + "".join(f"<th>{html.escape(m[-3:])}</th>" for m in members) + "</tr>"
    )
    rows = []
    for e in sorted(events, key=lambda e: (contigs.index(e["contig"]), e["pos"])):
        is_pass = e["filter"] == "PASS"
        cls = e["svtype"].lower() + ("" if is_pass else " nopass")
        gt_cells = "".join(
            f'<td class="c">{html.escape(gt)}<span class="muted">:{gq}:{sq}</span></td>'
            for gt, gq, sq in e["gt"]
        )
        rows.append(
            f'<tr class="{cls}" data-chr="{html.escape(e["contig"])}" '
            f'data-type="{e["svtype"]}" data-pass="{1 if is_pass else 0}">'
            f'<td>{html.escape(e["contig"])}</td><td>{e["pos"]:,}</td><td>{e["end"]:,}</td>'
            f'<td class="t">{e["svtype"]}</td><td>{e["size"]:,}</td>'
            f'<td>{html.escape(e["filter"])}</td><td>{html.escape(e["rescue"])}</td>'
            f'<td class="g">{html.escape(e["genes"] or "—")}</td>'
            f'<td class="g">{html.escape(e["impact"])}</td>'
            f'<td class="g">{html.escape(e.get("inherit") or "—")}</td>'
            f'<td class="c">{",".join(m[-3:] for m in e["carriers"]) or "—"}</td>{gt_cells}</tr>'
        )
    out_file.write_text(
        f"<!doctype html><title>{html.escape(title)} — événements</title>"
        f"<style>{_EVENTS_CSS}</style><h1>{html.escape(title)} — "
        f"{n_del} DEL + {n_dup} DUP</h1>"
        f'<div class="bar">'
        f'<label>chromosome <select id="chr">{chr_opts}</select></label>'
        f'<label>type <select id="type"><option value="all">tous</option>'
        f'<option>DEL</option><option>DUP</option></select></label>'
        f'<label><input type="checkbox" id="pass" checked> PASS uniquement</label>'
        f'<span class="muted"><b id="count"></b> lignes</span></div>'
        f'<p class="muted">GT<span class="muted">:GQ:SQ</span> par membre (suffixe du '
        f'barcode ; ordre : {html.escape(", ".join(members))}). Lignes grisées = non-PASS.</p>'
        # the header must live in <thead>: a loose <tr> lands in an implicit
        # tbody and the filter JS (tbody tr) hides/sorts it like a data row
        f'<table><thead>{head}</thead><tbody>{"".join(rows)}</tbody></table>'
        f"<script>{_EVENTS_JS}</script>"
    )


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


class _GeneIndex:
    """Overlap lookup of genes/exons in a tabix-indexed Ensembl GTF/GFF."""

    def __init__(self, gff_path: str | Path, ndd_tsv: str | Path | None = None):
        self.tbx = pysam.TabixFile(str(gff_path))
        self._contigs = set(self.tbx.contigs)
        # GeneTrek NDD tiers, keyed by ENSG and by symbol
        self.ndd: dict[str, int] = {}   # 2 = high-confidence, 1 = candidate
        if ndd_tsv:
            with open(ndd_tsv) as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    tier = 0
                    if row.get("hc_ndd_v7", "").strip().lower() == "true":
                        tier = 2
                    elif row.get("candidate_ndd_v5", "").strip().lower() == "true":
                        tier = 1
                    if tier:
                        for key in (row.get("ensg", ""), row.get("hgnc_symbol", "")):
                            if key:
                                self.ndd[key] = max(self.ndd.get(key, 0), tier)

    def _hits(self, contig: str, start: int, end: int):
        """[(name, ensg, coding, exon_overlap)] of genes overlapping the site."""
        if contig not in self._contigs:
            contig = contig.removeprefix("chr")
            if contig not in self._contigs:
                return []
        genes: dict[str, list] = {}
        try:
            rows = self.tbx.fetch(contig, max(start, 0), end)
        except (ValueError, OSError):
            return []
        for row in rows:
            f = row.split("\t")
            if len(f) < 9 or f[2] not in ("gene", "exon"):
                continue
            gid = re.search(r'gene_id="?([^";.]+)', f[8])
            name = re.search(r'gene_name="?([^";]+)"?', f[8])
            if not gid:
                continue
            g = genes.setdefault(
                gid.group(1),
                [name.group(1) if name else gid.group(1), gid.group(1),
                 "protein_coding" in f[8], False],
            )
            if f[2] == "exon":
                g[3] = True
            elif "protein_coding" in f[8]:
                g[2] = True
        return list(genes.values())

    def genes(self, contig: str, start: int, end: int, max_names: int = 6) -> str:
        hits = self._hits(contig, start, end)
        coding = [h[0] for h in hits if h[2]]
        names = coding or [h[0] for h in hits]
        out = ", ".join(names[:max_names])
        if len(names) > max_names:
            out += f" (+{len(names) - max_names})"
        return out

    def impact(self, contig: str, start: int, end: int) -> tuple[int, str]:
        """(score, label): 0 intergénique, 1 intronique, 2 exonique,
        +1 gène NDD candidat / +2 NDD haute confiance (GeneTrek)."""
        hits = [h for h in self._hits(contig, start, end) if h[2]]  # coding only
        if not hits:
            return 0, "0 — intergénique"
        exonic = any(h[3] for h in hits)
        score = 2 if exonic else 1
        label = "exonique" if exonic else "intronique"
        ndd_tier, ndd_genes = 0, []
        for name, ensg, _coding, _exon in hits:
            tier = max(self.ndd.get(ensg, 0), self.ndd.get(name, 0))
            if tier:
                ndd_genes.append(name)
                ndd_tier = max(ndd_tier, tier)
        if ndd_tier:
            score += ndd_tier
            label += " · NDD-" + ("HC" if ndd_tier == 2 else "cand")
            label += f" ({', '.join(sorted(set(ndd_genes))[:3])})"
        return score, f"{score} — {label}"


def denovo_outliers(
    denovo_by_family: dict[str, list[dict]], threshold: int
) -> dict[str, int]:
    """Children carrying more GQ-filtered de novo candidates than `threshold`.

    A genuine genome carries <1 de novo CNV on average; a child with dozens is
    a technical outlier (waviness under the QC bar, contamination, swap) whose
    candidates would drown the family tables — measured on dataset_paper: two
    such children held 39 of the 81 top-impact candidates.
    """
    counts: dict[str, int] = defaultdict(int)
    for rows in denovo_by_family.values():
        for r in rows:
            counts[r["child"]] += 1
    return {c: n for c, n in counts.items() if n > threshold}


def _read_cohort_vcf(path: Path):
    """One dict per site + compact genotype matrices.

    Dict-of-dicts genotypes cost ~80 bytes/genotype and OOM-killed the report
    on 238k sites x 2670 samples (635M genotypes); nalt/gq live in two numpy
    matrices (int8+int16 = 3 bytes/genotype, ~1.9 GB at that scale) and each
    site carries its row index.
    """
    sites = []
    nalt_rows, gq_rows = [], []
    with pysam.VariantFile(str(path)) as vcf:
        samples = list(vcf.header.samples)
        n = len(samples)
        for rec in vcf:
            af = rec.info.get("AF", (0.0,))
            af = af[0] if isinstance(af, tuple) else af
            nalt = np.full(n, -1, np.int8)   # -1 = missing
            gq = np.zeros(n, np.int16)
            for i, s in enumerate(samples):
                fmt = rec.samples[s]
                gt = fmt.get("GT")
                if gt is not None and None not in gt:
                    nalt[i] = sum(gt)
                    gq[i] = int(fmt.get("GQ") or 0)
            sites.append({
                "id": rec.id, "contig": rec.contig, "start": rec.start, "end": rec.stop,
                "svtype": rec.info.get("SVTYPE", "?"), "size": rec.stop - rec.start,
                "af": float(af), "ns": int(rec.info.get("NS", 0)),
                "ncarrier": int(rec.info.get("NCARRIER", 0)),
                "hwe": float(rec.info.get("HWE_P", 1.0)),
                "filter": ";".join(rec.filter.keys()) or "PASS",
                "idx": len(sites),
            })
            nalt_rows.append(nalt)
            gq_rows.append(gq)
    nalt_m = np.vstack(nalt_rows) if nalt_rows else np.zeros((0, len(samples)), np.int8)
    gq_m = np.vstack(gq_rows) if gq_rows else np.zeros((0, len(samples)), np.int16)
    return samples, sites, nalt_m, gq_m


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
    genes_gff: str | Path | None = None,
    ndd_tsv: str | Path | None = None,
    max_denovo_per_child: int = 20,
    family_events: bool = False,
) -> Path:
    run_dir = Path(run_dir)
    cohort_dir = run_dir / "cohort"
    vcf_path = cohort_dir / "cohort.vcf.gz"
    if not vcf_path.exists():
        raise FileNotFoundError(f"{vcf_path} absent — lancez le pipeline cohorte d'abord")
    out_html = Path(out_html) if out_html else run_dir / "report.html"
    cohort_name = cohort_name or run_dir.name

    gene_idx = _GeneIndex(genes_gff, ndd_tsv) if genes_gff else None

    def _genes(contig: str, start: int, end: int) -> str:
        return gene_idx.genes(contig, start, end) if gene_idx else ""

    def _impact(contig: str, start: int, end: int) -> tuple[int, str]:
        return gene_idx.impact(contig, start, end) if gene_idx else (0, "")

    samples, sites, nalt_m, gq_m = _read_cohort_vcf(vcf_path)
    passing = [s for s in sites if s["filter"] == "PASS"]
    n_artifact = sum(1 for s in sites if "CohortArtifact" in s["filter"])
    # "rare" by carrier count (an AF cutoff is meaningless in small cohorts:
    # one carrier out of 44 samples is already AF 1.1%)
    rare = sorted(
        (s for s in passing if 0 < s["ncarrier"] <= 2),
        key=lambda s: (s["ncarrier"], -s["size"]),
    )

    sweep_rows, summary_rows = [], []
    if (cohort_dir / "trio_sweep.tsv").exists():
        sweep_rows = list(csv.DictReader(open(cohort_dir / "trio_sweep.tsv"), delimiter="\t"))
    if (cohort_dir / "mendel_summary.tsv").exists():
        summary_rows = list(csv.DictReader(open(cohort_dir / "mendel_summary.tsv"), delimiter="\t"))

    # mendel_report.tsv scales with trios x sites (11 GB on dataset_paper):
    # stream it, keeping only GQ-filtered de novo rows and per-family
    # violation rows (the latter only when family workspaces are written)
    denovo_by_family: dict[str, list[dict]] = defaultdict(list)
    viol_by_family: dict[str, list[dict]] = defaultdict(list)
    viol_counts: dict[str, int] = defaultdict(int)
    # (family, site) -> [(child, inheritance)] for the events pages — family
    # runs only (scales with trios x carried sites, unbounded on big cohorts)
    inherit_by_site: dict[tuple[str, str], list[tuple[str, str]]] = defaultdict(list)
    collect_inherit = family_events and make_family_workspaces
    mendel_fields: list[str] = []
    if (cohort_dir / "mendel_report.tsv").exists():
        # de novo candidacy demands a PASS site: the rescue reservoir and
        # cohort-artifact sites otherwise flood the tables (dataset_paper:
        # ~490 "de novo"/child on all sites vs a median of 2 on PASS)
        pass_ids = {s["id"] for s in sites if s["filter"] == "PASS"}
        with open(cohort_dir / "mendel_report.tsv") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            mendel_fields = list(reader.fieldnames or [])
            for r in reader:
                if collect_inherit and r.get("inheritance"):
                    inherit_by_site[(r["family_id"], r["site_id"])].append(
                        (r["child"], r["inheritance"])
                    )
                if r["violation"] == "1":
                    viol_counts[r["family_id"]] += 1
                    if make_family_workspaces:
                        viol_by_family[r["family_id"]].append(r)
                if (
                    r["inheritance"] == "denovo"
                    and r["site_id"] in pass_ids
                    and min(int(r["child_gq"]), int(r["father_gq"]), int(r["mother_gq"]))
                    >= min_denovo_gq
                ):
                    denovo_by_family[r["family_id"]].append(r)

    # ---- de novo per-child guard: technical outliers leave the family tables
    dn_outliers = denovo_outliers(denovo_by_family, max_denovo_per_child)
    if dn_outliers:
        logger.warning(
            "garde-fou de novo: %d enfant(s) au-dessus de %d candidats exclu(s): %s",
            len(dn_outliers), max_denovo_per_child,
            ", ".join(f"{c} ({n})" for c, n in sorted(dn_outliers.items())),
        )

    # ---- per-sample stats (vectorized over the genotype matrices)
    genotyped_per_sample = (nalt_m >= 0).sum(axis=0)
    carrier_per_sample = (nalt_m > 0).sum(axis=0)
    masked_gq = np.ma.masked_array(gq_m, mask=(nalt_m < 0))
    median_gq_per_sample = np.ma.median(masked_gq, axis=0).filled(0)
    sample_stats = {
        s: {
            "genotyped": int(genotyped_per_sample[i]),
            "carrier": int(carrier_per_sample[i]),
            "median_gq": float(median_gq_per_sample[i]),
        }
        for i, s in enumerate(samples)
    }

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
            fam_viol = viol_by_family.get(fid, [])
            n_viol = viol_counts.get(fid, 0)
            denovo_all = denovo_by_family.get(fid, [])
            denovo = [r for r in denovo_all if r["child"] not in dn_outliers]
            fam_outliers = sorted(
                {r["child"] for r in denovo_all} & set(dn_outliers)
            )
            if make_family_workspaces:
                fam_dir.mkdir(parents=True, exist_ok=True)
                fam_events: list[dict] = []
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
                                if family_events:
                                    inh = inherit_by_site.get((fid, rec.id or ""), [])
                                    inh_s = " · ".join(
                                        f"{_INHERIT_FR.get(v, v)} ({c[-3:]})" for c, v in inh
                                    ) if len(inh) > 1 else (
                                        _INHERIT_FR.get(inh[0][1], inh[0][1]) if inh else ""
                                    )
                                    fam_events.append(
                                        _event_row(rec, members, _genes, _impact, inherit=inh_s)
                                    )
                pysam.tabix_index(str(fam_dir / f"{fid}.vcf.gz"), preset="vcf", force=True)
                if family_events:
                    _write_events_page(
                        fam_dir / "events.html", f"Famille {fid}", members, fam_events
                    )
                if (fam_viol or denovo) and mendel_fields:
                    with open(fam_dir / "mendel.tsv", "w", newline="") as fh:
                        writer = csv.DictWriter(fh, fieldnames=mendel_fields, delimiter="\t")
                        writer.writeheader()
                        writer.writerows(fam_viol + denovo)

            dn_rows = [
                [r["site_id"], f'{r["contig"]}:{int(r["start"]) + 1}-{r["end"]}',
                 r["svtype"], int(r["end"]) - int(r["start"]),
                 _genes(r["contig"], int(r["start"]), int(r["end"])),
                 _impact(r["contig"], int(r["start"]), int(r["end"]))[1],
                 r["child"], r["child_gq"], r["dnq"]]
                for r in denovo
            ]
            dn_rows.sort(key=lambda row: -int(str(row[5]).split(" ")[0] or 0))
            family_sections.append(f"""
<details><summary>Famille {html.escape(fid)} — {len(members)} membres,
{len(denovo)} de novo (PASS, GQ≥{min_denovo_gq}), {n_viol} violations</summary>
<p>Membres : {", ".join(map(html.escape, members))} —
workspace : <code>families/{html.escape(fid)}/</code>{f' · <a href="families/{html.escape(fid)}/events.html">tous les événements (filtres chr/type)</a>' if family_events and make_family_workspaces else ""}</p>
{f'<p class=warn>⚠ {", ".join(map(html.escape, fam_outliers))} : plus de {max_denovo_per_child} candidats de novo — profil technique atypique, candidats exclus du tableau.</p>' if fam_outliers else ""}
{_table(["site", "locus", "type", "taille", "gènes", "impact", "enfant", "GQ", "DNQ"], dn_rows, limit=20)
 if dn_rows else "<p>Aucun candidat de novo au seuil.</p>"}
</details>""")

    # ---- per-sample wave QC (evidence-store meta, absent on old stores)
    wave_scores: dict[str, float | None] = {}
    for s in samples:
        store_path = run_dir / "samples" / s / f"{s}.evidence.zarr"
        wave_scores[s] = None
        if store_path.exists():
            try:
                from .extract.store import EvidenceStore

                wave_scores[s] = EvidenceStore(store_path).meta.get("wave_score")
            except Exception:  # report survives a broken/foreign store
                pass
    n_wavy = sum(1 for v in wave_scores.values() if v is not None and v >= WAVE_WARN)

    # ---- HTML assembly
    svc = Counter(s["svtype"] for s in passing)
    sample_rows = sorted(
        (
            [s, st["genotyped"], st["carrier"],
             round(st["median_gq"], 1),
             ("—" if wave_scores[s] is None else
              f"{wave_scores[s]:.3f} ⚠" if wave_scores[s] >= WAVE_WARN
              else f"{wave_scores[s]:.3f}")]
            for s, st in sample_stats.items()
        ),
        key=lambda r: r[0],
    )
    def _carrier_ids(site: dict, cap: int = 6) -> str:
        ids = sorted(samples[i] for i in np.flatnonzero(nalt_m[site["idx"]] > 0))
        out = ", ".join(ids[:cap])
        if len(ids) > cap:
            out += f" (+{len(ids) - cap})"
        return out

    def _rare_rows(svtype: str) -> list[list]:
        rows = []
        for s in rare:
            if s["svtype"] != svtype:
                continue
            score, label = _impact(s["contig"], s["start"], s["end"])
            rows.append((score, [
                s["id"], f'{s["contig"]}:{s["start"] + 1}-{s["end"]}', s["size"],
                _genes(s["contig"], s["start"], s["end"]), label,
                f'{s["af"]:.4f}', _carrier_ids(s),
            ]))
        # triage: NDD/exonique d'abord, puis porteurs rares, puis taille
        rows.sort(key=lambda t: (-t[0], -t[1][2]))
        return [r for _score, r in rows]
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
<span class=kpi><b>{n_wavy}</b>échantillons wavy (score ≥ {WAVE_WARN})</span>
<span class=kpi><b>{len(dn_outliers)}</b>enfants outliers de novo (&gt;{max_denovo_per_child})</span>
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
<p>Impact : 0 intergénique, 1 intronique, 2 exonique ; +1 gène NDD candidat,
+2 gène NDD haute confiance (GeneTrek). Tri par impact puis taille.</p>
<h3>DEL ({sum(1 for s in rare if s["svtype"] == "DEL")})</h3>
{_table(["site", "locus", "taille", "gènes", "impact", "AF", "porteurs"], _rare_rows("DEL"), limit=40)}
<h3>DUP ({sum(1 for s in rare if s["svtype"] == "DUP")})</h3>
{_table(["site", "locus", "taille", "gènes", "impact", "AF", "porteurs"], _rare_rows("DUP"), limit=40)}

<h2>Échantillons</h2>
<p>Le score de wave mesure les oscillations mégabase de couverture survivant à la
correction GC (artefact de librairie) ; au-delà de {WAVE_WARN}, les calls depth-only
de l'échantillon sont à considérer avec prudence.</p>
{_table(["échantillon", "génotypés", "porteur", "GQ médian", "wave"], sample_rows)}

</body></html>"""

    out_html.write_text(page)
    logger.info(
        "rapport: %s (%d familles, %d sites rares)",
        out_html, len(family_sections), len(rare),
    )
    return out_html
