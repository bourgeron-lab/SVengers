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
td.val{text-align:left;white-space:nowrap} td.val label{margin-right:.4em}
a.igv{font-size:11px;padding:1px 6px;border:1px solid #9ab;border-radius:4px;text-decoration:none;color:#2a6f97}
.bar button{font-size:12px;padding:2px 8px;cursor:pointer}
"""

_EVENTS_JS = """
const rows = Array.from(document.querySelectorAll("tbody tr"));
const chrSel = document.getElementById("chr"), typeSel = document.getElementById("type"),
      valSel = document.getElementById("val"), ownSel = document.getElementById("own"),
      supSel = document.getElementById("sup"),
      passCb = document.getElementById("pass"),
      count = document.getElementById("count");
function apply(){
  const c = chrSel.value, t = typeSel.value, v = valSel.value, p = passCb.checked,
        o = ownSel.value, sp = supSel.value;
  let n = 0;
  for (const r of rows){
    const ok = (c === "all" || r.dataset.chr === c)
            && (t === "all" || r.dataset.type === t)
            && (v === "all" || (r.dataset.val || "none") === v)
            && (o === "all" || (o === "any" ? r.dataset.own !== "aucun" : r.dataset.own === o))
            && (sp === "all" || r.dataset.support === sp)
            && (!p || r.dataset.pass === "1");
    r.hidden = !ok;
    if (ok) n++;
  }
  count.textContent = n.toLocaleString("fr-FR") + " / " + rows.length.toLocaleString("fr-FR");
}
chrSel.onchange = typeSel.onchange = valSel.onchange = ownSel.onchange = supSel.onchange = passCb.onchange = apply;

// manual validation (IGV check): three exclusive radio buttons per row
// (vrai/faux/incertain — a forced binary choice pollutes labels on genuinely
// ambiguous calls: borderline coverage, fuzzy breakpoint, repeat region),
// remembered in this browser (localStorage) and exportable as TSV
const PAGE = "vedette-val:" + document.title + ":";
const VAL_LABEL = {true: "vrai", false: "faux", unsure: "incertain"};
function store(){ try { return window.localStorage; } catch (e) { return null; } }
function setVal(r, v){
  r.dataset.val = v || "";
  for (const opt of ["true", "false", "unsure"]) r.querySelector(".v-" + opt).checked = v === opt;
  const s = store();
  if (s) { try { v ? s.setItem(PAGE + r.dataset.site, v) : s.removeItem(PAGE + r.dataset.site); } catch (e) {} }
}
for (const r of rows){
  const s = store(); let v = "";
  if (s) { try { v = s.getItem(PAGE + r.dataset.site) || ""; } catch (e) {} }
  if (v) setVal(r, v);
  for (const opt of ["true", "false", "unsure"])
    r.querySelector(".v-" + opt).onchange = e => { setVal(r, e.target.checked ? opt : ""); apply(); };
  // send the locus (with 20% flanks) to a running IGV desktop (batch port).
  // A background fetch from a file:// page to localhost is blocked by some
  // browsers (private-network access): try it, and fall back to a real
  // navigation in a dedicated tab, which is never blocked (IGV answers "OK").
  r.querySelector(".igv").onclick = e => {
    e.preventDefault();
    const url = e.currentTarget.href;
    fetch(url, {mode: "no-cors"}).catch(() => window.open(url, "igv"));
  };
}
document.getElementById("export").onclick = () => {
  // resolve columns by header text: hardcoded cell indices break silently
  // every time a column is inserted
  const COL = {};
  document.querySelectorAll("thead th").forEach((th, i) => { COL[th.textContent.trim()] = i; });
  const EXPORT = [["filter", "FILTER"], ["sqmax", "SQmax"], ["carrier_frac", "% port."],
                  ["ref_carrier_frac", "fréq. réf."], ["lcr", "LCR"], ["segdup", "segdup"],
                  ["heritage", "héritage"], ["porteurs", "porteurs"], ["appel_propre", "appel propre"],
                  ["support", "support"]];
  const head = ["site_id", "chr", "pos", "end", "type", "taille",
                ...EXPORT.map(c => c[0]), "validation"];
  const lines = [head.join("\\t")];
  for (const r of rows){
    if (!r.dataset.val) continue;
    const td = r.children;
    lines.push([r.dataset.site, r.dataset.chr, r.dataset.pos, r.dataset.end, r.dataset.type,
                r.dataset.size,
                ...EXPORT.map(c => (td[COL[c[1]]] || {textContent: ""}).textContent),
                VAL_LABEL[r.dataset.val] || r.dataset.val].join("\\t"));
  }
  const a = document.createElement("a");
  a.href = URL.createObjectURL(new Blob([lines.join("\\n") + "\\n"], {type: "text/tab-separated-values"}));
  a.download = document.title.replace(/[^\\w.-]+/g, "_") + ".validation.tsv";
  a.click();
};
document.getElementById("import").onchange = e => {
  const f = e.target.files[0]; if (!f) return;
  f.text().then(txt => {
    // site ids are positional (vedette_site_N) and renumber when the cohort is
    // re-merged: match by id, then by exact locus, then by RO >= 0.9 same type
    const bySite = {}, byLocus = {};
    for (const r of rows){
      bySite[r.dataset.site] = r;
      byLocus[r.dataset.chr + ":" + r.dataset.pos + "-" + r.dataset.end + ":" + r.dataset.type] = r;
    }
    const L = txt.split(/\\r?\\n/); const h = L[0].split("\\t");
    const iS = h.indexOf("site_id"), iV = h.indexOf("validation"),
          iC = h.indexOf("chr"), iP = h.indexOf("pos"), iE = h.indexOf("end"), iT = h.indexOf("type");
    if (iS < 0 || iV < 0) { alert("TSV attendu avec colonnes site_id et validation"); return; }
    function findRow(c){
      let r = bySite[c[iS]];
      if (r || iC < 0 || iP < 0 || iE < 0 || iT < 0) return r;
      r = byLocus[c[iC] + ":" + c[iP] + "-" + c[iE] + ":" + c[iT]];
      if (r) return r;
      const p = parseInt(c[iP]), e = parseInt(c[iE]); let best = null, bestRo = 0;
      for (const q of rows){
        if (q.dataset.chr !== c[iC] || q.dataset.type !== c[iT]) continue;
        const qp = parseInt(q.dataset.pos), qe = parseInt(q.dataset.end);
        const ov = Math.min(e, qe) - Math.max(p, qp);
        if (ov <= 0) continue;
        const ro = ov / Math.max(e - p, qe - qp);
        if (ro >= 0.9 && ro > bestRo) { best = q; bestRo = ro; }
      }
      return best;
    }
    let n = 0;
    for (const line of L.slice(1)){
      const c = line.split("\\t"); const r = findRow(c); if (!r) continue;
      const v = (c[iV] || "").trim().toLowerCase();
      const norm = v === "vrai" ? "true" : v === "faux" ? "false"
        : v === "incertain" ? "unsure" : ["true", "false", "unsure"].includes(v) ? v : "";
      setVal(r, norm);
      n++;
    }
    apply(); alert(n + " validations importées");
  });
};
apply();

// column sorting: click a header to sort, click again to reverse
const tbody = document.querySelector("tbody");
const ths = Array.from(document.querySelectorAll("th"));
const NUM = /^-?[0-9]+([.][0-9]+)?$/;
function cellKey(r, i){
  const k = ths[i].dataset.key;
  if (k) return r.dataset[k] || "";
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


_NALT_GT = {-1: "./.", 0: "0/0", 1: "0/1", 2: "1/1"}


class _FamilySpill:
    """Buffered per-family append-only TSV spill files.

    The events pass is site-major (one sweep over the genotype matrices) but
    the pages are per family; holding every family's rows in memory scales
    with families x carried sites (~15M rows on dataset_paper), so rows are
    spilled to one file per family and re-read family by family at assembly.
    """

    def __init__(self, root: Path, flush_at: int = 500_000):
        self.root = root
        self.root.mkdir(parents=True, exist_ok=True)
        self.buf: dict[str, list[str]] = defaultdict(list)
        self.n = 0
        self.flush_at = flush_at

    def add(self, fid: str, line: str) -> None:
        self.buf[fid].append(line)
        self.n += 1
        if self.n >= self.flush_at:
            self.flush()

    def flush(self) -> None:
        for fid, lines in self.buf.items():
            with open(self.root / f"{fid}.tsv", "a") as fh:
                fh.write("\n".join(lines) + "\n")
        self.buf.clear()
        self.n = 0

    def read(self, fid: str) -> list[list[str]]:
        path = self.root / f"{fid}.tsv"
        if fid in self.buf or not path.exists():
            self.flush()
        if not path.exists():
            return []
        return [line.split("\t") for line in path.read_text().splitlines()]

    def cleanup(self) -> None:
        import shutil

        shutil.rmtree(self.root, ignore_errors=True)


_EVENTS_MAX_MEMBERS = 20  # no events page for catch-all pseudo-families


def own_call_index(run_dir: Path, members: list[str]) -> dict:
    """{(contig, svtype): [(start, stop, is_pass)]} of the family's OWN calls.

    A cohort site records a breakpoint cluster; re-genotyping then assigns a
    genotype to every sample, whether or not that sample ever produced a call
    there. On PMS, 80% of the PASS events a family "carries" have no call of
    their own behind them — and blind IGV review scored those 0/33, against
    3/9 for the ones backed by an own PASS call. The distinction has to be
    visible, so it is computed here from the published per-sample VCFs.
    """
    idx: dict[tuple[str, str], list[tuple[int, int, bool]]] = defaultdict(list)
    for m in members:
        vcf_path = run_dir / "samples" / m / f"{m}.filtered.vcf.gz"
        if not vcf_path.exists():
            continue
        try:
            with pysam.VariantFile(str(vcf_path)) as vf:
                for rec in vf:
                    svtype = rec.info.get("SVTYPE")
                    if svtype:
                        idx[(rec.contig, svtype)].append(
                            (rec.pos, rec.stop, list(rec.filter.keys()) == ["PASS"])
                        )
        except (OSError, ValueError):  # a report must survive a missing VCF
            continue
    for k in idx:
        idx[k].sort()
    return idx


def own_call_status(idx: dict, contig: str, pos: int, end: int, svtype: str,
                    min_ro: float = 0.5) -> str:
    """"PASS", "filtré" or "" for the family's own call behind a cohort site."""
    best = ""
    for op, oe, is_pass in idx.get((contig, svtype), ()):
        if op > end:
            break
        ov = min(end, oe) - max(pos, op)
        if ov > 0 and ov / max(end - pos, oe - op) >= min_ro:
            if is_pass:
                return "PASS"
            best = "filtré"
    return best


# A fragmented event leaves gaps where the HMM called CN=0 islands or where a
# segment failed a filter. Tolerating a gap proportional to the larger of the
# two neighbours closes them without chaining unrelated calls: measured over
# 25 PMS families it merges 0.9% of rows and the widest span stays 8.7 Mb
# (none above 15 Mb). The absolute cap bounds the worst case.
GAP_ABS = 2_000
GAP_FRAC = 0.02
GAP_MAX = 250_000


def collapse_events(events: list[dict], genes_fn=None, impact_fn=None,
                    segdup_iv=None, lcr_iv=None) -> list[dict]:
    """Merge overlapping same-type events carried by exactly the same members.

    A cohort site is one breakpoint cluster, and in a cohort ascertained for a
    locus (Phelan-McDermid: terminal 22q13 deletions of every size) each
    patient's deletion spans many of them, so one biological event shows up as
    a stack of nested rows. Grouping on the exact carrier set keeps genuinely
    different events apart: a de novo deletion in the child stays separate
    from an inherited one the parents also carry, even where they overlap.

    Neighbours separated by a small gap are merged too (see GAP_*): a terminal
    deletion comes back as several intervals with CN=0 islands or filtered
    segments between them — FAM-01's 8.9 Mb 22q13 deletion arrived as
    four pieces split by gaps of 3, 15 and 55 kb.
    """
    groups: dict[tuple, list[dict]] = defaultdict(list)
    for e in events:
        groups[(e["contig"], e["svtype"], tuple(sorted(e["carriers"])))].append(e)
    out: list[dict] = []
    for (contig, svtype, _carriers), evs in groups.items():
        evs.sort(key=lambda e: (e["pos"], e["end"]))
        cur: dict | None = None
        for e in evs:
            gap_tol = min(
                max(GAP_ABS, int(GAP_FRAC * max(
                    cur["end"] - cur["pos"], e["end"] - e["pos"]))), GAP_MAX,
            ) if cur is not None else 0
            if cur is not None and e["pos"] <= cur["end"] + gap_tol:
                cur["end"] = max(cur["end"], e["end"])
                cur["_ids"].append(e["id"])
                for k in ("sqmax", "cfrac", "refcf"):
                    a, b = cur.get(k) or "", e.get(k) or ""
                    if b and (not a or float(b) > float(a)):
                        cur[k] = b
                if e.get("support") == "fort":
                    cur["support"] = "fort"
                if e.get("own") == "PASS" or (e.get("own") and not cur.get("own")):
                    cur["own"] = e["own"]
                for k in ("rescue", "inherit"):
                    if e.get(k) and e[k] not in (cur.get(k) or ""):
                        cur[k] = " · ".join(x for x in (cur.get(k), e[k]) if x)
                continue
            if cur is not None:
                out.append(cur)
            cur = dict(e, _ids=[e["id"]])
        if cur is not None:
            out.append(cur)
    # re-annotate merged spans: genes, impact and masks describe the union
    for e in out:
        e["n_sites"] = len(e["_ids"])
        e["size"] = e["end"] - e["pos"]
        if e["n_sites"] > 1:
            start0 = e["pos"] - 1
            if genes_fn is not None:
                e["genes"] = genes_fn(e["contig"], start0, e["end"])
            if impact_fn is not None:
                e["impact"] = impact_fn(e["contig"], start0, e["end"])[1]
            if lcr_iv is not None:
                e["lcr"] = f"{lcr_iv.overlap_fraction(e['contig'], start0, e['end']):.2f}"
            if segdup_iv is not None:
                e["segdup"] = f"{segdup_iv.overlap_fraction(e['contig'], start0, e['end']):.2f}"
    out.sort(key=lambda e: (e["contig"], e["pos"]))
    return out


def _write_all_family_events(
    run_dir: Path,
    families: dict[str, list[str]],
    samples: list[str],
    sites: list[dict],
    nalt_m,
    gq_m,
    cn_m,
    genes_fn,
    impact_fn,
    inh_spill: _FamilySpill,
    pass_only: bool = False,
    segdup_iv=None,
    lcr_iv=None,
    collapse: bool = True,
    min_carrier_gq: int = 20,
) -> int:
    """families/<FID>/events.html for every family, in ONE sweep of the
    genotype matrices (no per-family rescan of the cohort VCF — that was the
    36 h trap on 1207-family cohorts). Returns the number of pages written."""
    sample_idx = {s: i for i, s in enumerate(samples)}
    fam_members: dict[str, list[tuple[str, int]]] = {
        fid: [(m, sample_idx[m]) for m in members if m in sample_idx]
        for fid, members in families.items()
        if 0 < len(members) <= _EVENTS_MAX_MEMBERS
    }
    # sample column -> family, for the carriers->families hop
    fam_of_col = np.full(len(samples), -1, np.int32)
    fids = sorted(fam_members)
    for k, fid in enumerate(fids):
        for _, i in fam_members[fid]:
            fam_of_col[i] = k

    ev_spill = _FamilySpill(run_dir / "families" / ".events_tmp")
    for site in sites:
        if pass_only and site["filter"] != "PASS":
            continue
        idx = site["idx"]
        # a carrier is a confident carrier: re-genotyping emits GT>0 with a
        # low GQ far too readily, and the trio signal says so plainly —
        # apparent de novo rate (a genotype-error proxy) on PMS is 16.7% for
        # DEL carriers at GQ<20 against 0.8% at GQ>=20, 13.9% vs 1.9% for DUP.
        carrier_cols = np.nonzero((nalt_m[idx] > 0) & (gq_m[idx] >= min_carrier_gq))[0]
        if not len(carrier_cols):
            continue
        fam_ks = np.unique(fam_of_col[carrier_cols])
        fam_ks = fam_ks[fam_ks >= 0]
        if not len(fam_ks):
            continue
        genes = genes_fn(site["contig"], site["start"], site["end"])
        impact = impact_fn(site["contig"], site["start"], site["end"])[1]
        def _num(v, fmt):
            return "" if v is None else format(float(v), fmt)
        # support du portage : le CN re-génotypé d'un porteur confirme-t-il le
        # type ? Mesuré sur les 77 trios PMS, le taux de de novo apparent (proxy
        # d'erreur de génotypage) tombe de 0,8 % à 0,3 % pour les DEL et de
        # 1,9 % à 0,3 % pour les DUP quand le CN est cohérent — mais la relecture
        # IGV montre que l'exiger coûte de vrais événements, d'où une colonne et
        # non un filtre par défaut.
        def _cn_ok(cn_val):
            if cn_m is None or cn_val < 0:
                return False
            return cn_val <= 1 if site["svtype"] == "DEL" else cn_val >= 3
        def _ovl(iv):
            return "" if iv is None else f"{iv.overlap_fraction(site['contig'], site['start'], site['end']):.2f}"
        base = "\t".join([
            site["id"] or "", site["contig"], str(site["start"] + 1), str(site["end"]),
            site["svtype"], str(site["size"]), site["filter"], site["rescue"],
            genes, impact,
            _num(site.get("sqmax"), ".1f"), _num(site.get("carrier_frac"), ".3f"),
            _num(site.get("ref_carrier_frac"), ".3f"),
            _ovl(lcr_iv), _ovl(segdup_iv),
        ])
        for k in fam_ks:
            fid = fids[int(k)]
            strong = any(
                nalt_m[idx, i] > 0 and gq_m[idx, i] >= min_carrier_gq
                and _cn_ok(int(cn_m[idx, i]) if cn_m is not None else -1)
                for _, i in fam_members[fid]
            )
            cells = "\t".join(
                f"{nalt_m[idx, i]}:{gq_m[idx, i]}:"
                + ("" if cn_m is None or cn_m[idx, i] < 0 else str(int(cn_m[idx, i])))
                for _, i in fam_members[fid]
            )
            ev_spill.add(fid, base + "\t" + ("fort" if strong else "faible") + "\t" + cells)
    ev_spill.flush()

    n_pages = 0
    for fid in fids:
        rows = ev_spill.read(fid)
        if not rows:
            continue
        inherit: dict[str, list[tuple[str, str]]] = defaultdict(list)
        for site_id, child, inh in inh_spill.read(fid):
            inherit[site_id].append((child, inh))
        members = [m for m, _ in fam_members[fid]]
        events = []
        for r in rows:
            site_id = r[0]
            inh = inherit.get(site_id, [])
            inh_s = " · ".join(
                f"{_INHERIT_FR.get(v, v)} ({c[-3:]})" for c, v in inh
            ) if len(inh) > 1 else (_INHERIT_FR.get(inh[0][1], inh[0][1]) if inh else "")
            gts, carriers = [], []
            for m, cell in zip(members, r[16:]):
                na, gq, cn = cell.split(":")
                gt_s = _NALT_GT.get(int(na), "./.")
                gts.append((gt_s, gq, cn))
                if gt_s not in ("0/0", "./.") and int(gq or 0) >= min_carrier_gq:
                    carriers.append(m)
            events.append({
                "id": site_id,
                "contig": r[1], "pos": int(r[2]), "end": int(r[3]), "svtype": r[4],
                "size": int(r[5]), "filter": r[6], "rescue": r[7], "genes": r[8],
                "impact": r[9], "sqmax": r[10], "cfrac": r[11], "refcf": r[12],
                "lcr": r[13], "segdup": r[14], "support": r[15],
                "inherit": inh_s, "carriers": carriers, "gt": gts,
            })
        fam_dir = run_dir / "families" / fid
        fam_dir.mkdir(parents=True, exist_ok=True)
        own_idx = own_call_index(run_dir, members)
        for e in events:
            e["own"] = own_call_status(own_idx, e["contig"], e["pos"], e["end"], e["svtype"])
        if collapse:
            events = collapse_events(events, genes_fn, impact_fn, segdup_iv, lcr_iv)
        _write_events_page(fam_dir / "events.html", f"Famille {fid}", members, events)
        n_pages += 1
    ev_spill.cleanup()
    return n_pages


def _write_events_page(
    out_file: Path, title: str, members: list[str], events: list[dict]
) -> None:
    """One interactive HTML page listing every DEL/DUP of a family workspace:
    chromosome and type selectors plus a PASS-only toggle filter the table
    client-side (no server, opens from disk). Genes/impact and per-member
    GT:GQ:SQ on each row.

    ``events`` rows: id, contig, pos (1-based), end, svtype, size, filter,
    rescue, genes, impact, inherit, carriers (list), gt (list of (GT, GQ, SQ)
    per member). Each row also carries an IGV button (locus sent to a running
    IGV desktop on its batch port) and a manual vrai/faux/incertain choice,
    remembered in the browser and exportable as TSV.
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
        "<th>rescue</th><th>SQmax</th><th>% port.</th><th>fréq. réf.</th>"
        "<th>LCR</th><th>segdup</th><th>sites</th><th data-key=\"own\">appel propre</th><th data-key=\"support\">support</th>"
        "<th>gènes</th><th>impact</th><th>héritage</th><th>porteurs</th>"
        '<th>IGV</th><th data-key="val">validation</th>'
        + "".join(f"<th>{html.escape(m[-3:])}</th>" for m in members) + "</tr>"
    )
    rows = []
    for e in sorted(events, key=lambda e: (contigs.index(e["contig"]), e["pos"])):
        is_pass = e["filter"] == "PASS"
        cls = e["svtype"].lower() + ("" if is_pass else " nopass")
        site = str(e.get("id") or f"{e['contig']}:{e['pos']}-{e['end']}")
        pad = max(500, round(e["size"] * 0.2))
        gt_cells = "".join(
            f'<td class="c">{html.escape(gt)}<span class="muted">:{gq}:{sq}</span></td>'
            for gt, gq, sq in e["gt"]
        )
        rows.append(
            f'<tr class="{cls}" data-chr="{html.escape(e["contig"])}" data-pos="{e["pos"]}" '
            f'data-end="{e["end"]}" data-size="{e["size"]}" data-site="{html.escape(site)}" '
            f'data-type="{e["svtype"]}" data-pass="{1 if is_pass else 0}" '
            f'data-own="{html.escape(e.get("own") or "aucun")}" '
            f'data-support="{html.escape(e.get("support") or "faible")}">'
            f'<td>{html.escape(e["contig"])}</td><td>{e["pos"]:,}</td><td>{e["end"]:,}</td>'
            f'<td class="t">{e["svtype"]}</td><td>{e["size"]:,}</td>'
            f'<td>{html.escape(e["filter"])}</td><td>{html.escape(e["rescue"])}</td>'
            f'<td>{html.escape(str(e.get("sqmax") or ""))}</td>'
            f'<td>{html.escape(str(e.get("cfrac") or ""))}</td>'
            f'<td>{html.escape(str(e.get("refcf") or ""))}</td>'
            f'<td>{html.escape(str(e.get("lcr") or ""))}</td>'
            f'<td>{html.escape(str(e.get("segdup") or ""))}</td>'
            f'<td title="{html.escape(" ".join(e.get("_ids", [])))}">{e.get("n_sites", 1)}</td>'
            f'<td>{html.escape(e.get("own") or "—")}</td>'
            f'<td>{html.escape(e.get("support") or "—")}</td>'
            f'<td class="g">{html.escape(e["genes"] or "—")}</td>'
            f'<td class="g">{html.escape(e["impact"])}</td>'
            f'<td class="g">{html.escape(e.get("inherit") or "—")}</td>'
            f'<td class="c">{",".join(m[-3:] for m in e["carriers"]) or "—"}</td>'
            f'<td><a class="igv" target="igv" title="Ouvrir dans IGV desktop (port 60151, ±20 %)" '
            f'href="http://localhost:60151/goto?locus={html.escape(e["contig"])}:{max(1, e["pos"] - pad)}-{e["end"] + pad}">IGV</a></td>'
            f'<td class="val"><label><input type="checkbox" class="v-true"> vrai</label> '
            f'<label><input type="checkbox" class="v-false"> faux</label> '
            f'<label><input type="checkbox" class="v-unsure"> incertain</label></td>{gt_cells}</tr>'
        )
    out_file.write_text(
        f"<!doctype html><title>{html.escape(title)} — événements</title>"
        f"<style>{_EVENTS_CSS}</style><h1>{html.escape(title)} — "
        f"{n_del} DEL + {n_dup} DUP</h1>"
        f'<div class="bar">'
        f'<label>chromosome <select id="chr">{chr_opts}</select></label>'
        f'<label>type <select id="type"><option value="all">tous</option>'
        f'<option>DEL</option><option>DUP</option></select></label>'
        f'<label>appel propre <select id="own">'
        f'<option value="any" selected>avec appel propre</option>'
        f'<option value="all">tous</option><option value="PASS">PASS seulement</option>'
        f'<option value="aucun">sans appel propre</option></select></label>'
        f'<label>support <select id="sup"><option value="all">tous</option>'
        f'<option value="fort">fort (CN cohérent)</option></select></label>'
        f'<label>validation <select id="val"><option value="all">toutes</option>'
        f'<option value="none">non validé</option><option value="true">vrai</option>'
        f'<option value="false">faux</option><option value="unsure">incertain</option></select></label>'
        f'<label><input type="checkbox" id="pass" checked> PASS uniquement</label>'
        f'<span class="muted"><b id="count"></b> lignes</span>'
        f'<button id="export" type="button">Exporter validations (TSV)</button>'
        f'<label class="muted">Importer <input type="file" id="import" accept=".tsv,.txt"></label></div>'
        f'<p class="muted">GT<span class="muted">:GQ:CN</span> par membre (suffixe du '
        f'barcode ; ordre : {html.escape(", ".join(members))}). Lignes grisées = non-PASS. '
        f'Porteur = génotype alternatif avec GQ suffisant (les génotypes moins sûrs restent '
        f'affichés mais ne comptent pas comme porteurs). Support « fort » = le CN re-génotypé '
        f"d'un porteur confirme le type (sur les trios PMS, taux de de novo apparent 0,3 % contre "
        f'0,8 % pour les DEL et 1,9 % pour les DUP sans cette condition) ; ce n\'est pas un filtre '
        f'par défaut car la relecture IGV montre qu\'il écarte aussi de vrais événements. '
        f'SQmax = meilleur SQ par échantillon au merge ; % port. = fraction de porteurs dans ce run ; '
        f'fréq. réf. = fraction de porteurs du site correspondant dans la cohorte de référence (--known-sites) ; '
        f'LCR / segdup = fraction du site chevauchant ces masques ; sites = nombre de sites cohorte '
        f'fusionnés dans la ligne (les points de cassure diffèrent d\'un patient à l\'autre, '
        f'survolez pour les identifiants) ; appel propre = un membre a-t-il son propre appel à ce '
        f'locus, ou le portage vient-il seulement du re-génotypage de cohorte (sur PMS, 80 % des '
        f'événements sont dans ce dernier cas et la relecture IGV les a tous jugés faux). '
        f'Le lien IGV envoie le locus (±20 %) à un IGV desktop ouvert sur cette machine '
        f'(View → Preferences → Advanced → « Enable port » 60151, CRAM de la famille chargés) ; '
        f'un onglet « OK » peut s\'ouvrir la première fois. Les cases vrai/faux/incertain sont '
        f'mémorisées dans ce navigateur et exportables en TSV (colonne <code>validation</code>).</p>'
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


def _read_cohort_vcf(path: Path, with_cn: bool = False):
    """One dict per site + compact genotype matrices.

    Dict-of-dicts genotypes cost ~80 bytes/genotype and OOM-killed the report
    on 238k sites x 2670 samples (635M genotypes); nalt/gq live in two numpy
    matrices (int8+int16 = 3 bytes/genotype, ~1.9 GB at that scale) and each
    site carries its row index. ``with_cn`` adds an int8 copy-number matrix
    (family events pages show GT:GQ:CN; 1 byte/genotype more).
    """
    sites = []
    nalt_rows, gq_rows, cn_rows = [], [], []
    with pysam.VariantFile(str(path)) as vcf:
        samples = list(vcf.header.samples)
        n = len(samples)
        # pysam raises on rec.info.get(key) when key is undeclared in the
        # header (rather than returning the default) — these three are only
        # present on cohorts that went through the current merge/cohort-filter
        # (SQMAX/SQMEAN from merge, CARRIER_FRAC/REF_CARRIER_FRAC from
        # cohort-filter); older cohorts re-filtered without a full re-merge
        # lack SQMAX entirely.
        has_sqmax = "SQMAX" in vcf.header.info
        has_carrier_frac = "CARRIER_FRAC" in vcf.header.info
        has_ref_carrier_frac = "REF_CARRIER_FRAC" in vcf.header.info
        for rec in vcf:
            af = rec.info.get("AF", (0.0,))
            af = af[0] if isinstance(af, tuple) else af
            nalt = np.full(n, -1, np.int8)   # -1 = missing
            gq = np.zeros(n, np.int16)
            cn = np.full(n, -1, np.int8) if with_cn else None
            for i, s in enumerate(samples):
                fmt = rec.samples[s]
                gt = fmt.get("GT")
                if gt is not None and None not in gt:
                    nalt[i] = sum(gt)
                    gq[i] = int(fmt.get("GQ") or 0)
                if with_cn:
                    v = fmt.get("CN")
                    if v is not None:
                        cn[i] = min(int(v), 127)
            sites.append({
                "id": rec.id, "contig": rec.contig, "start": rec.start, "end": rec.stop,
                "svtype": rec.info.get("SVTYPE", "?"), "size": rec.stop - rec.start,
                "af": float(af), "ns": int(rec.info.get("NS", 0)),
                "ncarrier": int(rec.info.get("NCARRIER", 0)),
                "hwe": float(rec.info.get("HWE_P", 1.0)),
                "filter": ";".join(rec.filter.keys()) or "PASS",
                "rescue": rec.info.get("RESCUE", "") or "",
                "sqmax": rec.info.get("SQMAX") if has_sqmax else None,
                "carrier_frac": rec.info.get("CARRIER_FRAC") if has_carrier_frac else None,
                "ref_carrier_frac": rec.info.get("REF_CARRIER_FRAC") if has_ref_carrier_frac else None,
                "idx": len(sites),
            })
            nalt_rows.append(nalt)
            gq_rows.append(gq)
            if with_cn:
                cn_rows.append(cn)
    nalt_m = np.vstack(nalt_rows) if nalt_rows else np.zeros((0, len(samples)), np.int8)
    gq_m = np.vstack(gq_rows) if gq_rows else np.zeros((0, len(samples)), np.int16)
    cn_m = (
        np.vstack(cn_rows) if cn_rows else np.zeros((0, len(samples)), np.int8)
    ) if with_cn else None
    return samples, sites, nalt_m, gq_m, cn_m


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


def pass_carrier_stats(
    nalt_m: np.ndarray, gq_m: np.ndarray, pass_mask: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """(PASS carrier count, fraction of them below GQ 20) per sample.

    Vectorized over the genotype matrices the report already holds — no extra
    cohort-VCF pass. The fraction is NaN for a sample with zero PASS carriers.
    """
    pass_nalt, pass_gq = nalt_m[pass_mask], gq_m[pass_mask]
    carrier_mask = pass_nalt > 0
    n_carrier = carrier_mask.sum(axis=0)
    n_low_gq = ((pass_gq < 20) & carrier_mask).sum(axis=0)
    low_gq_frac = np.divide(
        n_low_gq.astype(np.float64), n_carrier,
        out=np.full(nalt_m.shape[1], np.nan), where=n_carrier > 0,
    )
    return n_carrier, low_gq_frac


def carrier_excess(n_carrier: np.ndarray) -> np.ndarray:
    """Per-sample PASS carrier count over the cohort median (1.0 = typical).

    The count of sites a sample carries is what actually drives runaway de
    novo counts: a child carrying 1.6x the cohort median produces thousands
    of spurious candidates even with a normal genotype-quality profile
    (dataset_paper SMP-05/U/X: 487-497k carriers vs a 308k median, Spearman
    0.69 between carrier count and de novo count among the outliers).

    Deliberately NOT the low-GQ fraction: cohort-wide, 75% of PASS carrier
    genotypes sit below GQ 20 (median 0.75, range 0.50-0.89), and those three
    samples are at 0.61-0.69 — percentiles 2-6, i.e. among the *cleanest* on
    that axis. A first version of this guard thresholded the fraction and
    flagged 2670/2670 samples.
    """
    med = float(np.median(n_carrier)) if len(n_carrier) else 0.0
    if med <= 0:
        return np.full(len(n_carrier), np.nan)
    return n_carrier / med


def generate_report(
    run_dir: str | Path,
    pedigree: Pedigree | None,
    out_html: str | Path | None = None,
    cohort_name: str = "",
    make_family_workspaces: bool = True,
    min_denovo_gq: int = 20,
    genes_gff: str | Path | None = None,
    ndd_tsv: str | Path | None = None,
    segdup_bed: str | Path | None = None,
    lcr_bed: str | Path | None = None,
    max_denovo_per_child: int = 20,
    family_events: bool = False,
    events_pass_only: bool = False,
    collapse_family_events: bool = True,
    min_carrier_gq: int = 20,
) -> Path:
    run_dir = Path(run_dir)
    cohort_dir = run_dir / "cohort"
    vcf_path = cohort_dir / "cohort.vcf.gz"
    if not vcf_path.exists():
        raise FileNotFoundError(f"{vcf_path} absent — lancez le pipeline cohorte d'abord")
    out_html = Path(out_html) if out_html else run_dir / "report.html"
    cohort_name = cohort_name or run_dir.name

    gene_idx = _GeneIndex(genes_gff, ndd_tsv) if genes_gff else None
    # segdup/LCR are properties of the locus, so they are annotated here from
    # the site span rather than carried through the merge: that keeps them
    # available on already-merged cohorts. lcr_ovl was dead in the model up
    # to v5 (contig-naming bug), which is exactly why seeing it per event
    # matters for manual review.
    from .filter.annotate import IntervalSet

    segdup_iv = IntervalSet.from_bed(segdup_bed) if segdup_bed else None
    lcr_iv = IntervalSet.from_bed(lcr_bed) if lcr_bed else None

    def _genes(contig: str, start: int, end: int) -> str:
        return gene_idx.genes(contig, start, end) if gene_idx else ""

    def _impact(contig: str, start: int, end: int) -> tuple[int, str]:
        return gene_idx.impact(contig, start, end) if gene_idx else (0, "")

    samples, sites, nalt_m, gq_m, cn_m = _read_cohort_vcf(vcf_path, with_cn=family_events)
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
    # per-family (site, child, inheritance) spill for the events pages —
    # scales with trios x carried sites, so it goes to disk, not memory
    inh_spill = _FamilySpill(run_dir / "families" / ".inherit_tmp") if family_events else None
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
                if inh_spill is not None and r.get("inheritance"):
                    inh_spill.add(
                        r["family_id"], f'{r["site_id"]}\t{r["child"]}\t{r["inheritance"]}'
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
    # per-sample carrier-load guard: samples carrying far more PASS sites than
    # the cohort median are the ones whose children's de novo counts explode,
    # whatever the site-level filters do (see carrier_excess docstring). The
    # low-GQ fraction is kept as a reported column — informative about the
    # callset as a whole — but it does not discriminate between samples.
    pass_mask = np.array([site["filter"] == "PASS" for site in sites], dtype=bool)
    n_pass_carrier, low_gq_carrier_frac = pass_carrier_stats(nalt_m, gq_m, pass_mask)
    excess = carrier_excess(n_pass_carrier)
    sample_stats = {
        s: {
            "genotyped": int(genotyped_per_sample[i]),
            "carrier": int(carrier_per_sample[i]),
            "median_gq": float(median_gq_per_sample[i]),
            "low_gq_carrier_frac": (
                float(low_gq_carrier_frac[i]) if not np.isnan(low_gq_carrier_frac[i]) else None
            ),
            "carrier_excess": float(excess[i]) if not np.isnan(excess[i]) else None,
        }
        for i, s in enumerate(samples)
    }
    # Calibrated on dataset_paper (2670 samples). PASS-carrier load is wildly
    # dispersed — median 1.00 by construction, p90 4.4x, p95 6.4x, p99 9.6x,
    # max 23.5x — so a "slightly above median" threshold is meaningless here:
    #   1.3x -> 1037 flagged (38.8%), 3x -> 409 (15.3%),
    #   5x   ->  212 flagged  (7.9%), 8x ->  57  (2.1%)
    # 5x keeps 18 of the 20 worst de novo offenders while flagging under 8%
    # of the cohort (outliers' median load 1.41x vs 0.91x for the rest).
    # A cohort-relative criterion needs a real cohort behind the median: on a
    # handful of related samples everyone carries the same sites, so the
    # guard stays silent below CARRIER_EXCESS_MIN_SAMPLES.
    CARRIER_EXCESS_WARN = 5.0
    CARRIER_EXCESS_MIN_SAMPLES = 20
    overloaded = {}
    if len(samples) >= CARRIER_EXCESS_MIN_SAMPLES:
        overloaded = {
            s: v["carrier_excess"] for s, v in sample_stats.items()
            if v["carrier_excess"] is not None and v["carrier_excess"] >= CARRIER_EXCESS_WARN
        }
    if overloaded:
        logger.warning(
            "garde-fou charge de porteurs: %d échantillon(s) à >= %.1fx la médiane cohorte "
            "(%d sites PASS portés): %s",
            len(overloaded), CARRIER_EXCESS_WARN, int(np.median(n_pass_carrier)),
            ", ".join(
                f"{s} ({v:.2f}x)"
                for s, v in sorted(overloaded.items(), key=lambda kv: -kv[1])[:25]
            ),
        )

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
workspace : <code>families/{html.escape(fid)}/</code>{f' · <a href="families/{html.escape(fid)}/events.html">tous les événements (filtres chr/type)</a>' if family_events else ""}</p>
{f'<p class=warn>⚠ {", ".join(map(html.escape, fam_outliers))} : plus de {max_denovo_per_child} candidats de novo — profil technique atypique, candidats exclus du tableau.</p>' if fam_outliers else ""}
{_table(["site", "locus", "type", "taille", "gènes", "impact", "enfant", "GQ", "DNQ"], dn_rows, limit=20)
 if dn_rows else "<p>Aucun candidat de novo au seuil.</p>"}
</details>""")

    # ---- family events pages: one sweep of the matrices, never a per-family
    # rescan of the cohort VCF (the 36 h trap on 1207-family cohorts)
    if family_events and pedigree is not None and inh_spill is not None:
        inh_spill.flush()
        n_pages = _write_all_family_events(
            run_dir, dict(families), samples, sites, nalt_m, gq_m, cn_m,
            _genes, _impact, inh_spill, pass_only=events_pass_only,
            segdup_iv=segdup_iv, lcr_iv=lcr_iv, collapse=collapse_family_events,
            min_carrier_gq=min_carrier_gq,
        )
        inh_spill.cleanup()
        logger.info("events: %d pages familles écrites", n_pages)

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
             ("—" if st["low_gq_carrier_frac"] is None else f"{st['low_gq_carrier_frac']:.2f}"),
             ("—" if st["carrier_excess"] is None else
              f"{st['carrier_excess']:.2f} ⚠" if s in overloaded
              else f"{st['carrier_excess']:.2f}"),
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
<span class=kpi><b>{len(overloaded)}</b>échantillons surchargés (≥{CARRIER_EXCESS_WARN}x la médiane de porteurs PASS)</span>
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
{_table(["échantillon", "génotypés", "porteur", "GQ médian", "% porteurs GQ<20", "charge porteurs", "wave"], sample_rows)}

</body></html>"""

    out_html.write_text(page)
    logger.info(
        "rapport: %s (%d familles, %d sites rares)",
        out_html, len(family_sections), len(rare),
    )
    return out_html
