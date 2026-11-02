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
import bisect
import csv
import html
import io
import logging
import os
import re
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pysam

from . import __version__
from .karyotype import is_anomaly
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
td.qcw{color:#b3362b}
"""

_EVENTS_JS = """
const rows = Array.from(document.querySelectorAll("tbody tr"));
const chrSel = document.getElementById("chr"), typeSel = document.getElementById("type"),
      valSel = document.getElementById("val"), ownSel = document.getElementById("own"),
      supSel = document.getElementById("sup"), qcSel = document.getElementById("qc"),
      passCb = document.getElementById("pass"), subCb = document.getElementById("sub"),
      count = document.getElementById("count");
function apply(){
  const c = chrSel.value, t = typeSel.value, v = valSel.value, p = passCb.checked,
        o = ownSel.value, sp = supSel.value, q = qcSel.value;
  let n = 0;
  for (const r of rows){
    const ok = (c === "all" || r.dataset.chr === c)
            && (t === "all" || r.dataset.type === t)
            && (v === "all" || (r.dataset.val || "none") === v)
            && (o === "all" || (o === "any" ? r.dataset.own !== "aucun" : r.dataset.own === o))
            && (sp === "all" || r.dataset.support === sp)
            && (q === "all" || r.dataset.qc === q)
            && (subCb.checked || r.dataset.sub !== "1")
            && (!p || r.dataset.pass === "1" || (subCb.checked && r.dataset.sub === "1"));
    r.hidden = !ok;
    if (ok) n++;
  }
  count.textContent = n.toLocaleString("fr-FR") + " / " + rows.length.toLocaleString("fr-FR");
}
chrSel.onchange = typeSel.onchange = valSel.onchange = ownSel.onchange = supSel.onchange
  = qcSel.onchange = passCb.onchange = subCb.onchange = apply;

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
                  ["support", "support"], ["qc_porteur", "QC porteur"]];
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
  if (k) { const v = r.dataset[k] || ""; return NUM.test(v) ? parseFloat(v) : v; }
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


def _events_link(fid: str, n_members: int) -> str:
    """The family summary's link to its events page — or why there is none:
    dataset_paper's 38-member pedigree got a link to a page never written."""
    if n_members > _EVENTS_MAX_MEMBERS:
        return (f" · pas de page d'événements : {n_members} membres "
                f"(au-delà de {_EVENTS_MAX_MEMBERS} par page)")
    return (f' · <a href="families/{html.escape(fid)}/events.html">'
            f"tous les événements (filtres chr/type)</a>")


def own_call_index(run_dir: Path, members: list[str]) -> dict:
    """{(contig, svtype): [(start, stop, is_pass, member)]} of the family's OWN calls.

    A cohort site records a breakpoint cluster; re-genotyping then assigns a
    genotype to every sample, whether or not that sample ever produced a call
    there. On PMS, 80% of the PASS events a family "carries" have no call of
    their own behind them — and blind IGV review scored those 0/33, against
    3/9 for the ones backed by an own PASS call. The distinction has to be
    visible, so it is computed here from the published per-sample VCFs.
    """
    idx: dict[tuple[str, str], list[tuple[int, int, bool, str]]] = defaultdict(list)
    for m in members:
        vcf_path = run_dir / "samples" / m / f"{m}.filtered.vcf.gz"
        if not vcf_path.exists():
            continue
        try:
            with pysam.VariantFile(str(vcf_path)) as vf:
                lengths = {c.name: c.length for c in vf.header.contigs.values()}
                for rec in vf:
                    svtype = rec.info.get("SVTYPE")
                    if svtype:
                        # per-sample VCFs written before vcfout clamped END
                        # run up to one depth bin past the contig end, and
                        # own_span hands these ends straight to the rows
                        stop = rec.stop
                        if lengths.get(rec.contig):
                            stop = min(stop, lengths[rec.contig])
                        idx[(rec.contig, svtype)].append(
                            (rec.pos, stop, list(rec.filter.keys()) == ["PASS"], m)
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
    for rec in idx.get((contig, svtype), ()):
        op, oe, is_pass = rec[0], rec[1], rec[2]
        if op > end:
            break
        ov = min(end, oe) - max(pos, op)
        if ov > 0 and ov / max(end - pos, oe - op) >= min_ro:
            if is_pass:
                return "PASS"
            best = "filtré"
    return best


# Depth refinement of large rows (refine_bounds). Alexandre's IGV bounds for
# 15 rows (2026-10-13) showed row ends off by 0.15-10 Mb for three reasons the
# call chaining could not fix one by one: a cohort site from another
# proband's breakpoints taken as the row's start (46.07 vs 46.38 Mb), a chain
# cut where one piece's cohort site does not match (50.44 vs 49.88 Mb), a
# sub-threshold row showing one LowSQ piece of a ten-piece deletion (50.50 vs
# 49.28 Mb). The carriers' own depth says where the event is in all three.
REFINE_MIN_SIZE = 50_000
REFINE_WIN = 10_000            # window of the walk; bounds are then set per bin
REFINE_MAX_UNKNOWN = 250_000   # masked/no-data stretch the walk may cross
REFINE_DELTA = 0.6             # |CN - normal| for the event state
REFINE_NORMAL = 0.4            # |CN - normal| for the normal state
REFINE_MAX_NOISE = 2           # isolated off-state windows the walk may cross
# Second IGV review (2026-10-14, 39 adjusted rows): right 7/7 on carriers with
# wave < 0.128, 9/12 at 0.128-0.2, 9/18 above 0.2 (a coin toss); both DUPs
# (stretched from 112 kb to 23 Mb and from 134 kb to 72 Mb on a wave-0.127
# genome) were wrong. Deletions only, and not on very wavy carriers.
REFINE_MAX_WAVE = 0.2
REFINE_TYPES = ("DEL",)


class _DepthCache:
    """Per-sample, per-contig copy number (1 kb bins, NaN where masked) and
    the contig's normal copy number, read from the evidence stores."""

    def __init__(self, run_dir: Path):
        self.run_dir = Path(run_dir)
        self._cn: dict[tuple[str, str], tuple | None] = {}
        self._wave: dict[str, float | None] = {}

    def clear(self) -> None:
        self._cn.clear()
        self._wave.clear()

    def wave(self, sample: str) -> float | None:
        """The sample's wave QC score (store meta), None when unknown."""
        if sample not in self._wave:
            self._wave[sample] = None
            try:
                from .extract.store import EvidenceStore

                w = EvidenceStore(self.run_dir / "samples" / sample / f"{sample}.evidence.zarr").meta.get("wave_score")
                self._wave[sample] = None if w is None else float(w)
            except Exception:  # a report must survive a missing or foreign store
                pass
        return self._wave[sample]

    def get(self, sample: str, contig: str):
        key = (sample, contig)
        if key not in self._cn:
            self._cn[key] = None
            path = self.run_dir / "samples" / sample / f"{sample}.evidence.zarr"
            try:
                from .extract.store import EvidenceStore

                st = EvidenceStore(path)
                normal = int((st.meta.get("ploidy") or {}).get(contig, 2))
                mean_depth = float(st.meta.get("mean_depth") or 0)
                if normal > 0 and mean_depth > 0 and st.has_contig(contig):
                    _, corr, mask = st.depth(contig)
                    cn = 2.0 * np.asarray(corr, dtype=np.float64) / mean_depth
                    cn[np.asarray(mask, dtype=bool)] = np.nan
                    self._cn[key] = (cn, normal, st.bin_size)
            except Exception:  # a report must survive a missing or foreign store
                pass
        return self._cn[key]


def refine_bounds(depth: _DepthCache, contig: str, pos: int, end: int, svtype: str,
                  carriers, others):
    """(pos, end) of a large deletion row set on its carriers' depth, or None
    (duplications and carriers with wave >= REFINE_MAX_WAVE keep their bounds).

    Each window compares the carriers' copy number with the family's
    non-carriers (their median; the contig's normal copy number when there
    are none): a dip the whole family shares — subtelomeric chr22, a parent's
    coverage wave — is coverage, not the event. The walk goes out from each
    end while every carrier with data is REFINE_DELTA away from that
    reference (masked stretches up to REFINE_MAX_UNKNOWN are crossed), and an
    end where the carriers match the reference is trimmed inward. A non-carrier
    sharing the event leaves the window "normal" and stops the walk, so a row
    is not stretched over a neighbouring event other members carry."""
    if svtype not in REFINE_TYPES or end - pos + 1 < REFINE_MIN_SIZE:
        return None
    if any((depth.wave(s) or 0.0) >= REFINE_MAX_WAVE for s in carriers):
        return None
    tracks = [depth.get(s, contig) for s in carriers]
    if not tracks or any(t is None for t in tracks):
        return None
    bin_size = tracks[0][2]
    ref_tracks = [t for t in (depth.get(s, contig) for s in others) if t is not None]
    if not ref_tracks:
        # no family reference: in low-coverage stretches (acrocentric p arms)
        # an absolute CN 2 reference let singletons' rows grow by 0.7 Mb
        return None
    n_bins = min(len(t[0]) for t in tracks + ref_tracks)
    nb = max(1, REFINE_WIN // bin_size)
    sign = -1.0 if svtype == "DEL" else 1.0

    def value(cn, sl, min_n):
        x = cn[sl]
        x = x[np.isfinite(x)]
        return float(np.median(x)) if len(x) >= min_n else float("nan")

    def status_of(sl, min_n):
        ref_vals = [value(t[0], sl, min_n) for t in ref_tracks]
        ref_vals = [v for v in ref_vals if np.isfinite(v)]
        if not ref_vals or abs(float(np.median(ref_vals)) - tracks[0][1]) >= REFINE_DELTA:
            return "unknown"   # no family data, or the whole family is off: no information
        known = []
        for cn, _normal, _ in tracks:
            v = value(cn, sl, min_n)
            if not np.isfinite(v):
                continue
            ref = float(np.median(ref_vals))
            d = (v - ref) * sign
            known.append("event" if d >= REFINE_DELTA
                         else "normal" if abs(v - ref) <= REFINE_NORMAL else "other")
        if not known:
            return "unknown"
        if all(k == "event" for k in known):
            return "event"
        return "normal" if all(k == "normal" for k in known) else "other"

    def status(k):
        if k < 0 or (k + 1) * nb > n_bins:
            return "edge"
        return status_of(slice(k * nb, (k + 1) * nb), nb // 2)

    def walk(k0, step):
        """Last event window reached from k0 moving by `step`."""
        last, unknown, noise, k = k0, 0, 0, k0 + step
        while True:
            st = status(k)
            if st == "event":
                last, unknown, noise = k, 0, 0
            elif st == "unknown":
                unknown += 1
                if unknown * REFINE_WIN > REFINE_MAX_UNKNOWN:
                    break
            elif st == "other":   # one noisy window inside a deletion is not its end
                noise += 1
                if noise > REFINE_MAX_NOISE:
                    break
            else:
                break
            k += step
        return last

    ks, ke = (pos - 1) // REFINE_WIN, (end - 1) // REFINE_WIN
    mid = (ks + ke) // 2
    # an end window without information (a dip the whole family shares, a
    # partly masked bin) says nothing: decide from the next one inward
    for _ in range(2):
        if status(ks) in ("unknown", "other") and ks < mid:
            ks += 1
        if status(ke) in ("unknown", "other") and ke > mid:
            ke -= 1
    s_ks, s_ke = status(ks), status(ke)
    if s_ks == "event":
        new_ks = walk(ks, -1)
    elif s_ks == "normal":   # trim: first event window with two more behind it
        new_ks = next((k for k in range(ks, mid + 1)
                       if all(status(k + j) == "event" for j in range(3))), ks)
    else:
        new_ks = ks
    if s_ke == "event":
        new_ke = walk(ke, 1)
    elif s_ke == "normal":
        new_ke = next((k for k in range(ke, mid - 1, -1)
                       if all(status(k - j) == "event" for j in range(3))), ke)
    else:
        new_ke = ke

    def edge_bin(k, first, keep):
        """Outermost bin, in window k and its outer neighbour, from which the
        carriers' 5-bin median is at the event's copy number; `keep` when the
        window was not moved."""
        if keep is not None:
            return keep
        lo, hi = max(0, (k - 1) * nb), min(n_bins, (k + 2) * nb)
        rng = range(lo, hi) if first else range(hi - 1, lo - 1, -1)
        for b in rng:
            if status_of(slice(max(0, b - 2), b + 3), 3) == "event":
                return b
        return k * nb if first else min(n_bins, (k + 1) * nb) - 1

    b_start = edge_bin(new_ks, True, (pos - 1) // bin_size if new_ks == ks else None)
    b_end = edge_bin(new_ke, False, (end - 1) // bin_size if new_ke == ke else None)
    new_pos, new_end = b_start * bin_size + 1, (b_end + 1) * bin_size
    if new_end - new_pos + 1 < 0.5 * (end - pos + 1):
        return None   # never let the depth walk eat most of a called event
    if abs(new_pos - pos) <= REFINE_WIN and abs(new_end - end) <= REFINE_WIN:
        return None
    return (new_pos if new_ks != ks else pos), (new_end if new_ke != ke else end)



def own_span(idx: dict, contig: str, pos: int, end: int, svtype: str,
             carriers, min_ro: float = 0.5, min_call_frac: float = 0.1):
    """The carriers' own breakpoints for a cohort site, or None.

    A cohort site carries its cluster representative's coordinates, and in a
    cohort ascertained for a locus the cluster pools patients with different
    breakpoints. Re-genotyping then gives a carrier GQ 99 on a site that does
    not match its event, as long as most of the site has the right copy
    number in it: FAM-10 (own call from 45,146,001, site from
    44,116,001) and FAM-09 (own calls from 43,785,001, site from
    43,570,001) were drawn up to 1 Mb too long.

    The carriers' own calls overlapping the site are chained with the
    collapse rule (`_bridges`): a large depth event comes back in pieces —
    FAM-09's deletion in nine, none overlapping its site by half — so the
    chain, not a single call, stands for the carrier's event. The best chain
    replaces the site's span when it matches the site at reciprocal overlap
    >= `min_ro`, the rule of the "appel propre" column. The span can grow as
    well as shrink: FAM-11's duplication was merged into a site
    starting 141 kb inside it. The reciprocal overlap is the scale guard: a
    small site inside a large own call (a polymorphism the mother also
    carries) keeps its own bounds.

    Only PASS own calls set the chain's ends: allowing any own call moved 57k
    PMS rows on calls the filter had rejected. Filtered own calls still link
    PASS pieces inside a chain — the middle 138 kb of FAM-08's
    deletion is a LowSQ call, and without it the event split in two. A call
    joins only if at least `min_call_frac` of it lies in the site: a call
    touching the site's edge belongs to the neighbour (a 19 kb subtelomeric
    PASS call sharing 1 bp with its site pushed rows past the end of chr22).

    An own PASS call whose own cohort site failed the cohort filters (5th
    field False, see `_cohort_site_pass`) links but does not set ends either:
    FAM-14's 23 kb PASS call at 50,796,001-50,819,001 belongs to the
    subtelomeric site every PMS sample carries, and ran the row past chr22.
    Nor does a call nested in a longer one: FAM-06's 712 kb LowSQ
    duplication holds five PASS calls of 0.1-0.8 kb, which clipped the row to
    a sub-kb breakpoint inside the event and split the 5.6 Mb duplication.
    """
    calls = []
    for rec in idx.get((contig, svtype), ()):
        op, oe, is_pass = rec[0], rec[1], rec[2]
        if op > end:
            break
        who = rec[3] if len(rec) > 3 else None
        ov = min(end, oe) - max(pos, op)
        if who in carriers and ov > 0 and ov >= min_call_frac * (oe - op):
            anchors = is_pass and (len(rec) < 5 or rec[4] is not False)
            calls.append((op, oe, anchors))
    # longest first at equal starts, so a nested call always follows its host
    calls.sort(key=lambda c: (c[0], -c[1]))
    reach = None
    for k, (op, oe, _anchors) in enumerate(calls):
        if reach is not None and oe <= reach:
            calls[k] = (op, oe, False)          # nested: links, never an end
        reach = oe if reach is None else max(reach, oe)
    # chains: [start, end, first PASS start, last PASS end]
    chains = []
    for op, oe, anchors in calls:
        if chains and _bridges(chains[-1][0], chains[-1][1], op, oe):
            c = chains[-1]
            c[1] = max(c[1], oe)
        else:
            chains.append([op, oe, None, None])
            c = chains[-1]
        if anchors:
            c[2] = op if c[2] is None else min(c[2], op)
            c[3] = oe if c[3] is None else max(c[3], oe)
    spans = [(c[2], c[3]) for c in chains if c[2] is not None]
    if not spans:
        return None

    def ro(c):
        return (min(end, c[1]) - max(pos, c[0])) / max(end - pos, c[1] - c[0])

    best = max(spans, key=ro)
    return best if ro(best) >= min_ro else None


def _site_index(sites: list[dict]):
    """{(contig, svtype): (sorted [(pos, end, is_pass)], starts)} of cohort sites."""
    by: dict[tuple[str, str], list[tuple[int, int, bool]]] = defaultdict(list)
    for st in sites:
        by[(st["contig"], st["svtype"])].append((st["start"] + 1, st["end"], st["filter"] == "PASS"))
    out = {}
    for k, v in by.items():
        v.sort()
        out[k] = (v, [x[0] for x in v])
    return out


def _cohort_site_pass(site_idx: dict, contig: str, op: int, oe: int, svtype: str,
                      min_ro: float = 0.5):
    """FILTER==PASS of the cohort site an own call belongs to (best reciprocal
    overlap >= min_ro), or None when no site matches it."""
    entry = site_idx.get((contig, svtype))
    if entry is None:
        return None
    lst, starts = entry
    length = max(oe - op, 1)
    # RO >= 0.5 bounds a matching site to at most twice the call's length
    i = bisect.bisect_left(starts, op - length)
    j = bisect.bisect_right(starts, oe)
    best, best_ro = None, 0.0
    for sp, se, ok in lst[i:j]:
        ov = min(se, oe) - max(sp, op)
        if ov <= 0:
            continue
        r = ov / max(se - sp, length)
        if r >= min_ro and r > best_ro:
            best, best_ro = ok, r
    return best


def _impact_points(label: str) -> float:
    """Leading number of an impact label ("+1.90 — élevé · ..."), 0 if none."""
    m = re.match(r"\s*([+-]?\d+(?:\.\d+)?)", label or "")
    return float(m.group(1)) if m else 0.0


def _reannotate(e: dict, genes_fn=None, impact_fn=None, segdup_iv=None, lcr_iv=None) -> None:
    """Genes, impact and masks of an event whose span changed."""
    start0 = e["pos"] - 1
    if genes_fn is not None:
        e["genes"] = genes_fn(e["contig"], start0, e["end"])
    if impact_fn is not None:
        e["impact"] = impact_fn(e["contig"], start0, e["end"], e["svtype"])[1]
    if lcr_iv is not None:
        e["lcr"] = f"{lcr_iv.overlap_fraction(e['contig'], start0, e['end']):.2f}"
    if segdup_iv is not None:
        e["segdup"] = f"{segdup_iv.overlap_fraction(e['contig'], start0, e['end']):.2f}"


# A fragmented event leaves gaps where the HMM called CN=0 islands or where a
# segment failed a filter. Tolerating a gap proportional to the larger of the
# two neighbours closes them without chaining unrelated calls: measured over
# 25 PMS families it merges 0.9% of rows and the widest span stays 8.7 Mb
# (none above 15 Mb). The absolute cap bounds the worst case.
GAP_ABS = 2_000
GAP_FRAC = 0.02
GAP_MAX = 250_000


def _bridges(a_pos: int, a_end: int, b_pos: int, b_end: int) -> bool:
    """Whether b (starting at or after a) continues event a.

    Overlap always does. A gap is bridged when it is within the tolerance
    (GAP_*) AND the smaller neighbour is at least as long as the gap: a 0.6 kb
    polymorphism 21 kb off a 4.9 Mb deletion is its own event (FAM-10),
    whereas fragments of one event are large.
    """
    gap = b_pos - a_end
    if gap <= 0:
        return True
    tol = min(max(GAP_ABS, int(GAP_FRAC * max(a_end - a_pos, b_end - b_pos))), GAP_MAX)
    return gap <= tol and min(a_end - a_pos, b_end - b_pos) >= gap


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
    for (_contig, _svtype, _carriers), evs in groups.items():
        evs.sort(key=lambda e: (e["pos"], e["end"]))
        cur: dict | None = None
        for e in evs:
            if cur is not None and _bridges(cur["pos"], cur["end"], e["pos"], e["end"]):
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
                cur["_inh"].append(e.get("inherit") or "")
                if e.get("rescue"):
                    cur["rescue"] = e["rescue"]
                continue
            if cur is not None:
                out.append(cur)
            cur = dict(e, _ids=[e["id"]], _inh=[e.get("inherit") or ""])
        if cur is not None:
            out.append(cur)
    # re-annotate merged spans: genes, impact and masks describe the union
    for e in out:
        # inheritance of a merged row: the label most of its sites agree on,
        # not the concatenation of every variant seen. A terminal deletion
        # picks up a handful of sites a parent also carries, and listing them
        # all buried the fact that the event is de novo.
        labels = [x for x in e.pop("_inh", []) if x]
        if labels:
            counts = Counter(labels)
            top, n = counts.most_common(1)[0]
            e["inherit"] = top if n == len(labels) else f"{top} ({n}/{len(labels)})"
        e["n_sites"] = len(e["_ids"])
        e["size"] = e["end"] - e["pos"]
        if e["n_sites"] > 1:
            _reannotate(e, genes_fn, impact_fn, segdup_iv, lcr_iv)
    out.sort(key=lambda e: (e["contig"], e["pos"]))
    return out


# the page writer of the running _write_all_family_events, reached by its
# forked workers (a closure over the genotype matrices: shared, not pickled)
_FAMILY_PAGE = None


def _run_family_page(fid: str) -> int:
    return _FAMILY_PAGE(fid)


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
    sample_qc: dict[str, list[str]] | None = None,
    wave_scores: dict[str, float | None] | None = None,
    wave_by_contig: dict[str, dict] | None = None,
    mosaic_by_sample: dict[str, list[dict]] | None = None,
    subthreshold: SubCalls | None = None,
    denovo_families: dict[str, int] | None = None,
    low_quality_prefixes: tuple[str, ...] = (),
    jobs: int = 1,
) -> int:
    """families/<FID>/events.html for every family, in ONE sweep of the
    genotype matrices (no per-family rescan of the cohort VCF — that was the
    36 h trap on 1207-family cohorts). Returns the number of pages written.

    The pages are independent: with ``jobs`` > 1 forked workers write them
    (own calls, depth bounds and collapsing were ~40 % of the report on one
    core)."""
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
    depth = _DepthCache(run_dir)
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
        impact = impact_fn(site["contig"], site["start"], site["end"], site["svtype"])[1]
        def _num(v, fmt):
            return "" if v is None else format(float(v), fmt)
        # support du portage : le CN re-génotypé d'un porteur confirme-t-il le
        # type ? Mesuré sur les 77 trios PMS, le taux de de novo apparent (proxy
        # d'erreur de génotypage) tombe de 0,8 % à 0,3 % pour les DEL et de
        # 1,9 % à 0,3 % pour les DUP quand le CN est cohérent — mais la relecture
        # IGV montre que l'exiger coûte de vrais événements, d'où une colonne et
        # non un filtre par défaut.
        def _cn_ok(cn_val, site=site):
            if cn_m is None or cn_val < 0:
                return False
            return cn_val <= 1 if site["svtype"] == "DEL" else cn_val >= 3
        def _ovl(iv, site=site):
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

    site_idx = _site_index(sites)

    def _family_page(fid: str) -> int:
        rows = ev_spill.read(fid)
        if not rows and not any((subthreshold is not None and m in subthreshold)
                                or m in (mosaic_by_sample or {})
                                for m, _ in fam_members[fid]):
            return 0
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
            for m, cell in zip(members, r[16:], strict=False):
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
        own_idx = own_call_index(run_dir, members)
        own_idx = {
            k: [(op, oe, ok, who,
                 _cohort_site_pass(site_idx, k[0], op, oe, k[1]) if ok else None)
                for op, oe, ok, who in v]
            for k, v in own_idx.items()
        }
        for e in events:
            e["own"] = own_call_status(own_idx, e["contig"], e["pos"], e["end"], e["svtype"])
            span = own_span(own_idx, e["contig"], e["pos"], e["end"], e["svtype"],
                            set(e["carriers"]))
            if span is not None and span != (e["pos"], e["end"]):
                e["pos"], e["end"] = span
                e["size"] = e["end"] - e["pos"]
                _reannotate(e, genes_fn, impact_fn, segdup_iv, lcr_iv)
        if collapse:
            events = collapse_events(events, genes_fn, impact_fn, segdup_iv, lcr_iv)
        events += mosaic_events(members, mosaic_by_sample or {}, genes_fn, impact_fn,
                                segdup_iv, lcr_iv)
        sub_masked: dict[str, int] = {}
        if subthreshold:
            # IGV review (2026-10-13): sub-threshold rows were true 1/10 on old
            # IP genomes and 1/13 on wavy ones, 15/28 elsewhere
            skip = {m for m in members
                    if (low_quality_prefixes and m.startswith(tuple(low_quality_prefixes)))
                    or ((wave_scores or {}).get(m) or -1.0) >= WAVE_WARN}
            sub_rows, sub_masked = subthreshold_events(members, subthreshold, events,
                                                       genes_fn, impact_fn, segdup_iv, lcr_iv,
                                                       skip=skip)
            events += sub_rows
        if not events:
            return 0
        # large rows: bounds set on the carriers' depth (refine_bounds), then
        # rows that now cover the same event are merged again
        depth.clear()
        for e in events:
            if e.get("own") == "mosaïque" or not e["carriers"]:
                continue
            span = refine_bounds(depth, e["contig"], e["pos"], e["end"], e["svtype"],
                                 e["carriers"], [m for m in members if m not in e["carriers"]])
            if span is not None:
                e["pos"], e["end"] = span
                e["size"] = e["end"] - e["pos"]
                _reannotate(e, genes_fn, impact_fn, segdup_iv, lcr_iv)
        plain = [e for e in events if not e.get("sub") and e.get("own") != "mosaïque"]
        others_rows = [e for e in events if e.get("sub") or e.get("own") == "mosaïque"]
        if collapse:
            plain = collapse_events(plain, genes_fn, impact_fn, segdup_iv, lcr_iv)
        events = plain + [
            e for e in others_rows
            if not e.get("sub") or not any(
                p["contig"] == e["contig"] and p["svtype"] == e["svtype"]
                and min(p["end"], e["end"]) - max(p["pos"], e["pos"])
                > 0.5 * max(p["end"] - p["pos"], e["end"] - e["pos"])
                for p in plain)
        ]
        events += constitutional_events(members, mosaic_by_sample or {}, events,
                                        genes_fn, impact_fn, segdup_iv, lcr_iv)
        events = fold_mosaic_fragments(events, depth)
        big = big_cnv_counts(events)
        for e in events:
            alerts = []
            for m in e["carriers"]:
                fl = []
                w = row_wave_score((wave_by_contig or {}).get(m),
                                   (wave_scores or {}).get(m), e["contig"])
                if w is not None and w >= WAVE_WARN:
                    fl.append(f"wave {w:.3f}")
                fl += (sample_qc or {}).get(m, [])
                if big[m] >= BIG_CNV_WARN:
                    fl.append(f"grands CNV {big[m]}")
                fl += mosaic_flags((mosaic_by_sample or {}).get(m, []), e)
                if sub_masked.get(m):
                    fl.append(f"sous le seuil NDD {sub_masked[m]} masqués")
                if fl:
                    alerts.append(f"{m[-3:]} {', '.join(fl)}")
            if "de novo" in (e.get("inherit") or ""):
                ids = [e.get("id")] + list(e.get("_ids") or [])
                n_fam = max((denovo_families or {}).get(i, 0) for i in ids)
                alerts += denovo_alerts(e["size"], n_fam)
            e["qc"] = " · ".join(alerts)
        fam_dir = run_dir / "families" / fid
        fam_dir.mkdir(parents=True, exist_ok=True)
        _write_events_page(fam_dir / "events.html", f"Famille {fid}", members, events)
        return 1

    from .parallel import pmap

    global _FAMILY_PAGE
    _FAMILY_PAGE = _family_page
    try:
        n_pages = sum(pmap(_run_family_page, fids, jobs))
    finally:
        _FAMILY_PAGE = None
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
        "<th>gènes</th><th data-key=\"impact\">impact</th><th>héritage</th><th>porteurs</th>"
        '<th data-key="qc">QC porteur</th>'
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
            f'data-sub="{1 if e.get("sub") else 0}" '
            f'data-own="{html.escape(e.get("own") or "aucun")}" '
            f'data-impact="{_impact_points(e.get("impact") or "")}" '
            f'data-qc="{"alerte" if e.get("qc") else "ok"}" '
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
            + (f'<td class="c qcw">⚠ {html.escape(e["qc"])}</td>' if e.get("qc")
               else '<td class="c">—</td>')
            + f'<td><a class="igv" target="igv" title="Ouvrir dans IGV desktop (port 60151, ±20 %)" '
            f'href="http://localhost:60151/goto?locus={html.escape(e["contig"])}:{max(1, e["pos"] - pad)}-{e["end"] + pad}">IGV</a></td>'
            f'<td class="val"><label><input type="checkbox" class="v-true"> vrai</label> '
            f'<label><input type="checkbox" class="v-false"> faux</label> '
            f'<label><input type="checkbox" class="v-unsure"> incertain</label></td>{gt_cells}</tr>'
        )
    out_file.write_text(
        # without the charset a page opened from disk is decoded as
        # windows-1252 by Chrome ("filtrÃ©", "â€”" everywhere)
        f'<!doctype html><meta charset="utf-8"><title>{html.escape(title)} — événements</title>'
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
        f'<label>QC porteur <select id="qc"><option value="all">tous</option>'
        f'<option value="ok">sans alerte</option><option value="alerte">avec alerte</option>'
        f'</select></label>'
        f'<label>validation <select id="val"><option value="all">toutes</option>'
        f'<option value="none">non validé</option><option value="true">vrai</option>'
        f'<option value="false">faux</option><option value="unsure">incertain</option></select></label>'
        f'<label><input type="checkbox" id="pass" checked> PASS uniquement</label>'
        f'<label><input type="checkbox" id="sub" checked> + sous le seuil (gènes NDD)</label>'
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
        f'QC porteur = alertes QC de chaque porteur : wave au-delà du seuil (calculé sans le '
        f'chromosome de la ligne, qu\'un grand CNV fait osciller à lui seul), charge de sites PASS '
        f'≥ 5× la médiane de la cohorte, nombre de candidats de novo au-delà du garde-fou, ou '
        f'« grands CNV » : au moins 3 événements ≥ 1 Mb avec appel propre PASS chez ce porteur '
        f'(médiane de dataset_paper 0, p99 1), ou « mosaïque » : le porteur a un segment chromosomique '
        f'à nombre de copies fractionnaire confirmé par la BAF (vedette mosaic) — '
        f'un événement porté par un génome techniquement atypique est à relire avant d\'y croire, '
        f'quel que soit son impact. Sous le seuil (rescue « SOUS-SEUIL NDD ») = appel propre '
        f'rejeté par le seul score SQ, avec un CN franc (DEL CN ≤ 1, DUP CN ≥ 3, GQ ≥ 90, ≥ 1 kb), '
        f'touchant un gène NDD établi (impact ≥ 0,90) et absent de tout génome hors de la famille '
        f': une liste à relire dans IGV, pas des appels validés '
        f'(case « + sous le seuil » pour les masquer). Au-delà de 3 par génome, le génome est '
        f'jugé bruité : ses lignes sont masquées et ses autres lignes portent l\'alerte '
        f'« sous le seuil NDD n masqués ». Pas de ligne sous le seuil pour les génomes ondulés '
        f'(wave ≥ {WAVE_WARN}) ni pour les génomes anciens de moindre qualité (préfixes de barcode '
        f'configurés, IP au GHFC) : relecture IGV du 13/10, 1 vraie sur 10 chez les IP et 1 sur '
        f'13 chez les ondulés, contre 15 sur 28 ailleurs. Alertes de novo : « de novo < 1 kb » et '
        f'« de novo récurrent (n familles) » sous 10 kb — les 7 faux de novo de la relecture, '
        f'aucun des 56 vrais. Bornes des délétions d\'au moins 50 kb : ajustées sur la profondeur '
        f'des porteurs, comparée à celle des autres membres de la famille, comme on le fait à '
        f'l\'œil dans IGV (relecture du 13/10 : erreur au début de 0,5 Mb en médiane ramenée à '
        f'2 kb) ; pas pour les duplications ni sur les génomes très ondulés (wave ≥ {REFINE_MAX_WAVE}), '
        f'où la seconde relecture l\'a trouvé juste une fois sur deux. Les fragments d\'une mosaïque '
        f'(même porteur seul, dans le segment, au même niveau de copie) sont regroupés dans sa ligne '
        f'(colonne sites) ; un événement plus profond à l\'intérieur reste une ligne à part. '
        f'Rescue « CHROMOSOMIQUE » = gain ou perte d\'au moins 3 Mb à un nombre de copies entier '
        f'(aneuploïdie, bras, grande région), confirmé par la BAF (vedette mosaic) et qu\'aucun '
        f'appel PASS ne couvrait : à confirmer par caryotype ou ACPA. '
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
        self._gff_path = str(gff_path)
        self._tbx, self._tbx_pid = pysam.TabixFile(self._gff_path), os.getpid()
        self._contigs = set(self._tbx.contigs)
        # GeneTrek gene annotations, keyed by ENSG and by symbol
        self.gt: dict[str, _GeneTrek] = {}
        if ndd_tsv:
            with open(ndd_tsv) as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    g = _GeneTrek.from_row(row)
                    if g is None:
                        continue
                    for key in (row.get("ensg", ""), row.get("hgnc_symbol", "")):
                        if key:
                            self.gt[key] = g

    @property
    def tbx(self) -> pysam.TabixFile:
        """The GFF handle of THIS process: the family pages are written by
        forked workers, and an inherited htslib handle shares its file offset
        with the parent's and the siblings' (concurrent fetches would read
        each other's blocks)."""
        if self._tbx_pid != os.getpid():
            self._tbx, self._tbx_pid = pysam.TabixFile(self._gff_path), os.getpid()
        return self._tbx

    def _hits(self, contig: str, start: int, end: int):
        """[(name, ensg, coding, exon_overlap, gene_start, gene_end)] of genes
        overlapping the site (gene bounds 1-based inclusive, 0 if unknown)."""
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
                 "protein_coding" in f[8], False, 0, 0],
            )
            if f[2] == "exon":
                g[3] = True
            else:
                g[4], g[5] = int(f[3]), int(f[4])
                if "protein_coding" in f[8]:
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

    def impact(self, contig: str, start: int, end: int,
               svtype: str = "DEL") -> tuple[float, str]:
        """(points, label): gene-content score of a CNV, after the ACMG/ClinGen
        CNV scoring rules (Riggs et al. 2020, sections 1-3), with GeneTrek
        standing in for the ClinGen dosage map. Not an ACMG classification —
        inheritance and phenotype (sections 4-5) are left to the reader.

        The previous score (0 intergenic / 1 intronic / 2 exonic, +1/+2 for an
        NDD gene) had 5 values and ignored dosage, inheritance mode and the
        type of event: half the coding genome landed on the same level.
        """
        hits = [h for h in self._hits(contig, start, end) if h[2]]  # coding
        if not hits:
            return -0.60, "-0.60 — aucun gène codant"
        dup = svtype == "DUP"
        scored = []
        for name, ensg, _coding, exonic, gs, ge in hits:
            g = self.gt.get(ensg) or self.gt.get(name)
            whole = gs > 0 and start < gs and end >= ge
            inside = gs > 0 and start >= gs and end <= ge
            pts, rule = _gene_points(g, dup, whole, exonic, inside)
            scored.append((pts, name, rule, g))
        scored.sort(key=lambda t: -t[0])
        best_pts, best_name, best_rule, best_g = scored[0]
        n_coding = sum(1 for h in hits if h[3] or (h[4] > 0 and start < h[4] and end >= h[5]))
        # section 3: number of protein-coding genes (losses 25/35, gains 35/50)
        lo, hi = (35, 50) if dup else (25, 35)
        n_pts = 0.90 if n_coding >= hi else 0.45 if n_coding >= lo else 0.0
        total = round(best_pts + n_pts, 2)
        band = ("élevé" if total >= 0.90 else "modéré" if total >= 0.45
                else "faible" if total >= 0 else "aucun gène codant")
        parts = [f"{total:+.2f} — {band}"]
        if best_pts > 0:
            parts.append(f"{best_rule} {best_name} ({best_g.describe()})")
        else:
            exonic_any = any(h[3] for h in hits)
            parts.append("exonique" if exonic_any else "intronique")
        others = sum(1 for p, *_ in scored[1:] if p >= 0.45)
        if others:
            parts.append(f"+{others} gène{'s' if others > 1 else ''} NDD")
        if n_coding > 1:
            parts.append(f"{n_coding} gènes codants")
        return total, " · ".join(parts)


class _GeneTrek:
    """The GeneTrek columns the impact score reads, for one gene."""

    __slots__ = ("established", "candidate", "dominant", "recessive_only",
                 "constrained", "shet", "loeuf")

    _CANDIDATE_COLS = (
        "candidate_ndd_v5", "sfari_1_20260322", "autism_combined_v7",
        "panelapp_eng_id_20260331", "panelapp_aus_id_20260331",
        "panelapp_eng_epilepsy_20260331", "panelapp_aus_epilepsy_20260331",
        "panelapp_aus_autism_20260331",
    )
    _MODE_COLS = ("inheritance", "inheritance_panelapp_20260617",
                  "inheritance_sysndd_definitive_v0_11_2")

    @classmethod
    def from_row(cls, row: dict):
        def flag(c):
            return row.get(c, "").strip().lower() == "true"

        def num(c):
            try:
                return float(row.get(c, ""))
            except ValueError:
                return None

        g = cls()
        g.established = flag("hc_ndd_v7") or flag("sysndd_definitive_v0_11_2")
        g.candidate = any(flag(c) for c in cls._CANDIDATE_COLS)
        modes = " ".join(row.get(c, "") for c in cls._MODE_COLS).lower()
        g.dominant = "dominant" in modes
        g.recessive_only = "recessive" in modes and not g.dominant
        g.shet, g.loeuf = num("shet_mean"), num("loeuf")
        # haploinsufficiency predictors (ClinGen 2H): s_het or LOEUF
        g.constrained = (flag("is_shet_constrained")
                         or (g.shet is not None and g.shet >= 0.1)
                         or (g.loeuf is not None and g.loeuf < 0.35))
        if not (g.established or g.candidate or g.constrained):
            return None
        return g

    def describe(self) -> str:
        bits = ["NDD établi" if self.established else "NDD candidat" if self.candidate
                else "contraint"]
        if self.dominant:
            bits.append("AD")
        elif self.recessive_only:
            bits.append("AR")
        if self.shet is not None:
            bits.append(f"s_het {self.shet:.2f}")
        elif self.loeuf is not None:
            bits.append(f"LOEUF {self.loeuf:.2f}")
        return ", ".join(bits)


def _gene_points(g, dup: bool, whole: bool, exonic: bool, inside: bool) -> tuple[float, str]:
    """(points, rule) of one coding gene hit by a CNV; see _GeneIndex.impact.

    Losses (section 2): an established NDD gene with dominant inheritance is
    the stand-in for a ClinGen haploinsufficient gene — complete deletion
    1.00 (2A), exonic 0.90 (2C/2D, 5'/3' not distinguished). Unknown mode
    0.45, recessive only 0.15 (a heterozygous loss is a carrier state).
    Candidate NDD genes 0.30, 0.45 when also constrained; a constrained gene
    with no NDD evidence 0.15 (2H). Intronic: 0.

    Gains: triplosensitivity is not in GeneTrek, so a complete duplication
    scores 0.15 at most (dominant NDD gene, "TS non évaluée"). A duplication
    with both breakpoints inside a gene can disrupt it (2I) and is scored as
    a partial loss: 0.45 for a dominant NDD gene, 0.15 otherwise.
    """
    if g is None or not (whole or exonic):
        return 0.0, ""
    if not dup:
        rule = "DEL complète" if whole else "DEL exonique"
        if g.established and g.dominant:
            return (1.00 if whole else 0.90), rule
        if g.established:
            return (0.15 if g.recessive_only else 0.45), rule
        if g.candidate:
            return (0.45 if g.constrained else 0.30), rule
        return 0.15, rule
    if inside and exonic and not whole:
        rule = "DUP intragénique"
        if g.established and g.dominant:
            return 0.45, rule
        return 0.15, rule
    if whole and g.established and g.dominant:
        return 0.15, "DUP complète (TS non évaluée)"
    return 0.0, ""


def sample_qc_flags(overloaded: dict, dn_outliers: dict) -> dict[str, list[str]]:
    """Per-sample QC alerts shown next to the carriers of a family row.

    The impact score ranks events by gene content alone, so a technically
    noisy genome floats its artefacts to the top: on dataset_paper, family
    FAM-03 held a dozen of the 57 rows at >= +1.00 (PTEN, FGFR2...).
    Carrier load and de novo count are per sample; wave is per row (see
    `row_wave_score`).
    """
    flags: dict[str, list[str]] = defaultdict(list)
    for smp, x in overloaded.items():
        flags[smp].append(f"charge {x:.1f}×")
    for smp, n in dn_outliers.items():
        flags[smp].append(f"de novo {n}")
    return dict(flags)


# Large-CNV burden alert: events >= BIG_CNV_SIZE with the carrier's own PASS
# call, per sample. dataset_paper: median 0, p99 1, max 13 (SMP-02: 13
# non-transmitted deletions >= 1 Mb on 5 chromosomes, missed by wave, carrier
# load and the de novo guard alike); >= 3 flags 2 of 2673 samples.
BIG_CNV_SIZE = 1_000_000
BIG_CNV_WARN = 3


def big_cnv_counts(events: list[dict]) -> Counter:
    """{carrier: PASS events >= BIG_CNV_SIZE backed by an own PASS call} of
    one family page (a sample belongs to one family, so it is its count)."""
    n: Counter = Counter()
    for e in events:
        if (e["filter"] == "PASS" and e.get("own") == "PASS"
                and e["end"] - e["pos"] >= BIG_CNV_SIZE):
            for m in e["carriers"]:
                n[m] += 1
    return n


# A clone carrying an arm-level loss comes back from the integer CN-HMM as
# a string of CN 1 fragments; the family page showed them as rows beside the
# mosaic row (second IGV review, 2026-10-14: a clonal genome's 6p, 10q, 12p
# and 17p losses spread over 9 rows). A fragment carried by the mosaic's
# sample alone, lying in the segment and at the segment's copy number, is
# folded into the mosaic row; a deeper event inside it (a constitutional
# deletion within a mosaic one, CN 1.0 inside CN 1.3, also reviewed) stays.
FOLD_MIN_INSIDE = 0.5
FOLD_CN_TOL = 0.25


def fold_mosaic_fragments(events: list[dict], depth=None) -> list[dict]:
    """Rows folded into the mosaic row that contains them (see FOLD_*)."""
    mos = [e for e in events if e.get("own") in ("mosaïque", "chromosomique")]
    if not mos:
        return events
    folded: set[int] = set()
    for m in mos:
        carrier = m["carriers"][0]
        n_in = 0
        for k, e in enumerate(events):
            if (e is m or e.get("own") in ("mosaïque", "chromosomique") or k in folded
                    or e["carriers"] != [carrier] or e["contig"] != m["contig"]
                    or e["svtype"] != m["svtype"]):
                continue
            inside = min(e["end"], m["end"]) - max(e["pos"], m["pos"])
            if inside < FOLD_MIN_INSIDE * max(e["end"] - e["pos"], 1):
                continue
            if depth is not None and m.get("mosaic_cn") is not None:
                t = depth.get(carrier, e["contig"])
                if t is not None:
                    cn, _normal, bs = t
                    x = cn[(e["pos"] - 1) // bs:e["end"] // bs]
                    x = x[np.isfinite(x)]
                    if len(x) and abs(float(np.median(x)) - m["mosaic_cn"]) > FOLD_CN_TOL:
                        continue   # its own level: a distinct event inside the mosaic
            folded.add(k)
            n_in += 1
            m.setdefault("_ids", []).extend([e.get("id")] + list(e.get("_ids") or []))
        if n_in:
            m["n_sites"] = m.get("n_sites", 1) + n_in
            m["inherit"] += f" · {n_in} fragment{'s' if n_in > 1 else ''} regroupé{'s' if n_in > 1 else ''}"
    return [e for k, e in enumerate(events) if k not in folded]


def constitutional_events(members: list[str], mosaic_by_sample: dict[str, list[dict]],
                          shown: list[dict], genes_fn=None, impact_fn=None,
                          segdup_iv=None, lcr_iv=None) -> list[dict]:
    """Family-page rows for the scan's whole-copy segments (aneuploidies,
    arm-level and large constitutional gains or losses) that no PASS row of
    the carrier already covers by half: a trisomy 21 and a trisomy 13 came
    back from the integer HMM as dozens of LowSQ pieces and were on no page.
    A segment the PASS rows do cover (a PMS terminal deletion) keeps the
    calls' bounds, which are finer than the scan's 1 Mb windows."""
    rows = []
    for m in members:
        for g in mosaic_by_sample.get(m, []):
            if g["mosaic"]:
                continue
            pos, end = g["start"] + 1, g["end"]
            covered = sum(
                max(0, min(e["end"], end) - max(e["pos"], pos))
                for e in shown
                if e["contig"] == g["contig"] and e["svtype"] == g["svtype"]
                and e.get("filter") == "PASS" and m in e["carriers"]
                and e.get("own") not in ("mosaïque", "chromosomique")
            )
            if covered >= 0.5 * (end - pos):
                continue
            also = sorted(
                x[-3:] for x in members if x != m and any(
                    h["contig"] == g["contig"] and h["svtype"] == g["svtype"]
                    and min(h["end"], end) - max(h["start"] + 1, pos) > 0.5 * (end - pos)
                    for h in mosaic_by_sample.get(x, []))
            )
            kind = "gain" if g["svtype"] == "DUP" else "perte"
            e = {
                "id": f"chrom_{m}_{g['contig']}_{g['start']}",
                "contig": g["contig"], "pos": pos, "end": end,
                "svtype": g["svtype"], "size": end - g["start"],
                "filter": "PASS", "rescue": "CHROMOSOMIQUE", "genes": "", "impact": "",
                "sqmax": "", "cfrac": "", "refcf": "", "lcr": "", "segdup": "",
                "support": "fort", "own": "chromosomique",
                "inherit": f"{kind} constitutionnel (CN {g['cn']:.2f}, BAF)"
                           + (f" · aussi chez {', '.join(also)}" if also else ""),
                "mosaic_cn": float(g["cn"]),
                "carriers": [m],
                "gt": [("chr", "", f"{g['cn']:.2f}") if x == m else ("—", "", "") for x in members],
                "n_sites": 1, "_ids": [],
            }
            _reannotate(e, genes_fn, impact_fn, segdup_iv, lcr_iv)
            rows.append(e)
    return rows


def mosaic_events(members: list[str], mosaic_by_sample: dict[str, list[dict]],
                  genes_fn=None, impact_fn=None, segdup_iv=None, lcr_iv=None) -> list[dict]:
    """Family-page rows for the members' mosaic segments (`vedette mosaic`).

    A mosaic loss or gain has no cohort site of its own: the integer CN-HMM
    either fragments it or, close to CN 2, does not call it at all (a PMS
    proband's terminal 22q13 loss at CN 1.62 had no call; it was only in the
    report's mosaic table, so the family page showed nothing there). Each
    segment becomes one row, carried by its sample, visible by default.
    """
    rows = []
    for m in members:
        for g in mosaic_by_sample.get(m, []):
            if not g["mosaic"]:
                continue
            frac = round(100 * g["fraction"])
            e = {
                "id": f"mosaic_{m}_{g['contig']}_{g['start']}",
                "contig": g["contig"], "pos": g["start"] + 1, "end": g["end"],
                "svtype": g["svtype"], "size": g["end"] - g["start"],
                "filter": "PASS", "rescue": "MOSAÏQUE", "genes": "", "impact": "",
                "sqmax": "", "cfrac": "", "refcf": "", "lcr": "", "segdup": "",
                "support": "—", "own": "mosaïque",
                "inherit": f"mosaïque {frac} % (CN {g['cn']:.2f})",
                "mosaic_cn": float(g["cn"]),
                "carriers": [m],
                "gt": [("mos", f"{frac}%", f"{g['cn']:.2f}") if x == m else ("—", "", "")
                       for x in members],
                "n_sites": 1, "_ids": [],
            }
            _reannotate(e, genes_fn, impact_fn, segdup_iv, lcr_iv)
            rows.append(e)
    return rows


# Sub-threshold tier: a call the filter model under-scores but whose depth is
# unambiguous, private to its family and hitting an established dominant NDD
# gene. A PMS proband's 9 kb deletion of the 3' end of SHANK3 (CN 1, RD 0.50,
# GQ 99, confirmed in IGV) had SQ 0.38 for a threshold of 5: depth-only, no
# heterozygous site to read a BAF from in the GC-rich 22q13.33, noisy depth —
# v5's known weak spot on small depth-only deletions. A review list, not a
# PASS set.
SUB_MIN_GQ = 90
SUB_MIN_SIZE = 1_000
SUB_MIN_IMPACT = 0.90
# ... but they pile up in a few noisy genomes: without a cap, dataset_paper
# had 254 rows in 39 families, 71 in one genome (median 89 kb, no QC alert)
# and 53 in another; PMS 46 rows, 25 of them in one family's parent and
# child. A real sub-threshold event in an NDD gene is a one-off (per-genome
# counts on dataset_paper: 27 genomes with 1, 3 with 2, 2 with 3, then 4, 5,
# 8 ... 71); past SUB_MAX_PER_SAMPLE the genome's rows are masked and its
# other rows carry a QC alert instead. Left: PMS 8 rows, dataset_paper 39.
SUB_MAX_PER_SAMPLE = 3


class SubCalls:
    """Sub-threshold candidates of a cohort (see `subthreshold_candidates`) in
    compact arrays. dataset_paper has ~6M of them (2,240 per genome): held as
    one dict each, plus an up-front sharer map, they took the report past
    32 GB; here they cost ~25 bytes each and sharers are looked up on demand,
    for the candidates of the family being written only."""

    _DT = np.dtype([("ci", np.int32), ("dup", np.int8), ("start", np.int64),
                    ("end", np.int64), ("cn", np.int8), ("gq", np.int16), ("sq", np.float32)])

    def __init__(self):
        self.samples: list[str] = []
        self.contigs: list[str] = []
        self._ci: dict[str, int] = {}
        self._rows: dict[str, np.ndarray] = {}
        self._chunks: dict[tuple[int, int], list] = defaultdict(list)
        self._by: dict[tuple[int, int], tuple[np.ndarray, np.ndarray, np.ndarray]] = {}

    def add(self, sample: str, calls: list[tuple]) -> None:
        """calls: (contig, svtype, start, end, cn, gq, sq) of one sample."""
        si = len(self.samples)
        self.samples.append(sample)
        rows = np.array([(self._ci.setdefault(c, len(self._ci)), t == "DUP", a, b, cn, gq,
                          np.nan if sq is None else sq)
                         for c, t, a, b, cn, gq, sq in calls], dtype=self._DT)
        self._rows[sample] = rows
        for ci, dup in set(zip(rows["ci"].tolist(), rows["dup"].tolist(), strict=True)):
            sel = rows[(rows["ci"] == ci) & (rows["dup"] == dup)]
            self._chunks[(ci, dup)].append((sel["start"], sel["end"], np.full(len(sel), si, np.int32)))

    def freeze(self) -> SubCalls:
        self.contigs = sorted(self._ci, key=self._ci.get)
        for k, parts in self._chunks.items():
            st = np.concatenate([x[0] for x in parts])
            en = np.concatenate([x[1] for x in parts])
            sm = np.concatenate([x[2] for x in parts])
            o = np.argsort(st, kind="stable")
            self._by[k] = (st[o], en[o], sm[o])
        self._chunks = defaultdict(list)
        self._index = {s: i for i, s in enumerate(self.samples)}
        return self

    def __contains__(self, sample: str) -> bool:
        return sample in self._rows

    def __bool__(self) -> bool:
        return bool(self._rows)

    def __len__(self) -> int:
        return sum(len(r) for r in self._rows.values())

    def calls(self, sample: str) -> list[dict]:
        return [{"sample": sample, "contig": self.contigs[r["ci"]],
                 "svtype": "DUP" if r["dup"] else "DEL", "start": int(r["start"]),
                 "end": int(r["end"]), "cn": int(r["cn"]), "gq": int(r["gq"]),
                 "sq": None if np.isnan(r["sq"]) else float(r["sq"])}
                for r in self._rows.get(sample, ())]

    def sharers(self, sample: str, contig: str, start: int, end: int, svtype: str,
                min_ro: float = 0.5) -> set[str]:
        """Other samples with a candidate at reciprocal overlap >= min_ro."""
        entry = self._by.get((self._ci.get(contig, -1), svtype == "DUP"))
        if entry is None:
            return set()
        st, en, sm = entry
        length = end - start
        # RO >= 0.5 bounds a match to at most twice the call's length
        i = np.searchsorted(st, start - length, "left")
        j = np.searchsorted(st, end, "right")
        s_, e_, x_ = st[i:j], en[i:j], sm[i:j]
        ov = np.minimum(e_, end) - np.maximum(s_, start)
        hit = (ov > 0) & (ov >= min_ro * np.maximum(length, e_ - s_)) & (x_ != self._index[sample])
        return {self.samples[k] for k in np.unique(x_[hit])}


def _sample_subthreshold(vcf_path: Path) -> list[tuple] | None:
    """One genome's sub-threshold candidates (see subthreshold_candidates);
    None when its VCF is missing or unreadable."""
    if not vcf_path.exists():
        return None
    calls = []
    try:
        with pysam.VariantFile(str(vcf_path)) as vf:
            lengths = {c.name: c.length for c in vf.header.contigs.values()}
            for rec in vf:
                if list(rec.filter.keys()) != ["LowSQ"]:
                    continue
                svtype = rec.info.get("SVTYPE")
                fmt = rec.samples[0]
                cn, gq = fmt.get("CN"), fmt.get("GQ")
                if cn is None or gq is None or gq < SUB_MIN_GQ:
                    continue
                if not ((svtype == "DEL" and cn <= 1) or (svtype == "DUP" and cn >= 3)):
                    continue
                end = min(rec.stop, lengths.get(rec.contig) or rec.stop)
                if end - rec.start < SUB_MIN_SIZE:
                    continue
                calls.append((rec.contig, svtype, rec.start, end, min(int(cn), 127),
                              min(int(gq), 32767), fmt.get("SQ")))
    except (OSError, ValueError):  # a report must survive a missing VCF
        return None
    return calls


def subthreshold_candidates(run_dir: Path, samples: list[str], jobs: int = 1) -> SubCalls:
    """The calls rejected by the score alone (FILTER is exactly LowSQ) with a
    clear copy number (DEL CN <= 1, DUP CN >= 3), GQ >= SUB_MIN_GQ and >=
    SUB_MIN_SIZE, from the published per-sample VCFs (read in parallel with
    ``jobs`` > 1, added in sample order)."""
    from .parallel import pmap

    out = SubCalls()
    paths = [run_dir / "samples" / s / f"{s}.filtered.vcf.gz" for s in samples]
    for s, calls in zip(samples, pmap(_sample_subthreshold, paths, jobs), strict=True):
        if calls:
            out.add(s, calls)
    return out.freeze()


def subthreshold_events(members: list[str], cands: SubCalls,
                        shown: list[dict], genes_fn=None, impact_fn=None,
                        segdup_iv=None, lcr_iv=None,
                        skip: set[str] | frozenset = frozenset()) -> tuple[list[dict], dict[str, int]]:
    """(rows, masked): family-page rows for the members' sub-threshold calls
    (see SUB_*) that hit an established NDD gene (impact >= SUB_MIN_IMPACT),
    are found in no genome outside the family and are not already on the
    page; `masked` = {member: n} of the members past SUB_MAX_PER_SAMPLE, whose
    rows are left out."""
    fam = set(members)
    rows, masked = [], {}
    for m in members:
        if m in skip:
            continue
        mine = []
        for c in cands.calls(m):
            # cheapest test first: ~2,240 candidates per genome, most shared
            others = cands.sharers(m, c["contig"], c["start"], c["end"], c["svtype"])
            if others - fam:
                continue
            pos, end = c["start"] + 1, c["end"]
            if impact_fn is None or impact_fn(c["contig"], c["start"], end, c["svtype"])[0] < SUB_MIN_IMPACT:
                continue
            if any(e["contig"] == c["contig"] and e["svtype"] == c["svtype"]
                   and min(e["end"], end) - max(e["pos"], pos) > 0.5 * max(e["end"] - e["pos"], end - pos)
                   for e in shown):
                continue
            also = sorted(x[-3:] for x in others & fam)
            sq = "" if c["sq"] is None else f"{c['sq']:.2f}"
            e = {
                "id": f"sub_{m}_{c['contig']}_{c['start']}",
                "contig": c["contig"], "pos": pos, "end": end, "svtype": c["svtype"],
                "size": end - c["start"], "filter": "LowSQ", "rescue": "SOUS-SEUIL NDD",
                "genes": "", "impact": "", "sqmax": sq, "cfrac": "", "refcf": "",
                "lcr": "", "segdup": "", "support": "fort", "own": "LowSQ",
                "inherit": f"sous le seuil (SQ {sq})"
                           + (f" · appel aussi chez {', '.join(also)}" if also else ""),
                "carriers": [m],
                "gt": [("0/1" if c["svtype"] == "DUP" or c["cn"] == 1 else "1/1",
                        str(c["gq"]), str(c["cn"])) if x == m else ("—", "", "")
                       for x in members],
                "n_sites": 1, "_ids": [], "sub": True,
            }
            _reannotate(e, genes_fn, impact_fn, segdup_iv, lcr_iv)
            mine.append(e)
        if len(mine) > SUB_MAX_PER_SAMPLE:
            masked[m] = len(mine)
        else:
            rows += mine
    return rows, masked


def mosaic_flags(segments: list[dict], e: dict) -> list[str]:
    """Alerts of a carrier with mosaic segments (`vedette mosaic`): the row's
    own segment when the row lies in one, else a sample-level warning.

    A clone present in a fraction of the cells makes the integer CN-HMM emit
    fragments: a dataset_paper parent's 50 Mb 10q loss at CN 1.33 came back
    as 57 CN=1 rows. The row says which mosaic it belongs to."""
    mos = [g for g in segments if g["mosaic"]]
    if not mos:
        return []
    for g in mos:
        if (g["contig"] == e["contig"] and g["svtype"] == e.get("svtype", g["svtype"])
                and min(g["end"], e["end"]) - max(g["start"], e["pos"] - 1) > 0):
            return [f"mosaïque {round(100 * g['fraction'])} % "
                    f"{g['contig']}:{g['start'] / 1e6:.1f}-{g['end'] / 1e6:.1f} Mb"]
    return [f"porteur mosaïque ({len(mos)} segment{'s' if len(mos) > 1 else ''})"]


def row_wave_score(wave_by_contig: dict | None, sample_score: float | None,
                   contig: str, min_contigs: int = 10) -> float | None:
    """A carrier's wave score for one row: the mean of its per-contig scores
    without the row's own contig.

    The score is low-frequency spectral energy, and a large CNV is a step in
    coverage: on PMS, 26 of the 48 carriers of a >= 1 Mb chr22 deletion were
    "wavy" (54%, against 13% of parents), their chr22 alone at 0.4-0.6. Left
    out of the mean, 5 remain; a genuinely wavy genome oscillates on every
    contig and keeps its flag. Falls back to the sample score on stores
    without per-contig values.
    """
    if wave_by_contig:
        vals = [v for k, v in wave_by_contig.items() if k != contig]
        if len(vals) >= min_contigs:
            return sum(vals) / len(vals)
    return sample_score


def denovo_event_count(rows: list[dict]) -> int:
    """De novo candidate events of one child: candidate sites merged by type
    with the collapse rule (`_bridges`), so that one event spanning several
    cohort sites (patients' breakpoints differ) counts once. Measured on PMS
    the effect is small — the children's candidates are mostly distinct
    small events (SMP-01: 116 sites, 114 events).
    """
    by: dict[tuple[str, str], list[tuple[int, int]]] = defaultdict(list)
    for r in rows:
        by[(r["contig"], r["svtype"])].append((int(r["start"]), int(r["end"])))
    n = 0
    for ivs in by.values():
        ivs.sort(key=lambda iv: (iv[0], -iv[1]))
        cur = None
        for a, b in ivs:
            if cur is not None and _bridges(cur[0], cur[1], a, b):
                cur = (cur[0], max(cur[1], b))
            else:
                n += 1
                cur = (a, b)
    return n


# Per-child de novo guard. An absolute cutoff stopped meaning anything once
# the callset was measured: with > 20 events it flagged 664 of 665
# dataset_paper children and all 77 PMS children (and hid their de novo
# tables). Cohort-relative instead, as the carrier-load guard:
#   dataset_paper  median 75, p90 367, p99 3173 -> 5x median flags 63 (9.5%)
#   PMS            median 100, max 576          -> 5x median flags 1
# The absolute floor stays the cutoff below DENOVO_REL_MIN_CHILDREN children
# (family runs), where a median says nothing.
DENOVO_REL = 5.0
DENOVO_REL_MIN_CHILDREN = 20


# De novo row alerts, from Alexandre's IGV review of 63 high-impact de novo
# rows (2026-10-13): all 7 false ones were small (177 bp - 3.1 kb; one 3.1 kb
# site in 3 unrelated families), the smallest true one 10.5 kb. Recurrence
# alone is no sign of artefact: PMS's true 22q13 deletions share sites across
# up to 40 families, dataset_paper's SHANK2 and 16p11.2 across 2-3. Flagging
# < 1 kb, and recurrent below 10 kb, catches the 7 and none of the 56 true.
DENOVO_SMALL = 1_000
DENOVO_RECURRENT_MAX = 10_000


def denovo_alerts(size: int, n_families: int) -> list[str]:
    out = []
    if size < DENOVO_SMALL:
        out.append("de novo < 1 kb")
    if n_families >= 2 and size < DENOVO_RECURRENT_MAX:
        out.append(f"de novo récurrent ({n_families} familles)")
    return out


def collapse_denovo(rows: list[dict], own_idx_of=None, site_idx=None) -> list[dict]:
    """A family's de novo candidate sites as events, one per child and event.

    In a cohort ascertained for a locus, the merge makes one site per
    breakpoint cluster, and a child carrying the largest deletion is
    re-genotyped as carrier of every site nested in it: a PMS proband's
    8.2 Mb terminal 22q13 deletion filled the report's de novo table with 61
    sites. Each site first takes the child's own bounds (`own_span`, as on the
    family pages), then sites of the same child and type that overlap or
    bridge (`_bridges`, the family pages' rule) become one event, with the
    union of their bounds."""
    by: dict[tuple, list] = defaultdict(list)
    for r in rows:
        pos, end = int(r["start"]) + 1, int(r["end"])
        if own_idx_of is not None:
            idx = own_idx_of(r["child"])
            span = own_span(idx, r["contig"], pos, end, r["svtype"], {r["child"]})
            if span is not None:
                pos, end = span
        by[(r["child"], r["contig"], r["svtype"])].append((pos, end, r))
    out = []
    for (child, contig, svtype), items in by.items():
        items.sort(key=lambda t: t[0])
        cur = None
        for pos, end, r in items:
            if cur is not None and _bridges(cur["pos"], cur["end"], pos, end):
                cur["end"] = max(cur["end"], end)
                cur["sites"].append(r)
            else:
                cur = {"child": child, "contig": contig, "svtype": svtype,
                       "pos": pos, "end": end, "sites": [r]}
                out.append(cur)
    return out


def denovo_family_counts(denovo_by_family: dict[str, list[dict]]) -> dict[str, int]:
    """{site_id: families with a de novo candidate at that site}."""
    fams: dict[str, set] = defaultdict(set)
    for fid, rows in denovo_by_family.items():
        for r in rows:
            fams[r["site_id"]].add(fid)
    return {k: len(v) for k, v in fams.items()}


def denovo_counts(denovo_by_family: dict[str, list[dict]]) -> dict[str, int]:
    """{child: de novo candidate events}, across families."""
    per_child: dict[str, list[dict]] = defaultdict(list)
    for rows in denovo_by_family.values():
        for r in rows:
            per_child[r["child"]].append(r)
    return {c: denovo_event_count(rs) for c, rs in per_child.items()}


def denovo_cutoff(counts: dict[str, int], floor: int) -> float:
    """Count above which a child is a technical outlier: DENOVO_REL x the
    cohort median, never below `floor`."""
    if len(counts) < DENOVO_REL_MIN_CHILDREN:
        return float(floor)
    return max(float(floor), DENOVO_REL * float(np.median(list(counts.values()))))


def denovo_outliers(
    denovo_by_family: dict[str, list[dict]], threshold: int
) -> dict[str, int]:
    """Children with more de novo candidate events than `denovo_cutoff`.

    A genuine genome carries <1 de novo CNV on average; a child far above its
    cohort is a technical outlier (waviness under the QC bar, contamination,
    swap) whose candidates would drown the family tables — measured on
    dataset_paper: two such children held 39 of the 81 top-impact candidates.
    """
    counts = denovo_counts(denovo_by_family)
    cut = denovo_cutoff(counts, threshold)
    return {c: n for c, n in counts.items() if n > cut}


def _read_cohort_region(task):
    """Sites and genotype rows of one region (all of the file when the
    region is None) — see _read_cohort_vcf."""
    from .parallel import fetch_region

    path, region, with_cn = task
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
        lengths = {c.name: c.length for c in vcf.header.contigs.values()}
        records = fetch_region(vcf, region) if region is not None else vcf
        for rec in records:
            # sites merged before merge clamped END run up to one depth bin
            # past the contig end, and family rows without an own span keep
            # the site's ends (same bound as own_call_index)
            stop = rec.stop
            if lengths.get(rec.contig):
                stop = min(stop, lengths[rec.contig])
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
                "id": rec.id, "contig": rec.contig, "start": rec.start, "end": stop,
                "svtype": rec.info.get("SVTYPE", "?"), "size": stop - rec.start,
                "af": float(af), "ns": int(rec.info.get("NS", 0)),
                "ncarrier": int(rec.info.get("NCARRIER", 0)),
                "hwe": float(rec.info.get("HWE_P", 1.0)),
                "filter": ";".join(rec.filter.keys()) or "PASS",
                "rescue": rec.info.get("RESCUE", "") or "",
                "sqmax": rec.info.get("SQMAX") if has_sqmax else None,
                "carrier_frac": rec.info.get("CARRIER_FRAC") if has_carrier_frac else None,
                "ref_carrier_frac": rec.info.get("REF_CARRIER_FRAC") if has_ref_carrier_frac else None,
            })
            nalt_rows.append(nalt)
            gq_rows.append(gq)
            if with_cn:
                cn_rows.append(cn)
    nalt_m = np.vstack(nalt_rows) if nalt_rows else np.zeros((0, n), np.int8)
    gq_m = np.vstack(gq_rows) if gq_rows else np.zeros((0, n), np.int16)
    cn_m = (
        np.vstack(cn_rows) if cn_rows else np.zeros((0, n), np.int8)
    ) if with_cn else None
    return samples, sites, nalt_m, gq_m, cn_m


def _read_cohort_vcf(path: Path, with_cn: bool = False, jobs: int = 1):
    """One dict per site + compact genotype matrices.

    Dict-of-dicts genotypes cost ~80 bytes/genotype and OOM-killed the report
    on 238k sites x 2670 samples (635M genotypes); nalt/gq live in two numpy
    matrices (int8+int16 = 3 bytes/genotype, ~1.9 GB at that scale) and each
    site carries its row index. ``with_cn`` adds an int8 copy-number matrix
    (family events pages show GT:GQ:CN; 1 byte/genotype more).

    ``jobs`` > 1 (indexed VCF) decodes position-ordered regions in parallel
    — one core spent ~1.5 h in this loop on 2670 samples — and stacks them
    in file order.
    """
    from .parallel import pmap, split_regions

    regions = split_regions(path, jobs)
    parts = list(pmap(_read_cohort_region, [(str(path), r, with_cn) for r in regions], jobs))
    samples = parts[0][0]
    sites = [site for part in parts for site in part[1]]
    for i, site in enumerate(sites):
        site["idx"] = i
    nalt_m = np.vstack([p[2] for p in parts])
    gq_m = np.vstack([p[3] for p in parts])
    cn_m = np.vstack([p[4] for p in parts]) if with_cn else None
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


def _sample_karyotype(task) -> dict | None:
    """sex_chromosome_stats + classify for one sample, with the panel of
    normals its store was called with (none: no informative bins, None)."""
    import zarr

    from .extract.store import EvidenceStore
    from .karyotype import classify, sex_chromosome_stats

    store_path, ped_sex = task
    try:
        store = EvidenceStore(store_path)
        pon = store.meta.get("pon")
        if not pon or not Path(pon).exists():
            return None
        st = sex_chromosome_stats(store, zarr.open(pon, mode="r"))
    except (OSError, KeyError, ValueError):
        return None
    if st is None:
        return None
    return {**st, **classify(st, ped_sex)}


def sex_karyotypes(run_dir: Path, samples: list[str], ped_sex: dict[str, str | None],
                   jobs: int = 1) -> dict[str, dict]:
    """{sample: karyotype record} for the samples whose store carries a PoN."""
    from .parallel import pmap

    tasks = [(str(run_dir / "samples" / s / f"{s}.evidence.zarr"), ped_sex.get(s)) for s in samples]
    return {s: k for s, k in zip(samples, pmap(_sample_karyotype, tasks, jobs), strict=True) if k}


STORE_SECTIONS = ("« Chromosomes sexuels » (caryotype X/Y)", "score de vague par génome "
                  "(colonne QC, alertes des lignes famille)", "bornes des lignes ajustées sur la profondeur")


def check_evidence_stores(run_dir: Path, samples: list[str], require: bool = True) -> list[str]:
    """The samples without an evidence store under samples/<S>/ of the run.

    Karyotypes, wave QC and depth-based row bounds read these stores and
    degrade to nothing when absent: a re-filter run (--calls_dir) published
    only the filter outputs, and the 2670-sample report of 2026-11-02 lost
    its sex-chromosome section without a word. Some missing: a warning;
    all missing: an error unless ``require`` is False."""
    missing = [s for s in samples
               if not (Path(run_dir) / "samples" / s / f"{s}.evidence.zarr").exists()]
    if not missing:
        return missing
    where = f"{Path(run_dir) / 'samples'}/<S>/<S>.evidence.zarr"
    if require and len(missing) == len(samples):
        raise FileNotFoundError(
            f"aucun store de preuves sous {where} ({len(samples)} échantillons) — sections "
            f"perdues : {'; '.join(STORE_SECTIONS)}. Liez les stores des appels "
            "(samples/<S>/ du run d'appel) ou relancez avec --allow-missing-stores"
        )
    logger.warning(
        "%d/%d échantillon(s) sans store de preuves sous %s — absents de : %s. Exemples : %s",
        len(missing), len(samples), where, "; ".join(STORE_SECTIONS),
        ", ".join(missing[:10]) + (" …" if len(missing) > 10 else ""),
    )
    return missing


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
    mosaic_dir: str | Path | None = None,
    subthreshold: bool = True,
    low_quality_prefixes: tuple[str, ...] = (),
    jobs: int = 1,
    require_stores: bool = True,
) -> Path:
    run_dir = Path(run_dir)
    cohort_dir = run_dir / "cohort"
    vcf_path = cohort_dir / "cohort.vcf.gz"
    if not vcf_path.exists():
        raise FileNotFoundError(f"{vcf_path} absent — lancez le pipeline cohorte d'abord")
    # before the hours of VCF decoding: the header has the sample list
    with pysam.VariantFile(str(vcf_path)) as _hdr:
        check_evidence_stores(run_dir, list(_hdr.header.samples), require=require_stores)
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

    def _impact(contig: str, start: int, end: int, svtype: str = "DEL") -> tuple[float, str]:
        return gene_idx.impact(contig, start, end, svtype) if gene_idx else (0.0, "")

    samples, sites, nalt_m, gq_m, cn_m = _read_cohort_vcf(vcf_path, with_cn=family_events, jobs=jobs)
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

    dn_families = denovo_family_counts(denovo_by_family)
    # ---- de novo per-child guard: technical outliers leave the family tables
    dn_counts = denovo_counts(denovo_by_family)
    dn_cut = denovo_cutoff(dn_counts, max_denovo_per_child)
    dn_outliers = {c: n for c, n in dn_counts.items() if n > dn_cut}
    if dn_outliers:
        logger.warning(
            "garde-fou de novo: %d enfant(s) au-dessus de %.0f candidats exclu(s): %s",
            len(dn_outliers), dn_cut,
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

        # de novo table: each child's sites take its own bounds (own_span)
        site_idx = _site_index(sites)
        _own_cache: dict[str, dict] = {}
        dn_depth = _DepthCache(run_dir)

        def _own_idx_of(child: str) -> dict:
            if child not in _own_cache:
                raw = own_call_index(run_dir, [child])
                _own_cache[child] = {
                    k: [(op, oe, ok, who,
                         _cohort_site_pass(site_idx, k[0], op, oe, k[1]) if ok else None)
                        for op, oe, ok, who in v]
                    for k, v in raw.items()
                }
            return _own_cache[child]

        for fid in sorted(families):
            members = families[fid]
            _own_cache.clear()  # a child belongs to one family: keep one family's calls
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

            def _num_max(vals):
                nums = [float(v) for v in vals if v not in ("", None, ".")]
                return f"{max(nums):g}" if nums else ""

            dn_rows = []
            dn_depth.clear()
            for ev in collapse_denovo(denovo, _own_idx_of, site_idx):
                span = refine_bounds(dn_depth, ev["contig"], ev["pos"], ev["end"], ev["svtype"],
                                     [ev["child"]], [m for m in members if m != ev["child"]])
                if span is not None:
                    ev["pos"], ev["end"] = span
                ids = [r["site_id"] for r in ev["sites"]]
                size = ev["end"] - ev["pos"] + 1
                dn_rows.append([
                    ids[0] + (f" (+{len(ids) - 1} sites)" if len(ids) > 1 else ""),
                    f'{ev["contig"]}:{ev["pos"]}-{ev["end"]}', ev["svtype"], size,
                    _genes(ev["contig"], ev["pos"] - 1, ev["end"]),
                    _impact(ev["contig"], ev["pos"] - 1, ev["end"], ev["svtype"])[1],
                    ev["child"], _num_max(r["child_gq"] for r in ev["sites"]),
                    _num_max(r["dnq"] for r in ev["sites"]),
                    ", ".join(denovo_alerts(size, max(dn_families.get(i, 0) for i in ids))) or "—",
                ])
            dn_rows.sort(key=lambda row: -_impact_points(row[5]))
            family_sections.append(f"""
<details><summary>Famille {html.escape(fid)} — {len(members)} membres,
{len(dn_rows)} de novo (PASS, GQ≥{min_denovo_gq}; {len(denovo)} sites), {n_viol} violations</summary>
<p>Membres : {", ".join(map(html.escape, members))} —
workspace : <code>families/{html.escape(fid)}/</code>{_events_link(fid, len(members)) if family_events else ""}</p>
{f'<p class=warn>⚠ {", ".join(map(html.escape, fam_outliers))} : plus de {dn_cut:.0f} candidats de novo ({DENOVO_REL:g}× la médiane de la cohorte) — profil technique atypique, candidats exclus du tableau.</p>' if fam_outliers else ""}
{_table(["site", "locus", "type", "taille", "gènes", "impact", "enfant", "GQ", "DNQ", "alerte"], dn_rows, limit=20)
 if dn_rows else "<p>Aucun candidat de novo au seuil.</p>"}
</details>""")

    # ---- per-sample wave QC (evidence-store meta, absent on old stores)
    wave_scores: dict[str, float | None] = {}
    wave_by_contig: dict[str, dict] = {}
    for s in samples:
        store_path = run_dir / "samples" / s / f"{s}.evidence.zarr"
        wave_scores[s] = None
        if store_path.exists():
            try:
                from .extract.store import EvidenceStore

                meta = EvidenceStore(store_path).meta
                wave_scores[s] = meta.get("wave_score")
                wave_by_contig[s] = meta.get("wave_by_contig") or {}
            except Exception:  # report survives a broken/foreign store
                pass
    n_wavy = sum(1 for v in wave_scores.values() if v is not None and v >= WAVE_WARN)
    sample_qc = sample_qc_flags(overloaded, dn_outliers)
    # ---- mosaic segments (`vedette mosaic --run-dir`), when present
    from .mosaic import read_tsv as _read_mosaic

    mosaic_by_sample: dict[str, list[dict]] = {}
    for s in samples:
        tsv = (Path(mosaic_dir) / f"{s}.mosaic.tsv" if mosaic_dir
               else run_dir / "samples" / s / f"{s}.mosaic.tsv")
        if tsv.exists():
            try:
                mosaic_by_sample[s] = _read_mosaic(tsv)
            except (OSError, ValueError, KeyError):
                pass
    # segments are window-aligned: the last 1 kb bin of a contig runs past its
    # end (chr22 50,819,000 for 50,818,468) — bound them to the real length
    with pysam.VariantFile(str(vcf_path)) as _hdr_vcf:
        _clen = {c.name: c.length for c in _hdr_vcf.header.contigs.values()}
    for segs in mosaic_by_sample.values():
        for g in segs:
            if _clen.get(g["contig"]):
                g["end"] = min(g["end"], _clen[g["contig"]])
    mosaic_rows = sorted(
        ([s, g["contig"], f'{g["start"] / 1e6:.1f}-{g["end"] / 1e6:.1f}', g["svtype"],
          f'{g["cn"]:.2f}', f'{round(100 * g["fraction"])} %',
          f'{g["baf_dev"]:.3f} / {g["baf_expected"]:.3f}']
         for s, segs in mosaic_by_sample.items() for g in segs if g["mosaic"]),
        key=lambda r: (r[0], r[1]),
    )
    n_mosaic_samples = len({r[0] for r in mosaic_rows})
    if mosaic_rows:
        mosaic_section = (
            "<p>Segments ≥ 3 Mb dont le nombre de copies est fractionnaire (perte CN = 2 − f, "
            "gain CN = 2 + f, f = fraction de cellules), confirmés par le déséquilibre allélique "
            "(|BAF − 0,5| observé au-delà du bruit / attendu pour f). Un clone présent dans une "
            "partie des cellules n'est pas constitutionnel : artefact de culture (lignée "
            "lymphoblastoïde) ou processus clonal dans le prélèvement — à vérifier avant toute "
            "interprétation de ses CNV.</p>"
            + _table(["échantillon", "chr", "Mb", "type", "CN", "fraction", "BAF obs. / att."],
                     mosaic_rows)
        )
    elif mosaic_by_sample:
        mosaic_section = "<p>Aucun segment mosaïque.</p>"
    else:
        mosaic_section = "<p>Scan non lancé (<code>vedette mosaic --run-dir</code>).</p>"
    # whole-copy segments of the same scan: aneuploidies and large constitutional
    # gains/losses, which the integer HMM calls in pieces the filter may reject
    # (a trisomy 21 and a trisomy 13 were nowhere in the reports, 2026-10-15)
    constit_rows = sorted(
        ([s, g["contig"], f'{g["start"] / 1e6:.1f}-{g["end"] / 1e6:.1f}', g["svtype"],
          f'{g["cn"]:.2f}', f'{g["baf_dev"]:.3f} / {g["baf_expected"]:.3f}']
         for s, segs in mosaic_by_sample.items() for g in segs if not g["mosaic"]),
        key=lambda r: (r[0], r[1]),
    )
    if constit_rows:
        mosaic_section += (
            "<h3>Anomalies chromosomiques constitutionnelles</h3>"
            "<p>Segments ≥ 3 Mb à un nombre de copies entier (aneuploïdie, gain ou perte "
            "d'un bras ou d'une grande région), confirmés par le déséquilibre allélique "
            "(|BAF − 0,5| observé au-delà du bruit / attendu). Autosomes seulement. Constat "
            "technique à confirmer par une méthode clinique (caryotype, ACPA).</p>"
            + _table(["échantillon", "chr", "Mb", "type", "CN", "BAF obs. / att."], constit_rows)
        )
    # sex chromosomes: the mosaic scan is autosomal, and the caller's own sex
    # inference sets the X/Y baseline — an extra or missing X otherwise shows,
    # at best, as fragmented gains (an XXX genome of dataset_paper)
    ped_sex = ({i.sample_id: {1: "M", 2: "F"}.get(i.sex) for i in pedigree.individuals.values()}
               if pedigree is not None else {})
    karyo = sex_karyotypes(run_dir, samples, ped_sex, jobs=jobs)
    if karyo:
        kc = Counter(k["karyotype"] for k in karyo.values())
        kar_rows = sorted(
            [s_, {"M": "masculin", "F": "féminin"}.get(ped_sex.get(s_) or "", "—"),
             f'{k["x_cn"]:.2f}', "—" if k["y_cn"] is None else f'{k["y_cn"]:.2f}',
             "—" if k["x_het"] is None else f'{k["x_het"]:.2f}',
             "—" if k["x_baf"] is None else f'{k["x_baf"]:.3f}', k["karyotype"], k["note"] or "—"]
            for s_, k in karyo.items() if is_anomaly(k)
        )
        for s_, k in karyo.items():
            if is_anomaly(k):
                sample_qc.setdefault(s_, []).append(
                    k["karyotype"] if not k["note"] else f'{k["karyotype"]} ({k["note"]})')
        mosaic_section += (
            "<h3>Chromosomes sexuels</h3>"
            f"<p>{len(karyo)} génomes : " + ", ".join(f"{n} {k}" for k, n in kc.most_common())
            + ". Nombre de copies de chrX et chrY sur la profondeur brute (bins où le panel de "
            "normaux distingue les sexes : régions pseudo-autosomiques et répétitions X/Y "
            "exclues), hétérozygotie de chrX rapportée aux autosomes (un seul X : ~0 ; deux : "
            "~1), |BAF − 0,5| des sites hétérozygotes de chrX (XX ~0,06 ; XXX ~0,17). "
            "Signalés : caryotype autre que XX/XY, niveau intermédiaire (mosaïque, perte "
            "partielle du Y), un seul X hétérozygote (contamination ?), sexe du pedigree "
            "discordant. Constat technique à confirmer par une méthode clinique.</p>"
            + (_table(["échantillon", "sexe pedigree", "CN chrX", "CN chrY", "hétéroz. X",
                       "BAF X", "caryotype", "remarque"], kar_rows)
               if kar_rows else "<p>Aucune anomalie.</p>")
        )

    # ---- family events pages: one sweep of the matrices, never a per-family
    # rescan of the cohort VCF (the 36 h trap on 1207-family cohorts)
    if family_events and pedigree is not None and inh_spill is not None:
        inh_spill.flush()
        sub_cands = subthreshold_candidates(run_dir, samples, jobs=jobs) if subthreshold else None
        if sub_cands is not None:
            logger.info("events: %d appels sous le seuil candidats (%d échantillons)",
                        len(sub_cands), len(sub_cands.samples))
        n_pages = _write_all_family_events(
            run_dir, dict(families), samples, sites, nalt_m, gq_m, cn_m,
            _genes, _impact, inh_spill, pass_only=events_pass_only,
            segdup_iv=segdup_iv, lcr_iv=lcr_iv, collapse=collapse_family_events,
            min_carrier_gq=min_carrier_gq, sample_qc=sample_qc,
            wave_scores=wave_scores, wave_by_contig=wave_by_contig,
            mosaic_by_sample=mosaic_by_sample, subthreshold=sub_cands,
            denovo_families=dn_families, low_quality_prefixes=tuple(low_quality_prefixes),
            jobs=jobs,
        )
        inh_spill.cleanup()
        logger.info("events: %d pages familles écrites", n_pages)


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
            score, label = _impact(s["contig"], s["start"], s["end"], s["svtype"])
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
<span class=kpi><b>{n_mosaic_samples if mosaic_by_sample else "—"}</b>échantillons avec mosaïque chromosomique</span>
<span class=kpi><b>{len(overloaded)}</b>échantillons surchargés (≥{CARRIER_EXCESS_WARN}x la médiane de porteurs PASS)</span>
<span class=kpi><b>{len(dn_outliers)}</b>enfants outliers de novo (&gt;{dn_cut:.0f})</span>
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
<p>Impact : points de contenu génique d'après la grille ACMG/ClinGen des CNV
(Riggs et al. 2020, sections 1-3), GeneTrek tenant lieu de carte de dosage
ClinGen — ce n'est pas une classification ACMG (hérédité et phénotype non
comptés). DEL : gène NDD établi à transmission dominante délété en entier
+1,00, atteinte exonique +0,90 ; mode inconnu +0,45 ; récessif seul +0,15 ;
gène NDD candidat +0,30 (+0,45 si contraint) ; gène contraint seul
(s_het ≥ 0,1 ou LOEUF &lt; 0,35) +0,15 ; ≥ 25 / ≥ 35 gènes codants +0,45 / +0,90.
DUP : triplosensibilité non disponible — duplication intragénique d'un gène
NDD dominant +0,45 (perte de fonction possible), duplication complète +0,15 ;
≥ 35 / ≥ 50 gènes codants +0,45 / +0,90. Aucun gène codant −0,60. Niveaux :
élevé ≥ 0,90, modéré ≥ 0,45. Tri par impact puis taille.</p>
<h3>DEL ({sum(1 for s in rare if s["svtype"] == "DEL")})</h3>
{_table(["site", "locus", "taille", "gènes", "impact", "AF", "porteurs"], _rare_rows("DEL"), limit=40)}
<h3>DUP ({sum(1 for s in rare if s["svtype"] == "DUP")})</h3>
{_table(["site", "locus", "taille", "gènes", "impact", "AF", "porteurs"], _rare_rows("DUP"), limit=40)}

<h2>Mosaïques chromosomiques</h2>
{mosaic_section}

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
