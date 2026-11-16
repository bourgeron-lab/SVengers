"""Curators' validations from genetics-viz as labels (issue #11).

genetics-viz (the GHFC curation tool) appends each SVengers validation to
``<data_dir>/validations/vedette.tsv``::

    FID  Variant  Sample  User  Inheritance  Validation  CuratedStart  CuratedEnd
    Comment  Ignore  Timestamp  SiteId  CallerVersion  ModelVersion

- ``Variant`` = ``chr:start-end:del|dup``, the SVengers site's bounds, read as
  1-based inclusive (the family export's convention); matching is by overlap,
  so a base either way changes nothing;
- ``Validation``: ``present``, ``absent``, ``uncertain``, ``different`` (wrong
  type or breakpoints);
- ``CuratedStart``/``CuratedEnd``: breakpoints fixed by the curator, same
  convention;
- append-only: the last row of a (FID, Variant, Sample) not marked
  ``Ignore = 1`` wins.

Labels follow ``events.html``'s: ``present`` = vrai (1), ``absent`` = faux
(0), ``uncertain`` = incertain. ``different`` has no counterpart there: the
call is neither the event nor absent, so it is no training label.
"""

from __future__ import annotations

import csv
import logging
import re
from dataclasses import dataclass
from pathlib import Path

from .label_trio import _overlapping

logger = logging.getLogger(__name__)

LABELS = {"present": 1, "absent": 0, "uncertain": None, "different": None}
_VARIANT = re.compile(r"^(?P<contig>[^:]+):(?P<start>[\d,]+)-(?P<end>[\d,]+):(?P<svtype>del|dup)$", re.I)

LABEL_COLUMNS = [
    "fid", "sample", "contig", "start", "end", "svtype", "validation", "label",
    "curated_start", "curated_end", "start_offset", "end_offset",
    "inheritance", "user", "timestamp", "site_id", "caller_version", "model_version", "comment",
]


@dataclass
class Validation:
    fid: str
    sample: str
    contig: str
    start: int            # 0-based half-open, as everywhere inside SVengers
    end: int
    svtype: str
    validation: str
    curated_start: int | None
    curated_end: int | None
    inheritance: str
    user: str
    timestamp: str
    site_id: str
    caller_version: str
    model_version: str
    comment: str

    @property
    def label(self) -> int | None:
        return LABELS[self.validation]

    def row(self) -> dict:
        def off(a, b):
            return "" if a is None else a - b
        return {
            "fid": self.fid, "sample": self.sample, "contig": self.contig,
            "start": self.start, "end": self.end, "svtype": self.svtype,
            "validation": self.validation, "label": "" if self.label is None else self.label,
            "curated_start": "" if self.curated_start is None else self.curated_start,
            "curated_end": "" if self.curated_end is None else self.curated_end,
            "start_offset": off(self.curated_start, self.start),
            "end_offset": off(self.curated_end, self.end),
            "inheritance": self.inheritance, "user": self.user, "timestamp": self.timestamp,
            "site_id": self.site_id, "caller_version": self.caller_version,
            "model_version": self.model_version, "comment": self.comment,
        }


def _int(v: str | None) -> int | None:
    v = (v or "").strip().replace(",", "")
    return int(v) if v.lstrip("-").isdigit() else None


def read_genetics_viz(path: str | Path) -> tuple[list[Validation], dict[str, int]]:
    """(validations, counts): the current validation of every (FID, Variant,
    Sample) of a genetics-viz ``vedette.tsv``; counts of rows read, ignored,
    superseded and unreadable."""
    counts = {"rows": 0, "ignored": 0, "superseded": 0, "unreadable": 0}
    last: dict[tuple[str, str, str], Validation] = {}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            counts["rows"] += 1
            if (r.get("Ignore") or "").strip() in ("1", "true", "True"):
                counts["ignored"] += 1
                continue
            m = _VARIANT.match((r.get("Variant") or "").strip())
            val = (r.get("Validation") or "").strip().lower()
            if m is None or val not in LABELS:
                counts["unreadable"] += 1
                continue
            cs, ce = _int(r.get("CuratedStart")), _int(r.get("CuratedEnd"))
            v = Validation(
                fid=r.get("FID", ""), sample=r.get("Sample", ""), contig=m["contig"],
                start=_int(m["start"]) - 1, end=_int(m["end"]), svtype=m["svtype"].upper(),
                validation=val,
                curated_start=None if cs is None else cs - 1, curated_end=ce,
                inheritance=r.get("Inheritance", ""), user=r.get("User", ""),
                timestamp=r.get("Timestamp", ""), site_id=r.get("SiteId", "") or "",
                caller_version=r.get("CallerVersion", "") or "",
                model_version=r.get("ModelVersion", "") or "", comment=r.get("Comment", "") or "",
            )
            key = (v.fid, r["Variant"].strip(), v.sample)
            # append-only file: a later row supersedes; the timestamp decides
            # when rows were merged out of order
            if key in last:
                counts["superseded"] += 1
                if v.timestamp < last[key].timestamp:
                    continue
            last[key] = v
    return list(last.values()), counts


def write_labels(validations: list[Validation], out_tsv: str | Path) -> None:
    with open(out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=LABEL_COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for v in sorted(validations, key=lambda v: (v.sample, v.contig, v.start)):
            w.writerow(v.row())


def label_features(features_tsv: str | Path, validations: list[Validation], sample: str,
                   out_tsv: str | Path, min_ro: float = 0.5) -> dict[str, int]:
    """Label one sample's feature rows (``svengers filter --features-out``) with
    its validations: same type, reciprocal overlap >= ``min_ro``. A row
    matching only ``present`` gets 1, only ``absent`` 0; rows matching none,
    an uncertain or different validation, or contradictory ones, are left
    out. Output: the feature columns plus ``label``, as ``svengers label``."""
    sites: dict[tuple[str, str], list] = {}
    for v in validations:
        if v.sample == sample:
            sites.setdefault((v.contig, v.svtype), []).append((v.start, v.end, v.label, None))
    for k in sites:
        sites[k].sort(key=lambda s: (s[0], s[1]))
    counts = {"tp": 0, "fp": 0, "unlabeled": 0, "conflict": 0}
    with open(features_tsv) as fin, open(out_tsv, "w", newline="") as fout:
        reader = csv.DictReader(fin, delimiter="\t")
        writer = csv.DictWriter(fout, fieldnames=[*reader.fieldnames, "label"], delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        for r in reader:
            hits = {lab for lab, _ in _overlapping(sites.get((r["contig"], r["svtype"]), []),
                                                   int(r["start"]), int(r["end"]), min_ro)}
            if None in hits or not hits:
                counts["unlabeled"] += 1
                continue
            if len(hits) > 1:
                counts["conflict"] += 1
                continue
            r["label"] = hits.pop()
            counts["tp" if r["label"] == 1 else "fp"] += 1
            writer.writerow(r)
    logger.info("labels %s (%s): %d TP, %d FP, %d unlabeled, %d conflicting",
                features_tsv, sample, counts["tp"], counts["fp"], counts["unlabeled"],
                counts["conflict"])
    return counts


def breakpoint_offsets(validations: list[Validation]) -> dict[str, list[int]]:
    """{'start': [...], 'end': [...]}: curated minus called breakpoint of the
    validations present with curated bounds (bp; negative: the event starts or
    ends before the call)."""
    out: dict[str, list[int]] = {"start": [], "end": []}
    for v in validations:
        if v.validation != "present":
            continue
        if v.curated_start is not None:
            out["start"].append(v.curated_start - v.start)
        if v.curated_end is not None:
            out["end"].append(v.curated_end - v.end)
    return out
