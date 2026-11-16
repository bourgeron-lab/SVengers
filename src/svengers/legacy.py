"""Names written before v2.0, when SVengers was called Vedette.

Outputs now carry the new name; every reader also accepts what Vedette wrote,
so the genomes called before the rename stay usable without re-reading their
CRAMs (evidence-store contract). The evidence store itself never carried the
tool name: nothing to translate there.
"""

from __future__ import annotations

import sys
from collections.abc import Iterator, Mapping
from contextlib import contextmanager
from pathlib import Path

NAME = "svengers"
LEGACY_NAME = "vedette"
NAMES = (NAME, LEGACY_NAME)

# per-sample raw calls of `call`: <S>.svengers.vcf.gz (<S>.vedette.vcf.gz before 2.0)
RAW_CALLS_SUFFIXES = tuple(f".{n}.vcf.gz" for n in NAMES)


def raw_calls_vcf(directory: str | Path, sample: str) -> Path:
    """A sample's raw calls VCF in ``directory``: ``<S>.svengers.vcf.gz``, else
    ``<S>.vedette.vcf.gz`` when only that one exists (the new name when
    neither does, for the error message)."""
    return existing_raw_calls(Path(directory) / f"{sample}{RAW_CALLS_SUFFIXES[0]}")


def existing_raw_calls(path: str | Path) -> Path:
    """``path`` itself, unless it is a missing ``<S>.svengers.vcf.gz`` whose
    pre-2.0 sibling ``<S>.vedette.vcf.gz`` exists: then the sibling."""
    path = Path(path)
    new, old = RAW_CALLS_SUFFIXES
    if not path.exists() and path.name.endswith(new):
        legacy = path.with_name(path.name[: -len(new)] + old)
        if legacy.exists():
            return legacy
    return path


def caller_version(meta: Mapping) -> str | None:
    """Caller version recorded in an export or model manifest
    (``svengers_version``, ``vedette_version`` before 2.0)."""
    for key in (f"{n}_version" for n in NAMES):
        if meta.get(key):
            return str(meta[key])
    return None



@contextmanager
def legacy_modules() -> Iterator[None]:
    """While unpickling, ``vedette.<module>`` is the imported ``svengers.<module>``.

    The filter models trained before 2.0 are pickles that name their classes
    ``vedette.filter.model.CalibratedModel``; the aliases are removed (and any
    previous ``vedette.*`` entry restored) on exit."""
    absent = object()
    saved: dict[str, object] = {}
    for name, module in list(sys.modules.items()):
        if name == NAME or name.startswith(NAME + "."):
            old = LEGACY_NAME + name[len(NAME):]
            saved[old] = sys.modules.get(old, absent)
            sys.modules[old] = module
    try:
        yield
    finally:
        for old, module in saved.items():
            if module is absent:
                sys.modules.pop(old, None)
            else:
                sys.modules[old] = module
