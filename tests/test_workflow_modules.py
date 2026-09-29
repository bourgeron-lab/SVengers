"""Static checks on the Nextflow process modules.

These catch a class of bug that only shows up on a real cohort run: a process
that passes one of its inputs straight through to its outputs, while a
publishDir directive republishes everything. Nextflow then links the published
entry onto itself and the input becomes an unreadable symlink cycle -- which is
what happened to 2071 of the 2670 published evidence stores of dataset_paper.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

MODULE_DIR = Path(__file__).resolve().parents[1] / "workflows" / "nextflow" / "modules"

# `path(store)` / `path store`: a bare identifier, as opposed to `path("x.vcf")`
# or `path("${sample}.gt.tsv")`, which name a file the process creates itself.
_BARE_PATH = re.compile(r"\bpath[(\s]\s*([A-Za-z_][A-Za-z0-9_]*)\s*[),]")


def _modules() -> list[Path]:
    return sorted(MODULE_DIR.glob("*.nf"))


def _section(text: str, name: str) -> str:
    """The body of an `input:` / `output:` block, up to the next block."""
    m = re.search(rf"^\s*{name}:\s*$(.*?)^\s*(?:input|output|script|exec|when):",
                  text, re.M | re.S)
    return m.group(1) if m else ""


def test_modules_exist() -> None:
    assert _modules(), f"no Nextflow module found under {MODULE_DIR}"


@pytest.mark.parametrize("module", _modules(), ids=lambda p: p.name)
def test_passthrough_inputs_are_not_published(module: Path) -> None:
    text = module.read_text()
    if "publishDir" not in text:
        return

    inputs = set(_BARE_PATH.findall(_section(text, "input")))
    outputs = set(_BARE_PATH.findall(_section(text, "output")))
    passthrough = inputs & outputs
    if not passthrough:
        return

    publish = re.search(r"publishDir(.*?)(?=^\s*(?:input|output|tag|label|cpus|memory):)",
                        text, re.M | re.S)
    directive = publish.group(1) if publish else ""
    assert "saveAs" in directive, (
        f"{module.name} re-emits its input(s) {sorted(passthrough)} and publishes "
        "every output: publishDir would link the published entry onto itself. "
        "Exclude the passthrough with saveAs."
    )
