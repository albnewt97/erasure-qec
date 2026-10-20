"""README result tables equal what scripts/render_readme_tables.py renders from results/*.json.

For each table name, README.md may carry ``<!-- BEGIN:name -->`` ...
``<!-- END:name -->``; the text between the marker lines must equal the rendered
table exactly. A table whose JSON or markers are missing is skipped, naming why.
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path
from types import ModuleType

import pytest

_ROOT = Path(__file__).resolve().parents[1]
README = _ROOT / "README.md"


def _load(name: str) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, _ROOT / "scripts" / f"{name}.py")
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


render_readme_tables = _load("render_readme_tables")


def _readme_block(name: str) -> str | None:
    """Text between the BEGIN/END marker lines of ``name`` in README.md, or None."""
    text = README.read_text()
    begin, end = f"<!-- BEGIN:{name} -->\n", f"\n<!-- END:{name} -->"
    if begin not in text or end not in text:
        return None
    start = text.index(begin) + len(begin)
    return text[start : text.index(end, start)]


@pytest.mark.parametrize("name", list(render_readme_tables.TABLES))
def test_readme_table_matches_rendered(name: str) -> None:
    source = render_readme_tables.source(name)
    rendered = render_readme_tables.render(name)
    if rendered is None:
        pytest.skip(f"{source.relative_to(_ROOT)} does not exist")
    block = _readme_block(name)
    if block is None:
        pytest.skip(f"README.md has no <!-- BEGIN:{name} --> / <!-- END:{name} --> markers")
    assert block == rendered


def test_every_table_ends_with_the_numpy_footnote(tmp_path: Path) -> None:
    """Rendering works for every table, including ablation from a stub payload."""
    ablation_table = _load("ablation_table")
    shots = 10_000
    cells = {
        (p, r_e, d): {
            "h_err": h,
            "b_err": 3 * h,
            **ablation_table.ratio_and_bound(h, 3 * h, shots, d),
        }
        for p in (0.005, 0.010)
        for r_e in ablation_table.R_ES
        for d, h in zip(ablation_table.DISTANCES, (400, 200, 100, 40, 0), strict=True)
    }
    doc = ablation_table.payload(cells, (0.005, 0.010), shots, "stub")
    (tmp_path / "ablation.json").write_text(json.dumps(doc, allow_nan=False))
    for name in ("headline.json", "wu_validation.json"):
        src = _ROOT / "results" / name
        if src.exists():
            (tmp_path / name).write_text(src.read_text())

    ablation = render_readme_tables.render("ablation", tmp_path)
    assert ablation is not None
    rows = [line for line in ablation.splitlines() if line.startswith(("| 0.50 |", "| 1.00 |"))]
    assert len(rows) == len(cells)
    assert "> " in ablation and "†" in ablation  # gated cells: lower bound only
    for name in render_readme_tables.TABLES:
        table = render_readme_tables.render(name, tmp_path)
        if table is not None:
            assert table.endswith("\n\n" + render_readme_tables.FOOTNOTE)
