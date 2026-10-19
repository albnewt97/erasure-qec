"""Classification helpers of ``scripts/wu_validation.py``: edge rule, status, Wu-in-CI."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

from erasure_qec.analysis.statistics import SweepPoint

_SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "wu_validation.py"
_spec = importlib.util.spec_from_file_location("wu_validation", _SCRIPT)
assert _spec is not None and _spec.loader is not None
wu_validation = importlib.util.module_from_spec(_spec)
sys.modules["wu_validation"] = wu_validation
_spec.loader.exec_module(wu_validation)

GRID = [0.01, 0.02, 0.03, 0.04]


def test_edge_rule_bracket_at_either_end() -> None:
    assert wu_validation.edge_touched((0.01, 0.02), [-1.0, 1.0], GRID) == "lower"
    assert wu_validation.edge_touched((0.03, 0.04), [-1.0, 1.0], GRID) == "upper"
    assert wu_validation.edge_touched((0.02, 0.03), [-1.0, 1.0], GRID) == "none"


def test_edge_rule_without_crossing() -> None:
    assert wu_validation.edge_touched(None, [-0.3, -0.2, -0.1], GRID) == "upper"
    assert wu_validation.edge_touched(None, [0.1, 0.2, 0.3], GRID) == "lower"
    assert wu_validation.edge_touched(None, [0.1, -0.2, -0.3], GRID) == "lower"


def test_no_crossing_status() -> None:
    assert wu_validation.no_crossing_status([0.0, 0.2]) == "no crossing in grid (below)"
    assert wu_validation.no_crossing_status([-0.2, -0.1]) == "no crossing in grid (above)"


def test_sign_changes_count_both_directions() -> None:
    assert wu_validation.sign_changes([-1.0, 0.0, -0.5, 0.5]) == 3
    assert wu_validation.sign_changes([]) == 0


def test_wu_inside_is_closed_and_null_safe() -> None:
    assert wu_validation.wu_inside((0.01, 0.02), 0.02) is True
    assert wu_validation.wu_inside((0.01, 0.02), 0.021) is False
    assert wu_validation.wu_inside((0.01, 0.02), None) is None
    assert wu_validation.wu_inside((None, None), 0.015) is None


def _pt(d: int, p: float, errors: int) -> SweepPoint:
    return SweepPoint("wu_herald_mwpm", d, d, p, 0.0, 1000, errors)


def test_unsaturated_flag_uses_the_cap_at_both_distances() -> None:
    points = [_pt(9, 0.01, 100), _pt(13, 0.01, 50), _pt(9, 0.02, 300), _pt(13, 0.02, 450)]
    g = wu_validation.usable_g(points, 9, 13)
    assert [(p, unsat) for p, _, unsat in g] == [(0.01, True), (0.02, False)]
    assert g[0][1] < 0.0 < g[1][1]


def test_unsaturated_crossings_needs_both_bracket_points_unsaturated() -> None:
    g = [(0.01, -0.2, True), (0.02, 0.1, True), (0.03, -0.1, True), (0.04, 0.3, False)]
    assert wu_validation.unsaturated_crossings(g) == 1
    assert wu_validation.unsaturated_crossings([(0.01, -0.2, False), (0.02, 0.1, True)]) == 0
