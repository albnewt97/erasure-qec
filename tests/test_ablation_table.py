"""The suppression-table bound in ``scripts/ablation_table.py`` is a real bound.

The bound used to be ``b_pl / h_pl_hi`` -- blind's *point estimate* over herald's
Wilson upper limit -- which the README quoted as "a Wilson lower bound `> N×`".
It was not one: nothing bounded the numerator from below, so the quoted value
could exceed the true ratio whenever blind's rate was over-estimated. It is now
``b_pl_lo / h_pl_hi``. These tests pin that the replacement is genuinely no
larger than both the superseded expression and the point ratio.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

from erasure_qec.analysis.statistics import per_round_p_l, wilson_interval

_SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "ablation_table.py"
_spec = importlib.util.spec_from_file_location("ablation_table", _SCRIPT)
assert _spec is not None and _spec.loader is not None
ablation_table = importlib.util.module_from_spec(_spec)
sys.modules["ablation_table"] = ablation_table
_spec.loader.exec_module(ablation_table)


def _superseded_bound(h_err: int, b_err: int, shots: int, d: int) -> float:
    """The pre-fix expression: blind's POINT estimate over herald's Wilson upper."""
    h_pl_hi = per_round_p_l(wilson_interval(h_err, shots)[1], d)
    b_pl = per_round_p_l(b_err / shots, d)
    return b_pl / h_pl_hi if h_pl_hi > 0 else float("inf")


# (h_err, b_err, shots, d): the shapes the real table hits -- well-measured
# cells, the single-digit d=9/d=11 r_e=0.98 cells, and the zero-error edge.
_COUNTS = [
    (h, b, shots, d)
    for shots in (2_000, 100_000)
    for d in (3, 5, 7, 9, 11)
    for h, b in [
        (0, 0),
        (0, 7),
        (1, 40),
        (3, 900),
        (23, 2400),
        (50, 3000),
        (400, 1500),
        (2000, 2600),
        (5000, 5000),
        (6000, 4000),
    ]
    if h <= shots and b <= shots
]


@pytest.mark.parametrize(("h_err", "b_err", "shots", "d"), _COUNTS)
def test_bound_is_below_point_ratio_and_below_superseded(
    h_err: int, b_err: int, shots: int, d: int
) -> None:
    """``ratio_lb <= superseded bound <= ratio`` for every cell shape.

    The middle term is the old, invalid expression: the new bound must be no
    larger (it only shrinks the numerator), and the old one was already below
    the point ratio (it only inflated the denominator). Chaining them is what
    shows the fix is a strict tightening in the safe direction, not a rescale.
    """
    cell = ablation_table.ratio_and_bound(h_err, b_err, shots, d)
    old = _superseded_bound(h_err, b_err, shots, d)

    assert cell["ratio_lb"] <= old + 1e-12
    assert old <= cell["ratio"] + 1e-12
    assert cell["ratio_lb"] <= cell["ratio"] + 1e-12
    # And the ingredients bracket the point estimates the right way round.
    assert cell["b_pl_lo"] <= cell["b_pl"]
    assert cell["h_pl_hi"] >= cell["h_pl"]


def test_fix_actually_moves_the_number() -> None:
    """Guard against a no-op 'fix': with b_err in the thousands the numerator's
    Wilson lower limit is a few percent below its point estimate, so the bound
    must be strictly smaller than the superseded one -- not merely <=."""
    h_err, b_err, shots, d = 3, 2400, 100_000, 11
    cell = ablation_table.ratio_and_bound(h_err, b_err, shots, d)
    old = _superseded_bound(h_err, b_err, shots, d)
    assert cell["ratio_lb"] < old
    # A few percent, not an order of magnitude: b_err is well measured.
    assert 0.90 < cell["ratio_lb"] / old < 1.0


def test_ablation_cell_end_to_end_respects_the_invariants() -> None:
    """The same invariants hold on a real (small) decoded cell, so the bound
    survives the trip through ``ablation_cell``'s simulation and dict merge.

    ``ablation_cell`` samples from a seeded stim circuit, and a seeded stim draw
    reproduces only on the same machine and stim version (SIMD width changes it).
    So the count assertions below are a band, not an equality. Measured
    2026-09-06 on arm64 / stim 1.16.0 / pymatching 2.4.0: h_err 405, b_err 573,
    margin 168. The same (d, p, r_e) at n=2048 in
    ``tests/test_herald_decoder.py`` gives a margin sd of 21.7 (McNemar, shared
    shots) and a 12-seed margin spread of 15.2, so ``> 60`` is roughly 7 sd below
    the measured margin -- portable, but still fails if herald-conditioning stops
    helping on a real decoded cell.
    """
    cell = ablation_table.ablation_cell(d=3, p=0.03, r_e=0.9, shots=2_000)
    assert cell["h_err"] > 0
    assert cell["b_err"] - cell["h_err"] > 60
    old = _superseded_bound(int(cell["h_err"]), int(cell["b_err"]), 2_000, 3)
    assert cell["ratio_lb"] <= old <= cell["ratio"]
    assert cell["ratio_lb"] < cell["ratio"]
