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
import json
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


def test_lambda_is_gated_at_both_distances() -> None:
    assert ablation_table.lambda_d(49, 500, 100_000, 5) is None
    assert ablation_table.lambda_d(500, 49, 100_000, 5) is None
    lam = ablation_table.lambda_d(500, 50, 100_000, 5)
    expected = per_round_p_l(500 / 100_000, 5) / per_round_p_l(50 / 100_000, 7)
    assert lam == pytest.approx(expected)
    assert lam is not None and lam > 1.0


def _stub_cells(shots: int) -> dict[tuple[float, float, int], dict[str, float]]:
    """Cells built from counts (no decoding): herald errors fall with d, one cell at 0."""
    counts = {3: (900, 2000), 5: (300, 1500), 7: (60, 1200), 9: (10, 1000), 11: (0, 800)}
    return {
        (p, r_e, d): {"h_err": h, "b_err": b, **ablation_table.ratio_and_bound(h, b, shots, d)}
        for p in (0.005, 0.010)
        for r_e in ablation_table.R_ES
        for d, (h, b) in counts.items()
    }


def test_cell_rows_carry_each_decoders_lambda() -> None:
    shots = 100_000
    rows = {
        (r["p"], r["r_e"], r["d"]): r for r in ablation_table.cell_rows(_stub_cells(shots), shots)
    }
    assert len(rows) == 2 * len(ablation_table.R_ES) * len(ablation_table.DISTANCES)
    r3 = rows[(0.005, 0.5, 3)]
    assert r3["lambda_pair"] == [3, 5]
    assert r3["lambda_herald"] == ablation_table.lambda_d(900, 300, shots, 3)
    assert r3["lambda_blind"] == ablation_table.lambda_d(2000, 1500, shots, 3)
    assert rows[(0.005, 0.5, 7)]["lambda_herald"] is None  # herald d = 9 has 10 errors
    assert rows[(0.005, 0.5, 7)]["lambda_blind"] is not None
    r11 = rows[(0.010, 0.98, 11)]
    assert r11["lambda_pair"] is None and r11["lambda_herald"] is None
    assert r11["gated"] and r11["ratio"] is None  # h_err = 0: infinite ratio -> null


def test_payload_is_strict_json_with_provenance() -> None:
    shots = 100_000
    doc = ablation_table.payload(_stub_cells(shots), (0.005, 0.010), shots, "abc1234")
    text = json.dumps(doc, allow_nan=False)
    back = json.loads(text)
    assert back["git_describe"] == "abc1234"
    assert back["parameters"]["p_values"] == [0.005, 0.010]
    assert back["parameters"]["seed"] == ablation_table.SEED == 0
    assert back["parameters"]["shots"] == shots
    assert set(back["software"]) >= {"stim", "pymatching", "numpy"}


def test_printed_table_covers_each_p(capsys: pytest.CaptureFixture[str]) -> None:
    shots = 100_000
    cells = _stub_cells(shots)
    for p in (0.005, 0.010):
        ablation_table._print_table(cells, p, shots)
    out = capsys.readouterr().out
    assert "at p=0.5%" in out and "at p=1.0%" in out
    assert "> " in out  # gated cells print a lower bound, not a point ratio
    assert "d=3->5:" in out
