"""Pairwise-crossing estimator: synthetic ground truth, edge cases, bootstrap."""

import math
import random
from collections.abc import Callable

import pytest

from erasure_qec.analysis.crossing import bootstrap_pairwise_crossing, pairwise_crossing
from erasure_qec.analysis.statistics import SweepPoint, metric_value

TRUE_CROSSING = 0.01
SHOTS = 10**6
DISTANCES = (5, 7, 9, 11)
ADJACENT_PAIRS = ((5, 7), (7, 9), (9, 11))
# a fixes P(0.01, d) = 0.5 * sigmoid(a) = 0.1 for every d; b is the slope.
A = math.log(0.25)
B = 100.0
NU = 1.5


def _expected_shot_rate(p: float, d: int) -> float:
    """P(p, d) = 0.5 * sigmoid(a + b (p - 0.01) d^(1/nu)); always < 1/2."""
    z = A + B * (p - TRUE_CROSSING) * d ** (1.0 / NU)
    return 0.5 / (1.0 + math.exp(-z))


def _grid() -> list[float]:
    return [round(0.004 + 0.002 * i, 3) for i in range(9)]  # 0.004 .. 0.020


def _synthetic_points() -> list[SweepPoint]:
    """Per-shot curves that all pass through P = 0.1 at p = 0.01 (T = d)."""
    return [
        SweepPoint("herald_mwpm", d, d, p, 0.0, SHOTS, round(_expected_shot_rate(p, d) * SHOTS))
        for d in DISTANCES
        for p in _grid()
    ]


def _pair_points(d1: int, d2: int, e1: list[int], e2: list[int]) -> list[SweepPoint]:
    """Hand-built curves on p = 0.01, 0.02, ... with the given error counts."""
    out: list[SweepPoint] = []
    for i, (a, b) in enumerate(zip(e1, e2, strict=True)):
        p = round(0.01 * (i + 1), 2)
        out.append(SweepPoint("herald_mwpm", d1, d1, p, 0.0, 10_000, a))
        out.append(SweepPoint("herald_mwpm", d2, d2, p, 0.0, 10_000, b))
    return out


# --- synthetic ground truth ---------------------------------------------------


@pytest.mark.parametrize(("d1", "d2"), [(5, 7), (7, 9), (9, 11), (5, 9), (5, 11), (7, 11)])
def test_shot_crossing_recovers_true_crossing(d1: int, d2: int) -> None:
    est = pairwise_crossing(_synthetic_points(), d1, d2, metric="shot")
    assert est.value == pytest.approx(TRUE_CROSSING, abs=0.0005)
    assert (est.metric, est.d1, est.d2) == ("shot", d1, d2)


def test_round_negative_wherever_shot_nonpositive() -> None:
    """For d2 > d1 and 0 < P < 1/2: g_shot <= 0 implies g_round < 0.

    With r(P, d) = (1 - (1 - 2P)^(1/d)) / 2 and T = d, the ratio r(P, d) / P is
    increasing in P and decreasing in d on 0 < P < 1/2. So if
    P(d2) <= P(d1) then r(P(d2), d2) <= P(d2) r(P(d1), d2) / P(d1)
    < P(d2) r(P(d1), d1) / P(d1) <= r(P(d1), d1). Hence no per-round crossing
    can lie below the per-shot bracket.
    """
    by_key = {(pt.d, pt.p): pt for pt in _synthetic_points()}
    n_checked = 0
    for d1, d2 in ADJACENT_PAIRS:
        for p in _grid():
            lo, hi = by_key[(d1, p)], by_key[(d2, p)]
            g_shot = math.log(metric_value(hi, "shot")) - math.log(metric_value(lo, "shot"))
            if g_shot <= 0.0:
                g_round = math.log(metric_value(hi, "round")) - math.log(metric_value(lo, "round"))
                assert g_round < 0.0, (d1, d2, p, g_shot, g_round)
                n_checked += 1
    assert n_checked > 0


def _sigmoid(z: float) -> float:
    return 1.0 / (1.0 + math.exp(-z))


def _points_from(rate: Callable[[float, int], float], grid: list[float]) -> list[SweepPoint]:
    """Points with errors = round(P_shot * 10**6) for P_shot = rate(p, d), T = d."""
    return [
        SweepPoint("herald_mwpm", d, d, p, 0.0, SHOTS, round(rate(p, d) * SHOTS))
        for d in DISTANCES
        for p in grid
    ]


def _model_a_shot_rate(p: float, d: int) -> float:
    """Model A: the per-shot rate obeys the scaling form (nu = 1)."""
    return 0.5 * _sigmoid(math.log(0.2 / 0.8) + 100.0 * (p - 0.01) * d ** (1.0 / 1.0))


def _model_b_shot_rate(p: float, d: int) -> float:
    """Model B: the per-round rate obeys the scaling form (nu = 1.5)."""
    r = 0.5 * _sigmoid(math.log(0.02 / 0.98) + 60.0 * (p - 0.01) * d ** (1.0 / 1.5))
    return 0.5 * (1.0 - (1.0 - 2.0 * r) ** d)


def test_shot_collapse_gives_round_crossing_above_and_drifting_down() -> None:
    """Model A: P_shot collapses at p = 0.01, so the per-shot crossings sit there.

    Models A and B together show that the per-round crossing is never below
    the per-shot crossing. Which of the two drifts with d depends on which
    quantity obeys the scaling form; that is an empirical question for the
    data, not a property of the definitions. Here the per-shot rate obeys it,
    so the per-round crossings lie above 0.01 and move down as d grows.
    """
    pts = _points_from(_model_a_shot_rate, [round(0.004 + 0.001 * k, 4) for k in range(27)])
    shot = [pairwise_crossing(pts, d1, d2, metric="shot").value for d1, d2 in ADJACENT_PAIRS]
    rnd = [pairwise_crossing(pts, d1, d2, metric="round").value for d1, d2 in ADJACENT_PAIRS]
    assert shot == [pytest.approx(TRUE_CROSSING, abs=0.0005)] * 3, shot
    assert all(v > 0.0105 for v in rnd), rnd
    assert rnd[0] > rnd[1] > rnd[2], rnd


def test_round_collapse_gives_shot_crossing_below_and_drifting_up() -> None:
    """Model B: the per-round rate collapses at p = 0.01, so the per-round
    crossings sit there.

    Models A and B together show that the per-round crossing is never below
    the per-shot crossing. Which of the two drifts with d depends on which
    quantity obeys the scaling form; that is an empirical question for the
    data, not a property of the definitions. Here the per-round rate obeys it,
    so the per-shot crossings lie below 0.01 and move up as d grows.
    """
    pts = _points_from(_model_b_shot_rate, [round(0.001 + 0.0005 * k, 4) for k in range(39)])
    rnd = [pairwise_crossing(pts, d1, d2, metric="round").value for d1, d2 in ADJACENT_PAIRS]
    shot = [pairwise_crossing(pts, d1, d2, metric="shot").value for d1, d2 in ADJACENT_PAIRS]
    assert rnd == [pytest.approx(TRUE_CROSSING, abs=0.0005)] * 3, rnd
    assert all(v < 0.0095 for v in shot), shot
    assert shot[0] < shot[1] < shot[2], shot


# --- edge cases -----------------------------------------------------------------


def test_no_crossing_returns_nan() -> None:
    est = pairwise_crossing(_pair_points(5, 7, [1000] * 4, [500] * 4), 5, 7, metric="shot")
    assert math.isnan(est.value)
    assert est.bracket is None
    assert est.n_sign_changes == 0


def test_two_sign_changes_reports_first() -> None:
    """g: - + - -> two sign changes; the first (neg -> non-neg) crossing is kept."""
    pts = _pair_points(5, 7, [1000, 1000, 1000], [500, 2000, 500])
    est = pairwise_crossing(pts, 5, 7, metric="shot")
    assert est.n_sign_changes == 2
    assert est.bracket == (0.01, 0.02)
    # g = -ln 2 then +ln 2: linear in ln p puts g = 0 at the geometric mean.
    assert est.value == pytest.approx(math.sqrt(0.01 * 0.02))


def test_three_sign_changes_reports_first() -> None:
    pts = _pair_points(5, 7, [1000] * 4, [500, 2000, 500, 2000])
    est = pairwise_crossing(pts, 5, 7, metric="shot")
    assert est.n_sign_changes == 3
    assert est.bracket == (0.01, 0.02)


def test_zero_error_points_are_skipped() -> None:
    pts = _synthetic_points()
    target = next(pt for pt in pts if pt.d == 7 and pt.p == 0.008)
    zeroed = [
        SweepPoint(pt.decoder, pt.d, pt.rounds, pt.p, pt.r_e, pt.shots, 0) if pt == target else pt
        for pt in pts
    ]
    removed = [pt for pt in pts if pt != target]
    for metric in ("shot", "round"):
        a = pairwise_crossing(zeroed, 5, 7, metric=metric)
        b = pairwise_crossing(removed, 5, 7, metric=metric)
        assert (a.bracket, a.n_sign_changes) == (b.bracket, b.n_sign_changes)
        assert a.value == b.value or (math.isnan(a.value) and math.isnan(b.value))


def test_mixed_conditions_and_duplicates_raise() -> None:
    pts = _pair_points(5, 7, [1000] * 2, [500, 2000])
    mixed = [*pts, SweepPoint("blind_mwpm", 5, 5, 0.05, 0.0, 10_000, 10)]
    with pytest.raises(ValueError):
        pairwise_crossing(mixed, 5, 7, metric="shot")
    duplicated = [*pts, pts[0]]
    with pytest.raises(ValueError):
        pairwise_crossing(duplicated, 5, 7, metric="shot")
    with pytest.raises(ValueError):
        pairwise_crossing(pts, 5, 5, metric="shot")


# --- bootstrap ------------------------------------------------------------------


def test_bootstrap_is_seed_deterministic() -> None:
    pts = _synthetic_points()
    a = bootstrap_pairwise_crossing(pts, 7, 11, metric="shot", n_boot=50, seed=3)
    b = bootstrap_pairwise_crossing(pts, 7, 11, metric="shot", n_boot=50, seed=3)
    assert a == b
    c = bootstrap_pairwise_crossing(pts, 7, 11, metric="shot", n_boot=50, seed=4)
    assert c[1] != a[1]


def test_bootstrap_is_independent_of_input_order() -> None:
    pts = _synthetic_points()
    shuffled = list(pts)
    random.Random(0).shuffle(shuffled)
    assert bootstrap_pairwise_crossing(
        pts, 7, 11, metric="shot", n_boot=50, seed=0
    ) == bootstrap_pairwise_crossing(shuffled, 7, 11, metric="shot", n_boot=50, seed=0)


def test_bootstrap_ci_brackets_point_estimate() -> None:
    est, (lo, hi), n_no = bootstrap_pairwise_crossing(
        _synthetic_points(), 7, 11, metric="shot", n_boot=200, seed=0
    )
    assert n_no == 0
    assert lo <= est.value <= hi


# --- saturation cap ---------------------------------------------------------------

# Both curves sit in the saturated tail (P_L_shot 0.44 .. 0.49) and cross between
# p = 0.01 and p = 0.02 (g < 0, then g >= 0).
_SATURATED_E1 = [4500, 4700, 4800]
_SATURATED_E2 = [4400, 4800, 4900]
CAP = 0.4


def test_crossing_in_saturated_tail_is_classified_saturated() -> None:
    pts = _pair_points(5, 7, _SATURATED_E1, _SATURATED_E2)
    plain = pairwise_crossing(pts, 5, 7, metric="shot")
    capped = pairwise_crossing(pts, 5, 7, metric="shot", saturation_cap=CAP)
    assert math.isfinite(plain.value) and not plain.saturated
    assert capped.saturated
    assert math.isnan(capped.value)
    assert capped.bracket == plain.bracket == (0.01, 0.02)
    assert capped.bracket_max_p_l_shot == pytest.approx(0.48)


def test_saturated_bootstrap_replicates_count_as_no_crossing() -> None:
    pts = _pair_points(5, 7, _SATURATED_E1, _SATURATED_E2)
    _, _, n_no_plain = bootstrap_pairwise_crossing(pts, 5, 7, metric="shot", n_boot=200)
    _, ci, n_no = bootstrap_pairwise_crossing(
        pts, 5, 7, metric="shot", n_boot=200, saturation_cap=CAP
    )
    assert n_no_plain < 200
    assert n_no == 200
    assert all(math.isnan(x) for x in ci)


def test_unsaturated_crossing_is_unchanged_by_the_cap() -> None:
    pts = _synthetic_points()
    for d1, d2 in ADJACENT_PAIRS:
        plain = bootstrap_pairwise_crossing(pts, d1, d2, metric="shot", n_boot=200)
        capped = bootstrap_pairwise_crossing(
            pts, d1, d2, metric="shot", n_boot=200, saturation_cap=CAP
        )
        assert capped == plain
        assert not capped[0].saturated
        assert capped[0].bracket_max_p_l_shot < CAP
