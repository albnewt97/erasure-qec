"""Unit tests for per-round conversion and CIs (PLAN.md §10)."""

import math

import pytest

from erasure_qec.analysis.statistics import (
    Estimate,
    SweepPoint,
    bootstrap_per_round,
    lambda_factor,
    per_round_estimate,
    per_round_p_l,
    shot_p_l_from_per_round,
    wilson_interval,
)


def test_scrambled_fixed_point() -> None:
    """P_L_shot = 1/2 is the fully-scrambled fixed point: p_L = 1/2 for any T."""
    for rounds in (1, 3, 7, 25):
        assert per_round_p_l(0.5, rounds) == pytest.approx(0.5)


def test_single_round_is_identity() -> None:
    """T = 1: the per-round rate equals the shot-level rate."""
    for p in (0.0, 0.01, 0.123, 0.4, 0.5):
        assert per_round_p_l(p, 1) == pytest.approx(p)


def test_conversion_round_trips_with_inverse() -> None:
    # The round trip is well-conditioned away from the p_L -> 1/2 scrambling
    # limit. At p_L very close to 1/2 with large T the shot-level rate
    # saturates ((1-2 p_L)^T underflows below machine epsilon vs 1), so it is
    # physically non-invertible; that regime is covered by the fixed-point test.
    for p_l in (0.0, 1e-4, 0.02, 0.2, 0.4):
        for rounds in (1, 3, 5, 11):
            back = per_round_p_l(shot_p_l_from_per_round(p_l, rounds), rounds)
            assert back == pytest.approx(p_l, abs=1e-9)


def test_conversion_matches_closed_form() -> None:
    # p_L = 1/2 (1 - (1 - 2 P)^(1/T)) evaluated by hand.
    assert per_round_p_l(0.2, 5) == pytest.approx(0.5 * (1 - 0.6 ** (1 / 5)))


def test_conversion_clamps_above_half() -> None:
    """P_L_shot > 1/2 (noise-limited) clamps rather than going complex."""
    val = per_round_p_l(0.7, 5)
    assert math.isfinite(val) and val == pytest.approx(0.5)


def test_per_round_is_monotonic_in_shot_rate() -> None:
    xs = [per_round_p_l(p, 7) for p in (0.01, 0.05, 0.1, 0.2, 0.3)]
    assert all(a < b for a, b in zip(xs, xs[1:], strict=False))


def test_wilson_interval_brackets_estimate() -> None:
    lo, hi = wilson_interval(50, 1000)
    assert lo < 0.05 < hi
    assert 0.0 <= lo < hi <= 1.0
    # Wider interval for fewer shots at the same proportion.
    lo2, hi2 = wilson_interval(5, 100)
    assert (hi2 - lo2) > (hi - lo)


def test_per_round_estimate_orders_low_value_high() -> None:
    pt = SweepPoint("herald_mwpm", d=5, rounds=5, p=0.01, r_e=0.5, shots=10000, errors=284)
    est = per_round_estimate(pt)
    assert est.low < est.value < est.high
    assert est.value == pytest.approx(per_round_p_l(0.0284, 5))


def test_bootstrap_is_deterministic_and_brackets() -> None:
    pt = SweepPoint("herald_mwpm", d=3, rounds=3, p=0.01, r_e=0.5, shots=10000, errors=450)
    a = bootstrap_per_round(pt, seed=0)
    b = bootstrap_per_round(pt, seed=0)
    assert a == b  # same seed -> identical CI
    assert a.low < a.value < a.high
    assert isinstance(a, Estimate)


def test_lambda_factor_ratio_and_validation() -> None:
    p3 = SweepPoint("herald_mwpm", d=3, rounds=3, p=0.01, r_e=0.5, shots=10000, errors=450)
    p5 = SweepPoint("herald_mwpm", d=5, rounds=5, p=0.01, r_e=0.5, shots=10000, errors=284)
    lam = lambda_factor(p3, p5, seed=1)
    expected = per_round_p_l(0.045, 3) / per_round_p_l(0.0284, 5)
    assert lam.value == pytest.approx(expected)
    assert lam.value > 1.0  # bigger code suppresses errors below threshold
    with pytest.raises(ValueError):
        lambda_factor(p3, p3)  # not a d, d+2 pair
