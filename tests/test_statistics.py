"""Unit tests for per-round conversion and CIs (PLAN.md §10)."""

import math

import pytest

from erasure_qec.analysis.statistics import (
    Estimate,
    LambdaEstimate,
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


def test_lambda_factor_reports_zero_dropped_on_healthy_statistics() -> None:
    """With hundreds of errors in both distances, no replicate is discarded."""
    p3 = SweepPoint("herald_mwpm", d=3, rounds=3, p=0.01, r_e=0.5, shots=10000, errors=450)
    p5 = SweepPoint("herald_mwpm", d=5, rounds=5, p=0.01, r_e=0.5, shots=10000, errors=284)
    lam = lambda_factor(p3, p5, seed=1)
    assert isinstance(lam, LambdaEstimate)
    assert isinstance(lam, Estimate)  # existing callers keep working
    assert lam.n_boot == 2000
    assert lam.n_boot_dropped == 0
    assert lam.dropped_fraction == 0.0
    assert lam.ci_conditioned is False


def test_lambda_factor_surfaces_dropped_replicates_on_low_error_d2() -> None:
    """A d+2 point with a handful of errors silently conditions the CI.

    ``d + 2`` logs 3 errors in 100k shots, so a meaningful share of bootstrap
    replicates draw zero errors there, hit a zero denominator, and are excluded.
    Those are exactly the replicates whose Lambda is ``+inf``, so the discards
    bias ``high`` low -- the count must be visible on the returned value rather
    than swallowed by the loop.
    """
    p3 = SweepPoint("herald_mwpm", d=3, rounds=3, p=0.001, r_e=0.98, shots=100000, errors=200)
    p5 = SweepPoint("herald_mwpm", d=5, rounds=5, p=0.001, r_e=0.98, shots=100000, errors=3)
    lam = lambda_factor(p3, p5, seed=7)
    assert lam.n_boot == 2000
    assert lam.n_boot_dropped > 0, "zero-denominator replicates must be counted"
    # P(Binomial(1e5, 3e-5) == 0) = e^-3 ~ 0.0498; 2000 draws puts this near 100.
    assert lam.n_boot_dropped == pytest.approx(2000 * math.exp(-3.0), rel=0.25)
    assert 0.0 < lam.dropped_fraction < 1.0
    # 5% dropped vs the 2.5% upper tail the 95% CI claims to measure: the
    # unconditional upper quantile is +inf, so the reported `high` is not it.
    assert lam.ci_conditioned is True
    assert math.isfinite(lam.high)  # ...yet a finite number is still returned


def test_lambda_factor_dropped_count_is_seed_deterministic() -> None:
    p3 = SweepPoint("herald_mwpm", d=3, rounds=3, p=0.001, r_e=0.98, shots=100000, errors=200)
    p5 = SweepPoint("herald_mwpm", d=5, rounds=5, p=0.001, r_e=0.98, shots=100000, errors=3)
    assert lambda_factor(p3, p5, seed=7) == lambda_factor(p3, p5, seed=7)


def test_ci_conditioned_threshold_tracks_the_requested_ci_level() -> None:
    """The flag compares discards against *this* interval's tail mass."""
    lo_drop = LambdaEstimate(value=2.0, low=1.0, high=3.0, n_boot=1000, n_boot_dropped=30)
    assert lo_drop.dropped_fraction == pytest.approx(0.03)
    assert lo_drop.ci_conditioned is True  # 3% > 2.5% tail of a 95% CI
    # A 90% CI has a 5% upper tail, so the same 3% no longer dominates it.
    assert (
        LambdaEstimate(
            value=2.0, low=1.0, high=3.0, n_boot=1000, n_boot_dropped=30, ci=0.90
        ).ci_conditioned
        is False
    )
    # No replicates drawn at all -> nothing to condition on, and no false "0%".
    none_drawn = LambdaEstimate(value=2.0, low=1.0, high=3.0)
    assert math.isnan(none_drawn.dropped_fraction)
    assert none_drawn.ci_conditioned is False
