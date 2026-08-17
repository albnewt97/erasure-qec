"""Threshold-fit recovery and window rule (PLAN.md §10)."""

from pathlib import Path

import pytest

from erasure_qec.analysis.statistics import SweepPoint, shot_p_l_from_per_round
from erasure_qec.analysis.threshold_fit import estimate_crossing, fit_threshold

FIXTURES = Path(__file__).resolve().parent / "fixtures"

# Ground-truth quadratic ansatz p_L = A + B x + C x^2, x = (p - p_th) d^(1/nu).
TRUE_P_TH = 0.037
TRUE_NU = 1.5
TRUE_A, TRUE_B, TRUE_C = 0.10, 0.9, 0.4


def _quadratic_points(shots: int = 400_000) -> list[SweepPoint]:
    """Exact-ansatz points on a near-threshold grid (quadratic stays positive).

    Deterministic: ``errors`` is the rounded expected count from inverting the
    per-round rate over T = d rounds, so the fit recovers the injected p_th
    exactly up to rounding.
    """
    points: list[SweepPoint] = []
    # A tight band around p_th keeps the quadratic ansatz valid for all d.
    p_values = [TRUE_P_TH * f for f in (0.75, 0.85, 0.95, 1.05, 1.15, 1.25)]
    for d in (3, 5, 7, 9):
        for p in p_values:
            x = (p - TRUE_P_TH) * d ** (1.0 / TRUE_NU)
            p_l = TRUE_A + TRUE_B * x + TRUE_C * x * x
            errors = round(shots * shot_p_l_from_per_round(p_l, d))
            points.append(SweepPoint("herald_mwpm", d, d, p, 0.5, shots, errors))
    return points


def test_recovers_known_p_th_from_ansatz() -> None:
    points = _quadratic_points()
    result = fit_threshold(points, window_factor=1.4, n_boot=200, seed=0)
    assert result.converged, result.message
    # p_th recovered within a few bootstrap sigma of the injected value.
    assert abs(result.p_th - TRUE_P_TH) < 3 * result.p_th_err + 1e-4
    assert result.nu == pytest.approx(TRUE_NU, abs=0.3)
    assert result.p_th_err > 0


def test_crossing_estimate_near_true_threshold() -> None:
    assert estimate_crossing(_quadratic_points()) == pytest.approx(TRUE_P_TH, rel=0.35)


def test_window_rule_limits_points_and_is_recorded() -> None:
    points = _quadratic_points()
    narrow = fit_threshold(points, window_factor=1.1, n_boot=50, seed=0)
    wide = fit_threshold(points, window_factor=1.3, n_boot=50, seed=0)
    # A wider window admits at least as many points; the factor is recorded.
    assert wide.n_points >= narrow.n_points
    assert narrow.window_factor == 1.1 and wide.window_factor == 1.3


def test_fit_degrades_gracefully_on_insufficient_data() -> None:
    one_curve = [SweepPoint("herald_mwpm", 3, 3, p, 0.5, 10000, 100) for p in (0.01, 0.02)]
    result = fit_threshold(one_curve)
    assert not result.converged
    assert "insufficient" in result.message


@pytest.mark.parametrize("name", ["baseline_pauli", "erasure_r50", "erasure_r98"])
def test_fits_synthetic_fixtures(name: str) -> None:
    """The committed exponential-collapse fixtures fit near their injected p_th
    (0.010 / 0.024 / 0.045), confirming the fit works on figure data too."""
    from erasure_qec.analysis.statistics import load_sweep

    injected = {"baseline_pauli": 0.010, "erasure_r50": 0.024, "erasure_r98": 0.045}[name]
    points = [p for p in load_sweep(FIXTURES / f"{name}.csv") if p.decoder == "herald_mwpm"]
    result = fit_threshold(points, window_factor=1.6, n_boot=60, seed=0)
    assert result.converged, result.message
    assert result.p_th == pytest.approx(injected, rel=0.15)


# --- Real Monte-Carlo baseline sweep: crossing/threshold regression. ---


def _real_baseline_path() -> Path:
    """The committed snapshot of the real Monte-Carlo baseline sweep.

    Pinned to the committed fixture (not the live, gitignored
    data/baseline_pauli.csv) so these regressions are deterministic: the live
    sweep is re-collected over time (e.g. densified through the crossing),
    which shifts the fit; a regression test must assert against a fixed grid.
    """
    return FIXTURES / "real_baseline_pauli.csv"


def _real_baseline(decoder: str) -> list[SweepPoint]:
    from erasure_qec.analysis.statistics import load_sweep

    return [p for p in load_sweep(_real_baseline_path()) if p.decoder == decoder]


def _dense_quadratic_sweep(
    decoder: str, *, penalty: float = 1.0, shots: int = 300_000
) -> list[SweepPoint]:
    """A synthetic sweep with ~8 p-points *through* the crossing (a dense fit-
    resolution band), following the exact quadratic ansatz with nu = 1.5.

    ``penalty >= 1`` inflates a decoder's rate (to mimic a worse decoder) while
    keeping the same nu, so both decoders must still recover nu = 1.5.
    """
    p_th, nu, a, b, c = 0.015, 1.5, 0.05, 1.0, 2.0
    p_grid = [0.008, 0.011, 0.013, 0.015, 0.017, 0.019, 0.021, 0.024]
    points: list[SweepPoint] = []
    for d in (3, 5, 7, 9, 11):
        for p in p_grid:
            x = (p - p_th) * d ** (1.0 / nu)
            p_l = min(max(penalty * (a + b * x + c * x * x), 1e-6), 0.45)
            errors = round(shots * shot_p_l_from_per_round(p_l, d))
            points.append(SweepPoint(decoder, d, d, p, 0.0, shots, errors))
    return points


def test_dense_grid_recovers_nu_for_both_decoders() -> None:
    """The real acceptance bar the sparse grid could not meet: with ~8 points
    through the crossing, a DEFAULT fit_threshold (auto crossing, default
    window) recovers nu in [1.2, 1.8] for BOTH herald and blind decoders. This
    is why baseline_pauli.yaml adds the dense extra_p band."""
    for decoder, penalty in (("herald_mwpm", 1.0), ("blind_mwpm", 1.15)):
        result = fit_threshold(_dense_quadratic_sweep(decoder, penalty=penalty))
        assert result.converged, (decoder, result.message)
        assert 1.2 <= result.nu <= 1.8, (decoder, result.nu)
        assert result.p_th == pytest.approx(0.015, abs=0.002), (decoder, result.p_th)
        # The dense band is not filtered out by the window rule.
        assert result.n_points >= 25, (decoder, result.n_points)


def test_estimate_crossing_on_real_baseline_is_near_threshold() -> None:
    """estimate_crossing must find the true ~1.5% crossing, NOT the low edge of
    the p-grid where the deep-sub-threshold curves all vanish together."""
    crossing = estimate_crossing(_real_baseline("herald_mwpm"))
    assert 0.012 <= crossing <= 0.020, crossing


def test_default_fit_on_real_baseline_recovers_threshold() -> None:
    """A default fit_threshold call (crossing auto-estimated, no hand-supplied
    p_center) recovers an effective crossing p_th in [0.012, 0.020] for both
    decoders with a physical herald nu.

    On the densified baseline the all-d fit's bootstrap errors are tight
    (~0.0008), and the all-d effective crossing is finite-size-biased *by a
    different amount* for each decoder, so the two all-d p_th values no longer
    agree at 1 sigma. The all-d effective nu is itself an artefact of mixing
    small and large distances: it is inflated well above the physical value
    (~2.6 for herald here, after the coin-flip fit-window cut raised it from
    the earlier ~2.2), which is exactly why we report the asymptotic d>=7 fit
    instead (see :func:`test_asymptotic_d_min_fit_on_real_baseline`)."""
    herald = fit_threshold(_real_baseline("herald_mwpm"))
    blind = fit_threshold(_real_baseline("blind_mwpm"))

    assert herald.converged and blind.converged
    # Both all-d effective crossings sit in-band around the ~1.4-1.8% region.
    assert 0.012 <= herald.p_th <= 0.020, herald.p_th
    assert 0.012 <= blind.p_th <= 0.020, blind.p_th

    # Asymptotic (d>=7) thresholds agree within combined bootstrap error bars,
    # and the all-d effective nu is inflated above the asymptotic nu by the
    # finite-size corrections that d_min=7 removes (rather than asserting a
    # magic band on the physically meaningless all-d nu).
    herald7 = fit_threshold(_real_baseline("herald_mwpm"), d_min=7)
    blind7 = fit_threshold(_real_baseline("blind_mwpm"), d_min=7)
    assert herald.nu > herald7.nu, (herald.nu, herald7.nu)
    assert abs(herald7.p_th - blind7.p_th) <= herald7.p_th_err + blind7.p_th_err


def test_asymptotic_d_min_fit_on_real_baseline() -> None:
    """Dropping the small distances (largest finite-size corrections) via
    d_min=7 gives the asymptotic threshold ~1.4% with a physical nu, lower than
    the all-d effective crossing (~1.7%, finite-size drift)."""
    result = fit_threshold(_real_baseline("herald_mwpm"), d_min=7)
    assert result.converged, result.message
    assert result.d_min == 7
    assert result.distances == (7, 9, 11)
    assert 0.012 <= result.p_th <= 0.016, result.p_th
    assert 1.2 <= result.nu <= 1.8, result.nu


# --- Real Monte-Carlo r_e=0.5 sweep: estimate_crossing ordering-inversion guard. ---


def _real_erasure_r50(decoder: str) -> list[SweepPoint]:
    """Committed snapshot of the real r_e=0.5 sweep (pre-densification grid)."""
    from erasure_qec.analysis.statistics import load_sweep

    return [p for p in load_sweep(FIXTURES / "real_erasure_r50.csv") if p.decoder == decoder]


def test_estimate_crossing_r50_not_pulled_by_saturated_tail() -> None:
    """Regression for the ragged r_e=0.5 tail: the old estimator returned
    ~0.0316 for herald_mwpm (already above threshold, where the curves
    re-converge toward 1/2), instead of the true ~2.2% crossing where the
    d-ordering inverts. The ordering-inversion guard must reject the
    above-threshold minimum and land the estimate at the transition."""
    herald = estimate_crossing(_real_erasure_r50("herald_mwpm"))
    blind = estimate_crossing(_real_erasure_r50("blind_mwpm"))
    assert 0.018 <= herald <= 0.026, herald
    assert 0.015 <= blind <= 0.024, blind
