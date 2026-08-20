"""Threshold-fit recovery and window rule (PLAN.md §10)."""

import math
from pathlib import Path

import numpy as np
import pytest

from erasure_qec.analysis.statistics import SweepPoint, load_sweep, shot_p_l_from_per_round
from erasure_qec.analysis.threshold_fit import (
    bootstrap_threshold_difference,
    directional_sensitivity_ci,
    estimate_crossing,
    failed_vs_success_pth,
    failure_guard_breakdown,
    fit_threshold,
)

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


@pytest.mark.parametrize(
    "name", ["synthetic_baseline_pauli", "synthetic_erasure_r50", "synthetic_erasure_r98"]
)
def test_fits_synthetic_fixtures(name: str) -> None:
    """The committed exponential-collapse fixtures fit near their injected p_th
    (0.010 / 0.024 / 0.045). NB this confirms only that the fitter inverts its
    own generative model on figure data -- the fixtures are synthetic, not
    measurements (see analysis/synthetic.py provenance warning)."""
    from erasure_qec.analysis.statistics import load_sweep

    injected = {
        "synthetic_baseline_pauli": 0.010,
        "synthetic_erasure_r50": 0.024,
        "synthetic_erasure_r98": 0.045,
    }[name]
    points = [p for p in load_sweep(FIXTURES / f"{name}.csv") if p.decoder == "herald_mwpm"]
    result = fit_threshold(points, window_factor=1.6, n_boot=60, seed=0)
    assert result.converged, result.message
    assert result.p_th == pytest.approx(injected, rel=0.15)


def test_fitter_rejects_adversarial_saturated_fixture() -> None:
    """The adversarial fixture (inverted d-ordering + saturated ~1/2 tail, no
    real crossing) must be REJECTED, not fit: fit_threshold returns
    converged=False for both decoders and both the d>=7 and all-d windows,
    rather than returning a spurious p_th."""
    from erasure_qec.analysis.statistics import load_sweep

    points = load_sweep(FIXTURES / "synthetic_adversarial.csv")
    assert points  # fixture loads
    for decoder in ("herald_mwpm", "blind_mwpm"):
        sub = [p for p in points if p.decoder == decoder]
        for d_min in (7, None):
            result = fit_threshold(sub, d_min=d_min)
            assert not result.converged, (decoder, d_min, result.p_th, result.message)


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

    The all-d effective crossing is finite-size-biased *by a different amount*
    for each decoder, so the two all-d p_th values no longer agree at 1 sigma.
    The all-d effective nu is itself an artefact of mixing small and large
    distances: it is inflated above the physical value (~2.2 for herald here),
    which is exactly why we report the asymptotic d>=7 fit instead (see
    :func:`test_asymptotic_d_min_fit_on_real_baseline`)."""
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
    """Committed snapshot of the real r_e=0.5 sweep (fixed noise model).

    Re-collected under the corrected constant-budget model (herald rate
    (2/3)*p*r_e); the earlier snapshot used the un-fixed p*r_e/2 rate. The
    higher erasure rate shifts the herald crossing up (~2.2% -> ~2.6%).
    """
    from erasure_qec.analysis.statistics import load_sweep

    return [p for p in load_sweep(FIXTURES / "real_erasure_r50.csv") if p.decoder == decoder]


# --- Bootstrap CI calibration: coverage of a known p_th at ~nominal rate. ---

_COV_P_TH, _COV_NU = 0.015, 1.5
_COV_A, _COV_B, _COV_C = 0.05, 1.0, 2.0
_COV_GRID = [0.008, 0.011, 0.013, 0.015, 0.017, 0.019, 0.021, 0.024]
_COV_DIST = (3, 5, 7, 9, 11)
_COV_SHOTS = 30_000


def _cov_truth(p: float, d: int) -> float:
    x = (p - _COV_P_TH) * d ** (1.0 / _COV_NU)
    return min(max(_COV_A + _COV_B * x + _COV_C * x * x, 1e-6), 0.45)


def _noisy_realization(rng: np.random.Generator) -> list[SweepPoint]:
    """One sweep drawn with real Binomial noise around the known ansatz."""
    points: list[SweepPoint] = []
    for d in _COV_DIST:
        for p in _COV_GRID:
            p_shot = shot_p_l_from_per_round(_cov_truth(p, d), d)
            errors = int(rng.binomial(_COV_SHOTS, p_shot))
            points.append(SweepPoint("herald_mwpm", d, d, p, 0.0, _COV_SHOTS, errors))
    return points


@pytest.mark.slow
def test_bootstrap_ci_covers_known_p_th_at_nominal_rate() -> None:
    """The 95% bootstrap percentile CI must cover the injected p_th at roughly
    its nominal rate over many independent noisy realizations. The old bootstrap
    (unweighted refit from popt on a frozen window, failures dropped) reports a
    CI too narrow to hit this; the full-pipeline bootstrap does. Seeded, so the
    coverage number is deterministic.

    The CI must also stay *informative* (not trivially wide): a fit could reach
    100% coverage by returning an absurd interval, so the median CI width is
    bounded too.
    """
    rng = np.random.default_rng(12345)
    n_real = 40
    covered = 0
    n_converged = 0
    rel_widths: list[float] = []
    for i in range(n_real):
        fit = fit_threshold(_noisy_realization(rng), n_boot=100, seed=i)
        if not fit.converged:
            continue
        n_converged += 1
        lo, hi = fit.p_th_ci
        if lo <= _COV_P_TH <= hi:
            covered += 1
        rel_widths.append((hi - lo) / _COV_P_TH)

    assert n_converged >= 0.9 * n_real, n_converged
    coverage = covered / n_converged
    # Nominal 0.95; allow Monte-Carlo slack + the mild under-coverage typical of
    # bootstrap percentile intervals on a nonlinear fit. Well below 0.8 would
    # mean the CI is too narrow (the bug this fix removes).
    assert 0.80 <= coverage <= 1.0, (coverage, covered, n_converged)
    # Informative: the interval is a small fraction of p_th, not absurdly wide.
    assert float(np.median(rel_widths)) < 0.20, float(np.median(rel_widths))


def test_estimate_crossing_r50_not_pulled_by_saturated_tail() -> None:
    """Regression for the ragged r_e=0.5 tail: a naive estimator returns the
    saturated high-p tail (where curves re-converge toward 1/2 and the ordering
    inverts, e.g. p>=0.046 here has d9,d11 at the 1/2 coin-flip limit), instead
    of the true crossing. The ordering-inversion guard must reject the
    above-threshold region and land the estimate at the transition (~2.6%
    herald, ~2.0% blind under the fixed model)."""
    herald = estimate_crossing(_real_erasure_r50("herald_mwpm"))
    blind = estimate_crossing(_real_erasure_r50("blind_mwpm"))
    assert 0.022 <= herald <= 0.030, herald
    assert 0.016 <= blind <= 0.024, blind


# --- Paired/unpaired bootstrap of the threshold difference (Phase 4). ---

_PD_NU = 1.5
_PD_DIST = (3, 5, 7, 9, 11)
_PD_P = (0.008, 0.010, 0.012, 0.014, 0.016, 0.018, 0.020, 0.024, 0.028)


def _pd_rate(p: float, p_th: float, d: int) -> float:
    x = (p - p_th) * d ** (1.0 / _PD_NU)
    return min(max(0.08 * math.exp(30.0 * x), 1e-6), 0.49)


def _paired_synthetic(
    p_th_h: float, p_th_b: float, shots: int
) -> tuple[dict[tuple[float, int], tuple[int, int, int, int]], list[SweepPoint], list[SweepPoint]]:
    """Nested-failure joint (herald-wrong is a subset of blind-wrong, so the two
    are strongly positively correlated on shared shots) with a known
    ``Delta = p_th_b - p_th_h``. Returns (joint, herald_points, blind_points)."""
    joint: dict[tuple[float, int], tuple[int, int, int, int]] = {}
    herald: list[SweepPoint] = []
    blind: list[SweepPoint] = []
    for d in _PD_DIST:
        for p in _PD_P:
            ph = shot_p_l_from_per_round(_pd_rate(p, p_th_h, d), d)
            pb = max(shot_p_l_from_per_round(_pd_rate(p, p_th_b, d), d), ph)
            n11 = round(shots * ph)  # both wrong
            n01 = round(shots * (pb - ph))  # blind-only wrong
            n00 = shots - n11 - n01
            joint[(p, d)] = (n00, 0, n01, n11)  # (both ok, herald-only, blind-only, both wrong)
            herald.append(SweepPoint("herald_mwpm", d, d, p, 0.0, shots, n11))
            blind.append(SweepPoint("blind_mwpm", d, d, p, 0.0, shots, n01 + n11))
    return joint, herald, blind


def _decorrelate(
    joint: dict[tuple[float, int], tuple[int, int, int, int]],
) -> dict[tuple[float, int], tuple[int, int, int, int]]:
    """Same per-(p,d) MARGINALS but herald/blind made independent (product joint)
    -- a 'shuffled pairing' that destroys the correlation while preserving each
    decoder's error rate."""
    out: dict[tuple[float, int], tuple[int, int, int, int]] = {}
    for (p, d), (n00, n10, n01, n11) in joint.items():
        n = n00 + n10 + n01 + n11
        ph = (n10 + n11) / n
        pb = (n01 + n11) / n
        c = [
            round(n * (1 - ph) * (1 - pb)),
            round(n * ph * (1 - pb)),
            round(n * (1 - ph) * pb),
            round(n * ph * pb),
        ]
        c[0] += n - sum(c)  # absorb rounding so counts sum to n
        out[(p, d)] = (c[0], c[1], c[2], c[3])
    return out


def _marginals(
    joint: dict[tuple[float, int], tuple[int, int, int, int]],
) -> tuple[list[SweepPoint], list[SweepPoint]]:
    herald: list[SweepPoint] = []
    blind: list[SweepPoint] = []
    for (p, d), (n00, n10, n01, n11) in joint.items():
        n = n00 + n10 + n01 + n11
        herald.append(SweepPoint("herald_mwpm", d, d, p, 0.0, n, n10 + n11))
        blind.append(SweepPoint("blind_mwpm", d, d, p, 0.0, n, n01 + n11))
    return herald, blind


def _ci_width(fit: object) -> float:
    lo, hi = fit.p_th_ci  # type: ignore[attr-defined]
    return hi - lo


@pytest.mark.slow
def test_paired_bootstrap_ci_covers_known_delta() -> None:
    """Over many seeded noisy realizations of a known Delta, the paired 95% CI
    covers the truth at ~nominal rate, and stays informative (median width
    bounded, so a trivially-wide CI cannot pass on coverage alone)."""
    p_th_h, p_th_b = 0.020, 0.017
    true_delta = p_th_b - p_th_h  # -0.003
    expected, _, _ = _paired_synthetic(p_th_h, p_th_b, shots=40_000)
    rng = np.random.default_rng(7)
    covered = 0
    n_ok = 0
    widths: list[float] = []
    for i in range(30):
        realized = {
            k: tuple(int(x) for x in rng.multinomial(sum(c), np.array(c, float) / sum(c)))
            for k, c in expected.items()
        }
        herald, blind = _marginals(realized)  # type: ignore[arg-type]
        r = bootstrap_threshold_difference(
            herald, blind, joint_counts=realized, n_boot=150, seed=i  # type: ignore[arg-type]
        )
        if not r.converged:
            continue
        n_ok += 1
        lo, hi = r.delta_ci
        if lo <= true_delta <= hi:
            covered += 1
        widths.append(hi - lo)
    assert n_ok >= 27, n_ok
    coverage = covered / n_ok
    assert 0.80 <= coverage <= 1.0, (coverage, covered, n_ok)
    assert float(np.median(widths)) < 0.002, float(np.median(widths))


def test_paired_ci_narrower_than_naive_marginal_difference() -> None:
    """The property that justifies the branch: the difference bootstrap CI is
    narrower than the interval you'd get by naively differencing the two
    marginal CIs (whose width is the SUM of the marginal widths)."""
    joint, herald, blind = _paired_synthetic(0.020, 0.017, shots=40_000)
    diff = bootstrap_threshold_difference(
        herald, blind, joint_counts=joint, n_boot=400, seed=2
    )
    herald_fit = fit_threshold(herald, n_boot=400, seed=2)
    blind_fit = fit_threshold(blind, n_boot=400, seed=2)
    naive_width = _ci_width(herald_fit) + _ci_width(blind_fit)
    paired_width = diff.delta_ci[1] - diff.delta_ci[0]
    assert diff.converged
    assert paired_width < naive_width, (paired_width, naive_width)


def test_shuffled_pairing_widens_ci() -> None:
    """Pairing must actually be applied: decorrelating the shared-shot joint
    (same marginals, correlation removed) widens the Delta CI. If the code
    silently ignored the correlation, both CIs would match."""
    joint, herald, blind = _paired_synthetic(0.020, 0.017, shots=40_000)
    shuffled = _decorrelate(joint)
    paired = bootstrap_threshold_difference(
        herald, blind, joint_counts=joint, n_boot=400, seed=3
    )
    unpaired = bootstrap_threshold_difference(
        herald, blind, joint_counts=shuffled, n_boot=400, seed=3
    )
    w_paired = paired.delta_ci[1] - paired.delta_ci[0]
    w_shuffled = unpaired.delta_ci[1] - unpaired.delta_ci[0]
    assert paired.correlation > unpaired.correlation  # correlation was used
    assert w_shuffled > w_paired, (w_shuffled, w_paired)


@pytest.mark.slow
def test_real_baseline_delta_is_control_consistent_with_zero() -> None:
    """r_e = 0 control: herald and blind are the SAME decoder (no heralds), so
    Delta must be consistent with zero. If this fails the pairing/pipeline is
    broken and nothing else on the branch is trustworthy."""
    pts = load_sweep(FIXTURES / "real_baseline_pauli.csv")
    herald = [p for p in pts if p.decoder == "herald_mwpm"]
    blind = [p for p in pts if p.decoder == "blind_mwpm"]
    r = bootstrap_threshold_difference(herald, blind, d_min=7, n_boot=400, seed=0)
    assert r.converged, r.message
    assert not r.excludes_zero, (r.delta, r.delta_ci)
    assert r.delta_ci[0] <= 0.0 <= r.delta_ci[1], r.delta_ci


@pytest.mark.slow
def test_real_r50_delta_excludes_zero_significant_separation() -> None:
    """r_e = 0.5: the herald-vs-blind threshold separation IS significant under
    the correct statistic (the difference bootstrap CI excludes zero), even
    though the marginal CIs nearly touch. Blind threshold is below herald, so
    Delta = blind - herald < 0."""
    pts = load_sweep(FIXTURES / "real_erasure_r50.csv")
    herald = [p for p in pts if p.decoder == "herald_mwpm"]
    blind = [p for p in pts if p.decoder == "blind_mwpm"]
    r = bootstrap_threshold_difference(herald, blind, d_min=7, n_boot=400, seed=0)
    assert r.converged, r.message
    assert r.excludes_zero, (r.delta, r.delta_ci)
    assert r.delta < 0.0 and r.delta_ci[1] < 0.0, (r.delta, r.delta_ci)


@pytest.mark.slow
def test_r50_discards_not_mar_but_significance_survives_imputation() -> None:
    """Phase-4 Task 1: the r_e=0.5 Delta CI is conditioned on both decoders
    converging, and ~10% of replicates are discarded. Those discards are NOT
    missing-at-random -- they cluster at high blind p_th (almost all are
    bound-pinning of the bistable blind crossing) -- so the caveat is real.
    But the significance SURVIVES imputing the discards from the least-negative
    decile of Delta. Pins both the caveat and its resolution."""
    pts = load_sweep(FIXTURES / "real_erasure_r50.csv")
    herald = [p for p in pts if p.decoder == "herald_mwpm"]
    blind = [p for p in pts if p.decoder == "blind_mwpm"]
    r = bootstrap_threshold_difference(herald, blind, d_min=7, n_boot=1000, seed=0)
    assert r.n_paired_failed > 0 and r.failures

    # Almost all discards are bound-pinning (informative), not sparse windows,
    # so there is no mechanical cause to fix (Task 1c/1d).
    breakdown = failure_guard_breakdown(r)
    bound = sum(v for k, v in breakdown.items() if k.endswith("bound_pin"))
    assert bound >= 0.9 * sum(breakdown.values()), breakdown

    # NOT missing-at-random: failed-blind p_th sits above kept-blind p_th (1b).
    failed_blind, success_blind = failed_vs_success_pth(r, "blind")
    assert float(np.median(failed_blind)) > float(np.median(success_blind))

    # ... yet the significance survives the plausible-pessimistic imputation (1e).
    _, excludes_zero = directional_sensitivity_ci(r, seed=0, top_fraction=0.10)
    assert excludes_zero


@pytest.mark.slow
def test_baseline_control_discards_are_missing_at_random() -> None:
    """r_e=0 control: discards should be missing-at-random (blind p_th similar
    for failed and kept replicates), and Delta stays consistent with zero even
    under the directional imputation. If this fails the discard mechanism is
    decoder-asymmetric where it must not be."""
    pts = load_sweep(FIXTURES / "real_baseline_pauli.csv")
    herald = [p for p in pts if p.decoder == "herald_mwpm"]
    blind = [p for p in pts if p.decoder == "blind_mwpm"]
    r = bootstrap_threshold_difference(herald, blind, d_min=7, n_boot=1000, seed=0)
    failed_blind, success_blind = failed_vs_success_pth(r, "blind")
    assert abs(float(np.median(failed_blind)) - float(np.median(success_blind))) < 0.002
    _, excludes_zero = directional_sensitivity_ci(r, seed=0, top_fraction=0.10)
    assert not excludes_zero  # control: still consistent with zero
