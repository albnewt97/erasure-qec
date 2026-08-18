"""Finite-size scaling threshold fit (PLAN.md §10).

Near the threshold the per-round logical rate collapses onto a universal curve
of the single scaling variable ``x = (p - p_th) * d^(1/nu)``. We fit the
quadratic ansatz

    p_L(p, d) = A + B x + C x^2

for the five parameters ``(p_th, nu, A, B, C)`` by weighted least squares
(``scipy.optimize.curve_fit``), with per-point weights from the Wilson
half-width of the per-round rate.

Data-window selection rule
--------------------------
The quadratic ansatz is only valid *near* the crossing, so the fit uses only
points whose physical error rate lies within a multiplicative ``window_factor``
of an initial crossing estimate ``p_center``:

    p_center / window_factor  <=  p  <=  p_center * window_factor

``p_center`` is estimated data-drivenly by :func:`estimate_crossing` as the p
that minimizes the spread of per-round ``p_L`` across code distances (curves
for different d intersect at threshold). ``window_factor`` defaults to 1.5 and
is configurable; widen it when the sweep is coarse, narrow it when the ansatz
starts to bend. The exact window and the points used are recorded on the
returned :class:`FitResult`.

Confidence on ``p_th`` comes from a seeded parametric bootstrap that re-runs the
*entire* pipeline per replicate: resample ``errors ~ Binomial(shots, P_L_shot)``
over **all** points for the decoder (not just the previously-selected window),
then re-run :func:`estimate_crossing`, the window selection, and the *weighted*
``curve_fit`` with the same ``p0``/``sigma`` construction as the point estimate.
Replicates whose window is too thin or whose fit fails are counted
(``n_boot_failed``), not silently dropped, and the reported interval is a
percentile CI (the ``p_th`` distribution is skewed), not merely a std. Measuring
the same estimator the point fit uses is what keeps the interval honest: the
earlier bootstrap refit an *unweighted* model from the point estimate on a
*frozen* window and dropped failures, all of which biased the CI narrow.
"""

from collections.abc import Sequence
from dataclasses import dataclass, field, replace

import numpy as np
import numpy.typing as npt
from scipy.optimize import curve_fit

from erasure_qec.analysis.statistics import (
    Z_95,
    SweepPoint,
    per_round_estimate,
    per_round_p_l,
)

_FitParams = tuple[float, float, float, float, float]
_N_PARAMS = 5  # (p_th, nu, A, B, C)

# Minimum degrees of freedom (points - parameters) for a fit to be trusted.
# Fewer than this and the 5-parameter ansatz is fitting noise, not signal.
_MIN_DOF = 3

# Per-round p_L at or above this is saturated (approaching the 1/2 at-chance
# limit) and outside the local quadratic ansatz's validity, so it is excluded
# from both the crossing estimate and the fit window. The cut is on the
# *per-round* rate (the variable the finite-size collapse is fit in), NOT the
# shot-level rate: at large T a legitimate near-crossing point still has
# P_L_shot -> 1/2 (per-round ~0.1 over T~11 rounds gives P_L_shot ~0.45), so a
# shot-level cap wrongly deletes near-crossing large-d data and breaks
# high-threshold fits (r_e=0.98). Imprecise near-saturation points are instead
# downweighted by their large 1-sigma Wilson errors. In estimate_crossing a p
# is dropped if *any* distance saturates, keeping the saturated large-d points
# present so the ordering-inversion guard can reject an above-threshold p.
_SATURATION_CAP = 0.4

# nu optimiser bounds; a fit landing on either is unconstrained, not converged.
_NU_BOUNDS = (0.5, 6.0)

# Bootstrap CI percentiles (a 95% two-sided *percentile* interval).
_CI_PERCENTILES = (2.5, 97.5)


@dataclass(frozen=True)
class FitResult:
    """Outcome of a finite-size scaling fit (§10)."""

    converged: bool
    p_th: float
    p_th_err: float
    nu: float
    nu_err: float
    coeffs: tuple[float, float, float]  # (A, B, C)
    p_center: float
    window_factor: float
    distances: tuple[int, ...]
    n_points: int
    d_min: int | None = None  # if set, only distances >= d_min were fit
    chi2_dof: float = float("nan")  # reduced chi-squared of the fit
    # 95% bootstrap *percentile* CIs (skewed distribution -> asymmetric, and a
    # better summary than +/- one std); p_th_err/nu_err keep the bootstrap std.
    p_th_ci: tuple[float, float] = (float("nan"), float("nan"))
    nu_ci: tuple[float, float] = (float("nan"), float("nan"))
    n_boot: int = 0  # bootstrap replicates attempted
    n_boot_failed: int = 0  # replicates whose window was too thin or fit failed
    message: str = ""
    # Points actually used, as parallel arrays (for plotting the collapse).
    used_p: npt.NDArray[np.float64] = field(
        default_factory=lambda: np.array([], dtype=np.float64)
    )
    used_d: npt.NDArray[np.int64] = field(
        default_factory=lambda: np.array([], dtype=np.int64)
    )
    used_p_l: npt.NDArray[np.float64] = field(
        default_factory=lambda: np.array([], dtype=np.float64)
    )

    def scaling_x(self) -> npt.NDArray[np.float64]:
        """Scaling variable ``x = (p - p_th) d^(1/nu)`` for the used points."""
        return (self.used_p - self.p_th) * self.used_d.astype(float) ** (1.0 / self.nu)

    def collapse_curve(
        self, x: npt.NDArray[np.float64]
    ) -> npt.NDArray[np.float64]:
        """The fitted universal parabola ``A + B x + C x^2`` at ``x``."""
        a, b, c = self.coeffs
        return a + b * x + c * x * x


def _model(
    pd: tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]], *params: float
) -> npt.NDArray[np.float64]:
    p, d = pd
    p_th, nu, a, b, c = params
    x = (p - p_th) * d ** (1.0 / nu)
    result: npt.NDArray[np.float64] = a + b * x + c * x * x
    return result


# A candidate crossing must sit at or below the below->above transition. Above
# threshold the curves re-converge on their way to the 1/2 at-chance limit,
# producing a *second*, spurious relative-spread minimum. This tolerance is how
# much larger the big-d curves must be than the small-d curves before the
# ordering counts as inverted (above threshold); small enough to keep the
# crossing itself (halves ~equal), large enough to reject the saturating tail.
_INVERSION_TOL = 0.10


def _ordering_inverted(pl_by_d: dict[int, float]) -> bool:
    """True if per-round ``p_L`` clearly *increases* with distance (above thr).

    Below threshold larger codes have smaller ``p_L`` (curves fan downward);
    above threshold the ordering inverts (larger codes do worse). Compares the
    mean ``p_L`` of the largest half of distances to the smallest half; a
    ``> _INVERSION_TOL`` relative excess flags the inverted, above-threshold
    ordering. Robust to single-point noise via the half-means.
    """
    ds = sorted(pl_by_d)
    half = len(ds) // 2
    if half == 0:
        return False
    lower = sum(pl_by_d[d] for d in ds[:half]) / half
    upper = sum(pl_by_d[d] for d in ds[-half:]) / half
    return upper > lower * (1.0 + _INVERSION_TOL)


def estimate_crossing(points: Sequence[SweepPoint]) -> float:
    """Estimate the threshold crossing as the p where the p_L curves coincide.

    At threshold the per-round ``p_L`` becomes distance-independent, so the
    curves for different ``d`` cross and the spread across distances is
    minimized *there*. Three robustness fixes over a naive ``max - min`` search:

    1. **Exclude saturated points** (per-round ``p_L >= _SATURATION_CAP``): a
       ``p`` where any distance has reached the ~1/2 at-chance limit is already
       above threshold and is skipped entirely, keeping the saturated
       large-d values present so the ordering-inversion guard below can see
       them.
    2. **Exclude inverted (above-threshold) orderings** (:func:`_ordering_inverted`):
       the estimator must sit at the transition between the below-threshold
       regime (``p_L`` decreasing in d) and the above-threshold regime
       (increasing in d). Above threshold the curves re-converge toward 1/2,
       creating a spurious relative-spread minimum that pulls the estimate past
       the crossing (e.g. the ragged r_e=0.5 tail). Any ``p`` whose ordering
       has already inverted is rejected.
    3. **Normalize by the mean** (relative spread ``(max - min) / mean``): raw
       ``max - min`` is trivially minimized in the deep sub-threshold tail,
       where every curve vanishes toward zero and the absolute spread is tiny
       — an artifact at the low edge of the p-grid. Dividing by the mean
       measures spread on the scale of the curves themselves.

    Falls back to the geometric-mean p when no p survives as a candidate.
    """
    by_p: dict[float, dict[int, float]] = {}
    for pt in points:
        by_p.setdefault(pt.p, {})[pt.d] = per_round_p_l(pt.p_l_shot, pt.rounds)

    relative_spread: dict[float, float] = {}
    for p, vals in by_p.items():
        values = list(vals.values())
        if len(values) < 2 or any(v >= _SATURATION_CAP for v in values):
            continue  # need >= 2 distances, all unsaturated, to compare
        if _ordering_inverted(vals):
            continue  # above threshold: p_L increases with d
        mean = sum(values) / len(values)
        if mean > 0.0:
            relative_spread[p] = (max(values) - min(values)) / mean

    if not relative_spread:
        ps = sorted(by_p)
        return float(np.sqrt(ps[0] * ps[-1])) if ps else float("nan")
    return min(relative_spread, key=lambda p: relative_spread[p])


def _select_window(
    points: Sequence[SweepPoint], p_center: float, window_factor: float
) -> list[SweepPoint]:
    lo, hi = p_center / window_factor, p_center * window_factor
    return [
        pt
        for pt in points
        if lo <= pt.p <= hi
        and per_round_p_l(pt.p_l_shot, pt.rounds) < _SATURATION_CAP
    ]


def fit_threshold(
    points: Sequence[SweepPoint],
    *,
    window_factor: float = 1.5,
    p_center: float | None = None,
    nu_guess: float = 1.5,
    n_boot: int = 200,
    seed: int = 0,
    d_min: int | None = None,
) -> FitResult:
    """Fit the finite-size scaling ansatz to ``points`` (one decoder/r_e).

    Args:
        d_min: if set, keep only points with ``d >= d_min`` before crossing
            estimation and window selection. Dropping the smallest distances
            (largest finite-size corrections) gives the asymptotic threshold;
            the effective crossing drifts upward when small-d curves are
            included (§10).

    Returns a :class:`FitResult`; ``converged=False`` (with a message and NaN
    parameters) when there is not enough data or the optimizer fails, so
    callers can degrade gracefully on partial sweeps.
    """
    if d_min is not None:
        points = [pt for pt in points if pt.d >= d_min]
    center = estimate_crossing(points) if p_center is None else p_center
    used = _select_window(points, center, window_factor)
    distances = tuple(sorted({pt.d for pt in used}))

    def _fail(msg: str, chi2_dof: float = float("nan")) -> FitResult:
        return FitResult(
            converged=False,
            p_th=float("nan"),
            p_th_err=float("nan"),
            nu=float("nan"),
            nu_err=float("nan"),
            coeffs=(float("nan"), float("nan"), float("nan")),
            p_center=center,
            window_factor=window_factor,
            distances=distances,
            n_points=len(used),
            d_min=d_min,
            chi2_dof=chi2_dof,
            message=msg,
        )

    # 5 free parameters; require >= _MIN_DOF degrees of freedom and >= 2
    # distances so the ansatz is constrained by data rather than fitting noise.
    if not _sufficient(used):
        return _fail(
            f"insufficient data in window: {len(used)} points, "
            f"{len(distances)} distances (need >= {_N_PARAMS + _MIN_DOF} points "
            f"for {_MIN_DOF} dof, >= 2 distances)"
        )

    core = _weighted_fit_arrays(used, center, nu_guess)
    if core is None:
        return _fail("curve_fit failed on the point estimate")
    popt = core.popt

    resid = (core.y - _model((core.p_arr, core.d_arr), *popt)) / core.sigma
    dof = len(used) - _N_PARAMS  # > 0 by the min_points gate above
    chi2_dof = float(np.sum(resid**2) / dof)

    # A parameter pinned at its optimiser bound means the data did not
    # constrain it: the crossing sits at/outside the fit window (p_th at the
    # p-range edge) or nu ran to a bound. Report these as *not converged* even
    # though curve_fit returned, so a reviewer is not misled by a bound value.
    p_span = float(core.p_arr.max() - core.p_arr.min())
    edge_tol = max(1e-3 * p_span, 1e-9)
    pinned: list[str] = []
    if abs(popt[0] - core.p_arr.min()) <= edge_tol:
        pinned.append("p_th at window minimum")
    elif abs(popt[0] - core.p_arr.max()) <= edge_tol:
        pinned.append("p_th at window maximum")
    if abs(popt[1] - _NU_BOUNDS[0]) <= 1e-3:
        pinned.append("nu at lower bound")
    elif abs(popt[1] - _NU_BOUNDS[1]) <= 1e-3:
        pinned.append("nu at upper bound")

    # Parametric bootstrap that re-runs the *whole* pipeline (resample all
    # points -> estimate_crossing -> window -> weighted fit), so the CI reflects
    # the same estimator the point fit uses, including crossing/window variance.
    boot, n_boot_failed = _bootstrap_pipeline(
        points, window_factor=window_factor, p_center=p_center,
        nu_guess=nu_guess, seed=seed, n_boot=n_boot,
    )
    if boot.size:
        lo_pct, hi_pct = _CI_PERCENTILES
        p_th_err = float(np.std(boot[:, 0]))
        nu_err = float(np.std(boot[:, 1]))
        p_th_ci = (
            float(np.percentile(boot[:, 0], lo_pct)),
            float(np.percentile(boot[:, 0], hi_pct)),
        )
        nu_ci = (
            float(np.percentile(boot[:, 1], lo_pct)),
            float(np.percentile(boot[:, 1], hi_pct)),
        )
    else:
        p_th_err = nu_err = float("nan")
        p_th_ci = nu_ci = (float("nan"), float("nan"))

    converged = not pinned
    message = (
        "ok"
        if converged
        else "not converged: " + "; ".join(pinned) + " pinned at optimiser bound"
    )
    return FitResult(
        converged=converged,
        p_th=popt[0],
        p_th_err=p_th_err,
        nu=popt[1],
        nu_err=nu_err,
        coeffs=(popt[2], popt[3], popt[4]),
        p_center=center,
        window_factor=window_factor,
        distances=distances,
        n_points=len(used),
        d_min=d_min,
        chi2_dof=chi2_dof,
        p_th_ci=p_th_ci,
        nu_ci=nu_ci,
        n_boot=n_boot,
        n_boot_failed=n_boot_failed,
        message=message,
        used_p=core.p_arr,
        used_d=core.d_arr.astype(np.int64),
        used_p_l=core.y,
    )


@dataclass(frozen=True)
class _CoreFit:
    """One weighted curve_fit and the arrays it ran on (for chi^2 / plotting)."""

    popt: _FitParams
    p_arr: npt.NDArray[np.float64]
    d_arr: npt.NDArray[np.float64]
    y: npt.NDArray[np.float64]
    sigma: npt.NDArray[np.float64]


def _sufficient(used: Sequence[SweepPoint]) -> bool:
    """Whether ``used`` can constrain the 5-parameter ansatz with >= _MIN_DOF."""
    return (
        len(used) >= _N_PARAMS + _MIN_DOF
        and len({pt.d for pt in used}) >= 2
    )


def _weighted_fit_arrays(
    used: Sequence[SweepPoint], center: float, nu_guess: float
) -> _CoreFit | None:
    """Build the weighted-LSQ arrays for ``used`` and run ``curve_fit`` once.

    This is the *single* definition of the ``p0`` / ``sigma`` / ``bounds``
    construction, shared by the point estimate and every bootstrap replicate so
    they measure the same (weighted) estimator. Returns ``None`` if the fit
    raises. The caller must have already checked :func:`_sufficient`.
    """
    p_arr = np.array([pt.p for pt in used], dtype=float)
    d_arr = np.array([pt.d for pt in used], dtype=float)
    y = np.array([per_round_p_l(pt.p_l_shot, pt.rounds) for pt in used])
    # per_round_estimate is a 95% Wilson interval, so its half-width is ~Z_95
    # standard errors; divide by Z_95 to recover a ~1-sigma weight for the fit.
    sigma = np.array(
        [
            max((e.high - e.low) / (2.0 * Z_95), 1e-9)
            for e in map(per_round_estimate, used)
        ]
    )
    a0 = float(np.median(y))
    # The p_th initial guess must be inside the p-range of the windowed points
    # (curve_fit rejects an out-of-bounds p0). The data-driven ``center`` can
    # fall outside it after saturated points are dropped, so clamp it; if the
    # crossing truly sits above every retained point the fit then pins p_th at
    # the bound and is reported not-converged rather than crashing.
    p_th0 = float(min(max(center, p_arr.min()), p_arr.max()))
    p0: _FitParams = (p_th0, nu_guess, a0, 0.0, 0.0)
    bounds = (
        [p_arr.min(), _NU_BOUNDS[0], -1.0, -1e6, -1e6],
        [p_arr.max(), _NU_BOUNDS[1], 1.0, 1e6, 1e6],
    )
    try:
        popt, _ = curve_fit(
            _model, (p_arr, d_arr), y, p0=p0, sigma=sigma,
            absolute_sigma=True, bounds=bounds, maxfev=20000,
        )
    except (RuntimeError, ValueError):
        return None
    params: _FitParams = (
        float(popt[0]), float(popt[1]), float(popt[2]), float(popt[3]), float(popt[4]),
    )
    return _CoreFit(params, p_arr, d_arr, y, sigma)


def _bootstrap_pipeline(
    points: Sequence[SweepPoint],
    *,
    window_factor: float,
    p_center: float | None,
    nu_guess: float,
    seed: int,
    n_boot: int,
) -> tuple[npt.NDArray[np.float64], int]:
    """Full-pipeline parametric bootstrap of ``(p_th, nu)``.

    Each replicate resamples ``errors ~ Binomial(shots, P_L_shot)`` over **all**
    ``points`` (the same d-filtered set the point estimate saw, not the frozen
    window), then re-runs :func:`estimate_crossing`, :func:`_select_window`, and
    the weighted fit exactly as the point estimate does. Returns the
    ``(n_ok, 2)`` array of recovered ``(p_th, nu)`` and the number of replicates
    that failed (thin window or non-converging fit) -- reported, never dropped.
    """
    rng = np.random.default_rng(seed)
    shots = np.array([pt.shots for pt in points])
    p_shot = np.array([pt.p_l_shot for pt in points])
    out: list[tuple[float, float]] = []
    n_failed = 0
    for _ in range(n_boot):
        drawn = rng.binomial(shots, p_shot)
        resampled = [
            replace(pt, errors=int(e))
            for pt, e in zip(points, drawn, strict=True)
        ]
        center = (
            estimate_crossing(resampled) if p_center is None else p_center
        )
        used = _select_window(resampled, center, window_factor)
        if not _sufficient(used):
            n_failed += 1
            continue
        core = _weighted_fit_arrays(used, center, nu_guess)
        if core is None:
            n_failed += 1
            continue
        out.append((core.popt[0], core.popt[1]))
    arr = np.array(out, dtype=float) if out else np.empty((0, 2))
    return arr, n_failed
