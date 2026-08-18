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

Confidence on ``p_th`` comes from a seeded parametric bootstrap over the sinter
counts: resample ``errors ~ Binomial(shots, P_L_shot)`` per point, refit, and
report the mean and standard deviation of the recovered ``p_th``.
"""

from collections.abc import Sequence
from dataclasses import dataclass, field

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
    min_points = _N_PARAMS + _MIN_DOF
    if len(used) < min_points or len(distances) < 2:
        return _fail(
            f"insufficient data in window: {len(used)} points, "
            f"{len(distances)} distances (need >= {min_points} points "
            f"for {_MIN_DOF} dof, >= 2 distances)"
        )

    p_arr = np.array([pt.p for pt in used], dtype=float)
    d_arr = np.array([pt.d for pt in used], dtype=float)
    y = np.array([per_round_p_l(pt.p_l_shot, pt.rounds) for pt in used])
    # per_round_estimate is a 95% Wilson interval, so its half-width is ~Z_95
    # standard errors; divide by Z_95 to recover a ~1-sigma weight for the fit.
    # (The earlier code used the 95% half-width directly as 1 sigma, deflating
    # chi^2/dof by ~Z_95^2 ~ 3.84x; see docs/AUDIT.md.)
    sigma = np.array(
        [
            max((e.high - e.low) / (2.0 * Z_95), 1e-9)
            for e in map(per_round_estimate, used)
        ]
    )

    a0 = float(np.median(y))
    # The p_th initial guess must be inside the p-range of the *windowed* points
    # (curve_fit rejects an out-of-bounds p0). The data-driven ``center`` can
    # fall outside it after the coin-flip cut removes above-crossing points at
    # large d, so clamp it; if the crossing truly sits above every retained
    # point, the fit then pins p_th at the bound and is reported as not
    # converged rather than crashing.
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
    except (RuntimeError, ValueError) as exc:
        return _fail(f"curve_fit failed: {exc}")

    resid = (y - _model((p_arr, d_arr), *popt)) / sigma
    dof = len(used) - _N_PARAMS  # > 0 by the min_points gate above
    chi2_dof = float(np.sum(resid**2) / dof)

    # A parameter pinned at its optimiser bound means the data did not
    # constrain it: the crossing sits at/outside the fit window (p_th at the
    # p-range edge) or nu ran to a bound. Report these as *not converged* even
    # though curve_fit returned, so a reviewer is not misled by a bound value.
    p_span = float(p_arr.max() - p_arr.min())
    edge_tol = max(1e-3 * p_span, 1e-9)
    pinned: list[str] = []
    if abs(float(popt[0]) - p_arr.min()) <= edge_tol:
        pinned.append("p_th at window minimum")
    elif abs(float(popt[0]) - p_arr.max()) <= edge_tol:
        pinned.append("p_th at window maximum")
    if abs(float(popt[1]) - _NU_BOUNDS[0]) <= 1e-3:
        pinned.append("nu at lower bound")
    elif abs(float(popt[1]) - _NU_BOUNDS[1]) <= 1e-3:
        pinned.append("nu at upper bound")

    boot = _refit_bootstrap(used, popt, bounds, seed=seed, n_boot=n_boot)
    p_th_err = float(np.std(boot[:, 0])) if boot.size else float("nan")
    nu_err = float(np.std(boot[:, 1])) if boot.size else float("nan")

    converged = not pinned
    message = (
        "ok"
        if converged
        else "not converged: " + "; ".join(pinned) + " pinned at optimiser bound"
    )
    return FitResult(
        converged=converged,
        p_th=float(popt[0]),
        p_th_err=p_th_err,
        nu=float(popt[1]),
        nu_err=nu_err,
        coeffs=(float(popt[2]), float(popt[3]), float(popt[4])),
        p_center=center,
        window_factor=window_factor,
        distances=distances,
        n_points=len(used),
        d_min=d_min,
        chi2_dof=chi2_dof,
        message=message,
        used_p=p_arr,
        used_d=d_arr.astype(np.int64),
        used_p_l=y,
    )


def _refit_bootstrap(
    used: Sequence[SweepPoint],
    p0: npt.NDArray[np.float64],
    bounds: tuple[list[float], list[float]],
    *,
    seed: int,
    n_boot: int,
) -> npt.NDArray[np.float64]:
    """Parametric bootstrap: return the (n_boot, 5) array of refit parameters."""
    rng = np.random.default_rng(seed)
    p_arr = np.array([pt.p for pt in used], dtype=float)
    d_arr = np.array([pt.d for pt in used], dtype=float)
    shots = np.array([pt.shots for pt in used])
    p_shot = np.array([pt.p_l_shot for pt in used])
    rounds = np.array([pt.rounds for pt in used])
    out: list[npt.NDArray[np.float64]] = []
    for _ in range(n_boot):
        resampled = rng.binomial(shots, p_shot) / shots
        y = np.array(
            [per_round_p_l(float(v), int(t)) for v, t in zip(resampled, rounds, strict=True)]
        )
        try:
            popt, _ = curve_fit(
                _model, (p_arr, d_arr), y, p0=p0, bounds=bounds, maxfev=20000
            )
            out.append(popt)
        except (RuntimeError, ValueError):
            continue
    return np.array(out) if out else np.empty((0, 5))
