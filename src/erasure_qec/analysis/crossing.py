"""Pairwise crossing of two distances' logical-error curves.

An estimator of the threshold that is independent of the finite-size collapse
fit in :mod:`erasure_qec.analysis.threshold_fit`: for two distances ``d1`` and
``d2`` it locates the physical error rate ``p`` at which the two curves of the
chosen ``metric`` (see :data:`~erasure_qec.analysis.statistics.Metric`) meet.

With ``g(p) = ln m(p, d2) - ln m(p, d1)``, below threshold the larger code has
the smaller logical rate, so ``g < 0`` when ``d2 > d1``; the crossing is where
``g`` first becomes non-negative. Between the two bracketing grid points the
crossing is located by linear interpolation of ``g`` in ``ln p``. No model of
the curves is fitted and no saturation cut is applied.
"""

import math
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

from erasure_qec.analysis.statistics import Metric, SweepPoint, metric_value

_CI_PERCENTILES = (2.5, 97.5)


@dataclass(frozen=True)
class CrossingEstimate:
    """The first crossing of the ``d1`` and ``d2`` curves of ``metric``.

    ``value`` is NaN and ``bracket`` is ``None`` when ``g`` never goes from
    negative to non-negative. ``bracket`` is the pair of adjacent grid values of
    ``p`` the interpolation used. ``n_sign_changes`` counts every change of sign
    class (negative vs non-negative) between adjacent usable grid points, in
    either direction; a value above 1 means the curves cross more than once on
    the grid and only the first crossing is reported.
    """

    value: float
    bracket: tuple[float, float] | None
    n_sign_changes: int
    metric: str
    d1: int
    d2: int


def _validated_curves(
    points: Sequence[SweepPoint], d1: int, d2: int
) -> tuple[dict[float, SweepPoint], dict[float, SweepPoint]]:
    """Points at ``d1`` and ``d2`` keyed by ``p``; raises on ambiguous input."""
    if d1 == d2:
        raise ValueError(f"d1 and d2 must differ, got {d1}")
    conditions = {(pt.decoder, pt.r_e) for pt in points}
    if len(conditions) > 1:
        raise ValueError(f"points must share one (decoder, r_e), got {sorted(conditions)}")
    curves: dict[int, dict[float, SweepPoint]] = {d1: {}, d2: {}}
    for pt in points:
        if pt.d in curves:
            if pt.p in curves[pt.d]:
                raise ValueError(f"duplicate point at d={pt.d}, p={pt.p}")
            curves[pt.d][pt.p] = pt
    return curves[d1], curves[d2]


def _crossing_from_curves(
    c1: dict[float, SweepPoint],
    c2: dict[float, SweepPoint],
    d1: int,
    d2: int,
    metric: Metric,
) -> CrossingEstimate:
    ps: list[float] = []
    gs: list[float] = []
    for p in sorted(set(c1) & set(c2)):
        m1 = metric_value(c1[p], metric)
        m2 = metric_value(c2[p], metric)
        if m1 > 0.0 and m2 > 0.0:
            ps.append(p)
            gs.append(math.log(m2) - math.log(m1))

    nonneg = [g >= 0.0 for g in gs]
    n_sign_changes = sum(a != b for a, b in zip(nonneg, nonneg[1:], strict=False))
    for i in range(len(gs) - 1):
        if gs[i] < 0.0 <= gs[i + 1]:
            lp0, lp1 = math.log(ps[i]), math.log(ps[i + 1])
            lp = lp0 - gs[i] * (lp1 - lp0) / (gs[i + 1] - gs[i])
            return CrossingEstimate(
                value=math.exp(lp),
                bracket=(ps[i], ps[i + 1]),
                n_sign_changes=n_sign_changes,
                metric=metric,
                d1=d1,
                d2=d2,
            )
    return CrossingEstimate(
        value=float("nan"),
        bracket=None,
        n_sign_changes=n_sign_changes,
        metric=metric,
        d1=d1,
        d2=d2,
    )


def pairwise_crossing(
    points: Sequence[SweepPoint], d1: int, d2: int, *, metric: Metric
) -> CrossingEstimate:
    """First crossing of the ``d1`` and ``d2`` curves of ``metric``.

    ``points`` must belong to one ``(decoder, r_e)`` condition and hold at most
    one point per ``(d, p)``; points at other distances are ignored. Only the
    ``p`` values present at both distances with a non-zero metric value at both
    are used (``ln`` of a zero rate is undefined), sorted ascending. With
    ``g(p) = ln m(d2) - ln m(d1)``, the estimate is the first adjacent pair
    ``(p_i, p_{i+1})`` with ``g_i < 0 <= g_{i+1}``, interpolated linearly in
    ``(ln p, g)`` to ``g = 0``. Returns ``value = NaN`` and ``bracket = None``
    if there is no such pair.
    """
    c1, c2 = _validated_curves(points, d1, d2)
    return _crossing_from_curves(c1, c2, d1, d2, metric)


def bootstrap_pairwise_crossing(
    points: Sequence[SweepPoint],
    d1: int,
    d2: int,
    *,
    metric: Metric,
    n_boot: int = 1000,
    seed: int = 0,
) -> tuple[CrossingEstimate, tuple[float, float], int]:
    """:func:`pairwise_crossing` with a parametric-bootstrap 95% percentile CI.

    Each replicate resamples ``errors ~ Binomial(shots, errors / shots)``
    independently for every point at ``d1`` and ``d2`` (one
    ``numpy.random.default_rng(seed)`` stream; points are visited in ``(d, p)``
    order so the result does not depend on the input order) and recomputes the
    crossing. The CI is the 2.5th/97.5th percentile over the replicates that
    have a crossing; the others are counted and returned as ``n_no_crossing``.
    If no replicate has a crossing the CI is ``(nan, nan)``.

    The CI reflects binomial sampling noise only. It does NOT include the
    interpolation error inside the bracket: the true curves need not be linear
    in ``(ln p, g)`` between the two grid points, and that error is not
    estimated here.

    Returns ``(point_estimate, ci, n_no_crossing)``, where ``point_estimate`` is
    :func:`pairwise_crossing` on the unresampled data.
    """
    c1, c2 = _validated_curves(points, d1, d2)
    estimate = _crossing_from_curves(c1, c2, d1, d2, metric)

    ordered = sorted([*c1.values(), *c2.values()], key=lambda pt: (pt.d, pt.p))
    shots = np.array([pt.shots for pt in ordered], dtype=np.int64)
    rates = np.array([pt.p_l_shot for pt in ordered], dtype=float)

    rng = np.random.default_rng(seed)
    values: list[float] = []
    n_no_crossing = 0
    for _ in range(n_boot):
        drawn = rng.binomial(shots, rates)
        r1: dict[float, SweepPoint] = {}
        r2: dict[float, SweepPoint] = {}
        for pt, e in zip(ordered, drawn, strict=True):
            target = r1 if pt.d == d1 else r2
            target[pt.p] = SweepPoint(pt.decoder, pt.d, pt.rounds, pt.p, pt.r_e, pt.shots, int(e))
        rep = _crossing_from_curves(r1, r2, d1, d2, metric)
        if math.isnan(rep.value):
            n_no_crossing += 1
        else:
            values.append(rep.value)

    if values:
        lo_pct, hi_pct = _CI_PERCENTILES
        ci = (float(np.percentile(values, lo_pct)), float(np.percentile(values, hi_pct)))
    else:
        ci = (float("nan"), float("nan"))
    return estimate, ci, n_no_crossing
