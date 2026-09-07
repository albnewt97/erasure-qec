"""Shot statistics: per-round conversion, CIs, and the Lambda factor (PLAN.md §10).

Sinter reports a *shot-level* logical error probability ``P_L_shot`` (errors /
shots over a ``T``-round memory experiment). The physically comparable quantity
is the *per-round* logical error rate

    p_L = 1/2 * (1 - (1 - 2 * P_L_shot)^(1/T))          [PLAN.md §10]

which inverts the accumulation ``P_L_shot = 1/2 (1 - (1 - 2 p_L)^T)`` of an
error that must survive an odd number of the ``T`` rounds. Two fixed points
pin it: ``P_L_shot = 1/2 -> p_L = 1/2`` (fully scrambled) and ``T = 1 ->
p_L = P_L_shot`` (identity).

Confidence intervals come two ways, both exposed:
- **Wilson** score interval on the binomial ``P_L_shot``, its endpoints mapped
  through the (monotonic) per-round conversion — closed form, no sampling.
- **Bootstrap** — resample ``errors ~ Binomial(shots, P_L_shot)`` and push each
  draw through the conversion, then take percentiles. Seeded for determinism.

:func:`lambda_factor` bootstraps a *ratio*, so it additionally has to discard
replicates with a zero denominator. Those discards are counted and surfaced on
:class:`LambdaEstimate` rather than dropped silently — see that function's
docstring for why the conditioning biases the interval's upper endpoint low.
"""

from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import sinter

# Standard-normal 97.5th percentile: the z for a two-sided 95% interval.
Z_95 = 1.959963984540054


@dataclass(frozen=True)
class SweepPoint:
    """One decoded (decoder, d, p, r_e) task: raw shot counts from sinter."""

    decoder: str
    d: int
    rounds: int
    p: float
    r_e: float
    shots: int
    errors: int

    @property
    def p_l_shot(self) -> float:
        """Shot-level logical error probability, ``errors / shots``."""
        return self.errors / self.shots


@dataclass(frozen=True)
class Estimate:
    """A point estimate with a (low, high) confidence interval."""

    value: float
    low: float
    high: float


@dataclass(frozen=True)
class LambdaEstimate(Estimate):
    """:class:`Estimate` plus the bootstrap-conditioning bookkeeping of
    :func:`lambda_factor`.

    A separate type rather than three more fields on :class:`Estimate`: the
    Wilson interval from :func:`per_round_estimate` has no replicates at all, so
    an ``n_boot_dropped = 0`` there would read as "checked, none dropped" when
    nothing was ever checked. Subclassing keeps ``.value``/``.low``/``.high``
    and ``isinstance(x, Estimate)`` working for existing callers.

    Naming mirrors :class:`~erasure_qec.analysis.threshold_fit.FitResult`'s
    ``n_boot`` / ``n_boot_failed`` pair.
    """

    n_boot: int = 0  # replicates drawn
    n_boot_dropped: int = 0  # replicates with a zero denominator, excluded
    ci: float = 0.95  # the interval's nominal coverage, for `ci_conditioned`

    @property
    def dropped_fraction(self) -> float:
        """Fraction of replicates excluded from the CI (NaN if none drawn)."""
        return self.n_boot_dropped / self.n_boot if self.n_boot else float("nan")

    @property
    def ci_conditioned(self) -> bool:
        """True when the discards outweigh the tail the CI claims to measure.

        The dropped replicates are exactly those where the ``d+2`` code logged
        zero errors -- the ones whose unconditional ``Lambda`` is ``+inf``. Once
        more than ``(1 - ci) / 2`` of the draws are dropped, the *true* upper
        quantile is ``+inf``, and the reported ``high`` is not an approximation
        of it: it is the upper endpoint of a different distribution (``Lambda``
        given that ``d+2`` failed at least once). Treat such an interval the way
        ``FitResult.resolved = False`` is treated -- do not quote it.
        """
        return bool(self.n_boot) and self.dropped_fraction > (1.0 - self.ci) / 2.0


def per_round_p_l(p_l_shot: float, rounds: int) -> float:
    """Convert a shot-level ``P_L_shot`` to the per-round ``p_L`` (§10).

    ``P_L_shot`` is clamped to ``[0, 1/2]`` first: a decoder can only exceed
    1/2 on noise-limited samples, and ``1 - 2 P_L_shot`` must stay non-negative
    for the fractional power to be real. ``p_L`` is monotonically increasing in
    ``P_L_shot`` on this domain, which is what lets CI endpoints be mapped
    through directly.
    """
    if rounds < 1:
        raise ValueError(f"rounds must be >= 1, got {rounds}")
    x = min(max(p_l_shot, 0.0), 0.5)
    return float(0.5 * (1.0 - (1.0 - 2.0 * x) ** (1.0 / rounds)))


def shot_p_l_from_per_round(p_l: float, rounds: int) -> float:
    """Inverse of :func:`per_round_p_l`: accumulate a per-round rate over T rounds."""
    if rounds < 1:
        raise ValueError(f"rounds must be >= 1, got {rounds}")
    x = min(max(p_l, 0.0), 0.5)
    return 0.5 * (1.0 - (1.0 - 2.0 * x) ** rounds)


def wilson_interval(errors: int, shots: int, z: float = Z_95) -> tuple[float, float]:
    """Wilson score interval for a binomial proportion ``errors / shots``."""
    if shots <= 0:
        raise ValueError("shots must be positive")
    n = float(shots)
    phat = errors / n
    denom = 1.0 + z * z / n
    center = (phat + z * z / (2.0 * n)) / denom
    half = z * np.sqrt(phat * (1.0 - phat) / n + z * z / (4.0 * n * n)) / denom
    return (max(0.0, center - half), min(1.0, center + half))


def per_round_estimate(point: SweepPoint, z: float = Z_95) -> Estimate:
    """Per-round ``p_L`` with a Wilson CI mapped through the conversion."""
    lo, hi = wilson_interval(point.errors, point.shots, z)
    return Estimate(
        value=per_round_p_l(point.p_l_shot, point.rounds),
        low=per_round_p_l(lo, point.rounds),
        high=per_round_p_l(hi, point.rounds),
    )


def bootstrap_per_round(
    point: SweepPoint,
    *,
    n_boot: int = 2000,
    seed: int = 0,
    ci: float = 0.95,
) -> Estimate:
    """Per-round ``p_L`` with a bootstrap CI (seeded; §10)."""
    rng = np.random.default_rng(seed)
    draws = rng.binomial(point.shots, point.p_l_shot, size=n_boot) / point.shots
    vals = np.array([per_round_p_l(float(x), point.rounds) for x in draws])
    lo_q, hi_q = (1.0 - ci) / 2.0, 1.0 - (1.0 - ci) / 2.0
    return Estimate(
        value=per_round_p_l(point.p_l_shot, point.rounds),
        low=float(np.quantile(vals, lo_q)),
        high=float(np.quantile(vals, hi_q)),
    )


def lambda_factor(
    point_d: SweepPoint,
    point_d2: SweepPoint,
    *,
    n_boot: int = 2000,
    seed: int = 0,
    ci: float = 0.95,
) -> LambdaEstimate:
    """Suppression factor ``Lambda = p_L(d) / p_L(d+2)`` at fixed p, with a
    bootstrap CI (§10). ``point_d2`` must be the ``d + 2`` partner of
    ``point_d`` (same decoder, p, r_e).

    **The CI is conditioned, and the conditioning is not ignorable.** A
    replicate whose ``d+2`` draw logs zero errors has a zero denominator, so its
    ``Lambda`` is ``+inf`` (or ``0/0``); such replicates are excluded from the
    percentile CI. That is not missing-at-random -- it removes precisely the
    largest values, so ``high`` is biased *low*, never high. The count is
    reported as ``n_boot_dropped`` (cf. ``FitResult.n_boot_failed``, which the
    threshold fit reports for the same reason).

    Once the dropped fraction exceeds the CI's own tail mass ``(1 - ci) / 2``,
    the unconditional upper quantile is ``+inf`` and ``high`` is the endpoint of
    a *different* distribution -- ``Lambda`` given that ``d+2`` failed at least
    once. :attr:`LambdaEstimate.ci_conditioned` flags exactly that case; do not
    quote such an interval. Real example from ``data/erasure_r98.csv``:
    ``herald_mwpm``, ``r_e = 0.98``, ``d = 3 -> 5`` at ``p = 0.001`` has 1 error
    in 100k shots at ``d = 5``, drops 745/2000 replicates (37%), and reports
    ``[10, 56]`` for a quantity whose honest upper limit is unbounded.

    The dropped replicates are *excluded, not reweighted*: the estimator is left
    as-is so committed figures stay reproducible. The ``>= 50``-observed-error
    gate in :func:`~erasure_qec.analysis.plotting.figure_lambda` keeps every
    plotted point clear of this -- all 425 gated-in pairs in ``data/`` drop zero
    replicates -- but ``lambda_factor`` is public API and callers that skip the
    gate must check ``ci_conditioned`` themselves.
    """
    if point_d2.d != point_d.d + 2:
        raise ValueError(f"expected d+2={point_d.d + 2}, got d={point_d2.d}")
    num = per_round_p_l(point_d.p_l_shot, point_d.rounds)
    den = per_round_p_l(point_d2.p_l_shot, point_d2.rounds)
    rng = np.random.default_rng(seed)
    num_draws = rng.binomial(point_d.shots, point_d.p_l_shot, n_boot) / point_d.shots
    den_draws = rng.binomial(point_d2.shots, point_d2.p_l_shot, n_boot) / point_d2.shots
    ratios = []
    n_dropped = 0
    for a, b in zip(num_draws, den_draws, strict=True):
        pb = per_round_p_l(float(b), point_d2.rounds)
        if pb > 0.0:
            ratios.append(per_round_p_l(float(a), point_d.rounds) / pb)
        else:
            n_dropped += 1  # zero denominator: unconditional Lambda is +inf
    arr = np.array(ratios) if ratios else np.array([np.nan])
    lo_q, hi_q = (1.0 - ci) / 2.0, 1.0 - (1.0 - ci) / 2.0
    return LambdaEstimate(
        value=num / den if den > 0 else float("nan"),
        low=float(np.nanquantile(arr, lo_q)),
        high=float(np.nanquantile(arr, hi_q)),
        n_boot=n_boot,
        n_boot_dropped=n_dropped,
        ci=ci,
    )


def load_sweep(paths: str | Path | Sequence[str | Path]) -> list[SweepPoint]:
    """Load and merge sinter CSV(s) into ``SweepPoint`` records (§10 schema).

    Rows sharing a strong id (sinter's resume appends) are summed. Zero-shot
    rows and rows missing required metadata are skipped, so a partially-filled
    sweep loads cleanly.
    """
    if isinstance(paths, str | Path):
        paths = [paths]
    existing = [str(p) for p in paths if Path(p).exists()]
    if not existing:
        return []
    points: list[SweepPoint] = []
    for stat in sinter.read_stats_from_csv_files(*existing):
        meta = stat.json_metadata
        if stat.shots <= 0 or meta is None or "d" not in meta or "p" not in meta:
            continue
        d = int(meta["d"])
        points.append(
            SweepPoint(
                decoder=stat.decoder or "unknown",
                d=d,
                rounds=int(meta.get("rounds", d)),
                p=float(meta["p"]),
                r_e=float(meta.get("r_e", 0.0)),
                shots=int(stat.shots),
                errors=int(stat.errors),
            )
        )
    return points


def group_curves(
    points: Iterable[SweepPoint], decoder: str, r_e: float
) -> dict[int, list[SweepPoint]]:
    """Points for one (decoder, r_e), grouped by ``d`` and sorted by ``p``."""
    curves: dict[int, list[SweepPoint]] = {}
    for pt in points:
        if pt.decoder == decoder and pt.r_e == r_e:
            curves.setdefault(pt.d, []).append(pt)
    for d in curves:
        curves[d].sort(key=lambda pt: pt.p)
    return curves


def available_conditions(points: Iterable[SweepPoint]) -> list[tuple[str, float]]:
    """Sorted distinct ``(decoder, r_e)`` conditions present in the data."""
    return sorted({(pt.decoder, pt.r_e) for pt in points})
