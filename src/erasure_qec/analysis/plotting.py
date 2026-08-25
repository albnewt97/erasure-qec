"""Deterministic figure generation from sweep CSVs (PLAN.md §10).

Every figure is byte-reproducible: the Agg backend is forced, a fixed rcParams
block pins fonts/sizes, all bootstrap CIs use a fixed seed, and PNG metadata
(the ``Software``/``Creation Time`` chunks matplotlib writes by default) is
stripped. Regenerating from the same CSVs yields identical bytes — the M8 gate.

Figures (§10):
  (i)  threshold panels — per-round ``p_L`` vs ``p``, one curve per ``d``, one
       panel per ``r_e``, the fitted ``p_th`` marked with its CI band;
  (ii) a data-collapse inset inside each panel (``p_L`` vs the scaling variable);
  (iii) the Lambda factor ``p_L(d)/p_L(d+2)`` vs ``p``;
  (iv) the herald_mwpm vs blind_mwpm ablation on identical tasks;
  (v)  the hook-regression figure — ``shortest_graphlike_error`` length for the
       correct vs broken CX schedule vs ``d`` (computed directly from the M3
       fixtures; no Monte Carlo, no CSV).

All CSV-fed figures degrade gracefully on partial sweeps: missing conditions,
too few points to fit, or a single ``p`` all render what data exists rather
than crashing.
"""

import argparse
from collections.abc import Sequence
from contextlib import AbstractContextManager
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.axes import Axes  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402

from erasure_qec.analysis.statistics import (  # noqa: E402
    SweepPoint,
    available_conditions,
    group_curves,
    lambda_factor,
    load_sweep,
    per_round_estimate,
)
from erasure_qec.analysis.threshold_fit import FitResult, fit_threshold  # noqa: E402

BOOTSTRAP_SEED = 0
_PNG_METADATA = {"Software": None, "Creation Time": None}

_RC = {
    "figure.dpi": 110,
    "savefig.dpi": 110,
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans"],
    "font.size": 9,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "grid.linewidth": 0.5,
    "lines.linewidth": 1.3,
    "lines.markersize": 4,
    "legend.frameon": False,
    "axes.spines.top": False,
    "axes.spines.right": False,
}


def _styled() -> AbstractContextManager[None]:
    """The deterministic rcParams context (one place for the matplotlib cast)."""
    return plt.rc_context(_RC)  # type: ignore[arg-type]


def _color_for_distance(d: int, all_d: Sequence[int]) -> tuple[float, float, float, float]:
    """Stable viridis color keyed on a distance's rank among all distances."""
    order = sorted(set(all_d))
    frac = order.index(d) / max(len(order) - 1, 1)
    return matplotlib.colormaps["viridis"](0.1 + 0.8 * frac)


def _save(fig: Figure, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, format="png", metadata=_PNG_METADATA)
    plt.close(fig)
    return path


def _plot_curves_on(ax: Axes, points: Sequence[SweepPoint]) -> list[int]:
    """Plot per-round ``p_L`` vs ``p`` curves (one per d) with Wilson CIs."""
    curves = group_curves(points, points[0].decoder, points[0].r_e)
    all_d = sorted(curves)
    for d in all_d:
        pts = curves[d]
        ps = np.array([pt.p for pt in pts])
        est = [per_round_estimate(pt) for pt in pts]
        y = np.array([e.value for e in est])
        lo = np.array([e.low for e in est])
        hi = np.array([e.high for e in est])
        color = _color_for_distance(d, all_d)
        ax.plot(ps, y, marker="o", color=color, label=f"d={d}")
        ax.fill_between(ps, lo, hi, color=color, alpha=0.2, linewidth=0)
    ax.set_xscale("log")
    ax.set_yscale("log")
    return all_d


# Smallest-distance cutoff for the asymptotic threshold fit: dropping small d
# (which carry the largest finite-size corrections) gives the asymptotic p_th.
_ASYMPTOTIC_D_MIN = 7


def _choose_threshold_fit(
    cond_points: Sequence[SweepPoint], *, window_factor: float, seed: int
) -> FitResult:
    """Fit with ``d >= 7`` when >= 3 such distances exist; else fall back to all d."""
    n_large = len({pt.d for pt in cond_points if pt.d >= _ASYMPTOTIC_D_MIN})
    if n_large >= 3:
        return fit_threshold(
            cond_points, window_factor=window_factor, seed=seed, d_min=_ASYMPTOTIC_D_MIN
        )
    return fit_threshold(cond_points, window_factor=window_factor, seed=seed)


def _mark_threshold(ax: Axes, fit: FitResult) -> None:
    # Gate on `resolved`, not `converged`: a fit that converges but whose CI is
    # wider than its own threshold (r_e=0.98 blind) must not be drawn as a line.
    if not fit.resolved:
        ax.text(
            0.03, 0.03, f"no fit:\n{fit.message.splitlines()[0]}",
            transform=ax.transAxes, fontsize=6, va="bottom", ha="left", color="0.4",
        )
        return
    ax.axvline(fit.p_th, color="crimson", linestyle="--", linewidth=1.0)
    lo, hi = fit.p_th_ci
    has_ci = np.isfinite(lo) and np.isfinite(hi) and hi > lo
    if has_ci:
        # Asymmetric 95% bootstrap percentile band (the p_th distribution is
        # skewed, so this is the honest interval, not p_th +/- one std).
        ax.axvspan(lo, hi, color="crimson", alpha=0.12)
    scope = f"$d\\geq{fit.d_min}$ fit" if fit.d_min is not None else "all-$d$ fit"
    ci_line = (
        f"\n95% CI [{lo*100:.2f}, {hi*100:.2f}]" if has_ci else ""
    )
    ax.text(
        0.97, 0.03,
        f"{scope}\n$p_{{th}}={fit.p_th*100:.2f}\\%${ci_line}\n$\\nu={fit.nu:.2f}$",
        transform=ax.transAxes, fontsize=6.5, va="bottom", ha="right", color="crimson",
    )


def _add_collapse_inset(ax: Axes, fit: FitResult) -> None:
    if not fit.resolved:
        return
    inset = ax.inset_axes((0.60, 0.60, 0.38, 0.38))
    x = fit.scaling_x()
    inset.scatter(x, fit.used_p_l, s=6, c="0.3", zorder=3)
    grid = np.linspace(float(x.min()), float(x.max()), 100)
    inset.plot(grid, fit.collapse_curve(grid), color="crimson", linewidth=1.0)
    inset.set_title("collapse", fontsize=6)
    inset.set_xlabel(r"$(p-p_{th})\,d^{1/\nu}$", fontsize=6)
    inset.tick_params(labelsize=5)


def figure_threshold_panels(
    points: Sequence[SweepPoint],
    out_path: Path,
    *,
    decoder: str = "herald_mwpm",
    window_factor: float = 1.5,
    seed: int = BOOTSTRAP_SEED,
) -> Path | None:
    """(i)+(ii): one panel per r_e for ``decoder``, each with a collapse inset.

    ``window_factor`` matches ``fit_threshold``'s default so the panel
    annotations equal the thresholds quoted in the README.
    """
    conditions = [(dec, re) for dec, re in available_conditions(points) if dec == decoder]
    if not conditions:
        return None
    with _styled():
        fig, axes = plt.subplots(
            1, len(conditions), figsize=(4.2 * len(conditions), 3.6),
            squeeze=False, layout="constrained",
        )
        for ax, (dec, r_e) in zip(axes[0], conditions, strict=True):
            cond_points = [pt for pt in points if pt.decoder == dec and pt.r_e == r_e]
            _plot_curves_on(ax, cond_points)
            fit = _choose_threshold_fit(
                cond_points, window_factor=window_factor, seed=seed
            )
            _mark_threshold(ax, fit)
            _add_collapse_inset(ax, fit)
            ax.set_title(f"{dec}, $r_e={r_e:g}$")
            ax.set_xlabel("physical error rate $p$")
            ax.set_ylabel("per-round logical error rate $p_L$")
            ax.legend(loc="upper left", fontsize=7)
        fig.suptitle("Threshold scaling (per-round $p_L$ vs $p$)", fontsize=11)
        return _save(fig, out_path)


# A Lambda = p_L(d)/p_L(d+2) point is only trustworthy when BOTH distances have
# enough logical errors: at low p / large d / high r_e the ratio divides a
# near-zero by a near-zero (e.g. 2 errors per 1e5 shots), so the bootstrap CI
# explodes (bars reaching ~200 around means near 50). Require this many observed
# logical errors in each distance before plotting the point.
_LAMBDA_MIN_ERRORS = 50


def figure_lambda(
    points: Sequence[SweepPoint],
    out_path: Path,
    *,
    decoder: str = "herald_mwpm",
    seed: int = BOOTSTRAP_SEED,
) -> Path | None:
    """(iii): Lambda = p_L(d)/p_L(d+2) vs p for each consecutive distance pair.

    Only points where both distances have >= ``_LAMBDA_MIN_ERRORS`` observed
    logical errors are plotted; the low-statistics tail (whose bootstrap CIs
    swamp the signal) is dropped.
    """
    conditions = [(dec, re) for dec, re in available_conditions(points) if dec == decoder]
    if not conditions:
        return None
    with _styled():
        fig, ax = plt.subplots(figsize=(5.0, 3.8), layout="constrained")
        plotted = False
        for _, r_e in conditions:
            curves = group_curves(points, decoder, r_e)
            all_d = sorted(curves)
            for d in all_d:
                if d + 2 not in curves:
                    continue
                by_p = {pt.p: pt for pt in curves[d]}
                by_p2 = {pt.p: pt for pt in curves[d + 2]}
                shared = sorted(
                    p
                    for p in set(by_p) & set(by_p2)
                    if by_p[p].errors >= _LAMBDA_MIN_ERRORS
                    and by_p2[p].errors >= _LAMBDA_MIN_ERRORS
                )
                if not shared:
                    continue
                lam = [lambda_factor(by_p[p], by_p2[p], seed=seed) for p in shared]
                yerr = np.array([[max(la.value - la.low, 0) for la in lam],
                                 [max(la.high - la.value, 0) for la in lam]])
                ax.errorbar(
                    shared, [la.value for la in lam], yerr=yerr, marker="s",
                    capsize=2, label=f"$r_e={r_e:g}$, $d{{=}}{d}\\!\\to\\!{d + 2}$",
                )
                plotted = True
        if not plotted:
            return None
        ax.axhline(1.0, color="0.6", linewidth=0.8, linestyle=":")
        ax.set_xscale("log")
        ax.set_xlabel("physical error rate $p$")
        ax.set_ylabel(r"$\Lambda = p_L(d)/p_L(d{+}2)$")
        ax.set_title(f"Distance suppression factor ({decoder})")
        ax.text(
            0.02, 0.98,
            f"points require $\\geq {_LAMBDA_MIN_ERRORS}$ logical errors in both $d$",
            transform=ax.transAxes, fontsize=6.5, va="top", ha="left", color="0.4",
        )
        ax.legend(fontsize=7)
        return _save(fig, out_path)


def figure_ablation(
    points: Sequence[SweepPoint],
    out_path: Path,
    *,
    seed: int = BOOTSTRAP_SEED,
) -> Path | None:
    """(iv): herald_mwpm vs blind_mwpm on identical tasks (per-round p_L vs p)."""
    decoders = sorted({pt.decoder for pt in points})
    if "herald_mwpm" not in decoders or "blind_mwpm" not in decoders:
        return None
    r_es = sorted({pt.r_e for pt in points if pt.r_e > 0.0}) or sorted({pt.r_e for pt in points})
    r_e = r_es[-1]
    curves_h = group_curves(points, "herald_mwpm", r_e)
    curves_b = group_curves(points, "blind_mwpm", r_e)
    shared_d = sorted(set(curves_h) & set(curves_b))
    if not shared_d:
        return None
    with _styled():
        fig, ax = plt.subplots(figsize=(5.2, 3.9), layout="constrained")
        for d in shared_d:
            color = _color_for_distance(d, shared_d)
            for pts, style, lab in (
                (curves_h[d], "-", "herald"), (curves_b[d], "--", "blind")
            ):
                ps = np.array([pt.p for pt in pts])
                y = np.array([per_round_estimate(pt).value for pt in pts])
                ax.plot(
                    ps, y, linestyle=style, marker="o" if style == "-" else "x",
                    color=color, label=f"d={d} {lab}",
                )
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("physical error rate $p$")
        ax.set_ylabel("per-round logical error rate $p_L$")
        ax.set_title(f"Herald-aware vs blind ablation ($r_e={r_e:g}$)")
        ax.legend(fontsize=6.5, ncol=2)
        return _save(fig, out_path)


def figure_hook_regression(
    out_path: Path, *, distances: Sequence[int] = (3, 5, 7)
) -> Path:
    """(v): shortest_graphlike_error length, correct vs broken schedule vs d.

    Computed directly from the M3 fixtures (imported, never modified): a small
    uniform depolarizing sprinkle gives the circuit graphlike error mechanisms,
    then ``shortest_graphlike_error`` is read off. The correct hook-safe
    schedule tracks the code distance d; the broken schedule collapses to
    ~ceil((d+1)/2), proving the schedule is load-bearing. No Monte Carlo.
    """
    from erasure_qec.circuits.builder import build
    from erasure_qec.circuits.scheduling import (
        BROKEN_SCHEDULE,
        HOOK_SAFE_SCHEDULE,
        Schedule,
    )
    from erasure_qec.config import NoiseParams
    from erasure_qec.noise.injector import PauliOnlyInjector

    def _sge(d: int, schedule: Schedule) -> int:
        inj = PauliOnlyInjector(NoiseParams(p=1e-3))
        return len(build(d, d, inj, schedule=schedule).shortest_graphlike_error())

    ds = list(distances)
    correct = [_sge(d, HOOK_SAFE_SCHEDULE) for d in ds]
    broken = [_sge(d, BROKEN_SCHEDULE) for d in ds]
    with _styled():
        fig, ax = plt.subplots(figsize=(4.6, 3.6), layout="constrained")
        ax.plot(ds, ds, color="0.6", linestyle=":", label="$d$ (ideal)")
        ax.plot(ds, correct, marker="o", color="seagreen", label="hook-safe schedule")
        ax.plot(ds, broken, marker="s", color="crimson", label="broken schedule")
        ax.set_xticks(ds)
        ax.set_xlabel("code distance $d$")
        ax.set_ylabel("shortest graphlike error length")
        ax.set_title("Hook regression: schedule is load-bearing")
        ax.legend()
        return _save(fig, out_path)


def render_all(
    data_dir: str | Path,
    figures_dir: str | Path,
    *,
    config_names: Sequence[str] = ("baseline_pauli", "erasure_r50", "erasure_r98"),
    seed: int = BOOTSTRAP_SEED,
) -> list[Path]:
    """Regenerate every figure from the CSVs found under ``data_dir``.

    Reads ``<name>.csv`` (and any ``<name>_lambda.csv``) for each config name,
    pooling all points so the threshold panels show every available r_e. The
    hook-regression figure needs no CSV. Returns the written figure paths.
    """
    data = Path(data_dir)
    figures = Path(figures_dir)
    csvs = [data / f"{n}.csv" for n in config_names]
    csvs += [data / f"{n}_lambda.csv" for n in config_names]
    points = load_sweep([c for c in csvs if c.exists()])

    written: list[Path] = []
    if points:
        csv_figures = [
            ("threshold_panels.png",
             figure_threshold_panels(points, figures / "threshold_panels.png", seed=seed)),
            ("lambda_vs_p.png",
             figure_lambda(points, figures / "lambda_vs_p.png", seed=seed)),
            ("ablation.png",
             figure_ablation(points, figures / "ablation.png", seed=seed)),
        ]
        for name, result in csv_figures:
            if result is not None:
                written.append(result)
            else:
                print(f"skipped {name}: insufficient data")
    else:
        print(f"no sweep CSVs found under {data}; only the hook figure will render")
    written.append(figure_hook_regression(figures / "hook_regression.png"))
    return written


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=Path("data"))
    parser.add_argument("--figures-dir", type=Path, default=Path("figures"))
    parser.add_argument("--seed", type=int, default=BOOTSTRAP_SEED)
    args = parser.parse_args()
    paths = render_all(args.data_dir, args.figures_dir, seed=args.seed)
    for p in paths:
        print(f"wrote {p}")


if __name__ == "__main__":
    main()
