"""Deterministic figure generation from sweep CSVs (PLAN.md §10).

Every figure is byte-reproducible: the Agg backend is forced, a fixed rcParams
block pins fonts/sizes, all bootstrap CIs use a fixed seed, and PNG metadata
(the ``Software``/``Creation Time`` chunks matplotlib writes by default) is
stripped. Regenerating from the same inputs yields identical bytes — the M8 gate.

Figures (§10):
  (i)  threshold panels — the logical failure probability after ``T = d``
       rounds (per-shot ``P_L_shot``) vs ``p``, one panel per ``r_e``, one curve
       per ``(decoder, d)`` (herald solid, blind dashed), each decoder's
       headline per-shot (7,11) crossing
       (:func:`~erasure_qec.analysis.crossing.headline_crossing`) marked with
       its 95% bootstrap CI band; no mark where that crossing is missing or
       saturated;
  (ii) the per-round Lambda factor ``p_L(d)/p_L(d+2)`` vs ``p``;
  (iii) the herald_mwpm vs blind_mwpm suppression ratio at fixed ``p``, from
       ``results/ablation.json`` (``scripts/ablation_table.py``); skipped when
       that file does not exist;
  (iv) the hook-regression figure — ``shortest_graphlike_error`` length for the
       correct vs broken CX schedule vs ``d`` (computed directly from the M3
       fixtures; no Monte Carlo, no CSV).

All CSV-fed figures degrade gracefully on partial sweeps: missing conditions,
too few points for a crossing, or a single ``p`` all render what data exists
rather than crashing.
"""

import argparse
import json
import math
from collections.abc import Sequence
from contextlib import AbstractContextManager
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

# The imports below sit after matplotlib.use() deliberately -- the backend must be
# selected before pyplot is imported. They carry no E402 suppression: ruff exempts
# matplotlib.use() (along with sys.path.* and os.environ mutation) from E402, so a
# suppression here would be dead and RUF100 now rejects it. pycodestyle proper does
# not grant that exemption -- if the linter is ever switched, the suppression returns.
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure

from erasure_qec.analysis.crossing import HEADLINE_PAIR, headline_crossing
from erasure_qec.analysis.statistics import (
    SweepPoint,
    available_conditions,
    group_curves,
    lambda_factor,
    load_sweep,
    metric_estimate,
)

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


_DECODER_STYLE = {"herald_mwpm": ("-", "o", "crimson"), "blind_mwpm": ("--", "x", "0.25")}
_THRESHOLD_YLABEL = "logical failure probability after T = d rounds"


def _plot_shot_curves_on(ax: Axes, points: Sequence[SweepPoint], all_d: Sequence[int]) -> None:
    """Per-shot ``P_L_shot`` vs ``p`` with Wilson CIs, one curve per ``d``.

    Points with zero observed errors are not drawn (``P_L_shot = 0`` has no
    place on the log axis); the crossing estimator skips them as well.
    """
    linestyle, marker, _ = _DECODER_STYLE.get(points[0].decoder, ("-", "o", "k"))
    short = points[0].decoder.removesuffix("_mwpm")
    curves = group_curves(points, points[0].decoder, points[0].r_e)
    for d in sorted(curves):
        pts = [pt for pt in curves[d] if pt.errors > 0]
        if not pts:
            continue
        ps = np.array([pt.p for pt in pts])
        est = [metric_estimate(pt, "shot") for pt in pts]
        color = _color_for_distance(d, all_d)
        ax.plot(
            ps,
            [e.value for e in est],
            linestyle=linestyle,
            marker=marker,
            color=color,
            label=f"d={d} {short}",
        )
        ax.fill_between(
            ps, [e.low for e in est], [e.high for e in est], color=color, alpha=0.15, linewidth=0
        )


def _mark_headline_crossing(ax: Axes, points: Sequence[SweepPoint], *, seed: int) -> str:
    """Mark the headline (7,11) crossing of one decoder; return its annotation line."""
    est, (lo, hi), _ = headline_crossing(points, seed=seed)
    d1, d2 = HEADLINE_PAIR
    name = points[0].decoder.removesuffix("_mwpm")
    if est.saturated:
        return f"{name} ({d1},{d2}): saturated, not a threshold"
    if not math.isfinite(est.value):
        return f"{name} ({d1},{d2}): no crossing on the grid"
    linestyle, _, color = _DECODER_STYLE.get(points[0].decoder, ("-", "o", "k"))
    ax.axvline(est.value, color=color, linestyle=linestyle, linewidth=1.0)
    has_ci = math.isfinite(lo) and math.isfinite(hi) and hi > lo
    if has_ci:
        ax.axvspan(lo, hi, color=color, alpha=0.12)
    ci = f" [{lo * 100:.2f}, {hi * 100:.2f}]" if has_ci else ""
    return f"{name} ({d1},{d2}): {est.value * 100:.2f}%{ci}"


def figure_threshold_panels(
    points: Sequence[SweepPoint],
    out_path: Path,
    *,
    decoders: Sequence[str] = ("herald_mwpm", "blind_mwpm"),
    seed: int = BOOTSTRAP_SEED,
) -> Path | None:
    """(i): per-shot failure probability vs ``p``, one panel per ``r_e``.

    Each panel overlays ``decoders`` (herald solid, blind dashed) and marks each
    decoder's headline per-shot (7,11) crossing with its 95% bootstrap CI band
    (binomial noise only); the annotation in the panel corner gives the value or
    why there is none.
    """
    r_es = sorted({pt.r_e for pt in points if pt.decoder in decoders})
    if not r_es:
        return None
    all_d = sorted({pt.d for pt in points if pt.decoder in decoders})
    with _styled():
        fig, axes = plt.subplots(
            1, len(r_es), figsize=(4.4 * len(r_es), 3.8), squeeze=False, layout="constrained"
        )
        for ax, r_e in zip(axes[0], r_es, strict=True):
            notes = []
            for dec in decoders:
                cond = [pt for pt in points if pt.decoder == dec and pt.r_e == r_e]
                if not cond:
                    continue
                _plot_shot_curves_on(ax, cond, all_d)
                notes.append(_mark_headline_crossing(ax, cond, seed=seed))
            ax.text(
                0.97,
                0.03,
                "per-shot crossing, % [95% CI]\n" + "\n".join(notes),
                transform=ax.transAxes,
                fontsize=6,
                va="bottom",
                ha="right",
                color="0.2",
            )
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_title(f"$r_e={r_e:g}$")
            ax.set_xlabel("physical error rate $p$")
            ax.set_ylabel(_THRESHOLD_YLABEL)
            ax.legend(loc="upper left", fontsize=5.5, ncol=2)
        fig.suptitle(
            "Threshold (herald solid, blind dashed; lines: per-shot (7,11) crossing)",
            fontsize=10,
        )
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

    That gate is also what keeps :func:`lambda_factor`'s conditioned bootstrap
    out of this figure: replicates with a zero ``d+2`` denominator are excluded
    from its CI, which biases the upper bar low. On the committed sweeps no
    gated-in point drops a single replicate, but that is a property of the data,
    not a theorem -- so any point whose CI *is* conditioned
    (:attr:`LambdaEstimate.ci_conditioned`) is counted and called out in the
    corner note rather than drawn as if its error bar meant what it usually
    means.
    """
    conditions = [(dec, re) for dec, re in available_conditions(points) if dec == decoder]
    if not conditions:
        return None
    with _styled():
        fig, ax = plt.subplots(figsize=(5.0, 3.8), layout="constrained")
        plotted = False
        n_conditioned = 0  # plotted points whose bootstrap CI is conditioned
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
                n_conditioned += sum(la.ci_conditioned for la in lam)
                yerr = np.array(
                    [
                        [max(la.value - la.low, 0) for la in lam],
                        [max(la.high - la.value, 0) for la in lam],
                    ]
                )
                ax.errorbar(
                    shared,
                    [la.value for la in lam],
                    yerr=yerr,
                    marker="s",
                    capsize=2,
                    label=f"$r_e={r_e:g}$, $d{{=}}{d}\\!\\to\\!{d + 2}$",
                )
                plotted = True
        if not plotted:
            return None
        ax.axhline(1.0, color="0.6", linewidth=0.8, linestyle=":")
        ax.set_xscale("log")
        ax.set_xlabel("physical error rate $p$")
        ax.set_ylabel(r"per-round $\Lambda = p_L(d)/p_L(d{+}2)$")
        ax.set_title(f"Distance suppression factor ({decoder})")
        note = f"points require $\\geq {_LAMBDA_MIN_ERRORS}$ logical errors in both $d$"
        if n_conditioned:
            note += (
                f"\n{n_conditioned} point(s) have a conditioned bootstrap CI:"
                "\nupper bar is a lower bound on the true limit"
            )
        ax.text(
            0.02,
            0.98,
            note,
            transform=ax.transAxes,
            fontsize=6.5,
            va="top",
            ha="left",
            color="0.4",
        )
        ax.legend(fontsize=7)
        return _save(fig, out_path)


def figure_ablation(ablation_json: Path, out_path: Path) -> Path | None:
    """(iii): ``p_L(blind) / p_L(herald)`` vs ``d`` from ``results/ablation.json``.

    One line per ``(r_e, p)``. Cells at or above the 50-herald-error gate show
    the point ratio (filled markers); gated cells show the >= 95% one-sided
    lower confidence bound ``ratio_lb`` (open up-triangles), never a point
    ratio. Returns ``None`` when ``ablation_json`` does not exist.
    """
    if not ablation_json.exists():
        return None
    doc = json.loads(ablation_json.read_text())
    cells = doc["cells"]
    if not cells:
        return None
    r_es = sorted({c["r_e"] for c in cells})
    p_values = sorted({c["p"] for c in cells})
    shots = doc["parameters"]["shots"]
    with _styled():
        fig, ax = plt.subplots(figsize=(5.2, 3.9), layout="constrained")
        for k, r_e in enumerate(r_es):
            color = matplotlib.colormaps["viridis"](0.15 + 0.7 * k / max(len(r_es) - 1, 1))
            for p, style in zip(p_values, ("-", "--", ":", "-."), strict=False):
                row = sorted((c for c in cells if c["r_e"] == r_e and c["p"] == p), key=_by_d)
                ds = [c["d"] for c in row]
                ys = [c["ratio_lb"] if c["gated"] else c["ratio"] for c in row]
                ax.plot(ds, ys, linestyle=style, color=color, label=f"$r_e={r_e:g}$, p={p:.1%}")
                for c in row:
                    y = c["ratio_lb"] if c["gated"] else c["ratio"]
                    if c["gated"]:
                        ax.plot(c["d"], y, marker="^", markerfacecolor="none", color=color)
                    else:
                        ax.plot(c["d"], y, marker="o", color=color)
        ax.set_yscale("log")
        ax.set_xticks(sorted({c["d"] for c in cells}))
        ax.set_xlabel("code distance $d$")
        ax.set_ylabel(r"per-round $p_L(\mathrm{blind}) / p_L(\mathrm{herald})$")
        ax.set_title(f"Herald-aware vs blind on identical shots ({shots:,} per cell)")
        ax.text(
            0.02,
            0.98,
            "open triangles: < 50 herald errors; >= 95% one-sided lower bound",
            transform=ax.transAxes,
            fontsize=6.5,
            va="top",
            ha="left",
            color="0.4",
        )
        ax.legend(fontsize=6.5, loc="lower right")
        return _save(fig, out_path)


def _by_d(cell: dict[str, float]) -> float:
    return cell["d"]


def figure_hook_regression(out_path: Path, *, distances: Sequence[int] = (3, 5, 7)) -> Path:
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
    ablation_json: str | Path | None = None,
) -> list[Path]:
    """Regenerate every figure from the CSVs found under ``data_dir``.

    Reads ``<name>.csv`` (and any ``<name>_lambda.csv``) for each config name,
    pooling all points so the threshold panels show every available r_e. The
    ablation figure is drawn from ``ablation_json`` and skipped (no file
    written) when it is ``None`` or does not exist. The hook-regression figure
    needs no input. Returns the written figure paths.
    """
    data = Path(data_dir)
    figures = Path(figures_dir)
    csvs = [data / f"{n}.csv" for n in config_names]
    csvs += [data / f"{n}_lambda.csv" for n in config_names]
    points = load_sweep([c for c in csvs if c.exists()])

    written: list[Path] = []
    if points:
        csv_figures = [
            (
                "threshold_panels.png",
                figure_threshold_panels(points, figures / "threshold_panels.png", seed=seed),
            ),
            ("lambda_vs_p.png", figure_lambda(points, figures / "lambda_vs_p.png", seed=seed)),
        ]
        for name, result in csv_figures:
            if result is not None:
                written.append(result)
            else:
                print(f"skipped {name}: insufficient data")
    else:
        print(f"no sweep CSVs found under {data}; only the hook figure will render")
    ablation = None
    if ablation_json is not None:
        ablation = figure_ablation(Path(ablation_json), figures / "ablation.png")
    if ablation is not None:
        written.append(ablation)
    else:
        print(f"skipped ablation.png: no ablation results at {ablation_json}")
    written.append(figure_hook_regression(figures / "hook_regression.png"))
    return written


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=Path("data"))
    parser.add_argument("--figures-dir", type=Path, default=Path("figures"))
    parser.add_argument("--ablation-json", type=Path, default=Path("results/ablation.json"))
    parser.add_argument("--seed", type=int, default=BOOTSTRAP_SEED)
    args = parser.parse_args()
    paths = render_all(
        args.data_dir, args.figures_dir, seed=args.seed, ablation_json=args.ablation_json
    )
    for p in paths:
        print(f"wrote {p}")


if __name__ == "__main__":
    main()
