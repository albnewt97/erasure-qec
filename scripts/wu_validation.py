"""Wu et al. (2022) validation: crossings, collapse fits and the comparison with Wu.

Reads the seven data/wu/*.csv sweeps (convention wu2022; sinter merges rows that
share a strong_id). For each config and decoder (wu_herald_mwpm, wu_blind_mwpm):

- the (9,13) and (7,9) per-shot pairwise crossings
  (:func:`erasure_qec.analysis.crossing.bootstrap_pairwise_crossing`: first
  g < 0 -> g >= 0 sign change of g = ln P_L(d2) - ln P_L(d1), linear in
  (ln p, g); 95% percentile CI of an independent-binomial parametric bootstrap,
  n_boot = 1000, seed = 0; binomial noise only, interpolation error inside the
  bracket not included), with the bracket, n_sign_changes,
  n_sign_changes_unsaturated (CrossingEstimate.n_sign_changes_unsaturated: sign
  changes between adjacent usable p at which P_L_shot < SATURATION_CAP = 0.4 at
  both distances, the cap of :mod:`erasure_qec.analysis.threshold_fit`) and the
  number of bootstrap replicates without a crossing. Saturated-crossing rule
  (added after the P5b blind rows were inspected; every row and pair): a
  crossing whose bracket has P_L_shot >= SATURATION_CAP at either distance at
  either bracket point is reported as "saturated, not a threshold" with no
  value, and bootstrap replicates whose crossing is saturated count as
  no-crossing. Where the point
  estimate has no crossing the status is "no crossing in grid (below)"
  (g >= 0 at the first usable p) or "(above)" (g < 0 at every usable p). A
  wu_blind_mwpm pair with no unsaturated crossing anywhere on its grid is
  reported as BLIND_BELOW. No off-grid value or CI is given;
- fit_threshold(metric="shot", d_min=7, window_factor=1.5, n_boot=1000,
  seed=0) and the gap (collapse - (9,13) pairwise) in percentage points;
- for wu_herald_mwpm: the edge-rule check of docs/archive/remediation/WU_DPAIR_DECISION.md
  on the config's full grid (p_values plus extra_p_edge), Wu's reference value,
  ours - Wu in pp and whether Wu's value lies inside our (9,13) 95% CI.

Prints markdown tables (3 decimals; rates in %), writes the full-precision
results with git describe, data sha256 digests, software versions, n_boot and
seed to results/wu_validation.json, and draws figures/wu_validation.png
(byte-reproducible; plotting style of :mod:`erasure_qec.analysis.plotting`).
Refuses to run on a dirty work tree. Regenerate the table, JSON and figure with:

    uv run python scripts/wu_validation.py
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib
import matplotlib.ticker
import numpy as np
import pymatching
import scipy
import sinter
import stim

from erasure_qec.analysis.crossing import bootstrap_pairwise_crossing
from erasure_qec.analysis.plotting import _save, _styled
from erasure_qec.analysis.statistics import SweepPoint, metric_value, wilson_interval
from erasure_qec.analysis.threshold_fit import _SATURATION_CAP, fit_threshold

ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "data" / "wu"
CONFIG_DIR = ROOT / "experiments" / "configs" / "wu"
JSON_OUT = ROOT / "results" / "wu_validation.json"
FIG_OUT = ROOT / "figures" / "wu_validation.png"

HERALD, BLIND = "wu_herald_mwpm", "wu_blind_mwpm"
PAIRS = ((9, 13), (7, 9))
N_BOOT = 1000
SEED = 0
D_MIN = 7
WINDOW_FACTOR = 1.5
SATURATION_CAP = _SATURATION_CAP
WU_LINE_LABEL = "Wu et al. (UF, d = 11/15)"
SATURATED = "saturated, not a threshold"
BLIND_BELOW = "no unsaturated crossing in grid; blind threshold below the grid"


@dataclass(frozen=True)
class WuRef:
    """Wu et al.'s reported threshold for one config (``value`` None if not tabulated)."""

    value: float | None
    text: str
    source: str


CONFIGS: dict[str, WuRef] = {
    "wu_fig3_re000": WuRef(0.00937, "0.937(4)%", "Fig. 3a"),
    "wu_fig3b_re050": WuRef(None, "not tabulated", "not tabulated"),
    "wu_fig3b_re090": WuRef(None, "not tabulated", "not tabulated"),
    "wu_fig3_re098": WuRef(0.0415, "4.15(2)%", "Fig. 3a"),
    "wu_fig3b_re100": WuRef(0.0513, "5.13%", "Sec. III text"),
    "wu_si_pm005": WuRef(0.0380, "3.80(2)%", "SI S6"),
    "wu_si_pmp": WuRef(0.0285, "2.85(1)%", "SI S6"),
}
FIGURE_B = ("wu_fig3_re000", "wu_fig3b_re050", "wu_fig3b_re090", "wu_fig3_re098", "wu_fig3b_re100")


def _load_collect_wu() -> Any:
    path = ROOT / "experiments" / "collect_wu.py"
    spec = importlib.util.spec_from_file_location("collect_wu", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules["collect_wu"] = module
    spec.loader.exec_module(module)
    return module


def _finite(x: float) -> float | None:
    """JSON-safe float: NaN/inf become null."""
    return x if math.isfinite(x) else None


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _git_describe() -> str:
    out = subprocess.run(
        ["git", "describe", "--always", "--dirty"],
        capture_output=True,
        text=True,
        check=False,
        cwd=ROOT,
    )
    return out.stdout.strip() or "unknown"


def load_points(path: Path, config: str, decoder: str) -> list[SweepPoint]:
    """``SweepPoint`` per (d, p) of ``decoder`` in ``path``, with ``r_e`` from ``R_e``."""
    points: list[SweepPoint] = []
    for s in sinter.read_stats_from_csv_files(str(path)):
        meta = s.json_metadata
        if (
            s.decoder != decoder
            or s.shots <= 0
            or not isinstance(meta, dict)
            or meta.get("convention") != "wu2022"
            or meta.get("config") != config
        ):
            continue
        points.append(
            SweepPoint(
                decoder, int(meta["d"]), int(meta["rounds"]), float(meta["p"]),
                float(meta["R_e"]), int(s.shots), int(s.errors),
            )
        )  # fmt: skip
    return points


def usable_g(points: list[SweepPoint], d1: int, d2: int) -> list[tuple[float, float, bool]]:
    """``(p, g, unsaturated)`` at every p usable by the crossing estimator, ascending.

    Usable: present at both distances with a non-zero per-shot rate at both.
    ``unsaturated``: the per-shot rate is below ``SATURATION_CAP`` at both.
    """
    c1 = {pt.p: metric_value(pt, "shot") for pt in points if pt.d == d1}
    c2 = {pt.p: metric_value(pt, "shot") for pt in points if pt.d == d2}
    out = []
    for p in sorted(set(c1) & set(c2)):
        if c1[p] > 0.0 and c2[p] > 0.0:
            unsat = c1[p] < SATURATION_CAP and c2[p] < SATURATION_CAP
            out.append((p, math.log(c2[p]) - math.log(c1[p]), unsat))
    return out


def no_crossing_status(gs: list[float]) -> str:
    """Status when the curves have no g < 0 -> g >= 0 change on the grid."""
    if not gs:
        return "no usable points"
    # With no g < 0 -> g >= 0 change, g < 0 at the first point implies g < 0 at
    # every point (above). g >= 0 at the first point is classified "below", as
    # in the edge rule, even if g later drops below 0 and stays there.
    return "no crossing in grid (below)" if gs[0] >= 0.0 else "no crossing in grid (above)"


def unsaturated_crossings(g: list[tuple[float, float, bool]]) -> int:
    """Adjacent usable pairs with g < 0 -> g >= 0 whose two points are both unsaturated."""
    return sum(a[1] < 0.0 <= b[1] and a[2] and b[2] for a, b in zip(g, g[1:], strict=False))


def edge_touched(bracket: tuple[float, float] | None, gs: list[float], grid: list[float]) -> str:
    """The edge-rule outcome of WU_DPAIR_DECISION.md: "none", "lower" or "upper"."""
    if bracket is not None:
        if bracket[0] == min(grid):
            return "lower"
        return "upper" if bracket[1] == max(grid) else "none"
    if gs and all(g < 0.0 for g in gs):
        return "upper"
    if gs and gs[0] >= 0.0:
        return "lower"
    return "none"


def wu_inside(ci: tuple[float | None, float | None], wu: float | None) -> bool | None:
    """Whether ``wu`` lies in the closed interval ``ci``; None if either is missing."""
    lo, hi = ci
    if wu is None or lo is None or hi is None:
        return None
    return lo <= wu <= hi


def crossing_entry(points: list[SweepPoint], d1: int, d2: int, decoder: str) -> dict[str, Any]:
    est, ci, n_no = bootstrap_pairwise_crossing(
        points, d1, d2, metric="shot", n_boot=N_BOOT, seed=SEED, saturation_cap=SATURATION_CAP
    )
    g = usable_g(points, d1, d2)
    gs = [x[1] for x in g]
    crossed = est.bracket is not None and not est.saturated
    status = "crossing" if crossed else no_crossing_status(gs)
    if est.saturated:
        status = SATURATED
    n_unsat = unsaturated_crossings(g)
    if decoder == BLIND and n_unsat == 0:
        if status == "no crossing in grid (above)":
            raise SystemExit(f"{decoder} ({d1},{d2}): g < 0 at every usable p; not 'below'")
        status = BLIND_BELOW
    return {
        "d1": d1,
        "d2": d2,
        "status": status,
        "value": _finite(est.value) if crossed else None,
        "ci95": [_finite(ci[0]), _finite(ci[1])] if crossed else None,
        "bracket": list(est.bracket) if est.bracket is not None else None,
        "bracket_max_p_l_shot": _finite(est.bracket_max_p_l_shot),
        "saturated": est.saturated,
        "n_unsaturated_crossings": n_unsat,
        "n_sign_changes": est.n_sign_changes,
        "n_sign_changes_unsaturated": est.n_sign_changes_unsaturated,
        "n_no_crossing": n_no,
        "g": [[p, gv] for p, gv, _ in g],
    }


def fit_entry(points: list[SweepPoint]) -> dict[str, Any]:
    f = fit_threshold(
        points, d_min=D_MIN, n_boot=N_BOOT, seed=SEED, window_factor=WINDOW_FACTOR, metric="shot"
    )
    return {
        "p_th": _finite(f.p_th),
        "p_th_ci95": [_finite(f.p_th_ci[0]), _finite(f.p_th_ci[1])],
        "nu": _finite(f.nu),
        "chi2_dof": _finite(f.chi2_dof),
        "converged": f.converged,
        "resolved": f.resolved,
        "reason": f.reason,
        "message": f.message,
        "rel_ci_width": _finite(f.rel_ci_width),
        "n_points": f.n_points,
        "distances": list(f.distances),
        "p_center": _finite(f.p_center),
        "n_boot_failed": f.n_boot_failed,
    }


def analyse(collect_wu: Any) -> dict[str, Any]:
    results: dict[str, Any] = {}
    for name, ref in CONFIGS.items():
        path = DATA_DIR / f"{name}.csv"
        config = collect_wu.load_wu_config(CONFIG_DIR / f"{name}.yaml")
        grid = sorted(config.p_values)
        rows = []
        for decoder in (HERALD, BLIND):
            if decoder not in config.decoders:
                continue
            points = load_points(path, name, decoder)
            if sorted({pt.p for pt in points}) != grid:
                raise SystemExit(f"{name} {decoder}: data p values differ from the config grid")
            pairs = {f"{d1},{d2}": crossing_entry(points, d1, d2, decoder) for d1, d2 in PAIRS}
            fit = fit_entry(points)
            x913 = pairs["9,13"]["value"]
            row: dict[str, Any] = {
                "decoder": decoder,
                "n_points": len(points),
                "crossings": pairs,
                "fit_shot_dmin7": fit,
                "gap_collapse_minus_913_pp": (
                    100.0 * (fit["p_th"] - x913)
                    if fit["p_th"] is not None and x913 is not None
                    else None
                ),
            }
            if decoder == HERALD:
                c = pairs["9,13"]
                br = tuple(c["bracket"]) if c["bracket"] is not None else None
                row["edge_rule"] = edge_touched(br, [g for _, g in c["g"]], grid)
                row["wu"] = {"value": ref.value, "text": ref.text, "source": ref.source}
                both = x913 is not None and ref.value is not None
                row["ours_minus_wu_pp"] = 100.0 * (x913 - ref.value) if both else None
                ci = c["ci95"] or [None, None]
                row["wu_inside_ci"] = wu_inside((ci[0], ci[1]), ref.value)
                x79 = pairs["7,9"]["value"]
                row["shift_913_minus_79_pp"] = (
                    100.0 * (x913 - x79) if x913 is not None and x79 is not None else None
                )
            rows.append(row)
        results[name] = {
            "R_e": config.R_e,
            "p_m": "p" if config.p_m_follows_p else config.p_m,
            "grid": grid,
            "csv_sha256": _sha256(path),
            "rows": rows,
        }
    return results


def _pct(x: float | None) -> str:
    return "—" if x is None else f"{100.0 * x:.3f}"


def _num(x: float | None) -> str:
    """3 decimals; a value that rounds to zero prints as 0.000 (the JSON keeps the sign)."""
    if x is None:
        return "—"
    text = f"{x:.3f}"
    return "0.000" if text == "-0.000" else text


def _ci(ci: list[float | None] | None) -> str:
    return "—" if not ci or None in ci else f"[{_pct(ci[0])}, {_pct(ci[1])}]"


def _crossing_like(value: float | None, ci: list[float | None] | None) -> str:
    return "—" if value is None else f"{_pct(value)} {_ci(ci)}"


def _crossing_cell(c: dict[str, Any]) -> str:
    return c["status"] if c["value"] is None else _crossing_like(c["value"], c["ci95"])


def markdown(results: dict[str, Any]) -> str:
    lines = [
        "Per-shot pairwise crossings, % (95% bootstrap percentile CI, "
        f"n_boot = {N_BOOT}, seed = {SEED}); sc = n_sign_changes, "
        "sc_u = n_sign_changes_unsaturated, nx = replicates without a crossing.",
        "",
        "| config | R_e | p_m | decoder | (9,13) | bracket | sc / sc_u / nx "
        "| (7,9) | bracket | sc / sc_u / nx |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]
    for name, entry in results.items():
        for row in entry["rows"]:
            cells = []
            for key in ("9,13", "7,9"):
                c = row["crossings"][key]
                br = c["bracket"]
                cells += [
                    _crossing_cell(c),
                    "—" if br is None else f"[{_pct(br[0])}, {_pct(br[1])}]",
                    f"{c['n_sign_changes']} / {c['n_sign_changes_unsaturated']} "
                    f"/ {c['n_no_crossing']}",
                ]
            lines.append(
                f"| {name} | {entry['R_e']} | {entry['p_m']} | {row['decoder']} | "
                + " | ".join(cells)
                + " |"
            )
    lines += [
        "",
        f"fit_threshold(metric = shot, d_min = {D_MIN}), % (95% bootstrap percentile CI, "
        f"n_boot = {N_BOOT}); gap = collapse − (9,13) pairwise, pp.",
        "",
        "| config | decoder | p_th | ν | χ²/dof | converged | resolved | reason | gap (pp) |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for name, entry in results.items():
        for row in entry["rows"]:
            f = row["fit_shot_dmin7"]
            lines.append(
                f"| {name} | {row['decoder']} | {_crossing_like(f['p_th'], f['p_th_ci95'])} "
                f"| {_num(f['nu'])} | {_num(f['chi2_dof'])} | {f['converged']} "
                f"| {f['resolved']} | {f['reason']} | {_num(row['gap_collapse_minus_913_pp'])} |"
            )
    lines += [
        "",
        "Comparison with Wu et al. (wu_herald_mwpm, (9,13) per-shot crossing).",
        "",
        "| config | ours, % [95% CI] | Wu et al. | source | ours − Wu (pp) | Wu inside CI "
        "| edge rule |",
        "|---|---|---|---|---|---|---|",
    ]
    for name, entry in results.items():
        row = entry["rows"][0]
        c = row["crossings"]["9,13"]
        inside = row["wu_inside_ci"]
        lines.append(
            f"| {name} | {_crossing_cell(c)} | {row['wu']['text']} | {row['wu']['source']} "
            f"| {_num(row['ours_minus_wu_pp'])} | {'—' if inside is None else inside} "
            f"| {row['edge_rule']} |"
        )
    lines += [
        "",
        "Herald saturation diagnostics: sc_u = n_sign_changes_unsaturated; max P_L = largest "
        "per-shot P_L at the two distances and two bracket points.",
        "",
        "| config | (9,13) sc_u | (9,13) max P_L | (7,9) sc_u | (7,9) max P_L |",
        "|---|---|---|---|---|",
    ]
    for name, entry in results.items():
        cs = entry["rows"][0]["crossings"]
        lines.append(
            f"| {name} | "
            + " | ".join(
                f"{cs[k]['n_sign_changes_unsaturated']} | {_num(cs[k]['bracket_max_p_l_shot'])}"
                for k in ("9,13", "7,9")
            )
            + " |"
        )
    lines += [
        "",
        "Herald finite-size window: shift = (9,13) − (7,9) per-shot crossing, pp.",
        "",
        "| config | (7,9), % [95% CI] | (9,13), % [95% CI] | shift (pp) |",
        "|---|---|---|---|",
    ]
    for name, entry in results.items():
        row = entry["rows"][0]
        cs = row["crossings"]
        lines.append(
            f"| {name} | {_crossing_cell(cs['7,9'])} | {_crossing_cell(cs['9,13'])} "
            f"| {_num(row['shift_913_minus_79_pp'])} |"
        )
    return "\n".join(lines)


def _curve(points: list[SweepPoint], d: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    pts = sorted((pt for pt in points if pt.d == d), key=lambda pt: pt.p)
    p = np.array([pt.p for pt in pts])
    y = np.array([pt.p_l_shot for pt in pts])
    ci = np.array([wilson_interval(pt.errors, pt.shots) for pt in pts])
    return p, y, np.vstack([y - ci[:, 0], ci[:, 1] - y])


def figure(results: dict[str, Any], out: Path) -> Path:
    import matplotlib.pyplot as plt

    with _styled():
        fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(11.0, 4.4))
        colours = {"wu_fig3_re000": "tab:blue", "wu_fig3_re098": "tab:orange"}
        styles = {9: "--", 13: "-"}
        for name, colour in colours.items():
            points = load_points(DATA_DIR / f"{name}.csv", name, HERALD)
            r_e = results[name]["R_e"]
            for d, ls in styles.items():
                p, y, err = _curve(points, d)
                ax_a.errorbar(
                    100 * p, y, yerr=err, color=colour, linestyle=ls, marker="o",
                    markersize=3, capsize=0, label=f"R_e = {r_e}, d = {d}",
                )  # fmt: skip
            c = results[name]["rows"][0]["crossings"]["9,13"]
            if c["value"] is not None:
                p9, y9, _ = _curve(points, 9)
                yx = float(np.exp(np.interp(np.log(c["value"]), np.log(p9), np.log(y9))))
                lo, hi = c["ci95"]
                ax_a.errorbar(
                    [100 * c["value"]], [yx],
                    xerr=[[100 * (c["value"] - lo)], [100 * (hi - c["value"])]],
                    color=colour, marker="D", markersize=6, markeredgecolor="black",
                    capsize=3, linestyle="none", zorder=5,
                    label=f"(9,13) crossing, R_e = {r_e}",
                )  # fmt: skip
        for i, wu in enumerate((0.00937, 0.0415)):
            ax_a.axvline(
                100 * wu, color="0.3", linestyle="--", linewidth=1.0,
                label=WU_LINE_LABEL if i == 0 else None,
            )  # fmt: skip
        ax_a.set_xscale("log")
        ax_a.set_yscale("log")
        ticks = [0.8, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0]
        ax_a.set_xticks(ticks, [f"{t:g}" for t in ticks])
        ax_a.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax_a.set_xlabel("physical error rate p (%)")
        ax_a.set_ylabel("logical failure probability after T = d rounds")
        ax_a.set_title("(a) wu_herald_mwpm, d = 9 and 13", loc="left")
        ax_a.legend(fontsize=7, loc="lower center")

        xs, ys, los, his, wx, wy = [], [], [], [], [], []
        for name in FIGURE_B:
            entry = results[name]
            row = entry["rows"][0]
            c = row["crossings"]["9,13"]
            if c["value"] is not None:
                xs.append(entry["R_e"])
                ys.append(100 * c["value"])
                los.append(100 * (c["value"] - c["ci95"][0]))
                his.append(100 * (c["ci95"][1] - c["value"]))
            if row["wu"]["value"] is not None:
                wx.append(entry["R_e"])
                wy.append(100 * row["wu"]["value"])
        ax_b.errorbar(
            xs, ys, yerr=[los, his], color="tab:green", marker="o", linestyle="none",
            capsize=3, label="this work: wu_herald_mwpm, (9,13) per-shot crossing, 95% CI",
        )  # fmt: skip
        ax_b.plot(
            wx, wy, marker="s", markerfacecolor="none", color="black", linestyle="none",
            markersize=7, label="Wu et al. (tabulated values only)",
        )  # fmt: skip
        ax_b.set_xlabel("erasure fraction R_e")
        ax_b.set_ylabel("threshold p (%)\n(logical failure probability after T = d rounds)")
        ax_b.set_title("(b) threshold vs R_e (p_m = 0)", loc="left")
        ax_b.legend(fontsize=7, loc="upper left")
        fig.tight_layout()
        return _save(fig, out)


def main() -> None:
    describe = _git_describe()
    if describe.endswith("-dirty") or describe == "unknown":
        raise SystemExit(f"git describe is {describe!r}; commit or stash changes first")
    results = analyse(_load_collect_wu())
    table = markdown(results)
    print(table)
    payload = {
        "description": "Wu et al. (2022) validation: per-shot (9,13) and (7,9) crossings, "
        "d_min=7 shot collapse fits, comparison with Wu's reported thresholds",
        "generated_by": "scripts/wu_validation.py",
        "git_describe": describe,
        "data_sha256": {f"data/wu/{p.name}": _sha256(p) for p in sorted(DATA_DIR.glob("*.csv"))},
        "software": {
            "python": sys.version.split()[0],
            "stim": stim.__version__,
            "pymatching": pymatching.__version__,
            "sinter": sinter.__version__,
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "matplotlib": matplotlib.__version__,
        },
        "parameters": {
            "n_boot": N_BOOT,
            "seed": SEED,
            "d_min": D_MIN,
            "window_factor": WINDOW_FACTOR,
            "metric": "shot",
            "pairs": [list(p) for p in PAIRS],
            "ci_percentiles": [2.5, 97.5],
            "saturation_cap": SATURATION_CAP,
            "crossing_rule": "first g<0 -> g>=0, g = ln P_L_shot(d2) - ln P_L_shot(d1), "
            "linear in (ln p, g); points with zero rate at either distance skipped",
            "n_sign_changes_unsaturated": "sign changes between adjacent usable p with "
            "P_L_shot < saturation_cap at both distances",
            "wu_inside_ci": "Wu's central value in the closed (9,13) herald 95% CI; "
            "Wu's own uncertainty not used",
        },
        "configs": results,
        "markdown": table,
    }
    JSON_OUT.parent.mkdir(parents=True, exist_ok=True)
    JSON_OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print(f"\nwrote {JSON_OUT.relative_to(ROOT)}")
    print(f"wrote {figure(results, FIG_OUT).relative_to(ROOT)}")


if __name__ == "__main__":
    main()
