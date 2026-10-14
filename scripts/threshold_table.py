"""Pairwise crossings and d >= 7 collapse fits for the committed sweeps, both metrics.

For each committed sweep (data/baseline_pauli.csv, data/erasure_r50.csv,
data/erasure_r98.csv) and each decoder in it, reports:

- the pairwise crossing of the (5,7), (7,9), (9,11) and (7,11) curves in each
  metric ("shot" = P_L_shot, "round" = per-round p_L at T = d), with a 95%
  parametric-bootstrap percentile CI (binomial noise only; interpolation error
  inside the bracket is not included), the bracketing grid points, the number
  of sign changes of g = ln m(d2) - ln m(d1), and the number of bootstrap
  replicates without a crossing;
- fit_threshold(d_min=7) in each metric: p_th with its 95% bootstrap
  percentile CI, nu, chi2/dof, converged, resolved and the reason code.

Prints markdown tables to stdout and writes the full-precision results, with
software versions, data-file sha256 digests, n_boot and seed, to
results/threshold_table_committed.json. Run:

    uv run python scripts/threshold_table.py
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pymatching
import scipy
import sinter
import stim

from erasure_qec.analysis.crossing import bootstrap_pairwise_crossing
from erasure_qec.analysis.statistics import Metric, available_conditions, load_sweep
from erasure_qec.analysis.threshold_fit import fit_threshold

DATA_FILES = ("data/baseline_pauli.csv", "data/erasure_r50.csv", "data/erasure_r98.csv")
PAIRS = ((5, 7), (7, 9), (9, 11), (7, 11))
METRICS: tuple[Metric, ...] = ("shot", "round")
N_BOOT = 1000
SEED = 0
D_MIN = 7
WINDOW_FACTOR = 1.5  # fit_threshold's default, recorded explicitly


def _finite(x: float) -> float | None:
    """JSON-safe float: NaN/inf become null."""
    return x if math.isfinite(x) else None


def _sha256(path: str) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _git_describe() -> str:
    out = subprocess.run(
        ["git", "describe", "--always", "--dirty"], capture_output=True, text=True, check=False
    )
    return out.stdout.strip() or "unknown"


def _pct(x: float | None, digits: int = 3) -> str:
    return "—" if x is None else f"{100.0 * x:.{digits}f}"


def _ci(ci: tuple[float | None, float | None]) -> str:
    lo, hi = ci
    return "—" if lo is None or hi is None else f"[{_pct(lo, 1)}, {_pct(hi, 1)}]"


def _condition_results(points: list[Any]) -> dict[str, Any]:
    pairwise: list[dict[str, Any]] = []
    for metric in METRICS:
        for d1, d2 in PAIRS:
            est, ci, n_no = bootstrap_pairwise_crossing(
                points, d1, d2, metric=metric, n_boot=N_BOOT, seed=SEED
            )
            pairwise.append(
                {
                    "metric": metric,
                    "d1": d1,
                    "d2": d2,
                    "value": _finite(est.value),
                    "bracket": list(est.bracket) if est.bracket is not None else None,
                    "n_sign_changes": est.n_sign_changes,
                    "ci95": [_finite(ci[0]), _finite(ci[1])],
                    "n_no_crossing": n_no,
                }
            )
    fits: list[dict[str, Any]] = []
    for metric in METRICS:
        f = fit_threshold(
            points,
            d_min=D_MIN,
            n_boot=N_BOOT,
            seed=SEED,
            window_factor=WINDOW_FACTOR,
            metric=metric,
        )
        fits.append(
            {
                "metric": metric,
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
        )
    return {"pairwise": pairwise, "fits": fits}


def _print_tables(results: dict[str, Any]) -> None:
    print(f"## Pairwise crossings (bootstrap 95% percentile CI, n_boot={N_BOOT}, seed={SEED})\n")
    print(
        "| data | decoder | metric | pair | crossing (%) | 95% CI (%) | bracket (%) "
        "| n_sign_changes | n_no_crossing |"
    )
    print("|---|---|---|---|---|---|---|---|---|")
    for path, entry in results.items():
        for cond in entry["conditions"]:
            for row in cond["pairwise"]:
                br = row["bracket"]
                bracket = "—" if br is None else f"[{_pct(br[0])}, {_pct(br[1])}]"
                print(
                    f"| {Path(path).stem} | {cond['decoder']} | {row['metric']} "
                    f"| ({row['d1']},{row['d2']}) | {_pct(row['value'])} | {_ci(row['ci95'])} "
                    f"| {bracket} | {row['n_sign_changes']} | {row['n_no_crossing']} |"
                )
    print(f"\n## fit_threshold(d_min={D_MIN}) (bootstrap 95% percentile CI, n_boot={N_BOOT})\n")
    print(
        "| data | decoder | metric | p_th (%) | 95% CI (%) | nu | chi2/dof "
        "| converged | resolved | reason |"
    )
    print("|---|---|---|---|---|---|---|---|---|---|")
    for path, entry in results.items():
        for cond in entry["conditions"]:
            for f in cond["fits"]:
                nu = "—" if f["nu"] is None else f"{f['nu']:.2f}"
                chi = "—" if f["chi2_dof"] is None else f"{f['chi2_dof']:.2f}"
                print(
                    f"| {Path(path).stem} | {cond['decoder']} | {f['metric']} "
                    f"| {_pct(f['p_th'])} | {_ci(f['p_th_ci95'])} | {nu} | {chi} "
                    f"| {f['converged']} | {f['resolved']} | {f['reason']} |"
                )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", default="results/threshold_table_committed.json")
    args = parser.parse_args()

    results: dict[str, Any] = {}
    for path in DATA_FILES:
        points = load_sweep(path)
        if not points:
            raise SystemExit(f"no points loaded from {path}")
        conditions = []
        for decoder, r_e in available_conditions(points):
            cond_points = [pt for pt in points if pt.decoder == decoder and pt.r_e == r_e]
            conditions.append({"decoder": decoder, "r_e": r_e, **_condition_results(cond_points)})
        results[path] = {"sha256": _sha256(path), "conditions": conditions}

    _print_tables(results)

    payload = {
        "description": "pairwise crossings and d_min=7 collapse fits, metric shot and round",
        "generated_by": "scripts/threshold_table.py",
        "git_describe": _git_describe(),
        "software": {
            "python": sys.version.split()[0],
            "stim": stim.__version__,
            "pymatching": pymatching.__version__,
            "sinter": sinter.__version__,
            "numpy": np.__version__,
            "scipy": scipy.__version__,
        },
        "parameters": {
            "n_boot": N_BOOT,
            "seed": SEED,
            "d_min": D_MIN,
            "window_factor": WINDOW_FACTOR,
            "pairs": [list(p) for p in PAIRS],
            "metrics": list(METRICS),
            "ci_percentiles": [2.5, 97.5],
            "crossing_rule": "first g<0 -> g>=0, g = ln m(d2) - ln m(d1), linear in (ln p, g); "
            "points with zero metric value at either distance skipped",
        },
        "data": results,
    }
    out = Path(args.json)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
