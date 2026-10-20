"""Headline per-shot thresholds and the herald-minus-blind difference, dense sweeps.

For each committed sweep (data/baseline_pauli.csv, data/erasure_r50.csv,
data/erasure_r98.csv; r_e = 0, 0.5, 0.98) and each decoder in it, reports:

- headline: the per-shot (7,11) pairwise crossing
  (:func:`erasure_qec.analysis.crossing.headline_crossing`: first g < 0 -> g >= 0
  sign change of g = ln P_L_shot(11) - ln P_L_shot(7), linear in (ln p, g); 95%
  percentile CI of an independent-binomial parametric bootstrap, n_boot = 1000,
  seed = 0; binomial noise only, interpolation error inside the bracket not
  included), with the bracket, n_sign_changes, n_sign_changes_unsaturated, the
  number of replicates without a crossing and the saturated-crossing rule of
  scripts/wu_validation.py: a crossing whose bracket has P_L_shot >=
  SATURATION_CAP = 0.4 at either distance at either bracket point is reported as
  "saturated, not a threshold" with no value, and saturated bootstrap replicates
  count as no-crossing;
- drift context: the (5,7), (7,9) and (9,11) per-shot crossings, same estimator
  and fields;
- systematic: fit_threshold(metric="shot", d_min=7, window_factor=1.5,
  n_boot=1000, seed=0) and the gap (collapse - (7,11) pairwise) in percentage
  points;
- secondary: fit_threshold(metric="round", d_min=7), labelled
  "per-round (secondary)".

For each r_e it also reports the herald-minus-blind difference of the (7,11)
per-shot crossings (:func:`erasure_qec.analysis.crossing.bootstrap_crossing_difference`:
each replicate resamples the herald points, then the blind points, independently
from one numpy.random.default_rng(seed) stream, n_boot = 1000, seed = 0), its 95%
percentile CI and the fraction of replicates in which either crossing is missing
or saturated. r_e = 0 has no heralds and is the control.

Every headline, drift and fit value is asserted equal to the corresponding entry
of results/threshold_table_dense.json (same estimator, data, n_boot and seed;
that file applies no saturation cap, so a crossing or replicate the cap changes
fails the check). Refuses to run on a dirty work tree. Prints markdown tables (3
decimals; rates in %) and writes results/headline.json with git describe, data
sha256 digests, software versions, n_boot and seed:

    uv run python scripts/headline_results.py
"""

from __future__ import annotations

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

from erasure_qec.analysis.crossing import (
    HEADLINE_PAIR,
    CrossingEstimate,
    bootstrap_crossing_difference,
    bootstrap_pairwise_crossing,
    headline_crossing,
)
from erasure_qec.analysis.statistics import SweepPoint, available_conditions, load_sweep
from erasure_qec.analysis.threshold_fit import SATURATION_CAP, FitResult, fit_threshold

ROOT = Path(__file__).resolve().parent.parent
DATA_FILES = ("data/baseline_pauli.csv", "data/erasure_r50.csv", "data/erasure_r98.csv")
DENSE_JSON = ROOT / "results" / "threshold_table_dense.json"
JSON_OUT = ROOT / "results" / "headline.json"

HERALD, BLIND = "herald_mwpm", "blind_mwpm"
DRIFT_PAIRS = ((5, 7), (7, 9), (9, 11))
N_BOOT = 1000
SEED = 0
D_MIN = 7
WINDOW_FACTOR = 1.5
SATURATED = "saturated, not a threshold"
SECONDARY_LABEL = "per-round (secondary)"
# Keys of a pairwise entry and of a fit entry that results/threshold_table_dense.json
# also records; each must agree exactly.
CROSS_CHECK_PAIR_KEYS = (
    "value",
    "ci95",
    "bracket",
    "n_sign_changes",
    "n_sign_changes_unsaturated",
    "n_no_crossing",
)


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


def crossing_status(est: CrossingEstimate) -> str:
    """ "crossing", "saturated, not a threshold" or "no crossing in grid"."""
    if est.saturated:
        return SATURATED
    return "crossing" if est.bracket is not None else "no crossing in grid"


def crossing_entry(boot: tuple[CrossingEstimate, tuple[float, float], int]) -> dict[str, Any]:
    """JSON entry of a ``bootstrap_pairwise_crossing`` result."""
    est, ci, n_no = boot
    status = crossing_status(est)
    crossed = status == "crossing"
    return {
        "d1": est.d1,
        "d2": est.d2,
        "status": status,
        "value": _finite(est.value) if crossed else None,
        "ci95": [_finite(ci[0]), _finite(ci[1])] if crossed else None,
        "bracket": list(est.bracket) if est.bracket is not None else None,
        "bracket_max_p_l_shot": _finite(est.bracket_max_p_l_shot),
        "saturated": est.saturated,
        "n_sign_changes": est.n_sign_changes,
        "n_sign_changes_unsaturated": est.n_sign_changes_unsaturated,
        "n_no_crossing": n_no,
    }


def fit_entry(f: FitResult) -> dict[str, Any]:
    return {
        "metric": f.metric,
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


def _fit(points: list[SweepPoint], metric: str) -> dict[str, Any]:
    return fit_entry(
        fit_threshold(
            points,
            d_min=D_MIN,
            n_boot=N_BOOT,
            seed=SEED,
            window_factor=WINDOW_FACTOR,
            metric="shot" if metric == "shot" else "round",
        )
    )


def condition_entry(
    path: str, decoder: str, r_e: float, points: list[SweepPoint]
) -> dict[str, Any]:
    headline = crossing_entry(headline_crossing(points, n_boot=N_BOOT, seed=SEED))
    drift = [
        crossing_entry(
            bootstrap_pairwise_crossing(
                points,
                d1,
                d2,
                metric="shot",
                n_boot=N_BOOT,
                seed=SEED,
                saturation_cap=SATURATION_CAP,
            )
        )
        for d1, d2 in DRIFT_PAIRS
    ]
    collapse = _fit(points, "shot")
    secondary = _fit(points, "round")
    gap = (
        100.0 * (collapse["p_th"] - headline["value"])
        if collapse["p_th"] is not None and headline["value"] is not None
        else None
    )
    return {
        "data": path,
        "decoder": decoder,
        "r_e": r_e,
        "headline": headline,
        "drift": drift,
        "collapse_shot_dmin7": collapse,
        "gap_collapse_minus_pairwise_pp": gap,
        "per_round_secondary": {"label": SECONDARY_LABEL, **secondary},
    }


def difference_entry(
    path: str, r_e: float, herald: list[SweepPoint], blind: list[SweepPoint]
) -> dict[str, Any]:
    value, ci, n_missing = bootstrap_crossing_difference(
        herald,
        blind,
        *HEADLINE_PAIR,
        metric="shot",
        n_boot=N_BOOT,
        seed=SEED,
        saturation_cap=SATURATION_CAP,
    )
    return {
        "data": path,
        "r_e": r_e,
        "pair": list(HEADLINE_PAIR),
        "value": _finite(value),
        "ci95": [_finite(ci[0]), _finite(ci[1])],
        "n_boot": N_BOOT,
        "seed": SEED,
        "n_missing_or_saturated": n_missing,
        "fraction_missing_or_saturated": n_missing / N_BOOT,
        "control": r_e == 0.0,
    }


def cross_check(conditions: list[dict[str, Any]], dense: dict[str, Any]) -> None:
    """Assert every headline, drift and fit value equals its threshold_table_dense.json entry."""
    for cond in conditions:
        matches = [
            c
            for c in dense["data"][cond["data"]]["conditions"]
            if c["decoder"] == cond["decoder"] and c["r_e"] == cond["r_e"]
        ]
        assert len(matches) == 1, (cond["data"], cond["decoder"], cond["r_e"])
        ref = matches[0]
        for entry in [cond["headline"], *cond["drift"]]:
            (row,) = [
                r
                for r in ref["pairwise"]
                if r["metric"] == "shot" and (r["d1"], r["d2"]) == (entry["d1"], entry["d2"])
            ]
            for key in CROSS_CHECK_PAIR_KEYS:
                assert entry[key] == row[key], (cond["decoder"], cond["r_e"], entry["d1"], key)
        for ours in (cond["collapse_shot_dmin7"], cond["per_round_secondary"]):
            (fit,) = [f for f in ref["fits"] if f["metric"] == ours["metric"]]
            for key, val in fit.items():
                assert ours[key] == val, (cond["decoder"], cond["r_e"], ours["metric"], key)


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


def _crossing_cell(c: dict[str, Any]) -> str:
    if c["status"] != "crossing":
        return str(c["status"])
    return f"{_pct(c['value'])} {_ci(c['ci95'])}"


def _bracket(c: dict[str, Any]) -> str:
    br = c["bracket"]
    return "—" if br is None else f"[{_pct(br[0])}, {_pct(br[1])}]"


def _fit_cell(f: dict[str, Any]) -> str:
    text = f"{_pct(f['p_th'])} {_ci(f['p_th_ci95'])}"
    return text if f["resolved"] else f"{text} (not resolved: {f['reason']})"


def markdown(conditions: list[dict[str, Any]], differences: list[dict[str, Any]]) -> str:
    d1, d2 = HEADLINE_PAIR
    lines = [
        f"Headline: per-shot ({d1},{d2}) pairwise crossing, % [95% CI] (bootstrap, "
        f"n_boot = {N_BOOT}, seed = {SEED}); sc / sc_u = n_sign_changes / "
        "n_sign_changes_unsaturated; nx = replicates without a crossing; gap = collapse − "
        "pairwise (pp).",
        "",
        f"| r_e | decoder | ({d1},{d2}), % [95% CI] | bracket, % | sc / sc_u | nx "
        f"| collapse d ≥ {D_MIN}, % [95% CI] | gap (pp) | {SECONDARY_LABEL}, % [95% CI] |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for c in conditions:
        h = c["headline"]
        lines.append(
            f"| {c['r_e']:g} | {c['decoder']} | {_crossing_cell(h)} | {_bracket(h)} "
            f"| {h['n_sign_changes']} / {h['n_sign_changes_unsaturated']} | {h['n_no_crossing']} "
            f"| {_fit_cell(c['collapse_shot_dmin7'])} "
            f"| {_num(c['gap_collapse_minus_pairwise_pp'])} "
            f"| {_fit_cell(c['per_round_secondary'])} |"
        )
    pairs = " | ".join(f"({a},{b}), % [95% CI] (sc_u)" for a, b in DRIFT_PAIRS)
    lines += [
        "",
        "Drift context: per-shot pairwise crossings of adjacent distances.",
        "",
        f"| r_e | decoder | {pairs} |",
        "|---|---|" + "---|" * len(DRIFT_PAIRS),
    ]
    for c in conditions:
        cells = " | ".join(
            f"{_crossing_cell(x)} ({x['n_sign_changes_unsaturated']})" for x in c["drift"]
        )
        lines.append(f"| {c['r_e']:g} | {c['decoder']} | {cells} |")
    lines += [
        "",
        f"Herald − blind, per-shot ({d1},{d2}) crossing, pp [95% CI] (independent resampling "
        f"of each decoder, n_boot = {N_BOOT}, seed = {SEED}); missing = replicates in which "
        "either crossing is missing or saturated.",
        "",
        "| r_e | herald − blind (pp) | 95% CI (pp) | missing | note |",
        "|---|---|---|---|---|",
    ]
    for diff in differences:
        lo, hi = diff["ci95"]
        ci = "—" if lo is None or hi is None else f"[{_num(100 * lo)}, {_num(100 * hi)}]"
        value = None if diff["value"] is None else 100.0 * diff["value"]
        note = "control (no heralds)" if diff["control"] else ""
        lines.append(
            f"| {diff['r_e']:g} | {_num(value)} | {ci} "
            f"| {diff['n_missing_or_saturated']}/{diff['n_boot']} | {note} |"
        )
    return "\n".join(lines)


def main() -> None:
    describe = _git_describe()
    if describe.endswith("-dirty") or describe == "unknown":
        raise SystemExit(f"git describe is {describe!r}; commit or stash changes first")
    conditions: list[dict[str, Any]] = []
    differences: list[dict[str, Any]] = []
    for path in DATA_FILES:
        points = load_sweep(ROOT / path)
        if not points:
            raise SystemExit(f"no points loaded from {path}")
        by_decoder: dict[str, list[SweepPoint]] = {}
        r_es = set()
        for decoder, r_e in available_conditions(points):
            cond_points = [pt for pt in points if pt.decoder == decoder and pt.r_e == r_e]
            conditions.append(condition_entry(path, decoder, r_e, cond_points))
            by_decoder[decoder] = cond_points
            r_es.add(r_e)
        if len(r_es) != 1 or set(by_decoder) != {HERALD, BLIND}:
            raise SystemExit(f"{path}: expected one r_e and both decoders")
        differences.append(
            difference_entry(path, r_es.pop(), by_decoder[HERALD], by_decoder[BLIND])
        )

    cross_check(conditions, json.loads(DENSE_JSON.read_text()))
    table = markdown(conditions, differences)
    print(table)
    payload = {
        "description": "headline per-shot (7,11) crossings with drift context, d_min=7 "
        "collapse fits (shot; round as secondary) and the herald-minus-blind difference",
        "generated_by": "scripts/headline_results.py",
        "git_describe": describe,
        "data_sha256": {path: _sha256(ROOT / path) for path in DATA_FILES},
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
            "headline_pair": list(HEADLINE_PAIR),
            "drift_pairs": [list(p) for p in DRIFT_PAIRS],
            "metric": "shot",
            "ci_percentiles": [2.5, 97.5],
            "saturation_cap": SATURATION_CAP,
            "crossing_rule": "first g<0 -> g>=0, g = ln P_L_shot(d2) - ln P_L_shot(d1), "
            "linear in (ln p, g); points with zero rate at either distance skipped",
            "saturated_crossing_rule": "a crossing whose bracket has P_L_shot >= "
            "saturation_cap at either distance at either bracket point has no value; "
            "saturated bootstrap replicates count as no-crossing",
            "difference_rule": "herald minus blind; each replicate resamples herald then "
            "blind points independently from one default_rng(seed) stream; a replicate "
            "with either crossing missing or saturated is counted, not used",
            "cross_check": "headline, drift and fit values equal "
            "results/threshold_table_dense.json",
        },
        "conditions": conditions,
        "herald_minus_blind": differences,
        "markdown": table,
    }
    JSON_OUT.parent.mkdir(parents=True, exist_ok=True)
    JSON_OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print(f"\nwrote {JSON_OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
