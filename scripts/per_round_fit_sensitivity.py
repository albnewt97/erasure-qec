"""Sensitivity of the per-round d >= 7 collapse fit to its saturation cap and to pin retention.

Both checks run on the committed dense sweeps (data/baseline_pauli.csv,
data/erasure_r50.csv, data/erasure_r98.csv) with the estimator of
scripts/headline_results.py: fit_threshold(metric="round", d_min=7,
window_factor=1.5, n_boot=1000). Neither changes the published estimator.

- cap: for every (decoder, r_e) condition and every cap in CAPS, the fit with
  threshold_fit._SATURATION_CAP temporarily set to that cap (seed 0). The cap
  0.4 rows are asserted equal to the per-round entries of results/headline.json.
- pins: for r_e = 0.98 blind, seeds 0-29. Each bootstrap replicate is
  fit_threshold(n_boot=0) on that replicate's resampled points, drawn exactly as
  _bootstrap_pipeline draws them, so it re-runs estimate_crossing, the window
  and the weighted fit and applies the point fit's pin check. A replicate is
  "failed" if the fit returns no p_th, and "pinned" if it returns p_th but is
  not converged (p_th at a window edge or nu at a bound). The pins-retained CI
  uses every non-failed replicate, which is what fit_threshold does; it is
  asserted equal to fit_threshold's CI at every seed. The pins-excluded CI drops
  the pinned replicates, which is what bootstrap_threshold_difference does.
  "resolved" applies the 1.0 relative-width guard to each CI.

Refuses to run on a dirty work tree. Writes results/per_round_fit_sensitivity.json
with git describe, data sha256 digests, software versions and parameters:

    uv run python scripts/per_round_fit_sensitivity.py
"""

from __future__ import annotations

import hashlib
import json
import math
import subprocess
import sys
from dataclasses import replace
from pathlib import Path
from typing import Any
from unittest import mock

import numpy as np
import scipy

from erasure_qec.analysis import threshold_fit
from erasure_qec.analysis.statistics import SweepPoint, available_conditions, load_sweep
from erasure_qec.analysis.threshold_fit import FitResult, fit_threshold

ROOT = Path(__file__).resolve().parent.parent
DATA_FILES = ("data/baseline_pauli.csv", "data/erasure_r50.csv", "data/erasure_r98.csv")
HEADLINE_JSON = ROOT / "results" / "headline.json"
JSON_OUT = ROOT / "results" / "per_round_fit_sensitivity.json"

CAPS = (0.2, 0.25, 0.3, 0.4)
PUBLISHED_CAP = 0.4
PIN_CONDITION = ("blind_mwpm", 0.98)
PIN_SEEDS = tuple(range(30))
N_BOOT = 1000
SEED = 0
D_MIN = 7
WINDOW_FACTOR = 1.5
CI_PERCENTILES = threshold_fit._CI_PERCENTILES
MAX_REL_CI_WIDTH = threshold_fit._MAX_REL_CI_WIDTH


def _finite(x: float) -> float | None:
    return x if math.isfinite(x) else None


def _git_describe() -> str:
    out = subprocess.run(
        ["git", "describe", "--always", "--dirty"],
        capture_output=True,
        text=True,
        check=False,
        cwd=ROOT,
    )
    return out.stdout.strip() or "unknown"


def _fit(points: list[SweepPoint], *, n_boot: int, seed: int) -> FitResult:
    return fit_threshold(
        points,
        d_min=D_MIN,
        n_boot=n_boot,
        seed=seed,
        window_factor=WINDOW_FACTOR,
        metric="round",
    )


def cap_entry(cap: float, f: FitResult) -> dict[str, Any]:
    """JSON entry of one fit at one cap."""
    return {
        "cap": cap,
        "p_th": _finite(f.p_th),
        "p_th_ci95": [_finite(f.p_th_ci[0]), _finite(f.p_th_ci[1])],
        "nu": _finite(f.nu),
        "rel_ci_width": _finite(f.rel_ci_width),
        "converged": f.converged,
        "resolved": f.resolved,
        "reason": f.reason,
        "n_points": f.n_points,
    }


def cap_sensitivity(points: list[SweepPoint]) -> list[dict[str, Any]]:
    """The per-round fit of one condition at every cap in CAPS (seed 0)."""
    rows = []
    for cap in CAPS:
        with mock.patch.object(threshold_fit, "_SATURATION_CAP", cap):
            rows.append(cap_entry(cap, _fit(points, n_boot=N_BOOT, seed=SEED)))
    return rows


def _ci_summary(values: list[float], p_th: float) -> dict[str, Any]:
    if not values:
        return {"p_th_ci95": [None, None], "rel_ci_width": None, "resolved": False}
    lo, hi = (float(np.percentile(values, q)) for q in CI_PERCENTILES)
    width = (hi - lo) / p_th
    return {"p_th_ci95": [lo, hi], "rel_ci_width": width, "resolved": width < MAX_REL_CI_WIDTH}


def pin_retention(points: list[SweepPoint], seed: int) -> dict[str, Any]:
    """Pins-retained and pins-excluded per-round CIs of one condition at one seed."""
    kept = [pt for pt in points if pt.d >= D_MIN]
    reference = _fit(points, n_boot=N_BOOT, seed=seed)
    shots = np.array([pt.shots for pt in kept])
    rates = np.array([pt.p_l_shot for pt in kept])
    rng = np.random.default_rng(seed)
    retained: list[float] = []
    excluded: list[float] = []
    n_failed = 0
    for _ in range(N_BOOT):
        drawn = rng.binomial(shots, rates)
        resampled = [replace(pt, errors=int(e)) for pt, e in zip(kept, drawn, strict=True)]
        rep = _fit(resampled, n_boot=0, seed=seed)
        if not math.isfinite(rep.p_th):
            n_failed += 1
            continue
        retained.append(rep.p_th)
        if rep.converged:
            excluded.append(rep.p_th)
    with_pins = _ci_summary(retained, reference.p_th)
    if with_pins["p_th_ci95"] != [reference.p_th_ci[0], reference.p_th_ci[1]]:
        raise SystemExit(f"seed {seed}: pins-retained CI does not reproduce fit_threshold")
    if n_failed != reference.n_boot_failed:
        raise SystemExit(f"seed {seed}: failed-replicate count does not reproduce fit_threshold")
    return {
        "seed": seed,
        "p_th": reference.p_th,
        "n_boot": N_BOOT,
        "n_failed": n_failed,
        "n_pinned": len(retained) - len(excluded),
        "pins_retained": with_pins,
        "pins_excluded": _ci_summary(excluded, reference.p_th),
    }


def _check_published(conditions: list[dict[str, Any]], headline: dict[str, Any]) -> None:
    """Cap-0.4 rows must equal results/headline.json per_round_secondary."""
    published = {(c["decoder"], c["r_e"]): c["per_round_secondary"] for c in headline["conditions"]}
    for cond in conditions:
        row = next(r for r in cond["caps"] if r["cap"] == PUBLISHED_CAP)
        pub = published[(cond["decoder"], cond["r_e"])]
        for key in ("p_th", "p_th_ci95", "reason", "resolved"):
            if row[key] != pub[key]:
                raise SystemExit(
                    f"{cond['decoder']} r_e={cond['r_e']}: cap {PUBLISHED_CAP} {key} "
                    f"{row[key]!r} != headline.json {pub[key]!r}"
                )


def main() -> None:
    describe = _git_describe()
    if describe.endswith("-dirty") or describe == "unknown":
        raise SystemExit(f"git describe is {describe!r}; commit or stash changes first")
    conditions: list[dict[str, Any]] = []
    pins: list[dict[str, Any]] = []
    for path in DATA_FILES:
        points = load_sweep(ROOT / path)
        for decoder, r_e in available_conditions(points):
            cond_points = [pt for pt in points if pt.decoder == decoder and pt.r_e == r_e]
            caps = cap_sensitivity(cond_points)
            conditions.append({"data": path, "decoder": decoder, "r_e": r_e, "caps": caps})
            for row in caps:
                print(
                    f"{decoder} r_e={r_e:g} cap={row['cap']}: p_th={row['p_th']} "
                    f"ci={row['p_th_ci95']} reason={row['reason']}"
                )
            if (decoder, r_e) == PIN_CONDITION:
                for seed in PIN_SEEDS:
                    entry = pin_retention(cond_points, seed)
                    pins.append(entry)
                    print(
                        f"pins seed {seed}: retained {entry['pins_retained']['rel_ci_width']:.3f}"
                        f" excluded {entry['pins_excluded']['rel_ci_width']:.3f}"
                        f" pinned {entry['n_pinned']}/{N_BOOT}"
                    )
    _check_published(conditions, json.loads(HEADLINE_JSON.read_text()))
    doc = {
        "description": __doc__,
        "generated_by": "scripts/per_round_fit_sensitivity.py",
        "git_describe": describe,
        "data_sha256": {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in DATA_FILES},
        "software": {
            "python": sys.version.split()[0],
            "numpy": np.__version__,
            "scipy": scipy.__version__,
        },
        "parameters": {
            "metric": "round",
            "d_min": D_MIN,
            "window_factor": WINDOW_FACTOR,
            "n_boot": N_BOOT,
            "cap_seed": SEED,
            "caps": list(CAPS),
            "published_cap": PUBLISHED_CAP,
            "pin_condition": {"decoder": PIN_CONDITION[0], "r_e": PIN_CONDITION[1]},
            "pin_seeds": [PIN_SEEDS[0], PIN_SEEDS[-1]],
            "max_rel_ci_width": MAX_REL_CI_WIDTH,
        },
        "cap_sensitivity": conditions,
        "pin_retention": pins,
        "pin_summary": {
            "resolved_pins_retained": sum(e["pins_retained"]["resolved"] for e in pins),
            "resolved_pins_excluded": sum(e["pins_excluded"]["resolved"] for e in pins),
            "n_seeds": len(pins),
        },
    }
    JSON_OUT.write_text(json.dumps(doc, indent=2, allow_nan=False) + "\n")
    print(f"wrote {JSON_OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
