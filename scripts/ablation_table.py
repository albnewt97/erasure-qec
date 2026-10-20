"""Herald-vs-blind sub-threshold suppression table, direct from circuits.

The most reproducible result in the repo: for each fixed sub-threshold physical
error rate ``p`` in ``P_VALUES`` (0.5% and 1.0%; ``--p`` overrides), decode
identical shots with the herald-conditioned matcher and the blind matcher and
report the suppression ratio ``p_L(blind) / p_L(herald)`` per distance and
erasure fraction, and each decoder's own per-round distance-suppression factor
``Lambda_d = p_L(d) / p_L(d + 2)`` (reported only when that decoder has >= 50
observed errors at both ``d`` and ``d + 2``; see ``lambda_d``). This is
computed **directly from circuits** with a fixed seed — no CSV, no threshold
fit, and independent of how the noise budget is normalised across ``r_e`` — so
a reviewer can regenerate every number here from this one script. It is not
fast: one ``p`` at 100k shots took **~46 minutes** single-core (measured
2026-09-05 at p = 1.0%), dominated by the ``d = 9`` and ``d = 11`` cells; the
default two-``p`` run is projected at up to twice that. Writes
results/ablation.json (git describe, software versions, seed, shots and every
cell) and refuses to run on a dirty work tree.

    uv run python scripts/ablation_table.py                 # default 100k shots
    uv run python scripts/ablation_table.py --shots 200000 --p 0.01

Cells below the >= 50 observed-error gate (the same gate the Lambda figure uses)
are reported as a lower confidence bound ``> N×`` rather than a point ratio -- the
denominator is then a handful of errors and the point value is not trustworthy.
That bound combines two-sided-95% Wilson endpoints (one-sided 97.5% each) by
Bonferroni -- blind's lower limit over herald's upper limit -- giving a >= 95%
one-sided lower bound on the ratio *for that one cell*. See ``ablation_cell``.
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pymatching
import scipy
import stim

from erasure_qec.analysis.statistics import per_round_p_l, wilson_interval
from erasure_qec.circuits.builder import build
from erasure_qec.config import NoiseParams
from erasure_qec.decoding.herald_matching import (
    BlindMatchingDecoder,
    HeraldMatchingDecoder,
)
from erasure_qec.noise.injector import ErasureInjector

ROOT = Path(__file__).resolve().parent.parent
JSON_OUT = ROOT / "results" / "ablation.json"

# Fixed sub-threshold operating points and the grid the table sweeps.
P_VALUES = (0.005, 0.010)
R_ES = (0.5, 0.98)
DISTANCES = (3, 5, 7, 9, 11)
# One statistical standard across the repo: the same >= 50 observed-error gate
# the Lambda figure uses. Below it, the point ratio is not reported; only a
# Wilson-derived lower bound is (the denominator is a handful of errors).
_MIN_ERRORS = 50
# The stim detector-sampler seed of every cell.
SEED = 0


def _errors(pred: np.ndarray, obs: np.ndarray) -> int:
    return int((pred[:, 0] != obs[:, 0]).sum())


def ratio_and_bound(h_err: int, b_err: int, shots: int, d: int) -> dict[str, float]:
    """Point suppression ratio and its lower confidence bound, from raw counts.

    Pure arithmetic, split out of :func:`ablation_cell` so the bound can be
    tested over a grid of counts without re-running the Monte Carlo. See
    :func:`ablation_cell` for the confidence statement this bound supports.

    Two invariants hold for every ``(h_err, b_err, shots, d)`` and are pinned by
    ``tests/test_ablation_table.py``: ``ratio_lb <= ratio``, and ``ratio_lb`` is
    no larger than the superseded bound ``b_pl / h_pl_hi`` (which used blind's
    point estimate in the numerator and was therefore not a bound at all).
    """
    h_pl = per_round_p_l(h_err / shots, d)
    b_pl = per_round_p_l(b_err / shots, d)
    b_pl_lo = per_round_p_l(wilson_interval(b_err, shots)[0], d)
    h_pl_hi = per_round_p_l(wilson_interval(h_err, shots)[1], d)
    return {
        "ratio": b_pl / h_pl if h_pl > 0 else float("inf"),
        "ratio_lb": b_pl_lo / h_pl_hi if h_pl_hi > 0 else float("inf"),
        "h_pl": h_pl,
        "b_pl": b_pl,
        "b_pl_lo": b_pl_lo,
        "h_pl_hi": h_pl_hi,
    }


def ablation_cell(d: int, p: float, r_e: float, shots: int) -> dict[str, float]:
    """Decode ``shots`` identical shots with both decoders; return the point
    ratio and a conservative lower confidence bound on it.

    The bound puts blind's Wilson **lower** shot-rate limit in the numerator and
    herald's Wilson **upper** limit in the denominator, each mapped through the
    monotonically increasing :func:`per_round_p_l`::

        ratio_lb = per_round_p_l(wilson_lo(b_err, shots), d)
                   / per_round_p_l(wilson_hi(h_err, shots), d)

    **What confidence statement this supports.** ``wilson_interval`` is called at
    its default ``z = Z_95``, i.e. it returns a *two-sided 95%* interval, so each
    individual endpoint used here is a *one-sided 97.5%* limit. Bonferroni over
    the two endpoints gives a joint coverage of at least ``1 - 0.025 - 0.025 =
    95%``: with probability >= 95% both ``b_lo <= b_pl_true`` and ``h_hi >=
    h_pl_true`` hold simultaneously, and on that event ``ratio_lb <=
    b_pl_true / h_pl_true``. So the claim made is a **>= 95% one-sided lower
    confidence bound on the ratio for this one cell**. It is conservative in
    three ways:

    * Bonferroni is an inequality, not an equality -- the true joint coverage is
      somewhat above 95%, so the bound is looser than a 95% bound need be.
    * It makes no use of the pairing. The two decoders run on the *same* shots,
      so ``b_err`` and ``h_err`` are positively correlated and a paired analysis
      would be tighter; Bonferroni holds under *any* dependence, which is exactly
      why we may ignore that correlation rather than model it.
    * It is per-cell, **not simultaneous over the table**. Ten cells each at
      >= 95% do not jointly hold at 95%; a table-wide simultaneous claim would
      need a further Bonferroni split (>= 99.5% per endpoint).

    Two honest caveats. Combining two *one-sided 95%* limits -- what a naive
    reading of "95% Wilson" suggests -- would support only a >= 90% joint
    statement; it is the 97.5% endpoints that buy the 95% here. And Wilson
    coverage is nominal, not exact: at the single-digit error counts the
    high-``d``/high-``r_e`` cells produce, the "95%" is an approximation with no
    guaranteed actual coverage, so those cells' bounds are weak evidence and the
    caller must say so.
    """
    circ = build(d, d, ErasureInjector(NoiseParams(p=p, r_e=r_e)))
    dets, obs = circ.compile_detector_sampler(seed=SEED).sample(shots, separate_observables=True)
    h_err = _errors(HeraldMatchingDecoder.from_circuit(circ).decode_batch(dets), obs)
    b_err = _errors(BlindMatchingDecoder(circ).decode_batch(dets), obs)
    return {"h_err": h_err, "b_err": b_err, **ratio_and_bound(h_err, b_err, shots, d)}


def _fmt_ratio(cell: dict[str, float]) -> str:
    if cell["h_err"] >= _MIN_ERRORS:
        return f"{cell['ratio']:.1f}×"
    # Below the 50-error gate: a lower confidence bound, not a point ratio.
    return f"> {cell['ratio_lb']:.0f}× †"


def lambda_d(err_d: int, err_d2: int, shots: int, d: int) -> float | None:
    """One decoder's per-round ``Lambda_d = p_L(d) / p_L(d + 2)`` at a fixed ``p``.

    ``p_L`` is :func:`per_round_p_l` of the shot rate ``errors / shots`` with
    ``T = d`` (resp. ``d + 2``) rounds. ``None`` unless both distances have
    ``>= _MIN_ERRORS`` observed errors (the gate of the Lambda figure).
    """
    if err_d < _MIN_ERRORS or err_d2 < _MIN_ERRORS:
        return None
    return per_round_p_l(err_d / shots, d) / per_round_p_l(err_d2 / shots, d + 2)


def _finite(x: float) -> float | None:
    """JSON-safe float: NaN/inf become null."""
    return x if math.isfinite(x) else None


def cell_rows(
    cells: dict[tuple[float, float, int], dict[str, float]], shots: int
) -> list[dict[str, Any]]:
    """One JSON row per ``(p, r_e, d)`` cell, with both decoders' ``Lambda_d``."""
    rows: list[dict[str, Any]] = []
    for (p, r_e, d), c in sorted(cells.items()):
        partner = cells.get((p, r_e, d + 2))
        lam: dict[str, float | None] = {"herald": None, "blind": None}
        if partner is not None:
            for name, key in (("herald", "h_err"), ("blind", "b_err")):
                lam[name] = lambda_d(int(c[key]), int(partner[key]), shots, d)
        rows.append(
            {
                "p": p,
                "r_e": r_e,
                "d": d,
                "shots": shots,
                "h_err": int(c["h_err"]),
                "b_err": int(c["b_err"]),
                "gated": c["h_err"] < _MIN_ERRORS,
                **{k: _finite(c[k]) for k in ("ratio", "ratio_lb", "h_pl", "b_pl")},
                **{k: _finite(c[k]) for k in ("b_pl_lo", "h_pl_hi")},
                "lambda_pair": [d, d + 2] if partner is not None else None,
                "lambda_herald": lam["herald"],
                "lambda_blind": lam["blind"],
            }
        )
    return rows


def _git_describe() -> str:
    out = subprocess.run(
        ["git", "describe", "--always", "--dirty"],
        capture_output=True,
        text=True,
        check=False,
        cwd=ROOT,
    )
    return out.stdout.strip() or "unknown"


def payload(
    cells: dict[tuple[float, float, int], dict[str, float]],
    p_values: tuple[float, ...],
    shots: int,
    describe: str,
) -> dict[str, Any]:
    """The results/ablation.json document."""
    return {
        "description": "herald vs blind suppression ratio and per-decoder per-round "
        "Lambda_d at fixed sub-threshold p, decoded on identical shots",
        "generated_by": "scripts/ablation_table.py",
        "git_describe": describe,
        "software": {
            "python": sys.version.split()[0],
            "stim": stim.__version__,
            "pymatching": pymatching.__version__,
            "numpy": np.__version__,
            "scipy": scipy.__version__,
        },
        "parameters": {
            "p_values": list(p_values),
            "r_es": list(R_ES),
            "distances": list(DISTANCES),
            "shots": shots,
            "seed": SEED,
            "min_errors": _MIN_ERRORS,
            "ratio": "per-round p_L(blind) / p_L(herald); below the gate (herald errors "
            "< min_errors) only ratio_lb is quoted",
            "ratio_lb": "per-round blind Wilson lower / herald Wilson upper (two-sided 95% "
            "endpoints); >= 95% one-sided lower bound per cell (Bonferroni)",
            "lambda": "per-round p_L(d) / p_L(d+2) of one decoder; null unless that "
            "decoder has >= min_errors errors at both d and d+2",
        },
        "cells": cell_rows(cells, shots),
    }


def _print_table(
    cells: dict[tuple[float, float, int], dict[str, float]], p: float, shots: int
) -> None:
    header = f"p_L(blind)/p_L(herald) at p={p:.1%}, {shots:,} shots, seed={SEED}"
    print(header)
    print("-" * len(header))
    print("| `d` | " + " | ".join(f"`r_e = {r}`" for r in R_ES) + " |")
    print("|---|" + "---|" * len(R_ES))
    for d in DISTANCES:
        row = " | ".join(_fmt_ratio(cells[(p, r_e, d)]) for r_e in R_ES)
        print(f"| {d} | {row} |")
    print(
        f"\n† below the {_MIN_ERRORS}-observed-error gate (the same CONSTANT as "
        "the Lambda figure, but not the same rule: that figure requires "
        f"{_MIN_ERRORS} errors\n  at both d and d+2, while this table gates on the "
        "herald count alone -- the fragile denominator here): reported as a "
        "lower confidence bound `> N×`, not a point ratio.\n  Bound = per-round "
        "blind Wilson LOWER / herald Wilson UPPER "
        "(two-sided-95% endpoints, i.e. one-sided 97.5% each);\n  Bonferroni "
        "makes it a >= 95% one-sided lower bound on that cell's ratio -- "
        "per-cell, not simultaneous over the table."
    )

    print("\nraw counts (herald_err / blind_err):")
    for d in DISTANCES:
        parts = " | ".join(
            f"r_e={r_e}: {int(cells[(p, r_e, d)]['h_err'])}/{int(cells[(p, r_e, d)]['b_err'])}"
            for r_e in R_ES
        )
        print(f"  d={d:2d}  {parts}")

    print("\nper-cell detail (point ratio, bound, per-round rates):")
    for r_e in R_ES:
        for d in DISTANCES:
            c = cells[(p, r_e, d)]
            print(
                f"  r_e={r_e:<5} d={d:2d}  h_err={int(c['h_err']):6d} "
                f"b_err={int(c['b_err']):6d}  ratio={c['ratio']:8.2f}  "
                f"ratio_lb={c['ratio_lb']:8.2f}  "
                f"b_pl_lo={c['b_pl_lo']:.3e}  h_pl_hi={c['h_pl_hi']:.3e}"
            )
    # Wilson 95% CI on the herald shot-rate at the highest-d, highest-r_e cell,
    # to make the low-statistics caveat concrete.
    worst = cells[(p, R_ES[-1], DISTANCES[-1])]
    lo, hi = wilson_interval(int(worst["h_err"]), shots)
    print(
        f"\nherald shot-rate CI at d={DISTANCES[-1]}, r_e={R_ES[-1]}: "
        f"{worst['h_err'] / shots:.2e} in [{lo:.2e}, {hi:.2e}] (95% Wilson)"
    )

    print("\nper-round Lambda_d = p_L(d)/p_L(d+2) (herald / blind; — below the gate):")
    rows = {(r["r_e"], r["d"]): r for r in cell_rows(cells, shots) if r["p"] == p}
    for r_e in R_ES:
        parts = []
        for d in DISTANCES[:-1]:
            r = rows[(r_e, d)]
            h, b = r["lambda_herald"], r["lambda_blind"]
            parts.append(
                f"d={d}->{d + 2}: {'—' if h is None else f'{h:.2f}'} / "
                f"{'—' if b is None else f'{b:.2f}'}"
            )
        print(f"  r_e={r_e:<5} " + "  ".join(parts))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shots", type=int, default=100_000)
    parser.add_argument("--p", type=float, nargs="+", default=list(P_VALUES))
    parser.add_argument("--out", type=Path, default=JSON_OUT)
    args = parser.parse_args()
    describe = _git_describe()
    if describe.endswith("-dirty") or describe == "unknown":
        raise SystemExit(f"git describe is {describe!r}; commit or stash changes first")
    p_values = tuple(args.p)

    cells = {
        (p, r_e, d): ablation_cell(d, p, r_e, args.shots)
        for p in p_values
        for r_e in R_ES
        for d in DISTANCES
    }
    for i, p in enumerate(p_values):
        if i:
            print("\n")
        _print_table(cells, p, args.shots)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    doc = payload(cells, p_values, args.shots, describe)
    out.write_text(json.dumps(doc, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
