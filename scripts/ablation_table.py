"""Herald-vs-blind sub-threshold suppression table, direct from circuits.

The most reproducible result in the repo: for a fixed sub-threshold physical
error rate ``p``, decode identical shots with the herald-conditioned matcher and
the blind matcher and report the suppression ratio ``p_L(blind) / p_L(herald)``
per distance and erasure fraction. This is computed **directly from circuits**
with a fixed seed — no CSV, no threshold fit, and independent of how the noise
budget is normalised across ``r_e`` — so a reviewer can regenerate every number
here in minutes.

    uv run python scripts/ablation_table.py                 # default 100k shots
    uv run python scripts/ablation_table.py --shots 200000

Cells below the >= 50 observed-error gate (the same gate the Lambda figure uses)
are reported as a Wilson lower bound ``> N×`` rather than a point ratio -- the
denominator is then a handful of errors and the point value is not trustworthy.
"""

from __future__ import annotations

import argparse

import numpy as np

from erasure_qec.analysis.statistics import per_round_p_l, wilson_interval
from erasure_qec.circuits.builder import build
from erasure_qec.config import NoiseParams
from erasure_qec.decoding.herald_matching import (
    BlindMatchingDecoder,
    HeraldMatchingDecoder,
)
from erasure_qec.noise.injector import ErasureInjector

# Fixed sub-threshold operating point and the grid the table sweeps.
FIXED_P = 0.01
R_ES = (0.5, 0.98)
DISTANCES = (3, 5, 7, 9, 11)
# One statistical standard across the repo: the same >= 50 observed-error gate
# the Lambda figure uses. Below it, the point ratio is not reported; only a
# Wilson-derived lower bound is (the denominator is a handful of errors).
_MIN_ERRORS = 50


def _errors(pred: np.ndarray, obs: np.ndarray) -> int:
    return int((pred[:, 0] != obs[:, 0]).sum())


def ablation_cell(d: int, p: float, r_e: float, shots: int) -> dict[str, float]:
    """Decode ``shots`` identical shots with both decoders; return the ratio and
    a Wilson lower bound on it (herald's 95% *upper* shot-rate limit in the
    denominator, blind held at its well-measured point estimate)."""
    circ = build(d, d, ErasureInjector(NoiseParams(p=p, r_e=r_e)))
    dets, obs = circ.compile_detector_sampler(seed=0).sample(
        shots, separate_observables=True
    )
    h_err = _errors(HeraldMatchingDecoder.from_circuit(circ).decode_batch(dets), obs)
    b_err = _errors(BlindMatchingDecoder(circ).decode_batch(dets), obs)
    h_pl = per_round_p_l(h_err / shots, d)
    b_pl = per_round_p_l(b_err / shots, d)
    ratio = b_pl / h_pl if h_pl > 0 else float("inf")
    h_pl_hi = per_round_p_l(wilson_interval(h_err, shots)[1], d)
    ratio_lb = b_pl / h_pl_hi if h_pl_hi > 0 else float("inf")
    return {"ratio": ratio, "ratio_lb": ratio_lb, "h_err": h_err, "b_err": b_err,
            "h_pl": h_pl, "b_pl": b_pl}


def _fmt_ratio(cell: dict[str, float]) -> str:
    if cell["h_err"] >= _MIN_ERRORS:
        return f"{cell['ratio']:.1f}×"
    # Below the 50-error gate: a Wilson lower bound, not a point ratio.
    return f"> {cell['ratio_lb']:.0f}× †"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shots", type=int, default=100_000)
    args = parser.parse_args()

    cells = {
        (r_e, d): ablation_cell(d, FIXED_P, r_e, args.shots)
        for r_e in R_ES
        for d in DISTANCES
    }

    header = f"p_L(blind)/p_L(herald) at p={FIXED_P:.1%}, {args.shots:,} shots, seed=0"
    print(header)
    print("-" * len(header))
    print("| `d` | " + " | ".join(f"`r_e = {r}`" for r in R_ES) + " |")
    print("|---|" + "---|" * len(R_ES))
    for d in DISTANCES:
        row = " | ".join(_fmt_ratio(cells[(r_e, d)]) for r_e in R_ES)
        print(f"| {d} | {row} |")
    print(f"\n† below the {_MIN_ERRORS}-observed-error gate (same as the Lambda "
          "figure): reported as a Wilson lower bound `> N×`, not a point ratio.")

    print("\nraw counts (herald_err / blind_err):")
    for d in DISTANCES:
        parts = " | ".join(
            f"r_e={r_e}: {int(cells[(r_e, d)]['h_err'])}/{int(cells[(r_e, d)]['b_err'])}"
            for r_e in R_ES
        )
        print(f"  d={d:2d}  {parts}")
    # Wilson 95% CI on the herald shot-rate at the highest-d, highest-r_e cell,
    # to make the low-statistics caveat concrete.
    worst = cells[(R_ES[-1], DISTANCES[-1])]
    lo, hi = wilson_interval(int(worst["h_err"]), args.shots)
    print(
        f"\nherald shot-rate CI at d={DISTANCES[-1]}, r_e={R_ES[-1]}: "
        f"{worst['h_err'] / args.shots:.2e} in [{lo:.2e}, {hi:.2e}] (95% Wilson)"
    )


if __name__ == "__main__":
    main()
