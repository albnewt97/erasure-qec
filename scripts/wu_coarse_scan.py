"""Coarse scan locating the (5,7) crossing for the Wu R_e = 0.5 and 0.9 grids (P5a).

For each R_e in (0.5, 0.9), with p_m = 0, rounds = d and the ``wu_herald_mwpm``
decoder only: at d in {5, 7} and 9 log-spaced p in [0.008, 0.05], samples
``SHOTS`` shots of the Wu sampling circuit (sampler seed ``SEED + i``, i the
index of the (R_e, d, p) point in loop order) and counts logical failures.
``c`` is the per-shot (metric "shot") pairwise crossing of the d = 5 and d = 7
curves from :func:`erasure_qec.analysis.crossing.pairwise_crossing`; the P5a
grid is ``round(c * f, 4)`` for f in ``GRID_FACTORS``.

Planning-only scratch: nothing is written to data/ and the counts are not used
in any result. Prints the counts, c, the grid and the runtime; exits with
status 1 if either R_e has no crossing. Run:

    uv run python scripts/wu_coarse_scan.py
"""

from __future__ import annotations

import math
import sys
import time

import numpy as np
import pymatching
import sinter
import stim

from erasure_qec.analysis.crossing import pairwise_crossing
from erasure_qec.analysis.statistics import SweepPoint
from erasure_qec.decoding.sinter_adapter import WU_DECODERS, wu_tasks_parallel
from erasure_qec.noise.wu import WuNoiseParams

R_ES = (0.5, 0.9)
DISTANCES = (5, 7)
P_MIN, P_MAX, P_POINTS = 0.008, 0.05, 9
SHOTS = 3000
SEED = 2026
DECODER = "wu_herald_mwpm"
PROCESSES = 8
GRID_FACTORS = (0.800, 0.845, 0.890, 0.935, 0.980, 1.025, 1.070, 1.115, 1.160, 1.205, 1.250)


def p_values() -> list[float]:
    """The 9 log-spaced scan points, ``numpy.geomspace(P_MIN, P_MAX, P_POINTS)``."""
    return [float(p) for p in np.geomspace(P_MIN, P_MAX, P_POINTS)]


def grid_from_crossing(c: float) -> list[float]:
    """The P5a grid ``round(c * f, 4)`` for f in ``GRID_FACTORS``."""
    return [round(c * f, 4) for f in GRID_FACTORS]


def _failures(task: sinter.Task, seed: int) -> int:
    """Logical failures of ``DECODER`` on ``SHOTS`` shots sampled from ``task.circuit``."""
    assert task.circuit is not None and task.detector_error_model is not None
    sampler = task.circuit.compile_detector_sampler(seed=seed)
    dets, obs = sampler.sample(SHOTS, separate_observables=True, bit_packed=True)
    compiled = WU_DECODERS[DECODER].compile_decoder_for_dem(dem=task.detector_error_model)
    pred = compiled.decode_shots_bit_packed(bit_packed_detection_event_data=dets)
    return int(np.any(pred != obs, axis=1).sum())


def main() -> int:
    t0 = time.perf_counter()
    print(
        f"# versions (stim pymatching sinter): {stim.__version__} {pymatching.__version__} "
        f"{sinter.__version__}"
    )
    print(
        f"# decoder {DECODER}; p_m = 0; rounds = d; {SHOTS} shots per point; "
        f"sampler seed {SEED} + point index"
    )
    ps = p_values()
    specs = [
        (d, WuNoiseParams(p, r_e), {"d": d, "p": p, "R_e": r_e})
        for r_e in R_ES
        for d in DISTANCES
        for p in ps
    ]
    tasks = wu_tasks_parallel(specs, processes=PROCESSES)
    t_build = time.perf_counter() - t0

    points: dict[float, list[SweepPoint]] = {r_e: [] for r_e in R_ES}
    print("R_e d p shots errors P_L_shot")
    for i, ((d, params, _), task) in enumerate(zip(specs, tasks, strict=True)):
        errors = _failures(task, SEED + i)
        points[params.R_e].append(SweepPoint(DECODER, d, d, params.p, params.R_e, SHOTS, errors))
        print(f"{params.R_e} {d} {params.p:.6f} {SHOTS} {errors} {errors / SHOTS:.4f}")

    missing = False
    for r_e in R_ES:
        est = pairwise_crossing(points[r_e], 5, 7, metric="shot")
        print(f"\nR_e = {r_e}: (5,7) per-shot crossing c = {est.value!r}")
        print(f"  bracket = {est.bracket}, n_sign_changes = {est.n_sign_changes}")
        if math.isnan(est.value):
            print("  NO CROSSING FOUND")
            missing = True
            continue
        print(f"  grid = {grid_from_crossing(est.value)}")

    print(f"\nruntime: task build {t_build:.1f} s, total {time.perf_counter() - t0:.1f} s")
    return 1 if missing else 0


if __name__ == "__main__":
    sys.exit(main())
