"""DEM-level noise statistics (Phase 2.1).

The circuit-wide *heralded fraction* -- the share of the detector-error-model's
total error mass that fires at least one herald detector -- is the honest answer
to "how much of the noise is erasure-converted." It is NOT ``r_e``: ``r_e`` is
the fraction of the *two-qubit-gate* budget converted, but measurement and reset
errors are never heralded and idle errors only when ``convert_idle`` is set, so
the circuit-wide fraction is lower. This is the number the README must quote
wherever ``r_e`` is described, so ``r_e`` is not silently compared to Wu et
al.'s ``R_e`` (fraction of *all* errors converted).
"""

from __future__ import annotations

import stim

from erasure_qec.circuits.builder import build
from erasure_qec.config import NoiseParams
from erasure_qec.noise.injector import ErasureInjector


def _contract_dem(circuit: stim.Circuit) -> stim.DetectorErrorModel:
    return circuit.detector_error_model(
        decompose_errors=True, flatten_loops=True, approximate_disjoint_errors=True
    ).flattened()


def _herald_detector_ids(dem: stim.DetectorErrorModel) -> set[int]:
    """Detector ids carrying the 4th sentinel coordinate ``= 1`` (herald bits)."""
    coords = dem.get_detector_coordinates()
    return {i for i, c in coords.items() if len(c) > 3 and c[3] == 1.0}


def heralded_fraction(params: NoiseParams, *, d: int = 5, rounds: int | None = None) -> float:
    """Fraction of total DEM error mass that fires >= 1 herald detector.

    Deterministic (a property of the decomposed DEM, no sampling). Sums each
    error mechanism's probability, weighting by whether any of its detectors is
    a herald detector, and divides by the total. ``rounds`` defaults to ``d``.
    """
    circuit = build(d, d if rounds is None else rounds, ErasureInjector(params))
    dem = _contract_dem(circuit)
    heralds = _herald_detector_ids(dem)
    total = 0.0
    heralded = 0.0
    for inst in dem:
        if inst.type != "error":
            continue
        prob = inst.args_copy()[0]
        dets = {t.val for t in inst.targets_copy() if t.is_relative_detector_id()}
        total += prob
        if dets & heralds:
            heralded += prob
    return heralded / total if total else float("nan")
