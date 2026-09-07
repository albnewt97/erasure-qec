"""DEM-level noise statistics (Phase 2.1).

The circuit-wide *heralded fraction* -- the share of the detector-error-model's
total error mass that fires at least one herald detector -- is the honest answer
to "how much of the noise is erasure-converted." Error mass here excludes the
no-op herald branches: ``HERALDED_ERASE`` heralds on the identity branch too, and
a mechanism that flips only a herald detector causes no syndrome and no logical
error, so it is not part of the error budget on either side of the ratio. It is
NOT ``r_e``: ``r_e`` is
the fraction of the *two-qubit-gate* budget converted, but measurement and reset
errors are never heralded and idle errors only when ``convert_idle`` is set, so
the circuit-wide fraction is lower. This is the number the README must quote
wherever ``r_e`` is described.

Wu et al. (arXiv:2201.03540) define ``R_e`` on the *same* per-two-qubit-gate
axis as ``r_e`` -- each 2q gate takes a Pauli error with probability
``p(1 - R_e)`` or an erasure with probability ``p R_e``. Their ``R_e`` also
equals their circuit-wide erasure fraction, but only because they set
``p_m = 0`` and neglect idle errors, leaving 2q gates as the sole error
locations. Ours does not, because we add unheralded measurement, reset and idle
noise at rate ``p``. So ``r_e`` and ``R_e`` are directly comparable as
conversion axes; it is the *circuit-wide* fractions that are not.
"""

from __future__ import annotations

from typing import NamedTuple

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


class DemErrorMass(NamedTuple):
    """Total DEM error-mechanism probability, split by herald/no-op status.

    ``no_op`` is the pathological bucket: ``HERALDED_ERASE`` heralds on all four
    Pauli branches including identity, so the decomposed DEM carries one
    mechanism per herald detector whose *only* target is that herald detector --
    no syndrome detector, no logical observable. Those branches announce "this
    qubit was erased" and then flip nothing.
    """

    unheralded: float
    heralded: float
    no_op: float


def dem_error_mass(params: NoiseParams, *, d: int = 5, rounds: int | None = None) -> DemErrorMass:
    """Partition the flattened decomposed DEM's error mass into three buckets.

    Deterministic (a property of the DEM, no sampling). ``rounds`` defaults to
    ``d``. A mechanism is *heralded* if any of its detector targets is a herald
    detector, and a *no-op* if it is heralded and has no non-herald detector and
    no logical observable target.
    """
    circuit = build(d, d if rounds is None else rounds, ErasureInjector(params))
    dem = _contract_dem(circuit)
    heralds = _herald_detector_ids(dem)
    unheralded = heralded = no_op = 0.0
    for inst in dem:
        if inst.type != "error":
            continue
        prob = inst.args_copy()[0]
        dets: set[int] = set()
        observables: set[int] = set()
        for target in inst.targets_copy():
            if target.is_relative_detector_id():
                dets.add(target.val)
            elif target.is_logical_observable_id():
                observables.add(target.val)
        if not dets & heralds:
            unheralded += prob
        elif (dets - heralds) or observables:
            heralded += prob
        else:
            no_op += prob
    return DemErrorMass(unheralded=unheralded, heralded=heralded, no_op=no_op)


def heralded_fraction(
    params: NoiseParams,
    *,
    d: int = 5,
    rounds: int | None = None,
    count_no_op_heralds: bool = False,
) -> float:
    """Fraction of the DEM's *error* mass that fires >= 1 herald detector.

    Deterministic (a property of the decomposed DEM, no sampling). ``rounds``
    defaults to ``d``.

    No-op herald branches are excluded from **both** the numerator and the
    denominator by default. ``HERALDED_ERASE`` replaces the qubit with ``I/2``
    and heralds on all four Pauli branches, identity included, so the decomposed
    DEM contains one mechanism per herald detector whose only target is that
    herald detector. Such a mechanism causes no syndrome and no logical error:
    it is not part of the error budget, and counting it in both the heralded
    mass and the total inflates the reported fraction. See
    :func:`dem_error_mass`.

    Pass ``count_no_op_heralds=True`` for the original (pre-correction)
    definition, which counted those branches in both. It is retained only so the
    historical numbers in docs/AUDIT.md stay reproducible; it overstates how much
    of the circuit's error budget is erasure-converted.
    """
    mass = dem_error_mass(params, d=d, rounds=rounds)
    if count_no_op_heralds:
        numerator = mass.heralded + mass.no_op
        total = mass.unheralded + mass.heralded + mass.no_op
    else:
        numerator = mass.heralded
        total = mass.unheralded + mass.heralded
    return numerator / total if total else float("nan")
