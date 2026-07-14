"""Rotated surface-code Z-memory circuit builder (PLAN.md §3).

This is the only module that emits ``stim.Circuit``. It assembles the
noiseless skeleton — the 8-TICK round structure (§3.1) using the hook-safe
CX schedule (§3.2) and the exact detector algebra of §3.4 — and defers all
noise to a :class:`~erasure_qec.noise.injector.NoiseInjector`. Herald record
bits (from the injector's ``HERALDED_ERASE`` calls, or from the
``probe_erasures`` debug hook) each get an immediate sentinel
``DETECTOR(x, y, t, 1)`` per §3.4/§5.
"""

from collections.abc import Sequence

import stim

from erasure_qec.circuits.layout import (
    Basis,
    ancilla_coords,
    data_coords,
    logical_z_support,
    plaquette_neighbors,
)
from erasure_qec.circuits.scheduling import HOOK_SAFE_SCHEDULE, Schedule
from erasure_qec.noise.injector import NoiseInjector

# Fixed probability for the §7.1 DEM hand-verification probe hook. This is a
# debug/testing affordance, independent of any NoiseInjector.
PROBE_ERASURE_PROBABILITY = 0.25


def _sort_key(q: complex) -> tuple[float, float]:
    """Stable ordering for qubit coordinates (row-major: y then x)."""
    return (q.imag, q.real)


def num_detectors(d: int, rounds: int) -> int:
    """Closed-form syndrome-detector count for a ``rounds``-round Z-memory (§3.5).

    ``(d**2 - 1)/2`` round-0 Z detectors, ``(rounds - 1)(d**2 - 1)`` bulk
    detectors, and ``(d**2 - 1)/2`` time-closing Z detectors. Does not include
    herald detectors, which depend on the noise model and probe hooks.
    """
    if rounds < 1:
        raise ValueError(f"rounds must be >= 1, got {rounds}")
    n_anc = d**2 - 1
    return n_anc // 2 + (rounds - 1) * n_anc + n_anc // 2


def build(
    d: int,
    rounds: int,
    injector: NoiseInjector,
    *,
    schedule: Schedule = HOOK_SAFE_SCHEDULE,
    probe_erasures: Sequence[tuple[complex, int]] = (),
    mid_round_probe_erasures: Sequence[tuple[complex, int, int]] = (),
    probe_q: float = 0.25,
) -> stim.Circuit:
    """Build the distance-``d``, ``rounds``-round rotated surface-code Z-memory.

    Args:
        d: Code distance (odd, >= 1); determines the d x d data-qubit grid.
        rounds: Number of stabilizer-measurement rounds T (>= 1).
        injector: Noise injector called at each fault location; ``NullInjector``
            yields the noiseless circuit.
        schedule: Per-basis CX slot order into the canonical [NE, NW, SE, SW]
            neighbor order. Defaults to the hook-safe schedule (§3.2).
        probe_erasures: Debug hook (§7.1): ``(qubit, after_round)`` pairs, each
            hand-inserting one ``HERALDED_ERASE(probe_q)`` on ``qubit`` right
            after round ``after_round``'s ``MR`` tick (before the next round's
            reset, or before the closing measurement if ``after_round`` is the
            last round). Independent of ``injector``.
        mid_round_probe_erasures: Debug hook (§7.3 criterion 6): ``(qubit,
            round_index, after_cx_layer)`` triples, each hand-inserting one
            ``HERALDED_ERASE(probe_q)`` on ``qubit`` within round ``round_index``,
            immediately after CX layer ``after_cx_layer`` (0-indexed: 0..3;
            ``1`` means "between CX layers 2 and 3"). Independent of
            ``injector``.
        probe_q: Erasure probability for all probe ``HERALDED_ERASE``
            instructions (default 0.25 for §7.1 hand verification).

    Returns:
        The assembled ``stim.Circuit``.
    """
    if rounds < 1:
        raise ValueError(f"rounds must be >= 1, got {rounds}")

    data = sorted(data_coords(d), key=_sort_key)
    ancilla_basis = ancilla_coords(d)
    ancillas = sorted(ancilla_basis, key=_sort_key)
    x_ancillas = [a for a in ancillas if ancilla_basis[a] is Basis.X]
    z_ancillas = [a for a in ancillas if ancilla_basis[a] is Basis.Z]
    neighbors = {a: plaquette_neighbors(a, d) for a in ancillas}

    index: dict[complex, int] = {q: i for i, q in enumerate(data + ancillas)}
    rev_index: dict[int, complex] = {i: q for q, i in index.items()}
    circuit = stim.Circuit()

    probes_by_round: dict[int, list[complex]] = {}
    for qubit, after_round in probe_erasures:
        probes_by_round.setdefault(after_round, []).append(qubit)

    mid_round_probes_by_round: dict[int, dict[int, list[complex]]] = {}
    for qubit, round_index, after_cx_layer in mid_round_probe_erasures:
        by_layer = mid_round_probes_by_round.setdefault(round_index, {})
        by_layer.setdefault(after_cx_layer, []).append(qubit)

    # --- Preamble: coordinates + data reset in |0> (Z basis). ---
    for q in data + ancillas:
        circuit.append("QUBIT_COORDS", [index[q]], [q.real, q.imag])
    data_idx = [index[q] for q in data]
    circuit.append("R", data_idx)
    injector.after_reset(circuit, data_idx)
    circuit.append("TICK")

    meas_count = 0
    prev_round_meas: dict[complex, int] = {}

    for t in range(rounds):
        if t > 0:
            circuit.append("SHIFT_COORDS", [], (0, 0, 1))
        meas_count = _emit_round(
            circuit,
            injector,
            index,
            rev_index,
            ancillas,
            x_ancillas,
            z_ancillas,
            neighbors,
            schedule,
            meas_count,
            mid_round_probes_by_round.get(t, {}),
            probe_q,
        )
        # Assign measurement-record indices in the MR order (== `ancillas`).
        cur_round_meas: dict[complex, int] = {}
        for a in ancillas:
            cur_round_meas[a] = meas_count
            meas_count += 1
        _emit_syndrome_detectors(
            circuit, t, ancillas, z_ancillas, cur_round_meas, prev_round_meas, meas_count
        )
        prev_round_meas = cur_round_meas
        meas_count = _emit_probe_erasures(
            circuit, probes_by_round.get(t, []), index, rev_index, meas_count, probe_q
        )

    # --- Time close: measure data, closing Z detectors, logical observable. ---
    circuit.append("SHIFT_COORDS", [], (0, 0, 1))
    data_idx = [index[q] for q in data]
    injector.before_measure(circuit, data_idx)
    circuit.append("M", data_idx)
    data_meas: dict[complex, int] = {}
    for q in data:
        data_meas[q] = meas_count
        meas_count += 1
    _emit_closing_detectors(
        circuit, z_ancillas, neighbors, data_meas, prev_round_meas, meas_count
    )
    obs = [
        stim.target_rec(data_meas[q] - meas_count) for q in logical_z_support(d)
    ]
    circuit.append("OBSERVABLE_INCLUDE", obs, 0)
    return circuit


def _emit_round(
    circuit: stim.Circuit,
    injector: NoiseInjector,
    index: dict[complex, int],
    rev_index: dict[int, complex],
    ancillas: Sequence[complex],
    x_ancillas: Sequence[complex],
    z_ancillas: Sequence[complex],
    neighbors: dict[complex, list[complex | None]],
    schedule: Schedule,
    meas_count: int,
    mid_round_probes: dict[int, list[complex]],
    probe_q: float,
) -> int:
    """Emit one 8-TICK stabilizer round (§3.1).

    Returns the updated ``meas_count``, accounting for any herald record bits
    produced by the injector during the CX layers (§5), or by
    ``mid_round_probes`` (§7.3 criterion 6): qubits to hand-erase right after
    the CX layer keyed by index (0..3).
    """
    all_indices = list(index.values())
    anc_idx = [index[a] for a in ancillas]
    x_idx = [index[a] for a in x_ancillas]

    # TICK 1: reset ancillas.
    circuit.append("R", anc_idx)
    injector.after_reset(circuit, anc_idx)
    circuit.append("TICK")

    # TICK 2: H on X-ancillas.
    circuit.append("H", x_idx)
    circuit.append("TICK")

    # TICKs 3-6: four CX layers. X-ancillas are controls, Z-ancillas targets.
    for layer in range(4):
        pairs: list[tuple[int, int]] = []
        for a in x_ancillas:
            nb = neighbors[a][schedule[Basis.X][layer]]
            if nb is not None:
                pairs.append((index[a], index[nb]))
        for a in z_ancillas:
            nb = neighbors[a][schedule[Basis.Z][layer]]
            if nb is not None:
                pairs.append((index[nb], index[a]))
        circuit.append("CX", [q for pair in pairs for q in pair])
        heralded = injector.on_two_qubit_gate(circuit, pairs)
        meas_count = _emit_herald_detectors(circuit, heralded, rev_index, meas_count)
        active = {q for pair in pairs for q in pair}
        injector.on_idle(circuit, [i for i in all_indices if i not in active])
        circuit.append("TICK")
        meas_count = _emit_probe_erasures(
            circuit, mid_round_probes.get(layer, []), index, rev_index, meas_count, probe_q
        )

    # TICK 7: H on X-ancillas.
    circuit.append("H", x_idx)
    circuit.append("TICK")

    # TICK 8: measure + reset all ancillas (order fixes the record layout).
    injector.before_measure(circuit, anc_idx)
    circuit.append("MR", anc_idx)
    return meas_count


def _emit_herald_detectors(
    circuit: stim.Circuit,
    heralded_qubits: Sequence[int],
    rev_index: dict[int, complex],
    meas_count: int,
) -> int:
    """Append a sentinel ``DETECTOR(x, y, 0, 1)`` for each just-appended herald bit.

    ``heralded_qubits`` must list qubit indices in the exact order their
    herald bits landed in the measurement record (immediately after the
    ``HERALDED_ERASE`` instruction that produced them, before anything else
    appends to the record). Returns the updated ``meas_count``.
    """
    n = len(heralded_qubits)
    for k, q in enumerate(heralded_qubits):
        coord = rev_index[q]
        circuit.append(
            "DETECTOR",
            [stim.target_rec(k - n)],
            (coord.real, coord.imag, 0, 1),
        )
    return meas_count + n


def _emit_probe_erasures(
    circuit: stim.Circuit,
    qubits: Sequence[complex],
    index: dict[complex, int],
    rev_index: dict[int, complex],
    meas_count: int,
    probe_q: float,
) -> int:
    """Hand-inserted ``HERALDED_ERASE`` probes for the §7.1 worksheet.

    Independent of the injector: fires at probability ``probe_q``.
    """
    if not qubits:
        return meas_count
    targets = [index[q] for q in qubits]
    circuit.append("HERALDED_ERASE", targets, probe_q)
    return _emit_herald_detectors(circuit, targets, rev_index, meas_count)


def _emit_syndrome_detectors(
    circuit: stim.Circuit,
    round_index: int,
    ancillas: Sequence[complex],
    z_ancillas: Sequence[complex],
    cur: dict[complex, int],
    prev: dict[complex, int],
    meas_count: int,
) -> None:
    """Round-0 Z-only detectors, or bulk XOR detectors for rounds >= 1 (§3.4)."""
    if round_index == 0:
        for a in z_ancillas:
            circuit.append(
                "DETECTOR",
                [stim.target_rec(cur[a] - meas_count)],
                (a.real, a.imag, 0),
            )
    else:
        for a in ancillas:
            circuit.append(
                "DETECTOR",
                [
                    stim.target_rec(cur[a] - meas_count),
                    stim.target_rec(prev[a] - meas_count),
                ],
                (a.real, a.imag, 0),
            )


def _emit_closing_detectors(
    circuit: stim.Circuit,
    z_ancillas: Sequence[complex],
    neighbors: dict[complex, list[complex | None]],
    data_meas: dict[complex, int],
    prev_round_meas: dict[complex, int],
    meas_count: int,
) -> None:
    """Time-closing Z detectors: (XOR of adjacent data) XOR last syndrome (§3.4)."""
    for a in z_ancillas:
        recs = [
            stim.target_rec(data_meas[nb] - meas_count)
            for nb in neighbors[a]
            if nb is not None
        ]
        recs.append(stim.target_rec(prev_round_meas[a] - meas_count))
        circuit.append("DETECTOR", recs, (a.real, a.imag, 0))
