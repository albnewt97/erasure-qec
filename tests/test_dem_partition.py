"""Frozen hand-verified ground truth for the Design-2 DEM partition (PLAN.md §6/§7).

Every literal index below was verified by hand and reconciled against Stim's
emitted DEM in ``docs/dem_worksheet.md`` (Probe A run 2026-07-08; Probes B & C
run 2026-07-09; see the three "RECONCILIATION RECORD" blocks). Magic numbers are
a *feature* here: the worksheet's whole value is that these concrete detector
indices were checked against the physics. Each is commented with its coordinate
meaning.

Design 2: the partition consumes the DEM built with ``decompose_errors=True,
flatten_loops=True, approximate_disjoint_errors=True`` (this is what
``partition_dem`` does internally), splits each mechanism at ``^`` separators
into graphlike components, drops the herald component, and deduplicates the
remaining per-component edges. The Y mechanism therefore contributes no new
edge — it re-supplies the same components as X and Z.

If any assertion here fails, the *test* or the worksheet is wrong, not the
source: stop and re-verify by hand rather than editing ``dem_partition.py``.

Parity in this repo (from the worksheet): (2,2) and (4,4) are Z-checks;
(4,2) and (2,4) are X-checks; (2,0) is a top-boundary X-check.
"""

import pytest
import stim

from erasure_qec.circuits.builder import build
from erasure_qec.decoding.dem_partition import (
    PartitionedDEM,
    partition_dem,
    partition_flattened_dem,
)
from erasure_qec.noise.injector import NullInjector


def _edge_set(
    part: PartitionedDEM, herald: int
) -> set[tuple[tuple[int, ...], tuple[int, ...]]]:
    """Canonical {(sorted dets, obs_mask)} set for one herald's edges."""
    return {(e.dets, e.obs_mask) for e in part.herald_table[herald]}


# ---------------------------------------------------------------------------
# Probe A — erasure on the center data qubit (3,3), between rounds 0 and 1.
# Interior qubit: each Pauli component fires a pair of round-1 detectors.
# ---------------------------------------------------------------------------

# Indices VERIFIED BY HAND 2026-07-08 (docs/dem_worksheet.md, Probe A record).
A_HERALD = 4  # herald detector, coords (3, 3, 0, 1)
A_DZ_22, A_DZ_44 = 6, 11  # Z-checks (2,2) & (4,4), round-1 detectors -> X component
A_DX_42, A_DX_24 = 7, 10  # X-checks (4,2) & (2,4), round-1 detectors -> Z component


def test_probe_a_center_data_qubit() -> None:
    circuit = build(d=3, rounds=2, injector=NullInjector(), probe_erasures=[(3 + 3j, 0)])
    part = partition_dem(circuit)

    assert part.herald_indices.tolist() == [A_HERALD]
    assert len(part.syndrome_indices) == 16  # 17 detectors total, minus the 1 herald
    assert sum(1 for instr in part.dem_pauli if instr.type == "error") == 0

    # Exactly two graphlike edges, both with empty obs mask ((3,3) is NOT on
    # logical_z_support(3), so no component flips L0). The Y hyperedge is absent
    # by Design 2 (its components == the X and Z edges already present).
    assert _edge_set(part, A_HERALD) == {
        ((A_DZ_22, A_DZ_44), ()),  # X on (3,3) -> Z-check pair, no L0
        ((A_DX_42, A_DX_24), ()),  # Z on (3,3) -> X-check pair, no L0
    }


# ---------------------------------------------------------------------------
# Probe B — erasure on the corner data qubit (1,1), between rounds 0 and 1.
# Corner qubit touched by only two checks: each component is a SINGLE-detector
# boundary edge. (1,1) IS in logical_z_support(3), so the X component flips L0.
# ---------------------------------------------------------------------------

# Indices VERIFIED BY HAND 2026-07-09 (docs/dem_worksheet.md, Probe B record).
B_HERALD = 4  # herald detector, coords (1, 1, 0, 1)
B_DX_20 = 5  # top-boundary X-check (2,0), round-1 -> Z component, no L0
B_DZ_22 = 6  # Z-check (2,2), round-1 -> X component, flips L0
B_L0 = 0  # logical observable index


def test_probe_b_corner_data_qubit_boundary_edges() -> None:
    circuit = build(d=3, rounds=2, injector=NullInjector(), probe_erasures=[(1 + 1j, 0)])
    part = partition_dem(circuit)

    assert part.herald_indices.tolist() == [B_HERALD]
    assert len(part.syndrome_indices) == 16
    assert sum(1 for instr in part.dem_pauli if instr.type == "error") == 0

    assert _edge_set(part, B_HERALD) == {
        ((B_DX_20,), ()),  # Z on (1,1) -> single boundary detector, no L0
        ((B_DZ_22,), (B_L0,)),  # X on (1,1) -> single boundary detector, flips L0
    }

    # Boundary edges are genuine 1-tuples: not padded to length 2, not
    # backfilled with a sentinel boundary index.
    for edge in part.herald_table[B_HERALD]:
        assert len(edge.dets) == 1


# ---------------------------------------------------------------------------
# Probe C — erasure on the Z-ancilla (2,2), round 1, between CX layers 2 & 3.
# Mid-round on an ancilla: the X component corrupts the ancilla's own
# measurement across consecutive rounds -> a TIME-LIKE edge; the Z component
# propagates to one data neighbour -> a single space-like detector. No L0.
# ---------------------------------------------------------------------------

# Indices VERIFIED BY HAND 2026-07-09 (docs/dem_worksheet.md, Probe C record).
C_HERALD = 4  # herald detector, coords (2, 2, 1, 1)
C_DX_20 = 5  # X-check (2,0), round-1 -> Z component (data propagation), no L0
C_DZ_22_R1 = 6  # Z-check (2,2), round-1 detector  -> time-like edge endpoint (t=1)
C_DZ_22_R2 = 13  # Z-check (2,2), round-2 detector  -> time-like edge endpoint (t=2)


def test_probe_c_mid_round_ancilla_time_like_edge() -> None:
    circuit = build(
        d=3, rounds=2, injector=NullInjector(), mid_round_probe_erasures=[(2 + 2j, 1, 1)]
    )
    part = partition_dem(circuit)
    coords = circuit.get_detector_coordinates()

    # Herald carries the round-1 timestamp in its coordinate (found via sentinel).
    assert part.herald_indices.tolist() == [C_HERALD]
    assert coords[C_HERALD] == [2.0, 2.0, 1.0, 1.0]
    assert len(part.syndrome_indices) == 16
    assert sum(1 for instr in part.dem_pauli if instr.type == "error") == 0

    assert _edge_set(part, C_HERALD) == {
        ((C_DX_20,), ()),  # Z component -> single space-like detector, no L0
        ((C_DZ_22_R1, C_DZ_22_R2), ()),  # X component -> time-like pair, no L0
    }

    # The time-like edge connects the SAME ancilla across DIFFERENT time slices
    # (coordinate index 2 is the round), not two spatial neighbours (no hook).
    (time_like_edge,) = [
        e for e in part.herald_table[C_HERALD] if e.dets == (C_DZ_22_R1, C_DZ_22_R2)
    ]
    d0, d1 = time_like_edge.dets
    assert coords[d0][:2] == coords[d1][:2]  # same (x, y) plaquette
    assert coords[d0][2] != coords[d1][2]  # different time coordinate


# ---------------------------------------------------------------------------
# Negative tests — malformed mechanisms the real noise model never produces,
# built as synthetic DEMs (heralds declared via the 4th-coordinate sentinel).
# ---------------------------------------------------------------------------


def test_two_herald_detectors_in_one_mechanism_raises() -> None:
    # D0 and D1 are both heralds (4th coord == 1); one mechanism touches both.
    dem = stim.DetectorErrorModel(
        """
        detector(0, 0, 0, 1) D0
        detector(1, 1, 0, 1) D1
        detector(2, 2, 0) D2
        error(0.1) D0 D1 D2
        """
    )
    with pytest.raises(ValueError, match="herald detectors"):
        partition_flattened_dem(dem)


def test_component_with_more_than_two_syndrome_detectors_raises() -> None:
    # A non-herald component (D1 D2 D3) has three syndrome detectors -> not
    # graphlike; D0 is the herald in its own component.
    dem = stim.DetectorErrorModel(
        """
        detector(0, 0, 0, 1) D0
        detector(1, 1, 0) D1
        detector(2, 2, 0) D2
        detector(3, 3, 0) D3
        error(0.1) D1 D2 D3 ^ D0
        """
    )
    with pytest.raises(ValueError, match="not graphlike"):
        partition_flattened_dem(dem)
