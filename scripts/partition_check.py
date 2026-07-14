"""Design-2 DEM partition check on the PLAN.md §7 probe circuits. Not shipped.

Usage: uv run python scripts/partition_check.py [A|B|C]

Builds the selected probe, partitions its DEM, and prints/asserts the
hand-verified ground truth for that probe. Probe C's expectations are filled
in only after its worksheet (guide Part 4) is reconciled.
"""

import sys

from erasure_qec.circuits.builder import build
from erasure_qec.decoding.dem_partition import partition_flattened_dem
from erasure_qec.noise.injector import NullInjector

probe_name = sys.argv[1] if len(sys.argv) > 1 else "A"

if probe_name == "A":
    title = "Probe A: center data qubit (3,3), heralded erasure between rounds 0/1"
    circuit = build(d=3, rounds=2, injector=NullInjector(),
                    probe_erasures=[(3 + 3j, 0)])
    expected_herald = 4
    # Verified 2026-07-08: two bulk 2-detector edges, no L0 ((3,3) not in support).
    expected_edges = [((6, 11), ()), ((7, 10), ())]
elif probe_name == "B":
    title = "Probe B: corner data qubit (1,1), heralded erasure between rounds 0/1"
    circuit = build(d=3, rounds=2, injector=NullInjector(),
                    probe_erasures=[(1 + 1j, 0)])
    expected_herald = 4
    # Prediction (guide Step 3.3): two SINGLE-detector boundary edges;
    # L0 on the {6} edge ((1,1) IS in logical_z_support(3)).
    # NOTE: adjust the L0 representation below to match your ConditionedEdge
    # (e.g. (0,) if obs_mask is a tuple of observable indices).
    expected_edges = [((5,), ()), ((6,), (0,))]
elif probe_name == "C":
    title = "Probe C: Z-ancilla (2,2), round 1, between CX layers 2 and 3"
    circuit = build(d=3, rounds=2, injector=NullInjector(),
                    mid_round_probe_erasures=[(2 + 2j, 1, 1)])
    expected_herald = 4
    expected_edges = [((5,), ()), ((6, 13), ())]  
else:
    raise SystemExit(f"unknown probe {probe_name!r}")

dem = circuit.detector_error_model(
    decompose_errors=True,
    flatten_loops=True,
    approximate_disjoint_errors=True,  # required by stim for HERALDED_ERASE
).flattened()
result = partition_flattened_dem(dem)

herald_indices = list(result.herald_indices)
n_syndrome = len(result.syndrome_indices)
n_pauli_errors = sum(1 for instr in result.dem_pauli if instr.type == "error")

print("=" * 70)
print(title)
print("=" * 70)
print(f"herald_indices             : {herald_indices}")
print(f"# syndrome indices         : {n_syndrome}   (expect 16)")
print(f"# error instrs in dem_pauli: {n_pauli_errors}   (expect 0)")
for h, edges in sorted(result.herald_table.items()):
    print(f"herald_table[{h}]: {len(edges)} edge(s)")
    for dets, obs in sorted((tuple(sorted(e.dets)), e.obs_mask) for e in edges):
        print(f"    edge dets={set(dets)}  obs_mask={obs}")
print("=" * 70)

assert n_syndrome == 16, f"# syndrome {n_syndrome} != 16"
assert n_pauli_errors == 0, f"dem_pauli has {n_pauli_errors} errors, expected 0"
if expected_herald is None:
    print("Probe C: expectations not frozen yet — inspect output manually.")
else:
    assert herald_indices == [expected_herald], f"{herald_indices} != [{expected_herald}]"
    edges = result.herald_table.get(expected_herald, [])
    edge_summary = sorted((tuple(sorted(e.dets)), e.obs_mask) for e in edges)
    assert edge_summary == expected_edges, f"{edge_summary} != {expected_edges}"
    print("ALL CHECKS PASSED")