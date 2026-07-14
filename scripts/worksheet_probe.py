"""Dev tool for the PLAN.md §7 hand-verification worksheet. Not shipped.

This script writes **scratch** diagrams to ``figures/dev/`` (gitignored). The
canonical, committed copies referenced by ``docs/dem_worksheet.md`` live in
``docs/figures/`` (six per-probe SVGs). Re-running this script does NOT
regenerate those committed versions — it only refreshes the scratch copies;
promote by hand if you deliberately want to update the documented diagrams.
"""

import sys

from erasure_qec.circuits.builder import build
from erasure_qec.circuits.layout import ancilla_coords, logical_z_support
from erasure_qec.noise.injector import NullInjector

# ---- choose the probe on the command line: A, B, or C -------------------
# A: center data qubit (3,3), between rounds 0/1
# B: corner data qubit (1,1), between rounds 0/1
# C: Z-ancilla (4,2), round 1, after 0-indexed CX layer 1 (= between CX layers 2 and 3)
PROBES = {
    "A": dict(probe_erasures=[(3 + 3j, 0)]),
    "B": dict(probe_erasures=[(1 + 1j, 0)]),
    "C": dict(mid_round_probe_erasures=[(2 + 2j, 1, 1)]),
}
probe_name = sys.argv[1] if len(sys.argv) > 1 else "A"

if probe_name == "A":
    c = build(d=3, rounds=2, injector=NullInjector(),
              probe_erasures=[(3 + 3j, 0)])          # center data qubit, after round 0
elif probe_name == "B":
    c = build(d=3, rounds=2, injector=NullInjector(),
              probe_erasures=[(1 + 1j, 0)])          # corner data qubit, after round 0
elif probe_name == "C":
    # after 0-indexed CX layer 1 (= between CX layers 2 and 3)
    c = build(d=3, rounds=2, injector=NullInjector(),
              mid_round_probe_erasures=[(2 + 2j, 1, 1)])
else:
    raise SystemExit(f"unknown probe {probe_name!r}")

print("=" * 70)
print("num_detectors:", c.num_detectors, "(expect 17)")
print("=" * 70)
print("DETECTOR COORDINATES (index -> [x, y, t] or [x, y, t, 1] for herald):")
for idx, coords in sorted(c.get_detector_coordinates().items()):
    tag = "  <-- HERALD" if len(coords) >= 4 and coords[3] == 1 else ""
    print(f"  D{idx}: {coords}{tag}")
print("=" * 70)
print("ANCILLA TYPES:", ancilla_coords(3))
print("LOGICAL Z SUPPORT:", logical_z_support(3))
print("=" * 70)
print("FLATTENED DEM (with decomposition):")
print(
    c.detector_error_model(
        decompose_errors=True,
        flatten_loops=True,
        approximate_disjoint_errors=True,
    )
)
print("=" * 70)
print("EXPLAIN (which circuit fault caused each DEM line):")
dem_filter = c.detector_error_model(
    decompose_errors=True, flatten_loops=True, approximate_disjoint_errors=True
)
for e in c.explain_detector_error_model_errors(dem_filter=dem_filter):
    print(e)
    print("-" * 40)

# Diagrams for tracing propagation on paper (used mainly for Probe C):
with open("figures/dev/probe_timeline.svg", "w") as f:
    f.write(c.diagram("timeline-svg")._repr_svg_())
with open("figures/dev/probe_detslice.svg", "w") as f:
    f.write(c.diagram("detslice-with-ops-svg")._repr_svg_())
print("Saved figures/dev/probe_timeline.svg and probe_detslice.svg")