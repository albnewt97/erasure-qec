# DEM hand-verification — record

This is a record of the milestone-M5 hand-verification (PLAN.md §7). Before I trusted
`dem_partition.py`, I predicted on paper exactly what Stim's detector error model (DEM) had to contain
for three small probe circuits, ran Stim, and reconciled every line. Below is what I predicted and why,
what Stim actually emitted, how the two reconciled, and what I took away — the sharpest result being
Probe C, where a Z error copied onto two data qubits fired only one detector, because of the CX
scheduling order. The discipline was to write each prediction before printing Stim's answer, so I
couldn't rationalize backwards from the output.

All three probes reconciled — Probe A on 2026-07-08, Probes B and C on 2026-07-09 — the Design-2
partition was verified against each, and the reconciled indices are now frozen in
`tests/test_dem_partition.py` and `scripts/partition_check.py`. Every probe is d=3, T=2, Z-memory,
noiseless except for a single hand-placed `HERALDED_ERASE`.

## Verified constants for this repository

These came from the live runs and are the facts every prediction below uses.

| Fact | Value |
|---|---|
| Parity convention | interior Z-checks at **(2,2), (4,4)**; interior X-checks at (2,4), (4,2). Boundary: Z-type at (0,4), (6,2); X-type at (2,0), (4,6). (Corrected from my first draft, whose parity was mirrored; all tables here use the verified convention.) |
| `logical_z_support(3)` | **[(1,1), (3,1), (5,1)]** — the y=1 row. So (3,3) ∉ support → Probe A has **no L0 anywhere**; (1,1) ∈ support → Probe B has **L0 on the X and Y lines**. |
| Probe q | 0.25 (builder default `probe_q`) → every DEM line = q/4 = **0.0625** |
| Detector layout (Probes A/B) | Herald is **D4** (emitted after round-0 detectors); D0–D3 = round-0 Z (t=0); D5–D12 = round-1 (t=1); D13–D16 = closing (t=2) |
| Herald time coordinate | Stamped with the emitting round (Probe A/B: `[x, y, 0, 1]`). Heralds are identified only by the 4th coordinate = 1, never by t. |
| Mid-round probe indexing | `mid_round_probe_erasures=(qubit, round, after_cx_layer)`, **0-indexed layers**: value **1** = between CX layers 2 and 3 |
| DEM printing | Printed with `decompose_errors=True` (shows `^` separators) — see the DEM note below |

## Background I relied on

### The rotated surface code

A distance-3 rotated surface code stores one logical qubit in 9 data qubits arranged in a 3×3 grid,
with 8 ancilla qubits interleaved, each responsible for repeatedly measuring one stabilizer — a parity
check on the 2 or 4 data qubits around it. A Z-stabilizer measures Z⊗Z⊗Z⊗Z and detects X errors on
those neighbors (X anticommutes with Z); an X-stabilizer measures X⊗X⊗X⊗X and detects Z errors.

In this repo's coordinate system (matching Stim's generator and `layout.py`), the data qubits sit at
odd-odd coordinates — (1,1), (1,3), (1,5), (3,1), (3,3), (3,5), (5,1), (5,3), (5,5), so d² = 9 — and the
ancillas at even coordinates (plaquette centers and boundary bumps), d²−1 = 8: four interior weight-4
checks at (2,2), (2,4), (4,2), (4,4), and four weight-2 boundary checks.

I verified the checkerboard assignment against `ancilla_coords(3)` on 2026-07-08:

| Ancilla | Type | Weight | Notes |
|---|---|---|---|
| (2,2) | **Z** | 4 | interior; neighbors (1,1), (3,1), (1,3), (3,3) |
| (4,4) | **Z** | 4 | interior; neighbors (3,3), (5,3), (3,5), (5,5) |
| (2,4) | **X** | 4 | interior; neighbors (1,3), (3,3), (1,5), (3,5) |
| (4,2) | **X** | 4 | interior; neighbors (3,1), (5,1), (3,3), (5,3) |
| (0,4) | **Z** | 2 | left-edge boundary check |
| (6,2) | **Z** | 2 | right-edge boundary check |
| (2,0) | **X** | 2 | top-edge boundary check |
| (4,6) | **X** | 2 | bottom-edge boundary check |

The consistency checks held: 4 Z + 4 X = 8; Z-type boundary checks on the left/right edges and X-type
on top/bottom, per PLAN.md §2; and the four round-0 (t=0) detectors from the live run sat exactly on the
four Z-ancillas — so the builder and the layout agree about parity.

### Measurement outcomes, detectors, and the DEM

Each round every ancilla is measured to a bit, which I write s_a(t) for ancilla *a* in round *t*.
Individual outcomes can be random (the code state is a superposition), but certain parities are
deterministic in a noiseless circuit: measuring the same stabilizer twice in a row gives
s_a(t) ⊕ s_a(t−1) = 0; data qubits start in |0⟩, a +1 eigenstate of every Z-stabilizer, so the first
Z-stabilizer reading is s_a(0) = 0 (X-stabilizer round-0 outcomes are random, which is why round 0 has
Z-only detectors); and after the final transversal Z-basis measurement of all data qubits (outcomes
m_q), each Z-stabilizer recomputed offline satisfies (⊕_{q∈∂a} m_q) ⊕ s_a(T−1) = 0. A **detector** is
exactly one of those deterministic parities, declared in the circuit; when noise flips it, it "fires."

For the d=3, T=2 probe the full inventory is:

| Detector group | Formula | Count | Time coord |
|---|---|---|---|
| Round-0 Z checks | D = s_a(0) | 4 | t = 0 |
| Round-1 bulk (all 8 ancillas) | D = s_a(1) ⊕ s_a(0) | 8 | t = 1 |
| Closing Z checks | D = (⊕_{q∈∂a} m_q) ⊕ s_a(1) | 4 | t = 2 |

That is 16 syndrome detectors. The probe adds 1 herald detector (below), so `circuit.num_detectors`
prints 17. I write **D(x, y, t)** for the detector whose declared coordinates are (x, y, t) — e.g.
D(2,4,1) is the round-1 comparison detector of the Z-ancilla at (2,4).

`circuit.detector_error_model()` enumerates every independent noise mechanism as
`error(p) D_i D_j ... [L0]` — with probability p this mechanism fires detectors D_i, D_j, … and (if L0
is present) flips logical observable 0. The DEM is the decoder's entire view of the world; if it is
wrong, everything downstream is garbage, which is why I verified it by hand. Two flags mattered.
`flatten_loops=True` expands repeat blocks so detector indices are plain integers I could match against
`get_detector_coordinates()`. `decompose_errors=True` asks Stim to suggest how to split multi-detector
mechanisms into graph-like pieces, shown with `^` separators; the live run used it and it worked
cleanly — Stim isolated the herald as its own component (`... ^ D4`) and split the Y mechanism into its
two natural pairs (`D6 D11 ^ D7 D10 ^ D4`). I read `A ^ B ^ C` as one physical mechanism with a
suggested split into components A, B, C, and reconciled against the detector *set* of the whole line.
One consequence carried into the partition: a decomposed instruction's `targets_copy()` contains
`stim.target_separator()` entries between components, so the herald-counting logic (§6) has to recognize
separators rather than count them as detectors or read their `.val` — I confirmed that held (Probe A,
below). The undecomposed view is the same lines with the `^`s removed and Y as one flat 5-detector
mechanism.

### What HERALDED_ERASE does

`HERALDED_ERASE(q) target` is atomic: with probability 1−q nothing happens and a 0 is appended to the
measurement record; with probability q the qubit is replaced by the maximally mixed state — a uniformly
random Pauli from {I, X, Y, Z}, each with probability q/4 — and a 1 is appended. That appended bit is
the **herald**; the builder wraps it in a detector with a 4th sentinel coordinate `DETECTOR(x, y, t, 1)`,
which I call **D_h**.

Two consequences shaped every prediction. First, D_h fires for all four Pauli outcomes, including I —
the herald says "an erasure event occurred," not "a nontrivial error occurred" — so there is always a
DEM line `error(q/4) D_h` with no syndrome detectors, the I component. It is easy to forget. Second,
conditional on the herald each of X, Y, Z has probability 1/4, so the probability of "an error that
looks like X to the Z-checks" is P(X)+P(Y) = 1/2. That conditional 1/2 is why the M6 decoder later sets
these edges to weight ln((1−½)/½) = 0. The worksheet is where I watched that 1/2 emerge from the raw
q/4 lines.

### CX propagation rules (needed only for Probe C)

For a CX gate with control **c** and target **t**, conjugating Paulis through: X on the control copies
onto the target (X_c → X_c X_t); Z on the target copies onto the control (Z_t → Z_c Z_t); X on the
target stays put (X_t → X_t) and Z on the control stays put (Z_c → Z_c). X flows in the same direction
as the CX arrow, Z against it; Y = X·Z obeys both rules at once.

### The three probes, and why three

| Probe | Erasure location | What it tests | Status |
|---|---|---|---|
| **A** | Center data qubit (3,3), between rounds 0/1 | Bulk 2-detector edges, closing-detector cancellation, observable mask (no L0 here: (3,3) ∉ support) | ✅ reconciled |
| **B** | Corner data qubit (1,1), between rounds 0/1 | Single-detector **boundary** edges; L0 appears here ((1,1) ∈ support) | ✅ reconciled |
| **C** | **Z-ancilla (2,2)**, round 1, between CX layers 2 and 3 (`after_cx_layer=1`, 0-indexed) | Time-like (measurement-like) edges + propagation through remaining CX layers; confirms no dangerous hook | ✅ reconciled |

Each probe exercises a different code path in `dem_partition.py`. A is the easiest, because there is
zero propagation to trace; C is the involved one.

## The scratch script and the index dictionary

I drove all three probes from a small dev tool, `scripts/worksheet_probe.py`, which builds the selected
probe and dumps the detector count, the coordinates, the DEM, and the
`explain_detector_error_model_errors()` output:

```python
"""Dev tool for the PLAN.md §7 hand-verification worksheet. Not shipped."""
import sys
import stim

from erasure_qec.circuits.builder import build
from erasure_qec.circuits.layout import ancilla_coords, logical_z_support
from erasure_qec.noise.injector import NullInjector

# ---- choose the probe on the command line: A, B, or C -------------------
probe_name = sys.argv[1] if len(sys.argv) > 1 else "A"

if probe_name == "A":
    c = build(d=3, rounds=2, injector=NullInjector(),
              probe_erasures=[(3 + 3j, 0)])          # center data qubit, after round 0
elif probe_name == "B":
    c = build(d=3, rounds=2, injector=NullInjector(),
              probe_erasures=[(1 + 1j, 0)])          # corner data qubit, after round 0
elif probe_name == "C":
    c = build(d=3, rounds=2, injector=NullInjector(),
              mid_round_probe_erasures=[(2 + 2j, 1, 1)])  # Z-ancilla (2,2), round 1, after 0-indexed layer 1
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
print("FLATTENED DEM (no decomposition):")
print(c.detector_error_model(flatten_loops=True))
print("=" * 70)
print("EXPLAIN (which circuit fault caused each DEM line):")
for e in c.explain_detector_error_model_errors():
    print(e)
    print("-" * 40)

# Diagrams for tracing propagation on paper (used mainly for Probe C):
with open("figures/dev/probe_timeline.svg", "w") as f:
    f.write(str(c.diagram("timeline-svg")))
with open("figures/dev/probe_detslice.svg", "w") as f:
    f.write(str(c.diagram("detslice-with-ops-svg")))
print("Saved figures/dev/probe_timeline.svg and probe_detslice.svg")
```

The script writes scratch diagrams to `figures/dev/` (gitignored); the canonical copies embedded below
live in `docs/figures/`.

The parity-check run confirmed the constants table: `num_detectors` was 17; exactly one detector carried
the herald tag, **D4** at coords `[3, 3, 0, 1]`, identified by the sentinel coordinate; there were four
detectors at t = 0 (D0–D3), eight at t = 1 (D5–D12), four at t = 2 (D13–D16), with the t=0 detectors on
the four Z-ancillas (2,2), (6,2), (0,4), (4,4); and `logical_z_support(3)` = [(1,1), (3,1), (5,1)], with
(3,3) not in it (Probe A has no L0) and (1,1) in it (Probe B carries L0 on the X and Y lines).

Probes A and B share the same layout — both place the erasure between rounds 0/1, so the herald lands at
index 4 either way — while Probe C's herald index differs, and I re-derived it from Probe C's own dump.
The index dictionary, filled from the live Probe-A coordinate dump (2026-07-08):

| Coordinate name | Meaning (this repo's parity) | Index |
|---|---|---|
| D_h (Probes A/B) | herald | **D4** |
| D(2,0,1) | X-boundary check (2,0), round 1 | **D5** |
| D(2,2,1) | **Z**-check (2,2), round 1 | **D6** |
| D(4,2,1) | **X**-check (4,2), round 1 | **D7** |
| D(6,2,1) | Z-boundary check (6,2), round 1 | **D8** |
| D(0,4,1) | Z-boundary check (0,4), round 1 | **D9** |
| D(2,4,1) | **X**-check (2,4), round 1 | **D10** |
| D(4,4,1) | **Z**-check (4,4), round 1 | **D11** |
| D(4,6,1) | X-boundary check (4,6), round 1 | **D12** |
| D(2,2,2) | Z-check (2,2), closing | **D13** |
| D(6,2,2) | Z-boundary (6,2), closing | **D14** |
| D(0,4,2) | Z-boundary (0,4), closing | **D15** |
| D(4,4,2) | Z-check (4,4), closing | **D16** |
| D0–D3 | round-0 Z detectors: (2,2), (6,2), (0,4), (4,4) | D0, D1, D2, D3 |

## Probe A — erasure on the center data qubit (3,3)

![Probe A timeline](figures/probe_center_3_3_timeline.svg)
*Timeline: the single `HERALDED_ERASE(0.25)` sits in the idle gap after round 0's `MR`, before round 1's reset — no CX gates between injection and measurement, so there is zero propagation to trace.*

![Probe A detslice](figures/probe_center_3_3_detslice.svg)
*Detector slices: (3,3) is interior, so each non-`I` Pauli lights a **pair** of round-1 detectors; the closing (t=2) detectors stay dark because the two flips cancel inside each Z-check parity.*

The erasure sits in the idle gap after round 0's MR and before round 1's ancilla resets. Round 0
measured a clean state, so the four t=0 detectors can never fire; whatever Pauli the erasure injects is
sitting on (3,3) when round 1 starts, so round 1's stabilizer measurements see it in full; and there
are no gates between the injection and round 1, so there is zero propagation to trace — the injected
Pauli is exactly the Pauli the stabilizers see.

(3,3) is the center qubit, the only one touched by all four interior checks: the Z-checks (2,2) and
(4,4) [detectors D6, D11], sensitive to X on (3,3); and the X-checks (2,4) and (4,2) [detectors D10, D7],
sensitive to Z on (3,3).

Conditional on the herald, the injected Pauli is I, X, Y, or Z, each with probability q/4 = 0.0625.
Working each case:

- **I** (0.0625): the erasure fired (herald bit = 1) but the replacement Pauli was identity; no
  stabilizer notices. Fires {D4} and nothing else, no observable flip.
- **X** (0.0625): X on (3,3) anticommutes with the two Z-stabilizers containing it, so ancillas (2,2)
  and (4,4) get flipped round-1 outcomes — **D6 and D11** fire. Then the subtle part, the closing
  detectors. The X persists to the end and flips the final measurement m(3,3). In the closing detector
  of Z-check (2,2),

      D13 = m(1,1) ⊕ m(3,1) ⊕ m(1,3) ⊕ m(3,3) ⊕ s_{(2,2)}(1)

  the error flips m(3,3) and it already flipped s_{(2,2)}(1) (that is why D6 fired). Two flips inside one
  parity cancel, so D13 does not fire, and identically D16 (the (4,4) closing detector). This
  cancellation is the whole content of the §3.4 detector algebra: a persistent data error between round
  t and the end produces a syndrome at round t and nothing afterward. X on (3,3) would flip logical Z̄
  only if (3,3) were on `logical_z_support(3)`; the support is the y=1 row, and (3,3) is not on it.
  Fires {D4, D6, D11}, no L0.
- **Z** (0.0625): Z on (3,3) anticommutes with the two X-stabilizers — **D10 and D7** fire. A Z error is
  invisible to Z-basis measurements: it does not affect the final data measurements m_q, and there are
  no X-type closing detectors in a Z-memory experiment; it also cannot flip Z̄ (Z̄ is made of Z's, and Z
  commutes with Z). Fires {D4, D7, D10}, no L0.
- **Y** = X·Z (0.0625): detector firing is linear, so Y fires the union of the X-case and Z-case sets;
  the observable is as in the X case (the Z part of Y cannot touch Z̄). Fires {D4, D6, D7, D10, D11},
  no L0.

So the predicted DEM, order aside (Stim sorts its own way; no L0 anywhere since (3,3) ∉ support):

```
error(0.0625) D4                    # I  component
error(0.0625) D4 D6 D11             # X  component (Z-check pair)
error(0.0625) D4 D7 D10             # Z  component (X-check pair)
error(0.0625) D4 D6 D7 D10 D11      # Y  component (all four)
```

— these four lines and nothing else, since NullInjector adds no other noise. Two checks on the
prediction itself: all four detector sets are distinct, so Stim has nothing to merge and every
probability should be exactly 0.0625; and the marginals P(Z-check pair fires) = P(X)+P(Y) = 0.125 = q/2
and P(X-check pair fires) = P(Z)+P(Y) = 0.125 = q/2 give conditional probability 1/2 each given the
herald — the 1/2 that becomes the weight-0 edge in M6.

**Reconciliation record — Probe A, run 2026-07-08. Passed.** Stim emitted (with `decompose_errors=True`,
so `^` marks suggested components):

| Stim line | Detector set | Matched prediction | L0 |
|---|---|---|---|
| `error(0.0625) D4` | {D4} | I component ✅ | none ✅ |
| `error(0.0625) D6 D11 ^ D4` | {D4, D6, D11} | X component ✅ | none ✅ |
| `error(0.0625) D7 D10 ^ D4` | {D4, D7, D10} | Z component ✅ | none ✅ |
| `error(0.0625) D6 D11 ^ D7 D10 ^ D4` | {D4, D6, D7, D10, D11} | Y component ✅ | none ✅ |

A perfect bijection: four lines, all at exactly 0.0625, every line containing D4, no t=2 detectors
(closing-detector cancellation confirmed), no L0, and no extra mechanisms. Stim's decomposition isolated
the herald as its own `^` component and split the Y line into the two natural pairs — the structure the
partition consumes.

Running the Design-2 partition on this DEM gave `herald_indices` == [4]; the 16 syndrome indices as
{0,…,16} \ {4}; `dem_pauli` with zero error instructions; and `herald_table[4]` holding exactly the two
per-component edges {6,11} and {7,10}, the Y line re-contributing the same two (deduped). I chose
per-component edges (Design 2) over whole-mechanism hyperedges because PyMatching is graphlike: the
4-detector Y hyperedge could never enter the matching graph, and the Y case is covered automatically
once both component edges go to weight 0. The I-component line (herald only, no syndromes) contributes
no edge. Every edge's obs mask is empty here, since (3,3) ∉ support; the L0-carrying path is Probe B's
job. On the separator question I had flagged: `partition_check.py` and the tests build the DEM with
`decompose_errors=True, flatten_loops=True`; the target walk treats `stim.target_separator()` entries as
component boundaries and never as detectors; and the "exactly one herald per mechanism" check counts
heralds across the whole line, then drops them — the herald arrives as its own `^` component, so the
check does not require it to share a component with syndrome detectors.

## Probe B — erasure on the corner data qubit (1,1)

![Probe B timeline](figures/probe_boundary_1_1_timeline.svg)
*Timeline: same idle-gap placement as Probe A, but on the corner qubit (1,1), which only two checks touch.*

![Probe B detslice](figures/probe_boundary_1_1_detslice.svg)
*Detector slices: because (1,1) has only two neighboring checks, each Pauli component fires a **single** detector — a boundary edge to the virtual node — and the X component additionally flips `L0`, since (1,1) is on the logical support.*

(1,1) is a corner qubit touched by only two checks instead of four, so some Pauli components fire only
one syndrome detector. In matching-graph language those are boundary edges — an edge from a detector
node to the virtual boundary node — and Probe B is what shows the partition represents single-detector
conditioned edges faithfully. The classic bug there is code that assumes every edge has two endpoints
and either crashes, pads, or silently drops these.

In this repo's layout, (1,1) is touched by exactly two checks: the interior **Z**-check at (2,2)
[round-1 detector **D6**], sensitive to X on (1,1); and the weight-2 **X** boundary check at (2,0) (top
edge) [round-1 detector **D5**], sensitive to Z on (1,1). Both memberships are in
`plaquette_neighbors(2+2j, 3)` and `plaquette_neighbors(2+0j, 3)`. And since `logical_z_support(3)` =
[(1,1), (3,1), (5,1)], (1,1) is on the support, so the X and Y components carry L0 — the observable-mask
path that Probe A could not exercise.

The placement is the same as Probe A (between rounds 0/1, zero propagation), q = 0.25, and the same
detector-index layout (herald = D4, since the probe is emitted at the same point in the build):

| Case (prob 0.0625 each) | Fires | L0? |
|---|---|---|
| I | {D4} | no |
| X | {D4, **D6**} — **one** syndrome detector | **YES** |
| Z | {D4, **D5**} — **one** syndrome detector | no |
| Y | {D4, D5, D6} | **YES** |

Closing-detector cancellation works exactly as in Probe A: the X error flips m(1,1) and s_{(2,2)}(1),
which cancel inside D13, so again there are no t=2 detectors anywhere. (The X-boundary check (2,0) has no
closing detector at all; closing detectors exist only for Z-checks.) Four lines, all at 0.0625, all
containing D4, L0 on exactly the X and Y lines, nothing else in the DEM. In the partition these become
the two single-detector boundary edges; in M6 they turn into PyMatching boundary edges (weight 0 when
the herald fires), and the test frozen here is what guarantees M6 receives them intact.

**Reconciliation record — Probe B, run 2026-07-09. Passed.** Herald: D4 at coords `[1, 1, 0, 1]`. Stim
emitted:

| Stim line | Detector set | Matched prediction | L0 |
|---|---|---|---|
| `error(0.0625) D4` | {D4} | I ✅ | none ✅ |
| `error(0.0625) D5 ^ D4` | {D4, D5} | Z (single boundary detector) ✅ | none ✅ |
| `error(0.0625) D6 L0 ^ D4` | {D4, D6} | X (single boundary detector) ✅ | **L0 ✅** |
| `error(0.0625) D5 ^ D6 L0 ^ D4` | {D4, D5, D6} | Y ✅ | L0 ✅ |

The key finding: Stim scopes the L0 target inside the component it belongs to (`D6 L0`), which confirmed
that Design-2's per-component obs-mask extraction is the right convention. The partition output
(2026-07-09) was `herald_table[4]` == two edges, `dets={5} obs=()` and `dets={6} obs=(0,)` —
single-detector boundary edges intact, the observable flag on exactly the physics-correct edge, and the
Y-line components deduped cleanly against identical `(dets, obs_mask)` keys. §7.3 criterion 5: closed.

## Probe C — erasure on Z-ancilla (2,2), mid-round

![Probe C timeline](figures/probe_mid_round_ancilla_2_2_timeline.svg)
*Timeline: the erasure is injected **mid-round on ancilla (2,2), between CX layers 2 and 3** (`mid_round_probe_erasures=[(2+2j, 1, 1)]`). At that instant (2,2) has already read its NE/SE neighbors (3,1),(3,3) but not yet its NW/SW neighbors (1,1),(1,3) — that split is exactly what makes the Z-case result traceable by eye.*

![Probe C detslice](figures/probe_mid_round_ancilla_2_2_detslice.svg)
*Detector slices: the Z component copies onto **two** data qubits (1,1) and (1,3), yet fires **only one** detector — (1,1)'s X-check reads it in a later layer (round-1 detector fires), while (1,3)'s X-check already read it before injection (never seen, T=2 has no round 2). Two data errors, one detector: the concrete proof that hook analysis is CX-schedule-dependent.*

I retargeted this probe after the parity check: (4,2) turned out to be an X-ancilla in this repo, so the
probe moved to (2,2), an interior Z-ancilla. The script line is `mid_round_probe_erasures=[(2 + 2j, 1, 1)]`
— round 1, after 0-indexed CX layer 1, i.e. between the 2nd and 3rd CX layers.

Probes A and B injected errors in an idle gap, with no gates between injection and measurement. Probe C
injects inside round 1, between CX layers 2 and 3, on the Z-ancilla at (2,2), so the injected Pauli
propagates through CX layers 3 and 4 before the ancilla is measured. I traced that propagation by hand
(rules above). The probe checks two things: that heralded ancilla errors produce time-like edges (the
measurement-error lookalike, D(a,1) and D(a,2) firing together), and that the residual propagation onto
data qubits produces exactly the edges Pauli algebra predicts — in particular no hook aligned with the
logical operator's matching graph.

"Between layers 2 and 3" is defined by `scheduling.py`, not by counting TICKs. The Z-ancilla at (2,2) is
a CX target; its data neighbors are the controls. The canonical slot order is [NE, NW, SE, SW], and the
Z-schedule visits slots in the ᴎ order NE → SE → NW → SW (PLAN.md §3.2). Looking up which coordinate
sits in each slot from `plaquette_neighbors(2+2j, 3)`:

| CX layer (0-idx) | Slot visited (Z-schedule) | Data qubit (control) | Done before injection? |
|---|---|---|---|
| 0 | NE | **(3,1)** | ✔ done |
| 1 | SE | **(3,3)** | ✔ done |
| 2 | NW | **(1,1)** | ✘ still to come |
| 3 | SW | **(1,3)** | ✘ still to come |

Filled from `plaquette_neighbors(2+2j, 3)` = [(3,1), (1,1), (3,3), (1,3)] (canonical [NE, NW, SE, SW]
order), run 2026-07-09: q_NW = (1,1) and q_SW = (1,3). The committed timeline
[`figures/probe_mid_round_ancilla_2_2_timeline.svg`](figures/probe_mid_round_ancilla_2_2_timeline.svg)
(embedded above) shows the erasure instruction sitting after the layer-2 CX on (2,2) and before layer 3.
Probe C's herald is not D4 — the mid-round probe is emitted inside round 1, so the herald lands between
the round-0 and round-1 syndrome detectors at a different position, found by the `<-- HERALD` tag — with
17 detectors again.

The ancilla was reset to |0⟩ at the start of the round and is MR-measured in the Z basis at TICK 8. The
remaining gates touching it are CX(q_NW → anc) and CX(q_SW → anc), with the ancilla the target of both.

**Case X on the ancilla** (0.0625). X on a CX target does not propagate to the control, so the X just
sits on the ancilla through layers 3 and 4 and flips the Z-basis MR outcome: s_{(2,2)}(1) flips. A
flipped s(1) appears in exactly two detectors — D(2,2,1) = s(1) ⊕ s(0) fires, and
D(2,2,2) = (⊕ m_q) ⊕ s(1) fires (the data qubits are untouched, so nothing cancels this time, in
contrast to Probe A). No data qubit was harmed, so no other detectors and no observable flip:
{D_h, D(2,2,1), D(2,2,2)}, a pure time-like edge, indistinguishable from a measurement flip on that
ancilla — the measurement-like edge §7.3 criterion 6 asks for.

**Case Z on the ancilla** (0.0625). Z on a CX target copies onto the control, so the Z copies onto q_NW
(layer 3) and onto q_SW (layer 4), and the ancilla ends the round carrying a Z too. But Z on the ancilla
does not flip its Z-basis MR (Z|0⟩ = |0⟩, Z|1⟩ = −|1⟩ — a phase, invisible to measurement), so the
ancilla's own detectors do not fire. What remains is a Z⊗Z error on the two data qubits q_NW, q_SW,
created mid-round-1.

A Z on a data qubit fires the round-1-vs-round-0 comparisons of the X-checks containing it — but the
timing matters, because the Z lands on q_NW during round 1, after some of round 1's CX layers already
ran. For each X-check containing the qubit, the question is whether that X-check had already executed
its CX with the qubit before the moment of injection. If the CX came in layers 3–4 (after injection),
the Z is seen already in round 1 and D(X-check, 1) fires; if it came in layers 1–2 (before injection),
round 1 misses it and the Z is first seen in round 2 — but this circuit has no round 2, T=2 means the
next thing is the final data measurement, Z errors are invisible to Z-basis data measurements, and there
are no X-type closing detectors, so that X-check never sees it at all.

This is the genuinely fiddly part. For each of q_NW and q_SW, listing the X-checks containing it and
which layer (from the X-schedule, Z-order NE→NW→SE→SW) reads it:

| Data qubit | X-checks containing it | Layer each reads it | Detectors fired |
|---|---|---|---|
| q_NW = **(1,1)** | boundary (2,0) only | SW slot → layer **3** (after Z lands at layer 2) → **seen** | **D5** |
| q_SW = **(1,3)** | interior (2,4) only | NW slot → layer **1** (before Z lands at layer 3) → **never seen** | none |

Filled and verified 2026-07-09. The takeaway that made this probe worth doing: the Z landed on two data
qubits but fired only one detector, purely because of CX layer ordering — the concrete demonstration
that hook analysis is schedule-dependent. The final Z-case prediction is {D_h} ∪ (the fired X-check
detectors). A few structural facts held in the table: everything fired is an X-type round-1 detector, so
the edge lives entirely in the X-syndrome half of the matching graph; a Z⊗Z data error cannot flip
logical Z̄ (Z commutes with Z), so no L0; and therefore this "hook" is harmless to the Z-memory
experiment — the dangerous hook the schedule was designed against would contribute a short chain in the
Z̄-relevant (Z-check) graph, and the DEM shows no such edge from this mechanism, which is §7.3
criterion 6's "no hook-shaped 2-data-qubit edge in the dangerous graph," made concrete.

**Case Y on the ancilla** (0.0625). Y = X·Z, so the union of both pictures: the time-like pair
{D(2,2,1), D(2,2,2)} and the X-check detectors from the Z case. No L0. **Case I** (0.0625): {D_h} alone,
as always.

For a probe with propagation, `explain_detector_error_model_errors()` is the arbiter: it names the exact
circuit instruction behind each DEM line, so any disagreement between my hand propagation and Stim
localizes to either the Pauli algebra (recheck the rules against the detslice) or the probe's placement
(recheck the schedule table — a wrong layer placement is the usual culprit, and the prediction is then
correct for a different circuit than the one built).

**Reconciliation record — Probe C, run 2026-07-09. DEM passed.** Herald: **D4** at coords `[2, 2, 1, 1]`
(time-stamped round 1, found by the sentinel tag). Stim emitted:

| Stim line | Detector set | Matched prediction | L0 |
|---|---|---|---|
| `error(0.0625) D4` | {D4} | I ✅ | none ✅ |
| `error(0.0625) D6 D13 ^ D4` | {D4, D6, D13} | X → time-like pair D(2,2,1)+D(2,2,2) ✅ | none ✅ |
| `error(0.0625) D5 ^ D4` | {D4, D5} | Z → only (2,0)'s round-1 detector; the (1,3) copy never seen ✅ | none ✅ |
| `error(0.0625) D6 D13 ^ D5 ^ D4` | {D4, D5, D6, D13} | Y = union ✅ | none ✅ |

The §7.3 criterion-6 structural claims held: a pure time-like edge {6,13} exists; there is no L0
anywhere — even though (1,1) is in the logical support, it received a Z, which commutes with Z̄; and no
Z̄-graph hook edge exists, since the data-error component fires only X-type detectors. The frozen
partition expectation for `partition_check.py C` is `expected_herald = 4`,
`expected_edges = [((5,), ()), ((6, 13), ())]`.

## Freezing the results into tests

I froze the reconciled, hand-verified indices into `tests/test_dem_partition.py` as literal constants
with their coordinate meanings — a case where magic numbers are the point, since their whole value is
that they were checked against the physics:

```python
# Indices VERIFIED BY HAND on 2026-07-08: see docs/dem_worksheet.md, Probe A.
# Parity in this repo: (2,2)/(4,4) are Z-checks; (2,4)/(4,2) are X-checks.
D_H = 4                       # herald, coords (3, 3, 0, 1)
D_Z22_R1, D_Z44_R1 = 6, 11    # Z-checks (2,2), (4,4), round-1 detectors
D_X42_R1, D_X24_R1 = 7, 10    # X-checks (4,2), (2,4), round-1 detectors
D_XB20_R1 = 5                 # X boundary check (2,0), round-1 (Probe B)
```

The Probe-A test asserts `herald_indices` == [D_H], 16 syndrome indices, `dem_pauli.num_errors == 0`,
the herald-table edges, and an empty obs mask on every edge (since (3,3) ∉ support). The Probe-B test
asserts herald == [4], the two single-detector boundary edges {D_XB20_R1} and {D_Z22_R1} plus their
union, and the obs mask set on exactly the {6} and {5,6} edges ((1,1) ∈ support). The Probe-C test
asserts the time-like pair, the data-Z edge, and the absence of L0 on all of them. A negative test per
§6 constructs a DEM line with two herald detectors and asserts the partition raises.

## Quick reference

Pauli ↔ detector rules (between-rounds data error, Z-memory): X on data q fires the round-(t+1)
detectors of the Z-checks containing q, closing detectors cancel, L0 iff q ∈ logical support; Z on data
q fires the round-(t+1) detectors of the X-checks containing q, with no closing effect and never L0; Y
is the union.

CX propagation: X flows control→target, Z flows target→control, Y does both.

Ancilla mid-round error (Z-ancilla = CX target): X on the ancilla gives a time-like edge
{D(a,t), D(a,t+1)}; Z on the ancilla copies onto every not-yet-visited control, and the ancilla's own
measurement is unaffected.

Herald facts: D_h fires for I too, so there is always an `error(q/4) D_h` line; conditional on the
herald, each nontrivial Pauli-equivalence class has probability 1/2, which becomes a weight-0 edge.

Stim flags: `flatten_loops=True` always; `decompose_errors=True` shows the `^` component splits, and the
partition handles the decomposed form (recognizing separators, per the DEM note above);
`explain_detector_error_model_errors()` is the arbiter of any disagreement.

Verified constants (2026-07-08): Z-checks (2,2)/(4,4); X-checks (2,4)/(4,2); support = y=1 row; Probe
A/B herald = D4; D5=(2,0,1), D6=(2,2,1), D7=(4,2,1), D10=(2,4,1), D11=(4,4,1), D13=(2,2,2).
