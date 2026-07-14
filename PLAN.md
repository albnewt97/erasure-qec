# PLAN.md — Herald-Conditioned Decoding of the Erasure-Converted Rotated Surface Code

> **This document is the single source of truth for the build.** Every implementation
> task references a section number here. Build order is strictly milestone-by-milestone:
> the noiseless circuit and its determinism tests come first; noise, DEM partitioning,
> and decoding only begin once M3 is green.

---

## 0. Project Summary

We implement a distance-d rotated surface code memory experiment (Z basis) in Stim,
under a **biased-erasure noise model**: a fraction `R_e` of the physical error budget
per two-qubit gate is converted into *heralded erasures* (Stim `HERALDED_ERASE`),
the remainder stays as unheralded two-qubit depolarizing noise. A custom
**herald-conditioned matching decoder** (PyMatching, wrapped as a `sinter.Decoder`)
reads the herald bits per shot and sets the corresponding matching-graph edge
weights to zero (conditional error probability 1/2 ⇒ weight ln((1−p)/p) = 0).

**Headline deliverable:** threshold curves p_L vs p_phys for d ∈ {3,5,7,9,11},
one panel per R_e ∈ {0, 0.5, 0.9, 0.98}, showing the threshold sliding from
≈1% toward ≈4–5%, plus an ablation (herald-aware vs. blind decoder on identical shots).

Hardware context: ¹⁷¹Yb neutral-atom Rydberg gates (Wu–Kolkowitz–Puri–Thompson 2022),
dual-rail superconducting cavities.

---

## 1. Repository Architecture

```
erasure-surface-code/
├── pyproject.toml
├── PLAN.md                          # this file
├── README.md
├── LICENSE
├── .gitignore
├── .pre-commit-config.yaml
├── .github/workflows/ci.yml
├── src/erasure_qec/
│   ├── __init__.py
│   ├── config.py                    # frozen dataclasses: NoiseParams, CodeParams, ExperimentConfig (+YAML I/O)
│   ├── circuits/
│   │   ├── layout.py                # pure geometry, NO stim imports
│   │   ├── scheduling.py            # CX orderings as data (correct Z/ᴎ + broken fixture)
│   │   └── builder.py               # only module that emits stim.Circuit
│   ├── noise/
│   │   ├── model.py                 # NoiseParams -> concrete channel specs
│   │   └── injector.py              # NoiseInjector protocol: Null / PauliOnly / BiasedErasure
│   ├── decoding/
│   │   ├── dem_partition.py         # DEM -> (pauli sub-DEM, herald->edges table)
│   │   ├── herald_matching.py       # per-shot reweighted matcher, fast/slow path
│   │   └── sinter_adapter.py        # sinter.Decoder wrappers: herald_mwpm, blind_mwpm
│   └── analysis/
│       ├── statistics.py            # per-round conversion, Wilson/bootstrap CIs, Λ
│       ├── threshold_fit.py         # finite-size scaling collapse fit
│       └── plotting.py              # all figures -> figures/
├── experiments/
│   ├── collect_threshold_sweep.py
│   ├── collect_lambda_scan.py
│   └── configs/{baseline_pauli,erasure_r50,erasure_r98}.yaml
├── tests/
│   ├── conftest.py
│   ├── test_layout.py
│   ├── test_detectors_deterministic.py
│   ├── test_distance_invariants.py
│   ├── test_dem_partition.py
│   ├── test_herald_decoder.py
│   └── test_end_to_end_montecarlo.py   # @pytest.mark.slow
├── data/        # gitignored
└── figures/     # committed
```

Conventions: Python ≥3.11, `uv` for env + lockfile, `ruff` (line length 100),
`mypy --strict`, pytest with `slow` marker. All public functions typed and docstringed.

Dependencies: `stim>=1.14`, `pymatching>=2.2`, `sinter>=1.14`, `numpy`, `scipy`,
`matplotlib`, `pyyaml`; dev: `pytest`, `pytest-cov`, `ruff`, `mypy`, `pre-commit`.

---

## 2. Geometry & Coordinate Specification (layout.py)

Follow Stim's `surface_code:rotated_memory_z` conventions exactly so circuits can be
diffed against the generator.

- **Data qubits:** (2i+1, 2j+1) for 0 ≤ i, j < d. Total d².
- **Ancillas:** even-coordinate plaquette centers, checkerboard X/Z assignment,
  weight-2 boundary stabilizers on alternating edges. Total d²−1.
  Z-type boundary checks sit on the left/right edges (they terminate Z̄ strings);
  X-type on top/bottom.
- `layout.py` exposes:
  - `data_coords(d) -> set[complex]` (use complex numbers x + iy as coordinates)
  - `ancilla_coords(d) -> dict[complex, Basis]` where `Basis ∈ {X, Z}`
  - `plaquette_neighbors(anc, d) -> list[complex | None]` in a FIXED canonical
    slot order [NE, NW, SE, SW]; `None` for off-lattice slots (boundary weight-2 checks)
  - `logical_z_support(d) -> list[complex]` — one horizontal row of data qubits
- **Invariant tests (test_layout.py):** counts (d², d²−1), every data qubit touched by
  ≤4 ancillas, X/Z ancilla adjacency graphs commute (every pair of X- and Z-plaquettes
  shares 0 or 2 data qubits), boundary weight-2 checks appear only on correct edges.

---

## 3. Noiseless Circuit Specification (builder.py) — BUILD THIS FIRST

### 3.1 Round structure (8 TICKs per round)

```
TICK 1 : R   on all ancillas
TICK 2 : H   on X-ancillas
TICK 3 : CX layer 1
TICK 4 : CX layer 2
TICK 5 : CX layer 3
TICK 6 : CX layer 4
TICK 7 : H   on X-ancillas
TICK 8 : MR  on all ancillas
```

X-ancillas are CX **controls** onto data; Z-ancillas are CX **targets** from data.

### 3.2 Hook-safe CX schedule (scheduling.py)

- X-stabilizers touch neighbors in **Z-shaped order:  NE → NW → SE → SW**
- Z-stabilizers touch neighbors in **ᴎ-shaped order: NE → SE → NW → SW**

This orients every mid-window ancilla-fault hook (weight-2 data error)
**perpendicular** to the logical operator it could shorten. `scheduling.py` stores
these as per-basis lists of slot indices into the canonical [NE,NW,SE,SW] order,
plus a deliberately **broken** schedule (both bases using the same order) kept as a
test fixture.

### 3.3 Full memory experiment

1. `QUBIT_COORDS` for every qubit (data + ancilla).
2. `R` on all data qubits (Z basis), then round 0.
3. Rounds 1..T−1 identical, each followed by `SHIFT_COORDS(0, 0, 1)`.
4. Final: `M` all data qubits; closing detectors; `OBSERVABLE_INCLUDE(0)` on
   `logical_z_support(d)`.

### 3.4 Detector algebra (exact spec)

Let s_a(t) = outcome of ancilla a in round t; m_q = final data measurement.

| Where | Detector | Declared for |
|---|---|---|
| Round 0 | D_{a,0} = s_a(0) | **Z-ancillas only** (data starts in |0⟩⊗n, X outcomes random) |
| Rounds 1..T−1 | D_{a,t} = s_a(t) ⊕ s_a(t−1) | all ancillas |
| Time close | D_{a,T} = (⊕_{q∈∂a} m_q) ⊕ s_a(T−1) | **Z-ancillas only** |

Detector coordinates: `DETECTOR(x, y, t)` with t maintained by `SHIFT_COORDS(0,0,1)`.
**Herald detectors (added in M4) use a 4th sentinel coordinate: `DETECTOR(x, y, t, 1)`.**
Syndrome detectors always have 3 coordinates (or 4th = 0). This sentinel is the
partitioning key in §6.

### 3.5 Determinism tests (test_detectors_deterministic.py) — gate for M2

- `builder.build(d, rounds, NullInjector())` compiles for d ∈ {3,5,7}.
- `circuit.compile_detector_sampler().sample(1000)` → **all zeros** (every detector
  deterministic in the noiseless circuit).
- Number of detectors equals the closed-form count:
  `n_det = (d²−1)/2  +  (T−1)(d²−1)  +  (d²−1)/2` for T rounds.
- Diff against `stim.Circuit.generated("surface_code:rotated_memory_z", distance=d,
  rounds=T)` structurally: same detector count, same observable, and
  `shortest_graphlike_error()` lengths match (exact instruction equality NOT required).

---

## 4. Distance & Hook Invariants (M3 gate)

With single-depolarizing noise sprinkled uniformly (any small p, e.g. use
`circuit = builder.build(...).with_inserted_noise(...)` or the PauliOnly injector at
p = 1e−3):

- `len(circuit.shortest_graphlike_error()) == d` for the correct schedule, d ∈ {3,5,7}.
- `circuit.search_for_undetectable_logical_errors(dont_explore_detection_event_sets_with_size_above=4, ...)`
  finds nothing shorter than d.
- The **broken** schedule fixture yields shortest error `< d` (regression test that
  proves the scheduling logic is load-bearing). Record the observed length in the test.

---

## 5. Noise Model Specification (M4)

Parameters: `NoiseParams(p, r_e, p_meas, p_reset, p_idle)`; defaults
p_meas = p_reset = p_idle = p for the uniform sweep.

Per two-qubit gate on (a, b):
1. `DEPOLARIZE2(p · (1 − r_e))` on (a, b)  — residual Pauli component.
2. `HERALDED_ERASE(p · r_e / 2)` on a, and independently on b.
   Stim semantics: with prob q the qubit is replaced by I/2 (uniform Pauli
   {I,X,Y,Z} each q/4) AND a herald bit is appended to the measurement record.

Supporting channels:
- `X_ERROR(p_meas)` immediately before every `MR`/`M`.
- `X_ERROR(p_reset)` immediately after every `R`.
- `DEPOLARIZE1(p_idle)` on qubits idle during a TICK.

Every herald record bit gets `DETECTOR(x, y, t, 1)` (sentinel per §3.4).
The builder must append these herald detectors in the same round they occur.

Key fact: conditional on a herald firing, each nontrivial Pauli on that qubit has
probability 1/2 → matching edge weight ln((1−½)/½) = **0**.

---

## 6. DEM Partition Specification (dem_partition.py, M5)

Input: `circuit.detector_error_model(decompose_errors=True, flatten_loops=True)`.

Algorithm:
1. Read `dem.get_detector_coordinates()`; classify each detector index:
   **herald** if it has a 4th coordinate == 1, else **syndrome**.
2. Walk error instructions. For each `error(p) targets`:
   - If targets contain **no herald detector** → emit unchanged into `dem_pauli`.
   - If targets contain **exactly one herald detector h** → strip h; record the
     remaining (syndrome detectors, logical flips) as a conditioned edge in
     `herald_table[h] : list[ConditionedEdge(dets, obs_mask, cond_p = 0.5-per-component)]`.
   - Targets with ≥2 herald detectors, or a herald detector combined via `^`
     decomposition ambiguity → raise; must not occur for this noise model (test it).
3. Output: `PartitionedDEM(dem_pauli: stim.DetectorErrorModel,
   herald_table: dict[int, list[ConditionedEdge]],
   herald_indices: np.ndarray, syndrome_indices: np.ndarray)`.

---

## 7. ★ HAND-VERIFICATION PLAN FOR THE DEM PARTITION (do on paper BEFORE trusting M5)

This is the only genuinely error-prone component. Verify it by hand on the smallest
non-trivial instance before running any Monte Carlo.

### 7.1 The probe circuit

- d = 3, T = 2 rounds, Z-memory, **NullInjector** (fully noiseless), then insert
  **exactly one** `HERALDED_ERASE(0.25)` by hand on ONE data qubit — choose the
  **center data qubit (3,3)** — placed **between round 0 and round 1** (i.e., after
  round 0's MR TICK, before round 1's reset TICK). Expose this in the builder as a
  debug hook: `builder.build(..., probe_erasures=[(qubit, after_round)])`.
- Expected detector inventory: 4 (round-0 Z) + 8 (round-1 bulk) + 4 (closing Z)
  = 16 syndrome detectors + **1 herald detector** = 17 total. Verify with
  `circuit.num_detectors`.

### 7.2 Paper worksheet — expected DEM lines

Conditional on the herald, the qubit suffers X, Y, or Z each with prob 1/4
(and I with 1/4). Stim will emit herald-linked error mechanisms with total
probabilities q·(component fraction). Work out, on paper, which detectors each
Pauli fires:

| Injected Pauli on (3,3) between rounds | Z-syndrome detectors fired (round-1 comparisons of the two Z-plaquettes adjacent to (3,3)) | X-syndrome detectors fired (round-1 comparisons of the two X-plaquettes adjacent to (3,3)) | Logical Z̄ flip? |
|---|---|---|---|
| Z | none | both adjacent X-plaquette detectors | no |
| X | both adjacent Z-plaquette detectors | none | **yes iff (3,3) ∈ logical_z_support(3)** — check your row choice |
| Y | both Z-plaquette detectors | both X-plaquette detectors | same as X |

Fill in the CONCRETE detector indices by running:
```
print(circuit.detector_error_model(flatten_loops=True))
print(circuit.get_detector_coordinates())
```
and matching coordinates (x, y, t) to the plaquettes adjacent to (3,3) at t = 1.
Write the four expected flattened DEM lines (X-like, Z-like, and the Y line
decomposed into X⊗Z parts by `decompose_errors=True`), each including the herald
detector index D_h. Every line's probability must be a simple function of q = 0.25
(e.g. two-component combinations appear as q/4 + q/4 marginals after Stim merges
symmetric mechanisms — record what Stim actually emits and confirm the algebra:
P(X or Y) = q/2 fires the Z-plaquette pair, P(Z or Y) = q/2 fires the X-plaquette
pair, and the correlation is captured by the decomposition suggestion `^`).

### 7.3 Acceptance criteria (these become test_dem_partition.py)

1. Partition classifies exactly 1 herald detector, 16 syndrome detectors.
2. `dem_pauli` is EMPTY of error instructions (no residual noise in the probe).
3. `herald_table` has exactly one key; its ConditionedEdges cover exactly the
   detector pairs derived in §7.2, no extras, no missing.
4. The observable mask on the X-like edge matches the paper prediction for whether
   (3,3) lies on `logical_z_support(3)`.
5. Repeat the whole worksheet with the probe on a **boundary** data qubit (1,1)
   — now some Pauli components fire only ONE syndrome detector (edge to the spatial
   boundary). Confirm the partition represents boundary edges correctly
   (PyMatching boundary convention: single-detector edges).
6. Repeat with the probe on an **ancilla-adjacent timing**: erasure placed
   mid-round between CX layers 2 and 3 on an ancilla qubit — confirm heralded
   ancilla erasures map to measurement-like edges (time-like detector pairs
   s_a(t)⊕s_a(t±1)) and that no hook-shaped 2-data-qubit edge appears (the
   erasure happens at a known point; the DEM tells you what Stim propagates —
   verify it matches Pauli propagation through the remaining CX layers done
   by hand on the 8-TICK diagram).

Only after 1–6 hold on paper AND in the test does M5 close.

### 7.4 Tooling for the worksheet

- `stim.Circuit.diagram("timeline-svg")` and `"detslice-with-ops-svg"` for the
  8-TICK propagation picture (save into `figures/dev/` while verifying).
- `dem.flattened()` for stable indices; never hand-verify against an unflattened DEM.
- `stim.Circuit.explain_detector_error_model_errors()` to cross-check which circuit
  fault each DEM line came from.

---

## 8. Herald-Conditioned Decoder (herald_matching.py, M6)

- Build base `pymatching.Matching` from `dem_pauli` ONCE at compile time.
- Maintain edge arrays (endpoints, base weights). Per shot:
  - Extract herald bits via `herald_indices`; if **no herald fired** → fast path
    (batched `decode_batch` on the base matcher over all such shots at once).
  - Else slow path: copy base weights, set weight = 0 for every edge in
    `herald_table[h]` for each fired h (adding the edge to the graph if it exists
    only conditionally), decode that single shot.
- Two-tier dispatch is mandatory: at small p·r_e most shots are herald-free.
- Correctness test (test_herald_decoder.py): force specific erasure patterns with
  probability-1 probe erasures; assert herald-aware decoder corrects a pattern the
  blind decoder (heralds stripped, erasure treated as extra depolarizing) fails.

## 9. Sinter Integration & Experiments (M7)

- `sinter_adapter.py`: subclass `sinter.Decoder`; `compile_decoder_for_dem(dem)`
  runs the partition and returns a `CompiledDecoder` with
  `decode_shots_bit_packed`. Register `{"herald_mwpm": ..., "blind_mwpm": ...}`.
- `collect_threshold_sweep.py`: YAML config → list[sinter.Task]
  (d ∈ {3,5,7,9,11}, T = d, p log-spaced in [1e−3, 1e−1], fixed r_e per config)
  → `sinter.collect(custom_decoders=..., save_resume_filepath=data/...)`.

## 10. Analysis (M8)

- Per-round conversion: p_L = ½(1 − (1 − 2·P_L_shot)^{1/T}).
- Threshold: fit p_L near crossing to A + Bx + Cx², x = (p − p_th)·d^{1/ν};
  bootstrap CIs over sinter stats.
- Figures: (i) threshold panels per r_e, (ii) collapse inset, (iii) Λ-factor
  Λ = p_L(d)/p_L(d+2) at fixed p, (iv) herald-aware vs blind ablation,
  (v) hook regression figure (correct vs broken schedule shortest error length).

---

## 11. Milestones & Gates

| M | Deliverable | Gate (must be green before next M) |
|---|---|---|
| M0 | Repo scaffold, pyproject, CI, pre-commit | `uv sync`, ruff, mypy, empty pytest pass |
| M1 | layout.py + scheduling.py + test_layout.py | geometry invariants pass d ∈ {3,5,7} |
| M2 | **Noiseless** builder.py + NullInjector + determinism tests | §3.5 all green |
| M3 | Distance/hook invariants | §4 all green incl. broken-schedule regression |
| M4 | noise/model.py + injector.py + herald detectors | circuit compiles; herald count = expected; determinism still holds at p=0 |
| M5 | dem_partition.py + §7 hand verification | §7.3 criteria 1–6 in tests |
| M6 | herald_matching.py | forced-erasure correctness test |
| M7 | sinter_adapter + experiment scripts | 10⁴-shot smoke run end-to-end |
| M8 | statistics + threshold_fit + plotting | figures reproduce from CSVs deterministically |
| M9 | README with figures, throughput table, CI badge | done |
