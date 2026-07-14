"""DEM -> (Pauli sub-DEM, herald -> conditioned-edge table) partition (PLAN.md §6).

Input contract
--------------
The partition consumes ``circuit.detector_error_model(decompose_errors=True,
flatten_loops=True)``. (Stim additionally requires ``approximate_disjoint_errors=True``
whenever the circuit contains ``HERALDED_ERASE`` — a >2-outcome channel — or DEM
construction raises; :func:`partition_dem` therefore always passes it.) The DEM
must be flattened (:meth:`stim.DetectorErrorModel.flattened`) so detector indices
and coordinates are stable and absolute (§7.4 — never hand-verify against an
unflattened DEM).

Design 2
--------
Heralded mechanisms are split at ``^`` separators into graphlike components; each
non-herald component becomes one ConditionedEdge; the herald component is dropped
after identification. Rationale: PyMatching is graphlike, so a multi-detector
hyperedge (e.g. the Y component of an erasure) can never enter the matching graph
directly; setting weight 0 on its constituent component edges covers it
automatically.

Detector classification uses the coordinate sentinel (§3.4): a detector whose 4th
coordinate equals 1 is a herald; everything else is a syndrome detector. Each
error mechanism is routed as:

- **0 heralds** -> emitted unchanged (including its ``^`` decomposition) into
  ``dem_pauli``.
- **exactly 1 herald** -> the herald's (own) component is dropped and every
  remaining component becomes a :class:`ConditionedEdge` under
  ``herald_table[h]`` (deduplicated; see below).
- **>= 2 heralds**, or a herald sharing a component with any syndrome detector or
  logical target -> raise. Neither can occur for this noise model (independent
  per-qubit ``HERALDED_ERASE`` plus independent ``DEPOLARIZE2``).
"""

from dataclasses import dataclass, field

import numpy as np
import numpy.typing as npt
import stim

# 4th DETECTOR coordinate (0-indexed slot 3) marks a herald detector (§3.4).
_SENTINEL_SLOT = 3

# §5 key fact: conditional on a herald firing, every conditioned edge's
# probability is 1/2 -> matching weight ln((1-1/2)/(1/2)) = 0.
HERALD_CONDITIONAL_PROBABILITY = 0.5


@dataclass(frozen=True)
class ConditionedEdge:
    """One herald-conditioned, graphlike matching-graph edge (Design 2).

    ``dets``: the syndrome detectors this edge connects (1 for a boundary edge,
    2 for a bulk edge), sorted. ``obs_mask``: logical observables it flips,
    sorted. ``cond_p``: always :data:`HERALD_CONDITIONAL_PROBABILITY` (§5).
    """

    dets: tuple[int, ...]
    obs_mask: tuple[int, ...]
    cond_p: float = HERALD_CONDITIONAL_PROBABILITY


@dataclass(frozen=True)
class PartitionedDEM:
    """Output of :func:`partition_dem` (§6)."""

    dem_pauli: stim.DetectorErrorModel
    herald_table: dict[int, list[ConditionedEdge]] = field(default_factory=dict)
    herald_indices: npt.NDArray[np.int64] = field(
        default_factory=lambda: np.array([], dtype=np.int64)
    )
    syndrome_indices: npt.NDArray[np.int64] = field(
        default_factory=lambda: np.array([], dtype=np.int64)
    )


def partition_dem(circuit: stim.Circuit) -> PartitionedDEM:
    """Partition ``circuit``'s detector error model per §6 / Design 2.

    Builds the DEM per the module's input contract
    (``decompose_errors=True, flatten_loops=True``, plus the
    ``approximate_disjoint_errors=True`` that stim requires for
    ``HERALDED_ERASE``), flattens it for stable indices, and delegates to
    :func:`partition_flattened_dem`.
    """
    dem = circuit.detector_error_model(
        decompose_errors=True,
        flatten_loops=True,
        approximate_disjoint_errors=True,
    ).flattened()
    return partition_flattened_dem(dem)


def partition_flattened_dem(dem: stim.DetectorErrorModel) -> PartitionedDEM:
    """Partition an already-flattened detector error model per Design 2.

    ``dem`` must already be flattened so detector coordinates/indices are
    absolute. Exposed separately from :func:`partition_dem` so the routing
    logic can be tested against hand-built DEMs, including the deliberately
    malformed inputs that trigger §6's raise conditions but never arise from
    this project's actual noise model.
    """
    herald_set, herald_indices, syndrome_indices = _classify_detectors(dem)

    dem_pauli = stim.DetectorErrorModel()
    herald_table: dict[int, list[ConditionedEdge]] = {}
    # Per herald h: frozenset(dets) -> obs_mask, used both to deduplicate the
    # X/Y/Z mechanisms' shared components and to catch a same-dets/different-mask
    # inconsistency that this noise model must never produce.
    edge_index: dict[int, dict[frozenset[int], tuple[int, ...]]] = {}

    for instr in dem:
        if instr.type != "error":
            dem_pauli.append(instr)
            continue

        components = _split_components(instr.targets_copy())
        total_heralds = sum(
            1
            for comp in components
            for t in comp
            if t.is_relative_detector_id() and t.val in herald_set
        )

        if total_heralds == 0:
            dem_pauli.append(instr)  # keep the mechanism (and its decomposition)
            continue
        if total_heralds >= 2:
            raise ValueError(
                f"error mechanism touches {total_heralds} herald detectors "
                f"(at most 1 is supported by this noise model): {instr}"
            )

        _route_single_herald(
            instr, components, herald_set, herald_table, edge_index
        )

    return PartitionedDEM(
        dem_pauli=dem_pauli,
        herald_table=herald_table,
        herald_indices=herald_indices,
        syndrome_indices=syndrome_indices,
    )


def _classify_detectors(
    dem: stim.DetectorErrorModel,
) -> tuple[set[int], npt.NDArray[np.int64], npt.NDArray[np.int64]]:
    """Split every detector index into herald vs syndrome by the §3.4 sentinel."""
    coords = dem.get_detector_coordinates()
    herald_indices: list[int] = []
    syndrome_indices: list[int] = []
    for i in range(dem.num_detectors):
        xy = coords.get(i, [])
        is_herald = len(xy) > _SENTINEL_SLOT and xy[_SENTINEL_SLOT] == 1.0
        (herald_indices if is_herald else syndrome_indices).append(i)
    return (
        set(herald_indices),
        np.array(herald_indices, dtype=np.int64),
        np.array(syndrome_indices, dtype=np.int64),
    )


def _split_components(targets: list[stim.DemTarget]) -> list[list[stim.DemTarget]]:
    """Split a target list into ``^``-separated components.

    Separators are consumed as delimiters only: they never appear inside a
    component and their ``.val`` is never read (a separator has no integer
    value).
    """
    components: list[list[stim.DemTarget]] = [[]]
    for t in targets:
        if t.is_separator():
            components.append([])
        else:
            components[-1].append(t)
    return components


def _route_single_herald(
    instr: stim.DemInstruction,
    components: list[list[stim.DemTarget]],
    herald_set: set[int],
    herald_table: dict[int, list[ConditionedEdge]],
    edge_index: dict[int, dict[frozenset[int], tuple[int, ...]]],
) -> None:
    """Handle a one-herald mechanism: validate, drop herald component, add edges.

    Two passes because the herald component may appear anywhere in the target
    list (e.g. it is last in ``D6 D11 ^ D4``): first locate/validate the herald
    component, then turn each remaining component into an edge.
    """
    herald_val: int | None = None
    for comp in components:
        herald_dets = [
            t.val for t in comp if t.is_relative_detector_id() and t.val in herald_set
        ]
        if not herald_dets:
            continue
        # The herald must be its own component (D6 D11 ^ D4 is valid;
        # D6 D4 ^ D11 is not) -- otherwise attribution is ambiguous.
        others = [t for t in comp if not (t.is_relative_detector_id() and t.val in herald_set)]
        if others:
            raise ValueError(
                "herald detector shares a component with syndrome/logical "
                f"targets (ambiguous attribution): {instr}"
            )
        herald_val = herald_dets[0]
    assert herald_val is not None  # guaranteed: total_heralds == 1

    for comp in components:
        if any(t.is_relative_detector_id() and t.val in herald_set for t in comp):
            continue  # drop the herald component

        syndrome_dets = [t.val for t in comp if t.is_relative_detector_id()]
        obs = [t.val for t in comp if t.is_logical_observable_id()]

        # Remaining (non-herald) component -> one graphlike ConditionedEdge.
        if len(syndrome_dets) > 2:
            raise ValueError(
                "conditioned component spans more than 2 syndrome detectors "
                f"(not graphlike): {instr}"
            )
        if not syndrome_dets and not obs:
            continue  # empty component (does not occur for this noise model)

        # I-component policy: an instruction whose only content is the herald
        # component reaches here with no remaining components, so it adds
        # nothing -- the identity outcome of an erasure needs no correction; the
        # herald index is still registered in herald_indices.
        dets = tuple(sorted(syndrome_dets))
        obs_mask = tuple(sorted(set(obs)))
        _add_edge(herald_val, dets, obs_mask, instr, herald_table, edge_index)


def _add_edge(
    herald: int,
    dets: tuple[int, ...],
    obs_mask: tuple[int, ...],
    instr: stim.DemInstruction,
    herald_table: dict[int, list[ConditionedEdge]],
    edge_index: dict[int, dict[frozenset[int], tuple[int, ...]]],
) -> None:
    """Insert one edge under ``herald``, deduplicating on (dets, obs_mask)."""
    key = frozenset(dets)
    seen = edge_index.setdefault(herald, {})
    if key in seen:
        if seen[key] != obs_mask:
            raise ValueError(
                f"conflicting observable masks for detectors {sorted(dets)} under "
                f"herald D{herald}: {seen[key]} vs {obs_mask} (from {instr})"
            )
        return  # exact duplicate (e.g. Y re-contributing an X/Z component)
    seen[key] = obs_mask
    herald_table.setdefault(herald, []).append(
        ConditionedEdge(dets=dets, obs_mask=obs_mask, cond_p=HERALD_CONDITIONAL_PROBABILITY)
    )
