"""Herald-conditioned MWPM decoder (PLAN.md §8).

Compiled ONCE per circuit: the herald-free sub-DEM (``dem_pauli`` from the §6
partition) becomes the base ``pymatching.Matching``; its edge array
(endpoints, base weights, fault ids) is extracted, and every herald's
ConditionedEdges are resolved onto it. Edges with no unheralded counterpart
in ``dem_pauli`` are kept separately as *conditional-only* edges and enter
the graph only on shots where their herald fired.

Per shot, herald bits are extracted via ``herald_indices``. Two-tier
dispatch (mandatory per §8):

- **Fast path** — no herald fired: all such shots are decoded in a single
  batched ``decode_batch`` call on the base matcher. At small ``p * r_e``
  most shots take this path.
- **Slow path** — >= 1 herald fired: shots are grouped by their herald
  *signature* (the set of fired heralds); ONE reweighted matcher is built
  per distinct signature (M7 throughput requirement) and each group is
  decoded batched. Reweighting copies the base weight vector and sets
  weight 0 on every edge conditioned on a fired herald (§5 key fact:
  conditional error probability 1/2 => weight ln((1-p)/p) = 0). The
  reweighted matcher is built vectorized from a precompiled sparse check
  matrix, not via per-edge Python calls.

Conditional-only edges (no unheralded counterpart in ``dem_pauli``) are kept
as permanent check-matrix columns at a large "absent" weight and dropped to 0
while their herald fires. This is behaviorally identical to inserting them on
demand: their detectors can only fire when their herald fired (the erasure's
Pauli is applied only alongside its herald bit), so an unfired conditional
edge is never on any matched path.

``BlindMatchingDecoder`` is the ablation baseline (§8's "blind" decoder):
the same static graph built from the FULL DEM — erasure mechanisms folded
in at their unconditional marginal probabilities, i.e. erasure treated as
extra depolarizing — with herald bits stripped before decoding.
"""

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
import pymatching
import stim
from scipy.sparse import csc_matrix

from erasure_qec.decoding.dem_partition import (
    PartitionedDEM,
    partition_dem,
    partition_flattened_dem,
)

_BoolArray = npt.NDArray[np.bool_]

# (min_det, max_det) for a bulk edge; (det, None) for a boundary edge.
_EdgeKey = tuple[int, int | None]

# Weight given to a conditional-only edge while its herald is NOT fired.
# Must exceed any real min-weight path (max plausible: ~50 edges x weight
# ~14 at p=1e-6) so the edge is never chosen, while staying small enough
# that pymatching's integer weight discretization keeps resolving genuine
# weight differences.
_ABSENT_WEIGHT = 1000.0


@dataclass(frozen=True)
class _Edge:
    """One matching-graph edge: endpoints, base weight, logical fault ids."""

    u: int
    v: int | None  # None -> edge to the virtual boundary node
    weight: float
    fault_ids: frozenset[int]


def _edge_key(dets: Sequence[int]) -> _EdgeKey:
    """Canonical dictionary key for an edge's detector endpoints."""
    if len(dets) == 1:
        return (int(dets[0]), None)
    if len(dets) == 2:
        a, b = sorted(int(d) for d in dets)
        return (a, b)
    raise ValueError(f"matching edges span 1 or 2 detectors, got {len(dets)}")


def _validate(detection_events: _BoolArray, num_detectors: int) -> _BoolArray:
    dets = np.asarray(detection_events, dtype=bool)
    if dets.ndim != 2 or dets.shape[1] != num_detectors:
        raise ValueError(
            f"expected detection events of shape (n_shots, {num_detectors}), "
            f"got {dets.shape}"
        )
    return dets


def _pad_block(preds: npt.NDArray[np.uint8], num_observables: int) -> _BoolArray:
    """Widen a pymatching prediction block to exactly ``num_observables`` columns."""
    preds2 = np.atleast_2d(np.asarray(preds))
    out = np.zeros((preds2.shape[0], num_observables), dtype=bool)
    w = min(preds2.shape[1], num_observables)
    out[:, :w] = preds2[:, :w].astype(bool)
    return out


class HeraldMatchingDecoder:
    """Per-shot herald-reweighted MWPM decoder (§8)."""

    def __init__(self, partition: PartitionedDEM) -> None:
        dem = partition.dem_pauli
        self.num_detectors: int = dem.num_detectors
        self.num_observables: int = dem.num_observables
        self._herald_cols: npt.NDArray[np.int64] = partition.herald_indices

        # Base matcher: compiled once from the herald-free sub-DEM. An
        # error-free dem_pauli (probe circuits) yields a valid matcher with
        # zero edges; the fast path then never needs to decode (see below).
        self._base: pymatching.Matching = pymatching.Matching.from_detector_error_model(dem)

        # Canonical edge arrays (endpoints, base weights, fault ids).
        self._base_edges: dict[_EdgeKey, _Edge] = {}
        for u, v, attrs in self._base.edges():
            key = _edge_key([u] if v is None else [u, v])
            self._base_edges[key] = _Edge(
                u=key[0],
                v=key[1],
                weight=float(attrs["weight"]),
                fault_ids=frozenset(int(f) for f in attrs["fault_ids"]),
            )

        # Conditioned edges with no unheralded counterpart in dem_pauli become
        # extra columns at _ABSENT_WEIGHT (see module docstring). Their fault
        # ids come from the ConditionedEdge's obs_mask; for edges that DO
        # resolve onto a base edge, the base fault ids are authoritative
        # (zero-weighting never changes which logical the edge flips).
        columns: list[_Edge] = list(self._base_edges.values())
        col_of_key: dict[_EdgeKey, int] = {
            key: i for i, key in enumerate(self._base_edges)
        }
        conditioned_cols: dict[int, set[int]] = {}
        for h, cond_edges in partition.herald_table.items():
            cols: set[int] = set()
            for ce in cond_edges:
                key = _edge_key(ce.dets)
                if key not in col_of_key:
                    col_of_key[key] = len(columns)
                    columns.append(
                        _Edge(
                            u=key[0],
                            v=key[1],
                            weight=_ABSENT_WEIGHT,
                            fault_ids=frozenset(ce.obs_mask),
                        )
                    )
                cols.add(col_of_key[key])
            conditioned_cols[int(h)] = cols

        # Precompiled arrays for the slow path: sparse check matrix (detector
        # rows x edge columns), sparse faults matrix (observable rows), and the
        # base weight vector. Per signature, only the weight vector is copied.
        self._col_weights: npt.NDArray[np.float64] = np.array(
            [e.weight for e in columns], dtype=np.float64
        )
        self._cols_to_zero: dict[int, npt.NDArray[np.int64]] = {
            h: np.array(sorted(cols), dtype=np.int64)
            for h, cols in conditioned_cols.items()
        }
        n_cols = len(columns)
        self._n_cols = n_cols
        check_rows: list[int] = []
        check_cols: list[int] = []
        fault_rows: list[int] = []
        fault_cols: list[int] = []
        for j, edge in enumerate(columns):
            check_rows.append(edge.u)
            check_cols.append(j)
            if edge.v is not None:
                check_rows.append(edge.v)
                check_cols.append(j)
            for f in edge.fault_ids:
                fault_rows.append(f)
                fault_cols.append(j)
        self._check_matrix = csc_matrix(
            (np.ones(len(check_rows), dtype=np.uint8), (check_rows, check_cols)),
            shape=(self.num_detectors, n_cols),
        )
        self._faults_matrix = csc_matrix(
            (np.ones(len(fault_rows), dtype=np.uint8), (fault_rows, fault_cols)),
            shape=(self.num_observables, n_cols),
        )

    @classmethod
    def from_circuit(cls, circuit: stim.Circuit) -> "HeraldMatchingDecoder":
        """Partition ``circuit``'s DEM (§6) and compile the decoder from it."""
        return cls(partition_dem(circuit))

    def decode_batch(self, detection_events: _BoolArray) -> _BoolArray:
        """Decode full detector vectors (herald bits included) to observable flips.

        Args:
            detection_events: bool array of shape ``(n_shots, num_detectors)``,
                exactly as sampled from the circuit (herald detector columns in
                place).

        Returns:
            bool array of shape ``(n_shots, num_observables)`` of predicted
            logical observable flips.
        """
        dets = _validate(detection_events, self.num_detectors)
        n = dets.shape[0]
        syndromes = dets.copy()
        syndromes[:, self._herald_cols] = False
        heralds = dets[:, self._herald_cols]

        out = np.zeros((n, self.num_observables), dtype=bool)
        quiet = ~heralds.any(axis=1)

        # Fast path: one batched decode over every herald-free shot.
        if bool(quiet.any()):
            if self._base_edges:
                preds = self._base.decode_batch(syndromes[quiet].astype(np.uint8))
                out[np.flatnonzero(quiet)] = _pad_block(preds, self.num_observables)
            elif syndromes[quiet].any():
                raise RuntimeError(
                    "herald-free shot has detection events but dem_pauli has no "
                    "error mechanisms to explain them"
                )

        # Slow path with herald-signature grouping (M7): one reweighted
        # matcher per DISTINCT set of fired heralds, each group decoded in a
        # single batched call.
        slow_rows = np.flatnonzero(~quiet)
        if slow_rows.size:
            self._decode_slow_grouped(slow_rows, heralds, syndromes, out)
        return out

    def _decode_slow_grouped(
        self,
        slow_rows: npt.NDArray[np.int64],
        heralds: _BoolArray,
        syndromes: _BoolArray,
        out: _BoolArray,
    ) -> None:
        """Decode heralded shots, grouped by their fired-herald signature."""
        if self._n_cols == 0:
            if syndromes[slow_rows].any():
                raise RuntimeError(
                    "heralded shot has detection events but no edges exist for them"
                )
            return

        groups: dict[bytes, list[int]] = {}
        for i in slow_rows:
            groups.setdefault(heralds[i].tobytes(), []).append(int(i))

        for signature, rows in groups.items():
            fired_mask = np.frombuffer(signature, dtype=bool)
            weights = self._col_weights.copy()
            for h in self._herald_cols[fired_mask]:
                cols = self._cols_to_zero.get(int(h))
                if cols is not None:
                    weights[cols] = 0.0
            matcher = pymatching.Matching.from_check_matrix(
                self._check_matrix,
                weights=weights,
                faults_matrix=self._faults_matrix,
                use_virtual_boundary_node=True,
            )
            preds = matcher.decode_batch(syndromes[rows].astype(np.uint8))
            out[rows] = _pad_block(preds, self.num_observables)


class BlindMatchingDecoder:
    """Static MWPM baseline: heralds stripped, erasure as extra depolarizing (§8).

    Built from the FULL DEM, so every erasure mechanism contributes its
    unconditional marginal probability to the static edge weights. Herald
    detector columns are zeroed before decoding — the herald *information*
    is discarded, which is the entire point of the ablation.
    """

    def __init__(self, circuit: stim.Circuit) -> None:
        dem = circuit.detector_error_model(
            decompose_errors=True,
            flatten_loops=True,
            approximate_disjoint_errors=True,
        ).flattened()
        self.num_detectors: int = dem.num_detectors
        self.num_observables: int = dem.num_observables
        self._herald_cols: npt.NDArray[np.int64] = partition_flattened_dem(dem).herald_indices
        self._matching: pymatching.Matching = pymatching.Matching.from_detector_error_model(dem)

    def decode_batch(self, detection_events: _BoolArray) -> _BoolArray:
        """Decode full detector vectors, ignoring the herald bits entirely."""
        dets = _validate(detection_events, self.num_detectors)
        syndromes = dets.copy()
        syndromes[:, self._herald_cols] = False
        preds = self._matching.decode_batch(syndromes.astype(np.uint8))
        return _pad_block(preds, self.num_observables)
