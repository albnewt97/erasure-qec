"""Tagged decoding DEM -> :class:`PartitionedDEM` for the Wu et al. model.

Input contract
--------------
:func:`partition_tagged_dem` consumes the flattened DEM of a decoding circuit
built with :class:`~erasure_qec.noise.wu.WuDecodingInjector`, i.e.
``decoding.detector_error_model(decompose_errors=True,
flatten_loops=True).flattened()``. In that circuit the flag qubit is
noiseless, so no error mechanism touches a herald detector. Instead every
erasure generator carries the instruction tag ``wu_e{k}:{Xa|Za|Xb|Zb}``, and
stim propagates the tag into the DEM without merging mechanisms whose tags
differ (stim 1.16.0).

Routing
-------
Herald detectors are classified by the 4th-coordinate sentinel
(:func:`~erasure_qec.decoding.dem_partition._classify_detectors`); the ``k``-th
herald detector in ascending index order belongs to location ``k``.

- **untagged** error mechanism -> emitted unchanged into ``dem_pauli``.
- **tagged** ``wu_e{k}`` mechanism -> split at ``^`` separators; each component
  (1 or 2 syndrome detectors, no herald detector) becomes a
  :class:`~erasure_qec.decoding.dem_partition.ConditionedEdge` under herald
  ``k``, deduplicated with the same observable-mask conflict check as
  :func:`~erasure_qec.decoding.dem_partition.partition_flattened_dem`. Tagged
  mechanisms are NOT added to ``dem_pauli``.

The result is consumed unchanged by
:class:`~erasure_qec.decoding.herald_matching.HeraldMatchingDecoder`.
"""

import re

import stim

from erasure_qec.decoding.dem_partition import (
    ConditionedEdge,
    PartitionedDEM,
    _add_edge,
    _classify_detectors,
    _split_components,
)

_TAG = re.compile(r"^wu_e(\d+):(Xa|Za|Xb|Zb)$")


def partition_tagged(sampling: stim.Circuit, decoding: stim.Circuit) -> PartitionedDEM:
    """Cross-check the two circuits' detector layouts, then partition ``decoding``'s DEM.

    Raises ``ValueError`` if the circuits differ in detector count or detector
    coordinates, i.e. if the sampling circuit's detection events cannot be
    decoded against the decoding circuit's DEM.
    """
    if sampling.num_detectors != decoding.num_detectors:
        raise ValueError(
            f"detector counts differ: sampling {sampling.num_detectors}, "
            f"decoding {decoding.num_detectors}"
        )
    if sampling.get_detector_coordinates() != decoding.get_detector_coordinates():
        raise ValueError("sampling and decoding circuits have different detector coordinates")
    dem = decoding.detector_error_model(decompose_errors=True, flatten_loops=True).flattened()
    return partition_tagged_dem(dem)


def partition_tagged_dem(dem: stim.DetectorErrorModel) -> PartitionedDEM:
    """Partition a flattened, tagged Wu decoding DEM (see module docstring).

    Raises ``ValueError`` on: an untagged mechanism touching a herald; a tag
    not of the form ``wu_e{k}:{Xa|Za|Xb|Zb}`` or with ``k`` >= the number of
    herald detectors; a tagged component touching a herald, spanning more than
    2 detectors, or flipping an observable with no detector; conflicting
    observable masks for the same edge; or a number of distinct tagged
    locations different from the number of herald detectors.
    """
    herald_set, herald_indices, syndrome_indices = _classify_detectors(dem)

    dem_pauli = stim.DetectorErrorModel()
    herald_table: dict[int, list[ConditionedEdge]] = {}
    edge_index: dict[int, dict[frozenset[int], tuple[int, ...]]] = {}
    locations_seen: set[int] = set()

    for instr in dem:
        if instr.type != "error":
            dem_pauli.append(instr)
            continue

        targets = instr.targets_copy()
        if not instr.tag:
            if any(t.is_relative_detector_id() and t.val in herald_set for t in targets):
                raise ValueError(f"untagged error mechanism touches a herald detector: {instr}")
            dem_pauli.append(instr)
            continue

        match = _TAG.match(instr.tag)
        if match is None:
            raise ValueError(f"unrecognised error-mechanism tag {instr.tag!r}: {instr}")
        k = int(match.group(1))
        if k >= len(herald_indices):
            raise ValueError(
                f"tag location {k} out of range for {len(herald_indices)} herald detectors"
            )
        locations_seen.add(k)
        herald = int(herald_indices[k])

        for comp in _split_components(targets):
            dets = [t.val for t in comp if t.is_relative_detector_id()]
            obs = [t.val for t in comp if t.is_logical_observable_id()]
            if any(d in herald_set for d in dets):
                raise ValueError(f"tagged component touches a herald detector: {instr}")
            if len(dets) > 2:
                raise ValueError(
                    f"tagged component spans more than 2 detectors (not graphlike): {instr}"
                )
            if not dets:
                if obs:
                    raise ValueError(f"tagged component flips an observable undetected: {instr}")
                continue
            _add_edge(
                herald,
                tuple(sorted(dets)),
                tuple(sorted(set(obs))),
                instr,
                herald_table,
                edge_index,
            )

    if len(locations_seen) != len(herald_indices):
        raise ValueError(
            f"{len(locations_seen)} tagged erasure locations but "
            f"{len(herald_indices)} herald detectors"
        )
    # Declare every observable of the input DEM in dem_pauli. Without this, a
    # DEM whose only observable-flipping mechanisms are tagged (R_e = 1,
    # p_m = 0) yields a dem_pauli with num_observables == 0, and the decoder
    # would size its observable outputs from that.
    for i in range(dem_pauli.num_observables, dem.num_observables):
        dem_pauli.append("logical_observable", [], [stim.target_logical_observable_id(i)])
    if (dem_pauli.num_detectors, dem_pauli.num_observables) != (
        dem.num_detectors,
        dem.num_observables,
    ):
        raise ValueError("dem_pauli lost detectors or observables of the input DEM")

    return PartitionedDEM(
        dem_pauli=dem_pauli,
        herald_table=herald_table,
        herald_indices=herald_indices,
        syndrome_indices=syndrome_indices,
    )
