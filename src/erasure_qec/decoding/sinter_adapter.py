"""sinter integration: ``herald_mwpm`` and ``blind_mwpm`` (docs/archive/PLAN_v1.md §9).

Both decoders consume IDENTICAL bit-packed shot data — full detector vectors
with the herald columns in place, exactly as sinter samples them from the
circuit. ``herald_mwpm`` routes shots through the M6 two-tier
herald-conditioned matcher; ``blind_mwpm`` zeroes the herald columns
internally and decodes on the static full-DEM matcher (erasure mechanisms
folded in at their unconditional marginals — "erasure treated as extra
depolarizing").

Input contract: the DEM handed to ``compile_decoder_for_dem`` must have been
produced with ``decompose_errors=True, flatten_loops=True,
approximate_disjoint_errors=True`` (the Design-2 partition contract, §6).
sinter does not apply these flags itself for ``HERALDED_ERASE`` circuits, so
every ``sinter.Task`` must carry the DEM explicitly
(``collect_threshold_sweep.py`` does this).

Wu et al. model: ``wu_herald_mwpm`` and ``wu_blind_mwpm`` decode the tagged
DEM of the Wu *decoding* circuit, while sinter samples the Wu *sampling*
circuit (``task.circuit``; sinter 1.16.0 ``_stim_then_decode_sampler.py``).
The two circuits differ by design and share only their detector layout; build
tasks with :func:`wu_task`, which checks that layout. ``wu_herald_mwpm``
partitions the DEM with :func:`partition_tagged_dem`; ``wu_blind_mwpm`` matches
the full tagged DEM with the herald columns zeroed.

Register via ``sinter.collect(..., custom_decoders=CUSTOM_DECODERS)``, or
``WU_DECODERS`` for Wu tasks.
"""

import multiprocessing
from collections.abc import Sequence
from typing import Any

import numpy as np
import numpy.typing as npt
import pymatching
import sinter
import stim

from erasure_qec.circuits.builder import build
from erasure_qec.decoding.dem_partition import partition_flattened_dem
from erasure_qec.decoding.herald_matching import BlindMatchingDecoder, HeraldMatchingDecoder
from erasure_qec.decoding.tagged_partition import partition_tagged_dem
from erasure_qec.noise.wu import WuDecodingInjector, WuNoiseParams, WuSamplingInjector

_U8Array = npt.NDArray[np.uint8]


def _unpack(bit_packed: _U8Array, num_detectors: int) -> npt.NDArray[np.bool_]:
    """Bit-packed shots (n, ceil(dets/8)) -> bool (n, num_detectors)."""
    return np.unpackbits(bit_packed, axis=1, count=num_detectors, bitorder="little").astype(bool)


def _pack(predictions: npt.NDArray[np.bool_]) -> _U8Array:
    """Bool observable predictions -> bit-packed (n, ceil(obs/8)) uint8."""
    return np.packbits(predictions, axis=1, bitorder="little")


class _CompiledHeraldMwpm(sinter.CompiledDecoder):
    """Bit-packed shim around the M6 ``HeraldMatchingDecoder``."""

    def __init__(self, decoder: HeraldMatchingDecoder) -> None:
        self._decoder = decoder

    def decode_shots_bit_packed(self, *, bit_packed_detection_event_data: _U8Array) -> _U8Array:
        shots = _unpack(bit_packed_detection_event_data, self._decoder.num_detectors)
        return _pack(self._decoder.decode_batch(shots))


class _CompiledBlindMwpm(sinter.CompiledDecoder):
    """Static full-DEM matcher; herald columns zeroed before every decode."""

    def __init__(self, dem: stim.DetectorErrorModel) -> None:
        self._num_detectors: int = dem.num_detectors
        self._num_observables: int = dem.num_observables
        self._herald_cols = partition_flattened_dem(dem).herald_indices
        self._matching = pymatching.Matching.from_detector_error_model(dem)

    def decode_shots_bit_packed(self, *, bit_packed_detection_event_data: _U8Array) -> _U8Array:
        shots = _unpack(bit_packed_detection_event_data, self._num_detectors)
        shots[:, self._herald_cols] = False  # strip the herald information
        preds = np.asarray(self._matching.decode_batch(shots.astype(np.uint8)), dtype=np.uint8)
        out = np.zeros((preds.shape[0], self._num_observables), dtype=bool)
        w = min(preds.shape[1], self._num_observables)
        out[:, :w] = preds[:, :w].astype(bool)
        return _pack(out)


class _CompiledBlindFromDecoder(sinter.CompiledDecoder):
    """Bit-packed shim around a ``BlindMatchingDecoder``."""

    def __init__(self, decoder: BlindMatchingDecoder) -> None:
        self._decoder = decoder

    def decode_shots_bit_packed(self, *, bit_packed_detection_event_data: _U8Array) -> _U8Array:
        shots = _unpack(bit_packed_detection_event_data, self._decoder.num_detectors)
        return _pack(self._decoder.decode_batch(shots))


class HeraldMwpmDecoder(sinter.Decoder):
    """§9 ``herald_mwpm``: per-shot herald-reweighted MWPM (M6) under sinter."""

    def compile_decoder_for_dem(self, *, dem: stim.DetectorErrorModel) -> sinter.CompiledDecoder:
        flat = dem.flattened()
        return _CompiledHeraldMwpm(HeraldMatchingDecoder(partition_flattened_dem(flat)))

    def decode_via_files(self, **kwargs: Any) -> None:
        raise NotImplementedError("herald_mwpm only supports in-memory decoding")


class BlindMwpmDecoder(sinter.Decoder):
    """§9 ``blind_mwpm``: herald-stripped static MWPM ablation baseline."""

    def compile_decoder_for_dem(self, *, dem: stim.DetectorErrorModel) -> sinter.CompiledDecoder:
        return _CompiledBlindMwpm(dem.flattened())

    def decode_via_files(self, **kwargs: Any) -> None:
        raise NotImplementedError("blind_mwpm only supports in-memory decoding")


class WuHeraldMwpmDecoder(sinter.Decoder):
    """``wu_herald_mwpm``: herald-reweighted MWPM on the tagged Wu decoding DEM."""

    def compile_decoder_for_dem(self, *, dem: stim.DetectorErrorModel) -> sinter.CompiledDecoder:
        return _CompiledHeraldMwpm(HeraldMatchingDecoder(partition_tagged_dem(dem.flattened())))

    def decode_via_files(self, **kwargs: Any) -> None:
        raise NotImplementedError("wu_herald_mwpm only supports in-memory decoding")


class WuBlindMwpmDecoder(sinter.Decoder):
    """``wu_blind_mwpm``: static MWPM on the full tagged Wu DEM, heralds stripped."""

    def compile_decoder_for_dem(self, *, dem: stim.DetectorErrorModel) -> sinter.CompiledDecoder:
        flat = dem.flattened()
        herald_indices = partition_tagged_dem(flat).herald_indices
        return _CompiledBlindFromDecoder(BlindMatchingDecoder.from_dem(flat, herald_indices))

    def decode_via_files(self, **kwargs: Any) -> None:
        raise NotImplementedError("wu_blind_mwpm only supports in-memory decoding")


CUSTOM_DECODERS: dict[str, sinter.Decoder] = {
    "herald_mwpm": HeraldMwpmDecoder(),
    "blind_mwpm": BlindMwpmDecoder(),
}

# Wu et al. decoders, kept in their own registry: they require the tagged Wu
# decoding DEM (see wu_task) and must not be run on repo-model tasks.
WU_DECODERS: dict[str, sinter.Decoder] = {
    "wu_herald_mwpm": WuHeraldMwpmDecoder(),
    "wu_blind_mwpm": WuBlindMwpmDecoder(),
}


def contract_dem(circuit: stim.Circuit) -> stim.DetectorErrorModel:
    """The §6-contract DEM for ``circuit``, as every ``sinter.Task`` must carry."""
    return circuit.detector_error_model(
        decompose_errors=True,
        flatten_loops=True,
        approximate_disjoint_errors=True,
    )


def wu_task(d: int, params: WuNoiseParams, metadata: dict[str, Any]) -> sinter.Task:
    """One sinter task for the Wu model at distance ``d`` with ``d`` rounds.

    ``circuit`` is the sampling circuit; ``detector_error_model`` is the tagged
    DEM of the decoding circuit (``decompose_errors=True, flatten_loops=True``).
    Raises ``ValueError`` if the two circuits' detector layouts differ.
    """
    sampling = build(d, d, WuSamplingInjector(d, params))
    decoding = build(d, d, WuDecodingInjector(d, params))
    if (
        sampling.num_detectors != decoding.num_detectors
        or sampling.get_detector_coordinates() != decoding.get_detector_coordinates()
    ):
        raise ValueError("Wu sampling and decoding circuits have different detector layouts")
    return sinter.Task(
        circuit=sampling,
        detector_error_model=decoding.detector_error_model(
            decompose_errors=True, flatten_loops=True
        ),
        json_metadata=metadata,
    )


def wu_tasks_parallel(
    specs: Sequence[tuple[int, WuNoiseParams, dict[str, Any]]], processes: int = 8
) -> list[sinter.Task]:
    """:func:`wu_task` for every ``(d, params, metadata)`` spec, built in a process pool.

    Results come back in ``specs`` order. Tasks cross the process boundary by
    pickle, which stim implements via ``str()`` (about 6 significant digits);
    sinter pickles every task to its workers the same way, so workers receive
    identical circuits and DEMs whether a task was built here or serially, and
    ``strong_id`` (computed from ``str()``) is unchanged.

    The pool is closed and joined before returning, so no pool worker is alive
    when the caller starts ``sinter.collect``.
    """
    with multiprocessing.get_context("spawn").Pool(processes) as pool:
        tasks = pool.starmap(wu_task, specs)
        pool.close()
        pool.join()
    return tasks
