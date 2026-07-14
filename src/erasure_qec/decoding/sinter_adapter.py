"""sinter integration: ``herald_mwpm`` and ``blind_mwpm`` (PLAN.md §9).

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

Register via ``sinter.collect(..., custom_decoders=CUSTOM_DECODERS)``.
"""

from typing import Any

import numpy as np
import numpy.typing as npt
import pymatching
import sinter
import stim

from erasure_qec.decoding.dem_partition import partition_flattened_dem
from erasure_qec.decoding.herald_matching import HeraldMatchingDecoder

_U8Array = npt.NDArray[np.uint8]


def _unpack(bit_packed: _U8Array, num_detectors: int) -> npt.NDArray[np.bool_]:
    """Bit-packed shots (n, ceil(dets/8)) -> bool (n, num_detectors)."""
    return np.unpackbits(
        bit_packed, axis=1, count=num_detectors, bitorder="little"
    ).astype(bool)


def _pack(predictions: npt.NDArray[np.bool_]) -> _U8Array:
    """Bool observable predictions -> bit-packed (n, ceil(obs/8)) uint8."""
    return np.packbits(predictions, axis=1, bitorder="little")


class _CompiledHeraldMwpm(sinter.CompiledDecoder):
    """Bit-packed shim around the M6 ``HeraldMatchingDecoder``."""

    def __init__(self, decoder: HeraldMatchingDecoder) -> None:
        self._decoder = decoder

    def decode_shots_bit_packed(
        self, *, bit_packed_detection_event_data: _U8Array
    ) -> _U8Array:
        shots = _unpack(bit_packed_detection_event_data, self._decoder.num_detectors)
        return _pack(self._decoder.decode_batch(shots))


class _CompiledBlindMwpm(sinter.CompiledDecoder):
    """Static full-DEM matcher; herald columns zeroed before every decode."""

    def __init__(self, dem: stim.DetectorErrorModel) -> None:
        self._num_detectors: int = dem.num_detectors
        self._num_observables: int = dem.num_observables
        self._herald_cols = partition_flattened_dem(dem).herald_indices
        self._matching = pymatching.Matching.from_detector_error_model(dem)

    def decode_shots_bit_packed(
        self, *, bit_packed_detection_event_data: _U8Array
    ) -> _U8Array:
        shots = _unpack(bit_packed_detection_event_data, self._num_detectors)
        shots[:, self._herald_cols] = False  # strip the herald information
        preds = np.asarray(
            self._matching.decode_batch(shots.astype(np.uint8)), dtype=np.uint8
        )
        out = np.zeros((preds.shape[0], self._num_observables), dtype=bool)
        w = min(preds.shape[1], self._num_observables)
        out[:, :w] = preds[:, :w].astype(bool)
        return _pack(out)


class HeraldMwpmDecoder(sinter.Decoder):
    """§9 ``herald_mwpm``: per-shot herald-reweighted MWPM (M6) under sinter."""

    def compile_decoder_for_dem(
        self, *, dem: stim.DetectorErrorModel
    ) -> sinter.CompiledDecoder:
        flat = dem.flattened()
        return _CompiledHeraldMwpm(HeraldMatchingDecoder(partition_flattened_dem(flat)))

    def decode_via_files(self, **kwargs: Any) -> None:
        raise NotImplementedError("herald_mwpm only supports in-memory decoding")


class BlindMwpmDecoder(sinter.Decoder):
    """§9 ``blind_mwpm``: herald-stripped static MWPM ablation baseline."""

    def compile_decoder_for_dem(
        self, *, dem: stim.DetectorErrorModel
    ) -> sinter.CompiledDecoder:
        return _CompiledBlindMwpm(dem.flattened())

    def decode_via_files(self, **kwargs: Any) -> None:
        raise NotImplementedError("blind_mwpm only supports in-memory decoding")


CUSTOM_DECODERS: dict[str, sinter.Decoder] = {
    "herald_mwpm": HeraldMwpmDecoder(),
    "blind_mwpm": BlindMwpmDecoder(),
}


def contract_dem(circuit: stim.Circuit) -> stim.DetectorErrorModel:
    """The §6-contract DEM for ``circuit``, as every ``sinter.Task`` must carry."""
    return circuit.detector_error_model(
        decompose_errors=True,
        flatten_loops=True,
        approximate_disjoint_errors=True,
    )
