"""sinter adapter + experiment config units (PLAN.md §9, M7).

Fast checks: the bit-packed sinter path must be exactly equivalent to the M6
decoders on identical shot data, and the three sweep configs must satisfy the
§9 grid requirements.
"""

from pathlib import Path

import numpy as np
import pytest

from erasure_qec.circuits.builder import build
from erasure_qec.config import NoiseParams, load_experiment_config
from erasure_qec.decoding.herald_matching import (
    BlindMatchingDecoder,
    HeraldMatchingDecoder,
)
from erasure_qec.decoding.sinter_adapter import CUSTOM_DECODERS, contract_dem
from erasure_qec.noise.injector import BiasedErasureInjector

CONFIG_DIR = Path(__file__).resolve().parent.parent / "experiments" / "configs"
CONFIG_NAMES = ["baseline_pauli", "erasure_r50", "erasure_r98"]
EXPECTED_R_E = {"baseline_pauli": 0.0, "erasure_r50": 0.5, "erasure_r98": 0.98}


def test_adapters_match_direct_decoders_on_identical_shots() -> None:
    """Bit-packed round trip: herald_mwpm == HeraldMatchingDecoder and
    blind_mwpm == BlindMatchingDecoder, both fed the SAME sampled shots."""
    circuit = build(3, 3, BiasedErasureInjector(NoiseParams(p=0.02, r_e=0.5)))
    dem = contract_dem(circuit)

    n = 512
    packed_dets, _ = circuit.compile_detector_sampler(seed=5).sample(
        n, separate_observables=True, bit_packed=True
    )
    dets = np.unpackbits(
        packed_dets, axis=1, count=circuit.num_detectors, bitorder="little"
    ).astype(bool)

    herald = CUSTOM_DECODERS["herald_mwpm"].compile_decoder_for_dem(dem=dem)
    blind = CUSTOM_DECODERS["blind_mwpm"].compile_decoder_for_dem(dem=dem)
    herald_packed = herald.decode_shots_bit_packed(
        bit_packed_detection_event_data=packed_dets
    )
    blind_packed = blind.decode_shots_bit_packed(
        bit_packed_detection_event_data=packed_dets
    )

    ref_herald = np.packbits(
        HeraldMatchingDecoder.from_circuit(circuit).decode_batch(dets),
        axis=1,
        bitorder="little",
    )
    ref_blind = np.packbits(
        BlindMatchingDecoder(circuit).decode_batch(dets), axis=1, bitorder="little"
    )
    assert (herald_packed == ref_herald).all()
    assert (blind_packed == ref_blind).all()
    # The herald information changes at least one prediction at this noise level.
    assert (herald_packed != blind_packed).any()


def test_registry_names_match_plan() -> None:
    assert set(CUSTOM_DECODERS) == {"herald_mwpm", "blind_mwpm"}


@pytest.mark.parametrize("name", CONFIG_NAMES)
def test_config_satisfies_sweep_grid_requirements(name: str) -> None:
    config = load_experiment_config(CONFIG_DIR / f"{name}.yaml")
    assert config.name == name
    assert config.r_e == EXPECTED_R_E[name]
    assert config.distances == (3, 5, 7, 9, 11)
    # p in [1e-3, 1e-1], at least 10 points, endpoints exact, strictly sorted.
    assert len(config.p_values) >= 10
    assert config.p_values[0] == pytest.approx(1e-3)
    assert config.p_values[-1] == pytest.approx(1e-1)
    assert all(a < b for a, b in zip(config.p_values, config.p_values[1:], strict=False))

    # The 13-point log-spaced base grid must be present; additive extra_p bands
    # (a fit-resolution densification) may add interior points beyond it.
    base = [10.0 ** (-3.0 + i * (2.0 / 12.0)) for i in range(13)]
    for p in base:
        assert any(v == pytest.approx(p) for v in config.p_values)
    extra = [v for v in config.p_values if not any(v == pytest.approx(p) for p in base)]
    # Each config adds a dense fit-resolution band bracketing its crossing:
    # baseline ~1.4%, r50 ~2.6% (measured under the fixed model). The r98 band
    # was extended up to 0.085 when the corrected constant-budget model raised
    # the erasure rate ~33% and shifted the herald crossing higher (see the
    # config comment); the test tracks that config change.
    expected_band = {
        "baseline_pauli": [0.011, 0.013, 0.015, 0.017, 0.019],
        "erasure_r50": [0.018, 0.020, 0.022, 0.024, 0.026],
        "erasure_r98": [0.030, 0.035, 0.040, 0.045, 0.050,
                        0.052, 0.055, 0.058, 0.061, 0.064,
                        0.067, 0.070, 0.073, 0.076, 0.080, 0.085],
    }[name]
    assert extra == pytest.approx(expected_band)

    assert config.max_shots >= 1 and config.max_errors >= 1
    assert set(config.decoders) <= set(CUSTOM_DECODERS)
