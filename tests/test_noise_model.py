"""Noise model + herald detector checks (PLAN.md §5, §7.1 probe, M4 gate)."""

import pytest
import stim

from erasure_qec.circuits.builder import build, num_detectors
from erasure_qec.circuits.layout import ancilla_coords, plaquette_neighbors
from erasure_qec.config import NoiseParams
from erasure_qec.noise.injector import BiasedErasureInjector, NullInjector, PauliOnlyInjector
from erasure_qec.noise.model import channel_rates

DISTANCES = [3, 5, 7]


# --- §7.1 hand-verified probe: the "analytic expectation" for herald count. ---


def test_probe_herald_count_matches_hand_derivation() -> None:
    """d=3, T=2, NullInjector, one HERALDED_ERASE(0.25) on center qubit (3,3)
    after round 0: 4 (round-0 Z) + 8 (round-1 bulk) + 4 (closing Z) = 16
    syndrome detectors + 1 herald detector = 17 total (§7.1).
    """
    circuit = build(3, 2, NullInjector(), probe_erasures=[(complex(3, 3), 0)])
    assert circuit.num_detectors == 17
    assert circuit.num_detectors == num_detectors(3, 2) + 1

    coords = circuit.get_detector_coordinates()
    herald = {i: xy for i, xy in coords.items() if len(xy) == 4 and xy[3] == 1.0}
    syndrome = {i: xy for i, xy in coords.items() if i not in herald}
    assert len(herald) == 1
    assert len(syndrome) == 16
    assert all(len(xy) == 3 for xy in syndrome.values())
    (herald_coord,) = herald.values()
    assert herald_coord == [3.0, 3.0, 0.0, 1.0]


def test_probe_erasure_at_last_round_is_placed_before_closing_shift() -> None:
    """after_round == rounds - 1 still yields exactly one extra herald detector."""
    circuit = build(3, 2, NullInjector(), probe_erasures=[(complex(1, 1), 1)])
    assert circuit.num_detectors == num_detectors(3, 2) + 1
    coords = circuit.get_detector_coordinates()
    herald = {i: xy for i, xy in coords.items() if len(xy) == 4 and xy[3] == 1.0}
    (herald_coord,) = herald.values()
    assert herald_coord == [1.0, 1.0, 1.0, 1.0]


def test_probe_supports_multiple_simultaneous_erasures() -> None:
    circuit = build(
        3, 2, NullInjector(), probe_erasures=[(complex(3, 3), 0), (complex(1, 1), 0)]
    )
    assert circuit.num_detectors == num_detectors(3, 2) + 2
    coords = circuit.get_detector_coordinates()
    herald = {i: xy for i, xy in coords.items() if len(xy) == 4 and xy[3] == 1.0}
    assert len(herald) == 2


# --- General analytic herald count for a full BiasedErasureInjector circuit. ---


def _expected_herald_count(d: int, rounds: int) -> int:
    """Independently recompute the number of heralded qubits per §5: every
    CX-layer edge heralds both its endpoints, across all 4 layers and all
    rounds. Uses layout.py directly (not builder.py internals).
    """
    ancillas = ancilla_coords(d)
    n_edges_per_round = 0
    for a in ancillas:
        n_edges_per_round += sum(1 for nb in plaquette_neighbors(a, d) if nb is not None)
    # Each ancilla-data edge fires exactly once across the 4 CX layers
    # (the schedule assigns each ancilla's non-None neighbors to distinct
    # layers), each firing heralds 2 qubits (the ancilla and the data qubit).
    return n_edges_per_round * 2 * rounds


def test_full_circuit_herald_count_matches_analytic_expectation() -> None:
    for d in DISTANCES:
        rounds = d
        circuit = build(d, rounds, BiasedErasureInjector(NoiseParams(p=0.01, r_e=0.9)))
        coords = circuit.get_detector_coordinates()
        herald_count = sum(1 for xy in coords.values() if len(xy) == 4 and xy[3] == 1.0)
        assert herald_count == _expected_herald_count(d, rounds)
        assert circuit.num_detectors == num_detectors(d, rounds) + herald_count


# --- Gate: with p=0 the circuit is byte-identical to the noiseless one. ---


def test_biased_erasure_at_p_zero_is_byte_identical_to_noiseless() -> None:
    for d in DISTANCES:
        rounds = d
        reference = build(d, rounds, NullInjector())
        biased = build(d, rounds, BiasedErasureInjector(NoiseParams(p=0.0, r_e=0.7)))
        assert biased == reference
        assert str(biased) == str(reference)


def test_pauli_only_at_p_zero_is_byte_identical_to_noiseless() -> None:
    for d in DISTANCES:
        rounds = d
        reference = build(d, rounds, NullInjector())
        pauli = build(d, rounds, PauliOnlyInjector(NoiseParams(p=0.0)))
        assert pauli == reference
        assert str(pauli) == str(reference)


# --- Edge behavior: r_e = 0 and r_e = 1 special cases. ---


def test_r_e_zero_never_emits_heralded_erase_or_herald_detectors() -> None:
    circuit = build(3, 3, BiasedErasureInjector(NoiseParams(p=0.01, r_e=0.0)))
    assert "HERALDED_ERASE" not in str(circuit)
    coords = circuit.get_detector_coordinates()
    assert all(len(xy) == 3 for xy in coords.values())


def test_r_e_one_never_emits_depolarize2() -> None:
    circuit = build(3, 3, BiasedErasureInjector(NoiseParams(p=0.01, r_e=1.0)))
    assert "DEPOLARIZE2" not in str(circuit)
    assert "HERALDED_ERASE" in str(circuit)


def test_pauli_only_injector_never_heralds() -> None:
    circuit = build(3, 3, PauliOnlyInjector(NoiseParams(p=0.05)))
    assert "HERALDED_ERASE" not in str(circuit)
    assert "DEPOLARIZE2" in str(circuit)


# --- noise/model.py: NoiseParams -> ChannelRates arithmetic (§5). ---


def test_channel_rates_formula() -> None:
    params = NoiseParams(p=0.02, r_e=0.9, p_meas=0.01, p_reset=0.03, p_idle=0.005)
    rates = channel_rates(params)
    assert rates.depolarize2 == 0.02 * (1 - 0.9)
    assert rates.herald == 0.02 * 0.9 / 2
    assert rates.meas_flip == 0.01
    assert rates.reset_flip == 0.03
    assert rates.idle_depolarize == 0.005


def test_channel_rates_defaults_meas_reset_idle_to_p() -> None:
    params = NoiseParams(p=0.03, r_e=0.5)
    rates = channel_rates(params)
    assert rates.meas_flip == 0.03
    assert rates.reset_flip == 0.03
    assert rates.idle_depolarize == 0.03


def test_noise_params_rejects_out_of_range_values() -> None:
    with pytest.raises(ValueError):
        NoiseParams(p=1.5)
    with pytest.raises(ValueError):
        NoiseParams(p=0.1, r_e=-0.1)
    with pytest.raises(ValueError):
        NoiseParams(p=0.1, p_meas=2.0)


# --- Compiles and samples correctly at a nonzero noise rate (sanity). ---


def test_biased_circuit_compiles() -> None:
    circuit = build(5, 5, BiasedErasureInjector(NoiseParams(p=0.005, r_e=0.5)))
    sampler = circuit.compile_detector_sampler()
    dets, obs = sampler.sample(100, separate_observables=True)
    assert dets.shape[1] == circuit.num_detectors
    assert isinstance(circuit, stim.Circuit)
