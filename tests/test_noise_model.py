"""Noise model + herald detector checks (PLAN.md §5, §7.1 probe, M4 gate)."""

import pytest
import stim

from erasure_qec.circuits.builder import build, num_detectors
from erasure_qec.circuits.layout import ancilla_coords, plaquette_neighbors
from erasure_qec.config import NoiseParams
from erasure_qec.noise.injector import ErasureInjector, NullInjector, PauliOnlyInjector
from erasure_qec.noise.model import channel_rates, nonidentity_pauli_probability

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


# --- General analytic herald count for a full ErasureInjector circuit. ---


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
        circuit = build(d, rounds, ErasureInjector(NoiseParams(p=0.01, r_e=0.9)))
        coords = circuit.get_detector_coordinates()
        herald_count = sum(1 for xy in coords.values() if len(xy) == 4 and xy[3] == 1.0)
        assert herald_count == _expected_herald_count(d, rounds)
        assert circuit.num_detectors == num_detectors(d, rounds) + herald_count


# --- Gate: with p=0 the circuit is byte-identical to the noiseless one. ---


def test_erasure_injector_at_p_zero_is_byte_identical_to_noiseless() -> None:
    for d in DISTANCES:
        rounds = d
        reference = build(d, rounds, NullInjector())
        erasure = build(d, rounds, ErasureInjector(NoiseParams(p=0.0, r_e=0.7)))
        assert erasure == reference
        assert str(erasure) == str(reference)


def test_pauli_only_at_p_zero_is_byte_identical_to_noiseless() -> None:
    for d in DISTANCES:
        rounds = d
        reference = build(d, rounds, NullInjector())
        pauli = build(d, rounds, PauliOnlyInjector(NoiseParams(p=0.0)))
        assert pauli == reference
        assert str(pauli) == str(reference)


# --- Edge behavior: r_e = 0 and r_e = 1 special cases. ---


def test_r_e_zero_never_emits_heralded_erase_or_herald_detectors() -> None:
    circuit = build(3, 3, ErasureInjector(NoiseParams(p=0.01, r_e=0.0)))
    assert "HERALDED_ERASE" not in str(circuit)
    coords = circuit.get_detector_coordinates()
    assert all(len(xy) == 3 for xy in coords.values())


def test_r_e_one_never_emits_depolarize2() -> None:
    circuit = build(3, 3, ErasureInjector(NoiseParams(p=0.01, r_e=1.0)))
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
    # Expected herald rate changed from p*r_e/2 to (2/3)*p*r_e: the old rate let
    # the per-gate non-identity budget shrink with r_e because it ignored that
    # an erasure is a non-identity error only 3/4 of the time (see docs/AUDIT.md
    # and channel_rates docstring). q solves 2*q*(3/4) = p*r_e.
    assert rates.herald == 0.02 * 0.9 / (2 * 0.75)
    assert rates.meas_flip == 0.01
    assert rates.reset_flip == 0.03
    assert rates.idle_depolarize == 0.005
    assert rates.idle_herald == 0.0  # convert_idle defaults False


def test_channel_rates_hold_gate_budget_constant() -> None:
    """The per-two-qubit-gate non-identity error probability is ~p for every
    r_e (to first order), so r_e converts errors rather than removing them."""
    p = 0.02
    for r_e in (0.0, 0.5, 0.9, 0.98, 1.0):
        rates = channel_rates(NoiseParams(p=p, r_e=r_e))
        # DEPOLARIZE2 contributes its full probability; each of the two heralds
        # contributes 3/4 of its probability (non-identity fraction).
        budget = rates.depolarize2 + 2 * rates.herald * 0.75
        assert budget == pytest.approx(p)


def test_error_budget_invariance() -> None:
    """Encodes the chosen convention (a): the exact per-2q-gate non-identity
    Pauli probability is p-conserving across r_e -- it equals p minus only the
    O(p^2) inclusion-exclusion overlap for every r_e, NOT the linear p(1 - r_e/4)
    shrink of the old p*r_e/2 rate.

    The correction p^2*(r_e - 3/4 r_e^2) is not monotone: it is 0 at r_e=0,
    peaks at ~p^2/3 near r_e=2/3, and is exactly p^2/4 at r_e=1.
    """
    p = 0.02
    r_es = [0.0, 0.25, 0.5, 2.0 / 3.0, 0.9, 0.98, 1.0]
    probs = [nonidentity_pauli_probability(NoiseParams(p=p, r_e=r)) for r in r_es]

    # Endpoints are exact: p at r_e=0 (no erasure), p - p^2/4 at r_e=1.
    assert probs[0] == pytest.approx(p)
    assert probs[-1] == pytest.approx(p - p * p / 4.0)
    # Every r_e stays within one O(p^2) overlap of p (the whole spread is < p^2/3
    # ~ 1.3e-4 here), i.e. the budget is CONSERVED, not shrunk by ~p*r_e.
    for prob in probs:
        assert p - p * p / 3.0 - 1e-12 <= prob <= p + 1e-12
    assert 0 < max(probs) - min(probs) < p * p / 3.0 + 1e-12

    # The old p*r_e/2 rate would fail this: at r_e=0.98 it gave p*(1-r_e/4) =
    # 0.755p, a 24.5% (=O(p), not O(p^2)) shrink far below the invariance band.
    old_budget = p * (1.0 - 0.98 / 4.0)
    assert old_budget < p - p * p / 3.0


def test_channel_rates_defaults_meas_reset_idle_to_p() -> None:
    params = NoiseParams(p=0.03, r_e=0.5)
    rates = channel_rates(params)
    assert rates.meas_flip == 0.03
    assert rates.reset_flip == 0.03
    assert rates.idle_depolarize == 0.03
    assert rates.idle_herald == 0.0  # convert_idle defaults False


def test_convert_idle_holds_idle_budget_constant() -> None:
    """convert_idle splits the idle DEPOLARIZE1 budget like the 2q gate, but a
    single idle qubit carries one erasure, so idle_herald = (4/3) p_idle r_e and
    the idle budget dep + (3/4) herald stays at p_idle for every r_e."""
    p_idle = 0.01
    for r_e in (0.0, 0.5, 0.9, 1.0):
        rates = channel_rates(
            NoiseParams(p=0.02, r_e=r_e, p_idle=p_idle, convert_idle=True)
        )
        assert rates.idle_depolarize == pytest.approx(p_idle * (1 - r_e))
        assert rates.idle_herald == pytest.approx(p_idle * r_e / 0.75)
        assert rates.idle_depolarize + 0.75 * rates.idle_herald == pytest.approx(p_idle)


def _herald_count(circuit: stim.Circuit) -> int:
    return sum(1 for xy in circuit.get_detector_coordinates().values()
               if len(xy) == 4 and xy[3] == 1.0)


def test_convert_idle_emits_extra_heralds_and_is_noiseless_at_p_zero() -> None:
    gate_only = build(5, 5, ErasureInjector(NoiseParams(p=0.02, r_e=0.98)))
    with_idle = build(
        5, 5, ErasureInjector(NoiseParams(p=0.02, r_e=0.98, convert_idle=True))
    )
    assert _herald_count(with_idle) > _herald_count(gate_only)
    # p = 0 stays byte-identical to the noiseless circuit even with convert_idle.
    reference = build(5, 5, NullInjector())
    zero = build(5, 5, ErasureInjector(NoiseParams(p=0.0, r_e=0.98, convert_idle=True)))
    assert zero == reference


def test_noise_params_rejects_out_of_range_values() -> None:
    with pytest.raises(ValueError):
        NoiseParams(p=1.5)
    with pytest.raises(ValueError):
        NoiseParams(p=0.1, r_e=-0.1)
    with pytest.raises(ValueError):
        NoiseParams(p=0.1, p_meas=2.0)


# --- Compiles and samples correctly at a nonzero noise rate (sanity). ---


def test_erasure_circuit_compiles() -> None:
    circuit = build(5, 5, ErasureInjector(NoiseParams(p=0.005, r_e=0.5)))
    sampler = circuit.compile_detector_sampler()
    dets, obs = sampler.sample(100, separate_observables=True)
    assert dets.shape[1] == circuit.num_detectors
    assert isinstance(circuit, stim.Circuit)
