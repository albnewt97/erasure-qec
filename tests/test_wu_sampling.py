"""Wu et al. (2022) sampling injector: branch statistics and circuit structure."""

import itertools
import math

import numpy as np
import numpy.typing as npt
import pytest
import stim

from erasure_qec.circuits.builder import build
from erasure_qec.circuits.layout import num_code_qubits
from erasure_qec.config import NoiseParams
from erasure_qec.noise.injector import ErasureInjector
from erasure_qec.noise.wu import WuNoiseParams, WuSamplingInjector, append_wu_gate_chain

SHOTS = 200_000
SEED = 20220114
N_SIGMA = 5.0

# Toy circuits (not the surface code): qubits a=0, b=1, flag=2. Each entry is
# (name, basis of a, basis of b). A qubit prepared and measured in Z detects X
# flips; one prepared and measured in X detects Z flips.
TOY_CIRCUITS = [("x_a,z_b", "Z", "X"), ("z_a,x_b", "X", "Z")]


def _toy_circuit(basis_a: str, basis_b: str, p_e: float, p_p: float) -> stim.Circuit:
    """Prepare a, b in their bases, apply one Wu chain, then record (h, bit_a, bit_b)."""
    c = stim.Circuit()
    c.append("R", [0, 1])
    for q, basis in ((0, basis_a), (1, basis_b)):
        if basis == "X":
            c.append("H", [q])
    append_wu_gate_chain(c, 0, 1, 2, p_e, p_p)
    for q, basis in ((0, basis_a), (1, basis_b)):
        c.append("M" if basis == "Z" else "MX", [q])
    return c


def _flips(pauli: str, basis: str) -> bool:
    """Whether a single-qubit Pauli flips a measurement in ``basis``."""
    return pauli not in ("I", basis)


def _flip_fraction(basis_a: str, basis_b: str, *, include_identity: bool) -> float:
    """Exact fraction of two-qubit Paulis flipping at least one measured bit."""
    paulis = [
        (pa, pb)
        for pa, pb in itertools.product("IXYZ", "IXYZ")
        if include_identity or (pa, pb) != ("I", "I")
    ]
    hits = sum(_flips(pa, basis_a) or _flips(pb, basis_b) for pa, pb in paulis)
    return hits / len(paulis)


def _sample(circuit: stim.Circuit) -> npt.NDArray[np.bool_]:
    return np.asarray(circuit.compile_sampler(seed=SEED).sample(SHOTS), dtype=np.bool_)


def _assert_rate(hits: npt.NDArray[np.bool_], expected: float, what: str) -> None:
    n = hits.size
    sigma = math.sqrt(expected * (1.0 - expected) / n)
    observed = float(hits.mean())
    assert abs(observed - expected) <= N_SIGMA * sigma, (
        f"{what}: observed {observed:.5f}, expected {expected:.5f} +/- {N_SIGMA} * {sigma:.5f}"
    )


# --- (a) Branch statistics on toy circuits. ---


@pytest.mark.parametrize(("name", "basis_a", "basis_b"), TOY_CIRCUITS)
def test_exclusive_chain_branch_statistics(name: str, basis_a: str, basis_b: str) -> None:
    """p_e = 0.3, p_p = 0.2: herald rate, erasure randomization, unheralded Pauli rate."""
    p_e, p_p = 0.3, 0.2
    rec = _sample(_toy_circuit(basis_a, basis_b, p_e, p_p))
    h, bit_a, bit_b = rec[:, 0], rec[:, 1], rec[:, 2]
    bit_a_name, bit_b_name = name.split(",")

    _assert_rate(h, p_e, "P(h)")
    _assert_rate(bit_a[h], 0.5, f"P({bit_a_name} | h)")
    _assert_rate(bit_b[h], 0.5, f"P({bit_b_name} | h)")
    _assert_rate(bit_a[h] & bit_b[h], 0.25, f"P({name} | h)")

    erasure_flip = _flip_fraction(basis_a, basis_b, include_identity=True)
    assert erasure_flip == 12 / 16
    _assert_rate((bit_a | bit_b) & h, p_e * erasure_flip, "P(any bit and h)")

    k = _flip_fraction(basis_a, basis_b, include_identity=False)
    assert k == 12 / 15
    not_h = ~h
    _assert_rate(
        ~bit_a[not_h] & ~bit_b[not_h], 1.0 - p_p * k / (1.0 - p_e), "P(all bits 0 | not h)"
    )


@pytest.mark.parametrize(("name", "basis_a", "basis_b"), TOY_CIRCUITS)
def test_erasure_only_chain_never_flips_without_herald(
    name: str, basis_a: str, basis_b: str
) -> None:
    """p_e = 0.4, p_p = 0: herald rate, and no unheralded flip at all."""
    p_e = 0.4
    rec = _sample(_toy_circuit(basis_a, basis_b, p_e, 0.0))
    h, bit_a, bit_b = rec[:, 0], rec[:, 1], rec[:, 2]
    _assert_rate(h, p_e, "P(h)")
    _assert_rate(bit_a[h], 0.5, "P(bit_a | h)")
    _assert_rate(bit_b[h], 0.5, "P(bit_b | h)")
    _assert_rate(bit_a[h] & bit_b[h], 0.25, f"P({name} | h)")
    assert not (bit_a[~h] | bit_b[~h]).any()


# --- (b) R_e = 0 reproduces the existing model's circuit text. ---


@pytest.mark.parametrize("d", [3, 5])
def test_re_zero_matches_erasure_injector_text(d: int) -> None:
    p = 0.01
    wu = build(d, d, WuSamplingInjector(d, WuNoiseParams(p, 0.0)))
    existing = build(
        d, d, ErasureInjector(NoiseParams(p, 0.0, p_meas=0.0, p_reset=0.0, p_idle=0.0))
    )
    assert str(wu) == str(existing)


# --- (c) One herald detector per CX pair. ---


@pytest.mark.parametrize("d", [3, 5])
def test_herald_detector_count_equals_cx_pairs(d: int) -> None:
    injector = WuSamplingInjector(d, WuNoiseParams(0.01, 0.5))
    circuit = build(d, d, injector)

    n_pairs = sum(
        len(inst.targets_copy()) // 2
        for inst in circuit.flattened()
        if isinstance(inst, stim.CircuitInstruction) and inst.name == "CX"
    )
    coords = circuit.get_detector_coordinates()
    four_coord = [xy for xy in coords.values() if len(xy) == 4]

    assert n_pairs > 0
    assert len(four_coord) == n_pairs
    assert all(xy[3] == 1.0 for xy in four_coord)
    assert len(injector.locations) == n_pairs
    assert [k for _, _, k in injector.locations] == list(range(n_pairs))


# --- (d) The surface-code sampling circuit samples. ---


@pytest.mark.parametrize("r_e", [0.5, 0.98, 1.0])
def test_surface_code_circuit_samples(r_e: float) -> None:
    d = 3
    circuit = build(d, d, WuSamplingInjector(d, WuNoiseParams(0.01, r_e, p_m=0.01)))
    samples = circuit.compile_detector_sampler(seed=SEED).sample(100)
    assert samples.shape == (100, circuit.num_detectors)


# --- Parameter validation and helpers. ---


def test_num_code_qubits_matches_builder() -> None:
    for d in (1, 3, 5):
        assert num_code_qubits(d) == 2 * d * d - 1
        assert build(d, d, WuSamplingInjector(d, WuNoiseParams(0.0, 0.0))).num_qubits == (
            num_code_qubits(d)
        )
    with pytest.raises(ValueError):
        num_code_qubits(0)


def test_wu_noise_params_rates_and_validation() -> None:
    params = WuNoiseParams(0.02, 0.25)
    assert params.p_e == pytest.approx(0.005)
    assert params.p_p == pytest.approx(0.015)
    assert params.meas == 0.0
    assert WuNoiseParams(0.02, 0.25, p_m=0.1, p_m_follows_p=True).meas == 0.02
    for p, r_e, p_m in ((-0.1, 0.5, 0.0), (1.1, 0.5, 0.0), (0.1, 1.5, 0.0), (0.1, 0.5, -0.01)):
        with pytest.raises(ValueError):
            WuNoiseParams(p, r_e, p_m=p_m)
