"""``NoiseParams`` -> concrete per-instruction channel probabilities (PLAN.md §5).

Kept separate from ``injector.py`` so the arithmetic (what probability each
stim instruction gets) is unit-testable independent of circuit construction.
"""

from dataclasses import dataclass

from erasure_qec.config import NoiseParams

# A heralded erasure replaces the qubit with the maximally mixed state, i.e. a
# uniform Pauli in {I, X, Y, Z}; it is a *non-identity* error only 3/4 of the
# time (the I outcome still heralds but causes no syndrome).
_ERASURE_NONIDENTITY_FRACTION = 0.75


@dataclass(frozen=True)
class ChannelRates:
    """Concrete probabilities for each noisy instruction the builder emits."""

    depolarize2: float  # DEPOLARIZE2 on the two-qubit-gate pair
    herald: float  # HERALDED_ERASE on each qubit of the pair, independently
    meas_flip: float  # X_ERROR immediately before M/MR
    reset_flip: float  # X_ERROR immediately after R
    idle_depolarize: float  # DEPOLARIZE1 on qubits idle during a TICK


def channel_rates(params: NoiseParams) -> ChannelRates:
    """Convert ``NoiseParams`` into the concrete per-instruction probabilities.

    The two-qubit-gate error budget is held constant across ``r_e``: PLAN.md §0
    defines ``r_e`` as the fraction of the physical error budget per two-qubit
    gate that is converted into heralded erasures, so the total non-identity
    error probability per gate must equal ``p`` for every ``r_e``.

    - Residual Pauli: ``DEPOLARIZE2(p * (1 - r_e))`` — non-identity mass
      ``p * (1 - r_e)``.
    - Erasure: ``HERALDED_ERASE(q)`` independently on each of the two qubits.
      An erasure is a non-identity error only 3/4 of the time, so the two
      qubits contribute ``2 * q * (3/4)`` of non-identity mass to first order in
      ``p``. Setting that equal to the erasure share ``p * r_e`` gives
      ``q = p * r_e / (2 * 3/4) = (2/3) * p * r_e``.

    To first order in ``p`` the per-gate non-identity probability is then ``p``
    for every ``r_e``, and ``r_e`` is exactly the heralded fraction of the
    two-qubit-gate error budget. Measurement, reset, and idle errors are *not*
    erasure-converted (they model separate physical processes), so the
    circuit-wide heralded fraction is below ``r_e``.

    This corrects the earlier ``p * r_e / 2`` rate, under which the per-gate
    non-identity budget shrank to ``p * (1 - r_e/4)`` as ``r_e`` grew, inflating
    the apparent threshold (see docs/AUDIT.md).
    """
    return ChannelRates(
        depolarize2=params.p * (1.0 - params.r_e),
        herald=params.p * params.r_e / (2.0 * _ERASURE_NONIDENTITY_FRACTION),
        meas_flip=params.meas,
        reset_flip=params.reset,
        idle_depolarize=params.idle,
    )


def nonidentity_pauli_probability(params: NoiseParams) -> float:
    """Exact probability a two-qubit-gate location yields >=1 non-identity Pauli.

    The gate location applies three *independent* channels: one ``DEPOLARIZE2``
    on the pair and one ``HERALDED_ERASE`` on each qubit. ``DEPOLARIZE2(q)`` is
    non-identity with probability exactly ``q``; ``HERALDED_ERASE(h)`` is
    non-identity with probability ``(3/4) h`` -- the ``I/4`` branch heralds but
    causes no Pauli. Independence gives

        P(non-identity) = 1 - (1 - dep) * (1 - (3/4) h)^2 .

    This is the tested error-budget axis (tests/test_noise_model.py::
    test_error_budget_invariance). Under the constant-budget convention
    (``herald = (2/3) p r_e``) the *linear* part is
    exactly ``p`` for every ``r_e``; the only ``r_e`` dependence is the O(p^2)
    inclusion-exclusion overlap -- ``p`` at ``r_e = 0`` falling to ``p - p^2/4``
    at ``r_e = 1``. Contrast the earlier ``p r_e / 2`` rate, whose budget fell
    *linearly* to ``p (1 - r_e/4)`` (a 25% shrink at ``r_e = 0.98``); ``p`` was
    not an iso-noise axis then. See docs/AUDIT.md.
    """
    rates = channel_rates(params)
    return 1.0 - (1.0 - rates.depolarize2) * (
        1.0 - _ERASURE_NONIDENTITY_FRACTION * rates.herald
    ) ** 2
