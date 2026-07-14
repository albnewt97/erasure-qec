"""``NoiseParams`` -> concrete per-instruction channel probabilities (PLAN.md §5).

Kept separate from ``injector.py`` so the arithmetic (what probability each
stim instruction gets) is unit-testable independent of circuit construction.
"""

from dataclasses import dataclass

from erasure_qec.config import NoiseParams


@dataclass(frozen=True)
class ChannelRates:
    """Concrete probabilities for each noisy instruction the builder emits."""

    depolarize2: float  # DEPOLARIZE2 on the two-qubit-gate pair
    herald: float  # HERALDED_ERASE on each qubit of the pair, independently
    meas_flip: float  # X_ERROR immediately before M/MR
    reset_flip: float  # X_ERROR immediately after R
    idle_depolarize: float  # DEPOLARIZE1 on qubits idle during a TICK


def channel_rates(params: NoiseParams) -> ChannelRates:
    """Convert ``NoiseParams`` into the concrete channel probabilities of §5.

    Per two-qubit gate on (a, b): ``DEPOLARIZE2(p * (1 - r_e))`` is the
    residual Pauli component, and ``HERALDED_ERASE(p * r_e / 2)`` is applied
    independently to each of a and b.
    """
    return ChannelRates(
        depolarize2=params.p * (1.0 - params.r_e),
        herald=params.p * params.r_e / 2.0,
        meas_flip=params.meas,
        reset_flip=params.reset,
        idle_depolarize=params.idle,
    )
