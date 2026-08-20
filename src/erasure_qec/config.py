"""Shared configuration dataclasses (PLAN.md §1).

``NoiseParams`` (M4) parameterizes the biased-erasure noise model;
``ExperimentConfig`` + :func:`load_experiment_config` (M7) parse the
``experiments/configs/*.yaml`` sweep definitions. Rounds are always T = d
(PLAN §9), so the config carries distances only.
"""

from dataclasses import dataclass
from pathlib import Path

import yaml


def _check_probability(value: float, name: str) -> None:
    if not 0.0 <= value <= 1.0:
        raise ValueError(f"{name} must be in [0, 1], got {value}")


@dataclass(frozen=True)
class NoiseParams:
    """Physical noise parameters for the biased-erasure model (PLAN.md §5).

    ``p_meas``, ``p_reset``, and ``p_idle`` default to ``p`` for the uniform
    sweep; pass them explicitly to decouple measurement/reset/idle noise from
    the two-qubit gate error budget. Use the :attr:`meas`/:attr:`reset`/
    :attr:`idle` properties to read the resolved (never-``None``) rates.

    ``convert_idle`` (default ``False``, preserving the gate-only erasure model
    as a named baseline) also erasure-converts the *idle* budget: a fraction
    ``r_e`` of each idle qubit's ``DEPOLARIZE1`` becomes a heralded erasure. This
    is the physically-motivated case for metastable ¹⁷¹Yb, where decay during
    idling is exactly what fluorescence monitoring heralds, and it raises the
    circuit-wide heralded fraction toward ``r_e`` (measurement/reset stay
    unheralded either way -- see analysis.dem_stats.heralded_fraction).
    """

    p: float
    r_e: float = 0.0
    p_meas: float | None = None
    p_reset: float | None = None
    p_idle: float | None = None
    convert_idle: bool = False

    def __post_init__(self) -> None:
        _check_probability(self.p, "p")
        _check_probability(self.r_e, "r_e")
        for name, value in (
            ("p_meas", self.p_meas),
            ("p_reset", self.p_reset),
            ("p_idle", self.p_idle),
        ):
            if value is not None:
                _check_probability(value, name)

    @property
    def meas(self) -> float:
        """Resolved measurement-flip probability (``p_meas`` or ``p``)."""
        return self.p if self.p_meas is None else self.p_meas

    @property
    def reset(self) -> float:
        """Resolved reset-flip probability (``p_reset`` or ``p``)."""
        return self.p if self.p_reset is None else self.p_reset

    @property
    def idle(self) -> float:
        """Resolved idle-depolarizing probability (``p_idle`` or ``p``)."""
        return self.p if self.p_idle is None else self.p_idle


@dataclass(frozen=True)
class ExperimentConfig:
    """One threshold-sweep definition (PLAN.md §9), loaded from YAML.

    ``p_values`` is the fully resolved, log-spaced physical error rate grid;
    rounds are implicitly T = d for every distance in ``distances``.
    """

    name: str
    r_e: float
    distances: tuple[int, ...]
    p_values: tuple[float, ...]
    max_shots: int
    max_errors: int
    decoders: tuple[str, ...]
    lambda_p: float  # fixed sub-threshold p for the Lambda scan (§10)

    def __post_init__(self) -> None:
        _check_probability(self.r_e, "r_e")
        _check_probability(self.lambda_p, "lambda_p")
        if not self.distances or any(d < 3 or d % 2 == 0 for d in self.distances):
            raise ValueError(f"distances must be odd and >= 3, got {self.distances}")
        for p in self.p_values:
            _check_probability(p, "p")
        if self.max_shots < 1 or self.max_errors < 1:
            raise ValueError("max_shots and max_errors must be >= 1")
        if not self.decoders:
            raise ValueError("at least one decoder is required")


def _log_spaced(log10_min: float, log10_max: float, count: int) -> tuple[float, ...]:
    """``count`` log-spaced points in [10**log10_min, 10**log10_max]."""
    if count < 2:
        raise ValueError(f"need at least 2 p points, got {count}")
    step = (log10_max - log10_min) / (count - 1)
    return tuple(10.0 ** (log10_min + i * step) for i in range(count))


def load_experiment_config(path: str | Path) -> ExperimentConfig:
    """Parse an ``experiments/configs/*.yaml`` sweep definition."""
    with open(path, encoding="utf-8") as f:
        raw = yaml.safe_load(f)
    if not isinstance(raw, dict):
        raise ValueError(f"{path}: expected a YAML mapping at top level")
    base_grid = _log_spaced(
        float(raw["p_log10_min"]), float(raw["p_log10_max"]), int(raw["p_points"])
    )
    # Optional additive dense band (e.g. extra fit-resolution points through the
    # crossing); merged with the log grid, sorted, and de-duplicated.
    extra = tuple(float(x) for x in raw.get("extra_p", ()))
    p_values = tuple(sorted(set(base_grid) | set(extra)))
    return ExperimentConfig(
        name=str(raw["name"]),
        r_e=float(raw["r_e"]),
        distances=tuple(int(d) for d in raw["distances"]),
        p_values=p_values,
        max_shots=int(raw["max_shots"]),
        max_errors=int(raw["max_errors"]),
        decoders=tuple(str(dec) for dec in raw["decoders"]),
        lambda_p=float(raw["lambda_p"]),
    )
