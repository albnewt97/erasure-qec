"""Synthetic sweep generation for figure fixtures (PLAN.md §10 support).

Produces deterministic, sinter-format CSVs whose per-round logical rates follow
a physically-plausible data-collapse ``p_L = A * exp(B * x)`` in the scaling
variable ``x = (p - p_th) * d^(1/nu)``. This universal function has a clean
crossing at ``p = p_th`` (all distances meet at ``A``), fans out below
threshold and grows above, and stays monotone and positive everywhere — unlike
the bare quadratic ansatz, which is only valid in a narrow near-threshold
window. Near the crossing it reduces to that quadratic, so
:func:`erasure_qec.analysis.threshold_fit.fit_threshold` recovers ``p_th`` from
these fixtures within its window; the *exact-ansatz* recovery test constructs
quadratic points directly.

These fixtures stand in for the real Monte-Carlo sweeps (not yet collected) so
analysis figures render from committed data. The module never samples:
``errors`` is the rounded expected count, so the fixtures are byte-stable.
"""

import math
from dataclasses import dataclass
from pathlib import Path

import sinter

from erasure_qec.analysis.statistics import shot_p_l_from_per_round


@dataclass(frozen=True)
class AnsatzParams:
    """Ground-truth collapse parameters for a synthetic sweep.

    ``a`` is the per-round rate exactly at threshold (where all distances
    cross); ``b`` is the exponential slope in the scaling variable.
    """

    p_th: float
    nu: float
    a: float
    b: float

    def per_round_p_l(self, p: float, d: int) -> float:
        """Collapse per-round rate ``A exp(B x)``, clamped to ``[1e-6, 0.49]``."""
        x = (p - self.p_th) * d ** (1.0 / self.nu)
        return min(max(self.a * math.exp(self.b * x), 1e-6), 0.49)


def synthetic_stats(
    ansatz: AnsatzParams,
    *,
    decoder: str,
    distances: tuple[int, ...],
    p_values: tuple[float, ...],
    shots: int,
    r_e: float,
) -> list[tuple[int, int, int, float, float]]:
    """Return ``(d, p, errors, shots, r_e)`` rows with T = d rounds each."""
    rows: list[tuple[int, int, int, float, float]] = []
    for d in distances:
        rounds = d
        for p in p_values:
            p_round = ansatz.per_round_p_l(p, d)
            p_shot = shot_p_l_from_per_round(p_round, rounds)
            errors = int(round(shots * p_shot))
            rows.append((d, rounds, errors, p, r_e))
    return rows


def write_synthetic_csv(
    path: str | Path,
    ansatz: AnsatzParams,
    *,
    decoders: tuple[str, ...],
    distances: tuple[int, ...],
    p_values: tuple[float, ...],
    shots: int,
    r_e: float,
    blind_penalty: float = 1.0,
) -> Path:
    """Write a deterministic sinter-format CSV for the given decoders.

    ``blind_penalty`` (>= 1) inflates ``blind_mwpm``'s per-round rate so the
    ablation figure has a visible gap; ``herald_mwpm`` uses the ansatz as-is.
    """
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    lines = [sinter.CSV_HEADER]
    for decoder in decoders:
        penalty = blind_penalty if decoder == "blind_mwpm" else 1.0
        scaled = AnsatzParams(
            p_th=ansatz.p_th,
            nu=ansatz.nu,
            a=min(ansatz.a * penalty, 0.4),
            b=ansatz.b,
        )
        for d, rounds, errors, p, r_e_val in synthetic_stats(
            scaled, decoder=decoder, distances=distances, p_values=p_values,
            shots=shots, r_e=r_e,
        ):
            stat = sinter.TaskStats(
                strong_id=f"{decoder}-d{d}-p{p:.6g}-re{r_e_val:g}",
                decoder=decoder,
                json_metadata={"d": d, "p": p, "r_e": r_e_val, "rounds": rounds},
                shots=shots,
                errors=errors,
            )
            lines.append(stat.to_csv_line())
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out
