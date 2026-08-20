"""Regenerate committed synthetic sweep fixtures (PLAN.md §10 support).

Writes deterministic, ansatz-exact sinter CSVs to ``tests/fixtures/synthetic_*``
for the three r_e regimes, plus one adversarial fixture. These drive the
byte-stable plotting/analysis tests; they are NOT measurements (see the
provenance warning in ``analysis/synthetic.py``: they are generated from the
same collapse variable ``fit_threshold`` inverts). The real Monte-Carlo sweeps
are committed under ``data/`` and pinned under ``tests/fixtures/real_*.csv``.
Run:

    uv run python experiments/make_synthetic_fixtures.py
"""

from pathlib import Path

from erasure_qec.analysis.synthetic import (
    AnsatzParams,
    write_adversarial_csv,
    write_synthetic_csv,
)

FIXTURES = Path(__file__).resolve().parent.parent / "tests" / "fixtures"

# 13 log-spaced p points in [1e-3, 1e-1] (same grid the real configs use).
P_VALUES = tuple(10.0 ** (-3.0 + i * (2.0 / 12.0)) for i in range(13))
DISTANCES = (3, 5, 7, 9, 11)
SHOTS = 200_000

# (name, r_e, injected p_th, blind-vs-herald penalty) — p_th slides with r_e.
SWEEPS = [
    ("synthetic_baseline_pauli", 0.0, 0.010, 1.15),
    ("synthetic_erasure_r50", 0.5, 0.024, 1.35),
    ("synthetic_erasure_r98", 0.98, 0.045, 1.60),
]


def main() -> None:
    FIXTURES.mkdir(parents=True, exist_ok=True)
    for name, r_e, p_th, penalty in SWEEPS:
        ansatz = AnsatzParams(p_th=p_th, nu=1.5, a=0.09, b=22.0)
        path = write_synthetic_csv(
            FIXTURES / f"{name}.csv",
            ansatz,
            decoders=("herald_mwpm", "blind_mwpm"),
            distances=DISTANCES,
            p_values=P_VALUES,
            shots=SHOTS,
            r_e=r_e,
            blind_penalty=penalty,
        )
        print(f"wrote {path}")

    # Adversarial: inverted ordering + saturated tail; the fitter must reject it.
    adv = write_adversarial_csv(
        FIXTURES / "synthetic_adversarial.csv",
        decoders=("herald_mwpm", "blind_mwpm"),
        distances=DISTANCES,
        p_values=P_VALUES,
        shots=SHOTS,
        r_e=0.5,
    )
    print(f"wrote {adv}")


if __name__ == "__main__":
    main()
