"""Regenerate committed synthetic sweep fixtures (PLAN.md §10 support).

Writes deterministic, ansatz-exact sinter CSVs to ``tests/fixtures/`` for the
three r_e regimes. These stand in for the real Monte-Carlo sweeps (not yet
collected) so the README/analysis figures render from committed data. The
injected ``p_th`` slides with r_e (~1% Pauli baseline -> ~4-5% erasure-limited),
mirroring the physics the real sweeps should show. Run:

    uv run python experiments/make_synthetic_fixtures.py
"""

from pathlib import Path

from erasure_qec.analysis.synthetic import AnsatzParams, write_synthetic_csv

FIXTURES = Path(__file__).resolve().parent.parent / "tests" / "fixtures"

# 13 log-spaced p points in [1e-3, 1e-1] (same grid the real configs use).
P_VALUES = tuple(10.0 ** (-3.0 + i * (2.0 / 12.0)) for i in range(13))
DISTANCES = (3, 5, 7, 9, 11)
SHOTS = 200_000

# (name, r_e, injected p_th, blind-vs-herald penalty) — p_th slides with r_e.
SWEEPS = [
    ("baseline_pauli", 0.0, 0.010, 1.15),
    ("erasure_r50", 0.5, 0.024, 1.35),
    ("erasure_r98", 0.98, 0.045, 1.60),
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


if __name__ == "__main__":
    main()
