"""Phase-4 deliverable: the herald-vs-blind threshold-difference statistic.

Emits the numbers docs/AUDIT.md "Paired-decoder separation" reports:

1. Whether the sweeps are genuinely paired (evidence: shared shot counts).
2. Before/after marginal p_th CIs at the old (200) and new (1000) n_boot.
3. Delta = p_th(blind) - p_th(herald) and its CI for every r_e, with
   n_paired_failed and the herald/blind bootstrap correlation.
4. The r_e = 0 control (Delta must be consistent with zero).
5. Seed-stability spread of the Delta CI over several seeds.

Run:  uv run python scripts/paired_separation.py
"""

from __future__ import annotations

from erasure_qec.analysis.statistics import SweepPoint, load_sweep
from erasure_qec.analysis.threshold_fit import (
    bootstrap_threshold_difference,
    fit_threshold,
)

SWEEPS = [
    (0.0, "data/baseline_pauli.csv"),
    (0.5, "data/erasure_r50.csv"),
    (0.98, "data/erasure_r98.csv"),
]
D_MIN = 7


def _split(csv: str) -> tuple[list[SweepPoint], list[SweepPoint]]:
    pts = load_sweep(csv)
    return (
        [p for p in pts if p.decoder == "herald_mwpm"],
        [p for p in pts if p.decoder == "blind_mwpm"],
    )


def _pairing_evidence() -> None:
    print("== 1. Are the sweeps paired? (shared shot counts per (p, d)) ==")
    for r_e, csv in SWEEPS:
        h, b = _split(csv)
        by: dict[tuple[float, int], dict[str, int]] = {}
        for p in h + b:
            by.setdefault((round(p.p, 6), p.d), {})[p.decoder] = p.shots
        both = [v for v in by.values() if len(v) == 2]
        differ = sum(1 for v in both if v["herald_mwpm"] != v["blind_mwpm"])
        print(f"  r_e={r_e}: {differ}/{len(both)} (p,d) pairs have DIFFERING "
              f"herald/blind shot counts -> {'INDEPENDENT' if differ else 'shared?'}")
    print("  => sinter samples each decoder independently; CSV stores only marginal")
    print("     (shots, errors). Real data uses the UNPAIRED difference bootstrap.\n")


def _marginal_ci_before_after() -> None:
    print("== 2. Marginal p_th 95% CI: n_boot 200 (before) vs 1000 (after) ==")
    print(f"  {'r_e':>4} {'decoder':8} {'p_th':>7} {'CI@200':>16} {'CI@1000':>16}")
    for r_e, csv in SWEEPS:
        h, b = _split(csv)
        for dec, pts in (("herald", h), ("blind", b)):
            f200 = fit_threshold(pts, d_min=D_MIN, n_boot=200, seed=0)
            f1000 = fit_threshold(pts, d_min=D_MIN, n_boot=1000, seed=0)
            if not f200.converged:
                print(f"  {r_e:>4} {dec:8} not resolved")
                continue
            def ci(f: object) -> str:
                lo, hi = f.p_th_ci  # type: ignore[attr-defined]
                return f"[{lo*100:.2f},{hi*100:.2f}]"
            print(f"  {r_e:>4} {dec:8} {f1000.p_th*100:6.3f}% "
                  f"{ci(f200):>16} {ci(f1000):>16}")
    print()


def _delta_all() -> None:
    print("== 3+4. Delta = p_th(blind) - p_th(herald), 95% CI, n_boot=1000 ==")
    for r_e, csv in SWEEPS:
        h, b = _split(csv)
        r = bootstrap_threshold_difference(h, b, d_min=D_MIN, n_boot=1000, seed=0)
        tag = " (CONTROL: must include 0)" if r_e == 0.0 else ""
        if not r.converged:
            print(f"  r_e={r_e}: NON-RESULT ({r.message}){tag}")
            continue
        print(f"  r_e={r_e}: delta={r.delta*100:+.3f}%  "
              f"CI=[{r.delta_ci[0]*100:+.3f},{r.delta_ci[1]*100:+.3f}]%  "
              f"excludes_zero={r.excludes_zero}  corr={r.correlation:+.3f}  "
              f"failed={r.n_paired_failed}/{r.n_boot}{tag}")
    print()


def _seed_stability() -> None:
    print("== 5. Seed stability of the r_e=0.5 Delta CI (5 seeds, n_boot=1000) ==")
    h, b = _split("data/erasure_r50.csv")
    los, his = [], []
    for seed in range(5):
        r = bootstrap_threshold_difference(h, b, d_min=D_MIN, n_boot=1000, seed=seed)
        los.append(r.delta_ci[0])
        his.append(r.delta_ci[1])
        print(f"  seed={seed}: CI=[{r.delta_ci[0]*100:+.3f},{r.delta_ci[1]*100:+.3f}]%  "
              f"excludes_zero={r.excludes_zero}")
    print(f"  lower-endpoint spread: {(max(los)-min(los))*100:.3f}%  "
          f"upper-endpoint spread: {(max(his)-min(his))*100:.3f}%")


def main() -> None:
    _pairing_evidence()
    _marginal_ci_before_after()
    _delta_all()
    _seed_stability()


if __name__ == "__main__":
    main()
