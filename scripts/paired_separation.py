"""Phase-4 deliverable: the herald-vs-blind threshold-difference statistic.

Emits the numbers docs/AUDIT.md "Paired-decoder separation" reports:

1. Whether the sweeps are genuinely paired (evidence: shared shot counts).
2. Before/after marginal p_th CIs at the old (200) and new (1000) n_boot.
3. Delta = p_th(blind) - p_th(herald) and its CI for every r_e, with
   n_paired_failed and the herald/blind bootstrap correlation.
4. The r_e = 0 control (Delta must be consistent with zero).
5. Seed-stability spread of the Delta CI over several seeds.

Run:  uv run python scripts/paired_separation.py                    (~72 s)
      uv run python scripts/paired_separation.py --full-stability   (~4.7 min)

The default sweeps five seeds, which is enough to show the r_e = 0.5 Delta CI
is not seed-fragile.  `--full-stability` replaces that section with the full
ten-seed table -- marginal CI spreads, `rel_ci_width`, `n_paired_failed`, the
tipping point, the herald/blind correlation, and the r_e = 0.98 blind upper
endpoint with its grid-pin check.  Those are the quantities quoted in the
docs/AUDIT.md "Reproducibility and environment sensitivity" seed-spread column;
running the ten-seed table costs ~24 s per seed, which is why it is opt-in.
"""

from __future__ import annotations

import argparse

import numpy as np
from scipy.stats import ks_2samp

from erasure_qec.analysis.statistics import SweepPoint, load_sweep
from erasure_qec.analysis.threshold_fit import (
    FitResult,
    ThresholdDifference,
    bootstrap_threshold_difference,
    directional_sensitivity_ci,
    failed_vs_success_pth,
    failure_guard_breakdown,
    fit_threshold,
    partial_information_implied_delta,
    tipping_point_discards,
)

SWEEPS = [
    (0.0, "data/baseline_pauli.csv"),
    (0.5, "data/erasure_r50.csv"),
    (0.98, "data/erasure_r98.csv"),
]
D_MIN = 7
_FULL_SEEDS = 10  # the seed range docs/AUDIT.md's seed-spread column quotes


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
        print(
            f"  r_e={r_e}: {differ}/{len(both)} (p,d) pairs have DIFFERING "
            f"herald/blind shot counts -> {'INDEPENDENT' if differ else 'shared?'}"
        )
    print("  => sinter samples each decoder independently; CSV stores only marginal")
    print("     (shots, errors). Real data uses the UNPAIRED difference bootstrap.\n")


def _marginal_ci_before_after() -> None:
    print("== 2. Marginal p_th 95% CI: n_boot 200 (before) vs 1000 (after) ==")
    print(f"  {'r_e':>4} {'decoder':8} {'p_th':>7} {'CI@200':>16} {'CI@1000':>16}")

    def ci(f: FitResult) -> str:
        lo, hi = f.p_th_ci
        return f"[{lo * 100:.2f},{hi * 100:.2f}]" if f.converged else "-"

    for r_e, csv in SWEEPS:
        h, b = _split(csv)
        for dec, pts in (("herald", h), ("blind", b)):
            f200 = fit_threshold(pts, d_min=D_MIN, n_boot=200, seed=0)
            f1000 = fit_threshold(pts, d_min=D_MIN, n_boot=1000, seed=0)
            # Gate the quoted p_th on `resolved`, not `converged`: r_e=0.98 blind
            # converges but its CI is wider than p_th, so it is not a threshold.
            if not f1000.resolved:
                note = f"(rel_ci {f1000.rel_ci_width:.2f})" if f1000.converged else ""
                print(
                    f"  {r_e:>4} {dec:8} {'not resolved':>7} {ci(f200):>16} {ci(f1000):>16}  {note}"
                )
                continue
            print(f"  {r_e:>4} {dec:8} {f1000.p_th * 100:6.3f}% {ci(f200):>16} {ci(f1000):>16}")
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
        print(
            f"  r_e={r_e}: delta={r.delta * 100:+.3f}%  "
            f"CI=[{r.delta_ci[0] * 100:+.3f},{r.delta_ci[1] * 100:+.3f}]%  "
            f"excludes_zero={r.excludes_zero}  corr={r.correlation:+.3f}  "
            f"failed={r.n_paired_failed}/{r.n_boot}{tag}"
        )
    print()


def _seed_stability(n_seeds: int = 5) -> None:
    print(f"== 5. Seed stability of the r_e=0.5 Delta CI ({n_seeds} seeds, n_boot=1000) ==")
    h, b = _split("data/erasure_r50.csv")
    los, his = [], []
    for seed in range(n_seeds):
        r = bootstrap_threshold_difference(h, b, d_min=D_MIN, n_boot=1000, seed=seed)
        los.append(r.delta_ci[0])
        his.append(r.delta_ci[1])
        print(
            f"  seed={seed}: CI=[{r.delta_ci[0] * 100:+.3f},{r.delta_ci[1] * 100:+.3f}]%  "
            f"excludes_zero={r.excludes_zero}"
        )
    print(
        f"  lower-endpoint spread: {(max(los) - min(los)) * 100:.3f}%  "
        f"upper-endpoint spread: {(max(his) - min(his)) * 100:.3f}%"
    )
    print(
        f"  (default is 5 seeds; --full-stability runs the documented {_FULL_SEEDS}-seed table)\n"
    )


def _full_stability_table() -> None:
    """The ten-seed table behind docs/AUDIT.md's "seed spread" column.

    Recomputes every seed-dependent quantity that section quotes, so a reviewer
    can check each one rather than taking it on trust.  ~24 s per seed.
    """
    seeds = range(_FULL_SEEDS)
    print(f"== 5. FULL seed-stability table (seeds 0-{_FULL_SEEDS - 1}, n_boot=1000) ==")
    print("   Reproduces the 'seed spread' column of docs/AUDIT.md")
    print("   'Reproducibility and environment sensitivity'.\n")

    h0, b0 = _split("data/baseline_pauli.csv")
    h5, b5 = _split("data/erasure_r50.csv")
    _, b98 = _split("data/erasure_r98.csv")

    marginals = [
        ("r_e=0 herald", h0),
        ("r_e=0 blind", b0),
        ("r_e=0.5 herald", h5),
        ("r_e=0.5 blind", b5),
        ("r_e=0.98 blind", b98),
    ]
    los: dict[str, list[float]] = {n: [] for n, _ in marginals}
    his: dict[str, list[float]] = {n: [] for n, _ in marginals}
    rel98: list[float] = []
    for seed in seeds:
        for name, pts in marginals:
            f = fit_threshold(pts, d_min=D_MIN, n_boot=1000, seed=seed)
            los[name].append(f.p_th_ci[0])
            his[name].append(f.p_th_ci[1])
            if name == "r_e=0.98 blind":
                rel98.append(f.rel_ci_width)

    # Difference bootstraps: r_e=0.5 (the significance claim) and the r_e=0 control.
    d5: dict[str, list[float]] = {k: [] for k in ("lo", "hi", "failed", "tip", "corr")}
    d0: dict[str, list[float]] = {k: [] for k in ("lo", "hi")}
    excl5: list[bool] = []
    for seed in seeds:
        r = bootstrap_threshold_difference(h5, b5, d_min=D_MIN, n_boot=1000, seed=seed)
        d5["lo"].append(r.delta_ci[0])
        d5["hi"].append(r.delta_ci[1])
        d5["failed"].append(r.n_paired_failed)
        d5["tip"].append(tipping_point_discards(r))
        d5["corr"].append(r.correlation)
        excl5.append(r.excludes_zero)
        c = bootstrap_threshold_difference(h0, b0, d_min=D_MIN, n_boot=1000, seed=seed)
        d0["lo"].append(c.delta_ci[0])
        d0["hi"].append(c.delta_ci[1])

    def spread(v: list[float]) -> str:
        return f"{(max(v) - min(v)) * 100:.3f}"

    print(f"  {'quantity (%)':<32} {'seed-0':>9} {'min':>9} {'max':>9} {'spread':>8}")
    for name, _ in marginals:
        for half, vals in (("lo", los[name]), ("hi", his[name])):
            print(
                f"  {name + ' CI ' + half:<32} {vals[0] * 100:9.3f} "
                f"{min(vals) * 100:9.3f} {max(vals) * 100:9.3f} {spread(vals):>8}"
            )
    for label, vals in (
        ("Delta CI (r_e=0.5) lo", d5["lo"]),
        ("Delta CI (r_e=0.5) hi", d5["hi"]),
        ("Delta CI (r_e=0) lo", d0["lo"]),
        ("Delta CI (r_e=0) hi", d0["hi"]),
    ):
        print(
            f"  {label:<32} {vals[0] * 100:9.3f} {min(vals) * 100:9.3f} "
            f"{max(vals) * 100:9.3f} {spread(vals):>8}"
        )

    print()
    print(
        f"  n_paired_failed (r_e=0.5):      seed-0={int(d5['failed'][0])}  "
        f"range {int(min(d5['failed']))}-{int(max(d5['failed']))}"
    )
    print(
        f"  tipping-point discards (0.5):   seed-0={int(d5['tip'][0])}  "
        f"range {int(min(d5['tip']))}-{int(max(d5['tip']))}"
    )
    print(
        f"  herald/blind corr (r_e=0.5):    seed-0={d5['corr'][0]:+.3f}  "
        f"range {min(d5['corr']):+.3f} ... {max(d5['corr']):+.3f}"
    )
    print(
        f"  r_e=0.98 blind rel_ci_width:    seed-0={rel98[0]:.2f}  "
        f"range {min(rel98):.2f}-{max(rel98):.2f}"
    )
    print(
        f"  all {_FULL_SEEDS} seeds exclude zero (r_e=0.5 Delta): {all(excl5)} "
        f"({sum(excl5)}/{len(excl5)})"
    )

    # The r_e=0.98 blind upper endpoint: is it a fitted crossing or a grid pin?
    grid = sorted({round(pt.p, 6) for pt in b98 if pt.d >= D_MIN})
    ups = his["r_e=0.98 blind"]
    print(
        f"\n  r_e=0.98 blind CI upper by seed (%): "
        f"{' '.join(f'{u * 100:.3f}' for u in sorted(ups))}"
    )
    pinned = sorted(((i, u) for i, u in enumerate(ups) if round(u, 6) in grid), key=lambda t: t[1])
    hits = ", ".join(f"{u * 100:.3f}% (p={round(u, 6)}, seed {i})" for i, u in pinned)
    print(f"  sweep grid points (d>={D_MIN}) coinciding: {hits}")
    print(f"  => {len(pinned)}/{len(ups)} seeds put the upper endpoint EXACTLY on a sweep")
    print("     grid point: that endpoint is set by bound-pinned replicates, not by a")
    print("     fitted crossing. No digit of it is a measurement.\n")


def _discard_diagnostics() -> None:
    print("== 6. Discarded-replicate diagnostics (conditioning on convergence) ==")
    for r_e, csv in [(0.0, "data/baseline_pauli.csv"), (0.5, "data/erasure_r50.csv")]:
        h, b = _split(csv)
        r = bootstrap_threshold_difference(h, b, d_min=D_MIN, n_boot=1000, seed=0)
        if not r.converged:
            print(f"  r_e={r_e}: {r.message}")
            continue
        print(f"  r_e={r_e}: {r.n_paired_failed}/{r.n_boot} discarded")
        breakdown = failure_guard_breakdown(r)
        print("    by guard: " + ", ".join(f"{k}={v}" for k, v in sorted(breakdown.items())))
        # Missing-at-random check on blind p_th (the concern's direction).
        failed_b, success_b = failed_vs_success_pth(r, "blind")
        if len(failed_b) and len(success_b):
            ks = ks_2samp(failed_b, success_b)
            print(
                f"    blind p_th  success: median={np.median(success_b) * 100:.3f}% "
                f"q90={np.percentile(success_b, 90) * 100:.3f}%   "
                f"failed: median={np.median(failed_b) * 100:.3f}% "
                f"q10={np.percentile(failed_b, 10) * 100:.3f}%"
            )
            higher = np.median(failed_b) > np.median(success_b)
            print(
                f"    KS(failed vs success blind p_th): stat={ks.statistic:.3f} "
                f"p={ks.pvalue:.2e}  -> discards {'ARE' if ks.pvalue < 0.05 else 'are NOT'} "
                f"MAR-inconsistent; failed blind p_th {'higher' if higher else 'not higher'}"
            )
        _sensitivity_audit(r)
    print()


def _sensitivity_audit(r: ThresholdDifference) -> None:
    """Audit the old imputation (Task 1) and report the tipping-point bound that
    replaces it (Task 2)."""
    deltas = np.array(r.blind_pth_draws) - np.array(r.herald_pth_draws)
    upper = float(np.percentile(deltas, 97.5))
    # Task 1: what the imputation actually draws, and whether the CI was right.
    thr = float(np.percentile(deltas, 90.0))
    top = deltas[deltas >= thr]
    imputed = np.random.default_rng(0).choice(top, size=r.n_paired_failed, replace=True)
    ci, excl = directional_sensitivity_ci(r, seed=0, top_fraction=0.10)
    q = np.percentile(imputed, [0, 25, 50, 75, 100]) * 100
    print(
        f"    [T1] imputed {r.n_paired_failed} vals: min/25/50/75/max="
        f"{q[0]:+.3f}/{q[1]:+.3f}/{q[2]:+.3f}/{q[3]:+.3f}/{q[4]:+.3f}%  "
        f"distinct={len(np.unique(imputed))}  >{upper * 100:+.3f}%: "
        f"{int((imputed > upper).sum())}/{r.n_paired_failed}"
    )
    print(
        f"    [T1] imputation CI=[{ci[0] * 100:+.3f},{ci[1] * 100:+.3f}]% excl0={excl} "
        f"(reproduces the earlier number; WEAK/arbitrary -- see below)"
    )
    # Task 2: tipping-point + partial-information bound.
    tip = tipping_point_discards(r)
    implied, unbounded = partial_information_implied_delta(r)
    n_ge0 = int((implied >= 0.0).sum())
    herald_med = float(np.median(r.herald_pth_draws))
    fb = np.array([f.blind_p_th for f in r.failures])
    fb = fb[np.isfinite(fb)]
    n_reach = int((fb >= herald_med).sum())
    print(
        f"    [T2a] tipping point: claim survives unless >= {tip} of "
        f"{r.n_paired_failed} discards give delta>=0"
    )
    print(
        f"    [T2b] implied delta from recorded p_th: >=0 in {n_ge0}/{len(implied)}; "
        f"a discard needs blind p_th >= herald centre {herald_med * 100:.3f}%"
    )
    print(f"          (failed-blind median {np.median(fb) * 100:.3f}%, {n_reach} reach it)")
    print(f"    [T2c] discards with no converged decoder (unbounded): {unbounded}")
    if r.excludes_zero:
        verdict = "SURVIVES" if n_ge0 < tip else "DOES NOT SURVIVE"
        print(f"    => significance claim: tipping {tip} vs implied>=0 {n_ge0} -> {verdict}")
    else:
        print(
            "    => control: delta already consistent with zero; no significance "
            "to defend (tipping/implied are informational)"
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--full-stability",
        action="store_true",
        help=f"run the documented {_FULL_SEEDS}-seed stability table instead of the "
        "5-seed default (~4.7 min rather than ~72 s)",
    )
    args = parser.parse_args()
    _pairing_evidence()
    _marginal_ci_before_after()
    _delta_all()
    if args.full_stability:
        _full_stability_table()
    else:
        _seed_stability()
    _discard_diagnostics()


if __name__ == "__main__":
    main()
