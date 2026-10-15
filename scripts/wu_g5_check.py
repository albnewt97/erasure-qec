"""G5 gate: the (9,13) per-shot crossing of wu_fig3_re000 is at least 0.90%.

Reads data/wu/wu_fig3_re000.csv (rows with convention wu2022, config
wu_fig3_re000, decoder wu_herald_mwpm; sinter merges rows that share a
strong_id), computes the first crossing of the d = 9 and d = 13 curves of
P_L_shot = errors / shots with
:func:`erasure_qec.analysis.crossing.bootstrap_pairwise_crossing`
(parametric bootstrap, 95% percentile CI, n_boot = 1000, seed = 0; binomial
noise only, interpolation error inside the bracket not included) and prints
PASS if the point estimate is >= 0.0090, else FAIL (also FAIL when there is no
crossing). Exit status 0 on PASS, 1 on FAIL. Run:

    uv run python scripts/wu_g5_check.py
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import sinter

from erasure_qec.analysis.crossing import bootstrap_pairwise_crossing
from erasure_qec.analysis.statistics import SweepPoint

DATA = Path("data/wu/wu_fig3_re000.csv")
CONFIG = "wu_fig3_re000"
DECODER = "wu_herald_mwpm"
D1, D2 = 9, 13
N_BOOT = 1000
SEED = 0
THRESHOLD = 0.0090


def load_points(path: Path) -> list[SweepPoint]:
    """``SweepPoint`` per (d, p) of ``DECODER`` in ``path``, with ``r_e`` from ``R_e``."""
    points: list[SweepPoint] = []
    for s in sinter.read_stats_from_csv_files(str(path)):
        meta = s.json_metadata
        if (
            s.decoder != DECODER
            or s.shots <= 0
            or not isinstance(meta, dict)
            or meta.get("convention") != "wu2022"
            or meta.get("config") != CONFIG
        ):
            continue
        points.append(
            SweepPoint(
                decoder=DECODER,
                d=int(meta["d"]),
                rounds=int(meta["rounds"]),
                p=float(meta["p"]),
                r_e=float(meta["R_e"]),
                shots=int(s.shots),
                errors=int(s.errors),
            )
        )
    return points


def main() -> int:
    if not DATA.exists():
        raise SystemExit(f"{DATA} not found")
    points = load_points(DATA)
    est, ci, n_no = bootstrap_pairwise_crossing(
        points, D1, D2, metric="shot", n_boot=N_BOOT, seed=SEED
    )
    n_d = {d: sum(pt.d == d for pt in points) for d in (D1, D2)}
    print(f"{DATA}: {n_d[D1]} points at d={D1}, {n_d[D2]} at d={D2} ({DECODER})")
    print(f"({D1},{D2}) per-shot crossing: {est.value!r}")
    print(f"95% bootstrap percentile CI (n_boot={N_BOOT}, seed={SEED}): {ci}")
    print(f"bracket: {est.bracket}; n_sign_changes: {est.n_sign_changes}")
    print(f"replicates without a crossing: {n_no}")
    passed = not math.isnan(est.value) and est.value >= THRESHOLD
    print(f"G5 (crossing >= {THRESHOLD}): {'PASS' if passed else 'FAIL'}")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
