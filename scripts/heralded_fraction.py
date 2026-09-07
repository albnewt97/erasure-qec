"""Report the circuit-wide heralded fraction of the DEM error budget.

The honest "how much of the noise is erasure-converted" number for a config --
which is NOT r_e (measurement/reset stay unheralded; idle only with
convert_idle). Quote this wherever r_e is described.

Two columns are printed. "corrected" is the current definition, which excludes
the no-op herald branches (HERALDED_ERASE heralds on its identity branch too, so
the DEM carries mechanisms that flip a herald detector and nothing else) from
both numerator and denominator. "legacy" counted them in both, inflating the
fraction; it is kept only to reproduce the pre-correction docs/AUDIT.md numbers.
Run:

    uv run python scripts/heralded_fraction.py                 # standard table
    uv run python scripts/heralded_fraction.py --p 0.02 --r_e 0.98 --convert-idle
"""

from __future__ import annotations

import argparse

from erasure_qec.analysis.dem_stats import dem_error_mass, heralded_fraction
from erasure_qec.config import NoiseParams


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--p", type=float, default=0.02)
    parser.add_argument("--d", type=int, default=5)
    parser.add_argument(
        "--r_e", type=float, default=None, help="single r_e; omit for the standard 0.5/0.98 table"
    )
    parser.add_argument("--convert-idle", action="store_true")
    args = parser.parse_args()

    if args.r_e is not None:
        params = NoiseParams(p=args.p, r_e=args.r_e, convert_idle=args.convert_idle)
        mass = dem_error_mass(params, d=args.d)
        print(f"p={args.p} d={args.d} r_e={args.r_e} convert_idle={args.convert_idle}")
        print(
            f"  DEM error mass: unheralded={mass.unheralded:.4f} "
            f"heralded={mass.heralded:.4f} no-op-herald={mass.no_op:.4f}"
        )
        print(
            f"  heralded_fraction (corrected, no-ops excluded) = "
            f"{heralded_fraction(params, d=args.d):.4f}"
        )
        print(
            f"  heralded_fraction (legacy, no-ops counted)     = "
            f"{heralded_fraction(params, d=args.d, count_no_op_heralds=True):.4f}"
        )
        return

    print(f"heralded fraction of the DEM error budget  (p={args.p}, d={args.d})")
    print("corrected = no-op herald branches excluded from numerator and denominator")
    print("legacy    = pre-correction definition, counted them in both (inflated)")
    print(f"{'':>6} {'convert_idle=False':>25} {'convert_idle=True':>25}")
    print(f"{'r_e':>6} {'corrected':>12} {'legacy':>12} {'corrected':>12} {'legacy':>12}")
    for r_e in (0.5, 0.98):
        row = []
        for convert_idle in (False, True):
            params = NoiseParams(p=args.p, r_e=r_e, convert_idle=convert_idle)
            row.append(heralded_fraction(params, d=args.d))
            row.append(heralded_fraction(params, d=args.d, count_no_op_heralds=True))
        print(f"{r_e:>6} " + " ".join(f"{v:>12.4f}" for v in row))
    print(
        "\nNote: heralded fraction < r_e always -- meas/reset errors are never "
        "heralded, so r_e (the 2q-gate fraction) overstates circuit-wide conversion."
    )


if __name__ == "__main__":
    main()
