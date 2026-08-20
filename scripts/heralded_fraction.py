"""Report the circuit-wide heralded fraction of the DEM error budget.

The honest "how much of the noise is erasure-converted" number for a config --
which is NOT r_e (measurement/reset stay unheralded; idle only with
convert_idle). Quote this wherever r_e is described. Run:

    uv run python scripts/heralded_fraction.py                 # standard table
    uv run python scripts/heralded_fraction.py --p 0.02 --r_e 0.98 --convert-idle
"""

from __future__ import annotations

import argparse

from erasure_qec.analysis.dem_stats import heralded_fraction
from erasure_qec.config import NoiseParams


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--p", type=float, default=0.02)
    parser.add_argument("--d", type=int, default=5)
    parser.add_argument("--r_e", type=float, default=None,
                        help="single r_e; omit for the standard 0.5/0.98 table")
    parser.add_argument("--convert-idle", action="store_true")
    args = parser.parse_args()

    if args.r_e is not None:
        hf = heralded_fraction(
            NoiseParams(p=args.p, r_e=args.r_e, convert_idle=args.convert_idle), d=args.d
        )
        print(f"p={args.p} d={args.d} r_e={args.r_e} convert_idle={args.convert_idle}: "
              f"heralded_fraction={hf:.4f}")
        return

    print(f"heralded fraction of total DEM error mass  (p={args.p}, d={args.d})")
    print(f"{'r_e':>6} {'convert_idle=False':>20} {'convert_idle=True':>20}")
    for r_e in (0.5, 0.98):
        off = heralded_fraction(NoiseParams(p=args.p, r_e=r_e), d=args.d)
        on = heralded_fraction(
            NoiseParams(p=args.p, r_e=r_e, convert_idle=True), d=args.d
        )
        print(f"{r_e:>6} {off:>20.4f} {on:>20.4f}")
    print("\nNote: heralded fraction < r_e always -- meas/reset errors are never "
          "heralded, so r_e (the 2q-gate fraction) overstates circuit-wide conversion.")


if __name__ == "__main__":
    main()
