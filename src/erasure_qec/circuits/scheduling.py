"""Hook-safe CX schedules (PLAN.md §3.2).

Schedules are stored as per-basis tuples of slot indices into layout.py's
canonical [NE, NW, SE, SW] neighbor order — they say which plaquette
neighbor a CX layer touches for each basis, not raw coordinates.
"""

from erasure_qec.circuits.layout import Basis

# Aliases into the canonical [NE, NW, SE, SW] slot order.
NE, NW, SE, SW = 0, 1, 2, 3

Schedule = dict[Basis, tuple[int, int, int, int]]

# X-stabilizers touch neighbors in Z-shaped order (NE -> NW -> SE -> SW);
# Z-stabilizers touch neighbors in "n"-shaped order (NE -> SE -> NW -> SW).
# This orients every mid-window ancilla-fault hook (weight-2 data error)
# perpendicular to the logical operator it could shorten.
HOOK_SAFE_SCHEDULE: Schedule = {
    Basis.X: (NE, NW, SE, SW),
    Basis.Z: (NE, SE, NW, SW),
}

# Regression fixture (§4): both bases use the SAME order — the Z ("n") order
# NE -> SE -> NW -> SW. Forcing the X-stabilizers onto this order makes their
# last two CX targets (NW, SW) vertically aligned, so a single mid-window
# X-ancilla fault becomes a vertical weight-2 X error running PARALLEL to the
# logical X-chain of the Z-memory. That halves the effective distance:
# shortest_graphlike_error drops to ceil((d+1)/2) < d. Kept only to prove the
# hook-safe schedule above is load-bearing.
BROKEN_SCHEDULE: Schedule = {
    Basis.X: (NE, SE, NW, SW),
    Basis.Z: (NE, SE, NW, SW),
}
