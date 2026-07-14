"""Pure geometry for the rotated surface code (PLAN.md §2).

No stim imports here: this module only computes coordinates and adjacency so
it can be unit tested in isolation from circuit construction.
"""

from enum import Enum


class Basis(Enum):
    """Stabilizer basis measured by an ancilla qubit."""

    X = "X"
    Z = "Z"


# Canonical plaquette-neighbor slot order: index 0=NE, 1=NW, 2=SE, 3=SW.
# North is -y, south is +y, east is +x, west is -x.
_NEIGHBOR_OFFSETS: tuple[complex, complex, complex, complex] = (
    1 - 1j,  # NE
    -1 - 1j,  # NW
    1 + 1j,  # SE
    -1 + 1j,  # SW
)


def data_coords(d: int) -> set[complex]:
    """Coordinates of the d**2 data qubits: (2i+1, 2j+1) for 0 <= i, j < d."""
    if d < 1:
        raise ValueError(f"distance must be >= 1, got {d}")
    return {complex(2 * i + 1, 2 * j + 1) for i in range(d) for j in range(d)}


def ancilla_coords(d: int) -> dict[complex, Basis]:
    """Coordinates of the d**2 - 1 ancilla qubits, keyed to their stabilizer basis.

    Ancillas sit at even-coordinate plaquette centers (2a, 2b) for a, b in
    [0, d]. Interior centers (0 < a < d and 0 < b < d) always host one
    ancilla, checkerboarded by the parity of a + b (Z on even, X on odd).
    Boundary centers (on the edge of the bounding box) only survive as
    weight-2 checks when their checkerboard type matches the required
    boundary type: Z on the left/right edges, X on the top/bottom edges.
    The four bounding-box corners never host an ancilla.
    """
    if d < 1:
        raise ValueError(f"distance must be >= 1, got {d}")
    ancillas: dict[complex, Basis] = {}
    for a in range(d + 1):
        for b in range(d + 1):
            is_left_right = a in (0, d)
            is_top_bottom = b in (0, d)
            if is_left_right and is_top_bottom:
                continue  # bounding-box corner: never present
            basis = Basis.Z if (a + b) % 2 == 0 else Basis.X
            if is_left_right and basis is not Basis.Z:
                continue  # left/right boundary only hosts Z-type weight-2 checks
            if is_top_bottom and basis is not Basis.X:
                continue  # top/bottom boundary only hosts X-type weight-2 checks
            ancillas[complex(2 * a, 2 * b)] = basis
    return ancillas


def plaquette_neighbors(anc: complex, d: int) -> list[complex | None]:
    """Data-qubit neighbors of `anc` in the fixed canonical order [NE, NW, SE, SW].

    A slot is `None` when it falls off the data-qubit lattice: weight-2
    boundary stabilizers only have 2 of the 4 slots populated.
    """
    if anc not in ancilla_coords(d):
        raise ValueError(f"{anc!r} is not a valid ancilla coordinate for d={d}")
    valid = data_coords(d)
    neighbors: list[complex | None] = []
    for offset in _NEIGHBOR_OFFSETS:
        candidate = anc + offset
        neighbors.append(candidate if candidate in valid else None)
    return neighbors


def logical_z_support(d: int) -> list[complex]:
    """One horizontal row of data qubits, spanning the left/right Z boundaries."""
    if d < 1:
        raise ValueError(f"distance must be >= 1, got {d}")
    return [complex(2 * i + 1, 1) for i in range(d)]
