"""Geometry invariants for erasure_qec.circuits.layout (PLAN.md §2)."""

import pytest

from erasure_qec.circuits.layout import (
    Basis,
    ancilla_coords,
    data_coords,
    logical_z_support,
    plaquette_neighbors,
)

DISTANCES = [3, 5, 7]


def _neighbor_qubits(anc: complex, d: int) -> set[complex]:
    return {q for q in plaquette_neighbors(anc, d) if q is not None}


@pytest.mark.parametrize("d", DISTANCES)
def test_qubit_counts(d: int) -> None:
    assert len(data_coords(d)) == d**2
    assert len(ancilla_coords(d)) == d**2 - 1


@pytest.mark.parametrize("d", DISTANCES)
def test_every_data_qubit_touched_by_at_most_four_ancillas(d: int) -> None:
    touch_count: dict[complex, int] = dict.fromkeys(data_coords(d), 0)
    for anc in ancilla_coords(d):
        for q in _neighbor_qubits(anc, d):
            touch_count[q] += 1
    assert all(count <= 4 for count in touch_count.values())


@pytest.mark.parametrize("d", DISTANCES)
def test_x_and_z_plaquettes_share_zero_or_two_data_qubits(d: int) -> None:
    ancillas = ancilla_coords(d)
    x_ancs = [a for a, basis in ancillas.items() if basis is Basis.X]
    z_ancs = [a for a, basis in ancillas.items() if basis is Basis.Z]
    for x_anc in x_ancs:
        x_qubits = _neighbor_qubits(x_anc, d)
        for z_anc in z_ancs:
            shared = x_qubits & _neighbor_qubits(z_anc, d)
            assert len(shared) in (0, 2), (x_anc, z_anc, shared)


@pytest.mark.parametrize("d", DISTANCES)
def test_boundary_weight2_checks_appear_only_on_correct_edges(d: int) -> None:
    for anc, basis in ancilla_coords(d).items():
        degree = len(_neighbor_qubits(anc, d))
        on_left_right = anc.real in (0, 2 * d)
        on_top_bottom = anc.imag in (0, 2 * d)
        if degree == 4:
            assert not on_left_right and not on_top_bottom
            continue
        assert degree == 2, (anc, degree)
        assert on_left_right != on_top_bottom  # boundary check, not an interior/corner point
        if on_left_right:
            assert basis is Basis.Z, f"left/right boundary check {anc} must be Z-type"
        if on_top_bottom:
            assert basis is Basis.X, f"top/bottom boundary check {anc} must be X-type"


@pytest.mark.parametrize("d", DISTANCES)
def test_logical_z_support_is_one_horizontal_row(d: int) -> None:
    support = logical_z_support(d)
    assert len(support) == d
    assert set(support).issubset(data_coords(d))
    assert {q.imag for q in support} == {1.0}
    assert {q.real for q in support} == {2 * i + 1 for i in range(d)}


@pytest.mark.parametrize("d", DISTANCES)
def test_plaquette_neighbors_order_and_none_slots(d: int) -> None:
    valid = data_coords(d)
    for anc in ancilla_coords(d):
        neighbors = plaquette_neighbors(anc, d)
        assert len(neighbors) == 4
        for slot in neighbors:
            assert slot is None or slot in valid
        ne, nw, se, sw = neighbors
        if ne is not None and se is not None:
            assert ne.real == se.real  # NE, SE share the eastern x
        if nw is not None and sw is not None:
            assert nw.real == sw.real  # NW, SW share the western x
        if ne is not None and nw is not None:
            assert ne.imag == nw.imag  # NE, NW share the northern y


def test_plaquette_neighbors_rejects_non_ancilla_coordinate() -> None:
    with pytest.raises(ValueError):
        plaquette_neighbors(complex(0, 0), 3)
