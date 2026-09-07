"""The six committed worksheet SVGs must regenerate from their documented parameters.

``docs/dem_worksheet.md`` embeds six diagrams from ``docs/figures/``, but
``scripts/worksheet_probe.py`` writes only to the gitignored ``figures/dev/``:
its docstring says the committed copies are promoted **by hand**. That manual
step had no verification, and it silently failed — the two Probe C figures in
the tree were generated with ``mid_round_probe_erasures=[(2 + 2j, 0, 1)]``
(round **0**) while the worksheet prose, ``scripts/partition_check.py`` and
``tests/test_dem_partition.py`` all describe round **1**. An exhaustive search
over all 8 ancillas x 2 rounds x 4 CX layers matched the committed bytes at
exactly one point, ``(2 + 2j, 0, 1)``, so the mismatch was not ambiguous.
This test is the missing verification.

Comparison is deliberately not byte-for-byte: Stim's SVG coordinates can differ
in the last emitted decimal across platforms (e.g. 663.294 vs 663.293), so the
markup skeleton is compared exactly and the numbers are compared with a
tolerance. A trailing newline is also tolerated, because ``end-of-file-fixer``
appends one to any SVG it is allowed to touch. ``docs/figures/`` is excluded
from that hook in ``.pre-commit-config.yaml`` -- these files must stay exactly
what the generator emits -- but a contributor running without pre-commit
should get a clear diff about *content*, never about one whitespace byte.
"""

import re
from pathlib import Path

import pytest

from erasure_qec.circuits.builder import build
from erasure_qec.noise.injector import NullInjector

FIGURES = Path(__file__).resolve().parents[1] / "docs" / "figures"
PROBE_SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "worksheet_probe.py"

# Mirrors the PROBES table in scripts/worksheet_probe.py. The literals below are
# asserted to appear verbatim in that script by
# test_probe_parameters_match_the_generator_script, so the three-way tie between
# script, test and committed figure cannot drift silently.
_ProbeKwargs = dict[str, list[tuple[complex, int]] | list[tuple[complex, int, int]]]
PROBES: dict[str, tuple[str, _ProbeKwargs]] = {
    "probe_center_3_3": ("(3 + 3j, 0)", {"probe_erasures": [(3 + 3j, 0)]}),
    "probe_boundary_1_1": ("(1 + 1j, 0)", {"probe_erasures": [(1 + 1j, 0)]}),
    "probe_mid_round_ancilla_2_2": (
        "(2 + 2j, 1, 1)",
        {"mid_round_probe_erasures": [(2 + 2j, 1, 1)]},
    ),
}

DIAGRAMS = {"timeline": "timeline-svg", "detslice": "detslice-with-ops-svg"}

_FLOAT = re.compile(r"-?\d+\.\d+")

# Stim's SVG geometry is float-formatted to 3 decimals; the observed
# cross-platform drift is in that last digit, so 2e-3 covers it without
# masking a real coordinate change (the smallest structural shift seen when
# Probe C moves between rounds is many pixels).
FLOAT_TOL = 2e-3


def _skeleton_and_floats(text: str) -> tuple[str, list[float]]:
    """Split SVG text into (markup with every float blanked, the float values)."""
    floats = [float(m.group()) for m in _FLOAT.finditer(text)]
    return _FLOAT.sub("#", text.rstrip("\n")), floats


@pytest.mark.parametrize("stem", sorted(PROBES))
@pytest.mark.parametrize("kind", sorted(DIAGRAMS))
def test_committed_figure_regenerates_from_documented_parameters(stem: str, kind: str) -> None:
    _, kwargs = PROBES[stem]
    circuit = build(d=3, rounds=2, injector=NullInjector(), **kwargs)  # type: ignore[arg-type]
    generated = circuit.diagram(DIAGRAMS[kind])._repr_svg_()

    path = FIGURES / f"{stem}_{kind}.svg"
    assert path.is_file(), f"{path} is missing"
    committed = path.read_text()

    gen_skeleton, gen_floats = _skeleton_and_floats(generated)
    com_skeleton, com_floats = _skeleton_and_floats(committed)

    assert gen_skeleton == com_skeleton, (
        f"{path.name} does not regenerate from its documented parameters "
        f"({kwargs}); the committed copy is stale. Re-promote it from "
        f"figures/dev/ via scripts/worksheet_probe.py."
    )
    assert len(gen_floats) == len(com_floats)
    worst = max(
        (abs(a - b) for a, b in zip(gen_floats, com_floats, strict=True)),
        default=0.0,
    )
    assert worst <= FLOAT_TOL, f"{path.name}: largest coordinate drift {worst} exceeds {FLOAT_TOL}"


def test_probe_parameters_match_the_generator_script() -> None:
    """The parameters this test regenerates from are the ones the script uses."""
    source = PROBE_SCRIPT.read_text()
    for stem, (literal, _) in PROBES.items():
        assert literal in source, f"{stem}: {literal!r} not found in {PROBE_SCRIPT.name}"


def test_no_extra_or_missing_committed_figures() -> None:
    expected = {f"{stem}_{kind}.svg" for stem in PROBES for kind in DIAGRAMS}
    assert {p.name for p in FIGURES.glob("*.svg")} == expected
