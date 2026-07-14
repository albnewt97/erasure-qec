"""M0 scaffold gate: confirms the package imports and pytest collects/runs."""

import erasure_qec


def test_package_importable() -> None:
    assert erasure_qec is not None
