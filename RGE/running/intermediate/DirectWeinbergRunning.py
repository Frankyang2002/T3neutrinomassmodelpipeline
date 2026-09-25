"""Direct LLSS -> Weinberg running for the scalar-only T3 intermediate EFT."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from RGE.running.intermediate.LegacyEFT1Compatibility import (
    legacy_build_direct_weinberg_running,
    legacy_export_direct_weinberg_insertion,
)


def build_direct_weinberg_running(
    flavor_boundary_path: Path,
    component_rge_path: Path,
    mu_high: Any,
    mu_low: Any,
    output_path: Path | None = None,
) -> dict[str, Any]:
    """Build the verified full-flavor direct contribution to the Weinberg RGE."""
    return legacy_build_direct_weinberg_running(
        flavor_seed_path=flavor_boundary_path,
        component_rge_path=component_rge_path,
        mu_high=mu_high,
        mu_low=mu_low,
        output_path=output_path,
    )


def export_direct_weinberg_insertion(
    transport_path: Path,
    output_path: Path,
) -> dict[str, Any]:
    """Export the direct Weinberg running term for lower-threshold matching."""
    return legacy_export_direct_weinberg_insertion(
        transport_path=transport_path,
        output_path=output_path,
    )
