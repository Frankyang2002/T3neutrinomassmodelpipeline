"""Dimension-five Wilson running in the scalar-only T3 intermediate EFT."""

from __future__ import annotations

from pathlib import Path

from RGE.running.intermediate.LegacyEFT1Compatibility import (
    legacy_scalar_only_dimension_five_wilson_rge,
)


def run_scalar_only_dimension_five_wilson_rge(
    wilson_seed_path: Path,
    quartic_seed_path: Path,
    rgbeta_path: Path,
    output_path: Path | None = None,
    *,
    seed_support_only: bool = False,
) -> dict:
    """Run the verified one-loop ``psi^2 phi^2`` Wilson-tensor evolution.

    The numerical implementation and serialized payload are unchanged.  The
    historical ordinal-EFT entry point is isolated behind the compatibility
    boundary rather than leaking into production orchestration.
    """
    return legacy_scalar_only_dimension_five_wilson_rge(
        wilson_seed_path=wilson_seed_path,
        quartic_seed_path=quartic_seed_path,
        rgbeta_path=rgbeta_path,
        output_path=output_path,
        seed_support_only=seed_support_only,
    )
