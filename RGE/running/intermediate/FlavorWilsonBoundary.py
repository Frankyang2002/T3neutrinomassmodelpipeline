"""Full-flavor Wilson boundary construction for the fermion-first T3 EFT."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from RGE.running.intermediate.LegacyEFT1Compatibility import (
    legacy_export_full_flavor_wilson_boundary,
)


def export_full_flavor_wilson_boundary(
    wilson_seed_path: Path,
    output_path: Path | None = None,
) -> dict[str, Any]:
    """Convert the matched one-generation LLSS seed to its flavor-tensor form."""
    return legacy_export_full_flavor_wilson_boundary(
        wilson_seed_path=wilson_seed_path,
        output_path=output_path,
    )
