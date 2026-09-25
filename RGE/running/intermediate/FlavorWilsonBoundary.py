"""Full-flavor Wilson boundary construction after the fermion threshold."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from RGE.running.intermediate.ScalarOnlyWilsonFlow import (
    export_full_flavor_wilson_boundary as _export_full_flavor_wilson_boundary,
)


def export_full_flavor_wilson_boundary(
    wilson_seed_path: Path,
    output_path: Path | None = None,
) -> dict[str, Any]:
    """Convert the matched one-generation LLSS seed to full flavor tensors."""
    return _export_full_flavor_wilson_boundary(
        wilson_seed_path=wilson_seed_path,
        output_path=output_path,
    )
