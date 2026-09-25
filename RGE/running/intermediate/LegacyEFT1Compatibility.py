"""Single compatibility boundary to the historical ``EFT1`` implementation.

The production pipeline identifies intermediate theories by physical field
content.  The original implementation modules predate that architecture and
encode the fermion-first scalar-only EFT using the ordinal name ``EFT1``.

Only this module should import those historical entry points.  Descriptive
production modules and validation backends depend on the aliases below instead,
so the old names can later be retired without changing the pipeline or the
serialized output contract.
"""

from __future__ import annotations

from RGE.running.eft1.EFT1DirectWeinberg import (
    build_direct_weinberg_flavor_transport as legacy_build_direct_weinberg_running,
    export_direct_weinberg_matchete as legacy_export_direct_weinberg_insertion,
)
from RGE.running.eft1.EFT1ThresholdResume import (
    Threshold2Continuation as LegacyScalarThresholdContinuation,
    rerun_threshold2_with_running as legacy_resume_scalar_threshold_with_running,
)
from RGE.running.eft1.EFT1WilsonFlow import (
    run_eft1_wilson_transport as legacy_component_wilson_transport,
    run_flavor_seed_export as legacy_export_full_flavor_wilson_boundary,
)
from RGE.running.eft1.EFT1WilsonRGE import (
    run_eft1_wilson_rge as legacy_scalar_only_dimension_five_wilson_rge,
)
from RGE.running.rgbeta.RGBetaT3Running import (
    run_rgbeta_t3_eft1 as legacy_scalar_only_renormalisable_rge,
)


__all__ = [
    "LegacyScalarThresholdContinuation",
    "legacy_build_direct_weinberg_running",
    "legacy_component_wilson_transport",
    "legacy_export_direct_weinberg_insertion",
    "legacy_export_full_flavor_wilson_boundary",
    "legacy_resume_scalar_threshold_with_running",
    "legacy_scalar_only_dimension_five_wilson_rge",
    "legacy_scalar_only_renormalisable_rge",
]
