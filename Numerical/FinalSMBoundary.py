"""Compatibility surface for the historical final-SM boundary module.

The canonical scalar-threshold boundary implementation now lives in
``Numerical.ScalarThresholdBoundary``.
"""

from Numerical.ScalarThresholdBoundary import (
    EFT1State,
    FinalSMBoundaryState,
    IntermediateScalarState,
    ScalarThresholdDiagnostic,
    build_final_sm_boundary,
    build_sm_weinberg_initial_conditions,
    build_weinberg_initial_conditions,
    project_after_scalar_threshold,
    scalar_threshold_masses,
)

__all__ = [
    "EFT1State",
    "FinalSMBoundaryState",
    "IntermediateScalarState",
    "ScalarThresholdDiagnostic",
    "build_final_sm_boundary",
    "build_sm_weinberg_initial_conditions",
    "build_weinberg_initial_conditions",
    "project_after_scalar_threshold",
    "scalar_threshold_masses",
]
