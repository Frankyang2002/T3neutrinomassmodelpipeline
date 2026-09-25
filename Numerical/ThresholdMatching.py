"""Compatibility surface for the historical threshold-projection module.

The canonical fermion-threshold boundary implementation now lives in
``Numerical.FermionThresholdBoundary``.
"""

from Numerical.FermionThresholdBoundary import (
    EFT1State,
    FermionThresholdDiagnostic,
    IntermediateScalarState,
    UVState,
    build_intermediate_scalar_boundary,
    fermion_singular_masses,
    fermion_threshold_masses,
    project_after_fermion_threshold,
)

__all__ = [
    "EFT1State",
    "FermionThresholdDiagnostic",
    "IntermediateScalarState",
    "UVState",
    "build_intermediate_scalar_boundary",
    "fermion_singular_masses",
    "fermion_threshold_masses",
    "project_after_fermion_threshold",
]
