"""Numerical boundary construction at the heavy-fermion threshold.

For the verified fermion-first path the heavy fermion F is removed first:

    ordinary:      UV -> SM + S1 + S2
    shared scalar: UV -> SM + S

This module projects only the renormalisable numerical state. Higher-dimensional
Wilson matching remains in the Matchete/RGE matching pipeline, and no finite
threshold correction to renormalisable couplings is introduced here.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from Numerical.IntermediateScalarState import (
    SharedT3IntermediateScalarState,
    T3IntermediateScalarState,
)
from Numerical.State import (
    SharedT3UVState,
    T3UVState,
)


UVState = T3UVState | SharedT3UVState
IntermediateScalarState = T3IntermediateScalarState | SharedT3IntermediateScalarState


@dataclass(frozen=True)
class FermionThresholdDiagnostic:
    """Physical singular masses of MF at one running scale."""

    mu_gev: float
    masses_gev: np.ndarray

    @property
    def minimum_mass_gev(self) -> float:
        return float(np.min(self.masses_gev))

    @property
    def maximum_mass_gev(self) -> float:
        return float(np.max(self.masses_gev))


def fermion_threshold_masses(state: UVState) -> FermionThresholdDiagnostic:
    """Return non-negative singular values of the running heavy mass matrix."""

    state = state.validated()

    masses = np.linalg.svd(
        np.asarray(state.MF, dtype=complex),
        compute_uv=False,
    )
    masses = np.sort(np.asarray(masses, dtype=float))

    if not np.all(np.isfinite(masses)):
        raise ValueError("Heavy-fermion singular masses are non-finite.")

    if np.any(masses < 0.0):
        raise RuntimeError("Singular-value decomposition returned a negative mass.")

    return FermionThresholdDiagnostic(
        mu_gev=float(state.mu_gev),
        masses_gev=masses,
    )


def build_intermediate_scalar_boundary(
    uv_state: UVState,
    *,
    matching_scale_gev: float | None = None,
) -> IntermediateScalarState:
    """Build the scalar-only renormalisable boundary immediately below F."""

    uv_state = uv_state.validated()

    if matching_scale_gev is None:
        matching_scale = float(uv_state.mu_gev)
    else:
        matching_scale = float(matching_scale_gev)

    if not np.isfinite(matching_scale) or matching_scale <= 0.0:
        raise ValueError("matching_scale_gev must be finite and positive.")

    if not np.isclose(
        matching_scale,
        uv_state.mu_gev,
        rtol=1.0e-12,
        atol=0.0,
    ):
        raise ValueError(
            "The UV state must first be evolved to the requested matching "
            "scale before constructing the scalar-only intermediate EFT."
        )

    if isinstance(uv_state, SharedT3UVState):
        return SharedT3IntermediateScalarState(
            mu_gev=matching_scale,
            representation=uv_state.representation,
            sm=uv_state.sm,
            mSSq=uv_state.mSSq,
            lambdaS=uv_state.lambdaS,
            lambda3=uv_state.lambda3,
            lambda4=uv_state.lambda4,
            lambda5=uv_state.lambda5,
        ).validated()

    optional_names = (
        "lambdaH1Adj",
        "lambdaH2Adj",
        "lambdaS1Adj",
        "lambdaS2Adj",
        "lambda12Adj",
        "lambda12Cross",
        "lambdaHHdagS2S2",
        "lambdaHHdagS1barS1bar",
        "lambdaS1bar2S2bar2",
        "lambdaS1barS2S2bar2",
        "lambdaS1S1bar2S2bar",
        "lambdaHHdagS1barS2barCross",
    )

    optionals = {
        name: getattr(uv_state, name)
        for name in optional_names
    }

    return T3IntermediateScalarState(
        mu_gev=matching_scale,
        representation=uv_state.representation,
        sm=uv_state.sm,
        mS1Sq=uv_state.mS1Sq,
        mS2Sq=uv_state.mS2Sq,
        lambdaS1=uv_state.lambdaS1,
        lambdaS2=uv_state.lambdaS2,
        lambdaH1=uv_state.lambdaH1,
        lambdaH2=uv_state.lambdaH2,
        lambda12=uv_state.lambda12,
        lambdaT3=uv_state.lambdaT3,
        **optionals,
    ).validated()


# Historical function aliases.
fermion_singular_masses = fermion_threshold_masses
project_after_fermion_threshold = build_intermediate_scalar_boundary

# Compatibility type alias used by downstream parameter names/configuration.
EFT1State = IntermediateScalarState
