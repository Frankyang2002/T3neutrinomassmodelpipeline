"""Numerical threshold helpers for the first T3 threshold.

The current symbolic pipeline defines EFT1 specifically as the theory after
integrating out the heavy fermion F:

    ordinary:      UV -> SM + S1 + S2
    shared scalar: UV -> SM + S

This module implements only the renormalisable state projection at that
boundary.  It does *not* replace the existing Matchete threshold matching for
higher-dimensional Wilson coefficients, and it does not invent finite
threshold corrections to renormalisable couplings.

At tree level the surviving renormalisable parameters are carried continuously
from the UV endpoint into EFT1, while F, MF and the F Yukawas disappear from
the active renormalisable state.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from Numerical.EFT1State import (
    SharedT3EFT1State,
    T3EFT1State,
)
from Numerical.State import (
    SharedT3UVState,
    T3UVState,
)


UVState = T3UVState | SharedT3UVState
EFT1State = T3EFT1State | SharedT3EFT1State


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


def fermion_singular_masses(state: UVState) -> FermionThresholdDiagnostic:
    """Return non-negative singular values of the running heavy mass matrix.

    For the Majorana branch these coincide with the Takagi singular masses.
    For the vector-like branch they are the usual bi-unitary singular masses.

    No automatic single matching scale is chosen here because the current
    field-level threshold plan integrates ``F`` as one group while a general
    3x3 MF can contain three distinct physical masses.
    """

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


def project_after_fermion_threshold(
    uv_state: UVState,
    *,
    matching_scale_gev: float | None = None,
) -> EFT1State:
    """Project the UV endpoint onto the renormalisable EFT1 state.

    ``matching_scale_gev`` is normally the same scale as ``uv_state.mu_gev``.
    A different value is rejected: the caller should first run the UV ODE to
    the desired matching scale, then project the endpoint.  This prevents an
    accidental discontinuous jump in scale without RGE evolution.
    """

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
            "scale before constructing EFT1."
        )

    if isinstance(uv_state, SharedT3UVState):
        return SharedT3EFT1State(
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

    return T3EFT1State(
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
