"""Scalar-threshold helpers connecting the intermediate scalar EFT to SMEFT.

The verified sequential path is

    UV -> scalar-only EFT after F -> final EFT after the scalar threshold.

This module handles only the numerical renormalisable state. The established
Matchete/Wilson machinery remains authoritative for the matched Weinberg
coefficient C5.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from Numerical.IntermediateScalarState import (
    SharedT3IntermediateScalarState,
    T3IntermediateScalarState,
)
from Numerical.State import SMNumericalState
from Numerical.WeinbergRunning import SMInitialConditions


IntermediateScalarState = T3IntermediateScalarState | SharedT3IntermediateScalarState


@dataclass(frozen=True)
class ScalarThresholdDiagnostic:
    """Running scalar mass information at one intermediate-EFT scale."""

    mu_gev: float
    shared_scalar: bool
    mass_squared_gev2: tuple[float, ...]
    masses_gev: tuple[float, ...]

    @property
    def minimum_mass_gev(self) -> float:
        return float(min(self.masses_gev))

    @property
    def maximum_mass_gev(self) -> float:
        return float(max(self.masses_gev))

    @property
    def degenerate(self) -> bool:
        if len(self.masses_gev) < 2:
            return True
        return bool(
            np.allclose(
                self.masses_gev,
                self.masses_gev[0],
                rtol=1.0e-10,
                atol=0.0,
            )
        )


@dataclass(frozen=True)
class FinalSMBoundaryState:
    """Renormalisable SM state immediately below the scalar threshold."""

    mu_gev: float
    sm: SMNumericalState

    def validated(self) -> "FinalSMBoundaryState":
        mu = float(self.mu_gev)

        if not np.isfinite(mu) or mu <= 0.0:
            raise ValueError("mu_gev must be finite and positive.")

        return FinalSMBoundaryState(
            mu_gev=mu,
            sm=self.sm.validated(),
        )


def _positive_scalar_mass(
    name: str,
    mass_squared_gev2: float,
) -> float:
    value = float(mass_squared_gev2)

    if not np.isfinite(value):
        raise ValueError(f"{name} must be finite.")

    if value <= 0.0:
        raise ValueError(
            f"{name} must be positive to define a physical scalar threshold; "
            f"got {value} GeV^2."
        )

    return float(np.sqrt(value))


def scalar_threshold_masses(
    state: IntermediateScalarState,
) -> ScalarThresholdDiagnostic:
    """Return physical running scalar masses at the current threshold scale."""

    state = state.validated()

    if isinstance(state, SharedT3IntermediateScalarState):
        mass = _positive_scalar_mass("mSSq", state.mSSq)

        return ScalarThresholdDiagnostic(
            mu_gev=float(state.mu_gev),
            shared_scalar=True,
            mass_squared_gev2=(float(state.mSSq),),
            masses_gev=(mass,),
        )

    mass1 = _positive_scalar_mass("mS1Sq", state.mS1Sq)
    mass2 = _positive_scalar_mass("mS2Sq", state.mS2Sq)

    return ScalarThresholdDiagnostic(
        mu_gev=float(state.mu_gev),
        shared_scalar=False,
        mass_squared_gev2=(
            float(state.mS1Sq),
            float(state.mS2Sq),
        ),
        masses_gev=(mass1, mass2),
    )


def project_after_scalar_threshold(
    eft1_state: IntermediateScalarState,
    *,
    matching_scale_gev: float | None = None,
) -> FinalSMBoundaryState:
    """Project the scalar-only endpoint onto the renormalisable final SM state."""

    eft1_state = eft1_state.validated()

    if matching_scale_gev is None:
        matching_scale = float(eft1_state.mu_gev)
    else:
        matching_scale = float(matching_scale_gev)

    if not np.isfinite(matching_scale) or matching_scale <= 0.0:
        raise ValueError(
            "matching_scale_gev must be finite and positive."
        )

    if not np.isclose(
        matching_scale,
        eft1_state.mu_gev,
        rtol=1.0e-12,
        atol=0.0,
    ):
        raise ValueError(
            "The intermediate scalar state must first be evolved to the requested "
            "scalar matching scale before constructing the final SM boundary."
        )

    return FinalSMBoundaryState(
        mu_gev=matching_scale,
        sm=eft1_state.sm,
    ).validated()


def build_weinberg_initial_conditions(
    boundary: FinalSMBoundaryState,
    c5_matrix: np.ndarray,
) -> SMInitialConditions:
    """Join the authoritative matched C5 to the running SM boundary state."""

    boundary = boundary.validated()

    c5 = np.asarray(c5_matrix, dtype=complex)

    if c5.shape != (3, 3):
        raise ValueError("c5_matrix must be a 3x3 matrix.")

    if not np.all(np.isfinite(c5.real)) or not np.all(
        np.isfinite(c5.imag)
    ):
        raise ValueError("c5_matrix must contain only finite entries.")

    if not np.allclose(
        c5,
        c5.T,
        rtol=1.0e-10,
        atol=1.0e-14,
    ):
        raise ValueError(
            "c5_matrix must be symmetric for the Majorana Weinberg operator."
        )

    sm = boundary.sm

    return SMInitialConditions(
        gY=sm.gY,
        g2=sm.g2,
        g3=sm.g3,
        lambdaH=sm.lambdaH,
        ye=np.asarray(sm.ye, dtype=complex).copy(),
        yu=np.asarray(sm.yu, dtype=complex).copy(),
        yd=np.asarray(sm.yd, dtype=complex).copy(),
        K=c5,
    ).validated()


# Compatibility alias for public type annotations/configuration code.
EFT1State = IntermediateScalarState
