"""Second-threshold helpers connecting EFT1 to the final SM+Weinberg EFT.

The current sequential T3 pipeline uses

    UV -> EFT1 after F -> final EFT after the scalar threshold.

For the ordinary branch, the second threshold integrates out S1 and S2 as one
group.  For the shared-scalar branch, it integrates out S.

This module handles only the numerical *renormalisable* state at that boundary.
It deliberately does not replace the existing Matchete/Wilson machinery that
constructs the matched Weinberg coefficient C5.

No automatic grouped scalar matching scale is chosen.  In the ordinary branch
mS1Sq and mS2Sq generally run differently, so a single grouped threshold scale
must remain an explicit pipeline choice.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from Numerical.EFT1State import (
    SharedT3EFT1State,
    T3EFT1State,
)
from Numerical.State import SMNumericalState
from RGE.running.weinberg.WeinbergRunning import SMInitialConditions


EFT1State = T3EFT1State | SharedT3EFT1State


@dataclass(frozen=True)
class ScalarThresholdDiagnostic:
    """Running scalar mass information at one EFT1 scale."""

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
    state: EFT1State,
) -> ScalarThresholdDiagnostic:
    """Return the physical running scalar masses at the current EFT1 scale.

    The masses are interpreted directly from the RGBeta quadratic parameters:

        MS1(mu) = sqrt(mS1Sq(mu))
        MS2(mu) = sqrt(mS2Sq(mu))

    or, in shared-scalar mode,

        MS(mu) = sqrt(mSSq(mu)).

    This is a diagnostic only; it does not select a matching scale.
    """

    state = state.validated()

    if isinstance(state, SharedT3EFT1State):
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
    eft1_state: EFT1State,
    *,
    matching_scale_gev: float | None = None,
) -> FinalSMBoundaryState:
    """Project the EFT1 endpoint onto the renormalisable final-SM state.

    Surviving SM couplings are carried continuously across this numerical
    projection.  Finite threshold effects and C5 matching are handled by the
    existing symbolic Matchete pipeline, not invented here.

    The EFT1 state must already have been evolved to the requested matching
    scale.
    """

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
            "The EFT1 state must first be evolved to the requested scalar "
            "matching scale before constructing the final SM boundary."
        )

    return FinalSMBoundaryState(
        mu_gev=matching_scale,
        sm=eft1_state.sm,
    ).validated()


def build_weinberg_initial_conditions(
    boundary: FinalSMBoundaryState,
    c5_matrix: np.ndarray,
) -> SMInitialConditions:
    """Build full-flavour final-SM+Weinberg numerical initial conditions.

    ``c5_matrix`` must already be the physical symmetric 3x3 coefficient
    produced by the established final-C5 matching stage.  This function only
    joins that coefficient to the running SM couplings at the scalar boundary.
    """

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
