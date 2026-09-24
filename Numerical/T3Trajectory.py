"""Multi-segment numerical trajectory controller for the T3 pipeline.

This module composes the numerical layers already implemented:

    UV state
      -> UV RGBeta running
      -> F-threshold projection
      -> EFT1 RGBeta running
      -> scalar-threshold projection
      -> final SM boundary
      -> optional existing SM+Weinberg running

It does not perform Matchete matching and does not invent C5.  The physical
symmetric C5 matrix must come from the existing matching pipeline and can then
be attached explicitly for the final low-energy Weinberg evolution.

The controller uses explicit matching scales.  It does not infer a grouped
fermion or scalar matching scale from running mass matrices.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping

import numpy as np

from Numerical.EFT1Runner import EFT1RunningResult, run_eft1_segment
from Numerical.EFT1State import SharedT3EFT1State, T3EFT1State
from Numerical.FinalSMBoundary import (
    FinalSMBoundaryState,
    ScalarThresholdDiagnostic,
    build_weinberg_initial_conditions,
    project_after_scalar_threshold,
    scalar_threshold_masses,
)
from Numerical.State import SharedT3UVState, T3UVState
from Numerical.ThresholdMatching import (
    FermionThresholdDiagnostic,
    fermion_singular_masses,
    project_after_fermion_threshold,
)
from Numerical.UVRunner import UVRunningResult, run_uv_segment
from RGE.running.weinberg.WeinbergRunning import (
    NumericalRGEResult,
    evolve_weinberg,
)


UVState = T3UVState | SharedT3UVState
EFT1State = T3EFT1State | SharedT3EFT1State


@dataclass(frozen=True)
class T3RenormalisableTrajectory:
    """Complete numerical path through the UV and EFT1 renormalisable stages."""

    mu_uv_initial_gev: float
    mu_fermion_threshold_gev: float
    mu_scalar_threshold_gev: float

    uv: UVRunningResult
    fermion_threshold: FermionThresholdDiagnostic
    eft1_initial_state: EFT1State
    eft1: EFT1RunningResult
    scalar_threshold: ScalarThresholdDiagnostic
    final_sm_boundary: FinalSMBoundaryState

    @property
    def uv_initial_state(self) -> UVState:
        return self.uv.initial_state

    @property
    def uv_threshold_state(self) -> UVState:
        return self.uv.final_state

    @property
    def eft1_threshold_state(self) -> EFT1State:
        return self.eft1.final_state


@dataclass(frozen=True)
class T3FullNumericalTrajectory:
    """Renormalisable T3 trajectory plus the existing SM+Weinberg segment."""

    renormalisable: T3RenormalisableTrajectory
    c5_at_scalar_threshold: np.ndarray
    weinberg: NumericalRGEResult


def _positive_finite_scale(name: str, value: float) -> float:
    result = float(value)

    if not np.isfinite(result) or result <= 0.0:
        raise ValueError(f"{name} must be finite and positive.")

    return result


def _validate_downward_threshold_order(
    mu_uv_initial_gev: float,
    mu_fermion_threshold_gev: float,
    mu_scalar_threshold_gev: float,
) -> tuple[float, float, float]:
    mu_uv = _positive_finite_scale(
        "initial UV scale",
        mu_uv_initial_gev,
    )
    mu_f = _positive_finite_scale(
        "fermion matching scale",
        mu_fermion_threshold_gev,
    )
    mu_s = _positive_finite_scale(
        "scalar matching scale",
        mu_scalar_threshold_gev,
    )

    if not (mu_uv > mu_f > mu_s):
        raise ValueError(
            "Sequential T3 running requires descending explicit scales "
            "mu_UV > mu_F > mu_S."
        )

    return mu_uv, mu_f, mu_s


def run_renormalisable_t3_trajectory(
    initial_state: UVState,
    uv_rgbeta_payload: Mapping[str, Any],
    eft1_rgbeta_payload: Mapping[str, Any],
    *,
    mu_fermion_threshold_gev: float,
    mu_scalar_threshold_gev: float,
    uv_save_scales_gev: (
        np.ndarray | list[float] | tuple[float, ...] | None
    ) = None,
    eft1_save_scales_gev: (
        np.ndarray | list[float] | tuple[float, ...] | None
    ) = None,
    rtol: float = 1.0e-8,
    atol: float = 1.0e-11,
    method: str = "DOP853",
    max_step_log: float = np.inf,
) -> T3RenormalisableTrajectory:
    """Run the complete renormalisable UV -> EFT1 -> SM boundary path.

    Matching scales are explicit inputs.  The running mass diagnostics are
    stored in the returned result so the caller can compare the chosen scales
    with the physical running masses.
    """

    initial_state = initial_state.validated()

    mu_uv, mu_f, mu_s = _validate_downward_threshold_order(
        initial_state.mu_gev,
        mu_fermion_threshold_gev,
        mu_scalar_threshold_gev,
    )

    uv = run_uv_segment(
        initial_state,
        uv_rgbeta_payload,
        mu_f,
        save_scales_gev=uv_save_scales_gev,
        rtol=rtol,
        atol=atol,
        method=method,
        max_step_log=max_step_log,
    )

    uv_threshold_state = uv.final_state
    fermion_diagnostic = fermion_singular_masses(
        uv_threshold_state
    )

    eft1_initial = project_after_fermion_threshold(
        uv_threshold_state,
        matching_scale_gev=mu_f,
    )

    eft1 = run_eft1_segment(
        eft1_initial,
        eft1_rgbeta_payload,
        mu_s,
        save_scales_gev=eft1_save_scales_gev,
        rtol=rtol,
        atol=atol,
        method=method,
        max_step_log=max_step_log,
    )

    eft1_threshold_state = eft1.final_state
    scalar_diagnostic = scalar_threshold_masses(
        eft1_threshold_state
    )

    final_sm_boundary = project_after_scalar_threshold(
        eft1_threshold_state,
        matching_scale_gev=mu_s,
    )

    return T3RenormalisableTrajectory(
        mu_uv_initial_gev=mu_uv,
        mu_fermion_threshold_gev=mu_f,
        mu_scalar_threshold_gev=mu_s,
        uv=uv,
        fermion_threshold=fermion_diagnostic,
        eft1_initial_state=eft1_initial,
        eft1=eft1,
        scalar_threshold=scalar_diagnostic,
        final_sm_boundary=final_sm_boundary,
    )


def continue_with_weinberg_running(
    trajectory: T3RenormalisableTrajectory,
    c5_matrix: np.ndarray,
    mu_low_gev: float,
    *,
    rtol: float = 1.0e-8,
    atol: float = 1.0e-11,
) -> T3FullNumericalTrajectory:
    """Attach matched C5 and run the existing final SM+Weinberg stage.

    The supplied ``c5_matrix`` is authoritative matching output from the
    existing C5 pipeline.  This function does not construct or modify it.
    """

    mu_low = _positive_finite_scale(
        "low-energy scale",
        mu_low_gev,
    )
    mu_matching = trajectory.mu_scalar_threshold_gev

    if mu_low >= mu_matching:
        raise ValueError(
            "Low-energy Weinberg running requires mu_low < mu_S."
        )

    initial = build_weinberg_initial_conditions(
        trajectory.final_sm_boundary,
        c5_matrix,
    )

    result = evolve_weinberg(
        initial,
        mu_matching,
        mu_low,
        rtol=rtol,
        atol=atol,
    )

    return T3FullNumericalTrajectory(
        renormalisable=trajectory,
        c5_at_scalar_threshold=np.asarray(
            c5_matrix,
            dtype=complex,
        ).copy(),
        weinberg=result,
    )
