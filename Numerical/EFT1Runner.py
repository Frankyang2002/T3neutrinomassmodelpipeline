"""Numerical one-loop running of the renormalisable T3 intermediate EFT.

EFT1 is the theory after the heavy fermion F has been integrated out:

    ordinary T3:      SM + S1 + S2
    shared-scalar T3: SM + S

The independent variable is t = ln(mu).  The RGBeta payload is expected to use

    16*pi^2 dX/dln(mu) = beta_X^(1),

and the full accepted trajectory is retained.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping

import numpy as np
from scipy.integrate import solve_ivp

from Numerical.BetaVector import beta_values_to_derivative
from Numerical.EFT1RGBetaEvaluator import (
    eft1_derivative_from_payload,
    evaluate_eft1_rgbeta_payload,
)
from Numerical.EFT1State import (
    SharedT3EFT1State,
    T3EFT1State,
)
from Numerical.EFT1StateVector import (
    EFT1State,
    pack_eft1_state,
    unpack_eft1_state,
)
from Numerical.StateVector import StateVectorLayout


@dataclass(frozen=True)
class EFT1RunningPoint:
    """One saved point on an EFT1 RGE trajectory."""

    mu_gev: float
    vector: np.ndarray


@dataclass(frozen=True)
class EFT1RunningResult:
    """Complete numerical result for one EFT1 RGE segment."""

    mu_initial_gev: float
    mu_final_gev: float
    layout: StateVectorLayout
    mu_gev: np.ndarray
    t: np.ndarray
    y: np.ndarray
    solver_success: bool
    solver_message: str
    nfev: int
    njev: int
    nlu: int

    @property
    def n_points(self) -> int:
        return int(self.mu_gev.size)

    @property
    def initial_vector(self) -> np.ndarray:
        return self.y[:, 0].copy()

    @property
    def final_vector(self) -> np.ndarray:
        return self.y[:, -1].copy()

    def state_at_index(self, index: int) -> EFT1State:
        return unpack_eft1_state(
            self.y[:, index],
            self.layout,
            mu_gev=float(self.mu_gev[index]),
        )

    @property
    def initial_state(self) -> EFT1State:
        return self.state_at_index(0)

    @property
    def final_state(self) -> EFT1State:
        return self.state_at_index(-1)

    def points(self) -> tuple[EFT1RunningPoint, ...]:
        return tuple(
            EFT1RunningPoint(
                mu_gev=float(self.mu_gev[index]),
                vector=self.y[:, index].copy(),
            )
            for index in range(self.n_points)
        )


def _validate_scales(
    mu_initial_gev: float,
    mu_final_gev: float,
) -> tuple[float, float]:
    mu_initial = float(mu_initial_gev)
    mu_final = float(mu_final_gev)

    if not np.isfinite(mu_initial) or not np.isfinite(mu_final):
        raise ValueError("RGE scales must be finite.")

    if mu_initial <= 0.0 or mu_final <= 0.0:
        raise ValueError("RGE scales must be positive.")

    if mu_initial == mu_final:
        raise ValueError("Initial and final RGE scales must be different.")

    return mu_initial, mu_final


def _validate_save_scales(
    save_scales_gev: np.ndarray | list[float] | tuple[float, ...] | None,
    *,
    mu_initial_gev: float,
    mu_final_gev: float,
) -> tuple[np.ndarray | None, np.ndarray | None]:
    if save_scales_gev is None:
        return None, None

    scales = np.asarray(save_scales_gev, dtype=float)

    if scales.ndim != 1 or scales.size == 0:
        raise ValueError(
            "save_scales_gev must be a non-empty one-dimensional sequence."
        )

    if not np.all(np.isfinite(scales)) or np.any(scales <= 0.0):
        raise ValueError(
            "All save_scales_gev entries must be finite and positive."
        )

    lower = min(mu_initial_gev, mu_final_gev)
    upper = max(mu_initial_gev, mu_final_gev)

    if np.any(scales < lower) or np.any(scales > upper):
        raise ValueError(
            "Every save scale must lie inside the integration interval."
        )

    scales = np.unique(scales)

    if mu_final_gev < mu_initial_gev:
        scales = scales[::-1]

    if scales[0] != mu_initial_gev:
        scales = np.insert(scales, 0, mu_initial_gev)

    if scales[-1] != mu_final_gev:
        scales = np.append(scales, mu_final_gev)

    return np.log(scales), scales.copy()


def _state_from_vector(
    t: float,
    y: np.ndarray,
    layout: StateVectorLayout,
) -> EFT1State:
    return unpack_eft1_state(
        y,
        layout,
        mu_gev=float(np.exp(t)),
    )


def _rhs_factory(
    payload: Mapping[str, Any],
    layout: StateVectorLayout,
):
    def rhs(t: float, y: np.ndarray) -> np.ndarray:
        state = _state_from_vector(t, y, layout)
        evaluated = evaluate_eft1_rgbeta_payload(payload, state)

        return beta_values_to_derivative(
            evaluated,
            layout,
            divide_by_loop_factor=True,
        )

    return rhs


def run_eft1_segment(
    initial_state: EFT1State,
    rgbeta_payload: Mapping[str, Any],
    mu_final_gev: float,
    *,
    save_scales_gev: np.ndarray | list[float] | tuple[float, ...] | None = None,
    rtol: float = 1.0e-8,
    atol: float = 1.0e-11,
    method: str = "DOP853",
    max_step_log: float = np.inf,
) -> EFT1RunningResult:
    """Evolve one renormalisable EFT1 state between two positive scales."""

    state = initial_state.validated()

    mu_initial, mu_final = _validate_scales(
        state.mu_gev,
        mu_final_gev,
    )

    y0, layout = pack_eft1_state(state)

    # Validate the payload/state mapping and the initial derivative before
    # entering solve_ivp.
    derivative0, derivative_layout = eft1_derivative_from_payload(
        rgbeta_payload,
        state,
    )

    if derivative_layout != layout:
        raise RuntimeError(
            "RGBeta EFT1 derivative layout does not match the state layout."
        )

    if derivative0.shape != y0.shape:
        raise RuntimeError(
            "Initial RGBeta EFT1 derivative does not match the state vector."
        )

    t_eval, exact_save_scales = _validate_save_scales(
        save_scales_gev,
        mu_initial_gev=mu_initial,
        mu_final_gev=mu_final,
    )

    solution = solve_ivp(
        _rhs_factory(rgbeta_payload, layout),
        t_span=(np.log(mu_initial), np.log(mu_final)),
        y0=y0,
        method=method,
        t_eval=t_eval,
        rtol=float(rtol),
        atol=float(atol),
        max_step=float(max_step_log),
    )

    if not solution.success:
        raise RuntimeError(
            "EFT1 numerical RGE integration failed: "
            + str(solution.message)
        )

    if exact_save_scales is not None:
        # When the user supplied physical save scales, retain those exact
        # values in the public trajectory instead of reconstructing them via
        # exp(log(mu)), which introduces avoidable floating-point drift.
        mu_values = exact_save_scales.copy()
    else:
        mu_values = np.asarray(np.exp(solution.t), dtype=float)
        mu_values[0] = mu_initial
        mu_values[-1] = mu_final

    return EFT1RunningResult(
        mu_initial_gev=mu_initial,
        mu_final_gev=mu_final,
        layout=layout,
        mu_gev=mu_values,
        t=np.asarray(solution.t, dtype=float),
        y=np.asarray(solution.y, dtype=float),
        solver_success=bool(solution.success),
        solver_message=str(solution.message),
        nfev=int(solution.nfev),
        njev=int(getattr(solution, "njev", 0) or 0),
        nlu=int(getattr(solution, "nlu", 0) or 0),
    )
