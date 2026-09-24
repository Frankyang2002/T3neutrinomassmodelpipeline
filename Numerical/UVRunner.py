"""Numerical one-loop UV running for the renormalisable T3 theory.

This is the first genuine UV ODE layer.  It evolves the complete canonical
``T3UVState`` or ``SharedT3UVState`` using RGBeta's exported one-loop beta
functions.

The independent variable is

    t = ln(mu),

and the RGBeta reporting convention is

    16*pi^2 dX/dt = beta_X^(1).

Unlike the older final-Weinberg numerical stage, this runner retains the full
integration trajectory rather than only the endpoint.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping

import numpy as np
from scipy.integrate import solve_ivp

from Numerical.RGBetaEvaluator import (
    derivative_from_rgbeta_payload,
    evaluate_rgbeta_payload,
)
from Numerical.State import SharedT3UVState, T3UVState
from Numerical.StateVector import (
    StateVectorLayout,
    pack_uv_state,
    unpack_uv_state,
)


UVState = T3UVState | SharedT3UVState


@dataclass(frozen=True)
class UVRunningPoint:
    """One saved point along a UV RGE trajectory."""

    mu_gev: float
    vector: np.ndarray


@dataclass(frozen=True)
class UVRunningResult:
    """Complete numerical result for one UV RGE segment."""

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
    def final_vector(self) -> np.ndarray:
        return self.y[:, -1].copy()

    @property
    def initial_vector(self) -> np.ndarray:
        return self.y[:, 0].copy()

    def state_at_index(self, index: int) -> UVState:
        """Reconstruct the validated state stored at one trajectory point."""

        return unpack_uv_state(
            self.y[:, index],
            self.layout,
            mu_gev=float(self.mu_gev[index]),
        )

    @property
    def final_state(self) -> UVState:
        return self.state_at_index(-1)

    @property
    def initial_state(self) -> UVState:
        return self.state_at_index(0)

    def points(self) -> tuple[UVRunningPoint, ...]:
        """Return immutable point objects for simple downstream iteration."""

        return tuple(
            UVRunningPoint(
                mu_gev=float(self.mu_gev[index]),
                vector=self.y[:, index].copy(),
            )
            for index in range(self.n_points)
        )


def _validate_scales(mu_initial_gev: float, mu_final_gev: float) -> tuple[float, float]:
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
) -> np.ndarray | None:
    """Convert optional physical save scales into increasing solver t values."""

    if save_scales_gev is None:
        return None

    scales = np.asarray(save_scales_gev, dtype=float)

    if scales.ndim != 1 or scales.size == 0:
        raise ValueError("save_scales_gev must be a non-empty one-dimensional sequence.")

    if not np.all(np.isfinite(scales)) or np.any(scales <= 0.0):
        raise ValueError("All save_scales_gev entries must be finite and positive.")

    lower = min(mu_initial_gev, mu_final_gev)
    upper = max(mu_initial_gev, mu_final_gev)

    if np.any(scales < lower) or np.any(scales > upper):
        raise ValueError(
            "Every save scale must lie inside the integration interval."
        )

    # solve_ivp requires t_eval to follow the integration direction.
    descending = mu_final_gev < mu_initial_gev
    scales = np.unique(scales)

    if descending:
        scales = scales[::-1]

    # Always retain exact segment endpoints.
    if scales[0] != mu_initial_gev:
        scales = np.insert(scales, 0, mu_initial_gev)

    if scales[-1] != mu_final_gev:
        scales = np.append(scales, mu_final_gev)

    return np.log(scales)


def _state_from_solver_vector(
    t: float,
    y: np.ndarray,
    layout: StateVectorLayout,
) -> UVState:
    return unpack_uv_state(
        y,
        layout,
        mu_gev=float(np.exp(t)),
    )


def _rhs_factory(
    payload: Mapping[str, Any],
    layout: StateVectorLayout,
):
    """Build ``dy/dt`` for a fixed RGBeta theory definition."""

    def rhs(t: float, y: np.ndarray) -> np.ndarray:
        state = _state_from_solver_vector(t, y, layout)

        # Evaluate against the same fixed RGBeta payload.  The payload contains
        # the representation-dependent beta equations; only the numerical
        # couplings vary with t.
        evaluated = evaluate_rgbeta_payload(payload, state)

        # derivative_from_rgbeta_payload would repack the state and return the
        # same layout.  The explicit conversion below avoids an unnecessary
        # second state pack in the hot ODE loop.
        from Numerical.BetaVector import beta_values_to_derivative

        return beta_values_to_derivative(
            evaluated,
            layout,
            divide_by_loop_factor=True,
        )

    return rhs


def run_uv_segment(
    initial_state: UVState,
    rgbeta_payload: Mapping[str, Any],
    mu_final_gev: float,
    *,
    save_scales_gev: np.ndarray | list[float] | tuple[float, ...] | None = None,
    rtol: float = 1.0e-8,
    atol: float = 1.0e-11,
    method: str = "DOP853",
    max_step_log: float = np.inf,
) -> UVRunningResult:
    """Evolve one complete UV T3 state between two positive scales.

    Parameters
    ----------
    initial_state:
        Canonical numerical UV state.  Its ``mu_gev`` is the initial scale.

    rgbeta_payload:
        Successful UV RGBeta JSON payload for the same representation.

    mu_final_gev:
        Final physical scale in GeV.  Upward and downward running are both
        supported.

    save_scales_gev:
        Optional physical scales at which the state must be stored.  Endpoints
        are inserted automatically.  If omitted, the adaptive solver's accepted
        points are retained.

    rtol, atol, method, max_step_log:
        Passed to ``scipy.integrate.solve_ivp``.  ``max_step_log`` is measured
        in the logarithmic variable t = ln(mu).

    Returns
    -------
    UVRunningResult
        Full saved trajectory and reconstructible endpoint state.
    """

    state = initial_state.validated()
    mu_initial, mu_final = _validate_scales(
        state.mu_gev,
        mu_final_gev,
    )

    y0, layout = pack_uv_state(state)

    # Fail before entering solve_ivp if the payload/state mapping is incomplete,
    # has the wrong representation metadata, or contains unsupported syntax.
    derivative0, derivative_layout = derivative_from_rgbeta_payload(
        rgbeta_payload,
        state,
    )

    if derivative_layout != layout:
        raise RuntimeError("RGBeta derivative layout does not match the UV state layout.")

    if derivative0.shape != y0.shape:
        raise RuntimeError("Initial RGBeta derivative does not match the state vector.")

    t_eval = _validate_save_scales(
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
            "UV numerical RGE integration failed: "
            + str(solution.message)
        )

    mu_values = np.exp(solution.t)

    # exp(log(mu)) introduces tiny floating-point endpoint drift at large
    # scales.  The physical segment endpoints are user inputs, so retain them
    # exactly in the public trajectory metadata.
    mu_values = np.asarray(mu_values, dtype=float)
    mu_values[0] = mu_initial
    mu_values[-1] = mu_final

    return UVRunningResult(
        mu_initial_gev=mu_initial,
        mu_final_gev=mu_final,
        layout=layout,
        mu_gev=np.asarray(mu_values, dtype=float),
        t=np.asarray(solution.t, dtype=float),
        y=np.asarray(solution.y, dtype=float),
        solver_success=bool(solution.success),
        solver_message=str(solution.message),
        nfev=int(solution.nfev),
        njev=int(getattr(solution, "njev", 0) or 0),
        nlu=int(getattr(solution, "nlu", 0) or 0),
    )
