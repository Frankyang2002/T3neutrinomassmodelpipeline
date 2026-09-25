"""Scale-dependent SM+Weinberg numerical trajectory.

This module owns numerical trajectory integration and saved running states.
The beta-function model is defined under ``RGE.running.weinberg.WeinbergRGE``.

Physical interpretation of the saved trajectory -- C5 -> m_nu, charged-lepton
basis rotation, and neutrino observables -- lives in
``physics.NeutrinoTrajectory``.  Those helpers are re-exported here for
compatibility with existing callers.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

from Numerical.SMWeinbergEvolution import (
    SMWeinbergInitialConditions,
    _beta,
    _pack,
    _unpack,
)
from physics.NeutrinoTrajectory import (
    ScaleDependentNeutrinoPoint,
    charged_lepton_mass_basis_matrix,
    neutrino_mass_from_c5,
    scale_dependent_neutrino_observables,
)


@dataclass(frozen=True)
class WeinbergTrajectoryPoint:
    """One saved SM+Weinberg point."""

    mu_gev: float
    gY: float
    g2: float
    g3: float
    lambdaH: float
    ye: np.ndarray
    yu: np.ndarray
    yd: np.ndarray
    c5: np.ndarray


@dataclass(frozen=True)
class WeinbergTrajectoryResult:
    """Full saved final-SM+Weinberg trajectory."""

    mu_initial_gev: float
    mu_final_gev: float
    t: np.ndarray
    y: np.ndarray
    mu_gev: np.ndarray
    solver_success: bool
    solver_message: str
    nfev: int

    @property
    def n_points(self) -> int:
        return int(self.mu_gev.size)

    def point_at_index(self, index: int) -> WeinbergTrajectoryPoint:
        (
            gY,
            g2,
            g3,
            lambdaH,
            ye,
            yu,
            yd,
            c5,
        ) = _unpack(self.y[:, index])

        return WeinbergTrajectoryPoint(
            mu_gev=float(self.mu_gev[index]),
            gY=float(gY),
            g2=float(g2),
            g3=float(g3),
            lambdaH=float(lambdaH),
            ye=np.asarray(ye, dtype=complex).copy(),
            yu=np.asarray(yu, dtype=complex).copy(),
            yd=np.asarray(yd, dtype=complex).copy(),
            c5=np.asarray(c5, dtype=complex).copy(),
        )

    @property
    def initial_point(self) -> WeinbergTrajectoryPoint:
        return self.point_at_index(0)

    @property
    def final_point(self) -> WeinbergTrajectoryPoint:
        return self.point_at_index(-1)


def _validate_scales(
    mu_initial_gev: float,
    mu_final_gev: float,
) -> tuple[float, float]:
    mu_initial = float(mu_initial_gev)
    mu_final = float(mu_final_gev)

    if not np.isfinite(mu_initial) or not np.isfinite(mu_final):
        raise ValueError("Weinberg trajectory scales must be finite.")

    if mu_initial <= 0.0 or mu_final <= 0.0:
        raise ValueError("Weinberg trajectory scales must be positive.")

    if mu_initial == mu_final:
        raise ValueError(
            "Initial and final Weinberg trajectory scales must be different."
        )

    return mu_initial, mu_final


def _prepare_save_scales(
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
            "Every save scale must lie inside the Weinberg integration interval."
        )

    scales = np.unique(scales)

    if mu_final_gev < mu_initial_gev:
        scales = scales[::-1]

    if scales[0] != mu_initial_gev:
        scales = np.insert(scales, 0, mu_initial_gev)

    if scales[-1] != mu_final_gev:
        scales = np.append(scales, mu_final_gev)

    return np.log(scales), scales.copy()


def run_weinberg_trajectory(
    initial: SMWeinbergInitialConditions,
    mu_initial_gev: float,
    mu_final_gev: float,
    *,
    save_scales_gev: np.ndarray | list[float] | tuple[float, ...] | None = None,
    rtol: float = 1.0e-8,
    atol: float = 1.0e-11,
    method: str = "DOP853",
    max_step_log: float = np.inf,
) -> WeinbergTrajectoryResult:
    """Run the existing SM+Weinberg equations while retaining the trajectory."""

    mu_initial, mu_final = _validate_scales(
        mu_initial_gev,
        mu_final_gev,
    )

    initial = initial.validated()
    y0 = _pack(initial)

    t_eval, exact_save_scales = _prepare_save_scales(
        save_scales_gev,
        mu_initial_gev=mu_initial,
        mu_final_gev=mu_final,
    )

    solution = solve_ivp(
        _beta,
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
            "Numerical Weinberg trajectory integration failed: "
            + str(solution.message)
        )

    if exact_save_scales is not None:
        mu_values = exact_save_scales.copy()
    else:
        mu_values = np.asarray(np.exp(solution.t), dtype=float)
        mu_values[0] = mu_initial
        mu_values[-1] = mu_final

    y = np.asarray(solution.y, dtype=float).copy()

    for index in range(y.shape[1]):
        (
            gY,
            g2,
            g3,
            lambdaH,
            ye,
            yu,
            yd,
            c5,
        ) = _unpack(y[:, index])

        c5 = 0.5 * (c5 + c5.T)

        cleaned = SMWeinbergInitialConditions(
            gY=gY,
            g2=g2,
            g3=g3,
            lambdaH=lambdaH,
            ye=ye,
            yu=yu,
            yd=yd,
            K=c5,
        ).validated()

        y[:, index] = _pack(cleaned)

    return WeinbergTrajectoryResult(
        mu_initial_gev=mu_initial,
        mu_final_gev=mu_final,
        t=np.asarray(solution.t, dtype=float),
        y=y,
        mu_gev=mu_values,
        solver_success=bool(solution.success),
        solver_message=str(solution.message),
        nfev=int(solution.nfev),
    )
