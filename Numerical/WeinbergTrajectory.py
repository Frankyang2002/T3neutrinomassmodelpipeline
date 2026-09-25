"""Scale-dependent SM+Weinberg trajectory and neutrino observables.

This module reuses the numerical final-EFT state packing from
``Numerical.WeinbergRunning`` and retains the full accepted trajectory. The
actual beta-function model is defined under ``RGE.running.weinberg.WeinbergRGE``.

The project convention is

    L_EFT contains (1/2) C5 O5 + h.c.
    m_nu(mu) = -(v^2/2) C5(mu).

At every saved scale this module can therefore construct m_nu(mu) and run the
existing Takagi/observable machinery from
``RGE.phenomenology.NeutrinoObservables``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

from RGE.phenomenology.NeutrinoObservables import (
    NeutrinoObservables,
    calculate_neutrino_observables,
)
from Numerical.WeinbergRunning import (
    SMInitialConditions,
    _beta,
    _pack,
    _unpack,
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
class ScaleDependentNeutrinoPoint:
    """Neutrino mass matrix and observables at one physical scale."""

    mu_gev: float
    c5: np.ndarray
    mass_matrix_gev: np.ndarray
    observables: NeutrinoObservables


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
    initial: SMInitialConditions,
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

        cleaned = SMInitialConditions(
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


def neutrino_mass_from_c5(
    c5: np.ndarray,
    *,
    vev_gev: float = 246.22,
) -> np.ndarray:
    """Return m_nu = -(v^2/2) C5 in the project's physical convention."""

    c5 = np.asarray(c5, dtype=complex)

    if c5.shape != (3, 3):
        raise ValueError("C5 must be a 3x3 matrix.")

    if not np.allclose(
        c5,
        c5.T,
        rtol=1.0e-10,
        atol=1.0e-14,
    ):
        raise ValueError("C5 must be symmetric.")

    vev = float(vev_gev)

    if not np.isfinite(vev) or vev <= 0.0:
        raise ValueError("vev_gev must be finite and positive.")

    return -(vev**2 / 2.0) * c5


def charged_lepton_mass_basis_matrix(
    mass_matrix_gev: np.ndarray,
    ye: np.ndarray,
) -> np.ndarray:
    """Rotate a Majorana neutrino mass matrix to the charged-lepton mass basis."""

    mass_matrix = np.asarray(mass_matrix_gev, dtype=complex)
    ye = np.asarray(ye, dtype=complex)

    if mass_matrix.shape != (3, 3):
        raise ValueError("mass_matrix_gev must be a 3x3 matrix.")
    if ye.shape != (3, 3):
        raise ValueError("ye must be a 3x3 matrix.")

    left, singular_values, _ = np.linalg.svd(ye)
    order = np.argsort(singular_values)
    left = left[:, order]

    rotated = left.T @ mass_matrix @ left
    return 0.5 * (rotated + rotated.T)


def scale_dependent_neutrino_observables(
    trajectory: WeinbergTrajectoryResult,
    *,
    vev_gev: float = 246.22,
    ordering: str = "AUTO",
    max_takagi_residual: float = 1.0e-7,
) -> tuple[ScaleDependentNeutrinoPoint, ...]:
    """Evaluate m_nu(mu) and Takagi observables at every saved scale."""

    points: list[ScaleDependentNeutrinoPoint] = []

    for index in range(trajectory.n_points):
        running = trajectory.point_at_index(index)

        mass_matrix = neutrino_mass_from_c5(
            running.c5,
            vev_gev=vev_gev,
        )

        mass_matrix_flavor = charged_lepton_mass_basis_matrix(
            mass_matrix,
            running.ye,
        )

        observables = calculate_neutrino_observables(
            mass_matrix_flavor,
            ordering=ordering,
        )

        if observables.takagi_residual > max_takagi_residual:
            raise RuntimeError(
                "Takagi decomposition residual is too large at "
                f"mu={running.mu_gev} GeV: "
                f"{observables.takagi_residual}"
            )

        points.append(
            ScaleDependentNeutrinoPoint(
                mu_gev=running.mu_gev,
                c5=running.c5.copy(),
                mass_matrix_gev=mass_matrix,
                observables=observables,
            )
        )

    return tuple(points)
