"""Physical interpretation of a saved SM+Weinberg trajectory.

This module owns the map from a numerical running C5 trajectory to physical
neutrino quantities:

    C5(mu) -> m_nu(mu) -> charged-lepton mass basis -> observables.

It deliberately contains no ODE integration and depends only on the
trajectory protocol defined below; any object implementing the existing
``n_points``/``point_at_index`` interface can be consumed.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

import numpy as np

from physics.NeutrinoMass import (
    numerical_neutrino_mass_matrix,
)
from physics.NeutrinoObservables import (
    NeutrinoObservables,
    calculate_neutrino_observables,
)


class _RunningPoint(Protocol):
    mu_gev: float
    ye: np.ndarray
    c5: np.ndarray


class _Trajectory(Protocol):
    @property
    def n_points(self) -> int: ...

    def point_at_index(self, index: int) -> _RunningPoint: ...


@dataclass(frozen=True)
class ScaleDependentNeutrinoPoint:
    """Neutrino mass matrix and observables at one physical scale."""

    mu_gev: float
    c5: np.ndarray
    mass_matrix_gev: np.ndarray
    observables: NeutrinoObservables


def neutrino_mass_from_c5(
    c5: np.ndarray,
    *,
    vev_gev: float = 246.22,
) -> np.ndarray:
    """Compatibility-named physics wrapper for C5 -> m_nu."""

    return numerical_neutrino_mass_matrix(
        c5,
        vev_gev=vev_gev,
    )


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
    trajectory: _Trajectory,
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
