"""Bridge from a numerical T3 trajectory to the existing final-C5 evaluator.

The authoritative final hierarchical Weinberg coefficient remains the JSON
produced by

    RGE.running.weinberg.FinalWeinbergCoefficient

and is evaluated numerically by

    RGE.running.weinberg.NumericalWeinbergStage.evaluate_final_weinberg_json.

This module supplies the missing adapter from the running numerical trajectory
to the configuration expected by that existing evaluator.

Important basis convention
--------------------------
``RGE.matching.FlavorC5Matching`` evaluates the heavy-generation sum in the
heavy-fermion mass basis.

For a Majorana heavy fermion the UV state enforces a complex-symmetric MF.
At the F threshold this bridge therefore performs the Takagi change of basis

    U_F^T MF U_F = diag(M_1, M_2, M_3),  M_r >= 0,

and rotates the heavy index of the lepton-by-heavy Yukawa matrices as

    y1 -> y1 U_F,
    y2 -> y2 U_F.

This follows the matrix layout used by ``FlavorC5Matching``, where the heavy
index is the Yukawa column index.

For non-Majorana heavy fermions, an off-diagonal MF would require a validated
bi-unitary Dirac/vectorlike basis transformation.  That transformation is not
implemented here, so such states are still rejected.

Shared-scalar mode is also rejected here: the current final-C5 numerical
adapter is written in terms of the ordinary y1/y2, MS1/MS2, lambdaT3
parameterization, whereas the shared numerical state uses h and lambda5.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

import numpy as np

from Numerical.T3Trajectory import T3RenormalisableTrajectory
from Numerical.State import SharedT3UVState, T3UVState
from Numerical.EFT1State import SharedT3EFT1State, T3EFT1State
from RGE.running.weinberg.NumericalWeinbergStage import (
    evaluate_final_weinberg_json,
)


@dataclass(frozen=True)
class FinalC5NumericalInputs:
    """Numerical quantities supplied to the existing final-C5 evaluator."""

    heavy_masses_gev: np.ndarray
    y1: np.ndarray
    y2: np.ndarray
    scalar_mass1_gev: float
    scalar_mass2_gev: float
    grouped_scalar_matching_scale_gev: float
    lambda_t3: complex

    def as_config(self) -> dict:
        """Return the config subset expected by evaluate_final_weinberg_json."""

        return {
            "t3": {
                "MF": self.heavy_masses_gev.tolist(),
                "MS1": float(self.scalar_mass1_gev),
                "MS2": float(self.scalar_mass2_gev),
                # The hierarchical final-C5 adapter uses MS explicitly when
                # supplied, avoiding the old equal-MS1/MS2 fallback.
                "MS": float(self.grouped_scalar_matching_scale_gev),
                "lambdaT3": complex(self.lambda_t3),
                "y1_real": self.y1.real.tolist(),
                "y1_imag": self.y1.imag.tolist(),
                "y2_real": self.y2.real.tolist(),
                "y2_imag": self.y2.imag.tolist(),
                "hbar": 1.0 / (16.0 * np.pi**2),
            }
        }


@dataclass(frozen=True)
class HeavyMassBasisData:
    """Heavy masses and Yukawas expressed in the validated mass basis."""

    masses_gev: np.ndarray
    y1: np.ndarray
    y2: np.ndarray
    rotation: np.ndarray
    takagi_residual: float


def _already_diagonal_positive(
    matrix: np.ndarray,
    *,
    offdiagonal_atol: float,
    imaginary_atol: float,
) -> np.ndarray | None:
    """Return positive real diagonal entries when no basis change is needed."""

    diagonal_matrix = np.diag(np.diag(matrix))

    if not np.allclose(
        matrix,
        diagonal_matrix,
        rtol=0.0,
        atol=float(offdiagonal_atol),
    ):
        return None

    diagonal = np.diag(matrix)

    if not np.allclose(
        diagonal.imag,
        0.0,
        rtol=0.0,
        atol=float(imaginary_atol),
    ):
        return None

    masses = np.asarray(diagonal.real, dtype=float)

    if not np.all(np.isfinite(masses)):
        raise ValueError("Heavy-fermion masses must be finite.")

    if np.any(masses <= 0.0):
        return None

    return masses


def _majorana_takagi_mass_basis(
    state: T3UVState,
    *,
    symmetry_atol: float,
    residual_rtol: float,
    residual_atol: float,
) -> HeavyMassBasisData:
    """Takagi-diagonalize Majorana MF and rotate the heavy Yukawa index.

    The convention is

        U_F^T MF U_F = D,

    with D real, non-negative and sorted in ascending order.  Since y1 and y2
    have shape (lepton, heavy), the same field redefinition acts on their
    columns as y -> y U_F.
    """

    matrix = np.asarray(state.MF, dtype=complex)

    if matrix.shape != (3, 3):
        raise ValueError("MF must be a 3x3 matrix.")

    if not np.all(np.isfinite(matrix.real)) or not np.all(
        np.isfinite(matrix.imag)
    ):
        raise ValueError("MF must contain only finite entries.")

    scale = max(1.0, float(np.linalg.norm(matrix, ord=2)))
    symmetry_tolerance = max(
        float(symmetry_atol),
        1.0e-12 * scale,
    )

    if not np.allclose(
        matrix,
        matrix.T,
        rtol=0.0,
        atol=symmetry_tolerance,
    ):
        raise ValueError(
            "Majorana MF must be complex symmetric before Takagi "
            "diagonalization."
        )

    # Remove only numerical antisymmetric roundoff before factorization.
    matrix = 0.5 * (matrix + matrix.T)

    # For a complex-symmetric matrix M, the right singular vectors provide a
    # Takagi basis up to column phases when the singular values are resolved.
    # We explicitly phase-fix U so U^T M U has positive real diagonal entries.
    _, _, vh = np.linalg.svd(matrix)
    rotation = vh.conj().T

    takagi_matrix = rotation.T @ matrix @ rotation
    diagonal = np.diag(takagi_matrix)

    phase = np.ones(3, dtype=complex)
    nonzero = np.abs(diagonal) > max(float(residual_atol), 1.0e-30)
    phase[nonzero] = np.exp(
        -0.5j * np.angle(diagonal[nonzero])
    )
    rotation = rotation @ np.diag(phase)

    takagi_matrix = rotation.T @ matrix @ rotation

    # Deterministic ordering: pair each Yukawa column with the corresponding
    # ascending positive heavy mass.
    masses = np.asarray(
        np.real(np.diag(takagi_matrix)),
        dtype=float,
    )
    order = np.argsort(masses)
    rotation = rotation[:, order]
    takagi_matrix = rotation.T @ matrix @ rotation
    masses = np.asarray(
        np.real(np.diag(takagi_matrix)),
        dtype=float,
    )

    offdiagonal = takagi_matrix - np.diag(np.diag(takagi_matrix))
    imaginary_diagonal = np.imag(np.diag(takagi_matrix))

    residual = max(
        float(np.linalg.norm(offdiagonal, ord="fro")),
        float(np.linalg.norm(imaginary_diagonal)),
    )

    allowed = float(residual_atol) + float(residual_rtol) * max(
        1.0,
        float(np.linalg.norm(matrix, ord="fro")),
    )

    if residual > allowed:
        raise RuntimeError(
            "Takagi diagonalization of Majorana MF did not reach the "
            f"requested residual tolerance: residual={residual:.6e}, "
            f"allowed={allowed:.6e}. Near-degenerate heavy masses may require "
            "a dedicated degenerate-subspace treatment."
        )

    if not np.all(np.isfinite(masses)):
        raise ValueError("Takagi heavy masses must be finite.")

    # Tiny negative values can arise from roundoff after the phase fix.
    negative_tolerance = float(residual_atol) + float(residual_rtol) * max(
        1.0,
        float(np.max(np.abs(masses))),
    )
    if np.any(masses < -negative_tolerance):
        raise RuntimeError(
            "Takagi diagonalization produced a negative heavy mass beyond "
            "numerical tolerance."
        )
    masses = np.maximum(masses, 0.0)

    if np.any(masses <= 0.0):
        raise ValueError(
            "Final-C5 numerical evaluation requires strictly positive "
            "Majorana heavy masses."
        )

    y1 = np.asarray(state.y1, dtype=complex) @ rotation
    y2 = np.asarray(state.y2, dtype=complex) @ rotation

    return HeavyMassBasisData(
        masses_gev=masses.copy(),
        y1=y1.copy(),
        y2=y2.copy(),
        rotation=rotation.copy(),
        takagi_residual=residual,
    )


def _heavy_mass_basis_data(
    state: T3UVState,
    *,
    offdiagonal_atol: float,
    imaginary_atol: float,
    takagi_residual_rtol: float,
    takagi_residual_atol: float,
) -> HeavyMassBasisData:
    """Return heavy masses and Yukawas in the basis expected by final C5."""

    matrix = np.asarray(state.MF, dtype=complex)

    if matrix.shape != (3, 3):
        raise ValueError("MF must be a 3x3 matrix.")

    diagonal_positive = _already_diagonal_positive(
        matrix,
        offdiagonal_atol=offdiagonal_atol,
        imaginary_atol=imaginary_atol,
    )

    if diagonal_positive is not None:
        identity = np.eye(3, dtype=complex)
        return HeavyMassBasisData(
            masses_gev=diagonal_positive.copy(),
            y1=np.asarray(state.y1, dtype=complex).copy(),
            y2=np.asarray(state.y2, dtype=complex).copy(),
            rotation=identity,
            takagi_residual=0.0,
        )

    if state.representation.majorana_fermion:
        return _majorana_takagi_mass_basis(
            state,
            symmetry_atol=offdiagonal_atol,
            residual_rtol=takagi_residual_rtol,
            residual_atol=takagi_residual_atol,
        )

    raise ValueError(
        "Final-C5 numerical evaluation requires a diagonal positive heavy-"
        "fermion mass basis for non-Majorana F. MF is off-diagonal, complex, "
        "or has non-positive diagonal entries; no validated bi-unitary "
        "Dirac/vectorlike basis rotation is implemented."
    )


def final_c5_inputs_from_trajectory(
    trajectory: T3RenormalisableTrajectory,
    *,
    offdiagonal_mass_atol: float = 1.0e-10,
    imaginary_mass_atol: float = 1.0e-10,
    takagi_residual_rtol: float = 1.0e-10,
    takagi_residual_atol: float = 1.0e-10,
) -> FinalC5NumericalInputs:
    """Extract the existing final-C5 evaluator inputs from a trajectory.

    Scale assignment
    ----------------
    y1, y2, MF:
        taken at the UV endpoint mu_F, immediately before F is removed.

    MS1, MS2, lambdaT3:
        taken at the EFT1 endpoint mu_S, immediately before S1/S2 are removed.

    MS:
        the explicit grouped scalar matching scale used by the numerical
        trajectory.
    """

    uv_state = trajectory.uv_threshold_state
    eft1_state = trajectory.eft1_threshold_state

    if isinstance(uv_state, SharedT3UVState) or isinstance(
        eft1_state,
        SharedT3EFT1State,
    ):
        raise NotImplementedError(
            "The current final-C5 numerical evaluator uses the ordinary "
            "y1/y2, MS1/MS2, lambdaT3 parameterization. A separate validated "
            "shared-scalar h/lambda5 adapter is required."
        )

    if not isinstance(uv_state, T3UVState):
        raise TypeError("Expected an ordinary T3UVState at the F threshold.")

    if not isinstance(eft1_state, T3EFT1State):
        raise TypeError("Expected an ordinary T3EFT1State at the scalar threshold.")

    heavy_basis = _heavy_mass_basis_data(
        uv_state,
        offdiagonal_atol=offdiagonal_mass_atol,
        imaginary_atol=imaginary_mass_atol,
        takagi_residual_rtol=takagi_residual_rtol,
        takagi_residual_atol=takagi_residual_atol,
    )

    if eft1_state.mS1Sq <= 0.0 or eft1_state.mS2Sq <= 0.0:
        raise ValueError(
            "Positive scalar mass-squared parameters are required for final-C5 "
            "numerical evaluation."
        )

    scalar_mass1 = float(np.sqrt(eft1_state.mS1Sq))
    scalar_mass2 = float(np.sqrt(eft1_state.mS2Sq))

    if not np.isfinite(scalar_mass1) or not np.isfinite(scalar_mass2):
        raise ValueError("Running scalar masses are non-finite.")

    return FinalC5NumericalInputs(
        heavy_masses_gev=heavy_basis.masses_gev.copy(),
        y1=heavy_basis.y1.copy(),
        y2=heavy_basis.y2.copy(),
        scalar_mass1_gev=scalar_mass1,
        scalar_mass2_gev=scalar_mass2,
        grouped_scalar_matching_scale_gev=float(
            trajectory.mu_scalar_threshold_gev
        ),
        lambda_t3=complex(eft1_state.lambdaT3),
    )


def evaluate_final_c5_from_trajectory(
    final_weinberg_path: Path,
    trajectory: T3RenormalisableTrajectory,
) -> np.ndarray:
    """Evaluate the authoritative final-C5 JSON on a numerical trajectory."""

    inputs = final_c5_inputs_from_trajectory(trajectory)
    config = inputs.as_config()

    c5 = np.asarray(
        evaluate_final_weinberg_json(
            Path(final_weinberg_path),
            config,
        ),
        dtype=complex,
    )

    if c5.shape != (3, 3):
        raise RuntimeError(
            "Existing final-C5 evaluator returned a matrix with shape "
            f"{c5.shape}, expected (3, 3)."
        )

    if not np.all(np.isfinite(c5.real)) or not np.all(
        np.isfinite(c5.imag)
    ):
        raise RuntimeError(
            "Existing final-C5 evaluator returned non-finite entries."
        )

    if not np.allclose(
        c5,
        c5.T,
        rtol=1.0e-10,
        atol=1.0e-14,
    ):
        raise RuntimeError(
            "Existing final-C5 evaluator returned a non-symmetric matrix."
        )

    return c5


@dataclass(frozen=True)
class FinalC5TrajectoryBuilder:
    """Callable adapter compatible with ParameterScan.make_t3_point_evaluator."""

    final_weinberg_path: Path

    def __call__(
        self,
        _parameters: Mapping[str, float],
        trajectory: T3RenormalisableTrajectory,
    ) -> np.ndarray:
        return evaluate_final_c5_from_trajectory(
            self.final_weinberg_path,
            trajectory,
        )
