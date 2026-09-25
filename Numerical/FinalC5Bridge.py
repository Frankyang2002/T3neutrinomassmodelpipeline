"""Compatibility surface for the historical final-C5 bridge module.

The canonical numerical adapter lives in ``Numerical.FinalC5TrajectoryAdapter``.
Thin wrappers are retained here because older tests and callers patch symbols on
this historical module path.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

import numpy as np

from Numerical.FinalC5TrajectoryAdapter import (
    FinalC5EvaluationInputs,
    FinalC5NumericalInputs,
    FinalC5TrajectoryEvaluator,
    HeavyMassBasisData,
    _heavy_mass_basis_data,
    extract_final_c5_inputs,
)
from Numerical.SMWeinbergStage import evaluate_final_weinberg_json
from Numerical.T3Trajectory import T3RenormalisableTrajectory


def final_c5_inputs_from_trajectory(
    trajectory: T3RenormalisableTrajectory,
    *,
    offdiagonal_mass_atol: float = 1.0e-10,
    imaginary_mass_atol: float = 1.0e-10,
    takagi_residual_rtol: float = 1.0e-10,
    takagi_residual_atol: float = 1.0e-10,
) -> FinalC5EvaluationInputs:
    """Historical wrapper for ``extract_final_c5_inputs``."""

    return extract_final_c5_inputs(
        trajectory,
        offdiagonal_mass_atol=offdiagonal_mass_atol,
        imaginary_mass_atol=imaginary_mass_atol,
        takagi_residual_rtol=takagi_residual_rtol,
        takagi_residual_atol=takagi_residual_atol,
    )


def evaluate_final_c5_from_trajectory(
    final_weinberg_path: Path,
    trajectory: T3RenormalisableTrajectory,
) -> np.ndarray:
    """Historical wrapper preserving patchability of the old module path."""

    inputs = final_c5_inputs_from_trajectory(trajectory)
    config = inputs.as_config()

    c5 = np.asarray(
        evaluate_final_weinberg_json(Path(final_weinberg_path), config),
        dtype=complex,
    )

    if c5.shape != (3, 3):
        raise RuntimeError(
            "Existing final-C5 evaluator returned a matrix with shape "
            f"{c5.shape}, expected (3, 3)."
        )

    if not np.all(np.isfinite(c5.real)) or not np.all(np.isfinite(c5.imag)):
        raise RuntimeError(
            "Existing final-C5 evaluator returned non-finite entries."
        )

    if not np.allclose(c5, c5.T, rtol=1.0e-10, atol=1.0e-14):
        raise RuntimeError(
            "Existing final-C5 evaluator returned a non-symmetric matrix."
        )

    return c5


@dataclass(frozen=True)
class FinalC5TrajectoryBuilder:
    """Historical callable adapter preserving old patch/import behavior."""

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


__all__ = [
    "FinalC5EvaluationInputs",
    "FinalC5NumericalInputs",
    "FinalC5TrajectoryBuilder",
    "FinalC5TrajectoryEvaluator",
    "HeavyMassBasisData",
    "_heavy_mass_basis_data",
    "evaluate_final_c5_from_trajectory",
    "evaluate_final_weinberg_json",
    "extract_final_c5_inputs",
    "final_c5_inputs_from_trajectory",
]
