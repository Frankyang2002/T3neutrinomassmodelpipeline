"""Numerical decomposition of the hierarchical final C5 coefficient.

The authoritative hierarchical matching result is

    C5_final = C5_hard + C5_direct

at the grouped scalar threshold.  ``FinalWeinbergAdapter`` already constructs
the symbolic physical Majorana matrices ``K_hard``, ``K_running`` and ``K``.
This module evaluates those same objects numerically on a
``T3RenormalisableTrajectory`` for plotting and diagnostics.

No new matching formula is introduced here.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import sympy as sp

from Numerical.FinalC5TrajectoryAdapter import extract_final_c5_inputs
from Numerical.SMWeinbergStage import (
    _hierarchical_scale_substitutions,
    _sympy_substitutions_from_config,
)
from Numerical.T3Trajectory import T3RenormalisableTrajectory
from RGE.matching.FinalWeinbergAdapter import (
    load_hierarchical_majorana_c5,
)


@dataclass(frozen=True)
class FinalC5ContributionBreakdown:
    """Numerical hard, direct-running, and combined physical C5 matrices."""

    hard: np.ndarray
    direct_running: np.ndarray
    combined: np.ndarray

    def validated(self) -> "FinalC5ContributionBreakdown":
        hard = np.asarray(self.hard, dtype=complex)
        direct = np.asarray(self.direct_running, dtype=complex)
        combined = np.asarray(self.combined, dtype=complex)

        for name, matrix in (
            ("hard", hard),
            ("direct_running", direct),
            ("combined", combined),
        ):
            if matrix.shape != (3, 3):
                raise ValueError(f"{name} must be a 3x3 matrix.")
            if not np.all(np.isfinite(matrix.real)) or not np.all(
                np.isfinite(matrix.imag)
            ):
                raise ValueError(f"{name} contains non-finite entries.")
            if not np.allclose(
                matrix,
                matrix.T,
                rtol=1.0e-10,
                atol=1.0e-14,
            ):
                raise ValueError(f"{name} must be symmetric.")

        if not np.allclose(
            hard + direct,
            combined,
            rtol=1.0e-9,
            atol=1.0e-14,
        ):
            raise ValueError(
                "Combined C5 does not equal hard + direct-running contributions."
            )

        return FinalC5ContributionBreakdown(
            hard=hard.copy(),
            direct_running=direct.copy(),
            combined=combined.copy(),
        )

    def as_json_dict(self) -> dict:
        result = self.validated()

        def block(matrix: np.ndarray) -> dict:
            return {
                "real": matrix.real.tolist(),
                "imag": matrix.imag.tolist(),
                "abs": np.abs(matrix).tolist(),
                "frobenius_norm": float(np.linalg.norm(matrix, ord="fro")),
            }

        hard_norm = float(np.linalg.norm(result.hard, ord="fro"))
        direct_norm = float(np.linalg.norm(result.direct_running, ord="fro"))

        return {
            "hard": block(result.hard),
            "direct_running": block(result.direct_running),
            "combined": block(result.combined),
            "direct_to_hard_norm_ratio": (
                None if hard_norm == 0.0 else direct_norm / hard_norm
            ),
        }


def _evaluate_symbolic_matrix(
    matrix: sp.Matrix,
    substitutions: dict[sp.Expr, complex | float],
) -> np.ndarray:
    values = np.empty((matrix.rows, matrix.cols), dtype=complex)

    for i in range(matrix.rows):
        for j in range(matrix.cols):
            value = sp.N(matrix[i, j].subs(substitutions), 18)

            if value.free_symbols:
                raise ValueError(
                    "Unresolved symbols remain in final-C5 contribution "
                    f"({i},{j}): {sorted(map(str, value.free_symbols))}"
                )

            values[i, j] = complex(value)

    return values


def evaluate_final_c5_contributions(
    final_weinberg_path: Path,
    trajectory: T3RenormalisableTrajectory,
) -> FinalC5ContributionBreakdown:
    """Evaluate the hard, direct, and combined C5 matrices at the scalar threshold."""

    inputs = extract_final_c5_inputs(trajectory)
    config = inputs.as_config()

    symbolic = load_hierarchical_majorana_c5(
        Path(final_weinberg_path),
        n_lepton=3,
        n_heavy=len(inputs.heavy_masses_gev),
        split_heavy_masses=True,
    )

    substitutions = _sympy_substitutions_from_config(config)
    substitutions.update(_hierarchical_scale_substitutions(config))

    result = FinalC5ContributionBreakdown(
        hard=_evaluate_symbolic_matrix(
            symbolic["K_hard"],
            substitutions,
        ),
        direct_running=_evaluate_symbolic_matrix(
            symbolic["K_running"],
            substitutions,
        ),
        combined=_evaluate_symbolic_matrix(
            symbolic["K"],
            substitutions,
        ),
    )

    return result.validated()
