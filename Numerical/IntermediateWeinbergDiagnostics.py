"""Numerical sampling of the direct intermediate-EFT Weinberg contribution.

The authoritative hierarchical one-loop coefficient contains the direct

    LLSS -> O5

mixing accumulated between the fermion and scalar thresholds.  The final-C5
adapter stores that contribution symbolically in ``K_running``.  This module
samples the *same fixed-order expression* at intermediate scales by replacing
the lower matching scale ``MS`` with the requested running scale.

This is deliberately not the LLSS self-running correction.  LLSS self-running
is O(hbar) and feeding it through the scalar loop would first affect C5 at
O(hbar^2); the authoritative one-loop calculation excludes that insertion.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

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
class IntermediateWeinbergRunning:
    """Sampled direct LLSS -> Weinberg contribution."""

    mu_gev: np.ndarray
    delta_c5: np.ndarray

    @property
    def delta_c5_abs(self) -> np.ndarray:
        return np.abs(self.delta_c5)


def _evaluate_matrix(
    matrix: sp.Matrix,
    substitutions: dict[sp.Expr, complex | float],
) -> np.ndarray:
    result = np.empty(matrix.shape, dtype=complex)

    for i in range(matrix.rows):
        for j in range(matrix.cols):
            value = sp.N(matrix[i, j].subs(substitutions), 18)

            if value.free_symbols:
                raise ValueError(
                    "Unresolved symbols remain in direct intermediate "
                    f"Weinberg entry ({i},{j}): "
                    f"{sorted(map(str, value.free_symbols))}"
                )

            result[i, j] = complex(value)

    return result


def sample_intermediate_direct_weinberg(
    final_weinberg_path: Path,
    trajectory: T3RenormalisableTrajectory,
    scales_gev: Sequence[float],
) -> IntermediateWeinbergRunning:
    """Sample the implemented fixed-order direct C5 contribution.

    ``K_running`` in ``FinalWeinbergAdapter`` contains the direct one-loop
    intermediate-EFT contribution.  Its lower scale enters through ``MS``.
    Replacing ``MS`` by each requested ``mu`` therefore gives the accumulated
    contribution at that scale without introducing an additional RGE model.
    """

    scales = np.asarray(scales_gev, dtype=float)

    if scales.ndim != 1 or scales.size < 2:
        raise ValueError(
            "scales_gev must be a one-dimensional sequence with at least "
            "two entries."
        )
    if not np.all(np.isfinite(scales)) or np.any(scales <= 0.0):
        raise ValueError("All scales_gev entries must be finite and positive.")

    mu_f = float(trajectory.mu_fermion_threshold_gev)
    mu_s = float(trajectory.mu_scalar_threshold_gev)
    lower = min(mu_f, mu_s)
    upper = max(mu_f, mu_s)

    if np.any(scales < lower) or np.any(scales > upper):
        raise ValueError(
            "Every intermediate Weinberg save scale must lie between the "
            "fermion and scalar thresholds."
        )

    inputs = extract_final_c5_inputs(trajectory)
    config = inputs.as_config()

    symbolic = load_hierarchical_majorana_c5(
        Path(final_weinberg_path),
        n_lepton=3,
        n_heavy=len(inputs.heavy_masses_gev),
        split_heavy_masses=True,
    )
    k_running = symbolic["K_running"]

    substitutions = _sympy_substitutions_from_config(config)
    # Validate all non-running scale substitutions once.  The MS value is
    # overwritten below at every requested scale.
    substitutions.update(_hierarchical_scale_substitutions(config))
    ms_symbol = sp.Symbol("MS")

    matrices: list[np.ndarray] = []
    for mu in scales:
        point_substitutions = dict(substitutions)
        point_substitutions[ms_symbol] = float(mu)
        matrices.append(
            _evaluate_matrix(k_running, point_substitutions)
        )

    values = np.asarray(matrices, dtype=complex)

    if values.shape != (scales.size, 3, 3):
        raise RuntimeError(
            "Unexpected sampled direct-Weinberg trajectory shape: "
            f"{values.shape}."
        )

    if not np.allclose(
        values,
        np.swapaxes(values, 1, 2),
        rtol=1.0e-10,
        atol=1.0e-14,
    ):
        raise RuntimeError(
            "Sampled direct intermediate Weinberg matrices are not symmetric."
        )

    return IntermediateWeinbergRunning(
        mu_gev=scales.copy(),
        delta_c5=values,
    )
