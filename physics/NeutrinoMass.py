"""Symbolic neutrino-mass construction from the Weinberg coefficient.

This module owns the symbolic physics conversion

    m_nu = -(v^2/2) C5

for the project's convention H^0=(v+h)/sqrt(2), together with the existing
symbolic output stage. It does not define or run any RGE.
"""

from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

from RGE.matching.FinalWeinbergAdapter import (
    is_hierarchical_final_c5,
    load_hierarchical_majorana_c5,
)
from RGE.matching.WeinbergFlavorMatching import (
    build_flavor_c5_from_matchete,
)


v = sp.Symbol("v")


def numerical_neutrino_mass_matrix(
    c5_matrix,
    *,
    vev_gev: float = 246.22,
):
    """Return numerical m_nu = -(v^2/2) C5 in the project convention."""

    import numpy as np

    c5 = np.asarray(c5_matrix, dtype=complex)

    if c5.shape != (3, 3):
        raise ValueError("C5 must be a 3x3 matrix.")

    vev = float(vev_gev)
    if not np.isfinite(vev) or vev <= 0.0:
        raise ValueError("vev_gev must be finite and positive.")

    return -(vev**2 / 2.0) * c5


def build_neutrino_mass_matrix(
    c5_matrix: sp.MatrixBase,
    *,
    vev: sp.Expr = v,
) -> sp.Matrix:
    """Return m_nu = -(v^2/2) C5 for H^0 = (v+h)/sqrt(2)."""

    return -sp.Rational(1, 2) * vev**2 * c5_matrix


def run_symbolic_neutrino_mass_stage(
    c5_path: Path,
    output_dir: Path,
    *,
    n_heavy: int = 3,
) -> dict:
    """Build and save the full symbolic neutrino-mass matrix."""

    c5_path = Path(c5_path)
    output_dir = Path(output_dir)
    data_dir = output_dir / "data"
    data_dir.mkdir(parents=True, exist_ok=True)

    if is_hierarchical_final_c5(c5_path):
        flavor_result = load_hierarchical_majorana_c5(
            c5_path,
            n_lepton=3,
            n_heavy=n_heavy,
            split_heavy_masses=True,
        )
        c5_input_kind = "final_weinberg_json"
    else:
        flavor_result = build_flavor_c5_from_matchete(
            c5_path,
            n_lepton=3,
            n_heavy=n_heavy,
            split_heavy_masses=True,
        )
        c5_input_kind = "legacy_scalar_c5"

    c5_matrix = flavor_result["K"]
    mass_matrix = build_neutrino_mass_matrix(c5_matrix)

    symmetry_difference = mass_matrix - mass_matrix.T
    if any(sp.simplify(entry) != 0 for entry in symmetry_difference):
        raise RuntimeError("Neutrino mass matrix is not symmetric.")

    mass_path = data_dir / "neutrino_mass_matrix.txt"
    mass_path.write_text(
        sp.sstr(mass_matrix) + "\n",
        encoding="utf-8",
    )

    summary = {
        "NeutrinoMassStatus": "Success",
        "NeutrinoMassConvention": "m_nu = -(v^2/2) C5",
        "C5InputKind": c5_input_kind,
        "C5InputFile": c5_path.relative_to(output_dir).as_posix()
        if c5_path.is_relative_to(output_dir)
        else str(c5_path),
        "NeutrinoMassMatrixFile": mass_path.relative_to(output_dir).as_posix(),
        "NeutrinoMassMatrixSymmetric": True,
    }

    summary_path = data_dir / "neutrino_mass_summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )

    return summary
