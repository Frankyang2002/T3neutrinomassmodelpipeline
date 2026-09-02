from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

from RGE.matching.FlavorMatchedC5 import flavor_match_from_c5_file


v = sp.Symbol("v")


def build_neutrino_mass_matrix(
    c5_matrix: sp.MatrixBase,
    *,
    vev: sp.Expr = v,
) -> sp.Matrix:
    """Return m_nu = -v^2 C5 in the convention currently used by the code."""

    return (-vev**2 * c5_matrix)


def run_neutrino_mass_stage(
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

    flavor_result = flavor_match_from_c5_file(
        c5_path,
        n_lepton=3,
        n_heavy=n_heavy,
        split_heavy_masses=True,
    )

    c5_matrix = flavor_result["K"]
    mass_matrix = build_neutrino_mass_matrix(c5_matrix)

    symmetry_difference = mass_matrix - mass_matrix.T

    if any(sp.simplify(entry) != 0 for entry in symmetry_difference):
        raise RuntimeError(
            "Neutrino mass matrix is not symmetric."
        )

    mass_path = data_dir / "neutrino_mass_matrix.txt"
    mass_path.write_text(
        sp.sstr(mass_matrix) + "\n",
        encoding="utf-8",
    )

    summary = {
        "NeutrinoMassStatus": "Success",
        "NeutrinoMassConvention": "m_nu = -v^2 C5",
        "NeutrinoMassMatrixFile": mass_path.relative_to(output_dir).as_posix(),
        "NeutrinoMassMatrixSymmetric": True,
    }

    summary_path = data_dir / "neutrino_mass_summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )

    return summary
