from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

from RGE.matching.FlavorMatchedC5 import (
    flavor_match_from_c5_file,
    is_final_weinberg_json,
    load_final_weinberg_flavor_matrix,
    write_flavor_outputs,
)
from RGE.running.WeinbergRunning import beta_weinberg_matrix, symbolic_complex_matrix


def run_flavor_matched_rge(
    c5_path: Path,
    output_dir: Path,
    *,
    n_heavy: int = 3,
    debug_outputs: bool = False,
) -> dict:
    """Construct the full flavor C5 matrix and its SMEFT beta matrix."""

    c5_path = Path(c5_path)
    output_dir = Path(output_dir)
    data_dir = output_dir / "data"
    data_dir.mkdir(parents=True, exist_ok=True)

    if is_final_weinberg_json(c5_path):
        flavor_result = load_final_weinberg_flavor_matrix(
            c5_path,
            n_lepton=3,
            n_heavy=n_heavy,
            split_heavy_masses=True,
        )
        K = flavor_result["K"]

        flavor_files = {
            "C5FlavorMatrixFile": "",
            "C5LoopKernelFile": "",
        }

        if debug_outputs:
            debug_dir = output_dir / "debug"
            debug_dir.mkdir(parents=True, exist_ok=True)

            matrix_path = debug_dir / "c5_flavor_matrix.txt"
            hard_path = debug_dir / "c5_hard_kernel.txt"
            running_path = debug_dir / "c5_running_kernel.txt"

            matrix_path.write_text(sp.sstr(K) + "\n", encoding="utf-8")
            hard_path.write_text(
                sp.sstr(flavor_result["hard_kernel"]) + "\n",
                encoding="utf-8",
            )
            running_path.write_text(
                sp.sstr(flavor_result["running_kernel"]) + "\n",
                encoding="utf-8",
            )

            flavor_files = {
                "C5FlavorMatrixFile": matrix_path.relative_to(output_dir).as_posix(),
                "C5LoopKernelFile": hard_path.relative_to(output_dir).as_posix(),
                "C5RunningKernelFile": running_path.relative_to(output_dir).as_posix(),
            }

        input_kind = "final_weinberg_json"
        matching_assumption = (
            "Physical Majorana C5 from hierarchical MSbar hard matching "
            "plus EFT1 running"
        )
    else:
        flavor_result = flavor_match_from_c5_file(
            c5_path,
            n_lepton=3,
            n_heavy=n_heavy,
            split_heavy_masses=True,
        )

        flavor_files = write_flavor_outputs(
            output_dir,
            flavor_result,
            debug_outputs=debug_outputs,
        )

        K = flavor_result["K"]
        input_kind = "legacy_scalar_c5"
        matching_assumption = (
            "3 lepton generations, 3 diagonal heavy-fermion generations"
        )

    Ye = symbolic_complex_matrix("Ye", 3, 3)
    Yu = symbolic_complex_matrix("Yu", 3, 3)
    Yd = symbolic_complex_matrix("Yd", 3, 3)

    beta_K = beta_weinberg_matrix(
        K,
        Ye,
        Yu,
        Yd,
        expand_result=False,
    )

    beta_path = data_dir / "c5_flavor_beta_matrix.txt"

    beta_path.write_text(
        sp.sstr(beta_K) + "\n",
        encoding="utf-8",
    )

    summary = {
        "FlavorRGEStatus": "Success",
        "FlavorMatchingAssumption": matching_assumption,
        "C5InputKind": input_kind,
        "C5FlavorMatrixFile": flavor_files["C5FlavorMatrixFile"],
        "C5LoopKernelFile": flavor_files["C5LoopKernelFile"],
        "C5FlavorBetaMatrixFile": beta_path.relative_to(output_dir).as_posix(),
        "LeptonFlavorGenerations": 3,
        "HeavyFlavorGenerations": n_heavy,
        "FlavorC5Symmetric": True,
        "FlavorBetaC5Symmetric": True,
    }

    summary_path = data_dir / "flavor_rge_summary.json"

    summary_path.write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )

    return summary


# ---------------------------------------------------------------------------
# Symbolic neutrino-mass stage
# Consolidated from the former standalone neutrino-mass stage module.
# ---------------------------------------------------------------------------

v = sp.Symbol("v")


def build_neutrino_mass_matrix(
    c5_matrix: sp.MatrixBase,
    *,
    vev: sp.Expr = v,
) -> sp.Matrix:
    """Return m_nu = -(v^2/2) C5 for H^0 = (v+h)/sqrt(2)."""

    return (-sp.Rational(1, 2) * vev**2 * c5_matrix)


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

    if is_final_weinberg_json(c5_path):
        flavor_result = load_final_weinberg_flavor_matrix(
            c5_path,
            n_lepton=3,
            n_heavy=n_heavy,
            split_heavy_masses=True,
        )
        c5_input_kind = "final_weinberg_json"
    else:
        flavor_result = flavor_match_from_c5_file(
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
