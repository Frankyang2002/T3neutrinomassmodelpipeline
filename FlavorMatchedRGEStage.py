from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

from FlavorMatchedC5 import flavor_match_from_c5_file, write_flavor_outputs
from SMEFTWeinbergFlavorRGE import beta_weinberg_matrix, symbolic_complex_matrix


def run_flavor_matched_rge(
    c5_path: Path,
    output_dir: Path,
    *,
    n_heavy: int = 3,
) -> dict:
    """Construct the full flavor C5 matrix and its SMEFT beta matrix."""

    c5_path = Path(c5_path)
    output_dir = Path(output_dir)

    flavor_result = flavor_match_from_c5_file(
        c5_path,
        n_lepton=3,
        n_heavy=n_heavy,
        split_heavy_masses=True,
    )

    flavor_files = write_flavor_outputs(
        output_dir,
        flavor_result,
    )

    K = flavor_result["K"]

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

    beta_path = output_dir / "c5_flavor_beta_matrix.txt"

    beta_path.write_text(
        sp.sstr(beta_K) + "\n",
        encoding="utf-8",
    )

    summary = {
        "FlavorRGEStatus": "Success",
        "FlavorMatchingAssumption": (
            "3 lepton generations, 3 diagonal heavy-fermion generations"
        ),
        "C5FlavorMatrixFile": flavor_files["C5FlavorMatrixFile"],
        "C5LoopKernelFile": flavor_files["C5LoopKernelFile"],
        "C5FlavorBetaMatrixFile": beta_path.name,
        "LeptonFlavorGenerations": 3,
        "HeavyFlavorGenerations": n_heavy,
        "FlavorC5Symmetric": True,
        "FlavorBetaC5Symmetric": True,
    }

    summary_path = output_dir / "flavor_rge_summary.json"

    summary_path.write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )

    return summary
