"""Output stage for the one-generation SMEFT Weinberg RGE benchmark."""

from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

from RGE.matching.FinalWeinbergAdapter import (
    is_hierarchical_final_c5,
    load_hierarchical_majorana_c5,
)
from RGE.matching.MatcheteC5Parsing import parse_matchete_c5
from RGE.running.weinberg.OneGenerationWeinbergBenchmark import (
    calculate_one_generation_weinberg_benchmark,
)


def _output_relative_path(path: Path, output_dir: Path) -> str:
    """Return a portable report path relative to one model directory."""

    try:
        return path.relative_to(output_dir).as_posix()
    except ValueError:
        return str(path)


def run_one_generation_weinberg_benchmark(
    c5_path: Path,
    output_dir: Path,
    *,
    debug_outputs: bool = False,
) -> dict:
    """Read one matched C5 coefficient, calculate its EFT beta, and save outputs."""

    c5_path = Path(c5_path)
    output_dir = Path(output_dir)
    data_dir = output_dir / "data"
    debug_dir = output_dir / "debug"
    data_dir.mkdir(parents=True, exist_ok=True)

    if debug_outputs:
        debug_dir.mkdir(parents=True, exist_ok=True)

    # The hierarchical pipeline supplies the physical full-flavor JSON. This
    # stage is the universal one-generation SMEFT RGE benchmark, so reduce that
    # physical coefficient to a 1x1 flavor problem.
    if is_hierarchical_final_c5(c5_path):
        adapted = load_hierarchical_majorana_c5(
            c5_path,
            n_lepton=1,
            n_heavy=1,
            split_heavy_masses=False,
        )
        kappa = sp.simplify(adapted["K"][0, 0])
        c5_input_kind = "final_weinberg_json_one_generation_reduction"
        matching_assumption = (
            "Hierarchical physical Majorana C5 reduced to one generation "
            "for the universal SMEFT RGE benchmark"
        )
    else:
        kappa = parse_matchete_c5(
            c5_path.read_text(encoding="utf-8")
        )
        c5_input_kind = "legacy_scalar_c5"
        matching_assumption = "Common heavy T3 threshold"

    result = calculate_one_generation_weinberg_benchmark(kappa)

    if sp.simplify(result["difference"]) != 0:
        raise RuntimeError(
            "One-generation Weinberg benchmark failed its internal SM check: "
            f"{result['difference']}"
        )

    beta_path = data_dir / "c5_beta.txt"
    ratio_path = debug_dir / "c5_beta_over_c5.txt"
    summary_path = data_dir / "rge_summary.json"

    beta_path.write_text(
        str(result["dot_kappa"]) + "\n",
        encoding="utf-8",
    )

    if debug_outputs:
        ratio_path.write_text(
            str(result["beta_ratio"]) + "\n",
            encoding="utf-8",
        )

    summary = {
        "RGEStatus": "Success",
        "TheoryBelowThreshold": "SMEFT",
        "MatchingAssumption": matching_assumption,
        "C5InputKind": c5_input_kind,
        "C5InputFile": _output_relative_path(c5_path, output_dir),
        "C5BetaFile": _output_relative_path(beta_path, output_dir),
        "C5BetaOverC5File": (
            _output_relative_path(ratio_path, output_dir)
            if debug_outputs
            else ""
        ),
        "RGEComponent": list(result["component"]),
        "BetaOverC5": str(result["beta_ratio"]),
        "C5Beta": str(result["dot_kappa"]),
    }

    summary_path.write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )

    return summary
