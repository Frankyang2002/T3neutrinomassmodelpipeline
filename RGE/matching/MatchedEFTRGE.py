from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

from RGE.general.AnomalousDimensions import calculate_complete_master_rge
from RGE.general.FermionBasis import build_gauge_sectors
from RGE.general.GaugeGenerators import g2
from RGE.general.MasterWeinbergRGE import MasterRGEInputs
from RGE.general.WilsonTensors import SparseWilsonLookup
from RGE.matching.WeinbergWilsonAdapter import (
    build_sm_eft,
    build_sm_yukawa,
    yd,
    ye,
    yu,
    build_weinberg_wilson_tensor,
    validate_weinberg_tensor_symmetry,
)

lambdaH = sp.Symbol("lambdaH")



def _matching_square_bracket(text: str, open_index: int) -> int:
    if text[open_index] != "[":
        raise ValueError("Expected '[' at open_index.")

    depth = 0

    for index in range(open_index, len(text)):
        char = text[index]

        if char == "[":
            depth += 1
        elif char == "]":
            depth -= 1

            if depth == 0:
                return index

    raise ValueError("Unbalanced Mathematica square brackets.")


def _first_top_level_argument(arguments: str) -> str:
    depth = 0

    for pos, char in enumerate(arguments):
        if char in "[{(":
            depth += 1
        elif char in "]})":
            depth -= 1
        elif char == "," and depth == 0:
            return arguments[:pos].strip()

    return arguments.strip()


def _replace_matchete_barred_couplings(text: str) -> str:
    marker = "Bar[Coupling["
    result = text

    while marker in result:
        start = result.index(marker)
        bar_open = start + len("Bar")
        bar_close = _matching_square_bracket(result, bar_open)
        inner = result[bar_open + 1 : bar_close]

        if not inner.startswith("Coupling["):
            raise ValueError(
                f"Unexpected Bar contents while parsing Matchete C5: {inner}"
            )

        coupling_open = inner.index("[")
        coupling_close = _matching_square_bracket(inner, coupling_open)
        arguments = inner[coupling_open + 1 : coupling_close]
        name = _first_top_level_argument(arguments)

        result = result[:start] + f"conjugate({name})" + result[bar_close + 1 :]

    return result


def _replace_matchete_couplings(text: str) -> str:
    marker = "Coupling["
    result = text

    while marker in result:
        start = result.index(marker)
        open_index = start + len("Coupling")
        close_index = _matching_square_bracket(result, open_index)
        arguments = result[open_index + 1 : close_index]
        name = _first_top_level_argument(arguments)
        result = result[:start] + name + result[close_index + 1 :]

    return result


def parse_matchete_c5(text: str) -> sp.Expr:
    """Parse the Matchete InputForm subset used by c5_coefficient.txt."""

    stripped = text.strip()

    if not stripped:
        raise ValueError("The matched C5 coefficient file is empty.")

    cleaned = _replace_matchete_barred_couplings(stripped)
    cleaned = _replace_matchete_couplings(cleaned)
    cleaned = cleaned.replace("Conjugate[", "conjugate(")
    cleaned = cleaned.replace("Log[", "log(")
    cleaned = cleaned.replace("Sqrt[", "sqrt(")
    cleaned = cleaned.replace("]", ")")
    cleaned = cleaned.replace("^", "**")

    return sp.simplify(
        sp.sympify(
            cleaned,
            locals={
                "log": sp.log,
                "sqrt": sp.sqrt,
                "conjugate": sp.conjugate,
            },
        )
    )


def higgs_quartic(a: int, b: int, c: int, d: int) -> sp.Expr:
    """lambda_abcd for V=(lambdaH/2)(H^dagger H)^2."""

    delta = lambda x, y: sp.Integer(1 if x == y else 0)

    return lambdaH * (
        delta(a, b) * delta(c, d)
        + delta(a, c) * delta(b, d)
        + delta(a, d) * delta(b, c)
    )


def calculate_matched_eft_rge(kappa: sp.Expr) -> dict:
    """Evaluate the matched Weinberg coefficient in the active SM EFT."""

    scalar_model, fermion_basis = build_sm_eft()
    yukawa = build_sm_yukawa(scalar_model, fermion_basis)

    wilson = build_weinberg_wilson_tensor(
        scalar_model,
        fermion_basis,
        kappa,
    )
    validate_weinberg_tensor_symmetry(wilson)
    coefficient = SparseWilsonLookup(wilson)

    nu = fermion_basis.global_index("L", 1)
    h0_r = scalar_model.block("H").local_to_global(3)
    component = (nu, nu, h0_r, h0_r)
    c_value = coefficient[component]

    inputs = MasterRGEInputs(
        fermion_dimension=fermion_basis.dimension,
        yukawa=yukawa,
        quartic=higgs_quartic,
        gauge_sectors=build_gauge_sectors(
            scalar_model,
            fermion_basis,
        ),
    )

    contributions = calculate_complete_master_rge(
        model=scalar_model,
        inputs=inputs,
        output_component=component,
        coefficient=coefficient,
    )

    beta_ratio = sp.factor(
        sp.simplify(contributions["total"] / c_value)
    )

    expected_ratio = (
        -3 * g2**2
        + 2 * lambdaH
        + 6 * yu * sp.conjugate(yu)
        + 6 * yd * sp.conjugate(yd)
        - ye * sp.conjugate(ye)
    )

    return {
        "component": component,
        "contributions": contributions,
        "beta_ratio": beta_ratio,
        "expected_ratio": expected_ratio,
        "difference": sp.factor(
            sp.expand(beta_ratio - expected_ratio)
        ),
        "dot_kappa": sp.factor(
            sp.simplify(2 * contributions["total"])
        ),
    }


def run_matched_eft_rge(
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

    # The hierarchical pipeline supplies the physical full-flavor JSON.  The
    # purpose of this stage is only the universal one-generation SMEFT RGE
    # benchmark, so reduce that physical coefficient to a 1x1 flavor problem.
    # Import locally to avoid the module-level cycle: FinalWeinbergAdapter uses
    # parse_matchete_c5 from this module for the hard threshold expression.
    from RGE.matching.FinalWeinbergAdapter import (
        is_final_weinberg_json,
        load_final_weinberg_flavor_matrix,
    )

    if is_final_weinberg_json(c5_path):
        adapted = load_final_weinberg_flavor_matrix(
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

    result = calculate_matched_eft_rge(kappa)

    if sp.simplify(result["difference"]) != 0:
        raise RuntimeError(
            "Matched-EFT RGE failed its internal SM benchmark: "
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


def _output_relative_path(path: Path, output_dir: Path) -> str:
    """Return a portable report path relative to one model directory."""

    try:
        return path.relative_to(output_dir).as_posix()
    except ValueError:
        return str(path)
