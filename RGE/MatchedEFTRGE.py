from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

from RGE.GeneralWeinbergRGEGenerator import (
    ComplexScalar,
    MasterRGEInputs,
    RGEModel,
    calculate_complete_master_rge,
    g2,
)
from RGE.T3RGETensors import FermionBasis, WeylFermion, build_gauge_sectors
from RGE.WeinbergWilsonAdapter import (
    build_weinberg_wilson_tensor,
    validate_weinberg_tensor_symmetry,
)

lambdaH = sp.Symbol("lambdaH")
ye, yu, yd = sp.symbols("ye yu yd")


class SparseWilsonLookup:
    """Sparse C_ijab lookup supporting the [] interface used by Eq. (4.85)."""

    def __init__(self, components):
        self.components = components

    def __getitem__(self, key):
        return sp.sympify(self.components.get(tuple(key), sp.S.Zero))


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


def _add_symmetric_yukawa(entries, i, j, a, value) -> None:
    entries[(i, j, a)] = value
    entries[(j, i, a)] = value


def build_sm_yukawa(
    scalar_model: RGEModel,
    fermion_basis: FermionBasis,
):
    """Build the one-generation SM y_ija tensor including colour multiplicity."""

    entries = {}
    h = scalar_model.block("H")
    hp_r = h.local_to_global(1)
    hp_i = h.local_to_global(2)
    h0_r = h.local_to_global(3)
    h0_i = h.local_to_global(4)
    root2 = sp.sqrt(2)

    nu = fermion_basis.global_index("L", 1)
    e = fermion_basis.global_index("L", 2)
    ec = fermion_basis.global_index("eC", 1)

    _add_symmetric_yukawa(entries, nu, ec, hp_r, ye / root2)
    _add_symmetric_yukawa(entries, nu, ec, hp_i, -sp.I * ye / root2)
    _add_symmetric_yukawa(entries, e, ec, h0_r, ye / root2)
    _add_symmetric_yukawa(entries, e, ec, h0_i, -sp.I * ye / root2)

    for colour in range(1, 4):
        q_up = fermion_basis.global_index(f"Q{colour}", 1)
        q_down = fermion_basis.global_index(f"Q{colour}", 2)
        uc = fermion_basis.global_index(f"uC{colour}", 1)
        dc = fermion_basis.global_index(f"dC{colour}", 1)

        _add_symmetric_yukawa(entries, q_up, uc, h0_r, yu / root2)
        _add_symmetric_yukawa(entries, q_up, uc, h0_i, sp.I * yu / root2)
        _add_symmetric_yukawa(entries, q_down, uc, hp_r, -yu / root2)
        _add_symmetric_yukawa(entries, q_down, uc, hp_i, -sp.I * yu / root2)

        _add_symmetric_yukawa(entries, q_up, dc, hp_r, yd / root2)
        _add_symmetric_yukawa(entries, q_up, dc, hp_i, -sp.I * yd / root2)
        _add_symmetric_yukawa(entries, q_down, dc, h0_r, yd / root2)
        _add_symmetric_yukawa(entries, q_down, dc, h0_i, -sp.I * yd / root2)

    def yukawa(i: int, j: int, a: int) -> sp.Expr:
        return sp.sympify(entries.get((i, j, a), sp.S.Zero))

    return yukawa


def build_sm_eft() -> tuple[RGEModel, FermionBasis]:
    """Construct the active one-generation SM EFT basis below the T3 threshold."""

    scalar_model = RGEModel(
        (
            ComplexScalar(
                name="H",
                su2_dimension=2,
                hypercharge=sp.Rational(1, 2),
            ),
        )
    )

    fermions = [
        WeylFermion(
            name="L",
            su2_dimension=2,
            hypercharge=-sp.Rational(1, 2),
        ),
        WeylFermion(
            name="eC",
            su2_dimension=1,
            hypercharge=sp.Integer(1),
        ),
    ]

    for colour in range(1, 4):
        fermions.extend(
            [
                WeylFermion(
                    name=f"Q{colour}",
                    su2_dimension=2,
                    hypercharge=sp.Rational(1, 6),
                ),
                WeylFermion(
                    name=f"uC{colour}",
                    su2_dimension=1,
                    hypercharge=-sp.Rational(2, 3),
                ),
                WeylFermion(
                    name=f"dC{colour}",
                    su2_dimension=1,
                    hypercharge=sp.Rational(1, 3),
                ),
            ]
        )

    return scalar_model, FermionBasis(fermions=tuple(fermions))


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


def run_matched_eft_rge(c5_path: Path, output_dir: Path) -> dict:
    """Read one matched C5 coefficient, calculate its EFT beta, and save outputs."""

    c5_path = Path(c5_path)
    output_dir = Path(output_dir)

    kappa = parse_matchete_c5(
        c5_path.read_text(encoding="utf-8")
    )
    result = calculate_matched_eft_rge(kappa)

    if sp.simplify(result["difference"]) != 0:
        raise RuntimeError(
            "Matched-EFT RGE failed its internal SM benchmark: "
            f"{result['difference']}"
        )

    beta_path = output_dir / "c5_beta.txt"
    ratio_path = output_dir / "c5_beta_over_c5.txt"
    summary_path = output_dir / "rge_summary.json"

    beta_path.write_text(
        str(result["dot_kappa"]) + "\n",
        encoding="utf-8",
    )
    ratio_path.write_text(
        str(result["beta_ratio"]) + "\n",
        encoding="utf-8",
    )

    summary = {
        "RGEStatus": "Success",
        "TheoryBelowThreshold": "SMEFT",
        "MatchingAssumption": "Common heavy T3 threshold",
        "C5InputFile": c5_path.name,
        "C5BetaFile": beta_path.name,
        "C5BetaOverC5File": ratio_path.name,
        "RGEComponent": list(result["component"]),
        "BetaOverC5": str(result["beta_ratio"]),
        "C5Beta": str(result["dot_kappa"]),
    }

    summary_path.write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )

    return summary
