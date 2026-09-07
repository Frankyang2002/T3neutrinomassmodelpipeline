from __future__ import annotations

"""
First end-to-end T3 RGE smoke test.

Inputs
------
output/T3_B_alpha_m1/rge_tensor_exchange.json
    UV representation, scalar quartics and raw Yukawa CG components.

output/T3_B_alpha_m1/c5_coefficient.txt
    Matched Weinberg coefficient from the existing matching pipeline.

Pipeline exercised
------------------
raw T3 exchange
  -> scalar real basis
  -> global Weyl basis
  -> y_ija
  -> lambda_abcd
  -> gauge generators t^A, theta^A
  -> matched kappa
  -> C_ijab
  -> Eq. (4.85)
  -> dot C_(nu nu H0_R H0_R).

For the normalization used by WeinbergWilsonAdapter,

    C_(nu nu H0_R H0_R) = kappa/2,

so the corresponding beta function is

    dot{kappa} = 2 dot{C}_(nu nu H0_R H0_R).
"""

import argparse
import json
from pathlib import Path

import sympy as sp
from sympy.parsing.mathematica import parse_mathematica

from tests.python.GeneralWeinbergRGEGenerator import (
    calculate_complete_master_rge,
)
from tests.python.Weinberg.T3RGETensors import (
    ComplexQuarticComponent,
    ComplexScalarFactor,
    T3RGETensors,
    quartic_tensor_from_components,
    scalar_model_from_exchange,
    wilson_component_function,
)
from tests.python.Weinberg.T3YukawaAdapter import (
    T3WeylConvention,
    real_yukawa_tensor_from_exchange,
)
from tests.python.Weinberg.WeinbergWilsonAdapter import (
    build_weinberg_wilson_tensor,
    validate_weinberg_tensor_symmetry,
)


PROJECT_ROOT = Path(__file__).resolve().parent

DEFAULT_EXCHANGE = (
    PROJECT_ROOT
    / "wolfram"
    / "output"
    / "T3_B_alpha_m1"
    / "rge_tensor_exchange.json"
)

DEFAULT_C5 = (
    PROJECT_ROOT
    / "wolfram"
    / "output"
    / "T3_B_alpha_m1"
    / "c5_coefficient.txt"
)


# -----------------------------------------------------------------------------
# Exact expression parsing
# -----------------------------------------------------------------------------


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

        result = (
            result[:start]
            + f"conjugate({name})"
            + result[bar_close + 1 :]
        )

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


def parse_wolfram_inputform(text: str) -> sp.Expr:
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

    try:
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
    except Exception as exc:
        raise ValueError(
            "Could not parse canonicalised Matchete C5 coefficient.\n"
            f"Raw expression:\n{stripped}\n\n"
            f"Canonicalised expression:\n{cleaned}"
        ) from exc


def quartic_components_from_exchange(
    data: dict,
) -> list[ComplexQuarticComponent]:
    """Convert the Wolfram complex-basis quartic JSON into Python objects."""

    result: list[ComplexQuarticComponent] = []

    for item in data.get("quartic_components", []):
        coefficient = parse_wolfram_inputform(str(item["coefficient"]))

        raw_factors = item["factors"]

        if len(raw_factors) != 4:
            raise ValueError(
                "Each exported quartic term must contain exactly four factors."
            )

        factors = tuple(
            ComplexScalarFactor(
                scalar_name=str(factor["scalar_name"]),
                component=int(factor["component"]),
                conjugated=bool(factor.get("conjugated", False)),
            )
            for factor in raw_factors
        )

        result.append(
            ComplexQuarticComponent(
                coefficient=coefficient,
                factors=factors,
            )
        )

    if not result:
        raise ValueError("No quartic component terms were found in the exchange.")

    return result


# -----------------------------------------------------------------------------
# End-to-end construction
# -----------------------------------------------------------------------------


def build_t3b_rge_tensors(
    exchange_path: Path,
    c5_path: Path,
) -> tuple[T3RGETensors, sp.Expr]:
    data = json.loads(exchange_path.read_text(encoding="utf-8"))

    scalar_model = scalar_model_from_exchange(data)

    scalar_model_y, fermion_basis, yukawa_tensor = (
        real_yukawa_tensor_from_exchange(
            data,
            T3WeylConvention(
                lepton_multiplicity=1,
                heavy_fermion_multiplicity=1,
                symmetrize_fermion_indices=True,
            ),
        )
    )

    # The Yukawa adapter and the generic scalar exchange reader must agree on
    # exactly the same scalar real basis.
    if scalar_model.total_real_scalar_dimension != (
        scalar_model_y.total_real_scalar_dimension
    ):
        raise AssertionError("Scalar-basis dimensions disagree between adapters.")

    for name in scalar_model.blocks:
        left = scalar_model.block(name)
        right = scalar_model_y.block(name)

        if (
            left.first != right.first
            or left.last != right.last
            or left.scalar != right.scalar
        ):
            raise AssertionError(
                f"Scalar-basis block mismatch for {name}."
            )

    complex_quartics = quartic_components_from_exchange(data)

    print(
        "Converting",
        len(complex_quartics),
        "complex quartic component terms to lambda_abcd..."
    )

    quartic_tensor = quartic_tensor_from_components(
        scalar_model,
        complex_quartics,
    )

    kappa = parse_wolfram_inputform(
        c5_path.read_text(encoding="utf-8")
    )

    wilson_tensor = build_weinberg_wilson_tensor(
        scalar_model,
        fermion_basis,
        kappa,
    )

    validate_weinberg_tensor_symmetry(wilson_tensor)

    tensors = T3RGETensors(
        scalar_model=scalar_model,
        fermion_basis=fermion_basis,
        yukawa_components=yukawa_tensor,
        quartic_components=quartic_tensor,
        wilson_components=wilson_tensor,
    )

    return tensors, kappa


# -----------------------------------------------------------------------------
# Smoke RGE
# -----------------------------------------------------------------------------


def run_neutral_component(
    tensors: T3RGETensors,
) -> tuple[tuple[int, int, int, int], dict[str, sp.Expr]]:
    """Run Eq. (4.85) for C_(nu nu H0_R H0_R)."""

    L = tensors.fermion_basis.block("L")
    H = tensors.scalar_model.block("H")

    nu = tensors.fermion_basis.global_index("L", 1)

    # H component convention:
    #   complex component 1 = H+
    #   complex component 2 = H0
    # and each complex component is stored as (R,I).
    h0_r = H.local_to_global(3)

    component = (nu, nu, h0_r, h0_r)

    inputs = tensors.master_inputs()
    coefficient = wilson_component_function(tensors.wilson_components)

    result = calculate_complete_master_rge(
        model=tensors.scalar_model,
        inputs=inputs,
        output_component=component,
        coefficient=coefficient,
        simplify_each=True,
    )

    return component, result


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run the first complete T3-B Eq. (4.85) RGE smoke test."
    )
    parser.add_argument(
        "--exchange",
        type=Path,
        default=DEFAULT_EXCHANGE,
    )
    parser.add_argument(
        "--c5",
        type=Path,
        default=DEFAULT_C5,
    )

    args = parser.parse_args()

    if not args.exchange.exists():
        raise FileNotFoundError(
            f"Missing RGE exchange JSON: {args.exchange}\n"
            "Run TestT3RGEExchangeExport.wl first."
        )

    if not args.c5.exists():
        raise FileNotFoundError(
            f"Missing matched C5 coefficient: {args.c5}\n"
            "Run the normal matching pipeline first."
        )

    print("Exchange:", args.exchange)
    print("Matched C5:", args.c5)
    print()

    tensors, kappa = build_t3b_rge_tensors(
        args.exchange,
        args.c5,
    )

    print()
    print("Model dimensions")
    print(
        "  real scalars:",
        tensors.scalar_model.total_real_scalar_dimension,
    )
    print(
        "  Weyl fermions:",
        tensors.fermion_basis.dimension,
    )
    print(
        "  nonzero y_ija:",
        len(tensors.yukawa_components),
    )
    print(
        "  nonzero lambda_abcd:",
        len(tensors.quartic_components),
    )
    print(
        "  nonzero C_ijab:",
        len(tensors.wilson_components),
    )

    print()
    print("Matched kappa:")
    print(" ", kappa)

    component, contributions = run_neutral_component(tensors)

    print()
    print(
        "RGE component:",
        f"dot C{component}",
    )

    ordered_terms = (
        "yukawa_wavefunction",
        "scalar_pair",
        "mixed_yukawa_gauge",
        "crossed_yukawa",
        "conjugate_coefficient_yukawa",
        "scalar_anomalous_dimension",
        "fermion_anomalous_dimension",
    )

    for name in ordered_terms:
        print()
        print(name + ":")
        print(" ", contributions[name])

    print()
    print("dot C:")
    print(" ", contributions["total"])

    beta_kappa = sp.simplify(
        2 * contributions["total"]
    )

    print()
    print("dot kappa = 2 dot C_(nu nu H0_R H0_R):")
    print(" ", beta_kappa)

    print()
    print("PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
