from __future__ import annotations

"""
End-to-end matched-EFT Weinberg-RGE smoke test for T3-B.

Pipeline
--------
1. Read Matchete's matched Weinberg coefficient kappa(M).
2. Treat kappa(M) as the boundary condition of the EFT below the heavy T3
   threshold.
3. Build the active one-generation Standard Model EFT:
       H, L, eC, Q, uC, dC
   with three colour copies for the quark fields.
4. Construct the full Weinberg tensor C_ijab from kappa(M).
5. Evaluate the validated generic Eq. (4.85) engine.

The heavy T3 fields S1, S2, F and their Yukawas y1, y2 do NOT appear as
active RGE fields below the matching threshold.  Their dependence remains
inside the matched boundary condition kappa(M).

Expected one-generation EFT result:
    beta_kappa / kappa
      = -3 g2^2
        + 2 lambdaH
        + 6 |yu|^2
        + 6 |yd|^2
        - |ye|^2.
"""

import argparse
from pathlib import Path

import sympy as sp


from tests.python.Weinberg.GeneralWeinbergRGEGenerator import (
    ComplexScalar,
    MasterRGEInputs,
    RGEModel,
    calculate_complete_master_rge,
    g2,
)
from RGE.general.T3RGETensors import (
    FermionBasis,
    WeylFermion,
    build_gauge_sectors,
)
from tests.python.Weinberg.WeinbergWilsonAdapter import (
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
        return sp.sympify(
            self.components.get(tuple(key), sp.S.Zero)
        )


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
    """Parse the Matchete InputForm subset used by the matched C5 output."""

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


def higgs_quartic(
    a: int,
    b: int,
    c: int,
    d: int,
) -> sp.Expr:
    """lambda_abcd for V=(lambdaH/2)(H^\dagger H)^2."""

    delta = lambda x, y: sp.Integer(1 if x == y else 0)

    return lambdaH * (
        delta(a, b) * delta(c, d)
        + delta(a, c) * delta(b, d)
        + delta(a, d) * delta(b, c)
    )


def add_symmetric_yukawa(entries, i, j, a, value):
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

    # Charged lepton: eC H^\dagger L.
    add_symmetric_yukawa(entries, nu, ec, hp_r, ye / root2)
    add_symmetric_yukawa(entries, nu, ec, hp_i, -sp.I * ye / root2)
    add_symmetric_yukawa(entries, e, ec, h0_r, ye / root2)
    add_symmetric_yukawa(entries, e, ec, h0_i, -sp.I * ye / root2)

    for colour in range(1, 4):
        q_up = fermion_basis.global_index(f"Q{colour}", 1)
        q_down = fermion_basis.global_index(f"Q{colour}", 2)
        uc = fermion_basis.global_index(f"uC{colour}", 1)
        dc = fermion_basis.global_index(f"dC{colour}", 1)

        # Up-type: uC (Q . H).
        add_symmetric_yukawa(entries, q_up, uc, h0_r, yu / root2)
        add_symmetric_yukawa(entries, q_up, uc, h0_i, sp.I * yu / root2)
        add_symmetric_yukawa(entries, q_down, uc, hp_r, -yu / root2)
        add_symmetric_yukawa(entries, q_down, uc, hp_i, -sp.I * yu / root2)

        # Down-type: dC H^\dagger Q.
        add_symmetric_yukawa(entries, q_up, dc, hp_r, yd / root2)
        add_symmetric_yukawa(entries, q_up, dc, hp_i, -sp.I * yd / root2)
        add_symmetric_yukawa(entries, q_down, dc, h0_r, yd / root2)
        add_symmetric_yukawa(entries, q_down, dc, h0_i, -sp.I * yd / root2)

    def yukawa(i: int, j: int, a: int) -> sp.Expr:
        return sp.sympify(
            entries.get((i, j, a), sp.S.Zero)
        )

    return yukawa


def build_sm_eft():
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

    fermion_basis = FermionBasis(
        fermions=tuple(fermions)
    )

    return scalar_model, fermion_basis


def run_matched_eft_rge(kappa: sp.Expr):
    scalar_model, fermion_basis = build_sm_eft()

    yukawa = build_sm_yukawa(
        scalar_model,
        fermion_basis,
    )

    wilson = build_weinberg_wilson_tensor(
        scalar_model,
        fermion_basis,
        kappa,
    )
    validate_weinberg_tensor_symmetry(wilson)

    C = SparseWilsonLookup(wilson)

    nu = fermion_basis.global_index("L", 1)
    h0_r = scalar_model.block("H").local_to_global(3)

    component = (nu, nu, h0_r, h0_r)
    c_value = C[component]

    assert sp.simplify(c_value - kappa / 2) == 0

    gauge_sectors = build_gauge_sectors(
        scalar_model,
        fermion_basis,
    )

    inputs = MasterRGEInputs(
        fermion_dimension=fermion_basis.dimension,
        yukawa=yukawa,
        quartic=higgs_quartic,
        gauge_sectors=gauge_sectors,
    )

    result = calculate_complete_master_rge(
        model=scalar_model,
        inputs=inputs,
        output_component=component,
        coefficient=C,
    )

    return component, c_value, result


def default_c5_path() -> Path:
    return (
        Path(__file__).resolve().parent
        / "wolfram"
        / "output"
        / "T3_B_alpha_m1"
        / "c5_coefficient.txt"
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run the matched T3-B C5 coefficient through the SM EFT RGE."
    )
    parser.add_argument(
        "--c5",
        type=Path,
        default=default_c5_path(),
        help="Path to Matchete c5_coefficient.txt.",
    )
    args = parser.parse_args()

    c5_path = args.c5.resolve()

    if not c5_path.exists():
        raise FileNotFoundError(c5_path)

    print("Matched C5:", c5_path)
    print()

    kappa = parse_wolfram_inputform(
        c5_path.read_text(encoding="utf-8")
    )

    print("Matched boundary condition kappa(M):")
    print(" ", kappa)
    print()

    component, c_value, result = run_matched_eft_rge(kappa)

    print("Active EFT below threshold:")
    print("  Scalars: H")
    print("  Fermions: L, eC, Q, uC, dC")
    print("  Heavy T3 fields S1, S2, F: integrated out")
    print()

    print("RGE component:")
    print(f"  dot C{component}")
    print()

    order = (
        "yukawa_wavefunction",
        "scalar_pair",
        "mixed_yukawa_gauge",
        "crossed_yukawa",
        "conjugate_coefficient_yukawa",
        "scalar_anomalous_dimension",
        "fermion_anomalous_dimension",
        "total",
    )

    print("Contributions to (dot kappa)/kappa:")

    for name in order:
        ratio = sp.factor(
            sp.simplify(result[name] / c_value)
        )
        print(f"  {name:32s}: {ratio}")

    beta_ratio = sp.factor(
        sp.simplify(result["total"] / c_value)
    )

    expected_ratio = (
        -3 * g2**2
        + 2 * lambdaH
        + 6 * yu * sp.conjugate(yu)
        + 6 * yd * sp.conjugate(yd)
        - ye * sp.conjugate(ye)
    )

    difference = sp.factor(
        sp.expand(beta_ratio - expected_ratio)
    )

    dot_kappa = sp.factor(
        sp.simplify(2 * result["total"])
    )

    expected_dot_kappa = sp.factor(
        sp.simplify(kappa * expected_ratio)
    )

    print()
    print("beta_kappa/kappa:")
    print(" ", beta_ratio)

    print()
    print("Expected beta_kappa/kappa:")
    print(" ", expected_ratio)

    print()
    print("Difference:")
    print(" ", difference)

    print()
    print("dot kappa:")
    print(" ", dot_kappa)

    print()
    print("Expected dot kappa = kappa(M) * SM running factor:")
    print(" ", expected_dot_kappa)

    print()

    if sp.simplify(difference) != 0:
        print("FAIL")
        return 1

    if sp.simplify(dot_kappa - expected_dot_kappa) != 0:
        print("FAIL: dot kappa does not factorise as expected.")
        return 1

    print("PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
