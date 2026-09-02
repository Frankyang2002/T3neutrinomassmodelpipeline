from __future__ import annotations

"""
Convert raw Yukawa CG components exported by T3RGETensorExport.wl qinto the
y_ija tensor used by GeneralWeinbergRGEGenerator_complete.py.

This module is deliberately separate from both the Wolfram exporter and the
generic RGE engine.

Convention used here
--------------------
Eq. (4.85) is evaluated in one global basis of left-handed Weyl fermions.

For the T3 Yukawa sector we use two fermion blocks:

    L : the SM lepton doublet
    F : the new T3 fermion multiplet

Each raw Wolfram entry already specifies the SU(2) component of L, F and the
scalar, together with whether F is being used in the conjugate
representation and whether the scalar is conjugated.

The map performed here is therefore

    raw CG component
        -> global Weyl indices (i,j)
        -> complex scalar component
        -> real scalar index a
        -> y_ija.

No overall Yukawa sign or 1/2 is guessed.  The coupling factor exported by
Wolfram is multiplied by the CG coefficient exactly as supplied.

The scotogenic regression below checks the nontrivial SU(2) epsilon
orientation

    L . eta = nu eta^0 - e eta^+

and then checks the complex-to-real conversion

    eta = (R + i I)/sqrt(2).
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

import sympy as sp

from RGE.general.T3RGETensors import (
    ComplexYukawaComponent,
    FermionBasis,
    WeylFermion,
    build_real_yukawa_tensor,
    load_component_exchange,
    parse_exact_expression,
    scalar_model_from_exchange,
)


SQRT2 = sp.sqrt(2)


# -----------------------------------------------------------------------------
# Raw exchange parsing
# -----------------------------------------------------------------------------


def _parse_wolfram_exact(value) -> sp.Expr:
    """Parse the small exact-expression subset emitted by the Wolfram exporter."""

    if isinstance(value, (int, float)):
        return sp.sympify(value)

    if not isinstance(value, str):
        raise TypeError(f"Expected string/int/float, received {type(value)!r}.")

    text = value.strip()

    # Mathematica exact syntax used by the exporter.
    text = text.replace("^", "**")

    while "Sqrt[" in text:
        start = text.rfind("Sqrt[")
        depth = 0
        close = None

        for pos in range(start + 5, len(text)):
            char = text[pos]

            if char == "[":
                depth += 1
            elif char == "]":
                if depth == 0:
                    close = pos
                    break
                depth -= 1

        if close is None:
            raise ValueError(f"Unbalanced Sqrt expression: {value!r}")

        inside = text[start + 5 : close]
        text = text[:start] + f"sqrt({inside})" + text[close + 1 :]

    # Mathematica symbols are valid SymPy symbol names for the coupling names
    # currently exported by the T3 builder.
    return sp.sympify(
        text,
        locals={
            "I": sp.I,
            "sqrt": sp.sqrt,
        },
    )


@dataclass(frozen=True)
class T3WeylConvention:
    """Global Weyl-basis convention used for exported T3 Yukawa components."""

    lepton_multiplicity: int = 1
    heavy_fermion_multiplicity: int = 1
    symmetrize_fermion_indices: bool = True


def fermion_basis_from_raw_exchange(
    data: Mapping,
    convention: T3WeylConvention = T3WeylConvention(),
) -> FermionBasis:
    """Construct L and F blocks from the raw Wolfram exchange metadata."""

    raw = {
        item["name"]: item
        for item in data["raw_fermions"]
    }

    if "L" not in raw or "F" not in raw:
        raise ValueError("Exchange must contain raw fermion metadata for L and F.")

    lepton = raw["L"]
    heavy = raw["F"]

    # The basis contains the un-conjugated physical representation once.
    # Individual Yukawa terms may request its conjugate representation.  That
    # orientation is treated at component conversion rather than by duplicating
    # the entire block.
    return FermionBasis(
        fermions=(
            WeylFermion(
                name="L",
                su2_dimension=int(lepton["su2_dimension"]),
                hypercharge=_parse_wolfram_exact(lepton["hypercharge"]),
                multiplicity=convention.lepton_multiplicity,
                conjugated_representation=False,
            ),
            WeylFermion(
                name="F",
                su2_dimension=int(heavy["su2_dimension"]),
                hypercharge=_parse_wolfram_exact(heavy["hypercharge"]),
                multiplicity=convention.heavy_fermion_multiplicity,
                conjugated_representation=False,
            ),
        )
    )


# -----------------------------------------------------------------------------
# SU(2) conjugate-representation component map
# -----------------------------------------------------------------------------


def su2_duality_matrix(dimension: int) -> sp.Matrix:
    r"""Return the canonical SU(2) duality matrix J for spin j=(d-1)/2.

    In the |j,m> basis ordered m=j,j-1,...,-j,

        J_{m,m'} = (-1)^(j-m) delta_{m',-m}.

    This converts a conjugate-representation component into the ordinary
    representation basis.  For d=2,

        J = [[0, 1],
             [-1, 0]],

    i.e. the familiar epsilon tensor up to the chosen canonical orientation.
    """

    if dimension < 1:
        raise ValueError("SU(2) dimension must be positive.")

    result = sp.zeros(dimension)

    for row in range(dimension):
        column = dimension - 1 - row
        result[row, column] = (-1) ** row

    return result


def conjugate_component_expansion(
    dimension: int,
    component: int,
) -> dict[int, sp.Expr]:
    """Expand one conjugate-representation component in the ordinary basis."""

    if not 1 <= component <= dimension:
        raise IndexError("Component lies outside the SU(2) representation.")

    J = su2_duality_matrix(dimension)

    # column vector e_component transformed by J
    column = component - 1
    result: dict[int, sp.Expr] = {}

    for ordinary_row in range(dimension):
        value = sp.simplify(J[ordinary_row, column])

        if value != 0:
            result[ordinary_row + 1] = value

    return result


# -----------------------------------------------------------------------------
# raw_yukawa_components -> ComplexYukawaComponent
# -----------------------------------------------------------------------------


def complex_yukawa_components_from_exchange(
    data: Mapping,
    convention: T3WeylConvention = T3WeylConvention(),
) -> tuple[FermionBasis, list[ComplexYukawaComponent]]:
    """Convert Wolfram raw Yukawa entries to complex-basis tensor components."""

    basis = fermion_basis_from_raw_exchange(data, convention)
    result: list[ComplexYukawaComponent] = []

    for entry in data.get("raw_yukawa_components", []):
        lepton = entry["lepton"]
        heavy = entry["fermion"]
        scalar = entry["scalar"]

        coupling = sp.Symbol(str(entry["coupling"]))
        cg = _parse_wolfram_exact(entry["cg_coefficient"])
        base_coefficient = sp.simplify(coupling * cg)

        lepton_component = int(lepton["component"])
        heavy_component = int(heavy["component"])

        # L in the exporter is currently marked as conjugated because the
        # Matchete invariant is written with Bar[L].  The physical Weyl tensor
        # convention used here takes L itself as the basis field.  For SU(2),
        # use pseudoreality to map the conjugate doublet back through J.
        if bool(lepton.get("conjugated_representation", False)):
            lepton_expansion = conjugate_component_expansion(
                basis.block("L").fermion.su2_dimension,
                lepton_component,
            )
        else:
            lepton_expansion = {lepton_component: sp.S.One}

        if bool(heavy.get("conjugated_representation", False)):
            heavy_expansion = conjugate_component_expansion(
                basis.block("F").fermion.su2_dimension,
                heavy_component,
            )
        else:
            heavy_expansion = {heavy_component: sp.S.One}

        for l_component, l_factor in lepton_expansion.items():
            for f_component, f_factor in heavy_expansion.items():
                i = basis.global_index("L", l_component)
                j = basis.global_index("F", f_component)

                result.append(
                    ComplexYukawaComponent(
                        fermion_i=i,
                        fermion_j=j,
                        scalar_name=str(scalar["name"]),
                        scalar_component=int(scalar["component"]),
                        scalar_conjugated=bool(scalar.get("conjugated", False)),
                        coefficient=sp.simplify(
                            base_coefficient * l_factor * f_factor
                        ),
                    )
                )

    return basis, result


def real_yukawa_tensor_from_exchange(
    data: Mapping,
    convention: T3WeylConvention = T3WeylConvention(),
):
    """Return the global Weyl basis and sparse real-scalar y_ija tensor."""

    scalar_model = scalar_model_from_exchange(data)
    fermion_basis, complex_components = complex_yukawa_components_from_exchange(
        data,
        convention,
    )

    tensor = build_real_yukawa_tensor(
        scalar_model=scalar_model,
        fermion_basis=fermion_basis,
        components=complex_components,
        symmetrize_fermions=convention.symmetrize_fermion_indices,
    )

    return scalar_model, fermion_basis, tensor


# -----------------------------------------------------------------------------
# Regression tests
# -----------------------------------------------------------------------------


def _self_check_duality() -> None:
    J2 = su2_duality_matrix(2)

    assert J2 == sp.Matrix([
        [0, 1],
        [-1, 0],
    ])

    # J J* = -1 for half-integer isospin.
    assert sp.simplify(J2 * J2.conjugate()) == -sp.eye(2)

    J3 = su2_duality_matrix(3)

    # J J* = +1 for integer isospin.
    assert sp.simplify(J3 * J3.conjugate()) == sp.eye(3)


def _self_check_scotogenic_epsilon() -> None:
    r"""Check L.eta = nu eta0 - e eta+ and its real-field decomposition.

    Use the canonical epsilon contraction

        epsilon_{12}=+1,
        epsilon_{21}=-1.

    With L=(nu,e) and eta=(eta+,eta0),

        epsilon_ab L_a eta_b
          = nu eta0 - e eta+.
    """

    h = sp.Symbol("h")

    exchange = {
        "scalars": [
            {
                "name": "eta",
                "su2_dimension": 2,
                "hypercharge": "1/2",
            }
        ],
        "raw_fermions": [
            {
                "name": "L",
                "su2_dimension": 2,
                "hypercharge": "-1/2",
            },
            {
                "name": "F",
                "su2_dimension": 1,
                "hypercharge": "0",
            },
        ],
    }

    scalar_model = scalar_model_from_exchange(exchange)
    fermion_basis = fermion_basis_from_raw_exchange(exchange)

    # Direct physical complex-basis epsilon terms:
    #   + h nu N eta0
    #   - h e  N eta+
    components = [
        ComplexYukawaComponent(
            fermion_i=fermion_basis.global_index("L", 1),
            fermion_j=fermion_basis.global_index("F", 1),
            scalar_name="eta",
            scalar_component=2,
            scalar_conjugated=False,
            coefficient=h,
        ),
        ComplexYukawaComponent(
            fermion_i=fermion_basis.global_index("L", 2),
            fermion_j=fermion_basis.global_index("F", 1),
            scalar_name="eta",
            scalar_component=1,
            scalar_conjugated=False,
            coefficient=-h,
        ),
    ]

    tensor = build_real_yukawa_tensor(
        scalar_model,
        fermion_basis,
        components,
        symmetrize_fermions=True,
    )

    nu = fermion_basis.global_index("L", 1)
    e = fermion_basis.global_index("L", 2)
    N = fermion_basis.global_index("F", 1)

    eta_plus_R = scalar_model.block("eta").local_to_global(1)
    eta_plus_I = scalar_model.block("eta").local_to_global(2)
    eta_zero_R = scalar_model.block("eta").local_to_global(3)
    eta_zero_I = scalar_model.block("eta").local_to_global(4)

    assert sp.simplify(tensor[(nu, N, eta_zero_R)] - h / SQRT2) == 0
    assert sp.simplify(tensor[(nu, N, eta_zero_I)] - sp.I * h / SQRT2) == 0

    assert sp.simplify(tensor[(e, N, eta_plus_R)] + h / SQRT2) == 0
    assert sp.simplify(tensor[(e, N, eta_plus_I)] + sp.I * h / SQRT2) == 0

    # y_ija is symmetric in the two Weyl slots under the convention used by
    # the master equation.
    assert tensor[(N, nu, eta_zero_R)] == tensor[(nu, N, eta_zero_R)]
    assert tensor[(N, e, eta_plus_R)] == tensor[(e, N, eta_plus_R)]


def _self_check() -> None:
    _self_check_duality()
    _self_check_scotogenic_epsilon()


if __name__ == "__main__":
    _self_check()
    print("T3YukawaAdapter self-check passed.")
