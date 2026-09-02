from __future__ import annotations

"""
Adapter between generalized T3 model data and the tensors used by
GeneralWeinbergRGEGenerator_complete.py.

This file deliberately does not parse pretty-printed Matchete Lagrangians.
Instead, the Wolfram/model-building side should export machine-readable
component tensors for the accepted UV interactions.  This module then:

- constructs the global real-scalar basis;
- constructs a global left-handed Weyl-fermion basis;
- embeds scalar and fermion gauge generators;
- converts complex-scalar Yukawa components into y_ija;
- expands complex quartic interactions into the fully symmetric real tensor
  lambda_abcd;
- stores/makes available the matched Wilson tensor C_ijab.

Important normalization convention
----------------------------------
The quartic conversion assumes the real-scalar convention used by Eq. (4.85),

    V_4 = (1/4!) lambda_abcd phi_a phi_b phi_c phi_d,

so lambda_abcd is obtained by four derivatives of the quartic polynomial.

For Yukawas, this adapter only converts the complex scalar into its real
(R,I) components.  The supplied complex Yukawa coefficient must already use
the same fermion ordering and overall normalization as y_ija in Eq. (4.85).
"""

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import sympy as sp

from RGE.general.GeneralWeinbergRGEGenerator import (
    C,
    GaugeSector,
    MasterRGEInputs,
    RGEModel,
    ComplexScalar,
    g1,
    g2,
    global_su2_generators,
    global_u1_generator,
    su2_complex_generators,
)


SQRT2 = sp.sqrt(2)


# -----------------------------------------------------------------------------
# Exact-expression parsing
# -----------------------------------------------------------------------------


def parse_exact_expression(value) -> sp.Expr:
    """Convert a JSON number/string into a SymPy expression.

    The exporter should prefer simple exact strings such as:
      "1/2", "-sqrt(3)/2", "I/sqrt(2)", "y1[1,2]".
    """

    if isinstance(value, (int, float)):
        return sp.sympify(value)

    if not isinstance(value, str):
        raise TypeError(f"Expected an exact-expression string, got {type(value)!r}.")

    # Keep the parser intentionally small and transparent.
    translated = (
        value.replace("Sqrt[", "sqrt(")
        .replace("I", "I")
    )

    # Handle the common simple Mathematica Sqrt[...] form.
    if "Sqrt[" in value:
        translated = value.replace("Sqrt[", "sqrt(").replace("]", ")")

    return sp.sympify(
        translated,
        locals={
            "I": sp.I,
            "sqrt": sp.sqrt,
        },
    )


# -----------------------------------------------------------------------------
# Fermion basis
# -----------------------------------------------------------------------------


@dataclass(frozen=True)
class WeylFermion:
    """One left-handed Weyl multiplet used by the y_ija/C_ijab convention."""

    name: str
    su2_dimension: int
    hypercharge: sp.Expr
    multiplicity: int = 1
    conjugated_representation: bool = False

    def __post_init__(self) -> None:
        if self.su2_dimension < 1:
            raise ValueError("SU(2) representation dimension must be positive.")
        if self.multiplicity < 1:
            raise ValueError("Fermion multiplicity must be positive.")


@dataclass(frozen=True)
class FermionBasisBlock:
    """One fermion multiplet/flavour copy inside the global Weyl basis."""

    fermion: WeylFermion
    copy: int
    first: int
    last: int

    @property
    def indices(self) -> range:
        return range(self.first, self.last + 1)


@dataclass
class FermionBasis:
    """Global basis and gauge generators for all left-handed Weyl fermions."""

    fermions: tuple[WeylFermion, ...]

    def __post_init__(self) -> None:
        start = 1
        blocks: list[FermionBasisBlock] = []

        for fermion in self.fermions:
            for copy in range(1, fermion.multiplicity + 1):
                stop = start + fermion.su2_dimension - 1
                blocks.append(
                    FermionBasisBlock(
                        fermion=fermion,
                        copy=copy,
                        first=start,
                        last=stop,
                    )
                )
                start = stop + 1

        self.blocks = tuple(blocks)
        self.dimension = start - 1

    def blocks_named(self, name: str) -> tuple[FermionBasisBlock, ...]:
        result = tuple(block for block in self.blocks if block.fermion.name == name)
        if not result:
            raise KeyError(f"Unknown fermion multiplet {name!r}.")
        return result

    def block(self, name: str, copy: int = 1) -> FermionBasisBlock:
        for block in self.blocks:
            if block.fermion.name == name and block.copy == copy:
                return block
        raise KeyError(f"No fermion block {name!r}, copy={copy}.")

    def global_index(self, name: str, component: int, copy: int = 1) -> int:
        block = self.block(name, copy)
        if not 1 <= component <= block.fermion.su2_dimension:
            raise IndexError(
                f"{name} component must lie in "
                f"1,...,{block.fermion.su2_dimension}."
            )
        return block.first + component - 1


def _fermion_local_su2_generators(
    fermion: WeylFermion,
) -> tuple[sp.Matrix, sp.Matrix, sp.Matrix]:
    """Return generators for the requested Weyl representation."""

    generators = su2_complex_generators(fermion.su2_dimension)

    if fermion.conjugated_representation:
        # For the conjugate representation: T^A -> -(T^A)^*.
        generators = tuple(-generator.conjugate() for generator in generators)

    return generators


def fermion_global_su2_generators(
    basis: FermionBasis,
) -> tuple[sp.Matrix, sp.Matrix, sp.Matrix]:
    """Embed SU(2)_L generators into the full Weyl-fermion basis."""

    result = [sp.zeros(basis.dimension) for _ in range(3)]

    for block in basis.blocks:
        local = _fermion_local_su2_generators(block.fermion)
        offset = block.first - 1

        for generator_number in range(3):
            for row in range(block.fermion.su2_dimension):
                for column in range(block.fermion.su2_dimension):
                    result[generator_number][
                        offset + row,
                        offset + column,
                    ] = local[generator_number][row, column]

    return tuple(sp.simplify(generator) for generator in result)


def fermion_global_u1_generator(basis: FermionBasis) -> sp.Matrix:
    """Return the U(1)_Y generator in the global Weyl-fermion basis."""

    result = sp.zeros(basis.dimension)

    for block in basis.blocks:
        hypercharge = sp.sympify(block.fermion.hypercharge)
        if block.fermion.conjugated_representation:
            hypercharge = -hypercharge

        for index in block.indices:
            result[index - 1, index - 1] = hypercharge

    return result


def build_gauge_sectors(
    scalar_model: RGEModel,
    fermion_basis: FermionBasis,
) -> tuple[GaugeSector, GaugeSector]:
    """Build the SU(2)_L and U(1)_Y sectors required by Eq. (4.85)."""

    return (
        GaugeSector(
            coupling=g2,
            scalar_generators=global_su2_generators(scalar_model),
            fermion_generators=fermion_global_su2_generators(fermion_basis),
        ),
        GaugeSector(
            coupling=g1,
            scalar_generators=(global_u1_generator(scalar_model),),
            fermion_generators=(fermion_global_u1_generator(fermion_basis),),
        ),
    )


# -----------------------------------------------------------------------------
# Scalar real-basis variables
# -----------------------------------------------------------------------------


def scalar_real_symbols(model: RGEModel) -> tuple[sp.Symbol, ...]:
    """Create one real symbol phi_1,...,phi_N for the global scalar basis."""

    return tuple(
        sp.Symbol(f"phi_{index}", real=True)
        for index in range(1, model.total_real_scalar_dimension + 1)
    )


def complex_scalar_component_expression(
    model: RGEModel,
    scalar_name: str,
    component: int,
    conjugated: bool,
    real_symbols: Sequence[sp.Symbol] | None = None,
) -> sp.Expr:
    r"""Return z_m=(R_m+i I_m)/sqrt(2), or its conjugate, in the global basis."""

    block = model.block(scalar_name)

    if not 1 <= component <= block.scalar.su2_dimension:
        raise IndexError(
            f"{scalar_name} complex component must lie in "
            f"1,...,{block.scalar.su2_dimension}."
        )

    symbols = tuple(real_symbols or scalar_real_symbols(model))

    local_r = 2 * component - 1
    local_i = 2 * component
    global_r = block.local_to_global(local_r)
    global_i = block.local_to_global(local_i)

    sign = -1 if conjugated else 1

    return (
        symbols[global_r - 1]
        + sign * sp.I * symbols[global_i - 1]
    ) / SQRT2


# -----------------------------------------------------------------------------
# Yukawa tensor y_ija
# -----------------------------------------------------------------------------


@dataclass(frozen=True)
class ComplexYukawaComponent:
    """One component coefficient multiplying psi_i psi_j times a complex scalar."""

    fermion_i: int
    fermion_j: int
    scalar_name: str
    scalar_component: int
    scalar_conjugated: bool
    coefficient: sp.Expr


def build_real_yukawa_tensor(
    scalar_model: RGEModel,
    fermion_basis: FermionBasis,
    components: Iterable[ComplexYukawaComponent],
    symmetrize_fermions: bool = False,
) -> dict[tuple[int, int, int], sp.Expr]:
    """Convert complex-scalar Yukawa components into y_ija in the real basis.

    For z=(R+iI)/sqrt(2):
      y_R = Y/sqrt(2),
      y_I = +i Y/sqrt(2).

    For z^*=(R-iI)/sqrt(2):
      y_R = Y/sqrt(2),
      y_I = -i Y/sqrt(2).

    The function does not add any extra 1/2 or sign associated with a chosen
    two-component-fermion Lagrangian convention.  The exported complex
    coefficient must already match Eq. (4.85)'s y convention.
    """

    result: dict[tuple[int, int, int], sp.Expr] = {}

    for item in components:
        if not 1 <= item.fermion_i <= fermion_basis.dimension:
            raise IndexError("fermion_i is outside the global fermion basis.")
        if not 1 <= item.fermion_j <= fermion_basis.dimension:
            raise IndexError("fermion_j is outside the global fermion basis.")

        block = scalar_model.block(item.scalar_name)
        if not 1 <= item.scalar_component <= block.scalar.su2_dimension:
            raise IndexError("scalar_component is outside its SU(2) multiplet.")

        local_r = 2 * item.scalar_component - 1
        local_i = 2 * item.scalar_component
        r = block.local_to_global(local_r)
        im = block.local_to_global(local_i)

        coefficient = sp.sympify(item.coefficient)
        imaginary_sign = -1 if item.scalar_conjugated else 1

        contributions = {
            (item.fermion_i, item.fermion_j, r): coefficient / SQRT2,
            (
                item.fermion_i,
                item.fermion_j,
                im,
            ): imaginary_sign * sp.I * coefficient / SQRT2,
        }

        if symmetrize_fermions and item.fermion_i != item.fermion_j:
            contributions.update(
                {
                    (item.fermion_j, item.fermion_i, r): coefficient / SQRT2,
                    (
                        item.fermion_j,
                        item.fermion_i,
                        im,
                    ): imaginary_sign * sp.I * coefficient / SQRT2,
                }
            )

        for key, value in contributions.items():
            result[key] = sp.simplify(result.get(key, sp.S.Zero) + value)

    return result


def yukawa_component_function(
    components: Mapping[tuple[int, int, int], sp.Expr],
):
    """Return the y(i,j,a) callable expected by MasterRGEInputs."""

    def y(i: int, j: int, a: int) -> sp.Expr:
        return sp.sympify(components.get((i, j, a), sp.S.Zero))

    return y


# -----------------------------------------------------------------------------
# Quartic tensor lambda_abcd
# -----------------------------------------------------------------------------


@dataclass(frozen=True)
class ComplexScalarFactor:
    """One complex scalar field appearing in a quartic component term."""

    scalar_name: str
    component: int
    conjugated: bool = False


@dataclass(frozen=True)
class ComplexQuarticComponent:
    """One component of a quartic interaction before conversion to real fields."""

    coefficient: sp.Expr
    factors: tuple[
        ComplexScalarFactor,
        ComplexScalarFactor,
        ComplexScalarFactor,
        ComplexScalarFactor,
    ]


def quartic_polynomial_from_components(
    model: RGEModel,
    components: Iterable[ComplexQuarticComponent],
) -> tuple[sp.Expr, tuple[sp.Symbol, ...]]:
    """Expand exported complex quartic components into real scalar fields."""

    symbols = scalar_real_symbols(model)
    polynomial = sp.S.Zero

    for term in components:
        product = sp.sympify(term.coefficient)

        for factor in term.factors:
            product *= complex_scalar_component_expression(
                model=model,
                scalar_name=factor.scalar_name,
                component=factor.component,
                conjugated=factor.conjugated,
                real_symbols=symbols,
            )

        polynomial += product

    return sp.expand(polynomial), symbols


def quartic_tensor_from_components(
    model: RGEModel,
    components: Iterable[ComplexQuarticComponent],
    simplify: bool = True,
) -> dict[tuple[int, int, int, int], sp.Expr]:
    r"""Return the symmetric real tensor lambda_abcd.

    The tensor is defined by

        lambda_abcd =
          d^4 V_4 / (d phi_a d phi_b d phi_c d phi_d),

    which is equivalent to V_4=(1/4!) lambda_abcd phi_a phi_b phi_c phi_d.
    Only sorted index keys are stored because lambda_abcd is fully symmetric.
    """

    polynomial, symbols = quartic_polynomial_from_components(model, components)
    n = model.total_real_scalar_dimension
    result: dict[tuple[int, int, int, int], sp.Expr] = {}

    # combinations_with_replacement avoids N^4 duplicate differentiation.
    from itertools import combinations_with_replacement

    for key in combinations_with_replacement(range(1, n + 1), 4):
        derivative = polynomial

        for index in key:
            derivative = sp.diff(derivative, symbols[index - 1])

        value = sp.simplify(derivative) if simplify else derivative

        if value != 0:
            result[key] = value

    return result


def quartic_component_function(
    components: Mapping[tuple[int, int, int, int], sp.Expr],
):
    """Return lambda(a,b,c,d) for MasterRGEInputs."""

    normalized = {
        tuple(sorted(key)): sp.sympify(value)
        for key, value in components.items()
    }

    def quartic(a: int, b: int, c: int, d: int) -> sp.Expr:
        return normalized.get(tuple(sorted((a, b, c, d))), sp.S.Zero)

    return quartic


# -----------------------------------------------------------------------------
# Wilson tensor C_ijab
# -----------------------------------------------------------------------------


class SparseWilsonLookup:
    def __init__(
        self,
        components: Mapping[tuple[int, int, int, int], sp.Expr],
    ):
        self.components = components

    def __getitem__(
        self,
        key: tuple[int, int, int, int],
    ) -> sp.Expr:
        return sp.sympify(
            self.components.get(key, sp.S.Zero)
        )

    def __call__(
        self,
        i: int,
        j: int,
        a: int,
        b: int,
    ) -> sp.Expr:
        return self[i, j, a, b]


def wilson_component_function(
    components: Mapping[tuple[int, int, int, int], sp.Expr],
):
    """Return a sparse Wilson-coefficient lookup callable."""

    return SparseWilsonLookup(components)


# -----------------------------------------------------------------------------
# Final assembly
# -----------------------------------------------------------------------------


@dataclass
class T3RGETensors:
    """All tensors/bases needed to call the generalized Eq. (4.85) engine."""

    scalar_model: RGEModel
    fermion_basis: FermionBasis
    yukawa_components: dict[tuple[int, int, int], sp.Expr]
    quartic_components: dict[tuple[int, int, int, int], sp.Expr]
    wilson_components: dict[tuple[int, int, int, int], sp.Expr]

    def master_inputs(self) -> MasterRGEInputs:
        """Build MasterRGEInputs from the converted UV tensors."""

        return MasterRGEInputs(
            fermion_dimension=self.fermion_basis.dimension,
            yukawa=yukawa_component_function(self.yukawa_components),
            quartic=quartic_component_function(self.quartic_components),
            gauge_sectors=build_gauge_sectors(
                self.scalar_model,
                self.fermion_basis,
            ),
        )


# -----------------------------------------------------------------------------
# Machine-readable exchange format
# -----------------------------------------------------------------------------


def load_component_exchange(path: str | Path) -> dict:
    """Load the JSON component exchange produced by the Wolfram adapter."""

    return json.loads(Path(path).read_text(encoding="utf-8"))


def scalar_model_from_exchange(data: Mapping) -> RGEModel:
    """Construct the scalar RGE model from exchange JSON."""

    scalars = tuple(
        ComplexScalar(
            name=item["name"],
            su2_dimension=int(item["su2_dimension"]),
            hypercharge=parse_exact_expression(item["hypercharge"]),
        )
        for item in data["scalars"]
    )

    return RGEModel(scalars=scalars)


def fermion_basis_from_exchange(data: Mapping) -> FermionBasis:
    """Construct the global Weyl basis from exchange JSON."""

    fermions = tuple(
        WeylFermion(
            name=item["name"],
            su2_dimension=int(item["su2_dimension"]),
            hypercharge=parse_exact_expression(item["hypercharge"]),
            multiplicity=int(item.get("multiplicity", 1)),
            conjugated_representation=bool(
                item.get("conjugated_representation", False)
            ),
        )
        for item in data["fermions"]
    )

    return FermionBasis(fermions=fermions)


# -----------------------------------------------------------------------------
# Regression/self-checks
# -----------------------------------------------------------------------------


def _self_check_quartic_normalization() -> None:
    """Recover the legacy doublet self-quartic normalization.

    For V=(lambda/2)(Phi^dagger Phi)^2, the real tensor must satisfy
      lambda_aaaa = 3 lambda,
      lambda_aabb = lambda   (a != b).
    """

    lam = sp.Symbol("lam", real=True)
    model = RGEModel(
        scalars=(
            ComplexScalar(
                name="Phi",
                su2_dimension=2,
                hypercharge=sp.Rational(1, 2),
            ),
        )
    )

    terms: list[ComplexQuarticComponent] = []

    # (lambda/2) sum_mn Phi_m^* Phi_m Phi_n^* Phi_n
    for m in range(1, 3):
        for n in range(1, 3):
            terms.append(
                ComplexQuarticComponent(
                    coefficient=lam / 2,
                    factors=(
                        ComplexScalarFactor("Phi", m, True),
                        ComplexScalarFactor("Phi", m, False),
                        ComplexScalarFactor("Phi", n, True),
                        ComplexScalarFactor("Phi", n, False),
                    ),
                )
            )

    tensor = quartic_tensor_from_components(model, terms)

    assert sp.simplify(tensor[(1, 1, 1, 1)] - 3 * lam) == 0
    assert sp.simplify(tensor[(1, 1, 2, 2)] - lam) == 0
    assert sp.simplify(tensor[(1, 1, 3, 3)] - lam) == 0


def _self_check_yukawa_conversion() -> None:
    """Check z=(R+iI)/sqrt(2) conversion used by y_ija."""

    y0 = sp.Symbol("y0")
    scalar_model = RGEModel(
        scalars=(
            ComplexScalar("Phi", 1, 0),
        )
    )
    fermion_basis = FermionBasis(
        fermions=(
            WeylFermion("f", 1, 0, multiplicity=2),
        )
    )

    tensor = build_real_yukawa_tensor(
        scalar_model,
        fermion_basis,
        (
            ComplexYukawaComponent(
                fermion_i=1,
                fermion_j=2,
                scalar_name="Phi",
                scalar_component=1,
                scalar_conjugated=False,
                coefficient=y0,
            ),
        ),
    )

    assert sp.simplify(tensor[(1, 2, 1)] - y0 / SQRT2) == 0
    assert sp.simplify(tensor[(1, 2, 2)] - sp.I * y0 / SQRT2) == 0


def _self_check() -> None:
    _self_check_quartic_normalization()
    _self_check_yukawa_conversion()


if __name__ == "__main__":
    _self_check()
    print("T3RGETensors self-check passed.")
