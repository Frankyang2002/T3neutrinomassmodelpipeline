from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations_with_replacement
from typing import Iterable, Mapping, Sequence

import sympy as sp

from RGE.general.RGEModel import RGEModel

# Quartic scalar components and conversion to the real-scalar tensor.

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
    from tests.python.Weinberg.T3RGETensors import SQRT2
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


def quartic_components_from_exchange(
    data: Mapping,
) -> tuple[ComplexQuarticComponent, ...]:
    from tests.python.Weinberg.T3RGETensors import parse_exact_expression
    """Parse exact complex-basis quartic terms from the Wolfram exchange."""

    result: list[ComplexQuarticComponent] = []

    for entry in data.get("quartic_components", []):
        factors = tuple(
            ComplexScalarFactor(
                scalar_name=str(factor["scalar_name"]),
                component=int(factor["component"]),
                conjugated=bool(factor.get("conjugated", False)),
            )
            for factor in entry["factors"]
        )

        if len(factors) != 4:
            raise ValueError("Every quartic exchange term must have four factors.")

        result.append(
            ComplexQuarticComponent(
                coefficient=parse_exact_expression(entry["coefficient"]),
                factors=factors,
            )
        )

    return tuple(result)
