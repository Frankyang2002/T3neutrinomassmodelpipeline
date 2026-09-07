from __future__ import annotations

import numpy as np
import sympy as sp

from RGE.general.RGECommon import C
from RGE.general.RGEModel import (
    ComplexScalar,
    RGEModel,
)

# Conversion utilities between the real-scalar Wilson tensor C_{ijab}
# and the complex neutral-scalar coefficient kappa used for the Weinberg operator.

def su2_weights(dimension: int) -> tuple[sp.Rational, ...]:
    """Return T3 weights ordered consistently with su2_complex_generators."""

    j = sp.Rational(dimension - 1, 2)
    return tuple(j - position for position in range(dimension))


def neutral_complex_component_position(scalar: ComplexScalar) -> int | None:
    """Return the 1-based complex-component position with Q=T3+Y=0."""

    target_t3 = -scalar.hypercharge
    weights = su2_weights(scalar.su2_dimension)

    try:
        return weights.index(target_t3) + 1
    except ValueError:
        return None


def neutral_real_pair(model: RGEModel, scalar_name: str) -> tuple[int, int] | None:
    """Return the global (R,I) indices of a scalar's neutral complex component."""

    block = model.block(scalar_name)
    complex_position = neutral_complex_component_position(block.scalar)

    if complex_position is None:
        return None

    local_r = 2 * complex_position - 1
    local_i = 2 * complex_position

    return (
        block.local_to_global(local_r),
        block.local_to_global(local_i),
    )


def neutral_k_symbols(model: RGEModel) -> dict[str, sp.Symbol]:
    """Create one symbolic K coefficient for every scalar with a neutral component."""

    result: dict[str, sp.Symbol] = {}

    for scalar in model.scalars:
        if neutral_real_pair(model, scalar.name) is not None:
            result[scalar.name] = sp.Symbol(f"K_{scalar.name}")

    return result


def c_to_k_substitutions(
    model: RGEModel,
    fermion_indices: tuple[int, int],
    coefficient=C,
    k_symbols: Mapping[str, sp.Symbol] | None = None,
) -> dict[sp.Expr, sp.Expr]:
    """Return neutral-neutral C -> K substitutions for every scalar multiplet."""

    i, j = fermion_indices
    symbols = dict(k_symbols or neutral_k_symbols(model))
    substitutions: dict[sp.Expr, sp.Expr] = {}

    n = model.total_real_scalar_dimension

    for a in range(1, n + 1):
        for b in range(1, n + 1):
            substitutions[coefficient[i, j, a, b]] = sp.S.Zero

    for scalar_name, k_symbol in symbols.items():
        pair = neutral_real_pair(model, scalar_name)

        if pair is None:
            continue

        r, im = pair

        substitutions[coefficient[i, j, r, r]] = k_symbol / 2
        substitutions[coefficient[i, j, r, im]] = sp.I * k_symbol / 2
        substitutions[coefficient[i, j, im, r]] = sp.I * k_symbol / 2
        substitutions[coefficient[i, j, im, im]] = -k_symbol / 2

    return substitutions


def convert_c_expression_to_k(
    expression: sp.Expr,
    model: RGEModel,
    fermion_indices: tuple[int, int],
    coefficient=C,
    k_symbols: Mapping[str, sp.Symbol] | None = None,
) -> sp.Expr:
    """Convert supported neutral C components into model-generated K symbols."""

    substitutions = c_to_k_substitutions(
        model=model,
        fermion_indices=fermion_indices,
        coefficient=coefficient,
        k_symbols=k_symbols,
    )

    return sp.simplify(sp.expand(expression.xreplace(substitutions)))
