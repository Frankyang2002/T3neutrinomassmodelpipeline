from __future__ import annotations

from fractions import Fraction
from typing import Callable, Mapping

import sympy as sp

from RGE.general.GaugeGenerators import (
    global_su2_generators,
    global_u1_generator,
)
from RGE.general.RGECommon import (
    C,
    QuarticComponent,
    _validate_output_component,
)
from RGE.general.RGEModel import RGEModel

# Component-level and partial Weinberg-RGE utilities.

LAMBDA = sp.IndexedBase("lambda")


def symbolic_quartic_component(
    a: int,
    b: int,
    c: int,
    d: int,
    lambda_tensor=LAMBDA,
) -> sp.Expr:
    """Return a formal lambda_abcd component without assuming a UV coupling basis."""

    return lambda_tensor[a, b, c, d]


def mapping_quartic_component(
    components: Mapping[tuple[int, int, int, int], sp.Expr],
) -> QuarticComponent:
    """Build a symmetric lambda_abcd lookup from an explicit component mapping."""

    normalized = {
        tuple(sorted(key)): sp.sympify(value)
        for key, value in components.items()
    }

    def component(a: int, b: int, c: int, d: int) -> sp.Expr:
        return normalized.get(tuple(sorted((a, b, c, d))), sp.S.Zero)

    return component


def _gauge_generator_contribution(
    theta: sp.Matrix,
    gauge_coupling: sp.Symbol,
    i: int,
    j: int,
    a_index: int,
    b_index: int,
    real_dimension: int,
    coefficient,
    numerical_coefficient,
) -> sp.Expr:
    """Compute a generic symmetrized scalar-generator contribution."""

    result = sp.S.Zero

    for c_index in range(real_dimension):
        for d_index in range(real_dimension):
            c = c_index + 1
            d = d_index + 1

            symmetrized_generators = (
                theta[a_index, c_index] * theta[b_index, d_index]
                + theta[b_index, c_index] * theta[a_index, d_index]
            )

            result += (
                numerical_coefficient
                * gauge_coupling**2
                * symmetrized_generators
                * coefficient[i, j, c, d]
            )

    return sp.simplify(result)


def first_gauge_component(
    model: RGEModel,
    output_component: tuple[int, int, int, int],
    coefficient=C,
) -> sp.Expr:
    r"""Return the gauge-generator contribution already checked in the legacy code.

    Implements

      - sum_alpha g_alpha^2 sum_cd
        (theta^alpha_ac theta^alpha_bd + theta^alpha_bc theta^alpha_ad) C_ijcd.

    SU(3) is absent because the scalar multiplets considered here are colour singlets.
    """

    _validate_output_component(model, output_component)
    i, j, a, b = output_component
    a_index, b_index = a - 1, b - 1
    n = model.total_real_scalar_dimension

    result = sp.S.Zero

    for theta in global_su2_generators(model):
        result += _gauge_generator_contribution(
            theta=theta,
            gauge_coupling=g2,
            i=i,
            j=j,
            a_index=a_index,
            b_index=b_index,
            real_dimension=n,
            coefficient=coefficient,
            numerical_coefficient=-1,
        )

    result += _gauge_generator_contribution(
        theta=global_u1_generator(model),
        gauge_coupling=g1,
        i=i,
        j=j,
        a_index=a_index,
        b_index=b_index,
        real_dimension=n,
        coefficient=coefficient,
        numerical_coefficient=-1,
    )

    return sp.simplify(result)


def first_scalar_component(
    model: RGEModel,
    output_component: tuple[int, int, int, int],
    quartic_component: QuarticComponent = symbolic_quartic_component,
    coefficient=C,
) -> sp.Expr:
    r"""Return sum_cd lambda_abcd C_ijcd in the global real-scalar basis."""

    _validate_output_component(model, output_component)

    i, j, a, b = output_component
    n = model.total_real_scalar_dimension
    result = sp.S.Zero

    for c in range(1, n + 1):
        for d in range(1, n + 1):
            quartic = quartic_component(a, b, c, d)

            if quartic != 0:
                result += quartic * coefficient[i, j, c, d]

    return sp.simplify(result)


RGEContribution = Callable[..., sp.Expr]


def calculate_partial_rge(
    model: RGEModel,
    output_component: tuple[int, int, int, int],
    quartic_component: QuarticComponent = symbolic_quartic_component,
    coefficient=C,
) -> dict[str, sp.Expr]:
    """Evaluate the generalized contributions implemented so far.

    This is deliberately called *partial* because the full master RGE has not yet
    been implemented in this new file.
    """

    contributions = {
        "first_gauge_component": first_gauge_component(
            model=model,
            output_component=output_component,
            coefficient=coefficient,
        ),
        "first_scalar_component": first_scalar_component(
            model=model,
            output_component=output_component,
            quartic_component=quartic_component,
            coefficient=coefficient,
        ),
    }

    contributions["total"] = sp.simplify(
        sum(contributions.values(), sp.S.Zero)
    )

    return contributions
