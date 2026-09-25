"""One-generation SMEFT Weinberg RGE benchmark.

This module owns the analytic/tensor RGE calculation only. File parsing,
stage output, and pipeline bookkeeping are handled elsewhere.
"""

from __future__ import annotations

import sympy as sp

from RGE.general.FermionBasis import build_gauge_sectors
from RGE.general.GaugeGenerators import g2
from RGE.general.WilsonTensorRGE import (
    WilsonRGEInputs,
    calculate_complete_wilson_tensor_rge,
)
from RGE.matching.WeinbergTensorAdapter import (
    SparseWilsonLookup,
    build_sm_eft,
    build_sm_yukawa,
    build_weinberg_wilson_tensor,
    validate_weinberg_tensor_symmetry,
    yd,
    ye,
    yu,
)


lambdaH = sp.Symbol("lambdaH")


def higgs_quartic(a: int, b: int, c: int, d: int) -> sp.Expr:
    """Return lambda_abcd for V=(lambdaH/2)(H^dagger H)^2."""

    delta = lambda x, y: sp.Integer(1 if x == y else 0)

    return lambdaH * (
        delta(a, b) * delta(c, d)
        + delta(a, c) * delta(b, d)
        + delta(a, d) * delta(b, c)
    )


def calculate_one_generation_weinberg_benchmark(kappa: sp.Expr) -> dict:
    """Evaluate the one-generation Weinberg coefficient with the general tensor RGE."""

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

    inputs = WilsonRGEInputs(
        fermion_dimension=fermion_basis.dimension,
        yukawa=yukawa,
        quartic=higgs_quartic,
        gauge_sectors=build_gauge_sectors(
            scalar_model,
            fermion_basis,
        ),
    )

    contributions = calculate_complete_wilson_tensor_rge(
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
