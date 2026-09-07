from __future__ import annotations

"""
Pure-SM Higgs-quartic validation of the generic psi^2 phi^2 Eq. (4.85).

We use the inert-doublet/scotogenic quartic normalization

    V(H) = (lambdaH / 2) (H^\dagger H)^2.

Writing

    H = (H+, H0),

with each complex component split as

    Phi = (R + i I)/sqrt(2),

gives

    H^\dagger H = 1/2 sum_a phi_a^2.

The real-scalar quartic tensor is defined by

    V_4 = (1/4!) lambda_abcd phi_a phi_b phi_c phi_d.

For the convention above,

    lambda_abcd
      = lambdaH (
          delta_ab delta_cd
        + delta_ac delta_bd
        + delta_ad delta_bc
        ).

With gauge and Yukawa interactions switched off, the SM Weinberg operator
must therefore satisfy

    beta_kappa / kappa = 2 lambdaH.

This test checks the production Eq. (4.85) implementation directly.
"""

import sympy as sp

from tests.python.GeneralWeinbergRGEGenerator import (
    ComplexScalar,
    MasterRGEInputs,
    RGEModel,
    calculate_complete_master_rge,
)
from tests.python.Weinberg.T3RGETensors import (
    FermionBasis,
    WeylFermion,
)
from tests.python.Weinberg.WeinbergWilsonAdapter import (
    build_weinberg_wilson_tensor,
    validate_weinberg_tensor_symmetry,
)


lambdaH = sp.Symbol("lambdaH")
kappa = sp.Symbol("kappa")


class SparseWilsonLookup:
    def __init__(self, components):
        self.components = components

    def __getitem__(self, key):
        return sp.sympify(
            self.components.get(tuple(key), sp.S.Zero)
        )


def zero_yukawa(i: int, j: int, a: int) -> sp.Expr:
    return sp.S.Zero


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


def main() -> int:
    scalar_model = RGEModel(
        (
            ComplexScalar(
                name="H",
                su2_dimension=2,
                hypercharge=sp.Rational(1, 2),
            ),
        )
    )

    fermion_basis = FermionBasis(
        fermions=(
            WeylFermion(
                name="L",
                su2_dimension=2,
                hypercharge=-sp.Rational(1, 2),
            ),
        )
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

    inputs = MasterRGEInputs(
        fermion_dimension=fermion_basis.dimension,
        yukawa=zero_yukawa,
        quartic=higgs_quartic,
        gauge_sectors=(),
    )

    result = calculate_complete_master_rge(
        model=scalar_model,
        inputs=inputs,
        output_component=component,
        coefficient=C,
    )

    print("=" * 72)
    print("SM HIGGS-QUARTIC WEINBERG BENCHMARK")
    print("=" * 72)

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

    for name in order:
        ratio = sp.factor(
            sp.simplify(result[name] / c_value)
        )
        print(f"{name:34s}: {ratio}")

    beta_ratio = sp.factor(
        sp.simplify(result["total"] / c_value)
    )

    expected = 2 * lambdaH
    difference = sp.factor(
        sp.expand(beta_ratio - expected)
    )

    print()
    print("Expected beta_kappa/kappa :", expected)
    print("Engine beta_kappa/kappa   :", beta_ratio)
    print("Difference                :", difference)
    print()

    if sp.simplify(difference) != 0:
        print("FAIL")
        return 1

    print("PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
