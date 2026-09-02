from __future__ import annotations

"""
Combined one-generation SM Weinberg-RGE benchmark.

This is the final Standard-Model regression for the generic Eq. (4.85) engine.

Conventions
-----------
Higgs quartic:
    V(H) = (lambdaH / 2) (H^\dagger H)^2

One-generation Yukawas:
    ye : charged lepton
    yu : up-type quark
    yd : down-type quark

Colour is represented by three independent copies of uC/dC and of the quark
doublet components.  SU(3) gauge interactions are unnecessary for this test,
because the Weinberg operator contains only colour singlets and the known
one-loop beta_kappa has no direct g3 term.

Expected result:
    beta_kappa / kappa
      = -3 g2^2
        + 2 lambdaH
        + 6 |yu|^2
        + 6 |yd|^2
        - |ye|^2

The charged-lepton result includes both the trace contribution and the
flavour-dependent external-lepton contribution.
"""

import sympy as sp

from RGE.general.GeneralWeinbergRGEGenerator import (
    ComplexScalar,
    MasterRGEInputs,
    RGEModel,
    calculate_complete_master_rge,
    g1,
    g2,
)
from RGE.general.T3RGETensors import (
    FermionBasis,
    WeylFermion,
    build_gauge_sectors,
)
from RGE.general.WeinbergWilsonAdapter import (
    build_weinberg_wilson_tensor,
    validate_weinberg_tensor_symmetry,
)


kappa = sp.Symbol("kappa")
lambdaH = sp.Symbol("lambdaH")
ye, yu, yd = sp.symbols("ye yu yd")


class SparseWilsonLookup:
    def __init__(self, components):
        self.components = components

    def __getitem__(self, key):
        return sp.sympify(
            self.components.get(tuple(key), sp.S.Zero)
        )


def higgs_quartic(a: int, b: int, c: int, d: int) -> sp.Expr:
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
    entries = {}

    h = scalar_model.block("H")
    hp_r = h.local_to_global(1)
    hp_i = h.local_to_global(2)
    h0_r = h.local_to_global(3)
    h0_i = h.local_to_global(4)

    root2 = sp.sqrt(2)

    # Charged lepton: eC H^\dagger L.
    nu = fermion_basis.global_index("L", 1)
    e = fermion_basis.global_index("L", 2)
    ec = fermion_basis.global_index("eC", 1)

    add_symmetric_yukawa(entries, nu, ec, hp_r, ye / root2)
    add_symmetric_yukawa(entries, nu, ec, hp_i, -sp.I * ye / root2)
    add_symmetric_yukawa(entries, e, ec, h0_r, ye / root2)
    add_symmetric_yukawa(entries, e, ec, h0_i, -sp.I * ye / root2)

    # Three independent colour copies.  Only multiplicity matters here.
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
        return sp.sympify(entries.get((i, j, a), sp.S.Zero))

    return yukawa


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

    fermion_basis = FermionBasis(fermions=tuple(fermions))

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

    beta_ratio = sp.factor(
        sp.simplify(result["total"] / c_value)
    )

    expected = (
        -3 * g2**2
        + 2 * lambdaH
        + 6 * yu * sp.conjugate(yu)
        + 6 * yd * sp.conjugate(yd)
        - ye * sp.conjugate(ye)
    )

    difference = sp.factor(
        sp.expand(beta_ratio - expected)
    )

    print("=" * 72)
    print("COMBINED ONE-GENERATION SM WEINBERG BENCHMARK")
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
