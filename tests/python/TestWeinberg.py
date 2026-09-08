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

from tests.python.GeneralWeinbergRGEGenerator import (
    ComplexScalar,
    MasterRGEInputs,
    RGEModel,
    calculate_complete_master_rge,
    g1,
    g2,
)
from tests.python.Weinberg.T3RGETensors import (
    FermionBasis,
    WeylFermion,
    build_gauge_sectors,
)
from tests.python.Weinberg.WeinbergWilsonAdapter import (
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


def check_sm_complete() -> int:
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


# -----------------------------------------------------------------------------
# Additional Weinberg regressions consolidated from the former standalone files
# -----------------------------------------------------------------------------

import numpy as np

from RGE.running.SMEFTWeinbergFlavorRGE import (
    beta_weinberg_matrix,
    one_generation_reduction,
    symbolic_complex_matrix,
    symbolic_symmetric_matrix,
)
from RGE.running.NumericalWeinbergRGE import (
    SMInitialConditions,
    evolve_weinberg,
    neutrino_mass_matrix,
)
from RGE.stages.NeutrinoMassStage import build_neutrino_mass_matrix, v


def check_full_flavor() -> None:
    ye1, yu1, yd1 = sp.symbols("ye yu yd")

    expected_1g = (
        -3 * g2**2
        + 2 * lambdaH
        + 6 * yu1 * sp.conjugate(yu1)
        + 6 * yd1 * sp.conjugate(yd1)
        - ye1 * sp.conjugate(ye1)
    )
    assert sp.simplify(one_generation_reduction() - expected_1g) == 0

    K = symbolic_symmetric_matrix("k", 3)
    Ye = symbolic_complex_matrix("ye", 3, 3)
    Yu = symbolic_complex_matrix("yu", 3, 3)
    Yd = symbolic_complex_matrix("yd", 3, 3)

    beta = beta_weinberg_matrix(K, Ye, Yu, Yd)
    symmetry = (beta - beta.T).applyfunc(sp.simplify)
    assert all(entry == 0 for entry in symmetry)

    e1, e2, e3 = sp.symbols("e1 e2 e3")
    u1, u2, u3 = sp.symbols("u1 u2 u3")
    d1, d2, d3 = sp.symbols("d1 d2 d3")

    Ye_diag = sp.diag(e1, e2, e3)
    Yu_diag = sp.diag(u1, u2, u3)
    Yd_diag = sp.diag(d1, d2, d3)
    beta_diag = beta_weinberg_matrix(K, Ye_diag, Yu_diag, Yd_diag)

    trace = (
        sum(x * sp.conjugate(x) for x in (e1, e2, e3))
        + 3 * sum(x * sp.conjugate(x) for x in (u1, u2, u3))
        + 3 * sum(x * sp.conjugate(x) for x in (d1, d2, d3))
    )
    universal = 2 * lambdaH - 3 * g2**2 + 2 * trace
    charged = (e1, e2, e3)

    for p in range(3):
        for q in range(3):
            expected = (
                universal
                - sp.Rational(3, 2)
                * (
                    charged[p] * sp.conjugate(charged[p])
                    + charged[q] * sp.conjugate(charged[q])
                )
            ) * K[p, q]
            assert sp.simplify(beta_diag[p, q] - expected) == 0


def check_mass_normalization() -> None:
    k11, k12, k13, k22, k23, k33 = sp.symbols(
        "k11 k12 k13 k22 k23 k33"
    )
    K = sp.Matrix(
        [
            [k11, k12, k13],
            [k12, k22, k23],
            [k13, k23, k33],
        ]
    )
    M = build_neutrino_mass_matrix(K)
    expected = -sp.Rational(1, 2) * v**2 * K

    assert all(sp.simplify(x) == 0 for x in (M - expected))
    assert all(sp.simplify(x) == 0 for x in (M - M.T))


def check_numerical_running() -> None:
    K0 = np.array(
        [
            [1.0e-14, 2.0e-15 + 1.0e-15j, 3.0e-15],
            [2.0e-15 + 1.0e-15j, 1.2e-14, -1.0e-15j],
            [3.0e-15, -1.0e-15j, 8.0e-15],
        ],
        dtype=complex,
    )

    initial = SMInitialConditions(
        gY=0.36,
        g2=0.64,
        g3=1.00,
        lambdaH=0.26,
        ye=np.array([2.9e-6, 6.1e-4, 1.02e-2]),
        yu=np.array([7.0e-6, 3.5e-3, 0.85]),
        yd=np.array([1.5e-5, 3.0e-4, 1.6e-2]),
        K=K0,
    )

    high, low = 1.0e6, 1.0e2
    down = evolve_weinberg(initial, high, low)

    assert np.allclose(down.K, down.K.T, rtol=1e-10, atol=1e-20)

    back = evolve_weinberg(
        SMInitialConditions(
            gY=down.gY,
            g2=down.g2,
            g3=down.g3,
            lambdaH=down.lambdaH,
            ye=down.ye,
            yu=down.yu,
            yd=down.yd,
            K=down.K,
        ),
        low,
        high,
    )

    relative_error = np.max(
        np.abs(back.K - K0) / np.maximum(np.abs(K0), 1e-30)
    )
    assert relative_error <= 1e-5

    # Also exercise the production K -> M_nu conversion.
    mass = neutrino_mass_matrix(down.K)
    assert mass.shape == (3, 3)


def main() -> int:
    assert check_sm_complete() == 0
    check_full_flavor()
    check_mass_normalization()
    check_numerical_running()

    print("PASS: complete one-generation SM Weinberg benchmark")
    print("PASS: three-flavor matrix RGE")
    print("PASS: M_nu = -(v^2/2) K normalization")
    print("PASS: numerical Weinberg running")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
