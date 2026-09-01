from __future__ import annotations

"""
Pure-SM charged-lepton Yukawa validation of Eq. (4.85).

Gauge couplings and scalar quartics are switched off.  We keep one lepton
generation and the charged-lepton Yukawa only.

In an all-left-handed two-component basis,

    L = (nu, e),      eC = charge-conjugate right-handed electron,

and the Yukawa interaction is proportional to

    y_e eC (H^\dagger L).

Using

    H+ = (H+_R + i H+_I)/sqrt(2),
    H0 = (H0_R + i H0_I)/sqrt(2),

the complex-conjugated Higgs components appearing in H^\dagger give

    y(nu,eC,H+_R) =  y_e/sqrt(2)
    y(nu,eC,H+_I) = -i y_e/sqrt(2)
    y(e, eC,H0_R) =  y_e/sqrt(2)
    y(e, eC,H0_I) = -i y_e/sqrt(2),

with symmetric fermion slots.

For one generation, the known SM Weinberg-operator result is

    beta_kappa / kappa = - |y_e|^2

when all gauge and quartic terms are disabled.
"""

import sympy as sp

from GeneralWeinbergRGEGenerator import (
    ComplexScalar,
    MasterRGEInputs,
    RGEModel,
    calculate_complete_master_rge,
)
from T3RGETensors import FermionBasis, WeylFermion
from WeinbergWilsonAdapter import (
    build_weinberg_wilson_tensor,
    validate_weinberg_tensor_symmetry,
)


ye = sp.Symbol("ye")
kappa = sp.Symbol("kappa")


class SparseWilsonLookup:
    def __init__(self, components):
        self.components = components

    def __getitem__(self, key):
        return sp.sympify(
            self.components.get(tuple(key), sp.S.Zero)
        )


def zero_quartic(
    a: int,
    b: int,
    c: int,
    d: int,
) -> sp.Expr:
    return sp.S.Zero


def build_charged_lepton_yukawa(
    scalar_model: RGEModel,
    fermion_basis: FermionBasis,
):
    """Return the one-generation SM charged-lepton y_ija tensor."""

    nu = fermion_basis.global_index("L", 1)
    e = fermion_basis.global_index("L", 2)
    ec = fermion_basis.global_index("eC", 1)

    h_block = scalar_model.block("H")

    hp_r = h_block.local_to_global(1)
    hp_i = h_block.local_to_global(2)
    h0_r = h_block.local_to_global(3)
    h0_i = h_block.local_to_global(4)

    root2 = sp.sqrt(2)

    entries = {
        (nu, ec, hp_r): ye / root2,
        (nu, ec, hp_i): -sp.I * ye / root2,
        (e, ec, h0_r): ye / root2,
        (e, ec, h0_i): -sp.I * ye / root2,
    }

    # y_ija is symmetric in its two Weyl-fermion slots.
    symmetric_entries = dict(entries)

    for (i, j, a), value in entries.items():
        symmetric_entries[(j, i, a)] = value

    def yukawa(i: int, j: int, a: int) -> sp.Expr:
        return sp.sympify(
            symmetric_entries.get((i, j, a), sp.S.Zero)
        )

    return yukawa, symmetric_entries


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
            WeylFermion(
                name="eC",
                su2_dimension=1,
                hypercharge=sp.Integer(1),
            ),
        )
    )

    yukawa, entries = build_charged_lepton_yukawa(
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

    inputs = MasterRGEInputs(
        fermion_dimension=fermion_basis.dimension,
        yukawa=yukawa,
        quartic=zero_quartic,
        gauge_sectors=(),
    )

    result = calculate_complete_master_rge(
        model=scalar_model,
        inputs=inputs,
        output_component=component,
        coefficient=C,
    )

    print("=" * 72)
    print("SM CHARGED-LEPTON YUKAWA WEINBERG BENCHMARK")
    print("=" * 72)

    print()
    print("Nonzero y_ija entries:")

    for key in sorted(entries):
        print(f"  y{key} = {entries[key]}")

    print()
    print("Eq. (4.85) contributions to (dot kappa)/kappa:")

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

    expected = -ye * sp.conjugate(ye)
    difference = sp.factor(
        sp.expand(beta_ratio - expected)
    )

    print()
    print("Expected beta_kappa/kappa :", expected)
    print("Engine beta_kappa/kappa   :", beta_ratio)
    print("Difference                :", difference)
    print()

    if sp.simplify(difference) != 0:
        print(
            "FAIL: the Yukawa sector still has a normalization, "
            "orientation, or permutation mismatch."
        )
        return 1

    print("PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
