from __future__ import annotations

"""
Pure-SM validation of the generic psi^2 phi^2 Eq. (4.85) implementation.

The purpose of this file is diagnostic.  It removes every BSM field and tests
the Weinberg operator in the Standard Model alone.

Benchmarks
----------
Gauge-only:
    16 pi^2 beta_kappa |_gauge = -3 g2^2 kappa

in the usual SM Weinberg-operator convention.  In particular there is no
g1^2 term.

We evaluate the neutral real component

    C_(nu nu H0_R H0_R) = kappa/2,

so beta_kappa / kappa is identical to beta_C / C.

The test prints:
  * every gauge contribution using the currently implemented anomalous
    dimensions;
  * the result after inserting the physically expected g_alpha^2 factor in
    gamma_{c,f}, without modifying the production engine.

If the second result still disagrees with -3 g2^2, then Eq. (4.85)'s gauge
terms or generator conventions require further correction beyond A.3.
"""

import sympy as sp

from tests.python.GeneralWeinbergRGEGenerator import (
    ComplexScalar,
    MasterRGEInputs,
    RGEModel,
    calculate_complete_master_rge,
    calculate_master_rge,
    quadratic_casimir_matrix,
    scalar_collinear_anomalous_dimension,
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


class SparseWilsonLookup:
    """Sparse C_ijab lookup supporting the [] interface used by Eq. (4.85)."""

    def __init__(self, components):
        self.components = components

    def __getitem__(self, key):
        return sp.sympify(self.components.get(tuple(key), sp.S.Zero))


def zero_yukawa(i: int, j: int, a: int) -> sp.Expr:
    return sp.S.Zero


def zero_quartic(a: int, b: int, c: int, d: int) -> sp.Expr:
    return sp.S.Zero


def pure_sm_bases() -> tuple[RGEModel, FermionBasis]:
    """Only H and L are required for the gauge-only Weinberg benchmark."""

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

    return scalar_model, fermion_basis


def neutral_component(
    scalar_model: RGEModel,
    fermion_basis: FermionBasis,
) -> tuple[int, int, int, int]:
    nu = fermion_basis.global_index("L", 1)
    h0_r = scalar_model.block("H").local_to_global(3)
    return nu, nu, h0_r, h0_r


def ratio_to_coefficient(
    value: sp.Expr,
    coefficient_value: sp.Expr,
) -> sp.Expr:
    return sp.factor(sp.simplify(value / coefficient_value))


def corrected_gamma_fermion(
    inputs: MasterRGEInputs,
    i: int,
    j: int,
) -> sp.Expr:
    r"""Diagnostic A.3 with the expected gauge-coupling factors inserted.

    gamma_cf = -3 sum_alpha g_alpha^2 C2(F_alpha).

    Yukawas are zero in this benchmark.
    """

    result = sp.S.Zero

    for sector in inputs.gauge_sectors:
        casimir = quadratic_casimir_matrix(
            sector.fermion_generators
        )
        result += (
            -3
            * sector.coupling**2
            * casimir[i - 1, j - 1]
        )

    return sp.simplify(result)


def print_result(
    title: str,
    result: dict[str, sp.Expr],
    c_value: sp.Expr,
) -> sp.Expr:
    print()
    print("=" * 72)
    print(title)
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
        value = ratio_to_coefficient(result[name], c_value)
        print(f"{name:34s}: {value}")

    return ratio_to_coefficient(result["total"], c_value)


def main() -> int:
    scalar_model, fermion_basis = pure_sm_bases()

    wilson = build_weinberg_wilson_tensor(
        scalar_model,
        fermion_basis,
        kappa,
    )
    validate_weinberg_tensor_symmetry(wilson)

    C = SparseWilsonLookup(wilson)
    component = neutral_component(scalar_model, fermion_basis)
    c_value = C[component]

    assert sp.simplify(c_value - kappa / 2) == 0

    gauge_sectors = build_gauge_sectors(
        scalar_model,
        fermion_basis,
    )

    base_inputs = MasterRGEInputs(
        fermion_dimension=fermion_basis.dimension,
        yukawa=zero_yukawa,
        quartic=zero_quartic,
        gauge_sectors=gauge_sectors,
    )

    # ----------------------------------------------------------------------
    # 1. Current production implementation
    # ----------------------------------------------------------------------

    current = calculate_complete_master_rge(
        model=scalar_model,
        inputs=base_inputs,
        output_component=component,
        coefficient=C,
    )

    current_ratio = print_result(
        "CURRENT ENGINE: gauge-only SM Weinberg RGE",
        current,
        c_value,
    )

    # ----------------------------------------------------------------------
    # 2. Diagnostic A.3 correction only
    # ----------------------------------------------------------------------

    def gamma_scalar(a: int, b: int) -> sp.Expr:
        return scalar_collinear_anomalous_dimension(
            model=scalar_model,
            inputs=base_inputs,
            a=a,
            b=b,
        )

    def gamma_fermion(i: int, j: int) -> sp.Expr:
        return corrected_gamma_fermion(
            base_inputs,
            i,
            j,
        )

    corrected_inputs = MasterRGEInputs(
        fermion_dimension=fermion_basis.dimension,
        yukawa=zero_yukawa,
        quartic=zero_quartic,
        gauge_sectors=gauge_sectors,
        gamma_scalar=gamma_scalar,
        gamma_fermion=gamma_fermion,
    )

    corrected = calculate_master_rge(
        model=scalar_model,
        inputs=corrected_inputs,
        output_component=component,
        coefficient=C,
    )

    corrected_ratio = print_result(
        "DIAGNOSTIC: A.3 WITH g_alpha^2 INSERTED",
        corrected,
        c_value,
    )

    expected = -3 * g2**2

    print()
    print("=" * 72)
    print("SM GAUGE BENCHMARK")
    print("=" * 72)
    print("Expected beta_kappa/kappa :", expected)
    print("Current engine            :", current_ratio)
    print("A.3-only correction       :", corrected_ratio)

    current_difference = sp.factor(
        sp.expand(current_ratio - expected)
    )
    corrected_difference = sp.factor(
        sp.expand(corrected_ratio - expected)
    )

    print()
    print("Current - expected        :", current_difference)
    print("Corrected - expected      :", corrected_difference)

    print()

    if sp.simplify(corrected_difference) == 0:
        print("PASS: A.3 was the only gauge-sector problem.")
        return 0

    print(
        "EXPECTED FAIL: inserting g_alpha^2 into A.3 is not sufficient; "
        "the remaining Eq. (4.85) gauge implementation must be corrected."
    )
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
