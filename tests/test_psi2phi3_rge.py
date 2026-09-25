"""Unit regressions for the generic dimension-six psi^2 phi^3 RGE kernel."""

from __future__ import annotations

import sympy as sp

from RGE.general.GaugeGenerators import GaugeSector
from RGE.general.Psi2Phi3RGE import calculate_psi2phi3_rge
from RGE.general.ScalarBasis import ComplexScalar, ScalarBasis
from RGE.general.WilsonTensorRGE import WilsonRGEInputs


class _SparseC6:
    def __init__(self, value: sp.Expr):
        self.value = value

    def __getitem__(self, key: tuple[int, int, int, int, int]) -> sp.Expr:
        i, j, a, b, c = key
        if (i, j, a, b, c) == (1, 1, 1, 1, 1):
            return self.value
        return sp.S.Zero


def _one_component_inputs():
    y, lam, gamma_s, gamma_f = sp.symbols(
        "y lam gamma_s gamma_f",
        real=True,
    )
    model = ScalarBasis(
        (
            ComplexScalar(
                name="X",
                su2_dimension=1,
                hypercharge=sp.Integer(0),
            ),
        )
    )

    def yukawa(i: int, j: int, a: int) -> sp.Expr:
        return y if (i, j, a) == (1, 1, 1) else sp.S.Zero

    def quartic(a: int, b: int, c: int, d: int) -> sp.Expr:
        return lam if (a, b, c, d) == (1, 1, 1, 1) else sp.S.Zero

    def scalar_gamma(a: int, b: int) -> sp.Expr:
        return gamma_s if (a, b) == (1, 1) else sp.S.Zero

    def fermion_gamma(i: int, j: int) -> sp.Expr:
        return gamma_f if (i, j) == (1, 1) else sp.S.Zero

    inputs = WilsonRGEInputs(
        fermion_dimension=1,
        yukawa=yukawa,
        quartic=quartic,
        gamma_scalar=scalar_gamma,
        gamma_fermion=fermion_gamma,
    )
    return model, inputs, y, lam, gamma_s, gamma_f


def test_one_component_permutation_multiplicities_match_eq_425() -> None:
    model, inputs, y, lam, gamma_s, gamma_f = _one_component_inputs()
    c = sp.Symbol("C")

    result = calculate_psi2phi3_rge(
        model,
        inputs,
        (1, 1, 1, 1, 1),
        coefficient=_SparseC6(c),
    )

    assert sp.simplify(result["yukawa_wavefunction"] - y**2 * c) == 0
    assert sp.simplify(result["scalar_pair"] - 3 * lam * c) == 0
    assert sp.simplify(result["mixed_yukawa_gauge"] - 18 * y**2 * c) == 0
    assert sp.simplify(
        result["conjugate_coefficient_yukawa"]
        - 3 * y**2 * sp.conjugate(c)
    ) == 0
    assert sp.simplify(
        result["scalar_anomalous_dimension"] - 3 * gamma_s * c
    ) == 0
    assert sp.simplify(
        result["fermion_anomalous_dimension"] - 2 * gamma_f * c
    ) == 0


def test_zero_interactions_give_zero_self_running() -> None:
    model = ScalarBasis(
        (
            ComplexScalar(
                name="X",
                su2_dimension=1,
                hypercharge=sp.Integer(0),
            ),
        )
    )

    inputs = WilsonRGEInputs(
        fermion_dimension=1,
        yukawa=lambda i, j, a: sp.S.Zero,
        quartic=lambda a, b, c, d: sp.S.Zero,
    )

    result = calculate_psi2phi3_rge(
        model,
        inputs,
        (1, 1, 1, 1, 1),
        coefficient=_SparseC6(sp.Symbol("C")),
    )

    assert result["total"] == 0


def test_one_component_gauge_permutation_factors() -> None:
    g, theta_value, t_value = sp.symbols("g theta t", real=True)
    c = sp.Symbol("C")
    model = ScalarBasis(
        (
            ComplexScalar(
                name="X",
                su2_dimension=1,
                hypercharge=sp.Integer(0),
            ),
        )
    )
    theta = sp.zeros(2)
    theta[0, 0] = theta_value
    fermion_generator = sp.Matrix([[t_value]])
    sector = GaugeSector(
        coupling=g,
        scalar_generators=(theta,),
        fermion_generators=(fermion_generator,),
    )
    inputs = WilsonRGEInputs(
        fermion_dimension=1,
        yukawa=lambda i, j, a: sp.S.Zero,
        quartic=lambda a, b, c, d: sp.S.Zero,
        gauge_sectors=(sector,),
    )

    result = calculate_psi2phi3_rge(
        model,
        inputs,
        (1, 1, 1, 1, 1),
        coefficient=_SparseC6(c),
    )

    assert sp.simplify(
        result["scalar_pair"]
        + 6 * g**2 * theta_value**2 * c
    ) == 0
    assert sp.simplify(
        result["mixed_yukawa_gauge"]
        - 24 * g**2 * t_value * theta_value * c
    ) == 0
