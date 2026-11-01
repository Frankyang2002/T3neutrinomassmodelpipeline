"""Offline exact checks and deliberate failure cases; no native Wolfram runtime."""
from __future__ import annotations

from fractions import Fraction as F
import pytest

from scripts.T3GeneralGroupFactorProjection import (
    C2, T2, ProjectionError, field_data,
    expected_gauge_betas, expected_scalar_gauge_factors,
    expected_yukawa_gauge_factors, project_native_beta,
)


def test_exact_su2_and_field_hypercharges():
    assert [C2(i) for i in (1, 2, 3)] == [0, F(3, 4), 2]
    assert [T2(i) for i in (1, 2, 3)] == [0, F(1, 2), 2]
    assert field_data((1, 3, 2), -2)["S2"] == (3, F(0))


def test_gauge_coefficients_across_dirac_and_majorana():
    assert expected_gauge_betas((1, 3, 2), -2, "UV") == {
        "gY": F(55, 6), "g2": F(-1, 2), "g3": F(-7)}
    assert expected_gauge_betas((1, 3, 2), -2, "EFT1") == {
        "gY": F(43, 6), "g2": F(-5, 2), "g3": F(-7)}
    assert expected_gauge_betas((2, 2, 1), -1, "UV")["g2"] == F(-17, 6)
    assert expected_gauge_betas((2, 2, 3), -1, "UV")["g2"] == F(7, 6)


def test_scalar_quartic_and_bilinear_wavefunction_gauge_factors():
    assert expected_scalar_gauge_factors((1, 1, 1, 1, 0, 0), (1, 3, 2), -2) == {
        "g2": F(-9, 2), "gY": F(-15, 2)}
    # Two hypercharge-zero triplets: gauge SU2 wavefunction -3*(2+2)=-12.
    assert expected_scalar_gauge_factors((0, 0, 0, 0, 0, 2), (1, 3, 2), -2) == {
        "g2": F(-12), "gY": F(0)}
    assert expected_yukawa_gauge_factors((1, 3, 2), -2) == {
        "g2": F(-9, 2), "gY": F(-3, 2)}


def test_native_mathematica_exact_polynomial_extraction():
    assert project_native_beta("Rational[41,6]*gY^3", gauge="gY") == F(41, 6)
    assert project_native_beta("-7*g3^3", gauge="g3") == F(-7)
    assert project_native_beta("2*g2^2*lambdaH1-3*g2^2*lambdaH1+9*g2^2*lambdaH2",
                               gauge="g2", coupling="lambdaH1") == -1
    assert project_native_beta("-Rational[15,2]*gY^2*lambdaH1+lambdaH1^2+3*Tr[y1]*lambdaH1",
                               gauge="gY", coupling="lambdaH1") == F(-15, 2)


def test_yukawa_matrix_projector_preserves_only_target():
    raw = "-9*g2^2*Matrix[y1][$i,$j]/2 + 4*gY^2*Matrix[y1][$i,$j] + 3*g2^2*Matrix[y2][$i,$j]"
    assert project_native_beta(raw, gauge="g2", matrix_coupling="y1") == F(-9, 2)
    assert project_native_beta(raw, gauge="gY", matrix_coupling="y1") == 4


def test_negative_controls_fail_closed():
    for raw in ("", "lambdaH1^2", "Tr[y1]*g2^2*lambdaH1"):
        with pytest.raises(ProjectionError):
            project_native_beta(raw, gauge="g2", coupling="lambdaH1")
    with pytest.raises(ProjectionError):
        expected_scalar_gauge_factors((0, 0, 0, 0, 0, 0), (1, 3, 2), -2)
    with pytest.raises(ProjectionError):
        project_native_beta("g2^2*lambdaH1^2", gauge="g2", coupling="lambdaH1")
