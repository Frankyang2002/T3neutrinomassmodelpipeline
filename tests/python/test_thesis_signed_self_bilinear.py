"""Signed-real self-bilinear convention tests (no Wolfram required).

The Matchete source fixes S=(S_minus+i*S_plus)/sqrt(2). Changing the sign
of the real holomorphic bilinear exchanges the two *masses*, not the fields.
"""
import math

import pytest

from Numerical.matching.GeneralT3FlavourC5 import (
    signed_real_mass_basis_parameters, validate_real_self_bilinear,
)

SELF = ('A_-2', 'A_+0', 'D_-2', 'D_+0', 'E_-2', 'E_+0')


@pytest.mark.parametrize('model', SELF)
@pytest.mark.parametrize('value', [1000.0, -1000.0])
def test_sign_preserves_original_matchete_field_labels(model, value):
    species = '2' if model.endswith('-2') else '1'
    msq = 1_250_000.0
    plus = f'T3{species}PlusMass'
    minus = f'T3{species}MinusMass'
    raw = {plus: math.sqrt(msq + 2 * abs(value)),
           minus: math.sqrt(msq - 2 * abs(value)),
           'lambdaT3': 0.1}
    result = signed_real_mass_basis_parameters(model, raw, value)
    assert result[plus]**2 == pytest.approx(msq + 2 * value)
    assert result[minus]**2 == pytest.approx(msq - 2 * value)
    assert result['lambdaT3'] == raw['lambdaT3']
    assert raw[plus] > raw[minus]  # Input not mutated.


def test_sign_crossing_swaps_order_but_preserves_spectrum():
    raw = {'T31PlusMass': math.sqrt(1000**2+100),
           'T31MinusMass': math.sqrt(1000**2-100)}
    plus = signed_real_mass_basis_parameters('D_+0', raw, 50)
    minus = signed_real_mass_basis_parameters('D_+0', raw, -50)
    assert plus['T31PlusMass'] == minus['T31MinusMass']
    assert plus['T31MinusMass'] == minus['T31PlusMass']
    assert sorted(plus.values()) == sorted(minus.values())


def test_negative_real_bilinear_is_not_mistaken_for_unsupported_complex_phase():
    assert validate_real_self_bilinear('D_+0', -800.0) == -800.0 + 0j
    assert validate_real_self_bilinear('A_-2', 800.0) == 800.0 + 0j


@pytest.mark.parametrize('bad', [0, 1j, 123 + 0.5j, complex(float('nan'), 0)])
def test_zero_complex_or_nonfinite_bilinear_is_rejected(bad):
    with pytest.raises(ValueError):
        validate_real_self_bilinear('D_+0', bad)


def test_small_numerical_imaginary_roundoff_is_discarded():
    assert validate_real_self_bilinear('D_+0', -500 + 1e-9j) == -500.0 + 0j


def test_complex_phase_symbol_rejects_unverified_reinterpretation():
    params = {'T31PlusMass': 1100.0, 'T31MinusMass': 1090.0}
    with pytest.raises(ValueError, match='T3NeutralPhase'):
        signed_real_mass_basis_parameters(
            'D_+0', params, -50.0, coefficient_symbols={'T3NeutralPhase'})


def test_no_scalar_field_redefinition_applied_to_unrelated_parameters():
    params = {'T31PlusMass': 1100., 'T31MinusMass': 1090.,
              'y1Conj': {'real': [[1]], 'imag': [[0]]}}
    got = signed_real_mass_basis_parameters('D_+0', params, -50)
    assert got['y1Conj'] == params['y1Conj']
    assert got['T31PlusMass'] == 1090
    assert got['T31MinusMass'] == 1100
