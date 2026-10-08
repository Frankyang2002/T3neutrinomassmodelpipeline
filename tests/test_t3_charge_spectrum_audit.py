"""Exact charge-spectrum regression tests."""
from fractions import Fraction
from Numerical.diagnostics.AuditT3ChargeSpectrum import audit, charge_spectrum


def test_doublet_and_triplet_charges():
    assert [x['Q'] for x in charge_spectrum(2, Fraction(1, 2))] == ['1', '0']
    assert [x['Q'] for x in charge_spectrum(3, Fraction(0))] == ['1', '0', '-1']
    assert [x['Q'] for x in charge_spectrum(1, Fraction(-1))] == ['-1']


def test_audit_preserves_unknown_multiplicity():
    item = {'model': {'dS1': 2, 'dS2': 2, 'dF': 1, 'alpha': -1},
            'vertices': {'T3Y1CG': {'nonzero_component_count': 2}}}
    result = audit({'cases': [item]})
    assert result['models'][0]['fields']['F']['weights_and_charges'][0]['Q'] == '0'
    assert result['models'][0]['verified_loop_multiplicity'] is None
    assert result['models'][0]['compensation_factor'] is None
