"""Canonical export integration: reject missing model metadata and unknown couplings."""
import pytest
from Numerical.core.GeneralCouplingLayout import spec_for


def test_extra_fermion_yukawa_types():
    assert spec_for('y1Conj').kind == 'matrix'
    assert spec_for('y2Conj').kind == 'matrix'


def test_general_quartic_complex_type():
    assert spec_for('lambdaT3GeneralQ000004I1').kind == 'complex'


def test_mass_bilinear_complex_type():
    assert spec_for('muT3GeneralS2Sq').kind == 'complex'


def test_unknown_operator_rejected():
    with pytest.raises(ValueError, match='Unrecognised'):
        spec_for('unlistedNewOperator')
