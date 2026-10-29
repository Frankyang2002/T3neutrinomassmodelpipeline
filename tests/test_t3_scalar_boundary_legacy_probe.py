import numpy as np
import pytest
from scripts.ProbeT3ScalarBoundaryAndLegacyC5 import mass_rescale, comparison_metrics


def test_mass_rescale_preserves_other_parameters():
    c={'base_state':{'ordinary':{'mS1Sq':1e6,'mS2Sq':4e6,'lambdaT3':{'real':.5,'imag':0}}}}
    scaled=mass_rescale(c,0.5)
    assert scaled['base_state']['ordinary']['mS1Sq']==250000
    assert scaled['base_state']['ordinary']['mS2Sq']==1000000
    assert c['base_state']['ordinary']['mS1Sq']==1e6
    assert scaled['base_state']['ordinary']['lambdaT3']==c['base_state']['ordinary']['lambdaT3']


def test_invalid_factor():
    with pytest.raises(ValueError):mass_rescale({'base_state':{'ordinary':{}}},0)


def test_comparison_identical():
    a=np.eye(3,dtype=complex)*1e-12
    assert comparison_metrics(a,a)['relative_max_error']==0


def test_comparison_detects_global_sign():
    a=np.eye(3,dtype=complex)*1e-12
    result=comparison_metrics(a,-a)
    assert result['hypothesis_relative_errors']['minus_fast']==0
    assert result['relative_max_error']==2


def test_comparison_detects_constant_scale():
    a=np.eye(3,dtype=complex)*1e-12
    result=comparison_metrics(a,2*a)
    assert result['hypothesis_relative_errors']['best_complex_scale_fast']<1e-14
