import numpy as np
import pytest
from scripts.SearchT3JointScalarYukawa import configure, quartic_names, sample_parameters


def history():
    return {'cfg':{'base_state':{'ordinary':{'mS1Sq':100.0,'mS2Sq':400.0,
       'lambdaT3':{'real':5e-9,'imag':0.},'lambdaHS1':{'real':0.1,'imag':0.},
       'y1':{},'y2':{}}}}, 'seed':np.zeros(18),'lambda_fixed':5e-9}


def test_quartic_names_excludes_topology():
    assert quartic_names(history()['cfg'],4)==('lambdaHS1',)


def test_log_uniform_mass_coordinates():
    r=np.full(21,0.5)
    y,f,q=sample_parameters(r,np.zeros(18),np.ones(18),('lambdaHS1',),1.,.1,10.)
    assert np.allclose(f,[1,1])
    assert np.allclose(y,.5)
    assert np.allclose(q,[0])


def test_mass_bounds():
    r=np.zeros(20)
    _,f,_=sample_parameters(r,np.zeros(18),np.ones(18),(),1.,.1,10.)
    assert np.allclose(f,.1)


def test_invalid_factor_raises_before_yukawa_mutation():
    with pytest.raises(ValueError): configure(history(),np.zeros(18),[0,1],(),())
