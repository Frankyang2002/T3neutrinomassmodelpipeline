"""Pure contract tests: no expensive RG/matching or Wolfram required."""
from copy import deepcopy
import numpy as np
import pytest

from scripts.RunT3GeneralZeeSensitivity import (
    deform, leading_tensor, verify_nominal_invariance, gamma_grid,
    params_from, compensated_quartic_fraction, baseline_consistency,
)


def sample():
    return {'base_state': {'ordinary': {
        'lambdaT3': {'real': 0.1, 'imag': 0.},
        'y1': {'real': [[0.01,0.02,0.03], [0.05,0.04,0.06], [0.08,0.07,0.09]],
               'imag': np.zeros((3,3)).tolist()},
        'y2': {'real': [[0.03,0.08,0.05], [0.02,0.04,0.01], [0.07,0.06,0.09]],
               'imag': np.zeros((3,3)).tolist()},
        'lambdaT3GeneralQ020101I1': {'real': 0.01, 'imag': 0.},
        'muT3GeneralS2Sq': {'real': 100, 'imag': 0.},
    }}}


@pytest.mark.parametrize('direction',('T1','T2','T3'))
@pytest.mark.parametrize('gamma',(0.99,1.,1.01,1.2))
def test_nominal_tensor_and_independent_scalar_invariant(direction,gamma):
    base=sample(); copy=deepcopy(base)
    changed=deform(base,direction,gamma)
    verify_nominal_invariance(base,changed)
    np.testing.assert_allclose(leading_tensor(base),leading_tensor(changed),rtol=1e-13,atol=1e-18)
    assert base==copy  # no mutation
    old=base['base_state']['ordinary']
    new=changed['base_state']['ordinary']
    for name in ('lambdaT3GeneralQ020101I1','muT3GeneralS2Sq'):
        assert old[name] == new[name]


def test_rescaling_directions_distinct():
    b=sample()
    a=deform(b,'T1',4.)['base_state']['ordinary']
    assert a['lambdaT3']['real']==pytest.approx(.4)
    assert a['y1']['real'][0][0]==pytest.approx(.005)
    assert a['y2']['real'][0][0]==pytest.approx(.015)
    a=deform(b,'T2',4.)['base_state']['ordinary']
    assert a['y1']['real'][0][0]==pytest.approx(.01)
    assert a['y2']['real'][0][0]==pytest.approx(.0075)
    a=deform(b,'T3',4.)['base_state']['ordinary']
    assert a['lambdaT3']['real']==pytest.approx(.1)
    assert a['y1']['real'][0][0]==pytest.approx(.04)
    assert a['y2']['real'][0][0]==pytest.approx(.0075)


def test_gamma_grid_and_validation():
    vals=gamma_grid(.99,1.01,11)
    assert len(vals)==11 and vals[5]==1.
    assert all(x<y for x,y in zip(vals,vals[1:]))
    with pytest.raises(ValueError): gamma_grid(.9,1.1,10)
    with pytest.raises(ValueError): gamma_grid(1.,1.1,11)


def test_yukawa_params():
    p=params_from(sample())
    assert len(p)==18 and p['y1_01']==pytest.approx(.02)
    assert p['y2_20']==pytest.approx(.07)


def test_bc_quartic_compensating_rule():
    for gam in (.8,1.,1.2):
        q=compensated_quartic_fraction(.1,.5,.5*gam)
        assert q * abs(.5*gam)==pytest.approx(.05)
    with pytest.raises(ValueError): compensated_quartic_fraction(.5,.5,.1)


def test_baseline_check_rejects_wrong_observable():
    from scripts.RunT3GeneralZeeSensitivity import OBS
    baseline={'prediction':{key:1. for key in OBS}}
    rows=[{'status':'Success',**{key:1. for key in OBS}} for _ in range(3)]
    baseline_consistency(rows,baseline)
    rows[0][OBS[0]]=2.
    with pytest.raises(RuntimeError): baseline_consistency(rows,baseline)
