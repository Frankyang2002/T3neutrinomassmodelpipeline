"""Fast unit checks; no Wolfram installation or fitted campaign required."""
from __future__ import annotations

from copy import deepcopy
import numpy as np
import pytest

from Numerical.running.GeneralT3BCQuarticChannelLogAudit import (
    CHANNELS, masked_endpoint, validate_source_inventory,
    eft_quartic_channels, relative_error,
)


def endpoint():
    return {
        'mu_gev': 1000.0,
        'lambdaT3': {'real':0.11,'imag':0.0},
        'reused_quartics': {
            'lambdaHHdagS2S2': {'real':0.012,'imag':0.0},
            'lambdaHHdagS1barS1bar': {'real':-0.021,'imag':0.0},
            'unrelated': {'real':3.0,'imag':0.0}},
        'quadratic_mu2_gev2': {'real':0,'imag':800},
    }


@pytest.mark.parametrize('active', [x[0] for x in CHANNELS])
def test_masked_endpoint_keeps_exactly_one_quartic(active):
    ep=endpoint();before=deepcopy(ep)
    got=masked_endpoint(ep,active)
    assert ep==before and got is not ep
    assert got['quadratic_mu2_gev2']==ep['quadratic_mu2_gev2']
    assert got['reused_quartics']['unrelated']==ep['reused_quartics']['unrelated']
    for k in ('lambdaT3','lambdaHHdagS2S2','lambdaHHdagS1barS1bar'):
        node=got['lambdaT3'] if k=='lambdaT3' else got['reused_quartics'][k]
        old=ep['lambdaT3'] if k=='lambdaT3' else ep['reused_quartics'][k]
        assert node==(old if active==k else {'real':0.0,'imag':0.0})


def test_invalid_active_quartic_fails():
    with pytest.raises(ValueError,match='unknown one-scale'):
        masked_endpoint(endpoint(),'lambdaH1')


def test_missing_direct_quartic_fails():
    ep=endpoint();del ep['reused_quartics']['lambdaHHdagS2S2']
    with pytest.raises(ValueError,match='missing a B/C direct'):
        masked_endpoint(ep,'lambdaT3')


def test_source_quartic_inventory_is_exact():
    src=' '.join('Bar[Coupling['+s+', {}, 0]]' for _,s,_ in CHANNELS)
    found=validate_source_inventory(src)
    assert len(found)==3 and set(found)=={x[1] for x in CHANNELS}
    with pytest.raises(ValueError,match='precisely three'):
        validate_source_inventory(src+' Coupling[lambdaBogus, {}, 0]')
    with pytest.raises(ValueError,match='precisely three'):
        validate_source_inventory(src.replace('Coupling[lambdaT3, {}, 0]', ''))


def encode(a):
    arr=np.asarray(a,complex)
    return {'real':arr.real.tolist(),'imag':arr.imag.tolist()}


def signed_case(scale=1.0):
    a=np.diag([1.,2.,3.])*scale*1e-14
    b=np.diag([3.,2.,1.])*scale*1e-14
    c=np.array([[0,1.,0],[1.,0,0],[0,0,0]])*scale*1e-14
    pairs=(('native:lambdaT3',a),('general:lambdaHHdagS2S2',b),
           ('general:lambdaHHdagS1barS1bar',c))
    return {
        'case':'1111',
        'quartic_provenance': {
            'bc_direct_quartic_channel_signs': {
                'general:lambdaHHdagS2S2': -1,
                'general:lambdaHHdagS1barS1bar': -1},
            'bc_direct_quartic_map_physically_certified':False,
        },
        'per_general_interaction':[
            {'name':name,'category':'quartic_decomposition',
             'c5_leading_log':encode(x)} for name,x in pairs
        ]+[{'name':'F_tree_source:11','category':'fermion_tree_decomposition',
             'c5_leading_log':encode(a)}],
        'c5_induced_wilson_log':encode(a+b+c),
    }


def test_eft_quartic_channels_reproduce_saved_total():
    channels,total,info=eft_quartic_channels(signed_case())
    assert set(channels)=={r[2] for r in CHANNELS}
    np.testing.assert_allclose(sum(channels.values()),total)
    assert info['all_eft1_quartics_sum_relative_error']<1e-13
    assert info['three_selected_eft1_quartics_sum_relative_error']<1e-13


def test_unsigned_eft_case_fails():
    case=signed_case()
    case['quartic_provenance']['bc_direct_quartic_channel_signs'].pop('general:lambdaHHdagS2S2')
    with pytest.raises(ValueError,match='TWO required'):
        eft_quartic_channels(case)


def test_extra_nonzero_eft_quartic_fails():
    case=signed_case()
    other=np.diag([1.,0.,0.])*1e-14
    case['per_general_interaction'].append({'name':'portal:test',
        'category':'quartic_decomposition','c5_leading_log':encode(other)})
    case['c5_induced_wilson_log']=encode(
        np.asarray(case['c5_induced_wilson_log']['real'])+other)
    with pytest.raises(ValueError,match='additional nonzero'):
        eft_quartic_channels(case)


def test_relative_error_is_norm_invariant():
    m=np.eye(3)*1e-16
    assert relative_error(m,m)==0
    assert relative_error(2*m,m)==pytest.approx(1)


def test_refuses_incompatible_matrices():
    with pytest.raises(ValueError,match='wrong-shaped'):
        relative_error(np.eye(2),np.eye(3))


def test_three_channel_full_comparison_with_synthetic_evaluator(monkeypatch):
    """Exercise masks, exact three-channel closure and diagnostic factor of 2."""
    import sys
    import types
    from Numerical.running.GeneralT3BCQuarticChannelLogAudit import run_channel_comparison
    from Numerical.running.GeneralT3MatcheteLogCoefficientAudit import decode_matrix

    saved=signed_case()
    ep=endpoint()
    qnames=[x[0] for x in CHANNELS]
    ref_vals=np.array([ep['lambdaT3']['real'],
                       ep['reused_quartics'][qnames[1]]['real'],
                       ep['reused_quartics'][qnames[2]]['real']])
    mats=[np.eye(3)*0.25e-14,np.eye(3)*0.5e-14,
          np.array([[0,1,0],[1,0,0],[0,0,0]])*0.3e-14]
    eft_values={x[0]:decode_matrix(next(v['c5_leading_log']
                   for v in saved['per_general_interaction']
                   if v['name']==x[2]))/(2*np.log(100)) for x in CHANNELS}
    base=sum((r*m for r,m in zip(ref_vals,mats)),np.zeros((3,3),complex))
    saved['c5_matching_anchored_reference']=encode(base)
    prior=2*sum(eft_values.values(),np.zeros((3,3),complex))

    def fake(model, endpoint, heavy, compiled):
        q=np.array([endpoint['lambdaT3']['real'],
                endpoint['reused_quartics'][qnames[1]]['real'],
                endpoint['reused_quartics'][qnames[2]]['real']])
        result=sum((r*m for r,m in zip(q,mats)),np.zeros((3,3),complex))
        if compiled=='shifted':
            hierarchy=heavy[0][0]/1e5
            result=result+np.log(2)*sum((
                2*q[i]/ref_vals[i]*eft_values[qnames[i]]/hierarchy
                for i in range(3)),np.zeros((3,3),complex))
        return result

    stub=types.ModuleType('scripts.RunT3GeneralScalarBCExploratoryFit')
    stub.c5_at_endpoint=fake
    monkeypatch.setitem(sys.modules,'scripts.RunT3GeneralScalarBCExploratoryFit',stub)
    import scripts
    monkeypatch.setattr(scripts,'RunT3GeneralScalarBCExploratoryFit',stub,raising=False)
    source=' + '.join('Coupling['+x[1]+', {}, 0]' for x in CHANNELS)
    heavy=(np.array([1e5]*3), np.ones((3,3)),np.ones((3,3)),0.)
    result=run_channel_comparison('B_-1',case=saved,endpoint=ep,heavy=heavy,
            coefficient_text=source,original_compiled='original',
            shifted_compiled='shifted',prior_matchete=prior,
            hierarchy_factors=(1,2))
    assert result['status']=='Computed'
    assert len(result['hierarchy_scan'])==2
    for scan in result['hierarchy_scan']:
        assert scan['checks']['matchete_original_three_channel_closure']<1e-12
        assert scan['checks']['matchete_log_three_channel_closure']<1e-12
        for channel in scan['channels']:
            c=channel['comparison']
            assert c['best_fit_complex_factor']['real']==pytest.approx(2,rel=1e-12)
            assert c['relative_matrix_residual_after_best_fit_factor']<1e-12
