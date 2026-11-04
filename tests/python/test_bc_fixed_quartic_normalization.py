"""Fixed B/C topology/direct quartic normalizations; no fitted overall scale."""
from __future__ import annotations

import ast
import math
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from Numerical.running.GeneralT3SequentialLeadingLog import TensorBundle
from Numerical.running.GeneralT3ConditionalBCQuarticMap import (
    apply_conditional_bc_direct_quartic_map,
)
from Numerical.running.GeneralT3ConditionalBCQuarticNormalization import (
    CHANNEL_FACTORS, CONVENTION, apply_conditional_bc_channel_normalization,
)
from scripts.CheckT3BCFixedNormalization import (
    compare_fixed, fixed_spurion_from_signed_case,
)

NAMES = ('native:lambdaT3', 'general:lambdaHHdagS2S2',
         'general:lambdaHHdagS1barS1bar', 'portal:lambdaH1Inv1')


def make_bundle(model='B_-1'):
    n = len(NAMES)
    arr = np.zeros((n, 4, 4, 8, 8))
    for c in range(n):
        for i in range(4):
            arr[c,i,i,:,:] = (0.01+c*0.007)*np.eye(8)
    descriptions = [
        {'name':name.split(':')[1],
         'origin':'native Matchete source CG' if name.startswith('native:')
         else 'sector_order_candidate_RGBeta_to_SU2_reference_basis'}
        for name in NAMES]
    return TensorBundle(model,np.diag(np.arange(1,9,dtype=float)),arr.sum(axis=0),
                        arr,NAMES,{'quartic_channels':descriptions,
                                   'matchete_rgbeta_map_certified':False})


@pytest.mark.parametrize('model',['B_-1','C_-1'])
def test_fixed_map_after_sign_map_and_no_mutation(model):
    before = make_bundle(model)
    signed = apply_conditional_bc_direct_quartic_map(model,before)
    orig = before.channels.copy()
    signed_copy = signed.channels.copy()
    result = apply_conditional_bc_channel_normalization(model,signed)
    np.testing.assert_array_equal(before.channels,orig)
    np.testing.assert_array_equal(signed.channels,signed_copy)
    for idx,name in enumerate(NAMES):
        sign = -1 if name.startswith('general:lambdaHHdagS') else 1
        factor = CHANNEL_FACTORS.get(name,1)
        np.testing.assert_allclose(result.channels[idx], sign*factor*orig[idx],
                                   rtol=1e-15,atol=1e-16)
    np.testing.assert_allclose(result.quartic,result.channels.sum(axis=0),rtol=1e-14)
    np.testing.assert_array_equal(result.mass2,before.mass2)
    assert result.provenance['bc_channel_normalization_factors']==CHANNEL_FACTORS
    assert result.provenance['bc_channel_normalization_convention']==CONVENTION
    assert result.provenance['bc_channel_normalization_affects_rges'] is False
    assert result.provenance['bc_channel_normalization_factors_fitted'] is False
    assert result.provenance['bc_channel_normalization_physical_certification'] is False
    assert not result.provenance['matchete_rgbeta_map_certified']
    assert len([x for x in result.provenance['quartic_channels']
                if 'conditional_log_normalization_factor' in x])==3


@pytest.mark.parametrize('model',['A_-4','C_+1','E_-2'])
def test_nonmixed_rejected(model):
    with pytest.raises(ValueError,match='restricted to B_-1/C_-1'):
        apply_conditional_bc_channel_normalization(model,make_bundle(model))


def test_unsigned_rejected():
    with pytest.raises(ValueError,match='sign map must run first'):
        apply_conditional_bc_channel_normalization('B_-1',make_bundle())


def test_double_application_rejected():
    signed = apply_conditional_bc_direct_quartic_map('B_-1',make_bundle())
    done = apply_conditional_bc_channel_normalization('B_-1',signed)
    with pytest.raises(ValueError,match='already been applied'):
        apply_conditional_bc_channel_normalization('B_-1',done)


def test_missing_channel_rejected():
    signed = apply_conditional_bc_direct_quartic_map('B_-1',make_bundle())
    names = signed.names[:2]+signed.names[3:]
    arr = np.delete(signed.channels,2,axis=0)
    broken = replace(signed,names=names,channels=arr,quartic=arr.sum(axis=0))
    with pytest.raises(ValueError,match='missing or duplicated'):
        apply_conditional_bc_channel_normalization('B_-1',broken)


def test_inconsistent_quartic_total_rejected():
    signed=apply_conditional_bc_direct_quartic_map('C_-1',make_bundle('C_-1'))
    broken=replace(signed,quartic=signed.quartic+0.1)
    with pytest.raises(ValueError,match='channel sum disagree'):
        apply_conditional_bc_channel_normalization('C_-1',broken)


def test_controller_only_adds_normalization_in_mixed_after_sign():
    p=Path(__file__).resolve().parents[2]/'scripts'/'RunT3GeneralSequentialFitAblation.py'
    s=p.read_text(encoding='utf-8')
    tree=ast.parse(s)
    defs={n.name:n for n in tree.body if isinstance(n,ast.FunctionDef)}
    def call_count(name,node):
        return sum(isinstance(x,ast.Call) and isinstance(x.func,ast.Name)
                   and x.func.id==name for x in ast.walk(node))
    assert call_count('apply_conditional_bc_channel_normalization',defs['_mixed'])==1
    assert call_count('apply_conditional_bc_channel_normalization',defs['_nonmixed'])==0
    mixed=s[s.index('def _mixed('):s.index('def preflight(')]
    assert mixed.index('tensor=apply_conditional_bc_direct_quartic_map(') < mixed.index(
        'tensor=apply_conditional_bc_channel_normalization(') < mixed.index(
        'scalar=ScalarEndpoint(tensor,state)')


def matrix_record(A):
    A=np.asarray(A,complex)
    return {'real':A.real.tolist(),'imag':A.imag.tolist()}


def fake_signed_case():
    # Source input already has two direct channel signs flipped.
    fields={name:matrix_record(np.eye(3,dtype=complex)*(i+1)*0.1)
            for i,name in enumerate(NAMES)}
    total=sum((np.asarray(x['real']) for x in fields.values()),np.zeros((3,3)))
    return {'model':'B_-1','status':'Computed','cases':{'1111':{
        'quartic_provenance':{
           'bc_direct_quartic_channel_signs':{
             'general:lambdaHHdagS2S2':-1,
             'general:lambdaHHdagS1barS1bar':-1}},
        'per_general_interaction':[
          {'category':'quartic_decomposition','name':name,'c5_leading_log':value}
          for name,value in fields.items()],
        'c5_induced_wilson_log':matrix_record(total)
    }}}


def test_readonly_fixed_factor_checker_guard_for_other_channels():
    data=fake_signed_case()
    with pytest.raises(ValueError,match='additional nonzero HHSS channels'):
        fixed_spurion_from_signed_case(data,'B_-1')
    data['cases']['1111']['per_general_interaction'][-1]['c5_leading_log']=matrix_record(np.zeros((3,3)))
    data['cases']['1111']['c5_induced_wilson_log']=matrix_record(np.eye(3)*0.6)
    result=fixed_spurion_from_signed_case(data,'B_-1')
    expected=(math.sqrt(3)*.1 +2*.2+2*.3)*np.eye(3)/(2*math.log(100))
    np.testing.assert_allclose(result,expected,rtol=1e-14,atol=1e-15)
    compare=compare_fixed(result,result*1.0002)
    assert abs(compare['best_fit_real_diagnostic_only']-1.0002)<1e-10
    assert compare['best_fit_matrix_residual_diagnostic_only']<1e-12


def test_checker_requires_signed_unnormalized_provenance():
    data=fake_signed_case()
    data['cases']['1111']['quartic_provenance']['bc_channel_normalization_factors']=dict(CHANNEL_FACTORS)
    with pytest.raises(ValueError,match='already channel-normalized'):
        fixed_spurion_from_signed_case(data,'B_-1')
