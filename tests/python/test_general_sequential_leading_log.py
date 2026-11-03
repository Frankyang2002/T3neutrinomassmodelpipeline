"""Small independent algebraic tests; run in the full local repository."""
from __future__ import annotations

from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0,str(ROOT))

from Numerical.running.GeneralT3SequentialLeadingLog import (
    MODELS, anchor_to_one_scale, decode, encode, geometry_contribution,
    general_extra_channels, number,
)


def test_all_16_models_explicit():
    assert len(MODELS) == len(set(MODELS)) == 16
    assert {'B_-1','C_-1','E_+2','A_-4'} <= set(MODELS)


def test_real_imag_fourth_derivative_independent_contributions():
    rows_r=[{'H1_H2_S1_S2_0based':[0,1,1,2], 'value_exact':'-sqrt(2)'}]
    rows_i=[{'H1_H2_S1_S2_0based':[0,1,1,2], 'value_exact':'1/2'}]
    q=geometry_contribution(n=8,real_rows=rows_r,imag_rows=rows_i,
                            z=2+3j,label='test')
    expected=-2*np.sqrt(2)+1.5
    assert q.shape==(4,4,8,8)
    assert np.isclose(q[0,1,1,2],expected)
    assert np.isclose(q[1,0,2,1],expected)
    assert np.isclose(q[0,1,2,1],expected)
    assert np.count_nonzero(q)==4


def test_full_sector_extras_with_real_complex_and_exclusion():
    z={'counts_H_BarH_S1_BarS1_S2_BarS2':[1,1,1,1,0,0],
       'source_registry_symbols':[{'rgbeta_symbol':'lambdaT3GeneralQ1'},
                                  {'rgbeta_symbol':'lambdaT3GeneralQ2'}],
       'exact_reference_real_component_tensor_basis':[
          {'id':'r1','Re_c_unit_tensor_components':[
              {'H1_H2_S1_S2_0based':[0,0,0,0],'value_exact':'4'}],
           'Im_c_unit_tensor_components':[]},
          {'id':'r2','Re_c_unit_tensor_components':[],
           'Im_c_unit_tensor_components':[
              {'H1_H2_S1_S2_0based':[1,1,1,1],'value_exact':'-3'}]},
       ]}
    geometry={'monomial_sectors':[z], 'registry_HHSS_extra_symbols':2}
    registry={'models':{'B_-1':{'total_additional_quartics':2}}}
    state={'lambdaT3GeneralQ1':2,
           'lambdaT3GeneralQ2':{'real':0,'imag':5}}
    names,parts,info=general_extra_channels('B_-1',registry,state,geometry=geometry)
    assert names==['lambdaT3GeneralQ1','lambdaT3GeneralQ2']
    assert parts[0][0,0,0,0]==8
    assert parts[1][1,1,1,1]==-15
    names,parts,_=general_extra_channels('B_-1',registry,state,
        excluded={'lambdaT3GeneralQ2'},geometry=geometry)
    assert names==['lambdaT3GeneralQ1']
    with pytest.raises(ValueError,match='missing GENERAL scalar coupling'):
        general_extra_channels('B_-1',registry,{},geometry=geometry)


def test_one_scale_anchor_no_double_counting():
    full=np.array([[3e-12,1j*1e-12,0],[1j*1e-12,2e-12,0],[0,0,1e-12]])
    log=np.eye(3,dtype=complex)*0.3e-12
    on=anchor_to_one_scale(full,log,wilson_running_on=True)
    off=anchor_to_one_scale(full,log,wilson_running_on=False)
    np.testing.assert_allclose(decode(on['c5_at_scalar_threshold']),full)
    np.testing.assert_allclose(decode(off['c5_at_scalar_threshold']),full-log)
    np.testing.assert_allclose(decode(on['c5_residual_matching']),full-log)
    assert on['source'].startswith('one_scale_minus_log')


def test_decode_and_number_validation():
    a=np.eye(3,dtype=complex)*(0.3+0.7j)
    assert np.allclose(decode(encode(a)),a)
    assert number({'real':-2.,'imag':0.5},'x')==-2+0.5j
    with pytest.raises(ValueError,match='nonfinite'):
        number(float('nan'),'x')


def test_malformed_general_tensor_rejected():
    with pytest.raises(ValueError,match='bad derivative component'):
        geometry_contribution(n=8,real_rows=[
            {'H1_H2_S1_S2_0based':[0,3,8,0],'value_exact':'1'}],
            imag_rows=[],z=1,label='invalid')


def test_one_scale_wilson_ablation_exact_on_and_off(monkeypatch):
    from scripts import RunT3GeneralSequentialFitAblation as ab
    target=np.array([[4,0.1j,0],[0.1j,3,0],[0,0,2]],complex)*1e-12
    leading=np.eye(3,dtype=complex)*0.2e-12
    bubble=np.eye(3,dtype=complex)*0.4e-12
    fake_result={
       'c5_eft1_leading_log':encode(leading),
       'c5_scalar_finite_bubble':encode(bubble),
       'per_channel':[], 'quartic_provenance':{},
       'matching_scale_cancellation_residual':1e-14,
       'weinberg_beta_projection_residual':0.,
       'weinberg_bubble_projection_residual':0.,
    }
    from Numerical.running import GeneralT3SequentialLeadingLog as module
    monkeypatch.setattr(ab, 'response',lambda *a,**k:fake_result)
    monkeypatch.setattr(ab, 'low_energy',lambda state,c5,**k:{
        'c5_at_low_scale':encode(c5),
        'masses_ev_ascending':[1.,2.,3.],
        'prediction':{'sin2_theta13':0.02},'chi2':0.})
    fake=SimpleNamespace(couplings={'mS1Sq':1e6,'mS2Sq':1e6},
        provenance={'general_all':True})
    cases=ab._final_cases('A_-4',target,SimpleNamespace(),fake,None,
                          uv_on=True,eft_on=False,uv_nfev=20,eft_nfev=0,
                          spinor_factor=1.,vev=246.22)
    assert set(cases)=={'1000','1001','1010','1011'}
    for stage in ('0','1'):
        np.testing.assert_allclose(decode(cases['101'+stage]['c5_at_scalar_threshold']),target)
        np.testing.assert_allclose(decode(cases['100'+stage]['c5_at_scalar_threshold']),target-leading)
        np.testing.assert_allclose(decode(cases['101'+stage]['c5_bookkeeping_residual_matching']),
                                   target-leading-bubble)


def test_candidate_quartic_mapping_keeps_extra(monkeypatch):
    from Numerical.running.GeneralT3SequentialLeadingLog import native_aliases
    import types
    # The live registry generator provides this exact shape; fixture avoids
    # requesting local Matchete files in the pure algebra test.
    parent=types.ModuleType('Numerical.lo_certification')
    coupling=types.ModuleType('Numerical.lo_certification.couplings')
    child=types.ModuleType('Numerical.lo_certification.couplings.extra_quartic_source_v6')
    child.model_registry_quartics=lambda model,reg:[
        {'Matchete_generator_candidate_symbol':'lambdaT3ExtraQ110200I1',
         'rgbeta_symbol':'lambdaT3GeneralQ110200I1'}]
    monkeypatch.setitem(sys.modules,'Numerical.lo_certification',parent)
    monkeypatch.setitem(sys.modules,'Numerical.lo_certification.couplings',coupling)
    monkeypatch.setitem(sys.modules,'Numerical.lo_certification.couplings.extra_quartic_source_v6',child)
    names=['lambdaT3','lambdaT3ExtraQ110200I1']
    vals,alias=native_aliases('A_-2',names,{'lambdaT3':0.1,
            'lambdaT3GeneralQ110200I1':0.25}, {})
    assert vals['lambdaT3ExtraQ110200I1']==0.25
    assert alias['lambdaT3ExtraQ110200I1']=='lambdaT3GeneralQ110200I1'


def test_dynamic_tree_contracts_all_flavour_axes_without_dropping_terms(monkeypatch):
    import types
    from collections import defaultdict
    from Numerical.running.GeneralT3DynamicFermionWilson import _dense_operator
    p=types.ModuleType('Numerical.running.GeneralT3WilsonPolynomialBatch')
    def orbit(k):
        i,a,j,b,x,y=k
        return sorted({(i,a,j,b,x,y),(j,b,i,a,x,y),
                       (i,a,j,b,y,x),(j,b,i,a,y,x)})
    def polynomial(raw):
        grouped=defaultdict(complex)
        for key,z in raw.items():
            k=min(orbit(key))
            grouped[k]+=z
        return {k:v/len(orbit(k)) for k,v in grouped.items()},dict(grouped)
    p.source_term_alignment=lambda op,row:{'flavour_axis_permutation_for_lepton_order':[1,0]}
    p.canonical_polynomial=polynomial
    p.distinct_orbit=orbit
    b=types.ModuleType('Numerical.running.GeneralT3OrderedTensor')
    b.scalar_basis=lambda model:SimpleNamespace(total_dimension=12)
    monkeypatch.setitem(sys.modules,'Numerical.running.GeneralT3WilsonPolynomialBatch',p)
    monkeypatch.setitem(sys.modules,'Numerical.running.GeneralT3OrderedTensor',b)
    m=np.zeros((3,3),complex);m[0,1]=2j
    op={'source_index':7,'components':[{'lepton_su2':[1,2],
         'scalar_real':[5,6],
         'ordered_flavour_coefficient_gev_inv':encode(m)}]}
    by,indices,meta=_dense_operator([op],[{'source_index':7}],model='A_-4')
    assert indices==(7,)
    assert by.shape==(1,2,3,2,3,8,8)
    # Matchete lepton axes [1,0] map input flavour [p,q]=[1,2] -> [2,1].
    assert by[0,0,1,1,0,0,1]==0.5j
    assert by[0,1,0,0,1,1,0]==0.5j
    assert np.count_nonzero(by)==4
    assert meta['operator_symmetrisation'].startswith('exact')
