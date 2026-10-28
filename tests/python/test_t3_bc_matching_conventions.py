"""Offline arithmetic, signed-map, data-provenance and output-isolation tests."""
from __future__ import annotations

import copy
import json
from pathlib import Path

import numpy as np
import pytest

from scripts import AnalyzeT3BCMatchingConventions as m
from scripts import RunT3ZeroExtraBoundaryBenchmark as b7


def evidence(*,wrong_sign=False, certified=False):
    blocks=[]
    for model in m.b9.MODELS:
        for i,(rgb,mat) in enumerate(b7.DIRECT.items()):
            blocks.append({'model':model,'block':('020020','020200')[i],
                           'rgbeta_couplings':[rgb], 'matchete_couplings':[mat],
                           'r_in_m':[['1']],
                           'lambda_rgbeta_from_lambda_matchete_if_assumptions_hold':
                              [['1' if wrong_sign else '-1']],
                           'assumed_relative_potential_sign':'-1',
                           'component_prefactor_one_ordered_orientation':'1/4',
                           'geometry_class':'DirectGeometryCandidate',
                           'physical_mapping_certified':False})
    p23={'format':'stage23-exact-source-prefactors-conditional-signed-maps-v1',
         'source_prefactors_verified':True, 'global_su2_component_basis_verified':False,
         'signed_physical_coupling_map_certified':certified,
         'totals':{'slots':50,'blocks':46},'blocks':blocks}
    p26={'format':'stage26-conditional-quartic-action-polynomials-v1',
         'source_relative_action_sign_supported':True,
         'signed_polynomial_identity_conditional_verified':True,
         'physical_coupling_mapping_certified':False,
         'models_checked':16,'slots_checked':50}
    return p23,p26


def make_scene(m1=1.3e6,m2=1e6,mu=0,lam=.1,qr=.01, name='general_quartics'):
    return {'scenario':name,
            'eft1_end':{'mu_gev':1e3,'mS1Sq_gev2':m1,'mS2Sq_gev2':m2,
                        'quadratic_mu2_gev2':b7.complex_json(mu),
                        'lambdaT3':b7.complex_json(lam),
                        'reused_quartics':{k:b7.complex_json(qr) for k in b7.DIRECT}}}


def mock_eval(scene,quartics,*,sign,unmixed):
    # Symmetric, nonlinear test function with mass, mixing, quartic and cross terms.
    q=sum(sign*complex(v) for v in quartics.values())
    mix=0 if unmixed else scene['mu']/1e6
    dmass=(scene['m1']-1.3e6)/1e6
    c=(1+dmass+scene['lambdaT3']*.1+mix+q+3*q*mix)
    return np.array([[c,.1*c,0],[.1*c,2*c,0],[0,0,.5*c]],dtype=complex)


def test_map_both_direct_BC():
    a,b=evidence()
    result=m.check_convention_evidence(a,b,models=m.b9.MODELS)
    assert set(result['direct_blocks'])==set(m.b9.MODELS)
    assert result['physical_mapping_certified'] is False
    assert all(block['matchete_from_rgbeta_conditional']=='-1'
               for model in result['direct_blocks'].values() for block in model.values())


def test_map_wrong_sign_rejected():
    a,b=evidence(wrong_sign=True)
    with pytest.raises(ValueError,match='conditional geometry/sign'):
        m.check_convention_evidence(a,b,models=('B_-1',))


def test_map_certification_drift_rejected():
    a,b=evidence(certified=True)
    with pytest.raises(ValueError,match='Stage23'):
        m.check_convention_evidence(a,b,models=('B_-1',))


def test_map_wrong_orientation_rejected():
    a,b=evidence()
    a['blocks'][0]['r_in_m']=[['-1']]
    with pytest.raises(ValueError,match='conditional geometry/sign'):
        m.check_convention_evidence(a,b,models=('B_-1',))


def test_map_missing_coupling_rejected():
    a,b=evidence()
    a['blocks']=a['blocks'][1:]
    with pytest.raises(ValueError,match='missing/ambiguous'):
        m.check_convention_evidence(a,b,models=('B_-1',))


def test_map_block_duplicate_rejected():
    a,b=evidence()
    a['blocks'].append(copy.deepcopy(a['blocks'][0]))
    with pytest.raises(ValueError,match='missing/ambiguous'):
        m.check_convention_evidence(a,b,models=('B_-1',))


def test_map_stage26_sign_drift_rejected():
    a,b=evidence()
    b['source_relative_action_sign_supported']=False
    with pytest.raises(ValueError,match='Stage26'):
        m.check_convention_evidence(a,b,models=('B_-1',))


def test_map_stage26_physical_cert_drift_rejected():
    a,b=evidence()
    b['physical_coupling_mapping_certified']=True
    with pytest.raises(ValueError,match='Stage26'):
        m.check_convention_evidence(a,b,models=('C_-1',))


def test_zero_matrix_rejected():
    with pytest.raises(ValueError,match='zero'):
        m.relative(np.zeros((3,3)), np.eye(3))


def test_decomposition_exact_closure_and_interference():
    z=make_scene(name='restricted_zero',qr=0)
    q=make_scene(m1=1.31e6,mu=1e4,qr=.02)
    e0=m.endpoint(z);e1=m.endpoint(q)
    no={k:0j for k in b7.DIRECT}
    a=mock_eval(e0,no,sign=+1,unmixed=True)
    b=mock_eval(e1,no,sign=+1,unmixed=True)
    c=mock_eval(e1,no,sign=+1,unmixed=False)
    d=mock_eval(e1,e1['quartics'],sign=+1,unmixed=True)
    e=mock_eval(e1,e1['quartics'],sign=+1,unmixed=False)
    out=m.decompose_effects(a,b,c,d,e)
    assert out['telescoping_residual_relative']<1e-12
    assert out['relative_mixing_at_zero_direct_quartics']>0
    assert out['relative_quartic_mass_mixing_interference']>0
    assert out['relative_total']>0
    assert out['relative_terms_are_norms_not_additive'] is True


def test_decomposition_nonsymmetric_rejected():
    x=np.eye(3,dtype=complex)
    broken=x.copy();broken[0,1]=.1
    with pytest.raises(ValueError,match='symmetric'):
        m.decompose_effects(x,broken,x,x,x)


def test_endpoints_preserve_mixed_mass():
    scene=make_scene(mu=10000+20000j)
    e=m.endpoint(scene)
    assert e['mu']==10000+20000j
    assert e['quartics']['lambdaHHdagS2S2']==.01


def test_endpoints_bad_scale_rejected():
    scene=make_scene()
    scene['eft1_end']['mu_gev']=100
    with pytest.raises(ValueError,match='1 TeV'):
        m.endpoint(scene)


def test_model_output_nonzero_mixing_and_minus_map():
    zero=make_scene(qr=0,name='restricted_zero')
    quartic=make_scene(mu=20500,m1=1.301e6,qr=.02)
    source={'chi2':1,'fit_accepted':True,'uv_on':False}
    report={'scenarios':{'restricted_zero':zero,'general_quartics':quartic}}
    result=m.model_output('B_-1',report,source,mock_eval)
    assert result['general_endpoint_mu2_over_min_m2']>0
    assert result['mappings']['identity_prior_preview']['effects']['relative_total']!=\
           result['mappings']['stage23_conditional_minus']['effects']['relative_total']
    assert result['physical_mixed_C5_certified'] is False
    assert result['sequential_eft1_general_C5_performed'] is False


def test_model_output_previous_preview_crosscheck():
    zero=make_scene(qr=0,name='restricted_zero')
    quartic=make_scene(mu=100,m1=1.3e6,qr=.01)
    source={'chi2':1,'fit_accepted':True,'uv_on':False}
    report={'scenarios':{'restricted_zero':zero,'general_quartics':quartic},
            'C5_preview':{'epsilon_C5_frobenius_uncertified':0.777}}
    with pytest.raises(ValueError,match='did not reproduce'):
        m.model_output('C_-1',report,source,mock_eval)


def test_model_output_rejects_nonzero_reference_mixing():
    zero=make_scene(qr=0,mu=1,name='restricted_zero')
    quartic=make_scene(mu=100,m1=1.3e6,qr=.01)
    source={'chi2':1,'fit_accepted':True,'uv_on':False}
    with pytest.raises(ValueError,match='restricted mixed'):
        m.model_output('B_-1',{'scenarios':{'restricted_zero':zero,
                                     'general_quartics':quartic}},source,mock_eval)


def test_model_output_optional_input_mixing():
    zero=make_scene(qr=0,name='restricted_zero')
    quartic=make_scene(mu=100,qr=.01)
    explicit=make_scene(mu=10000,qr=.01,name='general_quartics_and_mixing')
    source={'chi2':1,'fit_accepted':True,'uv_on':False}
    result=m.model_output('C_-1',{'scenarios':{'restricted_zero':zero,
                    'general_quartics':quartic}},source,mock_eval,mixing_scenario=explicit)
    assert result['explicit_mixing_scenario']['mapping_results'][
        'stage23_conditional_minus']['relative_change_from_quartic_only_endpoint']>0


def test_new_run_ids_only(tmp_path):
    dest=m.check_output('small_uv_running_5e-9','clean_v1',root=tmp_path)
    dest.mkdir(parents=True)
    with pytest.raises(ValueError,match='already exists'):
        m.check_output('small_uv_running_5e-9','clean_v1',root=tmp_path)
    with pytest.raises(ValueError,match='run ID'):
        m.check_output('small_uv_running_5e-9','../unsafe',root=tmp_path)


def test_run_id_path_stays_in_general_scalar(tmp_path):
    p=m.pilot_path('small_uv_frozen_5e-9','saved_v1','B_-1',root=tmp_path)
    assert p==tmp_path/'output/studies/general_scalar_extended/small_uv_frozen_5e-9/saved_v1/B_m1_paired_pilot.json'


def test_summary_only_reports_conditional():
    zero=make_scene(qr=0,name='restricted_zero')
    quartic=make_scene(mu=100,qr=.01)
    result=m.model_output('B_-1',{'scenarios':{'restricted_zero':zero,
                    'general_quartics':quartic}},
                    {'chi2':11,'fit_accepted':False,'uv_on':False},mock_eval)
    s=m.format_summary('large_uv_frozen_0p5','new',{'B_-1':result},mixed_run=False)
    assert 'Not a physical' in s
    assert 'NO' in s
    assert 'conditional' in s.lower()


def test_mass_matrix_json_roundtrip():
    z=np.eye(3,dtype=complex)+2j*np.eye(3,dtype=complex)
    p=m.matrix_json(z)
    assert np.allclose(np.asarray(p['real'])+1j*np.asarray(p['imag']),z)


def test_one_scale_eval_uses_actual_complex_mass_rotation(tmp_path, monkeypatch):
    """Numerical plumbing integration with a simple test expression, not Matchete physics."""
    import sys
    import types
    from types import SimpleNamespace

    top=types.ModuleType('Numerical')
    top.__path__=[]
    exp=types.ModuleType('Numerical.experimental')
    exp.__path__=[]
    module=types.ModuleType('Numerical.experimental.GeneralT3MixedScalarMassBasis')

    class FakeMassBasis:
        @staticmethod
        def build(model,m1,m2,mu):
            matrix=np.array([[m1,mu],[np.conjugate(mu),m2]],dtype=complex)
            eig,U=np.linalg.eigh(matrix)
            return SimpleNamespace(eigenvalues_sq=eig,U=U,degenerate=False)

    module.MixedDoubletMassBasis=FakeMassBasis
    cp=types.ModuleType('scripts.RunGeneralT3NumericalC5')
    cp.compile_coefficient=lambda expression:compile(expression,'<test symbolic>','eval')
    monkeypatch.setitem(sys.modules,'Numerical',top)
    monkeypatch.setitem(sys.modules,'Numerical.experimental',exp)
    monkeypatch.setitem(sys.modules,'Numerical.experimental.GeneralT3MixedScalarMassBasis',module)
    monkeypatch.setitem(sys.modules,'scripts.RunGeneralT3NumericalC5',cp)
    coeff=tmp_path/'synthetic_coefficient.inputform.txt'
    coeff.write_text(
        'hbar*(y1a*y2b+y1b*y2a)*(lambdaT3 + '
        'lambdaT3ExtraQ020020I1 + lambdaT3ExtraQ020200I1 + '
        'T3MEU11*T3MEU22 + MSX/MSY)/MF',encoding='utf-8')
    y=np.eye(3,dtype=complex)
    monkeypatch.setattr(m,'historical_yukawa_fermion',lambda source:
                        (np.array([1e5,1.1e5,1.2e5]),y,y))
    ev=m.make_matchete_evaluator('B_-1',{},coeff)
    scene={'m1':1.3e6,'m2':1e6,'mu':1e4+2e4j,'lambdaT3':.1+0j}
    q={k:.01+0j for k in b7.DIRECT}
    mixed_plus=ev(scene,q,sign=+1,unmixed=False)
    unmixed=ev(scene,q,sign=+1,unmixed=True)
    mixed_minus=ev(scene,q,sign=-1,unmixed=False)
    assert mixed_plus.shape==(3,3)
    assert np.isfinite(mixed_plus).all()
    assert np.linalg.norm(mixed_plus-unmixed)>0
    assert np.linalg.norm(mixed_plus-mixed_minus)>0
    assert np.allclose(mixed_plus,mixed_plus.T)
