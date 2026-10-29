"""Regression/contract tests for the B/C selected extra Wilson source batch."""
from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
from types import ModuleType, SimpleNamespace
import sys

import numpy as np
import pytest
import sympy as sp

from scripts import RunT3BCExtraWilsonSequential as new


TEMPLATE = {(1,3):-sp.Rational(1,4),(1,4):-sp.I/4,
            (2,3):-sp.I/4,(2,4):sp.Rational(1,4)}


def source_fixture(model='B_-1'):
    stage = new.STAGE_NAMES[model]
    f = sp.Symbol('MF')
    lam1,lam2 = sp.symbols('lambdaHHdagS1barS1bar lambdaHHdagS2S2')
    y1,y2 = sp.symbols('y1 y2')
    factor=sp.Integer(2) if model=='B_-1' else 2*sp.sqrt(3)/3
    beta=factor*(sp.conjugate(lam1)*sp.conjugate(y1)**2+
                 sp.conjugate(lam2)*sp.conjugate(y2)**2)/f
    stage17={
        'model':stage,'status':'ExtraEFT1ScalarPairMixingDerived',
        'source_convention':'lambda_Matchete=-lambda_RGBeta, same named scalar field basis',
        'full_eft1_wilson_running_implemented':False,
        'nonzero_weinberg_real_components':4,
        'components':[{'higgs_real_indices':list(pair),
                       'beta_16pi2_expression':sp.sstr(beta*v)}
                      for pair,v in TEMPLATE.items()],
    }
    cert={'status':'FullExtraWeinbergTensorProjectionVerified',
          'model':stage,'independent_components_verified':30,
          'mixed_components_matched_to_stage17':4,
          'all_component_residuals_zero':True,
          'beta_kappa_16pi2_expression':sp.sstr(beta)}
    norm={'status':'ProductionDirectC5NormalizationVerified',
          'direct_running_log_coefficient_equals_beta_kappa':True,
          'beta_kappa_16pi2':'4*lambdaT3*conjugate(y1)*conjugate(y2)/(3*MF)'}
    stage19={'model':stage,'status':'SelectedExtraWeinbergSourceTransported',
             'full_coupled_eft1_wilson_evolution_implemented':False,
             'finite_scalar_matching_implemented':False,'complete_two_loop_result':False,
             'projection_certificate':cert,'normalization_certificate':norm}
    seed={'Level':1,'IntegratedFields':['F'],'ActiveHeavyFields':['S1','S2']}
    wilson={'status':'Success'}; matched={'status':'Success'}
    return stage17,stage19,seed,wilson,matched


@pytest.fixture
def projection_stubs(monkeypatch):
    mod=ModuleType('Numerical.experimental.GeneralT3WeinbergProjection')
    mod._analytic_mixed_template=lambda: {key:complex(v) for key,v in TEMPLATE.items()}
    mod.verify_production_c5_normalization=lambda wilson,match: {
       'status':'ProductionDirectC5NormalizationVerified',
       'beta_kappa_16pi2':'4*lambdaT3*conjugate(y1)*conjugate(y2)/(3*MF)'}
    monkeypatch.setitem(sys.modules,mod.__name__,mod)
    return mod


@pytest.mark.parametrize('model',['B_-1','C_-1'])
def test_exact_saved_stage17_matches_saved_stage19_and_current_normalization(projection_stubs,model):
    stage17,stage19,seed,rge,match=source_fixture(model)
    cert,norm=new.certify_saved_selected_source(model,stage17,stage19,seed,rge,match)
    assert cert['independent_components_verified']==30
    assert norm['status']=='ProductionDirectC5NormalizationVerified'


@pytest.mark.parametrize('pair',[(1,3),(1,4),(2,3),(2,4)])
def test_tampered_stage17_component_fails_closed(projection_stubs,pair):
    s,h,seed,rge,match=source_fixture()
    row=next(x for x in s['components'] if x['higgs_real_indices']==list(pair))
    row['beta_16pi2_expression']='2*('+row['beta_16pi2_expression']+')'
    with pytest.raises(ValueError,match='differs from certified full tensor'):
        new.certify_saved_selected_source('B_-1',s,h,seed,rge,match)


def test_duplicate_higgs_support_rejected(projection_stubs):
    s,h,seed,rge,match=source_fixture()
    s['components'][1]['higgs_real_indices']=[1,3]
    with pytest.raises(ValueError,match='Stage17 Higgs-real tensor support|repeated'):
        new.certify_saved_selected_source('B_-1',s,h,seed,rge,match)


def test_missing_full_tensor_certificate_rejected(projection_stubs):
    s,h,seed,rge,match=source_fixture()
    h['projection_certificate']['independent_components_verified']=4
    with pytest.raises(ValueError,match=r'full SU\(2\) projection certificate'):
        new.certify_saved_selected_source('B_-1',s,h,seed,rge,match)


def test_changed_source_convention_rejected(projection_stubs):
    s,h,seed,rge,match=source_fixture()
    s['source_convention']='different sign'
    with pytest.raises(ValueError,match='Stage17 scalar-pair source'):
        new.certify_saved_selected_source('B_-1',s,h,seed,rge,match)


def test_saved_direct_normalization_mismatch_rejected(projection_stubs):
    s,h,seed,rge,match=source_fixture()
    h['normalization_certificate']['beta_kappa_16pi2']='3*lambdaT3*conjugate(y1)*conjugate(y2)/MF'
    with pytest.raises(ValueError,match='direct normalization do not agree'):
        new.certify_saved_selected_source('B_-1',s,h,seed,rge,match)


def test_historical_stage19_completeness_guard(projection_stubs):
    s,h,seed,rge,match=source_fixture()
    h['complete_two_loop_result']=True
    with pytest.raises(ValueError,match='Stage19 source status'):
        new.certify_saved_selected_source('B_-1',s,h,seed,rge,match)


def test_source_paths_are_BC_only():
    assert new.sources('B_-1')['stage17'].as_posix().endswith('/B_m1_extra_eft1_scalar_pair_beta.json')
    assert new.sources('C_-1')['seed'].as_posix().endswith('T3_C_alpha_m1/data/eft1_after_F_wilson_seed.json')
    with pytest.raises(ValueError,match='only B_-1'):
        new.sources('A_+0')


def test_frozen_selected_source_has_negative_log_sign(monkeypatch):
    from Numerical.experimental import GeneralT3ExtraWilsonTransport as trans
    from Numerical.experimental import GeneralT3WeinbergProjection as proj
    class Comp:
        def __init__(self,key):self.higgs_pair=key
        def evaluate(self,l2,l1,f1,f2):return complex(TEMPLATE[self.higgs_pair])*(2*(np.conj(l1)*f1+np.conj(l2)*f2))
    monkeypatch.setattr(trans,'compile_stage17_sources',lambda *a,**k:tuple(Comp(key) for key in TEMPLATE))
    monkeypatch.setattr(trans,'frozen_seed_flavor_contractions',lambda m,a,b:(a,a))
    monkeypatch.setattr(proj,'_analytic_mixed_template',lambda:{k:complex(v) for k,v in TEMPLATE.items()},raising=False)
    view=new.EFT1TrajectoryView(np.array([1e5,1e3]),(SimpleNamespace(
        lambdaHHdagS2S2=.03,lambdaHHdagS1barS1bar=.02),None))
    a=np.eye(3,dtype=complex)
    got=new.frozen_selected_source('B_-1',{}, {},{},view,np.array([2e5]*3),a,a)
    expected=2*(.02+.03)*np.log(1e3/1e5)/new.LOOP_FACTOR*a
    assert np.allclose(got,expected,rtol=1e-13,atol=1e-17)
    assert np.diag(got).real.max()<0


def test_frozen_selected_source_vanishes_at_zero_portal(monkeypatch):
    from Numerical.experimental import GeneralT3ExtraWilsonTransport as trans
    from Numerical.experimental import GeneralT3WeinbergProjection as proj
    class Comp:
        def __init__(self,key):self.higgs_pair=key
        def evaluate(self,l2,l1,f1,f2):return complex(TEMPLATE[self.higgs_pair])*(l1*f1+l2*f2)
    monkeypatch.setattr(trans,'compile_stage17_sources',lambda *a,**k:tuple(Comp(k) for k in TEMPLATE))
    monkeypatch.setattr(trans,'frozen_seed_flavor_contractions',lambda m,a,b:(a,a))
    monkeypatch.setattr(proj,'_analytic_mixed_template',lambda:{k:complex(v) for k,v in TEMPLATE.items()},raising=False)
    view=new.EFT1TrajectoryView(np.array([1e5,1e3]),(SimpleNamespace(
        lambdaHHdagS2S2=0,lambdaHHdagS1barS1bar=0),None))
    assert np.linalg.norm(new.frozen_selected_source('B_-1',{},{},{},view,
                         np.ones(3),np.eye(3),np.eye(3)))==0


def test_summary_writes_only_selected_models(tmp_path,monkeypatch):
    for model in new.MODELS:
        x={'format':new.FORMAT,'status':new.STATUS,'model':model,
            'complete_general_scalar_sequential_matching':False,
            'relative_extra_C5_vs_direct_at_scalar_threshold':.01,
            'relative_extra_C5_vs_frozen_extra':.03,
            'direct_only_masses_ev_descending':[.05,.01,.001],
            'updated_masses_ev_descending':[.05001,.01,.001],
            'source_integral_projection':{'projection_relative_residual':0.0}}
        (tmp_path/f'{model.replace("+","p").replace("-","m")}_extra_wilson.json').write_text(json.dumps(x))
    def load(p):return json.loads(Path(p).read_text())
    monkeypatch.setattr(new.neutrino,'load',load)
    def save(p,x):p.parent.mkdir(exist_ok=True,parents=True);p.write_text(json.dumps(x))
    monkeypatch.setattr(new.neutrino,'atomic_json',save)
    result=new.summary(tmp_path,new.MODELS)
    assert result['model_count']==2
    assert result['complete_general_scalar_sequential_matching'] is False
    assert (tmp_path/'summary.csv').is_file()


def test_summary_rejects_false_full_matching_claim(tmp_path,monkeypatch):
    x={'format':new.FORMAT,'status':new.STATUS,'model':'B_-1',
       'complete_general_scalar_sequential_matching':True}
    (tmp_path/'B_m1_extra_wilson.json').write_text(json.dumps(x))
    monkeypatch.setattr(new.neutrino,'load',lambda p:json.loads(p.read_text()))
    with pytest.raises(ValueError,match='incompatible result'):
        new.summary(tmp_path,('B_-1',))


def test_trajectory_view_exposes_saved_states():
    fake=(SimpleNamespace(mu_gev=1e5),SimpleNamespace(mu_gev=1e3))
    x=new.EFT1TrajectoryView(np.array([1e5,1e3]),fake)
    assert x.state_at_index(1).mu_gev==1e3


@pytest.mark.parametrize('model',['B_-1','C_-1'])
def test_repository_source_artifacts_if_present(model):
    """Source-level verification is performed against local checkout when present."""
    ps=new.sources(model)
    if any(not ps[k].is_file() for k in ('stage17','stage19','seed','wilson_rge','matching')):
        pytest.skip('full repository source artifacts not present in isolated environment')
    from Numerical.experimental.GeneralT3WeinbergProjection import (
        verify_production_c5_normalization,)
    objects={key:json.loads(ps[key].read_text(encoding='utf-8-sig'))
             for key in ('stage17','stage19','seed','wilson_rge','matching')}
    cert,norm=new.certify_saved_selected_source(model,
        objects['stage17'],objects['stage19'],objects['seed'],
        objects['wilson_rge'],objects['matching'])
    assert cert['independent_components_verified']==30
    assert norm['status']=='ProductionDirectC5NormalizationVerified'


def fake_compute_setup(tmp_path,monkeypatch,*,samples=65):
    model='B_-1';tag='B_m1'
    folder=tmp_path/'continuous_v1';out=folder/new.DEFAULT_OUTPUT_ID
    (folder/tag).mkdir(parents=True)
    (folder/new.DIRECT_OUTPUT_ID).mkdir()
    (folder/tag/'continuous_trajectory.json').write_text('{}')
    (folder/tag/'sidecar.json').write_text('{}')
    normal=np.diag([1.e-15,2.e-15,3.e-15]).astype(complex)
    direct_report={
        'format':new.direct.FORMAT,'status':new.direct.STATUS,
        'model':model,'saved_points':samples,
        'full_sequential_wilson_matching_performed':False,
        'full_general_scalar_matching_certified':False,
        'lambdaT3_log_integral':{'real':-0.012,'imag':0},
        'c5_at_scalar_threshold_gev_inverse':{'real':normal.real.tolist(),
                    'imag':normal.imag.tolist()},
        'takagi_masses_ev_descending':[3.e-2,2.e-2,1.e-2],
    }
    (folder/new.DIRECT_OUTPUT_ID/f'{tag}_direct_wilson.json').write_text(json.dumps(direct_report))
    filepaths={k:tmp_path/f'{k}.json' for k in ('stage17','stage19','seed',
                 'matching','wilson_rge','eft1_beta')}
    for p in filepaths.values():p.write_text('{}')
    filepaths['stage17'].write_text(json.dumps({'source_convention':'lambda_Matchete=-lambda_RGBeta, same named scalar field basis'}))
    monkeypatch.setattr(new,'sources',lambda model:filepaths)
    monkeypatch.setattr(new,'digest',lambda p:'hash')
    monkeypatch.setattr(new.neutrino,'verify_continuous',lambda model,path:{
        'saved_candidate_eft1_endpoint':'sidecar.json','uv_fermion_threshold_state':{},
        'eft1_scalar_endpoint_state':{}})
    monkeypatch.setattr(new.neutrino,'load',lambda p:json.loads(Path(p).read_text()))
    monkeypatch.setattr(new.neutrino,'decode_uv',lambda r:object(),raising=False)
    monkeypatch.setattr(new.neutrino,'decode_eft',lambda r:object(),raising=False)
    monkeypatch.setattr(new.neutrino,'complex_pair',lambda x,label:complex(x['real'],x['imag']),raising=False)
    monkeypatch.setattr(new.neutrino,'complex_matrix',lambda x,label:np.array(x['real'])+1j*np.array(x['imag']),raising=False)
    monkeypatch.setattr(new.neutrino,'matrix_payload',lambda x:{'real':np.asarray(x).real.tolist(),
                                                          'imag':np.asarray(x).imag.tolist()},raising=False)
    monkeypatch.setattr(new.neutrino,'validated_c5',lambda x:np.asarray(x,dtype=complex))
    def save(path,payload):
        path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(json.dumps(payload))
    monkeypatch.setattr(new.neutrino,'atomic_json',save)
    monkeypatch.setattr(new.direct,'heavy_inputs',lambda uv:(np.array([1e5]*3),
            np.eye(3,dtype=complex),np.eye(3,dtype=complex),0.0))
    monkeypatch.setattr(new,'certify_saved_selected_source',lambda *x,**k:({'status':'source'},
                                                                               {'status':'norm'}))
    view=new.EFT1TrajectoryView(np.array([1e5,1e3]),(object(),object()))
    monkeypatch.setattr(new,'running_scalar_background',lambda *a,**k:(view,{
       'lambdaT3_integral':-0.012,'source_portals_at_F':{},
       'source_portals_at_S':{},'eft1_nfev':7,'endpoint_parity':{}}))
    delta=np.diag([1.e-17,2.e-17,3.e-17]).astype(complex)
    stage18={'status':'ExtraEFT1LLHHTensorSourceIntegrated',
             'higgs_real_components':[]}
    projection={'projection_relative_residual':0,'status':'IntegratedExtraWeinbergSourceProjected'}
    monkeypatch.setattr(new,'selected_source_integral',lambda *a,**k:(delta,stage18,projection))
    monkeypatch.setattr(new,'frozen_selected_source',lambda *a,**k:delta*.9)
    def low_energy(model,c5,eft):
        masses=[3.e-2,2.e-2,1.e-2] if np.linalg.norm(c5-normal)<1e-19 else [3.03e-2,2.02e-2,1.01e-2]
        return {'c5_at_100gev_gev_inverse':{'real':c5.real.tolist(),
                'imag':c5.imag.tolist()},
                'takagi_masses_ev_descending':masses}
    monkeypatch.setattr(new.neutrino,'low_energy',low_energy,raising=False)
    return folder,out,normal,delta


def test_compute_only_adds_selected_extra_source(monkeypatch,tmp_path):
    folder,out,normal,delta=fake_compute_setup(tmp_path,monkeypatch)
    report=new.compute('B_-1',folder=folder,out=out)
    assert report['status']==new.STATUS
    target=np.array(report['direct_plus_selected_extra_C5_at_scalar_threshold_gev_inverse']['real'])
    assert np.allclose(target,normal.real+delta.real,rtol=1e-14,atol=1e-32)
    assert report['full_coupled_eft1_wilson_evolution_implemented'] is False
    assert report['complete_general_scalar_sequential_matching'] is False
    assert report['uv_running_repeated'] is False
    assert report['updated_masses_ev_descending'][0]==pytest.approx(.0303)
    assert (out/'B_m1_extra_wilson.json').is_file()


def test_compute_resume_checks_saved_provenance(monkeypatch,tmp_path):
    folder,out,normal,delta=fake_compute_setup(tmp_path,monkeypatch)
    first=new.compute('B_-1',folder=folder,out=out)
    second=new.compute('B_-1',folder=folder,out=out,resume=True)
    assert second==first
    with pytest.raises(ValueError,match='output exists or inputs changed'):
        new.compute('B_-1',folder=folder,out=out,resume=False)


def test_compute_blocks_mismatched_sample_count(monkeypatch,tmp_path):
    folder,out,normal,delta=fake_compute_setup(tmp_path,monkeypatch,samples=65)
    with pytest.raises(ValueError,match='previous direct source used different'):
        new.compute('B_-1',folder=folder,out=out,samples=33)


def test_compute_blocks_direct_full_matching_claim(monkeypatch,tmp_path):
    folder,out,normal,delta=fake_compute_setup(tmp_path,monkeypatch)
    direct_file=folder/new.DIRECT_OUTPUT_ID/'B_m1_direct_wilson.json'
    data=json.loads(direct_file.read_text())
    data['full_sequential_wilson_matching_performed']=True
    direct_file.write_text(json.dumps(data))
    with pytest.raises(ValueError,match='previous direct-only matching result incompatible'):
        new.compute('B_-1',folder=folder,out=out)


def test_compute_blocks_inconsistent_direct_log_integral(monkeypatch,tmp_path):
    folder,out,normal,delta=fake_compute_setup(tmp_path,monkeypatch)
    direct_file=folder/new.DIRECT_OUTPUT_ID/'B_m1_direct_wilson.json'
    data=json.loads(direct_file.read_text())
    data['lambdaT3_log_integral']={'real':-0.014,'imag':0}
    direct_file.write_text(json.dumps(data))
    with pytest.raises(ValueError,match='EFT1 source integration differs'):
        new.compute('B_-1',folder=folder,out=out)
