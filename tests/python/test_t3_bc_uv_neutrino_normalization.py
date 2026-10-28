"""Batch 16 offline unit/contract tests; mock integrations do not certify physics."""
from __future__ import annotations

from copy import deepcopy
from pathlib import Path
import math
from types import SimpleNamespace, ModuleType
import sys

import numpy as np
import pytest

from scripts import ProbeT3BCUVNeutrinoNormalization as mod


def _source():
    return {'lambdaT3':5e-9, 'cfg':{'base_state':{'ordinary':{
        'lambdaT3':{'real':5e-9,'imag':0},'mS1Sq':1.2e6,'mS2Sq':1e6}}}}


def test_lambda_scales_both_representations_and_preserves_historical():
    original=_source(); same=deepcopy(original)
    changed,l=mod.scaled_lambda_source(original,yukawa_factor=0.01,lambda_multiplier=1)
    assert l==pytest.approx(5e-5)
    assert changed['lambdaT3']==pytest.approx(l)
    assert changed['cfg']['base_state']['ordinary']['lambdaT3']['real']==pytest.approx(l)
    assert original==same


@pytest.mark.parametrize('f,k', [(0,1),(-1,1),(float('nan'),1),(0.01,0),
                                  (0.01,float('inf')),(0.0001,100)])
def test_bad_lambda_factors_rejected(f,k):
    with pytest.raises(ValueError):
        mod.scaled_lambda_source(_source(),yukawa_factor=f,lambda_multiplier=k)


def test_yukawa_source_saved_is_ordered_and_scaled():
    names=mod.b11.NAMES
    saved={'best_parameters':{name:float(i+1) for i,name in enumerate(names[::-1])}}
    out=mod.seeded_yukawa(saved,{},source_name='saved_poor_fit',factor=.005)
    assert np.allclose(out,np.array([saved['best_parameters'][n] for n in names])*.005)


def test_historical_source_uses_original_seed(monkeypatch):
    historical={};saved={'best_parameters':{n:0. for n in mod.b11.NAMES}}
    monkeypatch.setattr(mod.b11,'seed_vector',lambda s:np.arange(18,dtype=float))
    out=mod.seeded_yukawa(saved,historical,source_name='historical_fit',factor=.01)
    assert np.allclose(out,np.arange(18)*.01)


@pytest.mark.parametrize('raw',['','1,1','0.1,0','-1,1','nan,1','inf,1','X,1'])
def test_bad_csv_rejected(raw):
    with pytest.raises(ValueError):mod.parse_positive_csv(raw,'x')


def test_csv_does_not_require_one_and_order_preserved():
    assert mod.parse_positive_csv('0.05,0.1,10','x')==(.05,.1,10.)


def test_screen_pass_and_fail_propagate():
    uv={'positive_scalar_masses':True,'scalar_masses_gev':[1000.,2500.],
        'min_fermion_mass_gev':1e5}
    eft={'positive_scalar_masses':True,'scalar_masses_gev':[1100.,2000.]}
    ok=mod.check_mass_screen(uv,eft,ratio_max=.5,
                       scalar_min_tev=.5,scalar_max_tev=5.)
    bad=mod.check_mass_screen(uv,{**eft,'scalar_masses_gev':[1000.,8000.]},
                       ratio_max=.5,scalar_min_tev=.5,scalar_max_tev=5.)
    assert ok['screen_pass'] and not bad['screen_pass']


def test_lambda_estimate_is_not_an_optimisation():
    r={'prediction':{'delta_m3l_sq_ev2':0.002511/16}}
    assert mod.estimate_lambda_multiplier(r)==pytest.approx(4.)
    assert mod.estimate_lambda_multiplier({'prediction':None}) is None
    assert mod.estimate_lambda_multiplier({'prediction':{'delta_m3l_sq_ev2':0}}) is None


def test_immutable_destination(tmp_path):
    out=mod.destination(tmp_path,'candidate_test')
    out.mkdir(parents=True)
    (out/'keep').write_text('safe')
    with pytest.raises(ValueError):mod.destination(tmp_path,'candidate_test')
    assert (out/'keep').read_text()=='safe'
    with pytest.raises(ValueError):mod.destination(tmp_path,'../unsafe')


def test_default_model_factors_match_batch15_passes():
    assert mod.FACTORS['B_-1']==(.005,.01,.015,.02)
    assert mod.FACTORS['C_-1']==(.005,.01,.015)


def test_summary_says_not_refit_and_has_three_statuses():
    r={'screen':{'rho_max':.5,'scalar_min_tev':.5,'scalar_max_tev':5.},
       'models':{'B_-1':{'best':None,'trials':[{'yukawa_source':'saved_poor_fit',
       'yukawa_factor':.01,'lambda_multiplier':1.,'lambdaT3_uv':5e-5,
       'rhoF':.03,'chi2':None,'status':'MassScreenFailAtMuS'}]}}}
    text=mod.summary_markdown(r)
    assert 'No Yukawa or lambda optimisation' in text
    assert 'MassScreenFailAtMuS' in text
    assert '\n' in text


def test_run_one_skips_c5_when_fermion_hierarchy_invalid(monkeypatch):
    # Minimal numerical mock: the first threshold alone is rejected.
    fake=ModuleType('Numerical.experimental.GeneralT3ScalarRGState')
    fake.unpack_general_scalar_state=lambda *a,**kw: SimpleNamespace(baseline=SimpleNamespace())
    fake.project_general_after_fermion=lambda *a,**kw: (_ for _ in ()).throw(AssertionError('must skip'))
    fake.pack_general_scalar_state=lambda *a,**kw: (_ for _ in ()).throw(AssertionError('must skip'))
    monkeypatch.setitem(sys.modules,fake.__name__,fake)
    mock=SimpleNamespace(size=3,representation=SimpleNamespace())
    monkeypatch.setattr(mod.b9,'prepare_initial',lambda *a,**kw:(np.zeros(5),mock,{},0))
    monkeypatch.setattr(mod.b7,'rhs_at',lambda *a,**kw:np.zeros(5))
    monkeypatch.setattr(mod.b7,'integrate_segment',lambda *a,**kw:
                SimpleNamespace(y=np.zeros((5,3)),nfev=10))
    monkeypatch.setattr(mod.b15,'scalar_spectrum',lambda *a,**kw:
                {'scalar_masses_gev':[100e3,140e3],
                 'min_fermion_mass_gev':100e3,'mu2_abs_gev2':0,
                 'positive_scalar_masses':True})
    output=mod.run_one('B_-1',{},None,{'UV':{}},None,None,
         quartic_fraction=.1,mu2_fraction=0,ratio_max=.5,
         scalar_min_tev=.5,scalar_max_tev=5.)
    assert output['status']=='MassScreenFailAtMuF'
    assert output['match_performed'] is False
    assert output['chi2'] is None


def test_reject_model_and_collision_before_any_integration(tmp_path):
    dest=tmp_path/'output/studies/general_scalar_extended/diagnostics/collision'
    dest.mkdir(parents=True)
    with pytest.raises(ValueError,match='models must'):
        mod.execute(run_id='one',models=('A_+0',),root=tmp_path)
    with pytest.raises(ValueError,match='run ID already exists'):
        mod.execute(run_id='collision',models=('B_-1',),root=tmp_path)


def test_run_one_success_path_calls_c5_only_after_two_threshold_screens(monkeypatch):
    # Full interface contract through fake numerical services; no physical RG.
    fake=ModuleType('Numerical.experimental.GeneralT3ScalarRGState')
    uvbaseline=SimpleNamespace(MF=np.diag([80e3,100e3,110e3]))
    eftbaseline=SimpleNamespace()
    fake.unpack_general_scalar_state=lambda v,layout,registry,**kw: (
        SimpleNamespace(baseline=uvbaseline) if kw['uv'] else
        SimpleNamespace(baseline=eftbaseline))
    fake.project_general_after_fermion=lambda state,registry: SimpleNamespace()
    eftlayout=SimpleNamespace(size=3,names=('baseline','extra'),
                              representation=SimpleNamespace())
    fake.pack_general_scalar_state=lambda state,registry:(np.zeros(3),eftlayout)
    monkeypatch.setitem(sys.modules,fake.__name__,fake)

    uvlayout=SimpleNamespace(size=3,names=('baseline','extra'),
                             representation=SimpleNamespace())
    monkeypatch.setattr(mod.b9,'prepare_initial',lambda *a,**kw:(np.zeros(5),uvlayout,{},0))
    monkeypatch.setattr(mod.b7,'rhs_at',lambda *a,**kw:np.zeros(5))
    nseg=[]
    def integrate(y0, rhs, start, stop, **kw):
        nseg.append((start,stop))
        return SimpleNamespace(y=np.zeros((5,3)),nfev=19)
    monkeypatch.setattr(mod.b7,'integrate_segment',integrate)
    monkeypatch.setattr(mod.b15,'scalar_spectrum',lambda vec,layout,reg,**kw: {
        'scalar_masses_gev':[1300.,2100.],
        'min_fermion_mass_gev':80e3,'mu2_abs_gev2':0,
        'positive_scalar_masses':True})
    monkeypatch.setattr(mod.b7,'check_payload',lambda *a,**kw:{'dummy':'beta'})
    monkeypatch.setattr(mod.b7,'snapshot',lambda *a,**kw:{'mS1Sq_gev2':1e6})
    monkeypatch.setattr(mod.b11,'heavy_flavour_at_muF',lambda uv:(
        np.array([80e3,100e3,110e3]),np.eye(3),np.eye(3),1e-10))
    calls=[]
    def c5(*args):
        calls.append('c5');return np.eye(3)*1e-15
    monkeypatch.setattr(mod.b11,'c5_at_endpoint',c5)
    fit=SimpleNamespace(observable_names=('delta_m21_sq_ev2','delta_m3l_sq_ev2',
        'sin2_theta12','sin2_theta13','sin2_theta23'),
        prediction=(7.537e-5,2.511e-3,.3088,.02248,.47),chi2=0.02)
    def low(c5,eft,target,**kwargs):
        calls.append('low');return fit,SimpleNamespace(K=c5)
    monkeypatch.setattr(mod.b11,'low_energy_fit',low)
    from scripts import AnalyzeT3BCMatchingConventions as b10
    monkeypatch.setattr(b10,'matrix_json',lambda arr:{'shape':list(arr.shape)})
    row=mod.run_one('B_-1',{'new':('extra',),'reused':(),
                    'cfg':{'numerical':{'vev_gev':246.22}}},
                    {},{'UV':{},'payloads':{'EFT1':{}}},'compile',None,
                    quartic_fraction=.1,mu2_fraction=0.,
                    ratio_max=.5,scalar_min_tev=.5,scalar_max_tev=5.)
    assert nseg==[(1e7,1e5),(1e5,1e3)]
    assert calls==['c5','low']
    assert row['screen_pass'] and row['status']=='CandidateOscillationEvaluated'
    assert row['chi2']==pytest.approx(.02)
    assert row['absolute_neutrino_fit_pass']
    assert row['match_performed']
