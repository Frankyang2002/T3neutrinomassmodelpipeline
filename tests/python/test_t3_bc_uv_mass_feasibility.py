"""Batch 15 offline contract tests. No synthetic result is called a live RGBeta run."""
from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace
import json
import math

import numpy as np
import pytest

from scripts import StudyT3BCUVScalarMassFeasibility as mod


@pytest.mark.parametrize('arg',["",'0,1','-1,1','nan,1','1,1','1,inf','1,2,3,4,5,6,7,8,9,10,11,12,13','bad,1'])
def test_bad_grid_fails(arg):
    with pytest.raises(ValueError):mod.parse_factors(arg,'test')


def test_grid_sorted_replay_present():
    assert mod.parse_factors('1,0.25,0.05','Y')==(0.05,0.25,1.)


def test_candidate_record_copies_and_scales(monkeypatch):
    saved={'best_parameters':{key:1+i/10 for i,key in enumerate(mod.b11.NAMES)}}
    ordinary={'mS1Sq':1.44e6,'mS2Sq':1.e6,'y1':{'real':[[0.]*3]*3},
              'y2':{'real':[[0.]*3]*3}}
    source={'cfg':{'base_state':{'ordinary':ordinary}},'lambdaT3':5e-9}
    def trial_record(src,values):
        row=deepcopy(src)
        row['fit_yukawa_vector']=values.copy()
        return row
    monkeypatch.setattr(mod.b11,'trial_record',trial_record)
    got=mod.candidate_record(source,saved,.25,2.)
    assert np.allclose(got['fit_yukawa_vector'],np.asarray(list(saved['best_parameters'].values()))*.25)
    assert got['cfg']['base_state']['ordinary']['mS1Sq']==1.44e6*4
    assert got['cfg']['base_state']['ordinary']['mS2Sq']==4.e6
    assert source['cfg']['base_state']['ordinary']['mS1Sq']==1.44e6
    assert source['cfg']['base_state']['ordinary']['mS2Sq']==1.e6


def _endpoint(s1,s2,f=100.,*,positive=True):
    return {'scalar_masses_gev':[s1*1000 if positive else None,
                                  s2*1000 if positive else None],
            'positive_scalar_masses':positive,'min_fermion_mass_gev':f*1000,
            'scalar_eigen_mass_squared_gev2':[(s1*1000)**2,(s2*1000)**2]}


def test_hierarchy_and_target_window():
    uv={'final':_endpoint(2.,3.,100.)}
    eft={'final':_endpoint(1.,2.,100.)}
    c=mod.classify(uv,eft,ratio_max=.5,scalar_min_tev=.5,scalar_max_tev=5.)
    assert c['screen_pass']
    assert c['ratio_muF']==pytest.approx(.03)
    assert c['reason']=='CandidatePassedScreenOnly'


def test_ordering_failure_not_screen_pass():
    c=mod.classify({'final':_endpoint(90.,120.,100.)},None,
                   ratio_max=.5,scalar_min_tev=.5,scalar_max_tev=5.)
    assert not c['ordering_at_muF'] and not c['screen_pass']
    assert c['reason']=='ScalarAtOrAboveFAtMuF'


def test_hierarchy_margin_distinct_from_ordering():
    c=mod.classify({'final':_endpoint(40.,60.,100.)},{'final':_endpoint(1.,2.)},
                   ratio_max=.5,scalar_min_tev=.5,scalar_max_tev=5.)
    assert c['ordering_at_muF'] and not c['strict_hierarchy_at_muF']
    assert not c['screen_pass'] and c['reason']=='HierarchyMarginNotReached'


def test_target_mass_window_independent():
    c=mod.classify({'final':_endpoint(10.,15.,100.)},{'final':_endpoint(10.,20.)},
                   ratio_max=.5,scalar_min_tev=.5,scalar_max_tev=5.)
    assert c['strict_hierarchy_at_muF'] and not c['near_1tev_at_muS']
    assert not c['screen_pass']


def test_scalar_tachyonic_fails_closed():
    c=mod.classify({'final':_endpoint(0.,1.,positive=False)},None,
                   ratio_max=.5,scalar_min_tev=.5,scalar_max_tev=5.)
    assert not c['ordering_at_muF']
    assert c['reason']=='NonPositiveScalarAtMuF'


def test_beta_inclusion_exclusion_exact(monkeypatch):
    """Use the real beta-sector code, with a synthetic algebraic beta evaluator."""
    import sys
    from types import ModuleType
    # Fake the three imports performed by the function without touching real RG.
    class FakeState:
        baseline=SimpleNamespace()
        extra={}
    mock_mod=ModuleType('Numerical.experimental.GeneralT3ScalarRGState')
    mock_mod.unpack_general_scalar_state=lambda *a,**kw: FakeState()
    eval_mod=ModuleType('Numerical.core.RGBetaEvaluator')
    eval_mod.state_environment=lambda x:{'y1':np.eye(3),'y2':np.eye(3)*2,'MF':np.eye(3)*1e5}
    def evaluate(expr,env):
        sy=np.sum(np.abs(env['y1'])**2)+np.sum(np.abs(env['y2'])**2)
        sm=np.sum(np.abs(env['MF'])**2)
        # Non-trivial background, pure Y, pure MF, and mixed Y×MF.
        return 3e6+sy*4e4+sm*.001+sy*sm*1e-5 if expr=='one' else 6e6+sm*.002+sy*sm*3e-5
    eval_mod.evaluate_inputform_expression=evaluate
    monkeypatch.setitem(sys.modules,'Numerical.experimental.GeneralT3ScalarRGState',mock_mod)
    monkeypatch.setitem(sys.modules,'Numerical.core.RGBetaEvaluator',eval_mod)
    layout=SimpleNamespace(size=4)
    state=np.asarray([0.,0.,0.,0.,0.,0.])
    result=mod.beta_sector_at_uv(state,layout,{},
                 {'mS1Sq':'one','mS2Sq':'two'},scale=1e7)
    for terms in result['terms'].values():
        summation=sum(terms[z] for z in ('no_y_no_MF','y_dependence_at_MF0',
                    'MF_dependence_at_y0','joint_y_MF_dependence'))
        assert summation==pytest.approx(terms['full_d_m2_d_lnmu_gev2'])
        assert terms['joint_y_MF_dependence']>0
    assert 'NOT unique' in result['interpretation']


def _saved():
    return {'chi2_final':17000.,'best_parameters':{n:1. for n in mod.b11.NAMES},
            'final_diagnostics':{'MF_muF_gev':[80e3,90e3,110e3]}}


def _settings(tmp_path):
    fit=tmp_path/'saved_fit.json';cfg=tmp_path/'saved_config.json'
    fit.write_text('{}',encoding='utf-8');cfg.write_text('{}',encoding='utf-8')
    return {'quartic_fraction':.1,'mu2_fraction':0.,'saved_fit_path':fit,
            'saved_config_path':cfg,'verified_source_files':7}


def test_mocked_scan_skips_eft_if_inverted(monkeypatch,tmp_path):
    monkeypatch.setattr(mod,'candidate_record',lambda src,saved,y,s:{'y':y,'s':s})
    monkeypatch.setattr(mod.b11,'beta_bundle',lambda *a,**kw:{'UV':{}})
    def uv(model,trial,registry,bundle,**kw):
        f=trial['y'];s=trial['s']
        heavy=140. if f==1 else 2.
        return {'initial':_endpoint(1.*s,1.2*s),
                'final':_endpoint(heavy*.7,heavy,100.),
                'uv_nfev':11,'first_vector':np.zeros(2),
                'final_vector':np.zeros(2),'layout':SimpleNamespace(size=0)}
    # Avoid exact saved-mass verification and beta baseline by intercepting the
    # baseline native services in the probe fixture only.
    class FakeState:
        baseline=SimpleNamespace(MF=np.diag([80e3,90e3,110e3]))
    import sys
    from types import ModuleType
    mock=ModuleType('Numerical.experimental.GeneralT3ScalarRGState')
    mock.unpack_general_scalar_state=lambda *a,**kw:FakeState()
    monkeypatch.setitem(sys.modules,'Numerical.experimental.GeneralT3ScalarRGState',mock)
    monkeypatch.setattr(mod,'beta_sector_at_uv',lambda *a,**kw:{'terms':{},'mu_gev':kw['scale']})
    monkeypatch.setattr(mod,'evaluate_uv',uv)
    called=[]
    def eft(model,uv,trial,registry,bundle,**kw):
        called.append(trial['y'])
        return {'final':_endpoint(1.,2.),'eft1_nfev':8}
    monkeypatch.setattr(mod,'evaluate_eft1',eft)
    scan=mod.model_scan(model='B_-1',source={},saved=_saved(),
         settings=_settings(tmp_path),registry={},yukawa_factors=(.1,1.),
         mass_factors=(1.,),ratio_max=.5,scalar_min_tev=.5,scalar_max_tev=5.,grid=3)
    assert len(scan['trials'])==2
    assert scan['trials'][0]['ratio_muF']==1.4
    assert scan['trials'][0]['eft1_nfev']==0
    assert scan['trials'][1]['screen_pass'] is True
    assert called==[.1]


def test_execute_new_path_and_readonly(monkeypatch,tmp_path):
    prior=tmp_path/'historical.json';prior.write_text('unchanged',encoding='utf-8')
    monkeypatch.setattr(mod.b7,'read_json',lambda p:{'model_count':16})
    monkeypatch.setattr(mod.b14,'fit_source',lambda root,model,run_id,registry:
                (_saved(),{},_settings(tmp_path)))
    monkeypatch.setattr(mod,'model_scan',lambda **k:{'model':k['model'],
              'source_fit_chi2':17000.,'n_trials':1,
              'n_candidate_screen_passes':0,'candidate_pass_coordinates':[],
              'beta_decomposition':{},
              'trials':[{'model':k['model'],'yukawa_factor':1.,
                   'uv_scalar_mass_factor':1.,'ratio_muF':1.5,'status':'CandidateScreenFail'}]})
    r=mod.execute(run_id='new_run',models=('B_-1',),
        yukawa_factors=(1.,),mass_factors=(1.,),root=tmp_path)
    dest=tmp_path/'output/studies/general_scalar_extended/diagnostics/new_run'
    assert (dest/'summary.md').exists()
    assert (dest/'scan.csv').exists()
    assert (dest/'screen.json').exists()
    assert r['historical_outputs_overwritten'] is False
    assert prior.read_text()=='unchanged'
    with pytest.raises(ValueError,match='output already exists'):
        mod.execute(run_id='new_run',root=tmp_path,models=('B_-1',))


def test_validate_before_registry_io(tmp_path):
    with pytest.raises(ValueError,match='model list'):
        mod.execute(run_id='run',models=('A_+0',),root=tmp_path)
    with pytest.raises(ValueError,match='invalid hierarchy'):
        mod.execute(run_id='run',models=('B_-1',),ratio_max=1.0,root=tmp_path)
    with pytest.raises(ValueError,match='grid must'):
        mod.execute(run_id='run',models=('B_-1',),grid=2,root=tmp_path)


def test_summary_is_not_fit_label():
    d={'screen':{'ratio_max':.5,'scalar_min_tev':.5,'scalar_max_tev':5.},
       'models':{'B_-1':{'beta_decomposition':{},'trials':[
           {'yukawa_factor':1.,'uv_scalar_mass_factor':1.,'ratio_muF':1.2,
            'status':'CandidateScreenFail'}]}}}
    t=mod.summary_markdown(d)
    assert 'not a refit' in t.lower()
    assert 'chi2 is NOT recomputed' in t
    assert 'EFT1 is deliberately skipped' in t
