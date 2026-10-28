"""Batch17 interface/optimizer tests; native Windows RGBeta running not mocked as physics."""
from __future__ import annotations
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace
import math
import numpy as np
import pytest
from scripts import FitT3BCUVHierarchyConstrained as mod


def source():
    return {'lambdaT3':5e-9,'cfg':{'base_state':{'ordinary':{
        'lambdaT3':{'real':5e-9,'imag':0},'mS1Sq':1e6,'mS2Sq':1e6}}}}


def options(**more):
    d=dict(factors=(.005,.01),sources=('historical_fit',),sobol_points=8,
      seed_count=2,max_nfev=20,local_fraction=.8,floor_fraction=.03,
      lambda_low=1e-9,lambda_high=.1,yukawa_abs_max=.15,
      ratio_max=.5,scalar_min=.5,scalar_max=5.)
    d.update(more);return d


def test_inventory_and_defaults():
    assert len(mod.NAMES)==19
    assert len(mod.b11.NAMES)==18
    assert mod.DEFAULT_FACTORS['B_-1']==(.005,.01,.015,.02)
    assert mod.DEFAULT_FACTORS['C_-1']==(.005,.01,.015)


@pytest.mark.parametrize('changes',[{'factors':(0,)},{'factors':(0.01,0.01)},
  {'sources':('nonsense',)},{'sobol_points':3},{'seed_count':0},
  {'max_nfev':0},{'local_fraction':0},{'floor_fraction':0},
  {'lambda_low':.2},{'lambda_high':1.2},
  {'yukawa_abs_max':0},{'ratio_max':1.0},{'scalar_max':.5}])
def test_invalid_options(changes):
    with pytest.raises(ValueError):mod.validate_options(**options(**changes))


def test_valid_options():
    mod.validate_options(**options())


def test_no_overwrite_or_traversal(tmp_path):
    dst=mod.destination(tmp_path,'my_run')
    dst.mkdir(parents=True)
    (dst/'original').write_text('preserve')
    with pytest.raises(ValueError):mod.destination(tmp_path,'my_run')
    assert (dst/'original').read_text()=='preserve'
    with pytest.raises(ValueError):mod.destination(tmp_path,'../bad')


def test_candidate_updates_both_lambda_representations_without_mutation(monkeypatch):
    s=source();prev=deepcopy(s)
    monkeypatch.setattr(mod.b11,'trial_record',lambda record,y:{**record,'yvalues':y.tolist()})
    x=np.r_[np.linspace(-.01,.01,18),-4.0]
    record,lam=mod.candidate_record(s,x,lambda_low=1e-9,lambda_high=.1)
    assert lam==pytest.approx(1e-4)
    assert record['lambdaT3']==pytest.approx(lam)
    assert record['cfg']['base_state']['ordinary']['lambdaT3']['real']==pytest.approx(lam)
    assert record['yvalues']==pytest.approx(x[:18])
    assert s==prev


@pytest.mark.parametrize('loglam',[-12.,0.,float('nan')])
def test_candidate_bounds_fail(loglam,monkeypatch):
    monkeypatch.setattr(mod.b11,'trial_record',lambda r,y:r)
    with pytest.raises(ValueError):mod.candidate_record(source(),np.r_[np.zeros(18),loglam],
                                           lambda_low=1e-9,lambda_high=.1)


def test_whitened_fit_and_correlated_covariance():
    target=SimpleNamespace(observable_names=mod.b11.FIELDS,
        central_values=np.array([2.,4.,6.,8.,10.]),
        covariance=np.diag([1.,4.,9.,16.,25.]))
    pred=dict(zip(mod.b11.FIELDS,[3.,6.,9.,12.,15.]))
    res=mod.whiten(target,pred)
    assert np.allclose(res,np.ones(5))
    assert float(res@res)==pytest.approx(5.)


def test_screen_penalty_nonzero_and_sensitive():
    inside=mod.screen_penalty({'rhoF':.25,'scalar_light_muS_tev':2,
         'scalar_heavy_muS_tev':3},ratio_max=.5,scalar_min=.5,scalar_max=5)
    outside=mod.screen_penalty({'rhoF':.7,'scalar_light_muS_tev':.3,
         'scalar_heavy_muS_tev':8},ratio_max=.5,scalar_min=.5,scalar_max=5)
    assert len(inside)==5 and np.all(inside>0)
    assert np.linalg.norm(outside)>np.linalg.norm(inside)
    assert np.isfinite(mod.screen_penalty({'rhoF':None},ratio_max=.5,
                                      scalar_min=.5,scalar_max=5)).all()


def target():
    return SimpleNamespace(observable_names=mod.b11.FIELDS,
       central_values=np.zeros(5),covariance=np.eye(5))


def test_screened_evaluator_wont_accept_mass_invalid_even_with_small_raw_chi2(monkeypatch):
    monkeypatch.setattr(mod,'candidate_record',lambda source,p,**kw:({},10**p[-1]))
    def eval_one(model,record,reg,bundle,compiled,target,**kw):
        assert kw['ratio_max']==.5
        if state['bad']:
            return dict(status='MassScreenFailAtMuF',screen_pass=False,chi2=None,
                rhoF=2.,scalar_light_muS_tev=None,scalar_heavy_muS_tev=None)
        return dict(status='CandidateOscillationEvaluated',screen_pass=True,
            chi2=5.,prediction=dict(zip(mod.b11.FIELDS,[1.]*5)),
            rhoF=.1,scalar_light_muS_tev=1.,scalar_heavy_muS_tev=3.)
    state={'bad':True}
    monkeypatch.setattr(mod.b16,'run_one',eval_one)
    obj=mod.ScreenedEvaluator('B_-1',{},None,None,None,target(),q=.1,mu=0,
              ratio_max=.5,scalar_min=.5,scalar_max=5.,lambda_low=1e-9,lambda_high=.1)
    x=np.r_[np.zeros(18),-4.]
    bad=obj.evaluate(x,label='bad')
    assert bad['chi2'] is None and obj.best is None and obj.match_calls==0
    state['bad']=False
    good=obj.evaluate(np.r_[np.full(18,.001),-4.],label='good')
    assert good['chi2']==pytest.approx(5.) and obj.best['chi2']==pytest.approx(5.)
    assert obj.screened==1 and obj.match_calls==1
    assert obj.calls==2
    assert obj.evaluate(x,label='cached')['chi2'] is None
    assert obj.calls==2


def test_screened_evaluator_chi2_inconsistency_fails_closed(monkeypatch):
    monkeypatch.setattr(mod,'candidate_record',lambda source,p,**kw:({},10**p[-1]))
    monkeypatch.setattr(mod.b16,'run_one',lambda *a,**kw:dict(
        status='CandidateOscillationEvaluated',screen_pass=True,
        chi2=0.,prediction=dict(zip(mod.b11.FIELDS,[1.]*5))))
    obj=mod.ScreenedEvaluator('C_-1',{},None,None,None,target(),q=.1,mu=0,
          ratio_max=.5,scalar_min=.5,scalar_max=5.,lambda_low=1e-9,lambda_high=.1)
    result=obj.evaluate(np.r_[np.zeros(18),-4],label='mismatch')
    assert result['status']=='NumericalTrialFailed' and obj.best is None
    assert 'does not agree' in result['reason']


def test_screened_evaluator_catches_exception_without_false_success(monkeypatch):
    monkeypatch.setattr(mod,'candidate_record',lambda *a,**kw:(_,)[1])
    obj=mod.ScreenedEvaluator('C_-1',{},None,None,None,target(),q=.1,mu=0,
          ratio_max=.5,scalar_min=.5,scalar_max=5.,lambda_low=1e-9,lambda_high=.1)
    out=obj.evaluate(np.r_[np.zeros(18),-4.],label='x')
    assert out['chi2'] is None and out['status']=='NumericalTrialFailed' and obj.failed==1


def test_seed_anchors_two_textures_and_tunable_loglambda(monkeypatch):
    monkeypatch.setattr(mod.b16,'seeded_yukawa',lambda *a,**kw:np.full(18,
       .025 if kw['source_name']=='historical_fit' else -.015))
    anchors=mod.seed_anchors('B_-1',source(),{'best_parameters':{}},
         factors=(.01,.02),sources=('historical_fit','saved_poor_fit'),
         lambda_low=1e-9,lambda_high=.1,yukawa_abs_max=.15,
         local_fraction=.8,floor_fraction=.03)
    assert len(anchors)==4
    for a in anchors:
        assert len(a['seed'])==19 and a['bounds'].shape==(2,19)
        assert np.all(a['bounds'][0]<a['seed']+1e-12)
        assert a['bounds'][0,-1]<a['bounds'][1,-1]
    assert anchors[0]['seed'][-1]==pytest.approx(math.log10(5e-5))
    assert anchors[1]['seed'][-1]==pytest.approx(math.log10(1.25e-5))


def test_summary_does_not_claim_certification():
    r={'model':'B_-1','run_id':'v1','status':'NoMassSafeNeutrinoFit',
       'mass_screen':{'rhoF_max':.5,'scalar_min_tev':.5,'scalar_max_tev':5.},
       'search':{'anchors':2,'sobol_points_per_anchor':8},
       'run_counts':{'mass_safe_c5_evaluations':0,'failed_numerical_evaluations':0},
       'best_mass_safe':None}
    s=mod.summary_markdown(r)
    assert 'No mass-safe point' in s and 'conditional' in s


def test_local_refinement_all_19_coordinates(monkeypatch):
    seen=[]
    def fake_eval(x,**kw):
        seen.append(x.copy())
        return {'status':'CandidateOscillationEvaluated','chi2':float(np.sum(x*x)),
          'whitened_residual':[float(x[0]),float(x[1]),float(x[2]),float(x[3]),float(x[4])]}
    evaluator=SimpleNamespace(evaluate=fake_eval,best={'chi2':2.})
    low=np.full(19,-1.);high=np.full(19,1.)
    anchor={'bounds':np.vstack((low,high)),'source':'historical_fit','factor':.01}
    entry={'parameters':dict(zip(mod.NAMES,np.full(19,.1))),'chi2':2.}
    out=mod.local_refinement(evaluator,entry,anchor,max_nfev=2)
    assert out['optimiser_nfev']>=1
    assert len(seen)>19  # finite-difference Jacobian of 19 variables


def test_full_run_contract_preserves_inputs_and_writes_candidate(tmp_path,monkeypatch):
    """Mock I/O and numerical services, not a native RGBeta correctness test."""
    import sys, json
    from types import ModuleType
    def file(path,body='{}'):
        path=tmp_path/path
        path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(body,encoding='utf8')
        return path
    original_paths=[file('old/fit.json','{"stay":true}'),
                    file('old/config.json','{"stay":true}'),
                    file('old/historic.json','{"stay":true}'),
                    file('old/historic_fit.json','{"stay":true}')]
    registry=file('RGE/general/generated/general_rgbeta_quartic_registry.json')
    beta_uv=file('RGE/UV.json');beta_eft=file('RGE/EFT.json')
    stage23=file('RGE/stage23.json');stage26=file('RGE/stage26.json')
    target_path=file('configs/target.json')
    matchete=file('c5/B_m1_final_one_loop_probe_weinberg_ordered_coefficient.inputform.txt','dummy')
    fake_state=ModuleType('Numerical.experimental.GeneralT3ScalarRGState')
    fake_state.load_registry=lambda p:{'test':'registry'}
    monkeypatch.setitem(sys.modules,fake_state.__name__,fake_state)
    fake_orch=ModuleType('Numerical.orchestration')
    fake_orch.RunSmallLambdaRefit=SimpleNamespace(TARGET=target_path)
    monkeypatch.setitem(sys.modules,'Numerical.orchestration',fake_orch)
    fake_fit=ModuleType('Numerical.fitting.OscillationFit')
    fake_fit.OscillationFitTarget=SimpleNamespace(from_json=lambda p:'target')
    monkeypatch.setitem(sys.modules,fake_fit.__name__,fake_fit)
    fake_c5=ModuleType('scripts.RunGeneralT3NumericalC5')
    fake_c5.DEFAULT_SOURCE=matchete.parent
    monkeypatch.setitem(sys.modules,fake_c5.__name__,fake_c5)
    historical=source()
    historical.update(cfg_path=original_paths[2],fit_path=original_paths[3])
    prior={'best_parameters':{n:.01 for n in mod.b11.NAMES}}
    settings={'quartic_fraction':.1,'mu2_fraction':0.,
              'saved_fit_path':original_paths[0],
              'saved_config_path':original_paths[1]}
    monkeypatch.setattr(mod.b14,'fit_source',lambda *a,**kw:(prior,historical,settings))
    monkeypatch.setattr(mod.b11,'beta_bundle',lambda *a,**kw:{
       'files':{'UV':beta_uv,'EFT1':beta_eft}})
    monkeypatch.setattr(mod.b11,'compile_source',lambda *a,**kw:'compiled')
    monkeypatch.setattr(mod.b10,'EVIDENCE23',stage23)
    monkeypatch.setattr(mod.b10,'EVIDENCE26',stage26)
    monkeypatch.setattr(mod.b10,'check_convention_evidence',lambda *a,**kw:{})
    monkeypatch.setattr(mod.b16,'seeded_yukawa',lambda *a,**kw:np.full(18,.01))
    monkeypatch.setattr(mod,'sobol_points',lambda anchor,n,**kw:[anchor['seed'].copy() for _ in range(n)])
    def fake_record(src,p,**kw):
        cfg=deepcopy(src['cfg'])
        cfg['base_state']['ordinary']['lambdaT3']={'real':float(10**p[-1]),'imag':0.}
        cfg['saved_yukawa']=p[:18].tolist()
        return {'cfg':cfg},float(10**p[-1])
    monkeypatch.setattr(mod,'candidate_record',fake_record)
    class FakeEvaluator:
        def __init__(self,*args,**kw):
            self.calls=0;self.screened=0;self.failed=0;self.best=None
        def evaluate(self,p,**kw):
            self.calls+=1;self.screened+=1
            r={'status':'CandidateOscillationEvaluated','screen_pass':True,'chi2':2.,
               'rhoF':.02,'scalar_light_muS_tev':2.,'scalar_heavy_muS_tev':3.,
               'lambdaT3_uv':float(10**p[-1]),
               'prediction':dict(zip(mod.b11.FIELDS,[1,2,3,4,5])),
               'parameters':dict(zip(mod.NAMES,p.astype(float).tolist()))}
            self.best=r
            return r
    monkeypatch.setattr(mod,'ScreenedEvaluator',FakeEvaluator)
    report=mod.run_model(model='B_-1',run_id='newfit',root=tmp_path,
         factors=(.01,),sources=('historical_fit',),sobol_n=2,max_nfev=2)
    assert report['status']=='CandidateConstrainedAccepted'
    dest=tmp_path/'output/studies/general_scalar_extended/small_uv_running_5e-9/newfit'
    assert (dest/'fit.json').is_file()
    assert (dest/'fitted_config.json').is_file()
    assert (dest/'summary.md').is_file()
    out=json.loads((dest/'fitted_config.json').read_text())
    assert out['constrained_general_scalar_uv_fit']['physical_matching_certified'] is False
    assert len(out['saved_yukawa'])==18
    assert all(p.read_text()=='{"stay":true}' for p in original_paths)
    with pytest.raises(ValueError,match='colliding'):
        mod.run_model(model='B_-1',run_id='newfit',root=tmp_path,
          factors=(.01,),sources=('historical_fit',),sobol_n=2,max_nfev=2)


def test_cli_main_wires_scalar_threshold_options_to_model(monkeypatch):
    """Regression: parser --scalar-*-tev must use Namespace.scalar_*_tev."""
    import sys
    captured=[]
    monkeypatch.setattr(mod,'run_model',lambda **kw: captured.append(kw))
    monkeypatch.setattr(sys,'argv',[
        'FitT3BCUVHierarchyConstrained',
        '--model','B_-1',
        '--run-id','b_uv_constrained_pilot_v1',
        '--seed-source','historical_fit',
        '--seed-factors','0.005,0.01,0.02',
        '--sobol-points','8',
        '--seed-count','1',
        '--max-nfev','20',
        '--rho-max','0.4',
        '--scalar-min-tev','0.75',
        '--scalar-max-tev','4.5',
    ])
    mod.main()
    assert len(captured)==1
    got=captured[0]
    assert got['model']=='B_-1'
    assert got['run_id']=='b_uv_constrained_pilot_v1'
    assert got['sources']==('historical_fit',)
    assert got['factors']==pytest.approx((.005,.01,.02))
    assert got['sobol_n']==8 and got['seed_count']==1 and got['max_nfev']==20
    assert got['ratio_max']==pytest.approx(.4)
    assert got['scalar_min']==pytest.approx(.75)
    assert got['scalar_max']==pytest.approx(4.5)


def test_cli_main_default_scalar_threshold_options(monkeypatch):
    import sys
    captured=[]
    monkeypatch.setattr(mod,'run_model',lambda **kw: captured.append(kw))
    monkeypatch.setattr(sys,'argv',[
        'FitT3BCUVHierarchyConstrained',
        '--model','C_-1',
        '--run-id','c_uv_constrained_pilot_v1',
    ])
    mod.main()
    assert len(captured)==1
    assert captured[0]['scalar_min']==pytest.approx(.5)
    assert captured[0]['scalar_max']==pytest.approx(5.)
    assert captured[0]['ratio_max']==pytest.approx(.5)
    assert captured[0]['sources']==mod.b16.SOURCES
