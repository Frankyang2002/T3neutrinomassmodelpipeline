"""Offline deterministic campaign logic tests; no Wolfram or numerical fit needed."""
from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from scripts import RunT3ReducedScalarSobolComparison as mod


def test_power_two_sobol_and_seed_first():
    seed=np.zeros(18)
    lo=np.full(18,-1.)
    hi=np.ones(18)
    first=mod.initial_trials(seed,lo,hi,8,123)
    again=mod.initial_trials(seed,lo,hi,8,123)
    assert len(first)==9 and first[0][0]=='historical_seed'
    assert np.array_equal(first[0][1],seed)
    assert all(np.array_equal(a[1],b[1]) for a,b in zip(first,again))
    assert all(np.all(np.abs(p)<=1) for _,p in first)
    with pytest.raises(ValueError,match='power of two'):
        mod.initial_trials(seed,lo,hi,7,123)


def test_only_mass_safe_results_used_for_local_starts():
    par={k:0. for k in mod.bulk.NAMES}
    rows=[{'status':'Rejected','chi2':None},
          {'status':'NeutrinoEvaluated','chi2':5.,'parameters':par},
          {'status':'NeutrinoEvaluated','chi2':2.,'parameters':par},
          {'status':'NeutrinoEvaluated','chi2':float('inf'),'parameters':par}]
    starts=mod.rank_valid_trials(rows,2)
    assert len(starts)==2 and [r['chi2'] for r in starts]==[2.,5.]


def test_optimisation_runs_even_when_historical_seed_already_passes(monkeypatch):
    class FakeEvaluator:
        def __init__(self):
            self.best=None
            self.count=0
            self.mass_safe=0
            self.failures=0
            self.problems=[]
        def evaluate(self,x,label):
            self.count+=1
            self.mass_safe+=1
            residual=np.asarray(x[:5])-.03
            row={'label':label,'status':'NeutrinoEvaluated',
                 'chi2':float(residual@residual),
                 'whitened_residual':residual.tolist(),
                 'parameters':dict(zip(mod.bulk.NAMES,map(float,x)))}
            if self.best is None or row['chi2']<self.best['chi2']:
                self.best=row
            return row
    monkeypatch.setattr(mod,'local_bounds',
           lambda seed,fraction:(np.full(18,-1.),np.full(18,1.)))
    ev=FakeEvaluator()
    hist={'seed':np.zeros(18)}
    result=mod.optimise(ev,hist,sobol_points=2,starts=1,
                         max_nfev=8,seed=5,local_fraction=1.)
    assert result['status']=='Accepted'
    assert result['local_optimisation_attempted']
    assert len(result['optimiser_starts'])==1
    assert result['counts']['evaluations']>3


def test_no_valid_seed_does_not_run_local_solver(monkeypatch):
    class Reject:
        def __init__(self):
            self.best=None;self.count=0;self.mass_safe=0
            self.failures=0;self.problems=[]
        def evaluate(self,x,label):
            self.count+=1
            return {'status':'Rejected','chi2':None,'reason':'ScalarMassWindow'}
    monkeypatch.setattr(mod,'local_bounds',
           lambda seed,fraction:(np.full(18,-1.),np.full(18,1.)))
    got=mod.optimise(Reject(),{'seed':np.zeros(18)},
                     sobol_points=4,starts=2,max_nfev=8,
                     seed=7,local_fraction=1.)
    assert got['status']=='NoMassSafeNeutrinoFit'
    assert not got['local_optimisation_attempted']
    assert got['counts']['evaluations']==5


def test_statuses_selection_and_no_overwrite(tmp_path):
    assert mod.fit_status(None)=='NoMassSafeNeutrinoFit'
    assert mod.fit_status({'chi2':6})=='PoorFit'
    assert mod.fit_status({'chi2':5.86})=='Accepted'
    assert mod.choose_csv_option('A_-2,B_-1',('A_-2','B_-1','C_-1'))==('A_-2','B_-1')
    with pytest.raises(ValueError,match='duplicate'):
        mod.choose_csv_option('B_-1,B_-1',('B_-1',))
    path=tmp_path/'new'/'fit.json'
    mod.write_json(path,{'status':'Accepted'})
    with pytest.raises(FileExistsError):
        mod.write_json(path,{'status':'Accepted'})


def test_saved_fit_config_explicitly_zeroes_reused_quartics(monkeypatch,tmp_path):
    names=mod.bulk.NAMES
    rec={'model':'B_-1','case':'small_uv_frozen_5e-9',
         'status':'Accepted','extra_zero_boundary_scale_gev':1e5,
         'new_additional_quartic_symbols':['lambdaNew'],
         'reused_additional_quartic_symbols_zeroed':['lambdaReuse'],
         'best_mass_safe':{'chi2':0.2,
             'parameters':dict(zip(names,np.zeros(18)))}}
    monkeypatch.setattr(mod.bulk,'REGISTRY',tmp_path/'registry.json',raising=False)
    monkeypatch.setattr(mod.bulk,'read',lambda p:{'models':{'B_-1':{}}},raising=False)
    monkeypatch.setattr(mod.bulk,'historical',
           lambda case,model,entry:{'cfg':{}},raising=False)
    monkeypatch.setattr(mod.bulk,'update_yukawas',
           lambda history,x:{'base_state':{'ordinary':{
               'lambdaReuse':{'real':.5,'imag':.1}}}},raising=False)
    mod.save_case(tmp_path,'B_-1','small_uv_frozen_5e-9',rec)
    import json
    root=tmp_path/'small_uv_frozen_5e-9'/'B_m1'
    result=json.loads((root/'fitted_config.json').read_text())
    assert result['base_state']['ordinary']['lambdaReuse']=={'real':0.,'imag':0.}
    assert result['reduced_scalar_zero_extra_boundary']['boundary_scale_gev']==1e5
    assert (root/'fit.json').is_file()


def test_numerical_matching_failure_is_not_called_mass_screen_failure(monkeypatch):
    class MatchingFailure:
        def __init__(self):
            self.best=None;self.count=0;self.mass_safe=0
            self.screen_pass=0;self.failures=0;self.problems=[]
        def evaluate(self,x,label):
            self.count+=1;self.screen_pass+=1;self.failures+=1
            return {'status':'NumericalFailure','chi2':None,'reason':'C5 mapping unavailable'}
    monkeypatch.setattr(mod,'local_bounds',
           lambda seed,fraction:(np.full(18,-1.),np.full(18,1.)))
    got=mod.optimise(MatchingFailure(),{'seed':np.zeros(18)},
                    sobol_points=2,starts=1,max_nfev=2,seed=19,
                    local_fraction=1.)
    assert got['status']=='MatchingOrEvaluationFailure'
    assert not got['local_optimisation_attempted']


def test_reference_expands_box_without_clipping(monkeypatch):
    seed=np.zeros(18)
    a2=np.zeros(18);a2[0]=2.75;a2[1]=-3.01
    monkeypatch.setattr(mod,'local_bounds',lambda seed,fraction:
        (np.full(18,-.5),np.full(18,.5)))
    lo,hi=mod.extended_bounds(seed,1.,(('A_-2_reference_fit',a2),))
    assert lo[1]<a2[1] and hi[0]>a2[0]
    assert np.all(lo>=-np.sqrt(4*np.pi))
    assert np.all(hi<=np.sqrt(4*np.pi))
    assert np.all(hi>lo)


def test_previous_accepted_fit_used_as_warm_start(monkeypatch,tmp_path):
    import json
    monkeypatch.setattr(mod,'ROOT',tmp_path)
    monkeypatch.setattr(mod,'OUT',tmp_path/'output/studies/reduced_scalar_sobol')
    monkeypatch.setattr(mod.bulk,'read',lambda path:json.loads(path.read_text()),raising=False)
    monkeypatch.setattr(mod.bulk,'sha',lambda path:__import__('hashlib').sha256(path.read_bytes()).hexdigest(),raising=False)
    fit_path=mod.OUT/'old_v1'/'small_uv_frozen_5e-9'/'A_m2'/'fit.json'
    fit_path.parent.mkdir(parents=True)
    ref=np.linspace(-.5,.5,18)
    fit={'model':'A_-2','case':'small_uv_frozen_5e-9',
         'status':'Accepted','best_mass_safe':{
         'chi2':.001,'parameters':dict(zip(mod.bulk.NAMES,ref.tolist()))}}
    fit_path.write_text(json.dumps(fit))
    y,meta,path=mod.previous_best_yukawas('small_uv_frozen_5e-9','A_-2','old_v1')
    assert np.array_equal(y,ref)
    assert meta['source']=='accepted_previous_reduced_fit'
    assert meta['previous_chi2']==.001
    assert path==fit_path
    fit['status']='PoorFit'
    fit_path.write_text(json.dumps(fit))
    with pytest.raises(ValueError,match='not an accepted'):
        mod.previous_best_yukawas('small_uv_frozen_5e-9','A_-2','old_v1')


def test_shared_a2_start_is_part_of_independent_search(monkeypatch):
    class FakeEvaluator:
        def __init__(self):
            self.count=0;self.mass_safe=0;self.failures=0
            self.screen_pass=0;self.best=None;self.problems=[]
        def evaluate(self,x,label):
            self.count+=1;self.mass_safe+=1;self.screen_pass+=1
            chi=float((float(x[0])-.8)**2)
            result={'label':label,'status':'NeutrinoEvaluated','chi2':chi,
                    'parameters':dict(zip(mod.bulk.NAMES,map(float,x))),
                    'whitened_residual':[float(x[0])-.8,0,0,0,0]}
            if self.best is None or chi<self.best['chi2']:
                self.best=result
            return result
    monkeypatch.setattr(mod,'local_bounds',lambda seed,fraction:
        (np.full(18,-.5),np.full(18,.5)))
    # Avoid expensive SciPy iterations; assert the explicit reference point
    # participates in the local solver choice, not only the shared evaluation.
    chosen=[]
    def fake_least_squares(fun,x0,**kwargs):
        chosen.append(np.array(x0))
        return SimpleNamespace(success=True,message='test',nfev=1,cost=0.)
    monkeypatch.setattr(mod,'least_squares',fake_least_squares)
    ev=FakeEvaluator()
    a2=np.zeros(18);a2[0]=.8
    got=mod.optimise(ev,{'seed':np.zeros(18)},sobol_points=2,starts=1,
                     max_nfev=5,seed=17,local_fraction=1.,
                     extra_trials=(('A_-2_reference_fit',a2),))
    assert got['best_mass_safe']['chi2']==0.
    assert got['best_mass_safe']['label']=='A_-2_reference_fit'
    assert got['optimiser_starts'][0]['start']=='A_-2_reference_fit'
    assert got['extra_trial_labels']==['A_-2_reference_fit']
    assert len(chosen)==1
