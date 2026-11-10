"""Deterministic optimiser tests only; NOT an EFT/Matchete integration test."""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
import numpy as np
import pytest

FILE=Path(__file__).resolve().parents[1]/'Numerical/parallel_t3/complete/all16_rough_yukawa_fit.py'
spec=importlib.util.spec_from_file_location('all16_sobol_fit_under_test',FILE)
m=importlib.util.module_from_spec(spec)
sys.modules[spec.name]=m
spec.loader.exec_module(m)


def cfg():
    a=np.arange(18,dtype=float).reshape(2,3,3)*.01
    b=-a-.1
    def pack(x):
        return {'y1_UV_times_1e4':x[0].tolist(),
                'y2_UV_times_1e4':x[1].tolist()}
    return {'fitting':{'input_Yukawa_scale':1e-4},
            'models':{'A_-4':pack(a),'D_-2':pack(b)}}


def constant_prediction(x,chi2=1.0):
    return m.Trial({'mixing_defined':True},[chi2**.5,0,0,0,0],
                   float(chi2),{'test_only':True})


def kwargs(**overrides):
    d={'model':'A_-4','target_chi2':.01,'max_nfev':2,
       'max_trials':50,'yukawa_bound':.005,
       'sobol_points':4,'screen_points':3,
       'active_schedule':(2,),'starts':1,
       'perturbation':1.0,'rng_seed':17}
    d.update(overrides)
    return d


def test_source_exact_two_seeds_are_distinct():
    a=m.seed_pair(cfg(),'A_-4')
    b=m.seed_pair(cfg(),'D_-2')
    assert a.shape==b.shape==(18,)
    assert not np.array_equal(a,b)
    y1,y2=m.physical_yukawas(a)
    assert y1[0,0]==0
    assert y2[-1,-1]==pytest.approx(.17e-4)


def test_search_box_respects_absolute_uv_bounds():
    seed=m.seed_pair(cfg(),'A_-4')
    lo,hi,c,w=m._search_box(seed,bound=.00005,perturbation=20.)
    assert np.max(hi)<=.5+1e-12
    assert np.min(lo)>=-.5-1e-12
    assert np.all(w>0)
    assert np.allclose(c-w,lo)
    assert np.allclose(c+w,hi)


def test_initial_benchmark_reaching_goal_stops_before_sobol(monkeypatch):
    monkeypatch.setattr(m,'all_models',lambda:('A_-4','D_-2'))
    count={'n':0}
    def evaluate(x):
        count['n']+=1
        return constant_prediction(x,chi2=0.)
    result=m.optimize_one(evaluate,cfg(),**kwargs())
    assert result['target_met']
    assert result['best_chi2']==0.
    assert result['unique_physics_evaluations']==1
    assert count['n']==2  # first benchmark + independent pipeline recheck
    assert result['sensitivity_probes']==[]
    assert result['local_least_squares_runs']==[]


def test_strict_goal_equality_is_not_accepted(monkeypatch):
    monkeypatch.setattr(m,'all_models',lambda:('A_-4','D_-2'))
    result=m.optimize_one(lambda x:constant_prediction(x,.01),cfg(),
        **kwargs(sobol_points=0,screen_points=0,active_schedule=(1,),max_nfev=1))
    assert result['verification']['passed']
    assert not result['target_met']


def test_sobol_sensitivity_and_bounded_least_squares_all_execute(monkeypatch):
    monkeypatch.setattr(m,'all_models',lambda:('A_-4','D_-2'))
    def evaluate(x):
        x=np.asarray(x,float)
        r=np.asarray([.8+.1*x[0],.5+.1*x[1],.4,0,0])
        return m.Trial({},r.tolist(),float(r@r),{'test_only':True})
    result=m.optimize_one(evaluate,cfg(),**kwargs())
    phases=[r['phase'] for r in result['trials']]
    assert len([k for k in phases if k.startswith('initial_seed')])==2
    assert phases.count('sobol')==4
    assert phases.count('sensitivity')>=1
    assert any(k.startswith('least_squares') for k in phases)
    assert result['search_method'].startswith('scrambled Sobol')
    assert result['active_parameter_stages']==[[0,1]] or len(result['active_parameter_stages'])==1
    assert not result['physical_C5_certified']
    assert result['unique_physics_evaluations']<=50


def test_sobol_is_seeded_deterministically(monkeypatch):
    monkeypatch.setattr(m,'all_models',lambda:('A_-4','D_-2'))
    run=lambda:m.optimize_one(lambda x:constant_prediction(x),cfg(),
                        **kwargs(screen_points=0,active_schedule=(1,),max_nfev=1))
    x=run()['trials']
    y=run()['trials']
    assert [(a['phase'],a.get('chi2')) for a in x]==[(a['phase'],a.get('chi2')) for a in y]


def test_trial_budget_counts_seeds_and_sobol(monkeypatch):
    monkeypatch.setattr(m,'all_models',lambda:('A_-4','D_-2'))
    result=m.optimize_one(lambda x:constant_prediction(x),cfg(),
                         **kwargs(max_trials=3))
    assert result['unique_physics_evaluations']==3
    assert result['optimizer_message'].startswith('Unique fresh-pipeline')
    assert result['verification']['passed']


def test_both_blocked_seeds_do_not_fabricate_goal(monkeypatch):
    monkeypatch.setattr(m,'all_models',lambda:('A_-4','D_-2'))
    def bad(x):raise RuntimeError('real pipeline input missing')
    result=m.optimize_one(bad,cfg(),**kwargs())
    assert result['status']=='BLOCKED_NO_VALID_SEED'
    assert not result['target_met']
    assert all(row['status']=='BLOCKED' for row in result['trials'])


def test_pipeline_recheck_can_reject_apparent_goal(monkeypatch):
    monkeypatch.setattr(m,'all_models',lambda:('A_-4','D_-2'))
    count={'n':0}
    def evaluator(x):
        count['n']+=1
        return constant_prediction(x,0.001 if count['n']==1 else 10.)
    result=m.optimize_one(evaluator,cfg(),**kwargs())
    assert not result['target_met']
    assert result['status']=='BLOCKED_NONREPRODUCIBLE_BEST'


def test_bad_active_schedule_is_rejected():
    with pytest.raises(m.PipelineBlocked,match='Invalid Sobol'):
        m.optimize_one(lambda x:constant_prediction(x),cfg(),
                       **kwargs(active_schedule=(0,)))


def test_recheck_does_not_count_as_trial(monkeypatch):
    monkeypatch.setattr(m,'all_models',lambda:('A_-4','D_-2'))
    count={'n':0}
    def evaluator(x):
        count['n']+=1
        return constant_prediction(x,0.001)
    out=m.optimize_one(evaluator,cfg(),**kwargs())
    assert count['n']==out['unique_physics_evaluations']+1
