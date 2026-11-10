"""Standalone safety/early-stop tests; do NOT claim physical Matchete verification."""
from __future__ import annotations

import sys
import types
from types import SimpleNamespace

import numpy as np
import pytest

from Numerical.parallel_t3.complete import all_model_conditional_optimizer as m


@pytest.fixture(autouse=True)
def fake_residual5(monkeypatch):
    """Use a small explicit five-observable Gaussian mock, never for production."""
    name='Numerical.parallel_t3.complete.two_model_conditional_uvfit_batch12'
    native=types.ModuleType(name)
    native.residual5=lambda obs,target:np.array([
        (obs[k]-target['central_values'][k])/target['one_sigma_errors'][k]
        for k in m.OBS],float)
    monkeypatch.setitem(sys.modules,name,native)


class FakeEngine:
    target={'central_values':{k:0. for k in m.OBS},
            'one_sigma_errors':{k:1. for k in m.OBS}}

    def __init__(self, *, offset=0.):
        self.offset=float(offset)
        self.calls=[]

    def evaluate(self,x,*,full=False):
        v=np.asarray(x,float)
        self.calls.append((v.copy(),full))
        r=v[:5]/1e-4-self.offset
        obs=dict(zip(m.OBS,r.tolist()))
        return m.ConditionalPrediction(obs,float(r@r),np.eye(3)*1e-16,
                           {'fake':'not physics'},{},None)


def test_exactly_sixteen_independent_model_fit_jobs(monkeypatch):
    inventory=tuple(['A_-4','D_-2']+[f'X_{i}' for i in range(14)])
    monkeypatch.setattr(m,'all_models',lambda:inventory)
    jj=m.jobs(inventory)
    assert len(jj)==16
    assert len({x['model'] for x in jj})==16
    assert sum(v['supported'] for v in jj)==2
    assert all(j['backend'] is None and j['reason'] for j in jj if not j['supported'])
    with pytest.raises(m.SourceBlocked,match='repeated'):
        m.jobs(("A_-4","A_-4"))


def test_two_stored_benchmarks_are_initial_guesses_not_two_fits():
    cfg={'fitting':{'input_Yukawa_scale':1e-4},'models':{
        'A_-4':{'y1_UV_times_1e4':(np.ones((3,3))*2).tolist(),
                'y2_UV_times_1e4':(np.ones((3,3))*2).tolist()},
        'D_-2':{'y1_UV_times_1e4':(np.zeros((3,3))).tolist(),
                'y2_UV_times_1e4':(np.zeros((3,3))).tolist()}}}
    engine=FakeEngine()
    label,x,screen=m.choose_start(engine.evaluate,cfg,policy='best-of-two')
    assert label=='D_-2'
    assert screen['A_-4']>screen['D_-2']
    np.testing.assert_array_equal(x,np.zeros(18))
    assert len(engine.calls)==2


def test_initial_seed_at_target_stops_without_scipy(monkeypatch):
    def should_not_optimize(*args,**kwargs):
        raise AssertionError('Do not invoke scipy when initial point meets goal')
    monkeypatch.setattr(m,'least_squares',should_not_optimize)
    engine=FakeEngine()
    x=np.zeros(18)
    out=m.optimize_job(engine.evaluate,x,max_nfev=100,bound=.005,target_chi2=.01)
    assert out['target_met_after_fresh_benchmark_validation'] is True
    assert out['best_chi2']==0.
    assert out['status']=='ConditionalTargetChi2Reached_NotPhysicalC5'
    assert out['optimizer_nfev'] is None
    assert len(engine.calls)==2
    assert engine.calls[-1][1] is True


def test_immediate_stop_on_first_successful_scipy_trial(monkeypatch):
    engine=FakeEngine()
    visited=[]
    def mock_scipy(f,x,**kw):
        visited.append('started')
        f(np.zeros(18))  # target met; must escape optimizer immediately
        visited.append('wrong: continued after target')
        return SimpleNamespace(success=False,nfev=1,message='no')
    monkeypatch.setattr(m,'least_squares',mock_scipy)
    out=m.optimize_job(engine.evaluate,np.ones(18)*1e-4,
                       max_nfev=10,bound=.005,target_chi2=.01)
    assert visited==['started']
    assert out['target_met_after_fresh_benchmark_validation']
    assert out['early_stop_triggered_inside_solver']
    assert out['stop_reason']=='verified_chi2_below_target'
    assert out['best_chi2']==0.


def test_above_goal_not_success_even_on_solver_convergence(monkeypatch):
    engine=FakeEngine()
    monkeypatch.setattr(m,'least_squares',lambda f,x,**kw:
                        SimpleNamespace(success=True,nfev=2,message='xtol'))
    out=m.optimize_job(engine.evaluate,np.ones(18)*1e-4,
                       max_nfev=3,bound=.005,target_chi2=.01)
    assert not out['target_met_after_fresh_benchmark_validation']
    assert out['status']=='ConditionalSolverConverged_AboveTarget_NotPhysicalC5'
    assert out['best_chi2']>0.01


def test_target_is_strictly_less_not_less_or_equal(monkeypatch):
    # Residual one component = 0.1, hence chi2 == 0.01 in exact binary float approx.
    engine=FakeEngine()
    monkeypatch.setattr(m,'least_squares',lambda f,x,**kw:
                        SimpleNamespace(success=False,nfev=1,message='max_nfev'))
    x=np.zeros(18)
    x[0]=1e-5
    target=float((x[0]/1e-4)**2)
    out=m.optimize_job(engine.evaluate,x,max_nfev=1,bound=.005,
                       target_chi2=target)
    assert out['best_chi2']==pytest.approx(target)
    assert not out['target_met_after_fresh_benchmark_validation']


def test_bad_bounds_and_target_are_rejected_before_first_pipeline_call():
    e=FakeEngine()
    with pytest.raises(m.SourceBlocked,match='outside'):
        m.optimize_job(e.evaluate,np.ones(18)*.006,max_nfev=3,bound=.005)
    with pytest.raises(m.SourceBlocked,match='target_chi2'):
        m.optimize_job(e.evaluate,np.ones(18)*1e-4,max_nfev=3,
                       bound=.005,target_chi2=0)
    assert not e.calls


def test_final_benchmark_must_reproduce_best_point(monkeypatch):
    class BrokenFinal(FakeEngine):
        def evaluate(self,x,*,full=False):
            v=super().evaluate(x,full=full)
            if full:
                return m.ConditionalPrediction(v.observables,1234.,v.c5_source_muS,
                                               v.matching,v.rg_nfev,v.result)
            return v
    e=BrokenFinal()
    with pytest.raises(m.SourceBlocked,match='did not reproduce'):
        m.optimize_job(e.evaluate,np.zeros(18),max_nfev=3,bound=.005)
