"""Offline unit tests: no RGEs or fitting are invoked."""
from __future__ import annotations

import importlib.util
import json
import sys
import types
from pathlib import Path

import numpy as np
import pytest

SCRIPT = Path(__file__).resolve().parents[2] / 'scripts' / 'BuildReducedScalarExistingInteractiveReports.py'
sp = importlib.util.spec_from_file_location('reduced_interactive_report_under_test', SCRIPT)
assert sp and sp.loader
m = importlib.util.module_from_spec(sp)
sp.loader.exec_module(m)


def test_grid_is_reciprocal_and_centred():
    grid = m.gamma_grid(.1, 10., 11)
    assert len(grid) == 11 and grid[5] == 1.0
    assert np.all(np.diff(grid) > 0)
    assert np.allclose(np.asarray(grid) * np.asarray(grid[::-1]), 1.)


@pytest.mark.parametrize('lo,hi,n',[(1.,10.,11),(.1,1.,11),(.1,10.,10),(.2,10.,11)])
def test_invalid_grid(lo, hi, n):
    with pytest.raises(ValueError):
        m.gamma_grid(lo, hi, n)


def test_t1_preserves_product():
    x = np.arange(1., 19.) / 13.
    y,lam=m.scale_inputs(x,.5,'T1',4.)
    assert np.allclose(y,x/2.) and lam == 2.
    assert np.allclose(lam*y[:9]*y[9:],.5*x[:9]*x[9:])


def test_t2_preserves_product():
    x=np.arange(1.,19.) / 13.
    y,lam=m.scale_inputs(x,.5,'T2',4.)
    assert np.allclose(y[:9],x[:9]) and np.allclose(y[9:],x[9:]/4.) and lam==2.
    assert np.allclose(lam*y[:9]*y[9:],.5*x[:9]*x[9:])


def test_t3_preserves_product():
    x=np.arange(1.,19.) / 13.
    y,lam=m.scale_inputs(x,.5,'T3',4.)
    assert np.allclose(y[:9],x[:9]*4) and np.allclose(y[9:],x[9:]/4.) and lam==.5
    assert np.allclose(lam*y[:9]*y[9:],.5*x[:9]*x[9:])


def test_shared_reference_selected_and_verified():
    a={'parameters':{f'p{i}':float(i) for i in range(18)},
       'chi2':1.,'prediction':{k:1. for k in m.OBS}}
    b={'parameters':{f'p{i}':float(i+25) for i in range(18)},
       'chi2':2.,'prediction':{k:2. for k in m.OBS}}
    ref={'chi2':12.,'prediction':{k:3. for k in m.OBS}}
    ctx={'fit':{'best_mass_safe':a,'shared_a2':ref},
         'donor':{'best_mass_safe':b},'bulk':types.SimpleNamespace(NAMES=tuple(a['parameters']))}
    x,obs=m.chosen_parameters(ctx,'optimised')
    assert x[0]==0 and obs is a
    x,obs=m.chosen_parameters(ctx,'shared_a_m2')
    assert x[0]==25 and obs is ref


def test_gamma_hist_does_not_mutate_original():
    old={'lambda_fixed':.5,'cfg':{'base_state':{'ordinary':{'lambdaT3':{'real':.5,'imag':0.}}}}}
    new=m.gamma_hist(old,2.)
    assert new['lambda_fixed']==2 and new['cfg']['base_state']['ordinary']['lambdaT3']['real']==2
    assert old['lambda_fixed']==.5 and old['cfg']['base_state']['ordinary']['lambdaT3']['real']==.5


def test_check_endpoint_rejects_different_fit():
    original={'chi2':1.,'prediction':{k:float(i) for i,k in enumerate(m.OBS)}}
    good={'status':'NeutrinoEvaluated','chi2':1.,'prediction':original['prediction']}
    m.check_endpoint(good,original,where='test')
    bad=dict(good, chi2=2.)
    with pytest.raises(ValueError):m.check_endpoint(bad,original,where='test')


def test_reuse_compatible_requires_fit_hash(tmp_path,monkeypatch):
    f=tmp_path/'fit.json';d=tmp_path/'donor.json'
    f.write_text('{}');d.write_text('{}')
    ctx={'fit_path':f,'donor_path':d,'case':m.CASES[0],'model':'A_-2'}
    p=tmp_path/'diagnostic.json'
    obj={'format':m.FORMAT,'status':'Success','mode':'optimised',
         'model':'A_-2','case':m.CASES[0],
         'source_fit_sha256':m.sha(f),'source_A_m2_fit_sha256':m.sha(d),
         'final_running':{'mu_gev':[1,2,3]}}
    m.write_json(p,obj)
    assert m.reuse_compatible(p,ctx,'optimised',n_points=3)
    f.write_text('CHANGED')
    assert not m.reuse_compatible(p,ctx,'optimised',n_points=3)


def test_report_calls_original_renderers(tmp_path, monkeypatch):
    calls=[]
    pkg=types.ModuleType('Numerical');plot=types.ModuleType('Numerical.plotting')
    rg=types.ModuleType('Numerical.plotting.InteractiveModelComparison')
    ev=types.ModuleType('Numerical.plotting.BuildFittedModelEnvelopes')
    rg.write_dashboard=lambda source,dest,title: (calls.append(('rge', source, dest)),dest.parent.mkdir(parents=True,exist_ok=True),dest.write_text('dashboard'))
    ev.build=lambda source,dest: (calls.append(('oval',source,dest)), dest.write_text('ovals'))
    for name,mod in [('Numerical',pkg),('Numerical.plotting',plot),
                     ('Numerical.plotting.InteractiveModelComparison',rg),
                     ('Numerical.plotting.BuildFittedModelEnvelopes',ev)]:
        monkeypatch.setitem(sys.modules,name,mod)
    src=tmp_path/'src'
    model=src/'sensitivity'/'T3_dS1_1_dS2_3_dF_2_alpha_m2'
    model.mkdir(parents=True)
    m.write_json(model/'sensitivity_diagnostics.json',{'n_points':3,'directions':{
        'T1':{'points':[{}, {'chi2':1.2},{}]}}})
    dest=m.render_existing(src,tmp_path/'display','optimised',m.CASES[0],'all')
    assert (dest/'interactive_rge.html').read_text()=='dashboard'
    assert (dest/'sensitivity_ovals.html').read_text()=='ovals'
    assert {call[0] for call in calls}=={'rge','oval'}
    assert m.read(src/'sensitivity'/'comparison_summary.json')[0]['verified_chi2']==1.2


def test_original_renderer_output_does_not_overlap_old_lambda_roots():
    assert 'display_reduced_scalar' in str(m.report_root('abc'))
    assert 'display_fermion_start' not in str(m.report_root('abc'))
    assert 'reduced_scalar_interactive_diagnostics' in str(m.diagnostics_root('abc'))


def test_frozen_eft_constructor_integration_failure_is_a_failed_gamma_point():
    import types
    calls = []

    class FakeReducedEvaluator:
        def __init__(self, case, model, entry, new, reused, hist, beta, kernel, target):
            lam = hist['lambda_fixed']
            calls.append(lam)
            if lam > 1:
                raise ValueError('T3 48-study campaign: Coupled EFT1 integration failed: '
                                 'Required step size is less than spacing between numbers.')

        def evaluate(self, vector, label):
            return {'status': 'NeutrinoEvaluated', 'chi2': 0.0,
                    'prediction': {key: 0.0 for key in m.OBS}}

    ctx = {'hist': {'lambda_fixed': .5,
                    'cfg': {'base_state': {'ordinary': {
                        'lambdaT3': {'real': .5, 'imag': 0.0}}}}},
           'reduced': types.SimpleNamespace(ReducedEvaluator=FakeReducedEvaluator),
           'case': m.CASES[1], 'model': 'A_-2', 'entry': {},
           'new': (), 'reused': (), 'beta': {}, 'kernel': None, 'target': None}
    cache = {}
    failed = m.evaluate_gamma(ctx, np.zeros(18), 10., 'T1', cache)
    assert failed['status'] == 'NumericalFailure'
    assert 'Required step size' in failed['reason']
    assert len(calls) == 1
    # T1 and T2 at the same gamma use the same failed frozen EFT; do not retry.
    assert m.evaluate_gamma(ctx, np.zeros(18), 10., 'T2', cache) == failed
    assert len(calls) == 1
    baseline = m.evaluate_gamma(ctx, np.zeros(18), 1., 'T1', cache)
    assert baseline['status'] == 'NeutrinoEvaluated'


def test_gamma_one_constructor_error_must_fail():
    import types

    def crash(*args):
        raise ValueError('T3 48-study campaign: Coupled EFT1 integration failed: underflow')

    ctx = {'hist': {'lambda_fixed': .5,
                    'cfg': {'base_state': {'ordinary': {
                        'lambdaT3': {'real': .5, 'imag': 0.0}}}}},
           'reduced': types.SimpleNamespace(ReducedEvaluator=crash),
           'case': m.CASES[1], 'model': 'A_-2', 'entry': {},
           'new': (), 'reused': (), 'beta': {}, 'kernel': None, 'target': None}
    with pytest.raises(ValueError, match='Coupled EFT1 integration failed'):
        m.evaluate_gamma(ctx, np.zeros(18), 1., 'T3', {})


def test_unexpected_constructor_error_remains_fatal():
    import types

    def crash(*args):
        raise ValueError('Live RGBeta beta inventory mismatch')

    ctx = {'hist': {'lambda_fixed': .5,
                    'cfg': {'base_state': {'ordinary': {
                        'lambdaT3': {'real': .5, 'imag': 0.0}}}}},
           'reduced': types.SimpleNamespace(ReducedEvaluator=crash),
           'case': m.CASES[1], 'model': 'A_-2', 'entry': {},
           'new': (), 'reused': (), 'beta': {}, 'kernel': None, 'target': None}
    with pytest.raises(ValueError, match='inventory mismatch'):
        m.evaluate_gamma(ctx, np.zeros(18), 10., 'T1', {})


def test_sensitivity_record_preserves_failed_gamma_point(tmp_path, monkeypatch):
    import types
    fit=tmp_path/'fit.json'
    donor=tmp_path/'donor.json'
    fit.write_text('{}')
    donor.write_text('{}')
    pred={key: float(i + 1) for i, key in enumerate(m.OBS)}
    monkeypatch.setattr(m, 'chosen_parameters',
                        lambda ctx, mode: (np.zeros(18), {'chi2':1.,'prediction':pred}))

    def eval_mock(ctx, vector, gamma, direction, cache):
        if gamma == 10. and direction in ('T1','T2'):
            return {'status': 'NumericalFailure', 'reason': 'Coupled EFT1 failed'}
        return {'status': 'NeutrinoEvaluated', 'chi2': 1., 'prediction': pred}

    monkeypatch.setattr(m, 'evaluate_gamma', eval_mock)
    target = types.SimpleNamespace(ordering='NO', observable_names=m.OBS,
         central_values=np.ones(5),covariance=np.eye(5),source='test')
    ctx={'fit_path':fit, 'donor_path':donor, 'case':m.CASES[1],
         'model':'A_-2', 'target':target}
    rec=m.sensitivity_record(ctx,'optimised',[.1,1.,10.])
    assert rec['status']=='Success' and rec['n_failed_points']==2
    assert rec['directions']['T1']['points'][-1]['status']=='NumericalFailure'
    assert rec['directions']['T1']['points'][-1]['message']=='Coupled EFT1 failed'
    assert all(rec['directions'][d]['points'][1]['status']=='NeutrinoEvaluated'
               for d in m.DIRECTIONS)


def test_expected_tachyonic_scalar_threshold_is_a_failed_trial():
    assert m.is_expected_scalar_rge_failure(ValueError(
        'mS1Sq must be positive to define a physical scalar threshold; got -123 GeV^2.'))
    assert not m.is_expected_scalar_rge_failure(ValueError(
        'Live RGBeta beta inventory mismatch: reused_quartic_symbols'))
