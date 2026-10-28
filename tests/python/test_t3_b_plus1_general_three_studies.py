"""Offline guards for Batch21. Full native RGBeta/Stage36 parity is live-only."""
from __future__ import annotations
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace
import math
import sys
import types
import numpy as np
import pytest

from scripts import FitT3BPlus1GeneralThreeStudies as f


def good_options(**kw):
    z=dict(sobol_points=8,local_nfev=25,seed_count=1,local_fraction=.8,
         floor_fraction=.03,abs_yukawa_max=2.,extra_fraction=None,extra_real=-.01,
         rho_max=.5,scalar_min_tev=.5,scalar_max_tev=5.)
    z.update(kw)
    return z


@pytest.mark.parametrize('case,expected,uv',[
    ('small_uv_frozen_5e-9',5e-9,False),
    ('large_uv_frozen_0p5',.5,False),
    ('small_uv_running_5e-9',5e-9,True),
])
def test_cases(case,expected,uv,tmp_path):
    c,fit,run,lam=f.input_files(case,tmp_path)
    assert c.name==f.STEM and fit.name==f.STEM
    assert '/fitted/' in c.as_posix() and '/fits/' in fit.as_posix()
    assert run==uv and lam==expected


@pytest.mark.parametrize('invalid', ['../hack','bad/name','x y','', '-abc','x'*65])
def test_run_id_rejects_invalid(invalid):
    with pytest.raises(ValueError):f.checked_run_id(invalid)


@pytest.mark.parametrize('value',['a','run_21','case-bp1','X'*64])
def test_run_id_allowed(value):
    assert f.checked_run_id(value)==value


@pytest.mark.parametrize('key,val',[
    ('sobol_points',3), ('sobol_points',0),('sobol_points',8192),
    ('local_nfev',-1),('local_nfev',1001),('seed_count',0),('seed_count',9),
    ('local_fraction',0),('local_fraction',3.5),('floor_fraction',0),
    ('abs_yukawa_max',0),('abs_yukawa_max',5),('extra_fraction',0),
    ('extra_fraction',1),('extra_real',math.inf),('extra_real',1),('rho_max',1),('scalar_min_tev',0),('scalar_max_tev',.1),
])
def test_bad_options(key,val):
    with pytest.raises(ValueError):f.options(**good_options(**{key:val}))


def test_good_options_and_screen():
    f.options(**good_options())
    seed=np.linspace(-.6,.7,18)
    lo,hi=f.search_box(seed,local_fraction=.8,floor_fraction=.03,abs_yukawa_max=2.)
    assert lo.shape==(18,) and hi.shape==(18,)
    assert np.all(lo<seed) and np.all(seed<hi)
    assert np.all(lo>=-2) and np.all(hi<=2)
    with pytest.raises(ValueError):f.search_box(np.full(18,4.),local_fraction=.8,
                                        floor_fraction=.03,abs_yukawa_max=2.)


@pytest.mark.parametrize('value', [{'real':1}, {'real':1,'imag':math.nan},
                     {'real':1,'imag':0,'bad':0}, None])
def test_bad_complex(value):
    with pytest.raises(ValueError):f.scalar_complex(value)


def test_scalar_complex():
    assert f.scalar_complex({'real':.3,'imag':-.2})==complex(.3,-.2)


def minimal_history(lam=5e-9,uv=False):
    one={'real':np.eye(3).tolist(),'imag':np.zeros((3,3)).tolist()}
    cfg={'base_state':{'ordinary':{'y1':deepcopy(one),'y2':deepcopy(one),
               'lambdaT3':{'real':lam,'imag':0.0}}},'numerical':{'vev_gev':246.22}}
    return {'cfg':cfg,'lambda_fixed':lam,'uv_on':uv,'seed':np.r_[np.eye(3).ravel(),np.eye(3).ravel()]}


@pytest.fixture
def mock_refit(monkeypatch):
    numerical=types.ModuleType('Numerical')
    orchestration=types.ModuleType('Numerical.orchestration')
    refit=types.ModuleType('Numerical.orchestration.RunSmallLambdaRefit')
    refit.matr=lambda d:np.asarray(d['real'])+1j*np.asarray(d['imag'])
    def setmat(d,x):
        d['real']=np.asarray(x.real).tolist()
        d['imag']=np.asarray(x.imag).tolist()
    refit.setmat=setmat
    orchestration.RunSmallLambdaRefit=refit
    monkeypatch.setitem(sys.modules,'Numerical',numerical)
    monkeypatch.setitem(sys.modules,'Numerical.orchestration',orchestration)
    monkeypatch.setitem(sys.modules,'Numerical.orchestration.RunSmallLambdaRefit',refit)
    return refit


def test_candidate_config_preserves_original_lambda_and_input(mock_refit):
    refit=mock_refit
    hist=minimal_history()
    before=deepcopy(hist)
    values=np.arange(18)/100
    c=f.candidate_config(hist,values)
    assert hist['cfg']==before['cfg']
    assert f.scalar_complex(c['base_state']['ordinary']['lambdaT3'])==complex(5e-9)
    assert np.allclose(refit.matr(c['base_state']['ordinary']['y1']).real,
                       values[:9].reshape(3,3))
    assert np.allclose(refit.matr(c['base_state']['ordinary']['y2']).real,
                       values[9:].reshape(3,3))
    assert c['benchmark_search']=={'enabled':False}
    with pytest.raises(ValueError):f.candidate_config(hist,np.ones(17))


def test_chi2_whitening():
    good=SimpleNamespace(whitened_residual=[1.,2.,-3.,0.,0.],chi2=14.)
    assert np.allclose(f.residual_from_fit(good),good.whitened_residual)
    with pytest.raises(ValueError):f.residual_from_fit(SimpleNamespace(whitened_residual=[1]*5,chi2=0))


def test_summary_includes_frozen_status_and_input_lambda():
    report={'case':'large_uv_frozen_0p5','status':'CandidateConstrainedPoorFit',
      'uv_running_per_trial':False,'historical_lambdaT3_fixed':.5,
      'extra_fraction':None,'extra_quartic_boundary_real':-.01,
      'extra_quartic_boundary_rule':'fixed real RGBeta boundary',
      'matching_kernel_stage36_parity_relative_error':0.,
      'counts':{'candidate_evaluations':1,'mass_safe_c5_evaluations':1,
        'numerical_failures':0},'best_mass_safe':{'chi2':7.,'rhoF':.01,
        'scalar_masses_muS_tev':[1.,1.2],
        'max_yukawa_singular_muF':.001,'max_scalar_quartic_abs':.5,
        'prediction':{'sin2_theta12':.3}}}
    md=f.summary(report)
    assert 'frozen UV' in md and '`0.5`' in md and '7' in md
    assert 'not certified' in md


def test_regression_batch20_output_reader(tmp_path,monkeypatch):
    class Faux:
        RESULT_FORMAT='ok'
        @staticmethod
        def validate_export(p,state):
            assert p['format']=='native'
            return ({'lambdaT3GeneralQ022000I1':'1'},'dummy_layout')
    monkeypatch.setattr(f,'b20',Faux)
    base=tmp_path/'output/studies/general_scalar_extended/b_p1_eft1_running/run1'
    base.mkdir(parents=True)
    native=base/'native_eft1_betas.json';native.write_text('{"format":"native"}')
    report={'format':'ok','status':'SuccessCoupledBPlus1GeneralEFT1Running',
            'model':f.MODEL,'new_stage29_quartic_evolved_in_eft1':True,
            'zero_control_vs_original_stage34':{'status':'ZeroExtraAgreesWithSavedStage34'},
            'input_sha256':{r'C:\\thesis\\native_eft1_betas.json':f.sha(native)}}
    (base/'report.json').write_text(__import__('json').dumps(report))
    payload,layout,path,reportpath=f.verify_beta_source(tmp_path,'run1',object())
    assert payload['format']=='native' and layout=='dummy_layout' and path==native
    report['input_sha256'][r'C:\\thesis\\native_eft1_betas.json']='0'*64
    (base/'report.json').write_text(__import__('json').dumps(report))
    with pytest.raises(ValueError,match='digest'):f.verify_beta_source(tmp_path,'run1',object())


def test_cached_fit_row_never_includes_chi2_on_error(mock_refit,monkeypatch):
    monkeypatch.setattr(f,'prepare_uv',lambda *a,**k:(_ for _ in ()).throw(ValueError('source missing')))
    fake=object.__new__(f.Evaluator)
    fake.historical=minimal_history()
    fake.betas={}
    fake.uv_betas=None
    fake.kernel=None
    fake.target=None
    fake.extra_fraction=None
    fake.extra_real=-.01
    fake.rho_max=.5
    fake.scalar_min_tev=.5
    fake.scalar_max_tev=5.
    fake.root=Path('/tmp/no-files')
    fake.extra_boundary=complex(-.01)
    fake.cache={};fake.best=None
    fake.calls=0;fake.failed=0;fake.mass_safe=0
    fake.frozen=None;fake.parity=None;fake.total_eft_nfev=0
    x=np.zeros(18)
    row=fake.evaluate(x,label='test')
    assert row['status']=='NumericalTrialFailed' and row['chi2'] is None
    assert fake.evaluate(x,label='same') is row and fake.calls==1
