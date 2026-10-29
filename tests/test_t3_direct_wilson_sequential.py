"""Contracts for the direct-only fermion-first Wilson continuation."""
from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import sympy as sp

from scripts import RunT3DirectWilsonSequential as seq


def test_model_specific_sources_cover_real_representation():
    hard, rge=seq.model_sources('E_-2')
    assert hard.as_posix().endswith('T3_E_alpha_m2/data/final_weinberg_coefficient.json')
    assert rge.as_posix().endswith('T3_E_alpha_m2/data/eft1_wilson_rge.json')


def test_pair_matrix_preserves_majorana_complex_orientation():
    y1=np.zeros((3,3),dtype=complex)
    y2=np.zeros((3,3),dtype=complex)
    y1[0,0]=1+2j
    y2[1,0]=3-4j
    m=seq.pair_matrix(np.array([5,0,0]),y1,y2)
    target=5*np.conj(y1[0,0])*np.conj(y2[1,0])
    assert m[0,1]==pytest.approx(target)
    assert m[1,0]==pytest.approx(target)
    assert m[0,0]==0


def test_pair_matrix_rejects_invalid_layout():
    with pytest.raises(ValueError,match='invalid heavy-flavour'):
        seq.pair_matrix(np.ones(2),np.ones((3,3)),np.ones((3,3)))


def test_pair_matrix_rejects_nonfinite_weights():
    with pytest.raises(ValueError,match='nonfinite'):
        seq.pair_matrix(np.array([1,np.nan,1]),np.ones((3,3)),np.ones((3,3)))


def test_direct_kernels_have_frozen_parity_by_construction():
    hard=lambda mf,lam,h:2*h*lam/mf
    source=lambda mf,lam,h:3*h*lam/mf
    obj=seq.DirectWilsonSource(hard,source,'hard','3','running')
    y1=np.eye(3,dtype=complex)
    y2=2*np.eye(3,dtype=complex)
    mf=np.array([200_000.,230_000.,260_000.])
    lam=1.2+0.35j
    full=obj.matrix(mf,y1,y2,lam)
    t=np.log(1000/100000)
    full+=t*obj.matrix(mf,y1,y2,lam,derivative=True)
    target=np.diag(8*seq.HBAR*lam/mf + 12*seq.HBAR*lam/mf*t)
    assert np.allclose(full,target,rtol=1e-12,atol=1e-28)


def test_log_derivative_and_source_beta_checked(monkeypatch,tmp_path):
    from RGE.matching import FinalWeinbergAdapter as adapter
    mf,ms,lam,hbar=sp.symbols('MF MS lambdaT3 hbar')
    hard=2*hbar*lam/mf
    running=3*hbar*lam*sp.log(ms/mf)/mf
    monkeypatch.setattr(adapter,'load_hierarchical_majorana_c5',
        lambda *a,**k:{'hard_kernel':hard,'running_kernel':running})
    a=tmp_path/'matching.json'; b=tmp_path/'rge.json'
    a.write_text(json.dumps({'status':'Success','scheme':'fixed_order_one_loop_MSbar_hierarchical',
        'combined':{'ready_for_physical_majorana_numerics':True}}))
    def rge(beta):
        b.write_text(json.dumps({'status':'Success','loop_order':1,
            'fixed_order_input':'tree-level EFT1 Wilson coefficients only',
            'weinberg_subspace_validation':{'matches_weinberg_subspace':True,
             'failure_count':0,'beta_kappa_16pi2':beta}}))
    rge('3*lambdaT3*conjugate(y1)*conjugate(y2)/MF')
    compiled=seq.compile_direct_source(a,b)
    assert complex(compiled.derivative(10,2j,1))==pytest.approx(0.6j)
    assert complex(compiled.hard(10,2j,1))==pytest.approx(0.4j)
    rge('4*lambdaT3*conjugate(y1)*conjugate(y2)/MF')
    with pytest.raises(ValueError,match='differs from hierarchical'):
        seq.compile_direct_source(a,b)


def test_nonlog_running_rejected(monkeypatch,tmp_path):
    from RGE.matching import FinalWeinbergAdapter as adapter
    mf,ms,lam,hbar=sp.symbols('MF MS lambdaT3 hbar')
    monkeypatch.setattr(adapter,'load_hierarchical_majorana_c5',
        lambda *a,**k:{'hard_kernel':2*hbar*lam/mf,
                       'running_kernel':hbar*lam*sp.log(ms/mf)**2/mf})
    a=tmp_path/'match.json';b=tmp_path/'beta.json'
    a.write_text(json.dumps({'status':'Success','scheme':'fixed_order_one_loop_MSbar_hierarchical',
        'combined':{'ready_for_physical_majorana_numerics':True}}))
    b.write_text(json.dumps({'status':'Success','loop_order':1,
        'fixed_order_input':'tree-level EFT1 Wilson coefficients only',
        'weinberg_subspace_validation':{'matches_weinberg_subspace':True,
         'failure_count':0,'beta_kappa_16pi2':'lambdaT3*conjugate(y1)*conjugate(y2)/MF'}}))
    with pytest.raises(ValueError,match='unhandled scalar|not a pure fixed-order log'):
        seq.compile_direct_source(a,b)


def test_resume_without_matching_source_is_explicit_control(monkeypatch,tmp_path):
    source_root=tmp_path/'continuous';target_root=tmp_path/'out'
    location=source_root/'A_m4'
    location.mkdir(parents=True)
    (location/'continuous_trajectory.json').write_text('{}')
    (location/'eft.json').write_text('{}')
    monkeypatch.setattr(seq.neutrino,'verify_continuous',lambda model,path:{
        'saved_candidate_eft1_endpoint':'eft.json'})
    monkeypatch.setattr(seq,'model_sources',lambda model:(tmp_path/'missing-match.json',tmp_path/'missing-rge.json'))
    monkeypatch.setattr(seq,'sha',lambda path:'digest')
    report=seq.calculate('A_-4',source_root=source_root,output_root=target_root)
    assert report['status']==seq.SKIP
    assert report['full_sequential_wilson_matching_performed'] is False
    assert 'takagi_masses_ev_descending' not in report
    resumed=seq.calculate('A_-4',source_root=source_root,output_root=target_root,resume=True)
    assert resumed==report


@pytest.mark.parametrize('model',['A_+0','C_-1','E_-2'])
def test_actual_repository_direct_source_if_present(model):
    """Exact Matchete/RGE source compatibility (runs in user's real checkout)."""
    match,beta=seq.model_sources(model)
    if not match.is_file() or not beta.is_file():
        pytest.skip('repository source artifacts are absent in isolated harness')
    result=seq.compile_direct_source(match,beta)
    assert result.beta_text and result.hard_text and result.running_text


def test_historical_frozen_parity_fails_on_wrong_normalisation(monkeypatch):
    import sys
    from dataclasses import dataclass
    from types import ModuleType
    @dataclass(frozen=True)
    class DummyEFT:
        lambdaT3: complex=2+0j
        def validated(self):
            return self
    fake=ModuleType('scripts.FitT3BPlus1GeneralThreeStudies')
    masses=np.array([2e5,2.1e5,2.2e5])
    y1=np.eye(3,dtype=complex)
    y2=2*np.eye(3,dtype=complex)
    direct=seq.DirectWilsonSource(lambda m,lam,h:2*h*lam/m,
                                  lambda m,lam,h:3*h*lam/m,'hard','3','run')
    lam=2+0j
    expected=np.diag(np.array([4*(2*seq.HBAR*lam/m)+4*(3*seq.HBAR*lam/m)*np.log(1e3/m) for m in masses]))
    class MockKernel:
        def __init__(self,path): pass
        def evaluate(self,uv,eft): return expected
    fake.MatchedKernel=MockKernel
    monkeypatch.setitem(sys.modules,fake.__name__,fake)
    assert seq.historical_frozen_parity(direct,Path('unused'),None,DummyEFT(),masses,y1,y2,lam)<1e-15
    fake.MatchedKernel=lambda path:SimpleNamespace(evaluate=lambda uv,eft:expected*2)
    with pytest.raises(ValueError,match='parity failed'):
        seq.historical_frozen_parity(direct,Path('unused'),None,DummyEFT(),masses,y1,y2,lam)


def test_summary_separates_supported_from_legacy(monkeypatch,tmp_path):
    reports={
       'A_+0':{'format':seq.FORMAT,'model':'A_+0','status':seq.STATUS,
                'full_sequential_wilson_matching_performed':False,
                'takagi_masses_ev_descending':[.05,.01,.001]},
       'A_-4':{'format':seq.FORMAT,'model':'A_-4','status':seq.SKIP,
                'full_sequential_wilson_matching_performed':False},
    }
    for model,report in reports.items():
        path=tmp_path/f'{seq.continuous.model_tag(model)}_direct_wilson.json'
        path.write_text(json.dumps(report))
    summary=seq.summarize(tmp_path,tmp_path,tuple(reports))
    assert summary['model_count']==2
    assert summary['direct_source_models']==1
    assert summary['missing_hierarchical_source_models']==1
    assert summary['full_sequential_matching_performed'] is False
    assert (tmp_path/'summary.csv').is_file()


def test_grouped_threshold_logs_reproduce_each_heavy_mass_scale():
    source=seq.DirectWilsonSource(lambda m,l,h:2*h*l/m,
        lambda m,l,h:3*h*l/m,'hard','3','run')
    mf=np.array([180e3,260e3,510e3])
    y1=np.eye(3,dtype=complex)
    y2=np.diag([1.3,2.1,3.7]).astype(complex)
    lam=0.7-0.2j
    hard,relocation,grouped=seq.grouped_hard_boundary(
        source,mf,y1,y2,lam,mu_f=100e3)
    assert np.linalg.norm(relocation)>0
    frozen=grouped+np.log(1e3/100e3)*source.matrix(mf,y1,y2,lam,derivative=True)
    expected=np.diag([2*y2[i,i]*(2*seq.HBAR*lam/mf[i]
                  +3*seq.HBAR*lam/mf[i]*np.log(1e3/mf[i]))
                  for i in range(3)])
    assert np.allclose(frozen,expected,rtol=1e-12,atol=1e-25)


def test_grouped_threshold_log_omission_would_fail_parity():
    source=seq.DirectWilsonSource(lambda m,l,h:2*h*l/m,
        lambda m,l,h:3*h*l/m,'hard','3','run')
    mf=np.array([250e3,300e3,350e3])
    a=b=np.eye(3,dtype=complex)
    hard,relocation,grouped=seq.grouped_hard_boundary(source,mf,a,b,1,mu_f=1e5)
    assert np.linalg.norm(relocation)>0
    naive=hard+np.log(1e3/1e5)*source.matrix(mf,a,b,1,derivative=True)
    fixed=grouped+np.log(1e3/1e5)*source.matrix(mf,a,b,1,derivative=True)
    assert np.linalg.norm(fixed-naive)>1e-9*np.linalg.norm(fixed)
