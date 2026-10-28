"""Offline safety/unit tests. Full Windows numerical integration is user-run."""
from __future__ import annotations

import ast
import copy
from pathlib import Path
import sys
from types import ModuleType,SimpleNamespace

import numpy as np
import pytest

from scripts import RunT3GeneralScalarBCExploratoryFit as fit


def test_parameter_names():
    assert len(fit.NAMES)==18 and len(set(fit.NAMES))==18
    assert fit.NAMES[0]=='y1_00' and fit.NAMES[-1]=='y2_22'


@pytest.mark.parametrize('case',tuple(fit.b9.CASES))
def test_destination_rejects_existing(case,tmp_path):
    dest=fit.fit_destination(case,'new_case',root=tmp_path)
    assert 'general_scalar_extended' in dest.parts
    dest.mkdir(parents=True)
    with pytest.raises(ValueError,match='already exists'):
        fit.fit_destination(case,'new_case',root=tmp_path)


@pytest.mark.parametrize('id',['../old','a/b','bad id',''])
def test_destination_rejects_traversal(id,tmp_path):
    with pytest.raises(ValueError):
        fit.fit_destination('small_uv_running_5e-9',id,root=tmp_path)


def test_reject_unknown_case(tmp_path):
    with pytest.raises(ValueError):
        fit.fit_destination('historical', 'okay',root=tmp_path)


def fake_base(monkeypatch):
    mod=ModuleType('Numerical.orchestration.RunSmallLambdaRefit')
    mod.matr=lambda d: np.asarray(d['real'])+1j*np.asarray(d['imag'])
    mod.setmat=lambda d,z: d.update(real=z.real.tolist(),imag=z.imag.tolist())
    orch=ModuleType('Numerical.orchestration')
    orch.RunSmallLambdaRefit=mod
    monkeypatch.setitem(sys.modules,'Numerical',ModuleType('Numerical'))
    monkeypatch.setitem(sys.modules,'Numerical.orchestration',orch)
    monkeypatch.setitem(sys.modules,'Numerical.orchestration.RunSmallLambdaRefit',mod)


def sample_source(v=0.01):
    raw={'real':np.full((3,3),v).tolist(),'imag':np.zeros((3,3)).tolist()}
    return {'cfg':{'base_state':{'ordinary':{'y1':copy.deepcopy(raw),
                                              'y2':copy.deepcopy(raw)}}}}


def test_trial_record_no_mutation(monkeypatch):
    fake_base(monkeypatch)
    src=sample_source()
    result=fit.trial_record(src,np.arange(18)*0.001)
    assert result['cfg']['base_state']['ordinary']['y1']['real'][0][1]==0.001
    assert result['cfg']['base_state']['ordinary']['y2']['real'][0][0]==pytest.approx(0.009)
    assert src['cfg']['base_state']['ordinary']['y1']['real'][0][1]==0.01


def test_trial_record_rejects_invalid(monkeypatch):
    fake_base(monkeypatch)
    with pytest.raises(ValueError):fit.trial_record(sample_source(),np.ones(17))
    with pytest.raises(ValueError):fit.trial_record(sample_source(),np.r_[np.ones(17),np.nan])


def test_seed_roundtrip(monkeypatch):
    fake_base(monkeypatch)
    src=sample_source(0.25)
    values=fit.seed_vector(src)
    assert values.shape==(18,) and np.all(values==0.25)
    changed=fit.trial_record(src,values)
    assert np.array_equal(fit.seed_vector(changed),values)


def test_seed_rejects_complex(monkeypatch):
    fake_base(monkeypatch)
    src=sample_source()
    src['cfg']['base_state']['ordinary']['y2']['imag'][1][1]=0.5
    with pytest.raises(ValueError,match='18 real Yukawa'):
        fit.seed_vector(src)


def test_real_finite():
    assert fit.real_finite('0.01','x')==0.01
    with pytest.raises(ValueError):fit.real_finite(float('inf'),'x')


def test_manifest_output_status():
    obj={'case':'small_uv_frozen_5e-9','model':'B_-1','run_id':'x',
         'status':'ExploratoryProbe','historical_fit_chi2':4.,
         'chi2_initial':25.,'chi2_final':20.,'passes_approx_chi2_5p86':False,
         'prediction':{'delta_m21_sq_ev2':7.5e-5},
         'candidate_evaluation_count':2,'sobol_failure_count':0}
    x=fit.summary_markdown(obj)
    assert 'uncertified' in x.lower()
    assert 'Historical chi2 = 4' in x
    assert 'general-scalar' in x


def test_scoring_adapter_calls_real_weinberg(monkeypatch):
    calls=[]
    def build_boundary(st,*,matching_scale_gev):
        calls.append(('threshold',matching_scale_gev));return 'boundary'
    def initial(boundary,c5):
        calls.append(('join',boundary));return 'start'
    def evolve(init,first,last,**kwargs):
        calls.append(('run',first,last));return SimpleNamespace(K=np.eye(3)*1e-15,
                                                      ye=np.eye(3))
    def mass(K,*,vev_gev):return K*vev_gev
    def basis(x,ye):return x
    def evaluate(matrix,target):return SimpleNamespace(chi2=1.1,ordering='NO',
      observable_names=('delta_m21_sq_ev2',),prediction=np.array([7.5e-5]),
      whitened_residual=np.array([0.]))
    root=ModuleType('Numerical')
    run=ModuleType('Numerical.running')
    threshold=ModuleType('Numerical.running.ScalarThresholdBoundary')
    threshold.build_final_sm_boundary=build_boundary
    threshold.build_sm_weinberg_initial_conditions=initial
    sm=ModuleType('Numerical.running.SMWeinbergEvolution');sm.evolve_sm_weinberg=evolve
    f=ModuleType('Numerical.fitting');ps=ModuleType('Numerical.fitting.ParameterScan')
    ps.neutrino_mass_from_c5=mass;ps.charged_lepton_mass_basis_matrix=basis
    obs=ModuleType('Numerical.fitting.OscillationFit');obs.evaluate_mass_matrix_fit=evaluate
    for n,mod in [('Numerical',root),('Numerical.running',run),
                  ('Numerical.running.ScalarThresholdBoundary',threshold),
                  ('Numerical.running.SMWeinbergEvolution',sm),
                  ('Numerical.fitting',f),('Numerical.fitting.ParameterScan',ps),
                  ('Numerical.fitting.OscillationFit',obs)]:
        monkeypatch.setitem(sys.modules,n,mod)
    result,low=fit.low_energy_fit(np.eye(3)*1e-15,'eft','target')
    assert result.chi2==1.1
    assert calls==[('threshold',1000.),('join','boundary'),('run',1000.,100.)]


def test_signed_map_and_complex_mass_basis_execution(monkeypatch):
    import mpmath as mp
    base=ModuleType('Numerical')
    experimental=ModuleType('Numerical.experimental')
    mb=ModuleType('Numerical.experimental.GeneralT3MixedScalarMassBasis')
    class Basis:
        def __init__(self,m1,m2,mu):
            matrix=np.array([[m1,mu],[np.conj(mu),m2]])
            self.eigenvalues_sq,self.U=np.linalg.eigh(matrix)
            self.tachyonic=False;self.degenerate=False
        @classmethod
        def build(cls,model,m1,m2,mu):return cls(m1,m2,mu)
    mb.MixedDoubletMassBasis=Basis
    for key,value in [('Numerical',base),('Numerical.experimental',experimental),
                      ('Numerical.experimental.GeneralT3MixedScalarMassBasis',mb)]:
        monkeypatch.setitem(sys.modules,key,value)
    # Exact evaluation is symmetric by ordered+A.T, and must reflect the signed
    # term (-2)*lambda_R in the negative-map reference.
    eft={'mu_gev':1000.,'mS1Sq_gev2':2e6,'mS2Sq_gev2':1e6,
         'quadratic_mu2_gev2':{'real':100.,'imag':50.},
         'lambdaT3':{'real':0.5,'imag':0.},
         'reused_quartics':{'lambdaHHdagS2S2':{'real':0.01,'imag':0.},
                            'lambdaHHdagS1barS1bar':{'real':0.02,'imag':0.}}}
    heavy=(np.array([1e5,1.1e5,1.2e5]),np.eye(3),np.eye(3),0.)
    q1=fit.b7.DIRECT['lambdaHHdagS2S2']
    q2=fit.b7.DIRECT['lambdaHHdagS1barS1bar']
    compiled=compile(f'({q1}+{q2}+lambdaT3)*y1a*y2b/MF+1e-7*T3MEU11','mock','eval')
    c=fit.c5_at_endpoint('B_-1',eft,heavy,compiled)
    assert c.shape==(3,3) and np.allclose(c,c.T)
    assert np.isfinite(c).all()
    # With no y mixing, the two direct quartics are subtracted in the diagonal.
    assert np.isclose(c[0,0].real,
          2*((0.5-0.01-0.02)/1e5+1e-7*Basis(2e6,1e6,100.+50.j).U[0,0].real))


def test_source_no_original_output_writes():
    tree=ast.parse(Path(fit.__file__).read_text(encoding='utf8'))
    text=ast.unparse(tree)
    assert 'dest.mkdir' in text
    assert 'b7.save_json(dest /' in text
    assert "historical_outputs_overwritten': False" in text
    assert 'refine_sobol_candidates' in text
    assert 'heavy_flavour_at_muF(uv_state)' in text
