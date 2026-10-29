"""Offline fail-closed tests for Stage 33; live benchmarks run on user machine."""
from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path
from types import SimpleNamespace
import sys
import types

import numpy as np
import pytest

from scripts import FermionThresholdHandoff as stage


@dataclass
class FakeSM:
    gY: float
    g2: float
    g3: float
    lambdaH: float
    ye: np.ndarray


@dataclass
class FakeUV:
    mu_gev: float
    sm: FakeSM
    MF: np.ndarray
    lambdaHHdagS2S2: complex | None = None


def test_16_models_and_input_contract():
    assert len(stage.MODELS) == 16
    assert stage.MODELS[0] == 'A_+0'
    assert len(set(stage.MODELS)) == 16
    assert stage.s30.START == 1e7 and stage.s30.STOP == 1e5


def test_complex_matrices_are_not_silently_projected_to_real():
    obj = FakeUV(1e5, FakeSM(0.37, 0.6, 0.9, 0.275,
                            np.array([[1+2j,3-4j]], dtype=complex)),
                 np.diag([1e5,2e5]).astype(complex), -0.015+0.001j)
    exported = stage.json_state(obj)
    assert exported['sm']['ye'] == {'real': [[1., 3.]], 'imag': [[2., -4.]]}
    assert exported['lambdaHHdagS2S2'] == {'real': -0.015, 'imag': 0.001}
    assert exported['MF']['real'] == [[1e5,0.], [0.,2e5]]
    assert stage.json_state(None) is None


def test_nonfinite_complex_rejected():
    with pytest.raises(ValueError, match='Nonfinite'):
        stage.json_state(complex(float('nan'),0))
    with pytest.raises(ValueError, match='Nonfinite'):
        stage.json_state(np.array([[complex(0,float('inf'))]]))


def test_strict_endpoint_check():
    stage._must_close(.276386, .276386000001, 'lambdaH')
    with pytest.raises(ValueError, match='Stage33 fail-closed'):
        stage._must_close(.276386, .278, 'lambdaH')
    with pytest.raises(ValueError, match='Stage33 fail-closed'):
        stage._must_close(.276386, float('nan'), 'lambdaH')


def test_source_hashes_windows_and_posix_paths(tmp_path, monkeypatch):
    (tmp_path / 'RGE' / 'general').mkdir(parents=True)
    path=tmp_path / 'RGE' / 'general' / 'registry.json'
    path.write_bytes(b'audited')
    other=[]
    for i in range(4):
        p=tmp_path/f'item{i}.txt';p.write_bytes(str(i).encode());other.append(p)
    mapping={
        'RGE\\general\\registry.json':sha256(path.read_bytes()).hexdigest(),
        **{f'item{i}.txt':sha256(p.read_bytes()).hexdigest() for i,p in enumerate(other)},
    }
    monkeypatch.setattr(stage, 'ROOT', tmp_path)
    stage.verify_sources({'input_sha256':mapping})
    path.write_bytes(b'tampered')
    with pytest.raises(ValueError,match='input changed'):
        stage.verify_sources({'input_sha256':mapping})


def test_report_and_inventory_are_fail_closed():
    obj={'model':'A_+0','status':'SuccessOneModelNumericalUVRunning',
         'one_loop_beta_evolution_performed':True,
         'full_coupled_uv_beta_system_evolved':True,
         'physical_coupling_mapping_certified':False,
         'mu_initial_gev':1e7,'mu_final_gev':1e5,
         'source_one_loop_beta_convention':stage.s30.BETA_CONVENTION,
         'new_quartic_couplings_evolved':5,'sidecar_extra_real_dimensions':10}
    stage._verify_report(obj,'A_+0',('n1','n2','n3','n4','n5'),())
    for change in [dict(physical_coupling_mapping_certified=True),
                   dict(full_coupled_uv_beta_system_evolved=False),
                   dict(sidecar_extra_real_dimensions=8)]:
        bad=dict(obj,**change)
        with pytest.raises(ValueError,match='Original UV report'):
            stage._verify_report(bad,'A_+0',('n1','n2','n3','n4','n5'),())


def test_full_compact_trajectory_audit_accepts_and_rejects_new_and_reused(monkeypatch):
    mod=types.ModuleType('Numerical.core.StateVector')
    def unpack(vec, layout, mu_gev):
        return SimpleNamespace(sm=SimpleNamespace(gY=vec[0],g2=.6,g3=.9,lambdaH=.276),
                               lambdaHHdagS2S2=complex(-.014,.002))
    mod.unpack_uv_state=unpack
    monkeypatch.setitem(sys.modules,'Numerical.core.StateVector',mod)
    layout=SimpleNamespace(size=1)
    solved=SimpleNamespace(y=np.array([[.38,.375,.37],
                                       [-.01,-.011,-.012],
                                       [0.,0.,0.]]),
                           t=np.log([1e7,1e6,1e5]),nfev=12)
    rows=[]
    for i,mu in enumerate((1e7,1e6,1e5)):
        rows.append({'mu_gev':mu,'gY':solved.y[0,i],'g2':.6,'g3':.9,'lambdaH':.276,
                     'extra_quartics':{'new':{'real':solved.y[1,i],'imag':0.}},
                     'reused_quartics':{'lambdaHHdagS2S2':{'real':-.014,'imag':.002}}})
    report={'trajectory':rows,'n_saved_points':3,'nfev':12}
    stage.verify_saved_trajectory(report,solved,layout,('new',),('lambdaHHdagS2S2',))
    broken=json.loads(json.dumps(report))
    broken['trajectory'][1]['extra_quartics']['new']['real']=-.5
    with pytest.raises(ValueError,match='replay'):
        stage.verify_saved_trajectory(broken,solved,layout,('new',),('lambdaHHdagS2S2',))
    broken=json.loads(json.dumps(report))
    del broken['trajectory'][1]['reused_quartics']
    with pytest.raises(ValueError,match='reused quartics'):
        stage.verify_saved_trajectory(broken,solved,layout,('new',),('lambdaHHdagS2S2',))


def test_zero_extra_compact_trajectory_audit(monkeypatch):
    mod=types.ModuleType('Numerical.core.StateVector')
    mod.unpack_uv_state=lambda v,l,mu_gev: SimpleNamespace(sm=SimpleNamespace(
        gY=.35,g2=.6,g3=.9,lambdaH=.27))
    monkeypatch.setitem(sys.modules,'Numerical.core.StateVector',mod)
    solved=SimpleNamespace(y=np.full((1,3),.35),t=np.log([1e7,1e6,1e5]),nfev=15)
    rows=[{'mu_gev':mu,'gY':.35,'g2':.6,'g3':.9,'lambdaH':.27,
           'extra_quartics':{}} for mu in [1e7,1e6,1e5]]
    stage.verify_saved_trajectory({'trajectory':rows,'n_saved_points':3,'nfev':15},
                                  solved,SimpleNamespace(size=1),(),())


def test_summary_sha256_must_match(tmp_path):
    source=tmp_path/'uv.json';source.write_text('saved',encoding='utf-8')
    digest=sha256(source.read_bytes()).hexdigest()
    summary={'format':'stage31-16-model-uv-summary-v1',
             'status':'Success16ModelConditionalUVRunning',
             'total_models':16,'new_quartic_slots_total':38,
             'reused_quartic_slots_total':12,
             'physical_coupling_mapping_certified':False,
             'baseline_stage30_sha256':digest,'stage31_sha256':{'A_-2':digest}}
    stage._verify_summary(summary,'A_+0',source)
    stage._verify_summary(summary,'A_-2',source)
    summary['stage31_sha256']['A_-2']='0'*64
    with pytest.raises(ValueError,match='hash differs'):
        stage._verify_summary(summary,'A_-2',source)


def test_reused_quartics_and_unmapped_sidecar_not_misrepresented():
    script=(Path(__file__).resolve().parents[2]/'scripts/FermionThresholdHandoff.py').read_text()
    assert 'build_intermediate_scalar_boundary(' in script
    assert 'fermion_threshold_masses(' in script
    assert "'wilson_coefficient_matching_performed': False" in script
    assert "'new_quartics_integrated_into_intermediate_eft': False" in script
    assert 's30.combined_rhs(' in script
    assert 'verify_saved_trajectory(record, solved' in script


def test_handoff_calls_existing_threshold_and_retains_unmapped_new_quartic(tmp_path,monkeypatch):
    original=tmp_path/'stage30_A_p0_coupled_uv_running.json'
    original.write_text('{"saved":true}', encoding='utf-8')
    combined=tmp_path/'stage31_16_model_uv_summary.json'
    combined.write_text('{"summary":true}', encoding='utf-8')
    m30={'model':'A_+0','status':'SuccessOneModelNumericalUVRunning'}
    s31={'format':'stage31-16-model-uv-summary-v1'}
    monkeypatch.setattr(stage.s30,'read_json',lambda p: m30 if p==original else s31)
    monkeypatch.setattr(stage,'_verify_summary',lambda *args: None)
    monkeypatch.setattr(stage,'_verify_report',lambda *args: None)
    monkeypatch.setattr(stage,'verify_sources',lambda *args: None)
    monkeypatch.setattr(stage.s30,'quartic_names',lambda reg: ('newQ',))
    monkeypatch.setattr(stage.s30,'REGISTRY',original)
    uv=FakeUV(1e5,FakeSM(.37,.59,.88,.276, np.eye(2,dtype=complex)),
              np.diag([1e5,1.1e5]).astype(complex))
    intermediate=SimpleNamespace(mu_gev=1e5)
    # SimpleNamespace cannot be exported as a canonical EFT dataclass;
    # use a local dataclass to enforce the real serializer contract.
    @dataclass
    class Intermediate:
        mu_gev: float
        sm: FakeSM
    interim=Intermediate(1e5,uv.sm)
    monkeypatch.setattr(stage,'_replay',lambda *args: (uv,{'newQ':-.012+0.003j},(),SimpleNamespace(nfev=212)))
    backend=types.ModuleType('Numerical.running.FermionThresholdBoundary')
    called=[]
    backend.fermion_threshold_masses=lambda st: SimpleNamespace(masses_gev=np.array([1e5,1.1e5]))
    def project(st,*,matching_scale_gev):
        called.append((st,matching_scale_gev))
        return interim
    backend.build_intermediate_scalar_boundary=project
    monkeypatch.setitem(sys.modules,'Numerical.running',types.ModuleType('Numerical.running'))
    monkeypatch.setitem(sys.modules,'Numerical.running.FermionThresholdBoundary',backend)
    payload=stage.handoff('A_+0',folder=tmp_path)
    assert called==[(uv,1e5)]
    assert payload['new_quartic_sidecar_at_threshold']['newQ']=={'real':-.012,'imag':.003}
    assert 'newQ' not in payload['intermediate_scalar_boundary']
    assert payload['wilson_coefficient_matching_performed'] is False
    assert payload['uv_threshold_state']['MF']['real'][0][0]==1e5
    assert (tmp_path/'stage33_A_p0_fermion_threshold_handoff.json').is_file()
