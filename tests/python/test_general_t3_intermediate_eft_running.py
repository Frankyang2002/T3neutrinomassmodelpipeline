"""Stage 34 offline contracts. Live RGBeta/production EFT1 are run on Windows."""
from __future__ import annotations

from dataclasses import dataclass, field
import json
from pathlib import Path
from hashlib import sha256
from types import SimpleNamespace
import sys
import types

import numpy as np
import pytest

from scripts import IntermediateEFTRunning as stage


@dataclass
class DummyRep:
    d_s1: int = 1
    d_s2: int = 3
    d_f: int = 2
    alpha: int = 0
    shared_scalar: bool = False


@dataclass
class DummySM:
    gY: float = 0.37
    g2: float = 0.60
    g3: float = 0.89
    lambdaH: float = 0.28


@dataclass
class DummyState:
    mu_gev: float = 1e5
    representation: DummyRep = field(default_factory=DummyRep)
    sm: DummySM = field(default_factory=DummySM)
    mS1Sq: float = 1.4e6
    mS2Sq: float = 1.1e6
    lambdaS1: float = 0.1
    lambdaS2: float = 0.12
    lambdaH1: float = 0.02
    lambdaH2: float = 0.03
    lambda12: float = 0.01
    lambdaT3: complex = 0.1+0j
    lambdaHHdagS2S2: complex | None = None

    def validated(self):
        return self


def _beta(keys=('gY', 'g2')):
    return {
        'format': stage.BETA_FORMAT,
        'status': 'Success',
        'physical_coupling_mapping_certified': False,
        'new_stage29_quartics_integrated': False,
        'metadata': {'SharedScalar': False, 'dS1': 1, 'dS2': 3,
                     'dF': 2, 'alpha': 0, 'IntegratedField': 'F',
                     'ActiveBSMFields': ['S1', 'S2'],
                     'ReportBetaConvention':
                     '16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g'},
        'report_betas': {k: '1' for k in keys},
        'beta_symbols': list(keys),
    }


def _stub_numerical(monkeypatch, *, lambda_error=False):
    vecmod = types.ModuleType('Numerical.running.IntermediateScalarStateVector')
    vecmod.pack_intermediate_scalar_state = lambda state: (np.array([0.4, 0.5]), SimpleNamespace(names=('gY', 'g2')))
    evalmod = types.ModuleType('Numerical.running.IntermediateScalarRGBetaEvaluator')
    evalmod.EXPECTED_CONVENTION = '16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g'
    monkeypatch.setitem(sys.modules, 'Numerical.running.IntermediateScalarStateVector', vecmod)
    monkeypatch.setitem(sys.modules, 'Numerical.running.IntermediateScalarRGBetaEvaluator', evalmod)


def test_stage34_inventory_and_scales():
    assert len(stage.MODELS) == 16
    assert set(stage.MODELS) == set(stage.s33.MODELS)
    assert stage.START == 100_000
    assert stage.STOP == 1_000
    assert stage.POINTS == 24
    assert stage.safe_model('A_+0') == 'A_p0'
    assert stage.safe_model('B_-1') == 'B_m1'
    with pytest.raises(ValueError, match='Unknown'):
        stage.safe_model('../A_+0')


def test_complex_scalar_and_matrix_validation():
    assert stage._complex({'real': -0.01, 'imag': 0.02}, 'x') == -.01+.02j
    assert np.all(stage._matrix({'real': np.eye(3).tolist(),
                                 'imag': np.zeros((3,3)).tolist()}, 'm') == np.eye(3))
    for raw in [{'real': 1}, {'real': float('nan'), 'imag': 0},
                {'real': 0, 'imag': float('inf')}]:
        with pytest.raises(ValueError, match='Stage34'):
            stage._complex(raw, 'bad')
    with pytest.raises(ValueError, match='expected finite 3x3'):
        stage._matrix({'real': [[1]], 'imag': [[0]]}, 'bad')
    with pytest.raises(ValueError, match='expected real'):
        stage._scalar(True, 'bad')


def test_beta_inventory_validation_passes(monkeypatch):
    _stub_numerical(monkeypatch)
    assert stage.validate_beta_payload(_beta(), DummyState()) == ('gY', 'g2')


@pytest.mark.parametrize('change', [
    {'status': 'Failed'},
    {'physical_coupling_mapping_certified': True},
    {'new_stage29_quartics_integrated': True},
    {'beta_symbols': ['gY']},
    {'report_betas': {'gY': '1'}},
])
def test_beta_inventory_rejects_inconsistent_payload(monkeypatch, change):
    _stub_numerical(monkeypatch)
    with pytest.raises(ValueError, match='Stage34'):
        stage.validate_beta_payload(dict(_beta(), **change), DummyState())


def test_bad_eft_beta_convention_and_representation(monkeypatch):
    _stub_numerical(monkeypatch)
    for key, value in [('alpha', -2), ('ReportBetaConvention', 'wrong'),
                       ('ActiveBSMFields', ['S1', 'S2', 'F'])]:
        obj = _beta()
        obj['metadata'][key] = value
        with pytest.raises(ValueError, match='metadata mismatch'):
            stage.validate_beta_payload(obj, DummyState())


def test_unresolved_beta_is_rejected(monkeypatch):
    _stub_numerical(monkeypatch)
    obj = _beta()
    obj['report_betas']['gY'] = 'BetaTerm[gY,1]'
    with pytest.raises(ValueError, match='beta symbol inventory mismatch'):
        stage.validate_beta_payload(obj, DummyState())


def test_boundary_hash_provenance_fail_closed(tmp_path, monkeypatch):
    # Use a deliberately fake but internally consistent 16-model summary.
    m='A_+0'; t='A_p0'
    uv=(tmp_path/'stage30_A_p0_coupled_uv_running.json'); uv.write_text('uv')
    us=(tmp_path/'stage31_16_model_uv_summary.json'); us.write_text('{}')
    uvsha=stage.sha(uv)
    raw={'format':stage.s33.FORMAT,
         'status':'SuccessRenormalisableFermionThresholdProjection',
         'model':m,'uv_running_replayed_and_trajectory_verified':True,
         'new_quartics_integrated_into_intermediate_eft':False,
         'wilson_coefficient_matching_performed':False,
         'intermediate_eft_running_performed':False,
         'physical_coupling_mapping_certified':False,
         'mu_fermion_threshold_gev':1e5,
         'source_report':'stage30_A_p0_coupled_uv_running.json',
         'source_report_sha256':uvsha,
         'source_summary_sha256':stage.sha(us),
         'heavy_fermion_singular_masses_gev':[1e5,1.1e5,1.2e5],
         'new_quartic_sidecar_at_threshold':{'lambdaT3GeneralQ000211I1':{'real':-.01,'imag':0}},
         'reused_quartics_in_intermediate_state':[],
         'intermediate_scalar_boundary':{'faked':True}}
    (tmp_path/f'stage33_{t}_fermion_threshold_handoff.json').write_text(json.dumps(raw))
    rows={n:{'status':raw['status'], 'report_sha256':uvsha,
             'heavy_fermion_singular_masses_gev':raw['heavy_fermion_singular_masses_gev'],
             'new_sidecar_count':1,'reused_count':0} for n in stage.MODELS}
    ss={'format':'stage33-16-model-fermion-threshold-handoff-v1',
        'status':'Success16RenormalisableThresholdProjections','models_checked':16,
        'physical_coupling_mapping_certified':False,
        'wilson_coefficient_matching_performed':False,'models':rows}
    (tmp_path/'stage33_16_model_fermion_threshold_summary.json').write_text(json.dumps(ss))
    monkeypatch.setattr(stage.s33, '_verify_summary', lambda *args: None)
    monkeypatch.setattr(stage, 'decode_intermediate_state', lambda x: DummyState())
    individual,state,extras,reused,hashes=stage.verify_boundary(m,tmp_path)
    assert len(extras) == 1 and reused == ()
    assert isinstance(state,DummyState)
    assert len(hashes) == 4
    uv.write_text('tampered')
    with pytest.raises(ValueError, match='SHA-256 mismatch'):
        stage.verify_boundary(m,tmp_path)
    uv.write_text('uv')
    raw['source_report'] = '../../untrusted.json'
    (tmp_path/f'stage33_{t}_fermion_threshold_handoff.json').write_text(json.dumps(raw))
    with pytest.raises(ValueError, match='source UV report path'):
        stage.verify_boundary(m,tmp_path)


def test_running_uses_native_intermediate_runner(monkeypatch,tmp_path):
    _stub_numerical(monkeypatch)
    state=DummyState()
    extras={'lambdaT3GeneralQ000211I1':{'real':-.01,'imag':0}}
    monkeypatch.setattr(stage, 'verify_boundary', lambda m, f: (
        {'new_quartic_sidecar_at_threshold':extras}, state,extras,(),{'src':'0'*64}))
    report_file=tmp_path/'stage34_A_p0_eft1_full_betas.json'
    report_file.write_text(json.dumps(_beta()))
    calls=[]
    runner=types.ModuleType('Numerical.running.IntermediateScalarRunner')
    class Solved:
        solver_success=True
        nfev=200
        mu_gev=np.geomspace(1e5,1e3,24)
        y=np.array([[.4]*24,[.5]*24])
        initial_vector=np.array([.4,.5])
        final_state=DummyState(mu_gev=1e3)
        def state_at_index(self, i):
            return DummyState(mu_gev=float(self.mu_gev[i]))
    def run(x,payload,mu_final_gev,**opts):
        calls.append((x,payload,mu_final_gev,opts))
        return Solved()
    runner.run_intermediate_scalar_segment=run
    monkeypatch.setitem(sys.modules,'Numerical.running.IntermediateScalarRunner',runner)
    mod=types.ModuleType('Numerical.running.ScalarThresholdBoundary')
    mod.scalar_threshold_masses=lambda st:SimpleNamespace(masses_gev=(1100.,1000.))
    monkeypatch.setitem(sys.modules,'Numerical.running.ScalarThresholdBoundary',mod)
    # s33.json_state accepts dataclasses but DummyRep, DummySM must remain simple.
    result=stage.run_model('A_+0',folder=tmp_path)
    assert calls and calls[0][2] == 1e3
    assert result['n_saved_points']==24
    assert result['new_stage29_quartics_evolved_in_eft1'] is False
    assert result['truncated_new_quartic_beta_system'] is True
    assert result['excluded_new_quartic_symbols'] == list(extras)
    assert result['active_wilson_operator_running_performed'] is False
    assert result['solver_nfev']==200
    assert (tmp_path/'stage34_A_p0_intermediate_eft_running.json').exists()


def test_wolfram_exporter_integrity():
    root=Path(__file__).resolve().parents[2]
    src=(root/'scripts/ExportGeneralT3EFT1OneLoop.wl').read_text()
    assert 'RGBeta`BetaTerm[RGBeta`Quartic,1]' in src
    assert 'T3RGBetaBuildEFT1[' in src
    assert 'new_stage29_quartics_integrated' in src
    assert 'Check[Catch' not in src
    assert 'T3RGBetaWriteJSON' in src
    runner=(root/'scripts/RunGeneralT3IntermediateEFTRunning.ps1').read_text()
    assert 'ExportGeneralT3EFT1OneLoop.wl' in runner
    assert 'IntermediateEFTRunning' in runner


def test_summary_rejects_missing_reports(tmp_path):
    with pytest.raises(FileNotFoundError):
        stage.write_summary(folder=tmp_path)
