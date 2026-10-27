"""Stage 35 offline contracts: no live RGBeta evaluation is required."""
from __future__ import annotations

from dataclasses import dataclass, field
from hashlib import sha256
import json
import sys
import types
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from scripts import Stage35ScalarThresholdHandoff as stage


@dataclass
class DummyRep:
    d_s1: int = 1
    d_s2: int = 3
    d_f: int = 2
    alpha: int = 0
    shared_scalar: bool = False


@dataclass
class DummySM:
    gY: float = .37
    g2: float = .59
    g3: float = .88
    lambdaH: float = .336
    yu: np.ndarray = field(default_factory=lambda: np.eye(3))
    yd: np.ndarray = field(default_factory=lambda: np.eye(3)*.01)
    ye: np.ndarray = field(default_factory=lambda: np.eye(3)*.001)


@dataclass
class DummyState:
    mu_gev: float = stage.START
    representation: DummyRep = field(default_factory=DummyRep)
    sm: DummySM = field(default_factory=DummySM)
    mS1Sq: float = 1.21e6
    mS2Sq: float = 1e6
    lambdaT3: complex = .1 + 0j

    def validated(self):
        return self


@dataclass
class DummyFinalSM:
    mu_gev: float
    sm: DummySM


def _snap(state, reused):
    return {'mu_gev': state.mu_gev, 'gY': state.sm.gY,
            'g2': state.sm.g2, 'g3': state.sm.g3,
            'lambdaH': state.sm.lambdaH,
            'reused_quartics': {x: {'real': .01, 'imag': 0} for x in reused}}


def _write(path, obj):
    path.write_text(json.dumps(obj, indent=2))


def _fixture(tmp_path, monkeypatch):
    # Mock only the canonical decoder and threshold service. All Stage34 JSON,
    # SHA-256, mass/snapshot and model-inventory logic in Stage35 is real.
    monkeypatch.setattr(stage.s34, 'decode_intermediate_state',
                        lambda raw: DummyState(mu_gev=float(raw['mu_gev'])))
    monkeypatch.setattr(stage.s34, '_snapshot', _snap)
    boundary_module = types.ModuleType('Numerical.running.ScalarThresholdBoundary')
    called = []
    def project(st, *, matching_scale_gev):
        called.append((st.mu_gev, matching_scale_gev))
        return DummyFinalSM(mu_gev=matching_scale_gev, sm=st.sm)
    boundary_module.build_final_sm_boundary = project
    boundary_module.scalar_threshold_masses = lambda st: SimpleNamespace(masses_gev=(1100., 1000.))
    monkeypatch.setitem(sys.modules, 'Numerical.running.ScalarThresholdBoundary', boundary_module)
    model, tag = 'A_+0', 'A_p0'
    init = DummyState(mu_gev=stage.START)
    end = DummyState(mu_gev=stage.STOP)
    extras = {'lambdaT3GeneralQ000211I1': {'real': -.01, 'imag': 0}}
    dummy_inputs = {}
    for i in range(5):
        p = tmp_path / f'dummy_source_{i}.txt'
        p.write_text(str(i))
        dummy_inputs[str(p)] = stage.sha(p)
    report = {
        'format': stage.s34.FORMAT,
        'status': 'SuccessSupportedRenormalisableIntermediateRunning',
        'model': model,
        'mu_initial_gev': stage.START,
        'mu_final_gev': stage.STOP,
        'one_loop_scalar_only_rgbeta_evolution_performed': True,
        'native_quartic_class_projected': True,
        'eft1_field_content': ['SM', 'S1', 'S2'],
        'active_wilson_operator_running_performed': False,
        'fermion_wilson_matching_performed': False,
        'physical_coupling_mapping_certified': False,
        'new_stage29_quartics_evolved_in_eft1': False,
        'truncated_new_quartic_beta_system': True,
        'excluded_new_quartic_symbols': list(extras),
        'reused_quartics_evolved_in_eft1': [],
        'new_quartic_sidecar_at_mu_f_only': extras,
        'solver_nfev': 47,
        'n_saved_points': 3,
        'intermediate_initial_state': {'mu_gev':stage.START,'sm':{'lambdaH':init.sm.lambdaH}},
        'intermediate_final_state': {'mu_gev':stage.STOP,'sm':{'lambdaH':end.sm.lambdaH}},
        'scalar_threshold_singular_masses_gev': [1100., 1000.],
        'trajectory': [_snap(init, ()), {'mu_gev': 1e4}, _snap(end, ())],
        'input_sha256': dummy_inputs,
    }
    path = tmp_path / f'stage34_{tag}_intermediate_eft_running.json'
    _write(path, report)
    rows = {m:{'report_sha256': 'f'*64, 'new_quartics_unmatched': 0,
               'reused_quartics_evolved': 0, 'lambdaH_initial': .276,
               'lambdaH_final': .336, 'solver_nfev': 47} for m in stage.MODELS}
    rows[model] = {'report_sha256': stage.sha(path), 'new_quartics_unmatched':1,
                   'reused_quartics_evolved':0,'lambdaH_initial': init.sm.lambdaH,
                   'lambdaH_final':end.sm.lambdaH,'solver_nfev':47}
    summary = {
        'format': 'stage34-16-model-supported-intermediate-running-v1',
        'status': 'Success16SupportedIntermediateSegments',
        'models_checked': 16,
        'mu_initial_gev': stage.START, 'mu_final_gev': stage.STOP,
        'new_quartics_not_in_intermediate_eft':38,'reused_quartics_evolved':12,
        'physical_coupling_mapping_certified': False,
        'wilson_coefficient_matching_performed':False, 'models':rows,
    }
    _write(tmp_path / 'stage34_16_model_intermediate_eft_summary.json', summary)
    return report, summary, path, dummy_inputs, called


def test_stage35_inventory():
    assert len(stage.MODELS) == 16
    assert stage.MODELS == stage.s34.MODELS
    assert stage.STOP == 1000
    assert stage.FORMAT.startswith('stage35-')
    assert stage.s34.safe_model('E_-4') == 'E_m4'


def test_snapshot_and_number_guards():
    stage._compare_json({'a': {'real': .3, 'imag': -.2}},
                        {'a': {'real': .3, 'imag': -.2}}, 'a')
    with pytest.raises(ValueError, match='snapshot key inventory'):
        stage._compare_json({'a': 1}, {'b': 1}, 'a')
    with pytest.raises(ValueError, match='state='):
        stage._compare_json({'a': .1}, {'a': .4}, 'a')
    with pytest.raises(ValueError, match='numeric'):
        stage._number(True, 'bad')
    with pytest.raises(ValueError, match='nonfinite'):
        stage._number(float('nan'), 'bad')


def test_scalar_projection_integrates_existing_threshold_api(tmp_path, monkeypatch):
    report, summary, path, inputs, called = _fixture(tmp_path, monkeypatch)
    result = stage.run_model('A_+0', folder=tmp_path)
    assert called == [(1000., 1000.)]
    assert result['status'] == stage.STATUS
    assert result['scalar_threshold_masses_gev'] == [1100., 1000.]
    assert result['final_sm_boundary_at_mu_s']['sm']['lambdaH'] == .336
    assert result['source_stage34_report_sha256'] == stage.sha(path)
    assert result['new_quartics_matched_or_evolved_below_mu_f'] is False
    assert result['scalar_wilson_matching_performed'] is False
    assert result['low_energy_weinberg_running_performed'] is False
    assert result['new_quartic_sidecar_at_mu_f_only'] == report['new_quartic_sidecar_at_mu_f_only']
    assert (tmp_path/'stage35_A_p0_scalar_threshold_handoff.json').is_file()


@pytest.mark.parametrize('mutate,pattern', [
    (lambda r: r.update(status='Failed'), 'individual report'),
    (lambda r: r.update(fermion_wilson_matching_performed=True), 'individual report'),
    (lambda r: r.update(new_stage29_quartics_evolved_in_eft1=True), 'individual report'),
    (lambda r: r.update(mu_final_gev=1100), 'individual report'),
    (lambda r: r.update(n_saved_points=5), 'saved trajectory'),
    (lambda r: r['trajectory'][-1].update(lambdaH=.4), 'state='),
    (lambda r: r.update(scalar_threshold_singular_masses_gev=[1100, 999]), 'state='),
    (lambda r: r['intermediate_final_state'].update(mu_gev=999), 'final-state scale'),
    (lambda r: r.update(excluded_new_quartic_symbols=[]), 'inventory mismatch'),
    (lambda r: r['trajectory'][1].update(mu_gev=2e5), 'not descending'),
])
def test_rejects_corrupted_report(tmp_path, monkeypatch, mutate, pattern):
    report, summary, path, inputs, called = _fixture(tmp_path, monkeypatch)
    mutate(report)
    _write(path, report)
    summary['models']['A_+0']['report_sha256'] = stage.sha(path)
    _write(tmp_path/'stage34_16_model_intermediate_eft_summary.json', summary)
    with pytest.raises(ValueError, match=pattern):
        stage.run_model('A_+0', folder=tmp_path)


def test_tampered_source_hash_fails(tmp_path, monkeypatch):
    report, summary, path, sources, called = _fixture(tmp_path, monkeypatch)
    first = Path(next(iter(sources)))
    first.write_text('tampered')
    with pytest.raises(ValueError, match='source file changed'):
        stage.run_model('A_+0', folder=tmp_path)


def test_stale_stage34_summary_hash_fails(tmp_path, monkeypatch):
    report, summary, path, inputs, called = _fixture(tmp_path, monkeypatch)
    path.write_text(path.read_text()+' ')
    with pytest.raises(ValueError, match='SHA-256 differs from summary'):
        stage.run_model('A_+0', folder=tmp_path)


def test_no_silent_missing_stage34_reports(tmp_path):
    with pytest.raises(FileNotFoundError):
        stage.write_summary(folder=tmp_path)


def test_invalid_model_rejected(tmp_path):
    with pytest.raises(ValueError, match='Unknown ordinary'):
        stage.run_model('../A_+0', folder=tmp_path)


def test_stage35_source_uses_verified_existing_threshold():
    source = (Path(__file__).resolve().parents[2]/'scripts/Stage35ScalarThresholdHandoff.py').read_text()
    assert 'build_final_sm_boundary(final_state, matching_scale_gev=STOP)' in source
    assert 'scalar_threshold_masses(final_state)' in source
    assert 'Stage34IntermediateEFTRunning' in source
    assert 'No C5/C12 Wilson matching' in source
    assert 'run_intermediate_scalar_segment' not in source
