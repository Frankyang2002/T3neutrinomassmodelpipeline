"""Stage 36 offline source-selection, matching and RGE handoff contracts."""
from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
import sys
import types
from types import SimpleNamespace

import numpy as np
import pytest

from scripts import Stage36FinalWeinbergRunning as stage


def _source_file(root, model, legacy=False):
    x, alpha = model.split('_')
    token = ('p' if alpha[0] == '+' else 'm') + alpha[1:]
    base = root / 'output'
    if legacy:
        base = base / 'hypercharge'
    else:
        base = base / 'full' / 'hypercharge'
    return base / f'T3_{x}_alpha_{token}' / 'data' / (
        'c5_coefficient.txt' if legacy else 'final_weinberg_coefficient.json')


def _mk_source(root, model):
    legacy = model in stage.LEGACY_MODELS
    path = _source_file(root, model, legacy)
    path.parent.mkdir(parents=True, exist_ok=True)
    if legacy:
        path.write_text('Matchete legacy coefficient')
    else:
        path.write_text(json.dumps({
            'status': 'Success', 'scheme': 'fixed_order_one_loop_MSbar_hierarchical',
            'renormalized_hard_threshold': {'status': 'Success'},
            'combined': {'ready_for_physical_majorana_numerics': True},
            'direct_running_contribution': {'manifestly_symmetric_under_pq': True},
        }))
    return path


def test_inventory_and_source_rule(tmp_path):
    assert len(stage.MODELS) == 16
    assert len(stage.LEGACY_MODELS) == 4
    assert set(stage.LEGACY_MODELS) < set(stage.MODELS)
    for model in stage.MODELS:
        p = _mk_source(tmp_path, model)
        if model in stage.LEGACY_MODELS:
            with pytest.raises(ValueError, match='allow-legacy'):
                stage.source_path(model, root=tmp_path)
        actual, kind = stage.source_path(model, root=tmp_path, allow_legacy=True)
        assert actual == p
        assert ('legacy' in kind) == (model in stage.LEGACY_MODELS)
    with pytest.raises(ValueError, match='Unsupported'):
        stage.source_path('UNKNOWN', root=tmp_path)
    assert stage.MU_LOW == 100.0


def test_matching_is_model_specific(tmp_path):
    _mk_source(tmp_path, 'A_+0')
    with pytest.raises(ValueError, match='model-specific matching artifact missing'):
        stage.source_path('D_+0', root=tmp_path)
    with pytest.raises(ValueError, match='model-specific matching artifact missing'):
        stage.source_path('A_-2', root=tmp_path)


def test_hierarchical_validation_fails_closed(tmp_path):
    path = _mk_source(tmp_path, 'A_+0')
    stage.verify_hierarchical_artifact(path)
    for key, value in [('status', 'Failed'), ('scheme', 'wrong')]:
        obj = json.loads(path.read_text())
        obj[key] = value
        path.write_text(json.dumps(obj))
        with pytest.raises(ValueError, match='not validated'):
            stage.verify_hierarchical_artifact(path)
        path = _mk_source(tmp_path, 'A_+0')
    obj = json.loads(path.read_text())
    obj['combined']['ready_for_physical_majorana_numerics'] = False
    path.write_text(json.dumps(obj))
    with pytest.raises(ValueError, match='not validated'):
        stage.verify_hierarchical_artifact(path)


@pytest.mark.parametrize('invalid', [
    np.ones((2, 2)),
    np.array([[1, 0, 0], [0.2, 1, 0], [0, 0, 1]], complex),
    np.zeros((3, 3)),
    np.full((3, 3), np.nan),
])
def test_c5_rejects_bad_matrices(invalid):
    with pytest.raises(ValueError, match='Stage36 fail-closed'):
        stage._valid_c5(invalid, 'C5')


def _install_module(monkeypatch, name, **entries):
    module = types.ModuleType(name)
    for k, v in entries.items():
        setattr(module, k, v)
    monkeypatch.setitem(sys.modules, name, module)


@pytest.mark.parametrize('model', ['A_+0', 'A_-4'])
def test_rge_uses_existing_services_with_correct_c5_source(tmp_path, monkeypatch, model):
    _mk_source(tmp_path, model)
    report33 = tmp_path/'stage33.json'
    report35 = tmp_path/'stage35.json'
    summary = tmp_path/'summary.json'
    for path in (report33, report35, summary):
        path.write_text('{}')
    uv = SimpleNamespace(MF=np.diag([100000., 110000., 120000.]))
    eft = SimpleNamespace(mu_gev=1000.)
    sm = SimpleNamespace(mu_gev=1000., sm=SimpleNamespace(lambdaH=.336))
    monkeypatch.setattr(stage, 'checked_inputs',
                        lambda m, f: (uv, eft, sm, report33, report35, summary))
    calls = []
    class MatchingInput:
        heavy_masses_gev = np.array([100000., 110000., 120000.])
        scalar_mass1_gev = 1130.
        scalar_mass2_gev = 1050.
        lambda_t3 = .1+0j
        def as_config(self):
            return {'t3': {'lambdaT3': self.lambda_t3}}
    def inputs(t):
        assert t.uv_threshold_state is uv and t.eft1_threshold_state is eft
        assert t.mu_scalar_threshold_gev == 1000.
        return MatchingInput()
    def eval_hier(path, config):
        calls.append(('hier', str(path)))
        return np.eye(3, dtype=complex)*1e-16
    def eval_legacy(path, config):
        calls.append(('legacy', str(path)))
        return np.eye(3, dtype=complex)*1.1e-16
    def build(boundary, K):
        assert boundary is sm
        calls.append(('boundary', float(np.linalg.norm(K))))
        return SimpleNamespace(K=K)
    def evolve(initial, start, stop, *, rtol, atol):
        assert start == 1000. and stop == 100. and rtol == 1e-8
        calls.append(('evolve', stop))
        return SimpleNamespace(K=initial.K*1.01, ye=np.eye(3)*.001,
                               gY=.35, g2=.6, g3=1.1, lambdaH=.37,
                               solver_success=True, mu_initial=1000.,
                               mu_final=100., nfev=35)
    def convert(K, *, vev_gev):
        assert vev_gev == stage.VEV
        return -(vev_gev**2/2.)*K
    def rotate(mass, ye):
        return mass
    def observables(mass, *, ordering):
        assert ordering == 'NO'
        return SimpleNamespace(masses_ev=np.array([.01, .012, .051]),
                               delta_m21_sq_ev2=4.4e-5,
                               delta_m31_sq_ev2=2.5e-3,
                               delta_m32_sq_ev2=2.46e-3,
                               pmns_abs=np.eye(3), takagi_residual=1e-12)
    _install_module(monkeypatch, 'Numerical.running.FinalC5TrajectoryAdapter',
                    extract_final_c5_inputs=inputs)
    _install_module(monkeypatch, 'Numerical.running.SMWeinbergStage',
                    evaluate_final_weinberg_json=eval_hier,
                    evaluate_symbolic_c5=eval_legacy)
    _install_module(monkeypatch, 'Numerical.running.ScalarThresholdBoundary',
                    build_sm_weinberg_initial_conditions=build)
    _install_module(monkeypatch, 'Numerical.running.SMWeinbergEvolution',
                    evolve_sm_weinberg=evolve)
    _install_module(monkeypatch, 'physics.NeutrinoMass',
                    numerical_neutrino_mass_matrix=convert)
    _install_module(monkeypatch, 'physics.NeutrinoTrajectory',
                    charged_lepton_mass_basis_matrix=rotate)
    _install_module(monkeypatch, 'physics.NeutrinoObservables',
                    calculate_neutrino_observables=observables)
    result = stage.run_model(model, folder=tmp_path, root=tmp_path, allow_legacy=True)
    assert result['status'] == stage.STATUS
    assert result['solver_nfev'] == 35
    assert result['physical_coupling_mapping_certified'] is False
    assert result['new_general_quartics_included_in_matching'] is False
    assert len(result['neutrino_observables']['masses_ev']) == 3
    assert calls[0][0] == ('legacy' if model in stage.LEGACY_MODELS else 'hier')
    assert calls[1][0] == 'boundary' and calls[2][0] == 'evolve'
    path = tmp_path / f'stage36_{stage.s35.s34.safe_model(model)}_low_energy_neutrinos.json'
    assert path.is_file()


def test_matching_does_not_allow_legacy_without_opt_in(tmp_path, monkeypatch):
    _mk_source(tmp_path, 'A_-4')
    dummy = tmp_path/'dummy'
    dummy.write_text('{}')
    monkeypatch.setattr(stage, 'checked_inputs',
                        lambda m, f: (None, None, None, dummy, dummy, dummy))
    with pytest.raises(ValueError, match='allow-legacy'):
        stage.run_model('A_-4', folder=tmp_path, root=tmp_path)


def test_source_code_uses_project_production_services():
    code = Path(stage.__file__).read_text()
    assert 'evaluate_final_weinberg_json' in code
    assert 'evaluate_symbolic_c5' in code
    assert 'build_sm_weinberg_initial_conditions' in code
    assert 'evolve_sm_weinberg' in code
    assert 'charged_lepton_mass_basis_matrix' in code
    assert 'numerical_neutrino_mass_matrix' in code
    assert 'sm_weinberg_running_performed' in code
