"""Fast tests for the thesis-input pilot (no Matchete/Wolfram rerun)."""
from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / 'scripts' / 'RunT3ThesisTwoModelPilot.py'
spec = importlib.util.spec_from_file_location('thesis_pilot', SCRIPT)
assert spec and spec.loader
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)


def test_actual_thesis_table_and_seeds():
    cfg = module.load_config(ROOT, module.CONFIG)
    assert cfg['scales_GeV'] == {'UV':1e7, 'F':1e5,'S':1e3,'low':100.0}
    assert cfg['lambdaT3_UV'] == 0.1
    assert cfg['class_E_extra_half_factor'] is False
    for model in module.MODELS:
        y = module.initial_parameters(cfg, model)
        assert y.shape == (18,)
        assert np.max(abs(y)) < 1e-3
    d = module.initial_parameters(cfg, 'D_-2')
    a = module.initial_parameters(cfg, 'A_-4')
    assert np.isclose(d[0], -0.13438183e-4)
    assert np.isclose(d[9], 1.36211245e-4)
    assert np.isclose(a[0], -0.26606112e-4)
    assert np.isclose(a[17], 0.67251936e-4)


def test_benchmark_input_overrides_no_hidden_class_E_half():
    cfg = module.load_config(ROOT, module.CONFIG)
    names = ['gY','g2','g3','lambdaH','lambdaS1','lambdaS2',
             'lambdaT3','lambdaT3GeneralQ001102I1','muT3GeneralS2Sq',
             'MF','mS1Sq','mS2Sq','y1','y2','y1Conj','yu','yd','ye']
    prior = {k: None for k in names}
    for model in module.MODELS:
        v = module.initial_parameters(cfg, model)
        result = module.overrides_for(cfg, model, names, v, prior)
        assert result['lambdaT3'] == 0.1
        assert result['lambdaT3GeneralQ001102I1'] == .01
        assert result['lambdaS1'] == .01
        assert result['muT3GeneralS2Sq'] == {'real':0.,'imag':0.}
        assert result['gY'] == 0.38265
        assert result['MF'][1][1]['real'] == 110000.0
        assert result['y1'][0][0]['real'] == v[0]
        assert result['y2'][0][0]['real'] == v[9]
        assert result['yu'] is None  # existing SM Yukawa initializer remains intact


def test_normalised_chi2_conventions():
    target = module.read_json(ROOT / 'configs/nufit_6_1_ic24_no.json')
    pred = dict(target['central_values'])
    assert np.allclose(module.residuals(pred, target), 0.)
    pred['sin2_theta13'] += target['one_sigma_errors']['sin2_theta13']
    assert np.isclose(module.residuals(pred, target) @ module.residuals(pred, target), 1.)


def test_sparse_workstream_b_prior_accepts_source_declared_scalar_couplings():
    """Regression: Workstream B _bench_overrides() is a sparse mapping."""
    cfg = module.load_config(ROOT, module.CONFIG)
    names = ['MF', 'gY', 'g2', 'g3', 'lambdaH', 'lambdaT3',
             'mS1Sq', 'mS2Sq', 'y1', 'y2', 'yu', 'yd', 'ye',
             'lambdaS1', 'lambdaS2', 'lambdaH1', 'lambdaH2', 'lambda12',
             'lambdaT3GeneralQ001102I1', 'muT3GeneralS2Sq']
    # The real runner does not explicitly override lambdaS1/H1/12 or the
    # optional bilinear, even when those coordinates are in the native beta.
    prior = {k: None for k in names if k not in (
        'lambdaS1', 'lambdaS2', 'lambdaH1', 'lambdaH2', 'lambda12',
        'muT3GeneralS2Sq')}
    for model in module.MODELS:
        got = module.overrides_for(cfg, model, names,
                                   module.initial_parameters(cfg, model), prior)
        assert set(got) <= set(names)
        assert set(prior) < set(got)
        assert all(got[k] == 0.01 for k in (
            'lambdaS1', 'lambdaS2', 'lambdaH1', 'lambdaH2', 'lambda12'))
        assert got['lambdaT3'] == 0.1
        assert got['muT3GeneralS2Sq'] == {'real':0., 'imag':0.}


def test_unknown_uv_source_key_is_rejected():
    cfg = module.load_config(ROOT, module.CONFIG)
    names = ['MF', 'gY', 'g2', 'g3', 'lambdaH', 'lambdaT3',
             'mS1Sq', 'mS2Sq', 'y1', 'y2']
    try:
        module.overrides_for(cfg, 'A_-4', names,
                             module.initial_parameters(cfg, 'A_-4'),
                             {'invented_parameter': 0.01})
    except ValueError as exc:
        assert 'unknown UV beta keys' in str(exc)
    else:
        raise AssertionError('Unrecognised beta key was incorrectly allowed')


def test_failed_evaluation_is_retried_but_completed_is_cached(tmp_path, monkeypatch):
    cfg = module.load_config(ROOT, module.CONFIG)
    tgt = module.read_json(ROOT / cfg['fitting']['target'])
    ev = module.Evaluator(tmp_path, cfg, 'D_-2', tgt, False, 3, [])
    x = module.initial_parameters(cfg, 'D_-2')
    count = {'n':0}

    def fake_chain(root, cfg, model, vec, digest, *, keep_work):
        count['n'] += 1
        if count['n'] == 1:
            raise ValueError('synthetic temporary input failure')
        return {'prediction': dict(tgt['central_values']),
                'masses_eV':[0.01, 0.012, 0.05], 'scenario':'analytic_triangle',
                'physics_status':module.SUCCESS_B2}

    monkeypatch.setattr(module, '_run_full_chain', fake_chain)
    first = ev.evaluate(x)
    assert first['status'] == 'BLOCKED'
    second = ev.evaluate(x)
    assert second['status'] == 'COMPLETE_CONDITIONAL'
    assert second['chi2_approx'] == 0
    third = ev.evaluate(x)
    assert third['status'] == 'COMPLETE_CONDITIONAL'
    assert count['n'] == 2
    assert ev.attempts == 2


def test_model_specific_exit_ignores_other_historic_blocker(tmp_path, monkeypatch):
    """D_-2 success must return zero even if A_-4 previously blocked."""
    import json
    import shutil

    (tmp_path / 'configs' / 'thesis_numerics').mkdir(parents=True)
    shutil.copy(ROOT / module.CONFIG, tmp_path / module.CONFIG)
    shutil.copy(ROOT / 'configs/nufit_6_1_ic24_no.json',
                tmp_path / 'configs/nufit_6_1_ic24_no.json')
    other = tmp_path / module.OUTPUT / 'A_m4' / 'pilot_summary.json'
    other.parent.mkdir(parents=True)
    other.write_text(json.dumps({'status': 'BASELINE_BLOCKED', 'best': None}), encoding='utf-8')
    monkeypatch.setattr(module, 'ROOT', tmp_path)

    def fake_evaluate(self, x):
        return {'status': 'COMPLETE_CONDITIONAL', 'chi2_approx': 85161.4,
                'candidate_digest': 'synthetic_digest',
                'candidate': np.asarray(x, float).tolist()}

    monkeypatch.setattr(module.Evaluator, 'evaluate', fake_evaluate)
    rc = module.main(['--baseline', '--model', 'D_-2'])
    assert rc == 0
    overview = (tmp_path / module.OUTPUT / 'two_model_overview.csv').read_text()
    assert 'A_-4,BASELINE_BLOCKED' in overview
    assert 'D_-2,BASELINE_COMPLETE' in overview


def test_sobol_refinement_calls_least_squares_and_improves():
    """Synthetic five-observable regression only; no physics success claimed."""
    import hashlib
    cfg = module.load_config(ROOT, module.CONFIG)
    thesis_seed = module.initial_parameters(cfg, 'D_-2')
    low, high, width = module._trial_box(thesis_seed, 1.0)
    optimum = thesis_seed.copy()
    optimum[:5] += .45 * width[:5]

    class FakeEvaluator:
        model = 'D_-2'
        limit = 80
        attempts = 0

        def evaluate(self, x):
            if self.attempts >= self.limit:
                raise RuntimeError('Pilot evaluation budget reached')
            self.attempts += 1
            x = np.asarray(x, float)
            errors = (x[:5]-optimum[:5]) / 1e-5
            return {'status': 'COMPLETE_CONDITIONAL', 'candidate': x.tolist(),
                    'candidate_digest': hashlib.sha256(x.tobytes()).hexdigest(),
                    'residuals': errors.tolist(),
                    'chi2_approx': float(errors @ errors)}

    fake = FakeEvaluator()
    result = module.fit_pilot(fake, thesis_seed, sobol_points=8,
                              screen_points=18, active_params=5, starts=1,
                              local_max_nfev=10, perturbation=1.0, rng_seed=17)
    assert result['status'] == 'PILOT_BEST_FOUND'
    assert result['least_squares_runs']
    assert result['least_squares_runs'][0]['local_solver'].startswith('scipy.optimize.least_squares')
    assert any(x['phase'] == 'sobol' for x in result['candidate_history'])
    assert any(x['phase'] == 'least_squares' for x in result['candidate_history'])
    assert result['best']['chi2_approx'] < result['seed']['chi2_approx']
    assert fake.attempts <= fake.limit
