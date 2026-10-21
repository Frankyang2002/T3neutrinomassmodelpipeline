"""Refine saved lambda=5e-9 Yukawa fits without rerunning Sobol.

Uses the authoritative fitted configuration and numerical evaluator, together with
sensitivity-ranked reduced least squares from OscillationOptimizer. Outputs are
isolated until --promote is requested. No physics implementation is modified.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path

import numpy as np

from Numerical.orchestration import RunSmallLambdaRefit as base
from Numerical.fitting.SobolBenchmarkSearch import BenchmarkParameter, _build_evaluator
from Numerical.fitting.OscillationOptimizer import (
    EvaluationRecorder, OptimizerParameter, estimate_local_sensitivities,
    optimise_reduced_least_squares, unit_vector_from_parameters,
)

OUT = base.ROOT / 'output' / 'studies' / 'lambda_5e-9_refit'
RESULTS = OUT / 'targeted_refinement'
FIELDS = ('delta_m21_sq_ev2', 'delta_m3l_sq_ev2', 'sin2_theta12', 'sin2_theta13', 'sin2_theta23')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    base.save(path, value)


def vector(cfg):
    return np.concatenate([base.matr(cfg['base_state']['ordinary'][k]).real.ravel() for k in ('y1', 'y2')])


def with_vector(cfg, x):
    payload = copy.deepcopy(cfg)
    for j, name in enumerate(('y1', 'y2')):
        spec = payload['base_state']['ordinary'][name]
        base.setmat(spec, np.asarray(x[9*j:9*j+9]).reshape(3, 3) + 1j * base.matr(spec).imag)
    return payload


def refine(row, args, mapping):
    model = row['model']
    fit_file = OUT / 'fits' / f'{model}.json'
    cfg_file = OUT / 'fitted' / f'{model}.json'
    if not fit_file.is_file() or not cfg_file.is_file():
        print(f'SKIP {model}: no saved fit and fitted config', flush=True)
        return
    prior = base.load(fit_file)
    score = float(prior.get('chi2_fit', np.inf))
    if not np.isfinite(score):
        print(f'SKIP {model}: no finite chi2', flush=True)
        return
    if score <= args.target:
        print(f'SKIP {model}: already chi2={score:.8g} <= {args.target:g}', flush=True)
        return
    cfg = base.load(cfg_file)
    expected = 2.5e-9 if base.cls(row) == 'E' else 5e-9
    actual = float(cfg['base_state']['ordinary']['lambdaT3']['real'])
    if not np.isclose(actual, expected, rtol=1e-10, atol=0):
        raise ValueError(f'{model}: lambdaT3 mismatch {actual} != {expected}')
    original = vector(cfg)
    if np.any(np.abs(original) >= args.bound):
        raise ValueError(f'{model}: existing Yukawas exceed --bound {args.bound}.')
    artifacts = base.locate_artifacts(model, mapping)
    names = [f'{y}_{i}{j}' for y in ('y1', 'y2') for i in range(3) for j in range(3)]
    defs = [BenchmarkParameter(n, f'ordinary.{y}.real.{i}.{j}', -args.bound, args.bound)
            for n, y, i, j in ((f'{y}_{i}{j}', y, i, j) for y in ('y1', 'y2') for i in range(3) for j in range(3))]
    evaluator = _build_evaluator(payload=cfg, parameters=defs,
        uv_rgbeta_path=artifacts[0], intermediate_rgbeta_path=artifacts[1],
        final_weinberg_path=artifacts[2], target_path=base.TARGET)
    params = tuple(OptimizerParameter(name=p.name, low=p.low, high=p.high, scale=p.scale).validated() for p in defs)
    config_sha = sha(cfg_file)
    manifest = {'model': model, 'original_chi2': score, 'original_fitted_sha256': config_sha,
                'artifacts_sha256': [sha(p) for p in artifacts], 'target_sha256': sha(base.TARGET),
                'bound': args.bound, 'target_chi2': args.target, 'budgets': [args.nfev12, args.nfev18],
                'method': 'saved_fit_sensitivity_ranked_12_then_18'}
    dest = RESULTS / model
    state_path = dest / 'result.json'
    if args.resume and state_path.is_file():
        cached = base.load(state_path)
        if cached.get('manifest') == manifest and cached.get('status') == 'Complete':
            print(f'REUSE {model}: chi2={cached["best_chi2"]:.8g}', flush=True)
            if args.promote:
                promote(model, fit_file, cfg_file, cached)
            return

    if not args.original_c5:
        from Numerical.orchestration.FastLambda5e9C5 import install
        install()
    seed_dict = dict(zip(names, original.tolist()))
    fit = evaluator(seed_dict)
    verified_start = float(fit.chi2)
    if not np.isclose(verified_start, score, atol=1e-5, rtol=1e-4):
        raise RuntimeError(f'{model}: stored chi2={score:.8g}, evaluator now={verified_start:.8g}; check numerical evaluator before refinement')
    best_score = verified_start
    best_x = original.copy()
    stage_log = []
    for count, budget in ((12, args.nfev12), (18, args.nfev18)):
        if best_score <= args.target:
            break
        recorder = EvaluationRecorder(evaluator=evaluator, parameters=params,
                                      observable_names=FIELDS, checkpoint_path=dest / f'progress_{count}.json')
        seed_unit = unit_vector_from_parameters(params, dict(zip(names, best_x.tolist())))
        baseline, sensitivity = estimate_local_sensitivities(recorder, seed_unit, step=0.02)
        result, active = optimise_reduced_least_squares(recorder, seed_unit, sensitivity,
                                                         active_count=count, max_nfev=budget)
        best = recorder.best
        if best is None:
            raise RuntimeError(f'{model}: no successful evaluation in {count}-parameter stage')
        trial = float(best['chi2'])
        if trial < best_score:
            best_score = trial
            best_x = np.asarray([best['parameters'][n] for n in names], dtype=float)
        stage_log.append({'active_count': count, 'max_nfev': budget, 'nfev': int(result.nfev),
                          'optimizer_success': bool(result.success), 'active_parameters': list(active),
                          'stage_best_chi2': trial, 'overall_best_chi2': best_score,
                          'numerical_evaluations': len(recorder.history)})
        print(f'  {model}: active={count} stage chi2={trial:.8g} best={best_score:.8g}', flush=True)
        save(dest / 'stage_summary.json', {'manifest': manifest, 'stages': stage_log, 'best_chi2': best_score})
    best_fit = evaluator(dict(zip(names, best_x.tolist())))
    checked = float(best_fit.chi2)
    if not np.isclose(best_score, checked, rtol=1e-6, atol=1e-6):
        raise RuntimeError(f'{model}: re-evaluation mismatch {best_score} vs {checked}')
    improved = checked < score
    result = {'status': 'Complete', 'manifest': manifest, 'best_chi2': checked,
              'improved': improved, 'passes_strict_5p86': checked <= 5.86,
              'reached_target': checked <= args.target, 'stages': stage_log,
              'parameters': dict(zip(names, best_x.tolist())),
              'prediction': dict(zip(best_fit.observable_names, map(float, best_fit.prediction)))}
    save(state_path, result)
    save(dest / 'fitted_candidate.json', with_vector(cfg, best_x))
    print(f'DONE {model}: {score:.8g} -> {checked:.8g}; target={args.target:g} achieved={checked<=args.target}', flush=True)
    if args.promote:
        promote(model, fit_file, cfg_file, result)


def promote(model, fit_file, cfg_file, result):
    if not result.get('improved'):
        print(f'NOT PROMOTED {model}: no improvement', flush=True)
        return
    dest = RESULTS / model
    orig = base.load(fit_file)
    current_cfg = base.load(cfg_file)
    if sha(cfg_file) != result['manifest']['original_fitted_sha256']:
        raise RuntimeError(f'{model}: original fitted config changed; refusing to promote stale result')
    archived = dest / 'before_promotion'
    archived.mkdir(parents=True, exist_ok=True)
    save(archived / 'fit.json', orig)
    save(archived / 'fitted.json', current_cfg)
    candidate = base.load(dest / 'fitted_candidate.json')
    candidate.setdefault('lambda_5e9_study', {})['refinement'] = {'from_chi2': orig['chi2_fit'],
                                                     'new_chi2': result['best_chi2'], 'target': result['manifest']['target_chi2']}
    newer = dict(orig)
    newer.update({'chi2_fit': result['best_chi2'],
                  'status': 'Success' if result['best_chi2'] <= 5.86 else 'PoorFit',
                  'passes_screening': result['best_chi2'] <= 5.86,
                  'fit_method': 'saved_fit_sensitivity_ranked_refinement',
                  'prediction': result['prediction'], 'refinement_result': str(dest / 'result.json')})
    save(cfg_file, candidate)
    save(fit_file, newer)
    print(f'PROMOTED {model}: chi2={result["best_chi2"]:.8g}; old files archived under {archived}', flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--model', default='T3_dS1_2_dS2_2_dF_3_alpha_m3', help='Exact model identifier; use --all-poor for every model above target')
    p.add_argument('--all-poor', action='store_true')
    p.add_argument('--target', type=float, default=1.0)
    p.add_argument('--bound', type=float, default=4.0)
    p.add_argument('--nfev12', type=int, default=150)
    p.add_argument('--nfev18', type=int, default=80)
    p.add_argument('--artifact-map', type=Path)
    p.add_argument('--resume', action='store_true')
    p.add_argument('--promote', action='store_true', help='Replace improved fitted and fit JSONs after archiving both')
    p.add_argument('--original-c5', action='store_true', help='Use original symbolic evaluator; default fast numerical C5')
    a = p.parse_args()
    if a.target <= 0 or a.bound <= 0 or a.nfev12 <= 0 or a.nfev18 <= 0:
        p.error('Positive --target, --bound and evaluation budgets required')
    mapping = base.load(a.artifact_map) if a.artifact_map else None
    selected = [r for r in base.original_rows() if a.all_poor or r['model'] == a.model]
    if not selected:
        p.error('Unknown model')
    for row in selected:
        refine(row, a, mapping)


if __name__ == '__main__':
    main()
