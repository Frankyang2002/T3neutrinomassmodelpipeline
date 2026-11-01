"""Per-model fixed-lambdaT3 Sobol + bounded local least-squares T3 fit.

Reuses the existing production numerical evaluator and fitting components.
Uses SAVED ordinary/restricted T3 RGE/C5 artifacts, not general-theory matching.
Never changes the frozen reference configurations or existing Sobol caches.

Preview:
  python -X utf8 -m scripts.RunThesisSobolModels --list
Single-model run:
  python -X utf8 -m scripts.RunThesisSobolModels --model A_-2 --run
All 16:
  python -X utf8 -m scripts.RunThesisSobolModels --run
Scalar alternative:
  python -X utf8 -m scripts.RunThesisSobolModels --lambda-t3 0.5 --run
"""
from __future__ import annotations

import argparse
import csv
from copy import deepcopy
import hashlib
import json
import math
from pathlib import Path
import re
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
TARGET = ROOT / 'configs/nufit_6_1_ic24_no.json'
REFERENCES = ROOT / 'Numerical/benchmarks/fitted_reference'
DEFAULT_SOURCE = ROOT / 'output/full/comparison'
REPS = {'A': (1, 3, 2), 'B': (2, 2, 1), 'C': (2, 2, 3),
        'D': (3, 1, 2), 'E': (3, 3, 2)}
ALPHAS = {'A': (-4, -2, 0), 'B': (-3, -1, 1), 'C': (-3, -1, 1),
          'D': (-2, 0, 2), 'E': (-4, -2, 0, 2)}
ALL_MODELS = tuple(f'{cls}_{alpha:+d}' for cls in REPS for alpha in ALPHAS[cls])
ARTIFACTS = ('uv_rgbeta_rge.json', 'eft1_rgbeta_rge.json', 'final_weinberg_coefficient.json')
OBS = ('delta_m21_sq_ev2', 'delta_m3l_sq_ev2', 'sin2_theta12', 'sin2_theta13', 'sin2_theta23')


def key_for(model: str) -> str:
    cls, sep, label = model.partition('_')
    if not sep or cls not in REPS:
        raise ValueError(f'Unknown T3 model {model!r}')
    try:
        alpha = int(label)
    except ValueError as exc:
        raise ValueError(f'Unknown T3 model {model!r}') from exc
    if alpha not in ALPHAS[cls]:
        raise ValueError(f'Unknown T3 model {model!r}')
    d1, d2, df = REPS[cls]
    return f'T3_dS1_{d1}_dS2_{d2}_dF_{df}_alpha_{"m" if alpha<0 else "p"}{abs(alpha)}'


def lambda_for(model: str, reference_lambda: float, requested: float | None,
               equal_class_e_multiplicity: bool) -> float:
    """Use repository's fitted scalar values unless a new common value is requested."""
    if requested is None:
        val = float(reference_lambda)
    else:
        val = float(requested)
        if equal_class_e_multiplicity and model.startswith('E_'):
            val *= 0.5
    if not math.isfinite(val) or val <= 0:
        raise ValueError('lambdaT3 must be finite and positive')
    return val


def scalar_slug(label: str, *, ref: bool) -> str:
    if ref:
        return 'reference_scalar'
    if not re.fullmatch(r'[0-9]+(?:\.[0-9]+)?(?:[eE][-+]?[0-9]+)?', label):
        raise ValueError('Invalid scalar label')
    return 'lambda_' + re.sub(r'[^A-Za-z0-9]', '_', label)


def _json(path: Path) -> dict:
    value = json.loads(path.read_text(encoding='utf-8-sig'))
    if not isinstance(value, dict):
        raise ValueError(f'Expected JSON object: {path}')
    return value


def build_search_payload(reference: dict, model: str, scalar: float, *, n_points: int,
                         target_chi2: float, seed: int, fit_lambda: bool = False) -> dict:
    """Rescale the existing seed and select fitted parameter directions.

    The square-root rescaling is a leading-order *initialisation heuristic*,
    not an equality of full RG-evolved observables.
    """
    if n_points < 8 or n_points & (n_points - 1):
        raise ValueError('n_points must be a power of two >= 8')
    if not math.isfinite(target_chi2) or target_chi2 <= 0:
        raise ValueError('target_chi2 must be finite and positive')
    p = deepcopy(reference)
    if p.get('kind') != 't3_pipeline_numerical_results_v1':
        raise ValueError('Not an integrated T3 numerical config')
    o = p['base_state']['ordinary']
    old = complex(float(o['lambdaT3']['real']), float(o['lambdaT3'].get('imag', 0.)))
    if abs(old.imag) > 1e-12 or old.real <= 0:
        raise ValueError('This real-input study requires a positive real reference lambdaT3')
    factor = math.sqrt(old.real / scalar)
    for name in ('y1', 'y2'):
        for part in ('real', 'imag'):
            if part in o[name]:
                o[name][part] = [[float(x) * factor for x in row] for row in o[name][part]]
    o['lambdaT3'] = {'real': scalar, 'imag': 0.0}
    params = []
    for name in ('y1', 'y2'):
        for i, row in enumerate(o[name]['real']):
            for j, item in enumerate(row):
                v = float(item)
                radius = max(2.0 * abs(v), 1.5e-4 * factor)
                params.append({'name': f'{name}_{i+1}{j+1}',
                               'binding': f'ordinary.{name}.real.{i}.{j}',
                               'low': v-radius, 'high': v+radius, 'scale': 'linear'})
    if fit_lambda:
        params.append({'name': 'lambdaT3_real', 'binding': 'ordinary.lambdaT3.real',
                       'low': 0.5*scalar, 'high': 2*scalar, 'scale': 'log'})
    p['ordering'] = 'NO'
    p['benchmark_search'] = {
        'enabled': True, 'method': 'sobol', 'use_current_model_representation': False,
        'n_points': int(n_points), 'target_chi2': float(target_chi2),
        'keep_best': 8, 'seed': int(seed),
        'oscillation_target': str(TARGET.resolve()),
        'local_refinement': {'enabled': True, 'seed_count': 3,
                             'active_counts': [5, 8, 12], 'max_nfev': [30, 40, 60],
                             'sensitivity_step': 0.02},
        'parameters': params,
    }
    p['sensitivity'] = {'enabled': True, 'gamma_min': 0.8, 'gamma_max': 1.2,
                        'n_points': 11, 'spacing': 'linear',
                        'oscillation_target': str(TARGET.resolve())}
    p['thesis_sobol'] = {
        'model': model, 'fixed_lambdaT3': not fit_lambda,
        'lambdaT3_seed': scalar, 'reference_lambdaT3': old.real,
        'yukawa_seed_rescale': factor,
        'note': 'Previous fitted-reference point used as centre. Rescaling is a seed heuristic. '
                'The complete numerical evaluator computes all fit scores.',
        'scope': 'Restricted ordinary T3 production matching, not full-general C5 matching.',
    }
    return p


def _signature(config: dict, files: list[Path], fitted_ref: Path) -> str:
    h = hashlib.sha256(json.dumps(config, sort_keys=True).encode())
    for path in [TARGET, fitted_ref, *files]:
        with path.open('rb') as fp:
            for chunk in iter(lambda: fp.read(1024*1024), b''):
                h.update(chunk)
    return h.hexdigest()


def _atomic(path: Path, obj: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.tmp')
    tmp.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    tmp.replace(path)


def fit_one(*, model: str, scenario: str, root: Path, output: Path, cfg_dir: Path,
            scalar_request: float | None, class_e_half: bool,
            n_points: int, target_chi2: float, seed: int, fit_lambda: bool,
            sensitivity: bool, force: bool, force_ls: bool) -> dict:
    # Delayed imports keep --list operational without installing Mathematica.
    from common.RunRecords import RunRecord
    from Numerical.orchestration.support.ResumePipelineNumericalResults import (
        discover_cases, parse_model_key)
    from Numerical.fitting import SobolBenchmarkSearch as sb
    from Numerical.fitting.OscillationFit import OscillationFitTarget
    from Numerical.orchestration.PipelineNumericalResults import _run_and_build_payload

    key = key_for(model)
    ref_path = REFERENCES / f'{key}.json'
    source = _json(ref_path)
    old_lam = float(source['base_state']['ordinary']['lambdaT3']['real'])
    lam = lambda_for(model, old_lam, scalar_request, class_e_half)
    config = build_search_payload(source, model, lam, n_points=n_points,
                                  target_chi2=target_chi2, seed=seed, fit_lambda=fit_lambda)
    cases = {k: c for sc, k, c in discover_cases(root) if sc == scenario}
    if key not in cases:
        raise FileNotFoundError(f'{key}: missing cached restricted matching at {root / scenario}. '
                                'Run the existing comparison production pipeline first.')
    case = cases[key]
    paths = [case/'data'/item for item in ARTIFACTS]
    for path in paths:
        if not path.is_file():
            raise FileNotFoundError(f'Missing saved matching/RGE artifact: {path}')
    sig = _signature(config, paths, ref_path)
    model_output = output / key
    cfg_file = cfg_dir / f'{key}.json'
    summary_file = model_output / 'fit_summary.json'
    diag_file = model_output / 'running_diagnostics.json'
    if summary_file.is_file() and diag_file.is_file() and not force:
        cached = _json(summary_file)
        if cached.get('input_signature') == sig and _json(diag_file).get('status') == 'Success':
            print(f'{model}: cached chi2={cached["chi2"]:.8g}', flush=True)
            return cached
    if (summary_file.exists() or diag_file.exists() or cfg_file.exists()) and not force:
        raise FileExistsError(f'{model}: previous incomplete or changed output at {model_output}. '
                              'Use a new --tag or --force to replace ONLY this study output.')
    if cfg_file.exists() and not force:
        raise FileExistsError(f'{cfg_file} already exists')

    d1, d2, df, alpha = parse_model_key(key)
    record = RunRecord(name=case.name.replace('_', '-', 1), alpha=alpha,
                       d_s1=d1, d_s2=d2, d_f=df, return_code=0, summary={}, output_dir=case)
    defs = sb._parameters_from_config(config)
    evaluator = sb._build_evaluator(payload=config, parameters=defs,
                                   uv_rgbeta_path=paths[0], intermediate_rgbeta_path=paths[1],
                                   final_weinberg_path=paths[2], target_path=TARGET)
    center = {}
    for item in config['benchmark_search']['parameters']:
        value = config['base_state']
        for part in item['binding'].split('.'):
            value = value[int(part)] if isinstance(value, list) else value[part]
        center[item['name']] = float(value)
    center_fit = evaluator(center)
    base_candidate = {'index': None, 'chi2': float(center_fit.chi2),
                      'parameters': center, 'ordering': str(center_fit.ordering),
                      'prediction': sb._prediction_dict(center_fit),
                      'whitened_residual': list(map(float, center_fit.whitened_residual)),
                      'source': 'fitted_reference_seed'}

    # Separate cache for every scalar/scenario. NEVER touch native global
    # output/benchmarks/<model>/sobol_benchmark.json.
    cache_path = model_output/'sobol_phase.json'
    existing = _json(cache_path) if cache_path.exists() else None
    if existing and existing.get('study_input_signature') == sig and not force:
        raw = existing['sobol']
    else:
        raw = sb._run_search(record=record, payload=config, definitions=defs,
                             evaluator=evaluator, cache_path=cache_path,
                             signature={'study_input_signature': sig}, n_points=n_points,
                             target_chi2=target_chi2, keep_best=8,
                             effective_seed=sb._effective_seed(seed, key))
        # Raw phase is intentionally separate from the refined output.
        _atomic(cache_path, {'study_input_signature': sig, 'sobol': raw})
    retained = list(raw['retained_best']) + [base_candidate]
    target = OscillationFitTarget.from_json(TARGET)
    local = config['benchmark_search']['local_refinement']
    best, refinement = sb.refine_sobol_candidates(
        evaluator=evaluator, observable_names=target.observable_names,
        definitions=defs, candidates=retained,
        target_chi2=target_chi2, seed_count=local['seed_count'],
        active_counts=local['active_counts'], max_nfev=local['max_nfev'],
        sensitivity_step=local['sensitivity_step'])

    # The native implementation deliberately skips LS when Sobol (or the
    # reference seed) already achieves target_chi2. For a thesis run explicitly
    # requesting Sobol + LS, try one bounded local confirmation even then.
    # Never worsen a fitted point or lose it if numerical derivatives fail.
    if force_ls and refinement.get('skipped', False):
        try:
            candidate, stage = sb._run_local_stage(
                evaluator=evaluator, observable_names=target.observable_names,
                parameters=sb._optimizer_parameters(defs), seed=best,
                active_count=min(5, len(defs)), max_nfev=20,
                sensitivity_step=local['sensitivity_step'],
                label='confirmation_least_squares')
            if float(candidate['chi2']) < float(best['chi2']):
                best = candidate
            refinement = {**refinement, 'skipped': False,
                          'stages': [stage], 'confirmation_only': True}
        except Exception as exc:
            refinement = {**refinement, 'least_squares_attempted': True,
                          'least_squares_error': f'{type(exc).__name__}: {exc}',
                          'confirmation_only': True}

    resolved = sb.apply_benchmark_parameters(config, best['parameters'], defs)
    resolved['benchmark_search']['enabled'] = False
    resolved['sensitivity']['enabled'] = bool(sensitivity)
    resolved['thesis_sobol']['best_chi2'] = float(best['chi2'])
    resolved['thesis_sobol']['method'] = 'Sobol plus reduced bounded least-squares (if needed)'
    resolved['thesis_sobol']['seed_used'] = best.get('source') or ('sobol_' + str(best.get('index')))
    # The report diagnostics are evaluated with this exact best input.
    _atomic(cfg_file, resolved)
    diagnostic = _run_and_build_payload(record=record, payload=resolved, config_path=cfg_file)
    diagnostic['thesis_sobol_fit'] = {
        'best_chi2': float(best['chi2']), 'seed_chi2': base_candidate['chi2'],
        'best': best, 'refinement': refinement, 'scalar_fixed': not fit_lambda,
        'scalar_input': lam, 'input_signature': sig,
    }
    _atomic(diag_file, diagnostic)
    summary = {
        'model': model, 'model_key': key, 'status': 'Success',
        'lambdaT3_requested': lam, 'lambdaT3_optimised': float(resolved['base_state']['ordinary']['lambdaT3']['real']),
        'lambdaT3_fixed': not fit_lambda, 'chi2': float(best['chi2']),
        'reference_seed_chi2': float(base_candidate['chi2']),
        'local_least_squares_ran': not bool(refinement.get('skipped', False)),
        'local_least_squares_attempted': bool(refinement.get('least_squares_attempted', False)) or not bool(refinement.get('skipped', False)),
        'n_sobol_evaluated': int(raw['evaluated_point_count']),
        'n_sobol_failed': int(raw['failed_point_count']),
        'input_signature': sig, 'source_artifacts': [str(p) for p in paths],
        'fitted_config': str(cfg_file), 'running_diagnostics': str(diag_file),
        'prediction': best['prediction'], 'parameters': best['parameters'],
    }
    _atomic(summary_file, summary)
    print(f'{model}: chi2={summary["chi2"]:.8g} scalar={lam:g} '
          f'LS={summary["local_least_squares_ran"]}', flush=True)
    return summary


def build_reports(output: Path, report: Path, summaries: list[dict]) -> None:
    if not summaries:
        return
    report.mkdir(parents=True, exist_ok=True)
    fields = ['model', 'lambdaT3_requested', 'lambdaT3_optimised', 'lambdaT3_fixed',
              'chi2', 'reference_seed_chi2', 'local_least_squares_ran', 'local_least_squares_attempted',
              'n_sobol_evaluated', 'n_sobol_failed', 'fitted_config']
    with (report/'fit_summary.csv').open('w', encoding='utf-8', newline='') as handle:
        writer=csv.DictWriter(handle, fieldnames=fields, extrasaction='ignore')
        writer.writeheader(); writer.writerows(summaries)
    from scripts.BuildThesisNumericHTML import write_outputs
    write_outputs(output, report)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--model', choices=ALL_MODELS, action='append', help='Repeat, otherwise all 16')
    p.add_argument('--list', action='store_true', help='Show the proposed benchmark list; no run')
    p.add_argument('--run', action='store_true', help='Actually run Sobol + least-squares')
    p.add_argument('--scenario', default='smallY_smallL',
                   help='Saved restricted matching artifact scenario')
    p.add_argument('--data-root', type=Path, default=DEFAULT_SOURCE)
    p.add_argument('--lambda-t3', type=float, default=None,
                   help='Common A-D lambdaT3; E uses half if --class-e-half')
    p.add_argument('--class-e-half', action=argparse.BooleanOptionalAction, default=True,
                   help='Use lambdaT3(E) = 0.5 * common value for charge multiplicity')
    p.add_argument('--fit-lambda', action='store_true',
                   help='Also vary lambdaT3 (default holds the scalar value fixed)')
    p.add_argument('--sobol-points', type=int, default=128)
    p.add_argument('--target-chi2', type=float, default=1e-6)
    p.add_argument('--seed', type=int, default=20261101)
    p.add_argument('--no-force-ls', action='store_true',
                   help='Allow native fitter to skip least squares when target already met')
    p.add_argument('--no-sensitivity', action='store_true',
                   help='Skip T1/T2/T3 gamma scans')
    p.add_argument('--tag', default=None,
                   help='Unique experiment label (default reference_scalar or lambda_<value>)')
    p.add_argument('--force', action='store_true',
                   help='Replace previous output for this experiment (does not change original refs)')
    a=p.parse_args(argv)
    selected = tuple(dict.fromkeys(a.model or ALL_MODELS))
    if a.lambda_t3 is not None and (not math.isfinite(a.lambda_t3) or a.lambda_t3 <= 0):
        p.error('--lambda-t3 must be finite and positive')
    if a.sobol_points < 8 or a.sobol_points & (a.sobol_points - 1):
        p.error('--sobol-points must be a power of two >= 8')
    if a.target_chi2 <= 0 or not math.isfinite(a.target_chi2):
        p.error('--target-chi2 must be finite and positive')
    tag = a.tag or scalar_slug(str(a.lambda_t3), ref=a.lambda_t3 is None)
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,60}', tag):
        p.error('--tag must use letters, numbers, underscore or hyphen')
    data_out = ROOT / 'output/thesis_sobol' / tag
    cfg_dir = ROOT / 'configs/thesis_sobol' / tag
    report = ROOT / 'Reports/output/thesis_sobol' / tag
    print(f'Study: {tag} | models: {len(selected)} | Sobol points: {a.sobol_points}', flush=True)
    for model in selected:
        path = REFERENCES / f'{key_for(model)}.json'
        if not path.is_file():
            print(f'  {model}: MISSING fitted-reference seed at {path}', flush=True)
            continue
        cfg = _json(path)
        ref = float(cfg['base_state']['ordinary']['lambdaT3']['real'])
        lam = lambda_for(model, ref, a.lambda_t3, a.class_e_half)
        print(f'  {model}: fixed scalar lambdaT3={lam:g} (reference={ref:g})', flush=True)
    if a.list or not a.run:
        print('Dry run only. Add --run to optimise and generate HTML reports.')
        return 0
    summaries=[]
    failures=[]
    for model in selected:
        try:
            summaries.append(fit_one(model=model, scenario=a.scenario, root=a.data_root,
                                     output=data_out, cfg_dir=cfg_dir, scalar_request=a.lambda_t3,
                                     class_e_half=a.class_e_half, n_points=a.sobol_points,
                                     target_chi2=a.target_chi2, seed=a.seed, fit_lambda=a.fit_lambda,
                                     sensitivity=not a.no_sensitivity, force=a.force,
                                     force_ls=not a.no_force_ls))
        except Exception as exc:
            msg = f'{type(exc).__name__}: {exc}'
            failures.append({'model': model, 'error': msg})
            print(f'  FAILED {model}: {msg}', flush=True)
    if summaries:
        try:
            build_reports(data_out, report, summaries)
            print('HTML reports:', report, flush=True)
        except Exception as exc:
            failures.append({'model': '<reports>', 'error': f'{type(exc).__name__}: {exc}'})
            print('HTML generation failed:', exc, flush=True)
    if failures:
        _atomic(report/'failures.json', failures)
    print(f'Complete: {len(summaries)}/{len(selected)} fits, failures={len(failures)}')
    return 1 if failures else 0


if __name__ == '__main__':
    raise SystemExit(main())
