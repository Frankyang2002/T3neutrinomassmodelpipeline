"""Continue only the unsuccessful 16-model GENERAL T3 Yukawa fits.

Read the existing general-RGE / general-Matchete C5 campaign, begin from each
saved best point, and deepen local optimisation from 12 to 18 Yukawa entries.
B(-1)/C(-1) use their existing conditional mixed-scalar evaluator; this script
never falls back to a restricted-Lagrangian matching kernel or refits accepted
models. The one-scale results are NOT a certified sequential EFT calculation.

PowerShell (project root):
    python -X utf8 -m scripts.ContinueT3GeneralPoorFits --self-test
    python -X utf8 -m scripts.ContinueT3GeneralPoorFits --run-id continue_general16_v1 --preflight
    python -u -X utf8 -m scripts.ContinueT3GeneralPoorFits --run-id continue_general16_v1
    python -u -X utf8 -m scripts.ContinueT3GeneralPoorFits --run-id continue_general16_v1 --resume
"""
from __future__ import annotations

import argparse
import csv
from copy import deepcopy
from hashlib import sha256
import json
import math
from pathlib import Path
import traceback
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CAMPAIGNS = (ROOT / 'output' / 'studies' / 'general_scalar_extended'
             / 'thesis_general_all16')
DEFAULT_SOURCE = 'thesis16_general_c5_full_v1'
DEFAULT_MODELS = ('A_-2', 'B_+1', 'C_-1', 'D_+0', 'D_+2')
SCREEN = 5.86
FORMAT = 't3-general-fit-continuation-18-yukawas-conditional-v1'
MIXED = {'B_-1', 'C_-1'}
NAMES = tuple(f'{name}_{i}{j}' for name in ('y1', 'y2')
              for i in range(3) for j in range(3))


def tag(model: str) -> str:
    return model.replace('-', 'm').replace('+', 'p')


def read(path: Path) -> dict:
    value = json.loads(path.read_text(encoding='utf-8-sig'))
    if not isinstance(value, dict):
        raise ValueError(f'Expected JSON object: {path}')
    return value


def save(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + '.tmp')
    temp.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    temp.replace(path)


def sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def model_parameters(config: dict) -> dict[str, float]:
    raw = config['base_state']['ordinary']
    out = {}
    for name in ('y1', 'y2'):
        real = np.asarray(raw[name]['real'], dtype=float)
        imag = np.asarray(raw[name]['imag'], dtype=float)
        if real.shape != (3, 3) or imag.shape != (3, 3):
            raise ValueError('Expected 3x3 Yukawa matrix')
        if not np.isfinite(real).all() or not np.isfinite(imag).all() or np.max(abs(imag)) > 1e-13:
            raise ValueError('Continuation currently supports only 18 REAL Yukawa entries')
        for i in range(3):
            for j in range(3):
                out[f'{name}_{i}{j}'] = float(real[i, j])
    return out


def parameters_from_result(result: dict) -> dict[str, float]:
    values = result.get('best_parameters')
    if not isinstance(values, dict) or set(values) != set(NAMES):
        raise ValueError('Saved result lacks the 18 fitted Yukawa parameters')
    parsed = {name: float(values[name]) for name in NAMES}
    if not all(math.isfinite(v) for v in parsed.values()):
        raise ValueError('Nonfinite fitted Yukawa')
    return parsed


def perturb(point: dict[str, float], rng: np.random.Generator, relative: float,
            absolute_floor: float) -> dict[str, float]:
    return {key: float(value + rng.normal() * relative * max(abs(value), absolute_floor))
            for key, value in point.items()}


def local_definitions(center: dict[str, float], *, window_factor: float,
                      window_floor: float, absolute_bound: float):
    from Numerical.fitting.SobolBenchmarkSearch import BenchmarkParameter
    defs = []
    for name in NAMES:
        val = center[name]
        half = max(abs(val) * window_factor, window_floor)
        lo, hi = max(-absolute_bound, val - half), min(absolute_bound, val + half)
        if not lo < val < hi:
            raise ValueError(f'No optimisation room for {name}: {lo}, {val}, {hi}')
        defs.append(BenchmarkParameter(
            name, f'ordinary.{name.split("_")[0]}.real.{name[-2]}.{name[-1]}',
            lo, hi))
    return defs


def candidate(fit: Any, point: dict[str, float], label: str) -> dict:
    return {'source': label, 'chi2': float(fit.chi2),
            'parameters': dict(point),
            'ordering': str(getattr(fit, 'ordering', 'NO')),
            'prediction': dict(zip(fit.observable_names, map(float, fit.prediction))),
            'whitened_residual': np.asarray(fit.whitened_residual, dtype=float).tolist()}


def source_info(model: str, source_root: Path) -> tuple[dict, dict | None]:
    summary = read(source_root / 'summary.json')
    rows = {x['model']: x for x in summary['rows']}
    if model not in rows:
        raise ValueError(f'{model} missing in original study')
    row = rows[model]
    path = source_root / tag(model) / 'result.json'
    old = read(path) if path.is_file() else None
    return row, old


def check_poor(row: dict) -> None:
    chi = row.get('chi2_final')
    if chi is not None and float(chi) <= SCREEN:
        raise ValueError('Already accepted: refusing to refit a successful model')
    if row['status'] not in ('ExploratoryPoorFit', 'Blocked'):
        raise ValueError('Unexpected original fit status: ' + str(row['status']))


def prepare_nonmixed(model: str, args: argparse.Namespace, source_settings: dict,
                     target: Any):
    from scripts import RunT3GeneralAll16FromReference as general
    uv, eft, provenance, kernel, coefficient, matching_kind = general.nonmixed_inputs(
        model, source_settings['canonical_run_id'],
        generate_general_c5=False, wolframscript=args.wolframscript)
    reference_path = general.ref_path(model)
    reference = read(reference_path)
    mu2 = (args.d0_mu2_fraction if model == 'D_+0'
           else float(source_settings['mu2_fraction']))
    fixed = general.fixed_overrides(
        reference, uv,
        quartic_fraction=float(source_settings['quartic_fraction']),
        mu2_fraction=mu2,
        extra_yukawa_fraction=float(source_settings['extra_yukawa_fraction']))
    evaluator = general.NonMixedEvaluator(model, uv, eft, kernel, fixed, target)

    def output_config(parameters: dict):
        cfg = deepcopy(reference)
        for name in ('y1', 'y2'):
            arr = np.array([[parameters[f'{name}_{i}{j}'] for j in range(3)]
                            for i in range(3)], dtype=float)
            cfg['base_state']['ordinary'][name] = general.matrix_json(arr)
        cfg['benchmark_search'] = {'enabled': False}
        cfg['sensitivity'] = {'enabled': False}
        return cfg

    def output_overrides(parameters: dict):
        return general.override_yukawas(fixed, parameters)

    return evaluator, output_config, output_overrides, {
        'method': matching_kind,
        'general_symbolic_provenance': provenance,
        'source_reference': str(reference_path.relative_to(ROOT)),
        'source_reference_sha256': sha(reference_path),
        'general_c5_coefficient': str(coefficient.relative_to(ROOT)),
        'general_c5_coefficient_sha256': sha(coefficient),
        'quartic_fraction': float(source_settings['quartic_fraction']),
        'mu2_fraction': mu2,
        'extra_yukawa_fraction': float(source_settings['extra_yukawa_fraction']),
        'modified_scalar_boundary': (model == 'D_+0' and
                                     mu2 != float(source_settings['mu2_fraction'])),
    }


def prepare_mixed(model: str, args: argparse.Namespace, source_settings: dict,
                  target: Any):
    from scripts import RunT3BCGeneralFromReference as bc_ref
    from scripts import RunT3GeneralScalarBCExploratoryFit as bc
    from scripts.RunGeneralT3NumericalC5 import DEFAULT_SOURCE
    source = bc_ref.build_reference_source(model)
    bc_ref.ensure_betas(model, source, export=False, wolframscript=args.wolframscript)
    bundle = bc.beta_bundle(model, source, source['registry'])
    coeff = (DEFAULT_SOURCE / ('B_m1_final' if model == 'B_-1' else 'C_m1_final'))
    coeff = coeff.with_name(coeff.name + '_one_loop_probe_weinberg_ordered_coefficient.inputform.txt')
    if not coeff.is_file():
        raise FileNotFoundError(str(coeff))
    compiled = bc.compile_source(coeff)
    evaluator = bc.CandidateEvaluator(
        'reference_scalar_general', model, source, source['registry'], bundle,
        compiled, target,
        quartic_fraction=float(source_settings['quartic_fraction']),
        mu2_fraction=float(source_settings['mu2_fraction']))

    def output_config(parameters: dict):
        arr = np.array([parameters[name] for name in NAMES], dtype=float)
        cfg = bc.trial_record(source, arr)['cfg']
        cfg['benchmark_search'] = {'enabled': False}
        cfg['sensitivity'] = {'enabled': False}
        return cfg

    return evaluator, output_config, None, {
        'method': 'mixed-scalar-conditional-one-scale',
        'source_reference': str(source['cfg_path'].relative_to(ROOT)),
        'source_reference_sha256': sha(source['cfg_path']),
        'general_c5_coefficient': str(coeff.relative_to(ROOT)),
        'general_c5_coefficient_sha256': sha(coeff),
        'quartic_fraction': float(source_settings['quartic_fraction']),
        'mu2_fraction': float(source_settings['mu2_fraction']),
        'extra_yukawa_fraction': 0.0,
        'modified_scalar_boundary': False,
    }


def c5_payload(details: dict, mixed: bool) -> dict:
    # Both evaluators return compact complex {real:[[...]],imag:[[...]]} matrices.
    keys = (('C5_muS_gev_inv', 'C5_low_gev_inv') if mixed else
            ('c5_at_scalar_threshold', 'c5_at_low_scale'))
    result = {}
    for dest, name in zip(('C5_muS_gev_inv', 'C5_low_gev_inv'), keys):
        raw = details[name]
        m = np.asarray(raw['real'], dtype=float) + 1j*np.asarray(raw['imag'], dtype=float)
        if m.shape != (3, 3) or not np.isfinite(m).all() or not np.allclose(m, m.T, rtol=1e-8, atol=1e-18):
            raise ValueError('Invalid conditional symmetric C5 matrix')
        result[dest] = raw
    return result


def run_one(model: str, args: argparse.Namespace, source_root: Path,
            source_settings: dict, dest: Path) -> dict:
    from Numerical.fitting.OscillationFit import OscillationFitTarget
    from Numerical.orchestration.RunSmallLambdaRefit import TARGET
    from Numerical.fitting.SobolBenchmarkSearch import refine_sobol_candidates

    row, old = source_info(model, source_root)
    check_poor(row)
    if old is None and model != 'D_+0':
        raise FileNotFoundError(f'Saved best point not available for {model}')
    target = OscillationFitTarget.from_json(TARGET)
    if model in MIXED:
        evaluator, config_builder, overrides_builder, meta = prepare_mixed(model, args, source_settings, target)
    else:
        evaluator, config_builder, overrides_builder, meta = prepare_nonmixed(model, args, source_settings, target)

    centers = []
    if old is not None:
        centers.append(('previous_best', parameters_from_result(old)))
    else:
        # No prior D(0) best: use its own general-reference seed, and a saved
        # general D(-2) flavour point as a distinct seed with the same 3x3 shape.
        from scripts import RunT3GeneralAll16FromReference as general
        centers.append(('D0_general_reference_seed', model_parameters(read(general.ref_path(model)))))
        alternate = source_root / 'D_m2' / 'result.json'
        if alternate.is_file():
            centers.append(('Dminus2_general_fitted_flavour_seed',
                            parameters_from_result(read(alternate))))

    original_score = float(row['chi2_final']) if row.get('chi2_final') is not None else None
    best = None
    details_log = []
    attempts = []
    rng = np.random.default_rng(args.seed + list(DEFAULT_MODELS).index(model))
    for label, center in centers:
        try:
            fit = evaluator(center)
            rec = candidate(fit, center, label)
            attempts.append(rec)
            print(f'  {label}: chi2={rec["chi2"]:.8g}', flush=True)
        except Exception as exc:
            details_log.append({'seed': label, 'status': 'SeedFailed',
                                'message': f'{type(exc).__name__}: {exc}'})
            print(f'  {label}: {type(exc).__name__}: {exc}', flush=True)

    if not attempts:
        raise RuntimeError('No valid general-matching seed: ' + json.dumps(details_log))
    best = min(attempts, key=lambda x: x['chi2'])
    for restart in range(args.restarts):
        if best['chi2'] <= args.target_chi2:
            break
        center = dict(best['parameters'])
        if restart:
            trial = perturb(center, rng, relative=args.perturbation,
                            absolute_floor=args.window_floor)
            try:
                trial_fit = evaluator(trial)
                trial_rec = candidate(trial_fit, trial, f'restart_{restart}_seed')
                if trial_rec['chi2'] < best['chi2']:
                    best = trial_rec
                    center = dict(trial)
            except Exception as exc:
                details_log.append({'restart': restart, 'status': 'PerturbationFailed',
                                    'message': f'{type(exc).__name__}: {exc}'})
        defs = local_definitions(center, window_factor=args.window_factor,
                                 window_floor=args.window_floor,
                                 absolute_bound=args.absolute_bound)
        try:
            candidate_best, history = refine_sobol_candidates(
                evaluator=evaluator, observable_names=target.observable_names,
                definitions=defs, candidates=[best],
                target_chi2=args.target_chi2, seed_count=1,
                active_counts=(12, 18),
                max_nfev=(args.stage12_nfev, args.stage18_nfev),
                sensitivity_step=args.sensitivity_step)
            details_log.append({'restart': restart, 'status': 'Completed',
                                'history': history})
            if candidate_best['chi2'] < best['chi2']:
                best = candidate_best
        except Exception as exc:
            details_log.append({'restart': restart, 'status': 'RefinementFailedRetainedBest',
                                'message': f'{type(exc).__name__}: {exc}'})
        print(f'  After restart {restart+1}/{args.restarts}: best chi2={best["chi2"]:.8g}', flush=True)

    final_fit, details = evaluator.evaluate(best['parameters'], detail=True)
    final_score = float(final_fit.chi2)
    # Preserve a previously better point; do not silently degrade the study.
    if original_score is not None and final_score >= original_score * (1 - 1e-9):
        status = 'NoImprovementRetainOriginal'
        use_for_merge = False
    else:
        status = 'ExploratoryAcceptedFit' if final_score <= SCREEN else 'ExploratoryPoorFit'
        use_for_merge = True
    cfg = config_builder(best['parameters'])
    cfg['general_fit_continuation'] = {
        'source_campaign': args.source_run_id, 'source_model': model,
        'chi2': final_score, 'fit_screen': SCREEN,
        'full_18_yukawas_allowed': True,
        'conditional_one_scale_c5': True,
        'complete_sequential_eft_certified': False,
        'source_scalar_boundary_modified': meta['modified_scalar_boundary']}
    out = {'format': FORMAT, 'model': model, 'status': status,
           'source_status': row['status'], 'source_chi2': original_score,
           'chi2_final': final_score, 'prediction': dict(zip(
               final_fit.observable_names, map(float, final_fit.prediction))),
           'whitened_residual': np.asarray(final_fit.whitened_residual).tolist(),
           'best_parameters': best['parameters'], 'best_source': best['source'],
           'local_refinement_continuation': details_log,
           'matching_and_rge_provenance': meta,
           'c5': c5_payload(details, mixed=model in MIXED),
           'numerical_details': details,
           'original_result_sha256': sha(source_root / tag(model) / 'result.json') if old else None,
           'fit_target_sha256': sha(Path(TARGET)),
           'use_for_merged_summary': use_for_merge,
           'physical_general_c5_certified': False,
           'full_general_sequential_matching_performed': False}
    dest.mkdir(parents=True, exist_ok=True)
    save(dest / 'result.json', out)
    save(dest / 'fitted_config.json', cfg)
    if overrides_builder is not None:
        save(dest / 'general_uv_overrides.json', overrides_builder(best['parameters']))
    save(dest / 'c5_conditional.json', out['c5'])
    return {'model': model, 'status': status, 'chi2_previous': original_score,
            'chi2_final': final_score, 'use_for_merged_summary': use_for_merge,
            'scalar_boundary_modified': meta['modified_scalar_boundary'],
            'result_path': str((dest / 'result.json').relative_to(ROOT))}


def merged_rows(source_root: Path, continuation: list[dict], dest: Path) -> None:
    originals = read(source_root / 'summary.json')['rows']
    newer = {x['model']: x for x in continuation}
    rows = []
    for entry in originals:
        model = entry['model']
        next_row = newer.get(model)
        updated = (next_row is not None and next_row.get('use_for_merged_summary') is True
                   and not next_row.get('scalar_boundary_modified', False))
        if updated:
            fitted = read(dest / tag(model) / 'result.json')
            rows.append({'model': model, 'status': fitted['status'],
                         'chi2_final': fitted['chi2_final'],
                         'prediction': fitted['prediction'],
                         'source': 'GeneralFitContinuation',
                         'result_path': next_row['result_path']})
        else:
            rows.append({'model': model, 'status': entry['status'],
                         'chi2_final': entry.get('chi2_final'),
                         'prediction': entry.get('prediction'),
                         'source': 'OriginalGeneralCampaign',
                         'result_path': str((source_root / tag(model) / 'result.json').relative_to(ROOT))
                         if (source_root / tag(model) / 'result.json').is_file() else None})
    save(dest / 'merged_all16.json', {'format': FORMAT, 'rows': rows,
                                     'accepted_count': sum(x['chi2_final'] is not None and
                                                           float(x['chi2_final']) <= SCREEN for x in rows),
                                     'blocked_count': sum(x['status'] == 'Blocked' for x in rows)})
    with (dest / 'merged_all16.csv').open('w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(['model', 'status', 'chi2_final', 'source', 'result_path',
                    'delta_m21_sq_ev2', 'delta_m3l_sq_ev2', 'sin2_theta12',
                    'sin2_theta13', 'sin2_theta23'])
        for x in rows:
            pred = x.get('prediction') or {}
            w.writerow([x['model'], x['status'], x['chi2_final'], x['source'],
                        x['result_path'], *(pred.get(k, '') for k in
                        ('delta_m21_sq_ev2', 'delta_m3l_sq_ev2', 'sin2_theta12',
                         'sin2_theta13', 'sin2_theta23'))])


def write_summary(dest: Path, rows: list[dict], source_root: Path) -> None:
    save(dest / 'summary.json', {'format': FORMAT, 'rows': rows})
    with (dest / 'summary.csv').open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=('model', 'status', 'chi2_previous',
                                          'chi2_final', 'use_for_merged_summary',
                                          'scalar_boundary_modified', 'result_path', 'message'),
                           extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)
    merged_rows(source_root, rows, dest)


def self_test() -> None:
    import tempfile
    m = {'base_state': {'ordinary': {
        'y1': {'real': [[(3*i+j)*1e-5 for j in range(3)] for i in range(3)],
               'imag': [[0]*3 for _ in range(3)]},
        'y2': {'real': [[(9+3*i+j)*1e-5 for j in range(3)] for i in range(3)],
               'imag': [[0]*3 for _ in range(3)]}}}}
    p = model_parameters(m)
    assert len(p) == 18 and p['y1_00'] == 0 and np.isclose(p['y2_22'], 17e-5)
    assert parameters_from_result({'best_parameters': p}) == p
    assert tag('A_-2') == 'A_m2' and tag('B_+1') == 'B_p1'
    a = perturb(p, np.random.default_rng(5), relative=0.1, absolute_floor=1e-5)
    assert any(a[k] != p[k] for k in p)
    assert set(DEFAULT_MODELS) == {'A_-2', 'B_+1', 'C_-1', 'D_+0', 'D_+2'}
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / 'x.json'
        save(path, m)
        assert read(path) == m
    print('Self-test passed: saved-parameter parsing, model tagging, seed variation, JSON IO')


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-id', default='continue_general16_v1')
    p.add_argument('--source-run-id', default=DEFAULT_SOURCE)
    p.add_argument('--models', default=','.join(DEFAULT_MODELS))
    p.add_argument('--preflight', action='store_true')
    p.add_argument('--self-test', action='store_true')
    p.add_argument('--resume', action='store_true')
    p.add_argument('--restarts', type=int, default=1)
    p.add_argument('--stage12-nfev', type=int, default=60)
    p.add_argument('--stage18-nfev', type=int, default=180)
    p.add_argument('--target-chi2', type=float, default=SCREEN)
    p.add_argument('--window-factor', type=float, default=2.0)
    p.add_argument('--window-floor', type=float, default=2e-4)
    p.add_argument('--absolute-bound', type=float, default=3.0)
    p.add_argument('--perturbation', type=float, default=0.15)
    p.add_argument('--sensitivity-step', type=float, default=0.02)
    p.add_argument('--seed', type=int, default=20261103)
    p.add_argument('--d0-mu2-fraction', type=float, default=0.01,
                   help='D(+0) extra bilinear UV boundary fraction. Values != 0.01 '
                        'change that model benchmark and are excluded from merged all16.')
    p.add_argument('--wolframscript', default='wolframscript')
    args = p.parse_args(argv)
    if args.self_test:
        self_test()
        return 0
    if (not args.run_id or '/' in args.run_id or '\\' in args.run_id or args.run_id in ('.', '..') or
        '/' in args.source_run_id or '\\' in args.source_run_id or args.source_run_id in ('.', '..')):
        p.error('run ID must be a simple directory name')
    if args.source_run_id == args.run_id:
        p.error('source and destination run IDs must differ')
    if (args.restarts < 1 or args.stage12_nfev < 2 or args.stage18_nfev < 2 or
        args.target_chi2 <= 0 or args.window_factor <= 0 or args.window_floor <= 0 or
        args.absolute_bound <= 0 or args.perturbation < 0 or not 0 < args.sensitivity_step <= 0.25 or
        not 0 <= args.d0_mu2_fraction <= 0.2):
        p.error('invalid optimisation settings')
    selected = [m.strip() for m in args.models.split(',') if m.strip()]
    if not selected or len(selected) != len(set(selected)) or set(selected) - set(DEFAULT_MODELS):
        p.error('Select only originally unsuccessful models: ' + ','.join(DEFAULT_MODELS))
    source_root = CAMPAIGNS / args.source_run_id
    dest = CAMPAIGNS / args.run_id
    if not (source_root / 'summary.json').is_file():
        p.error(f'Missing general fit campaign: {source_root}')
    source_settings = read(source_root / 'settings.json')
    rows0 = {r['model']: r for r in read(source_root / 'summary.json')['rows']}
    for m in selected:
        check_poor(rows0[m])
    sig = {'format': FORMAT, 'selected': selected,
           'source_run_id': args.source_run_id,
           'source_summary_sha256': sha(source_root / 'summary.json'),
           'source_settings_sha256': sha(source_root / 'settings.json'),
           'restarts': args.restarts, 'stage12_nfev': args.stage12_nfev,
           'stage18_nfev': args.stage18_nfev,
           'target_chi2': args.target_chi2,
           'window_factor': args.window_factor,
           'window_floor': args.window_floor,
           'absolute_bound': args.absolute_bound,
           'perturbation': args.perturbation,
           'sensitivity_step': args.sensitivity_step,
           'seed': args.seed, 'd0_mu2_fraction': args.d0_mu2_fraction}
    if args.preflight:
        for model in selected:
            row, old = source_info(model, source_root)
            print(model, row['status'], 'chi2=' + str(row.get('chi2_final')),
                  '18-entry seed=' + ('yes' if old else 'none (D+0 blocked)'), flush=True)
        print('Source:', source_root)
        print('Destination:', dest)
        print('No numerical integrations or file writes performed.')
        return 0
    if (dest / 'settings.json').exists():
        if not args.resume:
            p.error('Output exists; use --resume with identical settings or choose a new --run-id')
        if read(dest / 'settings.json') != sig:
            p.error('Resume signature mismatch; use another run ID')
    else:
        if dest.exists() and any(dest.iterdir()):
            p.error('Destination exists and is nonempty; choose a new --run-id')
        save(dest / 'settings.json', sig)
    rows: list[dict] = []
    for i, model in enumerate(selected, 1):
        result_file = dest / tag(model) / 'result.json'
        error_file = dest / tag(model) / 'error.json'
        if args.resume and result_file.is_file():
            old = read(result_file)
            row = {'model': model, 'status': old['status'],
                   'chi2_previous': old.get('source_chi2'),
                   'chi2_final': old.get('chi2_final'),
                   'use_for_merged_summary': old.get('use_for_merged_summary', False),
                   'scalar_boundary_modified': old['matching_and_rge_provenance'].get('modified_scalar_boundary', False),
                   'result_path': str(result_file.relative_to(ROOT))}
            print(f'[{i}/{len(selected)}] Reusing {model}: {row["status"]}', flush=True)
        else:
            print(f'[{i}/{len(selected)}] Continuing {model} (full-general; 18 Yukawas)', flush=True)
            try:
                row = run_one(model, args, source_root, source_settings, dest / tag(model))
                error_file.unlink(missing_ok=True)
            except Exception as exc:
                error = {'format': FORMAT, 'model': model, 'status': 'Blocked',
                         'error': f'{type(exc).__name__}: {exc}',
                         'traceback': traceback.format_exc(limit=14)}
                save(error_file, error)
                row = {'model': model, 'status': 'Blocked',
                       'chi2_previous': rows0[model].get('chi2_final'),
                       'chi2_final': None, 'use_for_merged_summary': False,
                       'scalar_boundary_modified': model == 'D_+0' and args.d0_mu2_fraction != float(source_settings['mu2_fraction']),
                       'message': error['error']}
                print('  BLOCKED:', error['error'], flush=True)
        rows.append(row)
        write_summary(dest, rows, source_root)
    print('New results:', dest / 'summary.csv')
    print('Thesis-ready merged overview:', dest / 'merged_all16.csv')
    return 2 if any(x['status'] == 'Blocked' for x in rows) else 0


if __name__ == '__main__':
    raise SystemExit(main())
