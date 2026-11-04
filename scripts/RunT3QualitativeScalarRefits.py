"""Three-point qualitative general-T3 scalar/Yukawa refit study.

Run as a module from the repository root. All physics is delegated to the
current *general* one-scale RGE/C5 evaluators; no restricted matching fallback.
The UV -> fermions -> scalars -> SMEFT sequential matching is NOT certified.

Two deliberately distinct scalar boundary policies:
 lambda-only: only lambdaT3 changes; additional physical UV quartics stay fixed.
 joint-quartics: extra general quartics change in proportion to |lambdaT3|;
                 ordinary scalar couplings and scalar bilinears stay fixed.

The B/C backend only permits 0 < quartic_fraction <= 1 and defines its seven
extra quartics as quartic_fraction * |lambdaT3|. Impossible lambda-only
compensations are explicitly UNSUPPORTED, never clamped.

Reoptimises all 18 *real* y1/y2 components at each allowed lambda value.
Conjugate Yukawas, where permitted, are retained at the saved UV boundary;
they are not independent optimiser parameters. No source files are changed.
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
import traceback
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'output/studies/thesis_qualitative/scalar_refits'
BENCHMARKS = (0.0005, 0.005, 0.5)
MODELS = tuple(f'{family}_{alpha:+d}' for family, alphas in {
    'A': (-4, -2, 0), 'B': (-3, -1, 1), 'C': (-3, -1, 1),
    'D': (-2, 0, 2), 'E': (-4, -2, 0, 2),
}.items() for alpha in alphas)
MIXED = frozenset(('B_-1', 'C_-1'))
NAMES = tuple(f'{name}_{i}{j}' for name in ('y1', 'y2') for i in range(3) for j in range(3))
OBS = ('delta_m21_sq_ev2', 'delta_m3l_sq_ev2',
       'sin2_theta12', 'sin2_theta13', 'sin2_theta23')
SCHEMA = 't3-thesis-qualitative-general-scalar-refit-v1'


class UnsupportedScalarBoundary(ValueError):
    """Backend cannot represent the requested independent scalar boundary."""


def tag(model: str) -> str:
    return model.replace('-', 'm').replace('+', 'p')


def lambda_tag(value: float) -> str:
    return 'lambda_' + format(value, '.12g').replace('.', 'p')


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def json_read(path: Path) -> dict:
    value = json.loads(path.read_text(encoding='utf-8-sig'))
    if not isinstance(value, dict):
        raise ValueError(f'Expected object: {path}')
    return value


def json_save(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.tmp')
    tmp.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    tmp.replace(path)


def real_lambda(config: dict) -> float:
    raw = config['base_state']['ordinary']['lambdaT3']
    z = complex(float(raw['real']), float(raw.get('imag', 0)))
    if not np.isfinite([z.real, z.imag]).all() or z.imag != 0 or z.real <= 0:
        raise ValueError('Current study requires positive real fitted lambdaT3')
    return float(z.real)


def fitted_vector(config: dict) -> np.ndarray:
    values = []
    for key in ('y1', 'y2'):
        spec = config['base_state']['ordinary'][key]
        real = np.asarray(spec['real'], dtype=float)
        imag = np.asarray(spec.get('imag', np.zeros((3, 3))), dtype=float)
        if real.shape != (3, 3) or imag.shape != (3, 3) or not np.isfinite(real).all() or not np.isfinite(imag).all():
            raise ValueError(f'Non-finite or incorrectly sized {key}')
        if np.max(abs(imag)) > 1e-12:
            raise ValueError('Current evaluator fits 18 real y1/y2 only')
        values.append(real.ravel())
    return np.concatenate(values)


def parameter_dict(x: np.ndarray) -> dict[str, float]:
    x = np.asarray(x, dtype=float)
    if x.shape != (18,) or not np.isfinite(x).all():
        raise ValueError('Expected 18 finite Yukawa parameters')
    return dict(zip(NAMES, map(float, x)))


def norms(x: np.ndarray) -> dict[str, float]:
    z = np.asarray(x, dtype=float).reshape(2, 3, 3)
    return {'y1_frobenius': float(np.linalg.norm(z[0])),
            'y2_frobenius': float(np.linalg.norm(z[1])),
            'combined_frobenius': float(np.linalg.norm(z)),
            'y1_max_abs': float(np.max(abs(z[0]))),
            'y2_max_abs': float(np.max(abs(z[1])))}


def yukawa_matrices(x: np.ndarray) -> dict:
    values = np.asarray(x, dtype=float).reshape(2, 3, 3)
    return {k: {'real': values[i].tolist(), 'imag': np.zeros((3, 3)).tolist()}
            for i, k in enumerate(('y1', 'y2'))}


def quartic_fraction(q0: float, old: float, new: float, mode: str) -> float:
    if not (np.isfinite([q0, old, new]).all() and 0 < q0 <= 1 and old > 0 and new > 0):
        raise ValueError('Invalid B/C quartic input')
    if mode == 'joint-quartics':
        return float(q0)
    if mode != 'lambda-only':
        raise ValueError('Unknown scalar policy')
    required = q0 * old / new
    if required > 1 or required <= 0:
        raise UnsupportedScalarBoundary(
            f'lambda-only cannot preserve B/C extra quartics: saved fraction={q0:g}, '
            f'lambda_old={old:g}, lambda_new={new:g}, required fraction={required:.9g} '
            '(backend requires 0 < quartic_fraction <= 1)')
    return float(required)


def scale_extra_quartics(frozen: dict, config: dict, factor: float) -> tuple[dict, list[str]]:
    """Only the *additional* UV quartics, not ordinary or SM couplings."""
    result = deepcopy(frozen)
    occupied = set(config['base_state']['ordinary']) | set(config['base_state']['sm'])
    altered = []
    for key, raw in frozen.items():
        if key in occupied or not key.startswith('lambda'):
            continue
        if isinstance(raw, (int, float)):
            result[key] = float(raw) * factor
        elif isinstance(raw, dict) and {'real', 'imag'} <= set(raw) and not isinstance(raw['real'], (dict, list)):
            result[key] = {'real': float(raw['real']) * factor,
                           'imag': float(raw['imag']) * factor}
        else:
            raise ValueError(f'Unrecognised extra quartic representation: {key}')
        altered.append(key)
    return result, sorted(altered)


def scalar_setup(model: str, config: dict, fit: dict, paths: dict,
                 value: float, mode: str, target: Any, canonical_run: str):
    """Construct a persistent evaluator for one lambda; original modules own physics."""
    from scripts import RunT3GeneralZeeSensitivity as zee
    from scripts import RunT3GeneralAll16FromReference as general
    old = real_lambda(config)
    changed = deepcopy(config)
    changed['base_state']['ordinary']['lambdaT3'] = {'real': value, 'imag': 0.0}
    if model not in MIXED:
        uv, eft, provenance, kernel, matched, scheme = general.nonmixed_inputs(
            model, canonical_run, generate_general_c5=False)
        if matched.resolve() != paths['matching'].resolve():
            raise ValueError('General Matchete matching source changed')
        frozen = json_read(paths['overrides'])
        if mode == 'lambda-only':
            fixed, scaled = deepcopy(frozen), []
        else:
            fixed, scaled = scale_extra_quartics(frozen, config, value / old)
        fixed['lambdaT3'] = {'real': value, 'imag': 0.0}
        evaluator = general.NonMixedEvaluator(model, uv, eft, kernel, fixed, target)
        return evaluator, {
            'matching': scheme, 'quartic_boundary': mode,
            'additional_quartics_scaled': scaled,
            'additional_quartics_held_fixed': mode == 'lambda-only',
            'conjugate_yukawas_held_fixed': True,
            'scalar_bilinears_held_fixed': True,
            'ordinary_scalar_inputs_held_fixed_except_lambdaT3': True}

    from scripts import RunT3BCGeneralFromReference as br
    from scripts import RunT3GeneralScalarBCExploratoryFit as bc
    q0, mu0 = zee._bc_input_fractions(fit)
    q = quartic_fraction(q0, old, value, mode)
    # Tracked reference supplies ONLY the registry/inventory metadata. The
    # physical state below is replaced by the merged general fitted_config.
    source = br.build_reference_source(model)
    source['cfg'] = changed
    source['lambdaT3'] = complex(value)
    bundle = bc.beta_bundle(model, source, source['registry'])
    coefficient = bc.compile_source(paths['matching'])
    evaluator = bc.CandidateEvaluator('reference_scalar_general', model, source,
        source['registry'], bundle, coefficient, target,
        quartic_fraction=q, mu2_fraction=mu0)
    return evaluator, {
        'matching': 'conditional-mixed-scalar-one-scale-signed-quartic-map',
        'quartic_boundary': mode,
        'saved_quartic_fraction': q0,
        'effective_quartic_fraction': q,
        'saved_extra_quartic_magnitude_at_uv': q0 * old,
        'new_extra_quartic_magnitude_at_uv': q * value,
        'extra_mixed_scalar_bilinear_fraction_fixed': mu0,
        'additional_quartics_held_fixed': mode == 'lambda-only',
        'additional_quartics_scaled': list(source['new'] + source['reused']) if mode == 'joint-quartics' else [],
        'conjugate_yukawas_held_fixed': True,
        'scalar_bilinears_held_fixed': True,
        'ordinary_scalar_inputs_held_fixed_except_lambdaT3': True,
        'signed_direct_quartic_map_physically_certified': False}


def fit_prediction(result: Any) -> dict[str, float]:
    out = dict(zip(result.observable_names, map(float, result.prediction)))
    if set(out) != set(OBS) or not np.isfinite(list(out.values())).all():
        raise ValueError('Bad five-observable prediction')
    return {k: out[k] for k in OBS}


def low_masses(details: dict, config: dict, ordering: str) -> list[float]:
    from physics.NeutrinoTrajectory import neutrino_mass_from_c5
    from physics.NeutrinoObservables import calculate_neutrino_observables
    from scripts import RunT3GeneralZeeSensitivity as zee
    raw = details.get('c5_at_low_scale', details.get('C5_low_gev_inv'))
    if raw is None:
        raise ValueError('Evaluator did not return low-scale C5')
    c5 = zee.matrix(raw)
    vev = float(config.get('numerical', {}).get('vev_gev', 246.22))
    matrix = neutrino_mass_from_c5(c5, vev_gev=vev)
    obs = calculate_neutrino_observables(matrix, ordering=ordering)
    return list(map(float, obs.masses_ev))


def baseline_check(evaluator, original: np.ndarray, fit: dict, tolerance: float = 1e-5) -> dict:
    baseline = evaluator(parameter_dict(original))
    prediction = fit_prediction(baseline)
    for key in OBS:
        expected = float(fit['prediction'][key])
        if abs(prediction[key] - expected) > max(1e-10, tolerance * abs(expected)):
            raise RuntimeError(f'Merged baseline mismatch for {key}: current={prediction[key]}, saved={expected}')
    return {'chi2': float(baseline.chi2), 'prediction': prediction}


def optimize(evaluator, original: np.ndarray, old: float, value: float, *,
             bound: float, max_nfev: int, regularization: float) -> tuple[Any, np.ndarray, dict]:
    """All 18 components are active; small UV-seed prior selects a local branch.

    Because five oscillation observables cannot identify 18 real Yukawas, the
    regularizer makes the seed-local comparison well-defined, not unique physics.
    """
    from scipy.optimize import least_squares
    scale_factor = math.sqrt(old / value)
    seed = original * scale_factor  # heuristic only, NOT a fitted solution
    if max(abs(seed)) >= bound:
        raise UnsupportedScalarBoundary('Rescaled Yukawa seed exceeds --bound; raise it explicitly')
    sizes = np.maximum(abs(seed), max(float(np.max(abs(seed))) * 0.1, 1e-7))
    best = {'score': math.inf, 'x': None, 'fit': None}
    failures: list[str] = []
    evaluations = 0
    failed_evaluations = 0

    def residual(x):
        nonlocal evaluations, failed_evaluations
        evaluations += 1
        reg = regularization * (x - seed) / sizes
        try:
            fit = evaluator(parameter_dict(x))
            physics = np.asarray(fit.whitened_residual, dtype=float)
            if physics.shape != (5,) or not np.isfinite(physics).all():
                raise ValueError('Invalid whitened residual vector')
            score = float(fit.chi2)
            if np.isfinite(score) and score < best['score']:
                best.update(score=score, x=np.array(x, copy=True), fit=fit)
        except Exception as exc:
            failed_evaluations += 1
            if len(failures) < 12:
                failures.append(f'{type(exc).__name__}: {exc}'[:400])
            physics = np.full(5, 1e4)
        return np.concatenate((physics, reg))

    # Evaluate the original best-fit-derived seed, then perform a NEW fit at
    # every lambda point; no result is generated by rescaling alone.
    residual(seed)
    outcome = least_squares(residual, seed, method='trf', jac='2-point',
                            bounds=(-bound, bound), x_scale=sizes,
                            max_nfev=max_nfev, ftol=1e-6, xtol=1e-6, gtol=1e-6)
    residual(outcome.x)
    if best['x'] is None:
        raise RuntimeError('No valid numerical evaluation during fit: ' + repr(failures[:3]))
    return best['fit'], best['x'], {
        'success': bool(outcome.success), 'message': str(outcome.message),
        'nfev_reported_by_scipy': int(outcome.nfev),
        'n_evaluator_attempts': int(evaluations),
        'n_failed_evaluations': int(failed_evaluations),
        'max_nfev_requested': int(max_nfev),
        'first_failures': failures,
        'all_18_real_yukawas_active': True,
        'seed_rescale_factor_heuristic': scale_factor,
        'branch_regularization_strength': regularization,
        'jacobian': '2-point', 'scipy_method': 'trf',
        'best_chi2': best['score']}


def point_run(model: str, config: dict, fit: dict, paths: dict, lam: float,
              args, target: Any, original: np.ndarray, baseline: dict) -> dict:
    old = real_lambda(config)
    result = {'model': model, 'lambdaT3': lam, 'mode': args.mode,
              'status': 'Failed', 'source_chi2': float(fit['chi2_final']),
              'baseline_recomputed_chi2': float(baseline['chi2'])}
    try:
        evaluator, boundary = scalar_setup(model, config, fit, paths, lam,
                                         args.mode, target, args.canonical_run_id)
        result['scalar_boundary'] = boundary
        _, optimum, history = optimize(evaluator, original, old, lam,
            bound=args.bound, max_nfev=args.max_nfev,
            regularization=args.regularization)
        final, details = evaluator.evaluate(parameter_dict(optimum), detail=True)
        pred = fit_prediction(final)
        mass = low_masses(details, config, str(final.ordering))
        result.update(status='AcceptedApproxFit' if float(final.chi2) <= args.chi2_cutoff
                      else 'PoorFit', chi2=float(final.chi2),
                      prediction=pred, neutrino_masses_ev=mass,
                      ordering=str(final.ordering),
                      whitened_residual=np.asarray(final.whitened_residual, dtype=float).tolist(),
                      fitted_yukawas=yukawa_matrices(optimum),
                      yukawa_norms=norms(optimum),
                      norm_ratio_to_reference=norms(optimum)['combined_frobenius'] / norms(original)['combined_frobenius'],
                      optimizer=history, numerical_details=details)
    except UnsupportedScalarBoundary as exc:
        result.update(status='UnsupportedScalarBoundary', reason=str(exc))
    except Exception as exc:
        result.update(status='Failed', reason=f'{type(exc).__name__}: {exc}'[:1500],
                      traceback=traceback.format_exc(limit=5)[-2500:])
    return result


def save_fitted_inputs(folder: Path, model: str, config: dict, paths: dict,
                       lam: float, mode: str, row: dict) -> None:
    """Save full local replay inputs; original accepted sources are never edited."""
    if row.get('status') not in ('AcceptedApproxFit', 'PoorFit'):
        return
    cfg = deepcopy(config)
    cfg['base_state']['ordinary']['lambdaT3'] = {'real':lam, 'imag':0.0}
    cfg['base_state']['ordinary'].update(deepcopy(row['fitted_yukawas']))
    cfg['thesis_qualitative_scalar_refit'] = {
        'schema': SCHEMA, 'mode':mode,
        'scalar_boundary': row['scalar_boundary'],
        'conditional_one_scale_only':True,
        'physics_modules_not_modified':True}
    json_save(folder / 'fitted_config.json', cfg)
    if model not in MIXED:
        from scripts import RunT3GeneralAll16FromReference as general
        frozen = json_read(paths['overrides'])
        if mode == 'joint-quartics':
            frozen, _ = scale_extra_quartics(frozen, config, lam/real_lambda(config))
        frozen['lambdaT3'] = {'real':lam, 'imag':0.0}
        vector = fitted_vector(cfg)
        overrides = general.override_yukawas(frozen, parameter_dict(vector))
        json_save(folder / 'general_uv_overrides.json', overrides)


def summarize(rows: list[dict]) -> dict:
    accepted = sorted((r for r in rows if r['status'] == 'AcceptedApproxFit'),
                      key=lambda r: r['lambdaT3'])
    slope = None
    if len(accepted) >= 2:
        a, b = accepted[0], accepted[-1]
        n0 = a['yukawa_norms']['combined_frobenius']
        n1 = b['yukawa_norms']['combined_frobenius']
        if n0 > 0 and n1 > 0:
            slope = math.log(n1 / n0) / math.log(b['lambdaT3'] / a['lambdaT3'])
    statuses = tuple(r['status'] for r in rows)
    complete = all(status in ('AcceptedApproxFit', 'UnsupportedScalarBoundary')
                   for status in statuses)
    # A model-level status is required by main() when tallying blocked models.
    # Explicitly unsupported B/C points are not failed numerical evaluations.
    if complete:
        model_status = ('CompleteWithUnsupported' if
                        'UnsupportedScalarBoundary' in statuses else 'Complete')
    else:
        model_status = 'Incomplete'
    return {'model': rows[0]['model'],
            'status': model_status,
            'benchmark_status': {format(r['lambdaT3'], '.12g'): r['status'] for r in rows},
            'accepted_count': len(accepted),
            'log_log_yukawa_norm_vs_lambda_slope_from_accepted_endpoints': slope,
            'naive_lambda_inverse_sqrt_slope': -0.5,
            'slope_is_diagnostic_not_an_exact_model_prediction': True}


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--models', default='pilot',
                   help='pilot = A_-4,B_-1; all = 16 models; or comma-separated')
    p.add_argument('--mode', choices=('lambda-only', 'joint-quartics'), default='lambda-only')
    p.add_argument('--canonical-run-id', default='full16_general_v1')
    p.add_argument('--run-id', default='three_benchmark_v1')
    p.add_argument('--bound', type=float, default=0.1)
    p.add_argument('--max-nfev', type=int, default=30)
    p.add_argument('--regularization', type=float, default=1e-3)
    p.add_argument('--chi2-cutoff', type=float, default=5.86)
    p.add_argument('--preflight', action='store_true', help='Read-only inputs and boundary checks')
    args = p.parse_args(argv)
    wanted = (('A_-4','B_-1') if args.models == 'pilot' else
              MODELS if args.models == 'all' else
              tuple(x.strip() for x in args.models.split(',') if x.strip()))
    if not wanted or len(wanted) != len(set(wanted)) or set(wanted) - set(MODELS):
        p.error('Invalid model list')
    if (not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,63}', args.run_id) or
        not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,63}', args.canonical_run_id)):
        p.error('Invalid run identifier')
    if not (args.bound > 0 and np.isfinite(args.bound) and args.max_nfev >= 1 and
            0 <= args.regularization <= 0.1 and args.chi2_cutoff > 0):
        p.error('Invalid optimiser options')
    # Import current source only after purely local CLI validation.
    from scripts import RunT3GeneralZeeSensitivity as zee
    selection = zee.fit_selections(list(wanted), 'merged')
    dest = OUT / args.run_id
    if not args.preflight and dest.exists():
        p.error(f'Refusing to overwrite existing outputs: {dest}')
    target = None if args.preflight else zee.target_object()
    all_rows = []
    model_summaries = []
    for model in wanted:
        print(f'{model}: checking accepted merged source', flush=True)
        try:
            source_args = argparse.Namespace(fit_run_id='merged',
                            canonical_run_id=args.canonical_run_id)
            config, fit, paths = zee.sources(model, source_args, selection[model])
            old = real_lambda(config)
            original = fitted_vector(config)
            allowed = {}
            for lam in BENCHMARKS:
                if model in MIXED:
                    q0, _ = zee._bc_input_fractions(fit)
                    try:
                        q = quartic_fraction(q0, old, lam, args.mode)
                        allowed[str(lam)] = {'status':'Supported', 'quartic_fraction':q}
                    except UnsupportedScalarBoundary as exc:
                        allowed[str(lam)] = {'status':'UnsupportedScalarBoundary', 'reason':str(exc)}
                else:
                    allowed[str(lam)] = {'status':'Supported'}
            if args.preflight:
                print(f'  source={selection[model]["source"]}, base_lambda={old:g}, '
                      f'points={allowed}', flush=True)
                continue
            # Require reproduction of accepted-fit source, never fit against a stale baseline.
            baseline_evaluator, baseline_boundary = scalar_setup(model,config,fit,paths,old,
                                           'lambda-only',target,args.canonical_run_id)
            baseline = baseline_check(baseline_evaluator, original, fit)
            model_rows = []
            for lam in BENCHMARKS:
                print(f'  lambdaT3={lam:g}, {args.mode}', flush=True)
                if allowed[str(lam)]['status'] == 'UnsupportedScalarBoundary':
                    row = {'model': model, 'lambdaT3': lam, 'mode': args.mode,
                           'status': 'UnsupportedScalarBoundary',
                           'reason': allowed[str(lam)]['reason']}
                else:
                    row = point_run(model, config, fit, paths, lam, args, target,
                                    original, baseline)
                row['source'] = selection[model]
                row['input_sha256'] = {k:sha(v) for k,v in paths.items()}
                row['input_paths'] = {k:v.relative_to(ROOT).as_posix() for k,v in paths.items()}
                row['baseline_lambdaT3'] = old
                row['original_yukawa_norms'] = norms(original)
                model_rows.append(row)
                point_dir = dest / tag(model) / lambda_tag(lam)
                json_save(point_dir / 'result.json', row)
                save_fitted_inputs(point_dir, model, config, paths, lam, args.mode, row)
                all_rows.append(row)
                print('    ' + row['status'] + (f', chi2={row["chi2"]:.6g}' if 'chi2' in row else
                                                ', '+row.get('reason','')[:180]), flush=True)
            model_summaries.append(summarize(model_rows))
        except Exception as exc:
            error = {'model':model,'status':'Blocked',
                     'reason':f'{type(exc).__name__}: {exc}'[:1500],
                     'source':selection[model]}
            model_summaries.append(error)
            print('  BLOCKED: ' + error['reason'],flush=True)
        if not args.preflight:
            json_save(dest/'summary.json', {'schema':SCHEMA,
                'mode':args.mode,'benchmarks':list(BENCHMARKS),
                'models':model_summaries,
                'conditional_one_scale_not_sequential_fermion_scalar_matching':True,
                'conjugate_yukawas_not_independently_optimised':True,
                'three_benchmarks_only':True,
                'merged_manifest_sha256':sha(zee.MERGED_FILE)})
    if args.preflight:
        print('Read-only preflight finished; no physical evaluations or outputs created.')
        return 0 if not model_summaries else 2
    fields = ('model','lambdaT3','mode','status','chi2',*OBS,'m1_ev','m2_ev','m3_ev',
              'y1_frobenius','y2_frobenius','combined_frobenius','norm_ratio_to_reference','reason')
    with (dest / 'comparison.csv').open('w',encoding='utf-8',newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in all_rows:
            flat = {k:row.get(k,'') for k in ('model','lambdaT3','mode','status','chi2','reason','norm_ratio_to_reference')}
            flat.update(row.get('prediction',{}))
            flat.update(row.get('yukawa_norms',{}))
            for i,m in enumerate(row.get('neutrino_masses_ev', []),1):
                flat[f'm{i}_ev'] = m
            writer.writerow({k: flat.get(k,'') for k in fields})
    accepted = sum(r['status']=='AcceptedApproxFit' for r in all_rows)
    unsupported = sum(r['status']=='UnsupportedScalarBoundary' for r in all_rows)
    blocked = sum(r.get('status') == 'Blocked' for r in model_summaries)
    print(f'Accepted points {accepted}/{len(all_rows)}; explicitly unsupported {unsupported}; '
          f'blocked models {blocked}; output {dest}')
    valid = (not blocked and len(all_rows)==3*len(wanted) and
             all(r['status'] in ('AcceptedApproxFit','UnsupportedScalarBoundary') for r in all_rows))
    return 0 if valid else 2


if __name__ == '__main__':
    raise SystemExit(main())
