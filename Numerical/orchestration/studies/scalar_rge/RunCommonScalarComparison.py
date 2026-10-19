"""Cross-model common-lambdaT3 sensitivity scan around independently fitted T3 benchmarks.

Run from repository root:
    python -m Numerical.orchestration.studies.scalar_rge.RunCommonScalarComparison           # inventory
    python -m Numerical.orchestration.studies.scalar_rge.RunCommonScalarComparison --run --limit 2
    python -m Numerical.orchestration.studies.scalar_rge.RunCommonScalarComparison --run

The shared lambda grid is calibrated using a *single reference model*: its
largest Yukawa singular value is set to each requested target at the UV scale.
All models are subsequently evaluated at exactly the same complex lambdaT3
values, but their Yukawa singular values generally differ. For each model,
y1 and y2 are multiplied by sqrt(lambda_old/lambda_new). Consequently the
leading scalar-Yukawa neutrino-mass prefactor is exactly preserved at the UV
input scale (including complex phases). Running and thresholds can still shift
low-energy observables. No benchmarks are refitted.

This requires ordinary split-scalar T3 configurations and their precomputed
RGBeta/matching artifacts. The original results are never modified.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from Numerical.orchestration.support.ResumePipelineNumericalResults import (
    PROJECT_ROOT, OUTPUT_ROOT, REQUIRED, discover_cases,
)
from Numerical.orchestration.studies.fitted_rge.RunFittedSensitivityComparison import (
    SCENARIO, OBS, best_candidate, model_class, read_json,
)
from Numerical.core.RGBetaEvaluator import load_rgbeta_payload
from Numerical.running.FinalC5TrajectoryAdapter import FinalC5TrajectoryEvaluator
from Numerical.fitting.OscillationFit import OscillationFitTarget
from Numerical.fitting.ParameterScan import make_t3_point_evaluator
from Numerical.fitting.ScanCLI import build_uv_state_from_config

TARGET_PATH = PROJECT_ROOT / 'configs' / 'nufit_6_1_ic24_no.json'
DEFAULT_TARGETS = (1e-4, 1e-3, 1e-2, 1e-1, 1.0)


def _get_complex_scalar(spec):
    if isinstance(spec, dict):
        return complex(float(spec.get('real', 0)), float(spec.get('imag', 0)))
    return complex(spec)


def _get_matrix(spec):
    real = np.asarray(spec['real'], dtype=float)
    imag = np.asarray(spec.get('imag', np.zeros_like(real)), dtype=float)
    if real.shape != (3, 3) or imag.shape != (3, 3):
        raise ValueError('Expected a 3x3 Yukawa matrix')
    return real + 1j * imag


def _set_matrix(spec, matrix):
    spec['real'] = np.asarray(matrix.real, dtype=float).tolist()
    spec['imag'] = np.asarray(matrix.imag, dtype=float).tolist()


def characteristic_yukawa(payload):
    ordinary = payload['base_state']['ordinary']
    return max(float(np.linalg.svd(_get_matrix(ordinary[k]), compute_uv=False)[0])
               for k in ('y1', 'y2'))


def deform_common_scalar(payload, new_lambda):
    """Preserve lambdaT3 * (y1*y2 + y2*y1) at the UV input scale."""
    result = deepcopy(payload)
    ordinary = result['base_state']['ordinary']
    old_lambda = _get_complex_scalar(ordinary['lambdaT3'])
    if abs(old_lambda) == 0 or abs(new_lambda) == 0:
        raise ValueError('Cannot compensate a zero lambdaT3')
    factor = np.sqrt(old_lambda / new_lambda)
    ordinary['lambdaT3'] = {'real': float(new_lambda.real), 'imag': float(new_lambda.imag)}
    for key in ('y1', 'y2'):
        _set_matrix(ordinary[key], _get_matrix(ordinary[key]) * factor)
    return result


def _view(payload):
    config = deepcopy(payload)
    config.setdefault('scan', {'bindings': {}})
    return SimpleNamespace(raw=config)


def evaluate(payload, paths, target):
    if bool(payload.get('representation', {}).get('shared_scalar', False)):
        raise ValueError('Shared-scalar field-content branch is not supported by this evaluator')
    uv = load_rgbeta_payload(paths[0])
    eft = load_rgbeta_payload(paths[1])
    c5 = FinalC5TrajectoryEvaluator(paths[2])
    scales = payload['scales']
    numerical = payload.get('numerical', {})
    point = make_t3_point_evaluator(
        uv_state_builder=lambda _parameters: build_uv_state_from_config(_view(payload), {}),
        c5_builder=c5, uv_rgbeta_payload=uv, eft1_rgbeta_payload=eft,
        target=target,
        mu_fermion_threshold_gev=float(scales['mu_fermion_threshold_gev']),
        mu_scalar_threshold_gev=float(scales['mu_scalar_threshold_gev']),
        mu_low_gev=float(scales['mu_low_gev']),
        vev_gev=float(numerical.get('vev_gev', 246.22)),
        rtol=float(numerical.get('rtol', 1e-8)),
        atol=float(numerical.get('atol', 1e-11)),
    )({})
    return {'chi2': float(point.chi2),
            'prediction': {k: float(v) for k, v in zip(point.observable_names, point.prediction)}}


def make_plots(rows, output, target):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print('matplotlib unavailable; data files are still complete')
        return
    from itertools import combinations
    labels = {
        'delta_m21_sq_ev2': r'$\Delta m_{21}^{2}\ (\mathrm{eV}^{2})$',
        'delta_m3l_sq_ev2': r'$\Delta m_{3\ell}^{2}\ (\mathrm{eV}^{2})$',
        'sin2_theta12': r'$\sin^{2}\theta_{12}$',
        'sin2_theta13': r'$\sin^{2}\theta_{13}$',
        'sin2_theta23': r'$\sin^{2}\theta_{23}$',
    }
    plotdir = output / 'plots'
    plotdir.mkdir(parents=True, exist_ok=True)
    valid = [r for r in rows if r['status'] in ('Success', 'PartialSuccess')]
    center = dict(zip(target.observable_names, target.central_values))
    palette = plt.get_cmap('tab20')
    for x, y in combinations(OBS, 2):
        fig, ax = plt.subplots(figsize=(8, 6))
        for i, row in enumerate(valid):
            result = read_json(Path(row['scan_json']))
            points = [p for p in result['points'] if p['status'] == 'Success']
            if not points:
                continue
            ax.plot([p['prediction'][x] for p in points],
                    [p['prediction'][y] for p in points],
                    '-o', lw=1.3, ms=3.5, color=palette(i % 20), label=row['model'])
        if x in center and y in center:
            ax.scatter(center[x], center[y], marker='x', s=70,
                       color='black', zorder=10, label='NuFIT centre')
        ax.set(xlabel=labels[x], ylabel=labels[y],
               title='Common scalar coupling, compensated Yukawas')
        ax.grid(alpha=0.2)
        ax.legend(fontsize=6, ncol=2)
        fig.tight_layout()
        fig.savefig(plotdir / f'{x}_vs_{y}.png', dpi=180)
        plt.close(fig)
    print('Wrote observable-pair figures to:', plotdir)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scenario', default=SCENARIO)
    parser.add_argument('--run', action='store_true')
    parser.add_argument('--limit', type=int, default=0)
    parser.add_argument('--reference-model', help='Model key for calibration; default is first eligible sorted model')
    parser.add_argument('--target-yukawas', nargs='+', type=float, default=list(DEFAULT_TARGETS))
    parser.add_argument('--max-yukawa', type=float, default=1.0,
                        help='Maximum allowed Yukawa singular value for every model and point')
    parser.add_argument('--force', action='store_true')
    args = parser.parse_args(argv)
    targets = [float(v) for v in args.target_yukawas]
    if not targets or any(not math.isfinite(v) or v <= 0 for v in targets):
        parser.error('target-yukawas must be finite and positive')
    if targets != sorted(set(targets)):
        parser.error('target-yukawas must be strictly increasing')
    if not math.isfinite(args.max_yukawa) or args.max_yukawa <= 0:
        parser.error('--max-yukawa must be positive')
    if args.limit < 0:
        parser.error('--limit must be >=0')

    cases = {}
    for scenario, key, case in discover_cases(OUTPUT_ROOT):
        if scenario == args.scenario:
            if key in cases:
                parser.error(f'Duplicate case: {key}')
            cases[key] = case
    selected = {}
    for key in sorted(cases):
        fit = best_candidate(args.scenario, key)
        if fit is None or any(not (cases[key] / 'data' / name).is_file() for name in REQUIRED):
            continue
        payload = read_json(fit[2])
        if payload.get('representation', {}).get('shared_scalar', False):
            continue
        selected[key] = (fit, payload, [cases[key] / 'data' / name for name in REQUIRED])
    if not selected:
        print('No eligible fitted ordinary-T3 cases with required numerical artifacts.')
        return 1
    reference = args.reference_model or next(iter(selected))
    if reference not in selected:
        parser.error('Reference model is not eligible; choose one of: ' + ', '.join(selected))
    reference_payload = selected[reference][1]
    lambda_ref = _get_complex_scalar(reference_payload['base_state']['ordinary']['lambdaT3'])
    y_ref = characteristic_yukawa(reference_payload)
    if abs(lambda_ref) == 0 or y_ref == 0:
        parser.error('Reference model must have nonzero lambdaT3 and Yukawas')
    # Common complex lambda values, calibrated to the reference model Yukawa norm.
    shared_lambdas = [lambda_ref * (y_ref / target_y)**2 for target_y in targets]
    print('Reference:', reference, 'lambdaT3=', lambda_ref, 'y_ref=', y_ref)
    for target_y, lam in zip(targets, shared_lambdas):
        print(f'  reference Yukawa target={target_y:g} => common lambdaT3={lam!s}')
    keys = list(selected)[:args.limit or None]
    print('Eligible cases:', len(selected), '; selected:', len(keys))
    if not args.run:
        print('Dry run; supply --run to calculate trajectories and plots.')
        return 0

    target = OscillationFitTarget.from_json(TARGET_PATH)
    output = PROJECT_ROOT / 'Reports' / 'output' / 'full' / 'common_scalar_comparison' / args.scenario
    output.mkdir(parents=True, exist_ok=True)
    rows = []
    for key in keys:
        (fit, base, paths) = selected[key]
        points = []
        original_lambda = _get_complex_scalar(base['base_state']['ordinary']['lambdaT3'])
        scan_path = output / key / 'common_scalar_scan.json'
        if scan_path.is_file() and not args.force:
            existing = read_json(scan_path)
            if existing.get('reference_model') != reference or existing.get('reference_yukawa_targets') != targets or existing.get('common_lambdaT3') != [[v.real, v.imag] for v in shared_lambdas] or existing.get('max_yukawa') != args.max_yukawa:
                raise RuntimeError(f'Existing scan settings differ: {scan_path}; use --force')
            points = existing['points']
        else:
            for target_y, lam in zip(targets, shared_lambdas):
                point = {'reference_yukawa_target': target_y,
                         'lambdaT3': {'real': float(lam.real), 'imag': float(lam.imag)},
                         'status': 'Pending', 'prediction': None, 'chi2': None}
                try:
                    deformed = deform_common_scalar(base, lam)
                    actual = characteristic_yukawa(deformed)
                    point['actual_yukawa_max_singular_value'] = actual
                    if actual > args.max_yukawa * (1 + 1e-12):
                        point.update(status='SkippedNonPerturbative', error=f'Max Yukawa singular value {actual:g} exceeds {args.max_yukawa:g}')
                    else:
                        point.update(evaluate(deformed, paths, target))
                        point['status'] = 'Success'
                except Exception as exc:
                    point.update(status='Failed', error=f'{type(exc).__name__}: {exc}')
                points.append(point)
                print(key, 'target=', target_y, point['status'], flush=True)
            scan_path.parent.mkdir(parents=True, exist_ok=True)
            scan_path.write_text(json.dumps({
                'model': key, 'reference_model': reference,
                'reference_yukawa_targets': targets,
                'reference_baseline_yukawa': y_ref,
                'reference_baseline_lambdaT3': [lambda_ref.real, lambda_ref.imag],
                'original_lambdaT3': [original_lambda.real, original_lambda.imag],
                'common_lambdaT3': [[v.real, v.imag] for v in shared_lambdas],
                'max_yukawa': args.max_yukawa,
                'definition': 'lambdaT3 set to common complex value; y1,y2 *= sqrt(old_lambda/new_lambda)',
                'points': points,
            }, indent=2) + '\n', encoding='utf-8')
        successful = sum(p['status'] == 'Success' for p in points)
        status = 'Success' if successful == len(points) else ('PartialSuccess' if successful else 'Failed')
        rows.append({'model': key, 'representation': model_class(key),
                     'fit_stage': fit[1], 'fit_chi2': fit[0],
                     'successful_points': successful, 'total_points': len(points),
                     'status': status, 'scan_json': str(scan_path)})
    with (output / 'comparison_summary.csv').open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (output / 'comparison_summary.json').write_text(json.dumps(rows, indent=2) + '\n', encoding='utf-8')
    make_plots(rows, output, target)
    print('Summary:', output / 'comparison_summary.csv')
    return 1 if any(r['status'] == 'Failed' for r in rows) else 0


if __name__ == '__main__':
    raise SystemExit(main())
