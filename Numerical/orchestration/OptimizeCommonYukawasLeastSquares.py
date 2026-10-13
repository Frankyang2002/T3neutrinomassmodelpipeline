"""Bounded local least-squares optimisation of common T3 Yukawas.

Read-only inputs: frozen 16 reference configs, existing trajectory artifacts, and
500-iteration checkpoint. Outputs only to a fresh common_yukawa_least_squares folder.

Usage:
  python -m Numerical.orchestration.OptimizeCommonYukawasLeastSquares
  python -m Numerical.orchestration.OptimizeCommonYukawasLeastSquares --run --max-nfev 40
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import time
from copy import deepcopy
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from Numerical.orchestration.RunCommonParameterComparisonV2 import (
    PROJECT_ROOT, SCENARIO, OBS, candidate_for_target, models_and_cases, load,
)
from Numerical.orchestration.OptimizeCommonYukawas import make_evaluator, BOUND
from Numerical.fitting.OscillationFit import OscillationFitTarget

SOURCE_DIR = PROJECT_ROOT / 'Reports' / 'output' / 'full' / 'common_yukawa_optimization' / SCENARIO
OUTPUT_DIR = PROJECT_ROOT / 'Reports' / 'output' / 'full' / 'common_yukawa_least_squares' / SCENARIO


def save_json_new(path, obj):
    with path.open('x', encoding='utf-8') as f:
        json.dump(obj, f, indent=2, allow_nan=False)
        f.write('\n')


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--run', action='store_true', help='Perform numerical optimisation')
    ap.add_argument('--max-nfev', type=int, default=40, help='Maximum solver residual-vector calls (finite-difference calls count)')
    ap.add_argument('--diff-step', type=float, default=1e-3, help='Relative finite-difference step for normalised Yukawas')
    args = ap.parse_args(argv)
    if args.max_nfev < 1 or args.diff_step <= 0:
        ap.error('Require --max-nfev >= 1 and --diff-step > 0')
    configs, cases = models_and_cases()
    names = sorted(configs)
    cp = load(SOURCE_DIR / 'checkpoint.json')
    common = load(SOURCE_DIR / 'optimized_common_input.json')
    source_name = common['source_model']
    if cp['source_candidate'] != source_name:
        raise ValueError('Checkpoint and common input source models disagree')
    x0 = np.asarray(cp['best_x'], dtype=float)
    if x0.shape != (18,) or not np.all(np.isfinite(x0)) or np.max(np.abs(x0)) > 1:
        raise ValueError('Invalid normalised 18-component checkpoint Yukawas')
    source = load(configs[source_name])
    target = OscillationFitTarget.from_json(PROJECT_ROOT / 'configs' / 'nufit_6_1_ic24_no.json')
    if tuple(target.observable_names) != OBS:
        raise ValueError(f'Unexpected observable order: {target.observable_names}')
    central = np.asarray(target.central_values, dtype=float)
    cov = np.asarray(target.covariance, dtype=float)
    whitening = np.linalg.cholesky(cov)
    payloads = [candidate_for_target(source, load(configs[n])) for n in names]
    evaluators = [make_evaluator(p, cases[n], target) for p, n in zip(payloads, names)]
    print('Models:', len(names), 'parameters:', len(x0), 'observables:', len(names)*len(OBS))
    print('Source checkpoint:', SOURCE_DIR / 'checkpoint.json')
    print('New output:', OUTPUT_DIR)
    print('Max residual evaluations:', args.max_nfev, '(a numerical Jacobian requires ~18 additional evaluations)')
    if not args.run:
        print('Dry run only. Add --run to begin; output directory is not created.')
        return 0
    if OUTPUT_DIR.exists():
        raise FileExistsError('Refusing to overwrite existing least-squares output: ' + str(OUTPUT_DIR))

    calls = 0
    start = time.monotonic()
    best = {'loss': math.inf, 'x': x0.copy(), 'rows': None}
    baseline_chi2 = float(cp['best_average_chi2'])
    # Each residual vector is whitened by the full NuFIT covariance matrix.
    # Division by sqrt(16) makes residual^2 equal mean chi2.
    def evaluate_vector(x):
        nonlocal calls
        calls += 1
        y1 = (np.asarray(x[:9])*BOUND).reshape(3, 3).tolist()
        y2 = (np.asarray(x[9:])*BOUND).reshape(3, 3).tolist()
        for p in payloads:
            p['base_state']['ordinary']['y1']['real'] = y1
            p['base_state']['ordinary']['y2']['real'] = y2
        residuals, rows = [], []
        for name, ev in zip(names, evaluators):
            fit = ev({})
            prediction = np.asarray(fit.prediction, dtype=float)
            if tuple(fit.observable_names) != OBS or prediction.shape != (5,) or not np.all(np.isfinite(prediction)):
                raise ValueError(f'Invalid prediction for {name}')
            z = np.linalg.solve(whitening, prediction - central)
            chi2 = float(z @ z)
            # Validate against pipeline's own fit chi2 to prevent silently changing objective.
            if not np.isclose(chi2, float(fit.chi2), rtol=2e-5, atol=1e-5):
                raise ValueError(f'Chi2 covariance mismatch for {name}: reconstructed={chi2}, pipeline={fit.chi2}')
            residuals.extend(z / math.sqrt(len(names)))
            rows.append({'model': name, 'status': 'Success', 'chi2': chi2,
                         **dict(zip(OBS, map(float, prediction)))})
        r = np.asarray(residuals, dtype=float)
        loss = float(r @ r)
        if loss < best['loss']:
            best.update(loss=loss, x=np.asarray(x).copy(), rows=deepcopy(rows))
            print(f'call {calls}: NEW BEST mean chi2={loss:.8f}; elapsed={time.monotonic()-start:.1f}s', flush=True)
        elif calls % 10 == 0:
            print(f'call {calls}: mean chi2={loss:.8f}; best={best["loss"]:.8f}; elapsed={time.monotonic()-start:.1f}s', flush=True)
        return r

    initial = evaluate_vector(x0)
    initial_loss = float(initial @ initial)
    if not np.isclose(initial_loss, baseline_chi2, rtol=2e-5, atol=1e-4):
        raise ValueError(f'Checkpoint mismatch: reconstructed mean chi2 {initial_loss} vs {baseline_chi2}')
    # Create new folder only after the initial objective is validated.
    OUTPUT_DIR.mkdir(parents=True, exist_ok=False)
    save_json_new(OUTPUT_DIR / 'initial_validation.json', {
        'checkpoint_mean_chi2': baseline_chi2, 'recomputed_mean_chi2': initial_loss,
        'model_count': len(names), 'yukawa_bound': BOUND})
    try:
        result = least_squares(evaluate_vector, x0, bounds=(-np.ones(18), np.ones(18)),
                               method='trf', jac='2-point', diff_step=args.diff_step,
                               max_nfev=args.max_nfev, verbose=1,
                               x_scale='jac', ftol=1e-7, xtol=1e-7, gtol=1e-7)
    except Exception as exc:
        # Preserve best point if solver encounters an unexpected numerical error.
        save_json_new(OUTPUT_DIR / 'solver_error.json', {'error': f'{type(exc).__name__}: {exc}',
                     'best_mean_chi2': best['loss'], 'best_x': best['x'].tolist(), 'calls': calls})
        raise
    # Prefer best observed point over solver's final iterate.
    best_x = best['x']
    common_out = deepcopy(common)
    common_out['shared_ordinary']['y1']['real'] = (best_x[:9]*BOUND).reshape(3,3).tolist()
    common_out['shared_ordinary']['y2']['real'] = (best_x[9:]*BOUND).reshape(3,3).tolist()
    common_out['average_chi2'] = best['loss']
    common_out['objective'] = 'bounded least squares of 80 NuFIT covariance-whitened residuals; fixed masses and lambdaT3'
    save_json_new(OUTPUT_DIR / 'optimized_common_input.json', common_out)
    save_json_new(OUTPUT_DIR / 'optimization_summary.json', {
        'initial_average_chi2': initial_loss, 'best_average_chi2': best['loss'],
        'solver_final_average_chi2': float(2*result.cost), 'success': bool(result.success),
        'message': str(result.message), 'solver_nfev': int(result.nfev),
        'actual_residual_calls': calls, 'elapsed_seconds': time.monotonic()-start,
        'model_count': len(names), 'parameter_count': 18,
        'max_nfev': args.max_nfev, 'diff_step': args.diff_step,
        'source_checkpoint': str(SOURCE_DIR / 'checkpoint.json')})
    save_json_new(OUTPUT_DIR / 'optimized_predictions.json', {'average_chi2': best['loss'], 'models': best['rows']})
    with (OUTPUT_DIR / 'optimized_predictions.csv').open('x', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=['model', 'status', 'chi2', *OBS])
        writer.writeheader(); writer.writerows(best['rows'])
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        sigmas = np.sqrt(np.diag(cov))
        for j, obs in enumerate(OBS):
            fig, ax = plt.subplots(figsize=(12, 5))
            ax.axhspan(central[j]-sigmas[j], central[j]+sigmas[j], alpha=.18, label='NuFIT Gaussian 1 sigma')
            ax.axhline(central[j], color='black', linestyle='--', label='NuFIT central')
            ax.scatter(range(len(names)), [r[obs] for r in best['rows']], s=30, label='Common Yukawa')
            ax.set_xticks(range(len(names)), names, rotation=90, fontsize=7)
            ax.set_ylabel(obs); ax.legend(); ax.grid(axis='y', alpha=.2)
            fig.tight_layout(); fig.savefig(OUTPUT_DIR / f'{obs}_least_squares.png', dpi=160); plt.close(fig)
    except ImportError:
        print('matplotlib unavailable; skipped plots', flush=True)
    print('Finished. Initial mean chi2:', initial_loss, 'Best mean chi2:', best['loss'])
    print('Output:', OUTPUT_DIR)
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
