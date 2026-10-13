"""Bounded common-Yukawa SPSA search across 16 T3 representations.

Read-only inputs. Writes ONLY to a new output directory; checkpoint after every iteration.
Usage: python -m Numerical.orchestration.OptimizeCommonYukawas --run --iterations 60
"""
from __future__ import annotations
import argparse
import csv
import json
import math
import time
from copy import deepcopy
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from Numerical.orchestration.RunCommonParameterComparisonV2 import (
    PROJECT_ROOT, REF, SCENARIO, SHARED, OBS, candidate_for_target,
    models_and_cases, load, evaluate,
)
from Numerical.fitting.OscillationFit import OscillationFitTarget
from Numerical.core.RGBetaEvaluator import load_rgbeta_payload
from Numerical.running.FinalC5TrajectoryAdapter import FinalC5TrajectoryEvaluator
from Numerical.fitting.ParameterScan import make_t3_point_evaluator
from Numerical.fitting.ScanCLI import build_uv_state_from_config
from types import SimpleNamespace

OUT = PROJECT_ROOT / 'Reports' / 'output' / 'full' / 'common_yukawa_optimization' / SCENARIO
SOURCE = PROJECT_ROOT / 'Reports' / 'output' / 'full' / 'common_parameter_comparison_v2' / SCENARIO / 'best_common_input.json'
BOUND = 1e-3


def make_evaluator(payload, case, target):
    data = case / 'data'
    uv = load_rgbeta_payload(data / 'uv_rgbeta_rge.json')
    eft = load_rgbeta_payload(data / 'eft1_rgbeta_rge.json')
    c5 = FinalC5TrajectoryEvaluator(data / 'final_weinberg_coefficient.json')
    scales = payload['scales']; numerical = payload.get('numerical', {})
    view = SimpleNamespace(raw={**payload, 'scan': {'bindings': {}}})
    evaluator = make_t3_point_evaluator(
        uv_state_builder=lambda _params: build_uv_state_from_config(view, {}),
        c5_builder=c5, uv_rgbeta_payload=uv, eft1_rgbeta_payload=eft,
        target=target, mu_fermion_threshold_gev=float(scales['mu_fermion_threshold_gev']),
        mu_scalar_threshold_gev=float(scales['mu_scalar_threshold_gev']),
        mu_low_gev=float(scales['mu_low_gev']),
        vev_gev=float(numerical.get('vev_gev', 246.22)),
        rtol=float(numerical.get('rtol', 1e-8)),
        atol=float(numerical.get('atol', 1e-11)),
    )
    return evaluator


def atomic_json(path, value):
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    tmp.replace(path)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--run', action='store_true')
    ap.add_argument('--iterations', type=int, default=60)
    ap.add_argument('--seed', type=int, default=20261013)
    ap.add_argument('--resume', action='store_true', help='Resume checkpoint in existing output directory')
    args = ap.parse_args(argv)
    if args.iterations < 1: ap.error('--iterations must be >= 1')
    configs, cases = models_and_cases()
    names = sorted(configs)
    common = load(SOURCE)
    source_name = common['source_model']
    source = load(configs[source_name])
    target = OscillationFitTarget.from_json(PROJECT_ROOT / 'configs' / 'nufit_6_1_ic24_no.json')
    targets = [candidate_for_target(source, load(configs[n])) for n in names]
    print('Models:', len(names), 'start:', source_name, 'output:', OUT, flush=True)
    print('Iterations:', args.iterations, 'approx evaluations:', 1 + 3 * args.iterations, 'each with 16 model trajectories', flush=True)
    if not args.run:
        print('Dry run. Use --run; use --resume to continue an existing run.')
        return 0
    if args.resume:
        if not (OUT / 'checkpoint.json').is_file():
            raise FileNotFoundError('No checkpoint to resume: ' + str(OUT))
    else:
        OUT.mkdir(parents=True, exist_ok=False)
    evaluators = [make_evaluator(p, cases[n], target) for p, n in zip(targets, names)]
    original = np.r_[np.asarray(source['base_state']['ordinary']['y1']['real'], float).ravel(),
                     np.asarray(source['base_state']['ordinary']['y2']['real'], float).ravel()] / BOUND
    if len(original) != 18: raise ValueError('Expected two real 3x3 Yukawa matrices')
    if np.max(np.abs(original)) > 1: raise ValueError('Starting Yukawas exceed bound')
    rng = np.random.default_rng(args.seed)
    def evaluate_all(x):
        x = np.clip(np.asarray(x, float), -1, 1)
        y1, y2 = (x[:9] * BOUND).reshape(3, 3), (x[9:] * BOUND).reshape(3, 3)
        for p in targets:
            p['base_state']['ordinary']['y1']['real'] = y1.tolist()
            p['base_state']['ordinary']['y2']['real'] = y2.tolist()
        rows = []
        for name, ev in zip(names, evaluators):
            try:
                fit = ev({})
                chi2 = float(fit.chi2)
                if not math.isfinite(chi2): raise ValueError('Non-finite chi2')
                rows.append({'model': name, 'status': 'Success', 'chi2': chi2,
                             'prediction': dict(zip(fit.observable_names, map(float, fit.prediction)))})
            except Exception as exc:
                rows.append({'model': name, 'status': 'Failed', 'chi2': None,
                             'error': f'{type(exc).__name__}: {exc}'})
        if any(r['status'] != 'Success' for r in rows):
            return 1e12, rows
        return sum(r['chi2'] for r in rows) / len(rows), rows
    if args.resume:
        cp = load(OUT / 'checkpoint.json')
        if cp['seed'] != args.seed: raise ValueError('Seed mismatch; use original --seed')
        x = np.asarray(cp['current_x'], float); best_x = np.asarray(cp['best_x'], float)
        best_loss = float(cp['best_average_chi2']); start = int(cp['completed_iterations'])
        rng.bit_generator.state = cp['rng_state']
    else:
        x = original.copy(); best_x = x.copy(); start = 0
        best_loss, rows = evaluate_all(x)
        if best_loss >= 1e12: raise RuntimeError('Initial common candidate failed; see source cross-evaluation')
        atomic_json(OUT / 'initial_evaluation.json', {'average_chi2': best_loss, 'rows': rows})
        print(f'Initial mean chi2 = {best_loss:.8g}', flush=True)
    start_time = time.monotonic()
    for k in range(start, args.iterations):
        # Simultaneous perturbation: three 16-model evaluations per iteration, independent of 18 dimensions.
        a = 0.025 / (1 + k / 20) ** 0.602
        c = 0.012 / (1 + k / 20) ** 0.101
        delta = rng.choice([-1.0, 1.0], size=18)
        plus, _ = evaluate_all(np.clip(x + c * delta, -1, 1))
        minus, _ = evaluate_all(np.clip(x - c * delta, -1, 1))
        # Normalize gradient to avoid giant jumps when chi2 is large.
        sign = np.sign(plus - minus)
        proposal = np.clip(x - a * sign * delta, -1, 1)
        # Reuse the better perturbed point if it beats the running optimum.
        if plus < best_loss:
            best_loss, best_x = plus, np.clip(x + c * delta, -1, 1)
        if minus < best_loss:
            best_loss, best_x = minus, np.clip(x - c * delta, -1, 1)
        proposal_loss, _ = evaluate_all(proposal)
        if proposal_loss < best_loss:
            best_loss, best_x = proposal_loss, proposal.copy()
        # Always retain the best point as the centre for the next iteration.
        x = best_x.copy()
        cp = {'seed': args.seed, 'completed_iterations': k + 1,
              'current_x': x.tolist(), 'best_x': best_x.tolist(),
              'best_average_chi2': best_loss, 'rng_state': rng.bit_generator.state,
              'source_candidate': source_name, 'yukawa_bound': BOUND,
              'elapsed_seconds_this_session': time.monotonic() - start_time}
        atomic_json(OUT / 'checkpoint.json', cp)
        print(f'iteration {k+1}/{args.iterations}: best mean chi2={best_loss:.8g}; elapsed={time.monotonic()-start_time:.1f}s', flush=True)
    final_loss, rows = evaluate_all(best_x)
    result = deepcopy(common)
    result['shared_ordinary']['y1']['real'] = (best_x[:9] * BOUND).reshape(3, 3).tolist()
    result['shared_ordinary']['y2']['real'] = (best_x[9:] * BOUND).reshape(3, 3).tolist()
    result['objective'] = 'minimum observed mean chi2 over 16 models; fixed masses and lambdaT3'
    result['average_chi2'] = final_loss
    atomic_json(OUT / 'optimized_common_input.json', result)
    atomic_json(OUT / 'optimized_predictions.json', {'average_chi2': final_loss, 'models': rows})
    with (OUT / f'optimized_predictions_iteration_{args.iterations}.csv').open('x', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=['model', 'status', 'chi2', *OBS, 'error'])
        writer.writeheader()
        for r in rows: writer.writerow({'model': r['model'], 'status': r['status'], 'chi2': r.get('chi2'),
            **r.get('prediction', {}), 'error': r.get('error', '')})
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        center = dict(zip(target.observable_names, map(float, target.central_values)))
        sigma = dict(zip(target.observable_names, np.sqrt(np.diag(target.covariance))))
        for obs in OBS:
            fig, ax = plt.subplots(figsize=(12, 5))
            ax.axhspan(center[obs]-sigma[obs], center[obs]+sigma[obs], alpha=.18, label='NuFIT Gaussian 1 sigma')
            ax.axhline(center[obs], color='black', linestyle='--', label='NuFIT central')
            for i, r in enumerate(rows):
                if r['status'] == 'Success': ax.scatter(i, r['prediction'][obs], color='C0')
            ax.set_xticks(range(16), names, rotation=90, fontsize=7)
            ax.set_ylabel(obs); ax.legend(); ax.grid(axis='y', alpha=.2)
            fig.tight_layout(); fig.savefig(OUT / f'{obs}_optimized.png', dpi=160); plt.close(fig)
    except ImportError:
        print('matplotlib missing; skipping plots')
    print('Finished. Best mean chi2:', final_loss, 'Output:', OUT, flush=True)
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
