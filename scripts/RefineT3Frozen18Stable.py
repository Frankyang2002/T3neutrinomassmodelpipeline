"""Re-fit saved frozen-UV general-quartic T3 Yukawas with a controlled 18D trust region.

This is an independent local refinement, not a replacement for any physics module.
Run from the repository root: python -m scripts.RefineT3Frozen18Stable ...
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.optimize import least_squares
from scripts import RunT3GeneralScalar48Campaign as bulk

CASES = ('small_uv_frozen_5e-9', 'large_uv_frozen_0p5')
FORMAT = 't3-frozen-stable-full18-v1'


def write_json(path: Path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    temporary.replace(path)


def eligible_candidate(row):
    if not isinstance(row, dict) or row.get('chi2') is None:
        return False
    try:
        return (math.isfinite(float(row['chi2'])) and
                all(math.isfinite(float(row['parameters'][name])) for name in bulk.NAMES))
    except (KeyError, TypeError, ValueError):
        return False


def load_seeds(model, case, source_campaign, from_runs):
    rows = []
    source = bulk.OUT / source_campaign / case / bulk.safe(model) / 'fit.json'
    if source.exists():
        saved = bulk.read(source)
        if saved.get('model') != model or saved.get('case') != case:
            raise ValueError('Source fit model/case mismatch')
        row = saved.get('best_mass_safe')
        if eligible_candidate(row):
            rows.append((row, str(source)))
    for run in from_runs:
        folder = bulk.OUT / run / case / bulk.safe(model)
        if not folder.is_dir():
            raise FileNotFoundError(f'Saved refinement folder not found: {folder}')
        for path in sorted(folder.glob('*.json')):
            saved = bulk.read(path)
            if saved.get('model') != model or saved.get('case') != case:
                continue
            candidates = list(saved.get('retained') or [])
            if isinstance(saved.get('best'), dict):
                candidates.append(saved['best'])
            for item in candidates:
                if eligible_candidate(item):
                    rows.append((item, str(path)))
    unique = {}
    for row, origin in rows:
        key = tuple(float(row['parameters'][name]) for name in bulk.NAMES)
        if key not in unique or row['chi2'] < unique[key][0]['chi2']:
            unique[key] = (row, origin)
    return sorted(unique.values(), key=lambda pair: float(pair[0]['chi2']))


def run(args):
    if args.case not in CASES:
        raise ValueError('Only frozen-UV cases are supported')
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,63}', args.run_id):
        raise ValueError('Invalid run ID')
    if args.nfev < 20 or not (0 < args.radius <= 0.1) or not (0 < args.diff_step <= 0.01):
        raise ValueError('Invalid refinement budget/radius/diff_step')
    if args.seeds < 1 or args.seeds > 20:
        raise ValueError('Invalid seed count')
    reg = bulk.read(bulk.REGISTRY)
    entry, new, reused = bulk.inventory(reg, args.model)
    if reused:
        raise ValueError('B(-1)/C(-1) need their special mixed-scalar evaluator')
    history = bulk.historical(args.case, args.model, entry)
    if history['uv_on']:
        raise ValueError('Frozen-only fitting requested for running-UV configuration')
    beta_path = bulk.OUT / args.source_campaign / 'betas' / (bulk.safe(args.model) + '.json')
    payload = bulk.read(beta_path)
    kernel, match_path, match_kind = bulk.matching(args.model)
    from Numerical.orchestration import RunSmallLambdaRefit as refit
    from Numerical.fitting.OscillationFit import OscillationFitTarget
    target = OscillationFitTarget.from_json(refit.TARGET)
    evaluator = bulk.Evaluator(args.case, args.model, entry, new, history, payload,
                                kernel, target, args.extra_real)
    lo, hi = bulk.b21.search_box(history['seed'], local_fraction=1.5,
          floor_fraction=.05, abs_yukawa_max=math.sqrt(4 * math.pi))
    width = hi - lo
    output = bulk.OUT / args.run_id / args.case / bulk.safe(args.model)
    output.mkdir(parents=True, exist_ok=True)
    state_path = output / 'state.json'
    trace_path = output / 'evaluations.jsonl'
    if state_path.exists():
        raise FileExistsError('Run already exists. Use a NEW --run-id to preserve all results.')
    seeds = load_seeds(args.model, args.case, args.source_campaign, args.from_run)
    if not seeds:
        raise ValueError('No valid saved seed found')
    used = seeds[:args.seeds]
    best = None
    stages = []
    errors = 0
    def record(row, name):
        nonlocal best
        with trace_path.open('a', encoding='utf-8') as f:
            f.write(json.dumps({'stage': name, 'status': row.get('status'),
                'chi2': row.get('chi2'), 'reason': row.get('reason')}, allow_nan=False) + '\n')
        if eligible_candidate(row) and (best is None or row['chi2'] < best['chi2']):
            best = dict(row)
            write_json(state_path, {'format': FORMAT, 'model': args.model, 'case': args.case,
                'status': 'Running', 'best': best, 'stages': stages, 'inputs': inputs})
            print(f'NEW BEST chi2={best["chi2"]:.9g}', flush=True)
    inputs = {'source_campaign': args.source_campaign, 'from_runs': args.from_run,
              'beta_sha256': bulk.sha(beta_path), 'matching_sha256': bulk.sha(match_path),
              'historical_sha256': bulk.sha(history['cfg_path']),
              'radius_unit': args.radius, 'diff_step_unit': args.diff_step,
              'max_nfev_per_stage': args.nfev, 'target_chi2': args.target,
              'seed_count': len(used)}
    for seed_index, (seed, origin) in enumerate(used, 1):
        if best is not None and best['chi2'] <= args.target:
            break
        seed_x = np.array([seed['parameters'][name] for name in bulk.NAMES])
        # Recompute each imported seed under the exact current solver inputs.
        baseline = evaluator.evaluate(seed_x, f'seed_{seed_index}')
        record(baseline, f'seed_{seed_index}')
        if baseline.get('chi2') is None:
            stages.append({'seed': seed_index, 'status': 'InvalidSeed', 'reason':baseline.get('reason')})
            continue
        seed_z = np.clip((seed_x - lo) / width, 0, 1)
        # Constrain the trust-region's available parameter range to the *local*
        # neighborhood of the seed.  This is different from the earlier global
        # [0,1] search box and prevents giant finite-difference excursions.
        lower = np.maximum(0, seed_z - args.radius)
        upper = np.minimum(1, seed_z + args.radius)
        lower = np.minimum(lower, seed_z - 1e-12)
        upper = np.maximum(upper, seed_z + 1e-12)
        lower = np.maximum(0, lower)
        upper = np.minimum(1, upper)
        # Strictly positive widths are required by scipy.
        if np.any(upper <= lower):
            raise ValueError('Degenerate local parameter bounds')
        def residual(z):
            nonlocal errors
            trial = evaluator.evaluate(lo + width * z, f'seed_{seed_index}_local')
            record(trial, f'seed_{seed_index}_local')
            if trial.get('chi2') is None:
                errors += 1
                return np.full(5, 1e5)
            v = np.asarray(trial.get('whitened_residual'), dtype=float)
            if v.shape != (5,) or not np.isfinite(v).all():
                errors += 1
                return np.full(5, 1e5)
            return v
        # NOTE: Explicit finite-difference increment in unit coordinates is
        # supplied through a user-controlled diff_step.  x_scale='jac' allows
        # a trust region matched to the residual sensitivities.
        result = least_squares(residual, seed_z, bounds=(lower, upper), method='trf',
                               jac='3-point', diff_step=args.diff_step,
                               x_scale='jac', max_nfev=args.nfev,
                               ftol=1e-12, xtol=1e-12, gtol=1e-12,
                               verbose=0)
        stage = {'seed': seed_index, 'source': origin,
                 'baseline_chi2': float(baseline['chi2']),
                 'optimizer_nfev': int(result.nfev), 'optimizer_njev':int(result.njev) if result.njev is not None else None,
                 'optimizer_status': int(result.status), 'optimizer_success':bool(result.success),
                 'optimizer_message':str(result.message), 'optimality': float(result.optimality),
                 'cost': float(result.cost), 'initial_radius_unit':args.radius,
                 'error_evaluations':errors,
                 'best_chi2':float(best['chi2']) if best else None}
        stages.append(stage)
        write_json(state_path, {'format':FORMAT, 'model':args.model,'case':args.case,
                                'status':'Running','best':best,'stages':stages,'inputs':inputs})
        print(f'SEED {seed_index}/{len(used)} nfev={result.nfev} status={result.status} '
              f'best_chi2={None if best is None else best["chi2"]} message={result.message}', flush=True)
    final_status = ('CandidateConstrainedAccepted' if best and best['chi2'] <= args.target
                    else 'CandidateConstrainedPoorFit' if best else 'NoMassSafeNeutrinoFit')
    if best:
        config = bulk.update_yukawas(history,np.array([best['parameters'][n] for n in bulk.NAMES]))
        config['conditional_general_scalar_fit'] = {
            'model':args.model, 'case':args.case,'chi2':best['chi2'],
            'source':'stable-local-full18', 'sequential_matching_certified':False}
        write_json(output/'fitted_config.json',config)
    write_json(state_path, {'format':FORMAT,'model':args.model,'case':args.case,
                            'status':final_status,'best':best,'stages':stages,'inputs':inputs,
                            'evaluations':evaluator.count,'numerical_failures':evaluator.failures})
    print(f'COMPLETE {args.model} {args.case} {final_status} '
          f'chi2={None if best is None else best["chi2"]}', flush=True)
    return 0 if final_status=='CandidateConstrainedAccepted' else 2


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-id',required=True)
    p.add_argument('--model',default='A_+0')
    p.add_argument('--case',default='small_uv_frozen_5e-9')
    p.add_argument('--source-campaign',default='overnight_six_v1_bulk48')
    p.add_argument('--from-run',action='append',default=[],help='May be specified repeatedly')
    p.add_argument('--seeds',type=int,default=3)
    p.add_argument('--nfev',type=int,default=600)
    p.add_argument('--radius',type=float,default=0.002,
                   help='Half-width of local optimization window in normalized coordinates')
    p.add_argument('--diff-step',type=float,default=1e-6)
    p.add_argument('--extra-real',type=float,default=-0.01)
    p.add_argument('--target',type=float,default=5.86)
    return run(p.parse_args(argv))


if __name__ == '__main__':
    raise SystemExit(main())
