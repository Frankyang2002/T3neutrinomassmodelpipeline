"""Prepare a reproducible thesis numerical input and optionally run the EXISTING pipeline.

Defaults to the repository's saved NuFIT-fitted reference point, *not* a new
physics model. --refit enables the pre-existing Sobol/local refinement path.

Run from project root:
  python -X utf8 -m scripts.PrepareThesisNumericPilot --model A_-2
  python -X utf8 -m scripts.PrepareThesisNumericPilot --model A_-2 --run
  python -X utf8 -m scripts.PrepareThesisNumericPilot --model A_-2 --refit --run
"""
from __future__ import annotations

import argparse
from copy import deepcopy
import json
import math
from pathlib import Path
import subprocess
import sys
import shutil

ROOT = Path(__file__).resolve().parents[1]
REPS = {'A': (1, 3, 2), 'B': (2, 2, 1), 'C': (2, 2, 3),
        'D': (3, 1, 2), 'E': (3, 3, 2)}
ALPHAS = {'A': {-4, -2, 0}, 'B': {-3, -1, 1}, 'C': {-3, -1, 1},
          'D': {-2, 0, 2}, 'E': {-4, -2, 0, 2}}
TARGET = ROOT / 'configs' / 'nufit_6_1_ic24_no.json'


def model_key(model: str) -> str:
    cls, atext = model.split('_', 1)
    alpha = int(atext)
    if cls not in REPS or alpha not in ALPHAS[cls]:
        raise ValueError(f'Not one of the 16 ordinary T3 models: {model}')
    d1, d2, df = REPS[cls]
    suffix = ('m' if alpha < 0 else 'p') + str(abs(alpha))
    return f'T3_dS1_{d1}_dS2_{d2}_dF_{df}_alpha_{suffix}'


def _finite_positive(value: float, name: str) -> float:
    value = float(value)
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f'{name} must be a finite positive number')
    return value


def make_config(reference: dict, *, lam: float | None = None,
                yukawa_factor: float = 1.0, gamma_min: float = 0.8,
                gamma_max: float = 1.2, gamma_points: int = 11,
                refit: bool = False, sobol_points: int = 64) -> dict:
    """Copy a complete known-good fitted config; never rewrite the reference."""
    factor = _finite_positive(yukawa_factor, 'yukawa_factor')
    gmin = _finite_positive(gamma_min, 'gamma_min')
    gmax = _finite_positive(gamma_max, 'gamma_max')
    if not (gmin < 1 < gmax) or gamma_points < 3 or gamma_points % 2 != 1:
        raise ValueError('Gamma range must straddle 1 and points must be odd >=3')
    if sobol_points < 8 or sobol_points & (sobol_points - 1):
        raise ValueError('Sobol points must be a power of two >= 8')
    data = deepcopy(reference)
    if data.get('kind') != 't3_pipeline_numerical_results_v1':
        raise ValueError('Source is not a pipeline numerical-results config')
    ordinary = data['base_state']['ordinary']
    for name in ('y1', 'y2'):
        for field in ('real', 'imag'):
            if field in ordinary[name]:
                ordinary[name][field] = [[float(v) * factor for v in row]
                                         for row in ordinary[name][field]]
    if lam is not None:
        ordinary['lambdaT3']['real'] = _finite_positive(lam, 'lambda_t3')
        ordinary['lambdaT3']['imag'] = 0.0
    data['ordering'] = 'NO'
    data['sensitivity'] = {
        'enabled': True, 'gamma_min': gmin, 'gamma_max': gmax,
        'n_points': int(gamma_points), 'spacing': 'linear',
        # An ABSOLUTE path avoids the current relative-path resolution issue
        # in Numerical/orchestration/PipelineNumericalResults.py.
        'oscillation_target': str(TARGET.resolve()),
    }
    if refit:
        parameters = [dict(name='lambdaT3_real', binding='ordinary.lambdaT3.real',
                           low=0.5 * max(1e-5, ordinary['lambdaT3']['real']),
                           high=2.0 * max(1e-5, ordinary['lambdaT3']['real']), scale='log')]
        for name in ('y1', 'y2'):
            for i, row in enumerate(ordinary[name]['real']):
                for j, value in enumerate(row):
                    v = float(value)
                    radius = max(abs(v) * 0.8, 2.0e-5)
                    parameters.append(dict(name=f'{name}_{i+1}{j+1}',
                                           binding=f'ordinary.{name}.real.{i}.{j}',
                                           low=v-radius, high=v+radius, scale='linear'))
        data['benchmark_search'] = dict(
            enabled=True, method='sobol', use_current_model_representation=True,
            n_points=int(sobol_points), target_chi2=5.0, keep_best=8, seed=20261101,
            oscillation_target=str(TARGET.resolve()),
            local_refinement=dict(enabled=True, seed_count=3,
                                  active_counts=[5, 8, 12],
                                  max_nfev=[25, 35, 45], sensitivity_step=0.02),
            parameters=parameters)
    else:
        data['benchmark_search'] = {'enabled': False}
    data['_thesis_numeric_pilot'] = {
        'source': 'Numerical/benchmarks/fitted_reference',
        'purpose': 'exploratory numerical running; reference targets are fitted, not predictions',
        'yukawa_factor_relative_to_reference': factor,
        'refit': bool(refit),
        'note': 'Uses original restricted matching/numerical pipeline, not full-general C5 matching',
    }
    return data


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--model', default='A_-2')
    ap.add_argument('--lambda-t3', type=float, default=None)
    ap.add_argument('--yukawa-factor', type=float, default=1.)
    ap.add_argument('--gamma-min', type=float, default=0.8)
    ap.add_argument('--gamma-max', type=float, default=1.2)
    ap.add_argument('--gamma-points', type=int, default=11)
    ap.add_argument('--refit', action='store_true', help='Use the EXISTING Sobol/local fitter; can be expensive')
    ap.add_argument('--sobol-points', type=int, default=64)
    ap.add_argument('--run', action='store_true', help='Run the original pipeline after creating the config')
    ap.add_argument('--output', type=Path)
    ap.add_argument('--study', default=None)
    args = ap.parse_args(argv)
    key = model_key(args.model)
    cls, aval = args.model.split('_', 1)
    alpha = int(aval)
    source = ROOT / 'Numerical' / 'benchmarks' / 'fitted_reference' / f'{key}.json'
    if not source.is_file():
        ap.error(f'Missing fitted reference: {source}')
    if not TARGET.is_file():
        ap.error(f'Missing NuFIT target: {TARGET}')
    ref = json.loads(source.read_text(encoding='utf-8-sig'))
    try:
        cfg = make_config(ref, lam=args.lambda_t3, yukawa_factor=args.yukawa_factor,
                          gamma_min=args.gamma_min, gamma_max=args.gamma_max,
                          gamma_points=args.gamma_points, refit=args.refit,
                          sobol_points=args.sobol_points)
    except (ValueError, KeyError, TypeError) as exc:
        ap.error(str(exc))
    suffix = '_refit' if args.refit else '_reference'
    dest = args.output or ROOT / 'configs' / 'thesis_numerics' / (key + suffix + '.json')
    dest = dest.resolve()
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.is_file():
        ap.error(f'Refusing to overwrite user config: {dest}. Use --output with a new name.')
    dest.write_text(json.dumps(cfg, indent=2) + '\n', encoding='utf-8')
    print('Prepared numerical input:', dest)
    print('Reference:', source)
    study = args.study or ('thesis_' + cls + '_' + ('m' if alpha < 0 else 'p') + str(abs(alpha)))
    d1, d2, df = REPS[cls]
    command = [sys.executable, '-X', 'utf8', 'pipeline.py',
               '--dims', str(d1), str(d2), str(df), '--alpha', str(alpha),
               '--restricted-production', '--threshold', 'F', '--threshold', 'S1', 'S2',
               '--study', study, '--numerical', str(dest), '--reset-numerical-configs']
    print('Existing pipeline command:')
    print(subprocess.list2cmdline(command))
    if args.run:
        # The existing pipeline's --reset-numerical-configs replaces one
        # configs/generated_models file; save its previous content first.
        generated = ROOT / 'configs' / 'generated_models' / f'{key}.json'
        if generated.is_file():
            backup = dest.parent / f'{key}_before_{study}_generated_backup.json'
            if backup.exists():
                ap.error(f'Backup already exists; refusing to overwrite: {backup}')
            shutil.copy2(generated, backup)
            print('Preserved previous generated model config:', backup)
        status = subprocess.run(command, cwd=ROOT, check=False).returncode
        if status:
            return status
        print('Build interactive reports:')
        print(f'{sys.executable} -X utf8 -m scripts.BuildThesisNumericHTML --study {study}')
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
