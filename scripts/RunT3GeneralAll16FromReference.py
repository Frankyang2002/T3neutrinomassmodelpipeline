"""Exploratory, resumable 16-model general-RGE / GENERAL-MATCHING Yukawa fitter.

For B(-1)/C(-1) reuse the existing mixed-scalar one-scale evaluator.
For the other 14, evolve canonical full-general UV/EFT1 RGEs and evaluate
full-general Matchete C5, retaining each allowed Yukawa-pair flavour channel.
This remains one-scale conditional matching, not sequential EFT matching.

No historical 1e-7/5e-9 fitted directories or old 48-study campaign are required.
All outputs are new; no GitHub or existing physical source files are modified.

Examples:
  python -X utf8 -m scripts.RunT3GeneralAll16FromReference --run-id thesis16_v1 --preflight
  python -X utf8 -m scripts.RunT3GeneralAll16FromReference --run-id thesis16_v1 --sobol-points 128 --export-bc-betas
  python -X utf8 -m scripts.RunT3GeneralAll16FromReference --run-id thesis16_v1 --resume --sobol-points 128 --export-bc-betas
"""
from __future__ import annotations

import argparse
import csv
from copy import deepcopy
from hashlib import sha256
import json
import math
from pathlib import Path
import re
import traceback
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT / 'Numerical' / 'benchmarks' / 'fitted_reference'
OUTPUT = ROOT / 'output' / 'studies' / 'general_scalar_extended' / 'thesis_general_all16'
REP = {'A': (1, 3, 2), 'B': (2, 2, 1), 'C': (2, 2, 3),
       'D': (3, 1, 2), 'E': (3, 3, 2)}
ALPHAS = {'A': (-4, -2, 0), 'B': (-3, -1, 1), 'C': (-3, -1, 1),
          'D': (-2, 0, 2), 'E': (-4, -2, 0, 2)}
MODELS = tuple(f'{family}_{alpha:+d}' for family in ALPHAS for alpha in ALPHAS[family])
MIXED = {'B_-1', 'C_-1'}
SELF = {'A_-2', 'A_+0', 'D_-2', 'D_+0', 'E_-2', 'E_+0'}
NAMES = tuple(f'{name}_{i}{j}' for name in ('y1', 'y2') for i in range(3) for j in range(3))
OBS = ('delta_m21_sq_ev2', 'delta_m3l_sq_ev2', 'sin2_theta12', 'sin2_theta13', 'sin2_theta23')
FORMAT = 't3-thesis-canonical-general-rge-conditional-c5-fit-v1'


def tag(model: str) -> str:
    return model.replace('-', 'm').replace('+', 'p')


def ref_path(model: str) -> Path:
    family, token = model.split('_')
    alpha = int(token)
    d1, d2, df = REP[family]
    return REFERENCE / f'T3_dS1_{d1}_dS2_{d2}_dF_{df}_alpha_{"m" if alpha < 0 else "p"}{abs(alpha)}.json'


def read(path: Path) -> dict:
    if not path.is_file():
        raise FileNotFoundError(str(path))
    obj = json.loads(path.read_text(encoding='utf-8-sig'))
    if not isinstance(obj, dict):
        raise ValueError(f'Expected JSON object in {path}')
    return obj


def save(path: Path, obj: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + '.tmp')
    temp.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    temp.replace(path)


def fingerprint(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def complex_json(value: complex) -> dict:
    z = complex(value)
    return {'real': float(z.real), 'imag': float(z.imag)}


def matrix(raw: dict) -> np.ndarray:
    return np.asarray(raw['real'], dtype=float) + 1j*np.asarray(raw['imag'], dtype=float)


def matrix_json(values: np.ndarray) -> dict:
    a = np.asarray(values, dtype=complex)
    return {'real': a.real.tolist(), 'imag': a.imag.tolist()}


def entry_matrix(values: np.ndarray) -> list:
    # GeneralRunningBenchmark._override_value expects per-entry JSON complex numbers,
    # not the compact {real: matrix, imag: matrix} config representation.
    a = np.asarray(values, dtype=complex)
    if a.shape != (3, 3):
        raise ValueError('Expected 3x3 matrix')
    return [[complex_json(x) for x in row] for row in a]


def seed_from_config(config: dict) -> np.ndarray:
    o = config['base_state']['ordinary']
    arrays = (matrix(o['y1']), matrix(o['y2']))
    if any(a.shape != (3, 3) or np.max(np.abs(a.imag)) > 1e-12 for a in arrays):
        raise ValueError('The first pilot fits 18 real y1/y2 entries only')
    return np.concatenate([a.real.ravel() for a in arrays])


def config_with_yukawas(config: dict, params: dict) -> dict:
    result = deepcopy(config)
    for k, name in enumerate(('y1', 'y2')):
        arr = np.array([float(params[f'{name}_{i}{j}']) for i in range(3)
                        for j in range(3)], dtype=float).reshape(3, 3)
        result['base_state']['ordinary'][name] = matrix_json(arr)
    result['benchmark_search'] = {'enabled': False}
    result['sensitivity'] = {'enabled': False}
    return result


def fixed_overrides(config: dict, uv_betas: dict, *, quartic_fraction: float,
                    mu2_fraction: float, extra_yukawa_fraction: float) -> dict:
    """All and only UV beta keys; preserve tracked ordinary and SM boundaries."""
    from Numerical.core.GeneralCouplingLayout import spec_for
    raw = config['base_state']
    ordinary = raw['ordinary']
    old_lambda = complex(float(ordinary['lambdaT3']['real']),
                         float(ordinary['lambdaT3']['imag']))
    if abs(old_lambda) == 0:
        raise ValueError('Reference lambdaT3 must be nonzero')
    extras = {}
    for name in uv_betas['report_betas']:
        if name in raw['sm']:
            value = raw['sm'][name]
        elif name in ordinary:
            value = ordinary[name]
        elif name in ('y1Conj', 'y2Conj'):
            parent = 'y1' if name == 'y1Conj' else 'y2'
            value = matrix_json(matrix(ordinary[parent]) * extra_yukawa_fraction)
        elif name.startswith('muT3General'):
            scale = min(float(ordinary['mS1Sq']), float(ordinary['mS2Sq']))
            value = complex_json(mu2_fraction * scale)
        elif name.startswith('lambda'):
            value = complex_json(quartic_fraction * abs(old_lambda))
        else:
            raise ValueError(f'Unrecognised full-general UV coupling key: {name}')
        kind = spec_for(name).kind
        if kind == 'matrix':
            if isinstance(value, dict):
                value = entry_matrix(matrix(value))
        elif kind == 'complex':
            if not isinstance(value, dict):
                value = complex_json(value)
        elif kind == 'real':
            if isinstance(value, dict):
                if abs(float(value.get('imag', 0))) > 1e-12:
                    raise ValueError(f'Imaginary real-valued coupling: {name}')
                value = float(value['real'])
        extras[name] = value
    return extras


def override_yukawas(overrides: dict, params: dict) -> dict:
    result = dict(overrides)
    for name in ('y1', 'y2'):
        vals = np.array([[float(params[f'{name}_{i}{j}']) for j in range(3)]
                         for i in range(3)])
        result[name] = entry_matrix(vals)
    return result


def general_coefficient_path(model: str) -> Path:
    return (ROOT / 'output' / 'general_matchete_matching' / tag(model) /
            'c5_coefficient.txt')


def general_coefficient_status(model: str) -> Path:
    return general_coefficient_path(model).with_name('general_matching_status.json')


def ensure_general_matching(model: str, *, generate: bool, wolframscript: str) -> Path:
    """Never substitute restricted production's final_weinberg_coefficient.json."""
    path = general_coefficient_path(model)
    state = general_coefficient_status(model)
    def ready():
        if not path.is_file() or path.stat().st_size == 0 or not state.is_file():
            return False
        status = read(state)
        return (status.get('model') == model and status.get('status') == 'Matched'
                and status.get('model_selection_verified') is True)
    if ready():
        return path
    if not generate:
        raise FileNotFoundError(f'{model}: full-GENERAL Matchete C5 absent: {path}. '
                                'Use --generate-general-c5, not --restricted-production.')
    from scripts.RunT3GeneralMatcheteMatching import run_model
    row = run_model(model, ROOT / 'output' / 'general_matchete_matching',
                    wolframscript, reuse=False)
    if not ready():
        raise RuntimeError(f'{model}: general Matchete matching unavailable: '
                           f'{row}. Inspect {path.parent / "general_matching.log"}')
    return path


def nonmixed_inputs(model: str, canonical_run_id: str, *,
                    generate_general_c5: bool = False,
                    wolframscript: str = 'wolframscript') -> tuple[dict, dict, dict, Any, Path, str]:
    from Numerical.orchestration.GeneralRunningBenchmark import load_canonical_general_exports
    from Numerical.matching.GeneralT3FlavourC5 import GeneralT3FlavourC5
    from scripts.RunT3GeneralSymbolicStudy import expected_yukawa
    uv, eft, provenance = load_canonical_general_exports(model, canonical_run_id)
    source = ensure_general_matching(model, generate=generate_general_c5,
                                     wolframscript=wolframscript)
    kernel = GeneralT3FlavourC5(model, source, wolframscript=wolframscript)
    expected = expected_yukawa(model)
    if expected and expected not in kernel.required_yukawas:
        raise ValueError(f'{model}: general Matchete C5 is missing allowed '
                         f'Yukawa flavour channel {expected}; refusing partial matching')
    return uv, eft, provenance, kernel, source, 'full-general-matchete-one-scale-channels'


class NonMixedEvaluator:
    """Reuse existing canonical general UV/EFT1 solvers and C5/mass/fit routines."""

    def __init__(self, model: str, uv: dict, eft: dict, kernel: Any,
                 fixed: dict, target: Any, *, include_diagnostics: bool = False):
        self.model, self.uv, self.eft, self.kernel = model, uv, eft, kernel
        self.fixed, self.target = fixed, target
        self.calls, self.failures = 0, 0

    def evaluate(self, params: dict, *, detail: bool = False):
        from Numerical.orchestration.GeneralRunningBenchmark import initial_state, _state_values
        from Numerical.running.UVRunner import run_uv_segment
        from Numerical.running.IntermediateScalarRunner import run_intermediate_scalar_segment
        from Numerical.running.FermionThresholdBoundary import build_intermediate_scalar_boundary
        from Numerical.running.IntermediateScalarStateVector import pack_intermediate_scalar_state
        from scripts.RunT3GeneralScalarBCExploratoryFit import low_energy_fit

        _, initial = initial_state(self.uv, mu_uv=1e7, mu_f=1e5, mu_s=1e3,
                                   overrides=override_yukawas(self.fixed, params))
        uv = run_uv_segment(initial, self.uv, 1e5, rtol=1e-7, atol=1e-10)
        boundary = build_intermediate_scalar_boundary(uv.final_state)
        _, eft_layout = pack_intermediate_scalar_state(boundary)
        keys = set(self.eft.get('report_betas', {}))
        if keys != set(eft_layout.names):
            raise ValueError('Full-general EFT1 RGE/state mismatch: ' +
                             f'missing={sorted(keys-set(eft_layout.names))}; ' +
                             f'extra={sorted(set(eft_layout.names)-keys)}')
        eft = run_intermediate_scalar_segment(boundary, self.eft, 1e3,
                                               rtol=1e-7, atol=1e-10)
        if min(eft.final_state.mS1Sq, eft.final_state.mS2Sq) <= 0:
            raise ValueError('Nonpositive scalar mass-squared at EFT1 endpoint')
        c5 = self.kernel.evaluate(uv.final_state, eft.final_state)
        c5 = np.asarray(c5, dtype=complex)
        if c5.shape != (3, 3) or not np.isfinite(c5).all():
            raise ValueError('Nonfinite/incorrect C5 from existing kernel')
        c5 = (c5 + c5.T) / 2
        fit, low = low_energy_fit(c5, eft.final_state, self.target)
        self.calls += 1
        if detail:
            return fit, {
                'uv_at_fermion_threshold': _state_values(uv.final_state),
                'eft1_at_scalar_threshold': _state_values(eft.final_state),
                'uv_nfev': int(uv.nfev), 'eft1_nfev': int(eft.nfev),
                'c5_at_scalar_threshold': matrix_json(c5),
                'c5_at_low_scale': matrix_json(low.K),
                'matching_scope': 'full-general Matchete quadratic Yukawa channels, one-scale',
                'general_c5_diagnostics': self.kernel.last_diagnostics,
            }
        return fit

    def __call__(self, params: dict):
        return self.evaluate(params)


def candidates(model: str, evaluator: Any, seed: np.ndarray, target: Any,
               *, sobol_points: int, rng_seed: int, bound: float,
               sobol_fraction: float, sobol_floor: float, keep: int) -> tuple[list, list, float | None]:
    from Numerical.fitting.SobolBenchmarkSearch import BenchmarkParameter, generate_sobol_parameter_points
    from scripts.RunT3GeneralScalarBCExploratoryFit import result_record
    seed_point = {k: float(v) for k, v in zip(NAMES, seed)}
    widths = np.maximum(np.abs(seed) * sobol_fraction, sobol_floor)
    low = np.maximum(-bound, seed - widths)
    high = np.minimum(bound, seed + widths)
    defs = tuple(BenchmarkParameter(key, f'ordinary.{key.split("_")[0]}.real.{key[-2]}.{key[-1]}',
                                    float(lo), float(hi))
                 for key, lo, hi in zip(NAMES, low, high))
    rows, failures = [], []
    chi_seed = None

    def evaluate(params: dict, source: str) -> None:
        nonlocal chi_seed
        try:
            fit = evaluator(params)
            result = result_record(fit, params, source)
            if result['chi2'] < 0 or not math.isfinite(result['chi2']):
                raise ValueError('Nonfinite or negative chi2')
            if source == 'reference_seed':
                chi_seed = result['chi2']
            rows.append(result)
            rows.sort(key=lambda x: x['chi2'])
            del rows[keep:]
            if source == 'reference_seed' or rows[0]['source'] == source:
                print(f'   {model}: {source}, best chi2={rows[0]["chi2"]:.8g}', flush=True)
        except Exception as exc:
            if len(failures) < 25:
                failures.append({'source': source, 'error': f'{type(exc).__name__}: {exc}'})

    evaluate(seed_point, 'reference_seed')
    points = generate_sobol_parameter_points(defs, sobol_points, seed=rng_seed)
    for i, point in enumerate(points, 1):
        evaluate(point, f'sobol_{i}')
    if not rows:
        raise RuntimeError(f'No valid numerical Sobol/seed points. First errors: {failures[:3]}')
    return rows, failures, chi_seed


def fit_nonmixed(model: str, args: argparse.Namespace, dest: Path) -> dict:
    from Numerical.fitting.OscillationFit import OscillationFitTarget
    from Numerical.fitting.SobolBenchmarkSearch import BenchmarkParameter, refine_sobol_candidates
    from Numerical.orchestration.RunSmallLambdaRefit import TARGET

    reference = read(ref_path(model))
    uv, eft, provenance, kernel, matching_path, matching_kind = nonmixed_inputs(
        model, args.canonical_run_id,
        generate_general_c5=args.generate_general_c5,
        wolframscript=args.wolframscript)
    fixed = fixed_overrides(reference, uv, quartic_fraction=args.quartic_fraction,
                            mu2_fraction=args.mu2_fraction,
                            extra_yukawa_fraction=args.extra_yukawa_fraction)
    target = OscillationFitTarget.from_json(TARGET)
    seed = seed_from_config(reference)
    if np.any(abs(seed) >= args.bound):
        raise ValueError('Reference Yukawa seed exceeds fit bound')
    evaluator = NonMixedEvaluator(model, uv, eft, kernel, fixed, target)
    rows, failures, seed_chi = candidates(
        model, evaluator, seed, target, sobol_points=args.sobol_points,
        rng_seed=args.seed + MODELS.index(model), bound=args.bound,
        sobol_fraction=args.sobol_fraction, sobol_floor=args.sobol_floor,
        keep=max(8, args.seed_count),
    )
    best = rows[0]
    history = {'enabled': False, 'reason': '--no-local-refinement'}
    if not args.no_local_refinement and best['chi2'] > args.target_chi2:
        defs = []
        for k, val in zip(NAMES, seed):
            half = max(abs(val)*args.sobol_fraction, args.sobol_floor)
            defs.append(BenchmarkParameter(k, f'ordinary.{k.split("_")[0]}.real.{k[-2]}.{k[-1]}',
                                           max(-args.bound, float(val-half)),
                                           min(args.bound, float(val+half))))
        try:
            trial, history = refine_sobol_candidates(
                evaluator=evaluator, observable_names=target.observable_names,
                definitions=defs, candidates=rows, target_chi2=args.target_chi2,
                seed_count=min(args.seed_count, len(rows)), active_counts=(5,8,12),
                max_nfev=(20,30,args.max_nfev), sensitivity_step=0.02)
            if trial['chi2'] < best['chi2']:
                best = trial
        except Exception as exc:
            history = {'enabled': True, 'status': 'FailedRetainedBestPoint',
                       'error': f'{type(exc).__name__}: {exc}'}
    elif not args.no_local_refinement:
        history = {'enabled': True, 'skipped': True, 'reason': 'chi2 target already reached'}
    fit, detail = evaluator.evaluate(best['parameters'], detail=True)
    score = float(fit.chi2)
    fitted = config_with_yukawas(reference, best['parameters'])
    fitted['exploratory_general_fit'] = {
        'model': model, 'chi2': score, 'scope': 'canonical full-general renormalisable RGEs and full-general '
        'Matchete one-scale Yukawa-channel C5 (conditional)',
        'physical_sequential_c5_calculated': False,
    }
    dest.mkdir(parents=True, exist_ok=True)
    save(dest / 'fitted_config.json', fitted)
    save(dest / 'general_uv_overrides.json', override_yukawas(fixed, best['parameters']))
    result = {
        'format': FORMAT, 'model': model,
        'status': 'ExploratoryAcceptedFit' if score <= 5.86 else 'ExploratoryPoorFit',
        'chi2_seed': seed_chi, 'chi2_final': score,
        'prediction': dict(zip(fit.observable_names, map(float, fit.prediction))),
        'target': dict(zip(fit.observable_names, map(float, fit.target))),
        'whitened_residual': np.asarray(fit.whitened_residual).tolist(),
        'best_source': best['source'], 'best_parameters': best['parameters'],
        'sobol_attempts': args.sobol_points, 'sobol_failures': failures,
        'local_refinement': history, 'evaluator_calls': evaluator.calls,
        'matching_kind': matching_kind, 'matching_file': str(matching_path),
        'general_symbolic_provenance': provenance,
        'fixed_general_scalar_boundary': {'quartic_fraction_of_lambdaT3': args.quartic_fraction,
                                          'bilinear_fraction_min_mass_squared': args.mu2_fraction,
                                          'conjugate_yukawa_fraction_of_ordinary': args.extra_yukawa_fraction},
        'general_rges_used': True, 'general_matchete_matching_used': True,
        'full_general_sequential_matching_performed': False,
        'self_bilinear_in_c5': 'T31PlusMass' in kernel.symbols or 'T32PlusMass' in kernel.symbols,
        'general_extra_yukawa_direct_c5_channels': bool(kernel.required_yukawas & {'y1Conj','y2Conj'}),
        'matching_channels': kernel.last_diagnostics.get('channels', []),
        'full_general_c5_flavour_pairings_evaluated': True,
        'numerical_details': detail,
        'input_sha256': {'reference': fingerprint(ref_path(model)),
                         'c5_kernel': fingerprint(matching_path)},
    }
    save(dest / 'result.json', result)
    with (dest / 'observables.csv').open('w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(['observable','prediction','target'])
        for k in target.observable_names:
            w.writerow([k, result['prediction'][k], result['target'][k]])
    return result


def preflight(model: str, args: argparse.Namespace) -> dict:
    path = ref_path(model)
    checks = {'reference_config': path.is_file()}
    if model in MIXED:
        from scripts.RunGeneralT3NumericalC5 import DEFAULT_SOURCE
        from scripts.RunT3ZeroExtraBoundaryBenchmark import beta_path
        coeff = DEFAULT_SOURCE / (tag(model) + '_final_one_loop_probe_weinberg_ordered_coefficient.inputform.txt')
        external = args.reuse_b_run if model == 'B_-1' else args.reuse_c_run
        if external:
            from scripts.RunT3BCGeneralFromReference import OUT as bc_out
            saved = bc_out / external / tag(model)
            checks['reusable_previous_fit'] = ((saved / 'result.json').is_file() and
                                                (saved / 'fitted_config.json').is_file())
        else:
            checks['matchete_one_scale_c5'] = coeff.is_file()
            for stage in ('UV','EFT1'):
                checks[f'general_{stage.lower()}_rgbeta'] = (args.export_bc_betas or
                                                                beta_path(model,stage).is_file())
    else:
        from Numerical.orchestration.GeneralRunningBenchmark import load_canonical_general_exports
        try:
            load_canonical_general_exports(model, args.canonical_run_id)
            checks['canonical_full_general_rg_exports'] = True
        except Exception as exc:
            checks['canonical_full_general_rg_exports'] = False
            checks['canonical_error'] = f'{type(exc).__name__}: {exc}'[:400]
        source = general_coefficient_path(model)
        status = general_coefficient_status(model)
        checks['full_general_matchete_c5'] = (source.is_file() and status.is_file()
            and read(status).get('status') == 'Matched'
            and read(status).get('model') == model)
        if not checks['full_general_matchete_c5']:
            checks['expected_general_c5_path'] = str(source)
            checks['generation_command'] = ('python -X utf8 -m '
                'scripts.RunT3GeneralMatcheteMatching --model ' + model)
    passed = all(v for k, v in checks.items() if isinstance(v, bool))
    return {'model': model, 'status': 'Ready' if passed else 'MissingInputs', 'checks': checks}


def summarise(root: Path, rows: list[dict]) -> None:
    root.mkdir(parents=True, exist_ok=True)
    header = ('model','status','chi2_initial','chi2_final','delta_m21_sq_ev2',
              'delta_m3l_sq_ev2','sin2_theta12','sin2_theta13','sin2_theta23',
              'method','output_folder','message')
    with (root / 'summary.csv').open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=header)
        w.writeheader()
        for row in rows:
            r = dict(row)
            pred = r.pop('prediction', {}) or {}
            for k in OBS:
                r[k] = pred.get(k)
            w.writerow({k: r.get(k,'') for k in header})
    save(root / 'summary.json', {'format': FORMAT, 'model_count': len(rows),
                                 'successful_fit_count': sum(r['status'] in
                                    ('ExploratoryAcceptedFit','ExploratoryPoorFit','ReusedFit')
                                    for r in rows),
                                 'rows': rows})


def validate_args(args: argparse.Namespace) -> None:
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,63}', args.run_id):
        raise ValueError('Unsafe run-id')
    if args.sobol_points < 2 or args.sobol_points & (args.sobol_points-1):
        raise ValueError('sobol-points must be a power of 2 >= 2')
    if not (args.max_nfev > 0 and 1 <= args.seed_count <= 8 and
            args.target_chi2 > 0 and args.bound > 0 and
            args.sobol_floor > 0 and args.sobol_fraction > 0 and
            0 <= args.quartic_fraction <= 1 and 0 <= args.mu2_fraction <= .2 and
            0 <= args.extra_yukawa_fraction <= 1):
        raise ValueError('Invalid fit/scalar settings')


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-id', required=True)
    p.add_argument('--models', help='Comma-separated model identifiers (default all 16)')
    p.add_argument('--canonical-run-id', default='full16_general_v1')
    p.add_argument('--preflight', action='store_true', help='Report data prerequisites; no optimisation')
    p.add_argument('--resume', action='store_true')
    p.add_argument('--export-bc-betas', action='store_true')
    p.add_argument('--generate-general-c5', action='store_true',
                   help='Generate 14 full-general Matchete C5 expressions if missing')
    p.add_argument('--wolframscript', default='wolframscript')
    p.add_argument('--reuse-b-run', help='Previously completed B(-1) direct run ID, e.g. thesis_b_m1_128')
    p.add_argument('--reuse-c-run', help='Previously completed C(-1) direct run ID')
    p.add_argument('--sobol-points', type=int, default=128)
    p.add_argument('--max-nfev', type=int, default=60)
    p.add_argument('--seed', type=int, default=20261101)
    p.add_argument('--seed-count', type=int, default=3)
    p.add_argument('--target-chi2', type=float, default=1e-8)
    p.add_argument('--bound', type=float, default=3.0)
    p.add_argument('--sobol-fraction', type=float, default=0.5)
    p.add_argument('--sobol-floor', type=float, default=1e-5)
    p.add_argument('--quartic-fraction', type=float, default=0.1)
    p.add_argument('--mu2-fraction', type=float, default=0.01)
    p.add_argument('--extra-yukawa-fraction', type=float, default=0.05)
    p.add_argument('--no-local-refinement', action='store_true')
    args = p.parse_args(argv)
    validate_args(args)
    wanted = [x.strip() for x in args.models.split(',') if x.strip()] if args.models else list(MODELS)
    if not wanted or len(wanted) != len(set(wanted)) or set(wanted) - set(MODELS):
        p.error('Invalid or repeated --models selection; use A_-2,B_-1,...')
    selected = [x for x in MODELS if x in wanted]
    out = OUTPUT / args.run_id
    if out.exists() and not (args.resume or args.preflight) and (out / 'settings.json').exists():
        p.error('Campaign folder already exists; use --resume or a fresh --run-id')
    out.mkdir(parents=True, exist_ok=True)
    signature = {'models': selected, 'canonical_run_id': args.canonical_run_id,
                 'sobol_points': args.sobol_points, 'max_nfev': args.max_nfev,
                 'seed': args.seed, 'seed_count': args.seed_count,
                 'target_chi2': args.target_chi2, 'bound': args.bound,
                 'sobol_fraction': args.sobol_fraction,
                 'sobol_floor': args.sobol_floor, 'quartic_fraction': args.quartic_fraction,
                 'mu2_fraction': args.mu2_fraction,
                 'extra_yukawa_fraction': args.extra_yukawa_fraction,
                 'no_local_refinement': args.no_local_refinement,
                 'matching_scheme': 'full-general-Matchete-Yukawa-channel-v1',
                 'reuse_b_run': args.reuse_b_run, 'reuse_c_run': args.reuse_c_run}
    if args.preflight:
        rows = []
        for i, model in enumerate(selected,1):
            try:
                rec = preflight(model, args)
            except Exception as exc:
                rec = {'model': model, 'status': 'PreflightFailed', 'message': str(exc)}
            rows.append(rec)
            print(f'[{i}/{len(selected)}] {model}: {rec["status"]}', flush=True)
            if rec['status'] != 'Ready':
                print('   ',rec.get('checks', rec.get('message')),flush=True)
        save(out / 'preflight.json', {'models': rows})
        print('Preflight:',out / 'preflight.json')
        return 0 if all(r['status'] == 'Ready' for r in rows) else 2

    signature_path = out / 'settings.json'
    if signature_path.is_file():
        if read(signature_path) != signature:
            p.error('Resume settings differ from existing campaign; choose another --run-id')
    else:
        save(signature_path, signature)

    rows = []
    for index, model in enumerate(selected, 1):
        dest = out / tag(model)
        print(f'[{index}/{len(selected)}] Fitting {model}', flush=True)
        try:
            prior = dest / 'result.json'
            if args.resume and prior.is_file():
                data = read(prior)
                if 'chi2_final' not in data or data.get('model') != model:
                    raise ValueError('Previously saved model result invalid')
                result = data
                status = 'ReusedFit'
            elif model in MIXED:
                from scripts import RunT3BCGeneralFromReference as bc
                external = args.reuse_b_run if model == 'B_-1' else args.reuse_c_run
                if external:
                    found = bc.OUT / external / tag(model) / 'result.json'
                    result = read(found)
                    if result.get('model') != model:
                        raise ValueError('Wrong reused mixed-scalar model')
                    if (not math.isclose(float(result['input_lambda_t3']),
                                         float(read(ref_path(model))['base_state']['ordinary']['lambdaT3']['real'])) or
                        not math.isclose(float(result['extra_quartic_fraction']), args.quartic_fraction) or
                        not math.isclose(float(result['extra_mixed_scalar_bilinear_fraction']), args.mu2_fraction)):
                        raise ValueError('Reused B/C benchmark has different fixed scalar inputs')
                    external_dir = found.parent
                else:
                    # Call the exact already-working direct B/C fitting implementation.
                    child = argparse.Namespace(
                        model=model, run_id=args.run_id, lambda_t3=None,
                        quartic_fraction=args.quartic_fraction, mu2_fraction=args.mu2_fraction,
                        export_betas=args.export_bc_betas, wolframscript=args.wolframscript,
                        sobol_points=args.sobol_points, seed=args.seed + MODELS.index(model),
                        bound=args.bound, sobol_fraction=args.sobol_fraction,
                        sobol_floor=args.sobol_floor, seed_count=args.seed_count,
                        target_chi2=args.target_chi2, max_nfev=args.max_nfev,
                        no_local_refinement=args.no_local_refinement)
                    external_dir = bc.run(child)
                    result = read(external_dir / 'result.json')
                dest.mkdir(parents=True, exist_ok=True)
                save(dest / 'result.json', dict(result, source_result_folder=str(external_dir),
                                                orchestration=FORMAT))
                save(dest / 'fitted_config.json', read(external_dir / 'fitted_config.json'))
                src_csv = external_dir / 'observables.csv'
                if src_csv.is_file():
                    (dest / 'observables.csv').write_bytes(src_csv.read_bytes())
                status = 'ExploratoryAcceptedFit' if float(result['chi2_final']) <= 5.86 else 'ExploratoryPoorFit'
            else:
                result = fit_nonmixed(model, args, dest)
                status = result['status']
            row = {'model': model, 'status': status,
                   'chi2_initial': result.get('chi2_seed'),
                   'chi2_final': result['chi2_final'],
                   'prediction': result.get('prediction'),
                   'method': ('mixed-scalar-conditional-one-scale' if model in MIXED
                              else 'canonical-full-general-rges-plus-general-Matchete-C5'),
                   'output_folder': str(dest)}
        except Exception as exc:
            dest.mkdir(parents=True, exist_ok=True)
            error = {'model': model, 'status': 'Blocked',
                     'error': f'{type(exc).__name__}: {exc}',
                     'traceback': traceback.format_exc(limit=9)}
            save(dest / 'error.json', error)
            row = {'model': model, 'status': 'Blocked',
                   'message': error['error'], 'output_folder': str(dest)}
            print(f'   BLOCKED: {error["error"]}', flush=True)
        if row['status'] != 'Blocked':
            (dest / 'error.json').unlink(missing_ok=True)
        rows.append(row)
        summarise(out, rows)
        if 'chi2_final' in row:
            print(f'   chi2={row["chi2_final"]:.8g}', flush=True)
    print(f'Campaign: {out / "summary.csv"}',flush=True)
    passed = sum('chi2_final' in row for row in rows)
    print(f'Completed {passed}/{len(selected)} model results.',flush=True)
    return 0 if passed == len(selected) else 2


if __name__ == '__main__':
    raise SystemExit(main())
