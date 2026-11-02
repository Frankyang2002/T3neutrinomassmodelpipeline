"""A(-2) general-Lagrangian Yukawa benchmark and y2Conj sensitivity.

This is an additive diagnostic runner: it does NOT replace the fitter.
All cases use the existing full-general UV/EFT1 running and Matchete
self-bilinear mass-eigenstate C5 evaluator with the same nonzero bilinear.
Only the independent UV boundary y2Conj is varied, not the C5 formula.

PowerShell (project root):
  python -X utf8 -m scripts.RunT3Aminus2YukawaBenchmark --benchmark-id a2_reference_v1
  python -X utf8 -m scripts.RunT3Aminus2YukawaBenchmark --benchmark-id a2_fit_from16_v1 \
      --fit-run-id thesis16_general_c5_full_v1
"""
from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
import math
from pathlib import Path
import re
import sys
from typing import Any

import numpy as np

from scripts import RunT3GeneralAll16FromReference as campaign

ROOT = Path(__file__).resolve().parents[1]
MODEL = 'A_-2'
OUT = ROOT / 'output' / 'studies' / 'general_scalar_extended' / 'a_m2_yukawa_benchmarks'
CAMPAIGN_OUT = campaign.OUTPUT
DEFAULT_FRACTIONS = (0.0, 0.05, 0.10)


def require(test: bool, explanation: str) -> None:
    if not test:
        raise ValueError('A(-2) Yukawa benchmark: ' + explanation)


def parse_fractions(raw: str) -> tuple[float, ...]:
    try:
        values = tuple(float(v.strip()) for v in raw.split(','))
    except ValueError as exc:
        raise ValueError('Fractions must be comma-separated numbers') from exc
    require(bool(values) and len(values) <= 12 and
            all(math.isfinite(v) and 0 <= v <= 1 for v in values),
            'fractions must contain 1–12 values between 0 and 1')
    require(len(set(values)) == len(values), 'repeated fractions')
    return values


def complex_value(value: Any) -> complex:
    if isinstance(value, dict):
        require('real' in value and 'imag' in value, 'malformed complex JSON')
        return complex(float(value['real']), float(value['imag']))
    return complex(value)


def complex_matrix(raw: Any) -> np.ndarray:
    """Read either reference {real,imag} matrices or per-entry matrix JSON."""
    if isinstance(raw, dict) and 'real' in raw and 'imag' in raw:
        arr = np.asarray(raw['real'], dtype=float) + 1j * np.asarray(raw['imag'], dtype=float)
    else:
        arr = np.asarray([[complex_value(x) for x in row] for row in raw], dtype=complex)
    require(arr.shape == (3, 3) and np.isfinite(arr.real).all() and
            np.isfinite(arr.imag).all(), 'invalid 3x3 Yukawa')
    return arr


def encoded_matrix(a: np.ndarray) -> dict[str, Any]:
    return campaign.matrix_json(a)


def yaml_free_json(value: Any) -> Any:
    """Ensure JSON-friendly primitive values for NumPy derived diagnostics."""
    if isinstance(value, np.ndarray):
        return yaml_free_json(value.tolist())
    if isinstance(value, np.generic):
        return yaml_free_json(value.item())
    if isinstance(value, complex):
        return {'real': value.real, 'imag': value.imag}
    if isinstance(value, dict):
        return {str(k): yaml_free_json(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [yaml_free_json(x) for x in value]
    return value


def spectral_masses(endpoint: dict[str, Any]) -> dict[str, Any]:
    """Use the existing positive-real bilinear mass-splitting convention."""
    require('muT3GeneralS2Sq' in endpoint, 'EFT1 missing running S2 self-bilinear')
    mu = complex_value(endpoint['muT3GeneralS2Sq'])
    m2sq = float(endpoint['mS2Sq'])
    m1sq = float(endpoint['mS1Sq'])
    require(np.isfinite([mu.real, mu.imag, m1sq, m2sq]).all(),
            'nonfinite scalar masses')
    require(mu.real > 0 and abs(mu.imag) <= 1e-10 * max(1.0, abs(mu)),
            'current matched mass basis requires positive-real bilinear')
    minus2, plus2 = m2sq - 2 * abs(mu), m2sq + 2 * abs(mu)
    require(minus2 > 0 and plus2 > 0 and m1sq > 0, 'tachyonic scalar masses')
    return {'mu2_gev2': {'real': mu.real, 'imag': mu.imag},
            'mS1_gev': math.sqrt(m1sq),
            'mS2_unsplit_gev': math.sqrt(m2sq),
            'mS2_minus_gev': math.sqrt(minus2),
            'mS2_plus_gev': math.sqrt(plus2),
            'mS2_minus_sq_gev2': minus2, 'mS2_plus_sq_gev2': plus2}


def params_from_config(config: dict[str, Any]) -> dict[str, float]:
    ordinary = config['base_state']['ordinary']
    arrays = (complex_matrix(ordinary['y1']), complex_matrix(ordinary['y2']))
    require(all(np.max(np.abs(a.imag)) < 1e-12 for a in arrays),
            'current benchmark fitter supports real ordinary Yukawas')
    return {f'{name}_{i}{j}': float(a[i, j].real)
            for name, a in zip(('y1', 'y2'), arrays)
            for i in range(3) for j in range(3)}


def load_fit(run_id: str, reference: dict[str, Any]) -> tuple[dict, dict, dict]:
    """Read A(-2) from a *single- or multi-model* full-general campaign.

    The fitted configuration and the numerical optimizer result must agree.
    A partial/in-progress campaign cannot be mistaken for a completed A(-2) fit.
    """
    require(bool(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,63}', run_id)),
            'unsafe fit-run-id')
    run_root = CAMPAIGN_OUT / run_id
    result_path = run_root / 'A_m2' / 'result.json'
    config_path = run_root / 'A_m2' / 'fitted_config.json'
    settings_path = run_root / 'settings.json'
    for needed in (settings_path, result_path, config_path):
        require(needed.is_file(), 'missing '+str(needed)+
                '; A(-2) must have completed in the 16-model campaign')
    result = campaign.read(result_path)
    fitted = campaign.read(config_path)
    settings = campaign.read(settings_path)

    require(result.get('model') == MODEL, 'incorrect source fit model')
    require(result.get('status') in ('ExploratoryAcceptedFit', 'ExploratoryPoorFit'),
            'A(-2) result is not a finished general fit')
    require('chi2_final' in result and isinstance(result.get('best_parameters'), dict),
            'incomplete original fit result')
    require(math.isfinite(float(result['chi2_final'])), 'nonfinite original fit chi2')
    require(settings.get('matching_scheme') == 'full-general-Matchete-Yukawa-channel-massbasis-v2',
            'source run is not the general-matching mass-basis-v2 campaign: '
            + repr(settings.get('matching_scheme')))
    models = settings.get('models')
    require(isinstance(models, list) and MODEL in models,
            'source campaign does not include A(-2)')
    require(fitted.get('representation') == reference.get('representation'),
            'reference/fitted representation mismatch')
    source_sha = result.get('input_sha256', {}).get('c5_kernel')
    source = campaign.general_coefficient_path(MODEL)
    require(source.is_file() and source_sha == sha256(source.read_bytes()).hexdigest(),
            'matching coefficient differs from original fitted run')

    params = {key: float(result['best_parameters'][key]) for key in campaign.NAMES}
    require(all(math.isfinite(x) for x in params.values()), 'nonfinite fitted Yukawas')
    config_params = params_from_config(fitted)
    require(set(config_params) == set(params), 'fitted configuration Yukawa keys differ')
    for key, val in params.items():
        require(math.isclose(val, config_params[key], rel_tol=1e-10, abs_tol=1e-14),
                f'fitted_config Yukawa {key} differs from optimizer best_parameters')

    fixed = result.get('fixed_general_scalar_boundary', {})
    require(math.isclose(float(settings['mu2_fraction']),
                         float(fixed['bilinear_fraction_min_mass_squared'])),
            'fit bilinear metadata mismatch')
    require(math.isclose(float(settings['quartic_fraction']),
                         float(fixed['quartic_fraction_of_lambdaT3'])),
            'fit quartic metadata mismatch')
    return params, settings, {
        'source': 'completed_A_m2_result_from_general_campaign',
        'fit_run_id': run_id,
        'source_campaign_model_count': len(models),
        'fit_chi2_original': float(result['chi2_final']),
        'fit_best_source': result['best_source'],
        'fit_result_file': str(result_path),
        'fit_config_file': str(config_path),
        'fit_result_sha256': sha256(result_path.read_bytes()).hexdigest(),
        'fit_config_sha256': sha256(config_path.read_bytes()).hexdigest(),
    }


def write_csv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    with path.open('w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def run(args: argparse.Namespace) -> int:
    from Numerical.fitting.OscillationFit import OscillationFitTarget
    from Numerical.orchestration.RunSmallLambdaRefit import TARGET

    require(bool(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,63}', args.benchmark_id)),
            'unsafe benchmark-id')
    fractions = parse_fractions(args.y2conj_fractions)
    output = OUT / args.benchmark_id
    require(not output.exists(), 'output exists; choose a new --benchmark-id')

    reference = campaign.read(campaign.ref_path(MODEL))
    params = params_from_config(reference)
    config = {'quartic_fraction': args.quartic_fraction,
              'mu2_fraction': args.mu2_fraction,
              'canonical_run_id': args.canonical_run_id}
    metadata = {'source': 'repository_reference'}
    if args.fit_run_id:
        params, settings, metadata = load_fit(args.fit_run_id, reference)
        config['quartic_fraction'] = float(settings['quartic_fraction'])
        config['mu2_fraction'] = float(settings['mu2_fraction'])
        config['canonical_run_id'] = settings['canonical_run_id']

    uv, eft, provenance, kernel, coefficient_file, kind = campaign.nonmixed_inputs(
        MODEL, config['canonical_run_id'], generate_general_c5=False,
        wolframscript=args.wolframscript)
    require({'y1', 'y2', 'y2Conj'} <= kernel.required_yukawas,
            'general Matchete C5 does not retain the expected Yukawas')
    require({'T32PlusMass', 'T32MinusMass'} <= kernel.symbols,
            'general C5 does not depend on both split scalar masses')
    status = campaign.read(campaign.general_coefficient_status(MODEL))
    require(status.get('matching_scheme') ==
            'FullGeneralSelfBilinearRealMassEigenfieldsOneScaleCandidate',
            'wrong matching scheme: expected full-general self-bilinear mass basis')

    output.mkdir(parents=True)
    target = OscillationFitTarget.from_json(TARGET)
    variants = []
    matrices = []
    summary = []
    y1 = np.array([[params[f'y1_{i}{j}'] for j in range(3)] for i in range(3)])
    y2 = np.array([[params[f'y2_{i}{j}'] for j in range(3)] for i in range(3)])
    print('A(-2) benchmark: '+('fitted run '+args.fit_run_id if args.fit_run_id else 'reference Yukawas'),
          flush=True)
    print('Y1 at UV boundary:\n', np.array2string(y1, precision=8), flush=True)
    print('Y2 at UV boundary:\n', np.array2string(y2, precision=8), flush=True)
    reference_y2 = complex_matrix(reference['base_state']['ordinary']['y2'])
    require(np.max(np.abs(reference_y2.imag)) < 1e-12, 'complex reference y2 unsupported')

    for fraction in fractions:
        label = 'y2conj_' + f'{fraction:.6g}'.replace('.', 'p')
        fixed = campaign.fixed_overrides(
            reference, uv, quartic_fraction=config['quartic_fraction'],
            mu2_fraction=config['mu2_fraction'], extra_yukawa_fraction=fraction)
        y2conj = complex_matrix(fixed['y2Conj'])
        require(np.allclose(y2conj, fraction * reference_y2, rtol=1e-12, atol=1e-14),
                'unexpected y2Conj boundary: fit uses frozen reference y2')
        case = {'case': label, 'y2conj_fraction_of_REFERENCE_y2': fraction,
                'uv_boundary_yukawas': {'y1': encoded_matrix(y1),
                                         'y2': encoded_matrix(y2),
                                         'y2Conj': encoded_matrix(y2conj)},
                'status': 'Blocked'}
        for yuk_name, matrix_value in (('y1', y1), ('y2', y2), ('y2Conj', y2conj)):
            for i in range(3):
                for j in range(3):
                    matrices.append({'case': label, 'scale': 'UV_boundary_1e7GeV',
                                     'yukawa': yuk_name, 'row': i, 'column': j,
                                     'real': float(matrix_value[i, j].real),
                                     'imag': float(matrix_value[i, j].imag)})
        try:
            evaluator = campaign.NonMixedEvaluator(MODEL, uv, eft, kernel, fixed, target)
            fit, detail = evaluator.evaluate(params, detail=True)
            spectral = spectral_masses(detail['eft1_at_scalar_threshold'])
            case.update(status='Evaluated', chi2=float(fit.chi2),
                        prediction=dict(zip(fit.observable_names,
                                            [float(v) for v in fit.prediction])),
                        target=dict(zip(fit.observable_names,
                                        [float(v) for v in fit.target])),
                        whitened_residual=np.asarray(fit.whitened_residual).tolist(),
                        scalar_masses=spectral,
                        c5_at_scalar_threshold=detail['c5_at_scalar_threshold'],
                        c5_at_low_scale=detail['c5_at_low_scale'],
                        running_uv_endpoint=detail['uv_at_fermion_threshold'],
                        running_eft1_endpoint=detail['eft1_at_scalar_threshold'],
                        matching_diagnostics=detail['general_c5_diagnostics'],
                        uv_nfev=detail['uv_nfev'], eft1_nfev=detail['eft1_nfev'])
            for yuk_name in ('y1', 'y2', 'y2Conj'):
                endpoint = complex_matrix(detail['uv_at_fermion_threshold'][yuk_name])
                for i in range(3):
                    for j in range(3):
                        matrices.append({'case': label, 'scale': 'fermion_threshold_1e5GeV',
                                         'yukawa': yuk_name, 'row': i, 'column': j,
                                         'real': float(endpoint[i, j].real),
                                         'imag': float(endpoint[i, j].imag)})
            row = {'case': label, 'fraction': fraction, 'status': 'Evaluated',
                   'chi2': float(fit.chi2),
                   'mS2_minus_gev': spectral['mS2_minus_gev'],
                   'mS2_plus_gev': spectral['mS2_plus_gev'],
                   'mu2_gev2': spectral['mu2_gev2']['real'],
                   **case['prediction']}
            summary.append(row)
            print(f'  {label}: chi2={fit.chi2:.8g} '
                  f'mS2minus={spectral["mS2_minus_gev"]:.6g} GeV '
                  f'mS2plus={spectral["mS2_plus_gev"]:.6g} GeV', flush=True)
        except Exception as exc:
            case['error'] = f'{type(exc).__name__}: {exc}'
            summary.append({'case': label, 'fraction': fraction,
                            'status': 'Blocked', 'error': case['error']})
            print('  '+label+': BLOCKED '+case['error'], flush=True)
        variants.append(case)
        campaign.save(output / (label + '.json'), yaml_free_json(case))

    write_csv(output / 'summary.csv', ['case', 'fraction', 'status', 'chi2',
              'delta_m21_sq_ev2', 'delta_m3l_sq_ev2', 'sin2_theta12',
              'sin2_theta13', 'sin2_theta23', 'mS2_minus_gev',
              'mS2_plus_gev', 'mu2_gev2', 'error'], summary)
    write_csv(output / 'yukawa_matrices.csv', ['case','scale','yukawa',
              'row','column','real','imag'], matrices)
    manifest = {
        'format': 't3-general-A-m2-yukawa-sensitivity-v1', 'model': MODEL,
        'scope': 'General UV/EFT1 RGEs + conditional mass-basis one-scale Matchete C5',
        'physical_sequential_matching_certified': False,
        'source': metadata,
        'benchmark_config': config,
        'y2Conj_held_fixed_per_case_during_original_y1_y2_fit': True,
        'y2Conj_fraction_multiplies_reference_y2_not_fitted_y2': True,
        'bilinear_nonzero_in_all_cases': True,
        'matchete_c5_source': str(coefficient_file),
        'matchete_c5_sha256': sha256(coefficient_file.read_bytes()).hexdigest(),
        'source_provenance': provenance,
        'matching_kind': kind,
        'matched_channels': sorted(kernel.required_yukawas),
        'cases': [{'name': row['case'], 'status': row['status'],
                   'chi2': row.get('chi2'), 'error': row.get('error')} for row in summary]
    }
    campaign.save(output / 'manifest.json', yaml_free_json(manifest))
    print('Output:', output)
    return 0 if all(x['status'] == 'Evaluated' for x in summary) else 2


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--benchmark-id', required=True)
    parser.add_argument('--fit-run-id',
                        help='Use completed A(-2) fit from a single- or 16-model general campaign; otherwise reference seed')
    parser.add_argument('--canonical-run-id', default='full16_general_v1')
    parser.add_argument('--quartic-fraction', type=float, default=0.1)
    parser.add_argument('--mu2-fraction', type=float, default=0.01)
    parser.add_argument('--y2conj-fractions', default='0,0.05,0.1',
                        help='UV fractions multiplying fixed REFERENCE y2')
    parser.add_argument('--wolframscript', default='wolframscript')
    args = parser.parse_args(argv)
    try:
        require(0 <= args.quartic_fraction <= 1 and 0 < args.mu2_fraction <= 0.2,
                'invalid quartic/bilinear fractions')
        return run(args)
    except (ValueError, FileNotFoundError, KeyError) as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
