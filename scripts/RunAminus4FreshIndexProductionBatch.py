"""A(-4) EFT1: validate SOURCE-LEVEL fresh-index fix, regenerate 17 native betas,
check 55 independent scalar-loop polynomials, and integrate seven-scale
renormalisable EFT1 ODE. Read-only on source/archives after ZIP extraction.

The corrected COMPLETE RGE/running/rgbeta/T3RGBetaModel.wl must already be
installed. Never modifies installed RGBeta/Matchete or historical results.

Run from repository root:
  python scripts/RunAminus4FreshIndexProductionBatch.py \
    --output-root output/a_minus4_next_batch/fresh_index_production

This does NOT certify the native Matchete quartic coupling dictionary,
GammaCC/Weyl phase or the physical Wilson coefficient C(muF).
"""
from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from Numerical.core.RGBetaEvaluator import evaluate_inputform_expression
from Numerical.running import Aminus4EFT1Couplings as eft
from Numerical.running.Aminus4QuarticHessianIdentities import beta_monomials_exact
from scripts.ExportAminus4EFT1Couplings import make_state
from scripts.RunAminus4FreshIndexCausalBatch import inspect_polynomial

MODEL = Path('RGE/running/rgbeta/T3RGBetaModel.wl')
RUNNER = Path('RGE/running/rgbeta/RunGeneralT3RGBeta.wl')
OLD_RAW = Path('output/a_minus4_eft1_batch2/A_m4_projected_raw.json')
OLD_BASIS = Path('output/a_minus4_eft1_batch2/A_m4_eft1_tensor_basis.json')
CAUSAL_BOTH = Path('output/a_minus4_next_batch/fresh_index_causal_v2/native/both.json')
PRIOR_ODE = Path('output/a_minus4_next_batch/conditional_quartic_ode')
ORIGINAL_SOURCE_SHA = 'a58c79f18b2abdf9d476588552dfa52e2f131f322094e5b72d1198162aa20abc'
FIXED_SOURCE_SHA = '147b25f172ebcc3644e3ae1be8a7686fe66e1b7c1673ca17dd10e444ce0fff0d'
QUARTICS = frozenset({'lambdaH', 'lambdaS1', 'lambdaS2', 'lambdaH1',
                       'lambdaH2', 'lambda12', 'lambdaT3',
                       'lambdaH2Adj', 'lambdaS2Adj'})
NONQUARTICS = eft.BETA_NAMES - QUARTICS
SCALAR_COORDS = ('lambdaH', 'lambdaS1', 'lambdaS2', 'lambdaS2Adj',
                 'lambdaH1', 'lambdaH2', 'lambdaH2Adj', 'lambda12',
                 'lambdaT3.real', 'lambdaT3.imag')


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def read(path: Path) -> dict:
    require(path.is_file(), 'Missing required input: ' + str(path))
    obj = json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(obj, dict), 'Expected JSON object: ' + str(path))
    return obj


def digest(path: Path) -> str:
    require(path.is_file(), 'Missing source: ' + str(path))
    return sha256(path.read_bytes()).hexdigest()


def check_model_source(root: Path) -> dict:
    p = root / MODEL
    actual = digest(p)
    require(actual == FIXED_SOURCE_SHA,
            f'Wrong A(-4) corrected source SHA-256: {actual}; expected {FIXED_SOURCE_SHA}. '
            'Do not proceed if local source differs from reviewed replacement.')
    src = p.read_text(encoding='utf8')
    require('T3RGBetaPortalAdjInvariant[3] :=' in src and
            'lambdaS2Adj,\n            {Bar @ S2, S2, Bar @ S2, S2}' in src and
            src.count('Module[{adj},') == 2,
            'Expected exactly two newly localised A(-4) adjoint invariants')
    return {'path': str(MODEL), 'sha256': actual,
            'original_sha256': ORIGINAL_SOURCE_SHA,
            'source_fix_numerical_scope': 'lambdaS2Adj and lambdaH2Adj in A(-4)',
            'other_models_not_numerically_certified': True}


def compatible_native(raw: dict) -> dict:
    require(raw.get('status') == 'Success' and
            raw.get('metadata', {}).get('GeneralScalarModel') == 'A_-4' and
            raw['metadata'].get('GeneralScalarStage') == 'EFT1' and
            raw['metadata'].get('GeneralScalarQuarticBetaBasis') == 'RGBetaFullQuarticProjection',
            'Regenerated betas are not the full native A(-4) EFT1 projection')
    require(raw['metadata'].get('GeneralScalarAdditionalQuartics') == 0 and
            raw['metadata'].get('GeneralScalarExistingStructuralCount') == 0 and
            set(raw['metadata'].get('GeneralScalarExistingQuarticBetasProjected', [])) == QUARTICS,
            'Unexpected A(-4) operator inventory')
    require(set(raw.get('report_betas', {})) == eft.BETA_NAMES,
            'Regenerated EFT1 must have exactly 17 running parameters')
    require(all(isinstance(v, str) and v.strip() and '$Failed' not in v and
                'BetaTerm[' not in v for v in raw['report_betas'].values()),
            'One or more native beta expressions is unresolved')
    payload = dict(raw)
    payload['metadata'] = dict(raw['metadata'])
    payload['metadata']['IntegratedField'] = 'F'
    payload['metadata']['ActiveBSMFields'] = ['S1', 'S2']
    eft.validate_beta_payload(payload)
    return payload


def value_maxabs(x: Any) -> float:
    return float(np.max(np.abs(np.asarray(x))))


def numerically_compare_beta_sets(expected: dict, actual: dict, points: list[dict],
                                  names: set[str] | frozenset[str], *, tolerance: float=2e-9) -> dict:
    """Compare complete beta expressions at nonzero gauge/Yukawa/scalar points."""
    require(set(names) <= set(expected) and set(names) <= set(actual),
            'Missing beta channel for full-expression comparison')
    max_normalized = 0.0
    worst = None
    for i, point in enumerate(points):
        env = eft._point(point)
        for name in sorted(names):
            a = evaluate_inputform_expression(expected[name], env)
            b = evaluate_inputform_expression(actual[name], env)
            diff = value_maxabs(np.asarray(a) - np.asarray(b))
            scaled = diff / max(1.0, value_maxabs(a), value_maxabs(b))
            if scaled > max_normalized:
                max_normalized, worst = scaled, {'point': i, 'name': name,
                                                'absolute': diff, 'normalized': scaled}
    require(max_normalized < tolerance,
            'Full-expression beta comparison failed: ' + str(worst))
    return {'samples': len(points), 'channels': sorted(names),
            'maximum_normalized_difference': max_normalized, 'worst': worst,
            'passes': True}


def native_export(root: Path, destination: Path, log: Path) -> None:
    exe = shutil.which('wolframscript')
    require(exe is not None, 'wolframscript is not on PATH')
    command = [exe, '-file', str(root / RUNNER), 'EFT1',
               '1', '3', '2', 'm4', 'A_-4', str(destination)]
    result = subprocess.run(command, cwd=root, capture_output=True,
                            text=True, check=False)
    log.write_text('Command: ' + subprocess.list2cmdline(command) + '\n' +
                   f'Exit code: {result.returncode}\nSTDOUT:\n{result.stdout}\n'
                   f'STDERR:\n{result.stderr}\n', encoding='utf8')
    require(result.returncode == 0 and destination.is_file(),
            f'Native beta regeneration failed (exit {result.returncode}); see {log}')


def coordinates(point: dict) -> np.ndarray:
    p = eft._point(point)
    return np.array([p[k] if not k.startswith('lambdaT3.') else
                     (complex(p['lambdaT3']).real if k.endswith('.real')
                      else complex(p['lambdaT3']).imag)
                     for k in SCALAR_COORDS], dtype=float)


def compare_saved(root: Path, new_info: dict, old_basis: Path) -> dict:
    old = read(root / old_basis)
    require(old.get('model') == 'A_-4' and len(old.get('running_couplings', [])) == 7,
            'Previous seven-scale reference data not found')
    new_grid = np.asarray(new_info['sampled_scales_gev'], float)
    old_grid = np.asarray(old['sampled_scales_gev'], float)
    require(np.allclose(new_grid, old_grid, rtol=1e-12, atol=1e-8),
            'Historical reference and new ODE scale grids differ')
    rows = []
    for newp, oldp, mu in zip(new_info['running_couplings'], old['running_couplings'], new_grid):
        diff = coordinates(newp) - coordinates(oldp)
        rows.append({'mu_gev': float(mu), **{
            'fixed_minus_previous_' + key: float(val)
            for key, val in zip(SCALAR_COORDS, diff)}})
    return {'endpoint': rows[-1], 'all_scales': rows,
            'max_scalar_quartic_absolute_grid_difference': max(
                abs(v) for row in rows for k, v in row.items() if k != 'mu_gev')}


def evaluate(root: Path, output_root: Path, *, samples: int=7) -> dict:
    root = root.resolve()
    output_root = (root / output_root).resolve()
    require(not output_root.exists(), 'Refusing to overwrite output: ' + str(output_root))
    require(samples == 7, 'Historical A(-4) comparison requires seven samples')
    source = check_model_source(root)
    require((root / CAUSAL_BOTH).is_file(),
            'Previous successful four-mode native experiment required: ' + str(CAUSAL_BOTH))
    output_root.mkdir(parents=True)
    native_path = output_root / 'A_m4_fixed_native_rgbeta.json'
    log = output_root / 'wolfram_export.log'
    native_export(root, native_path, log)
    fresh = compatible_native(read(native_path))
    old = compatible_native(read(root / OLD_RAW))
    witness = read(root / CAUSAL_BOTH)
    require(witness.get('mode') == 'both' and
            witness.get('fresh_triplet_adjoint_dummy') is True and
            witness.get('fresh_portal_adjoint_dummy') is True and
            witness.get('source_model_SHA256') == ORIGINAL_SOURCE_SHA and
            set(witness.get('report_betas', {})) == QUARTICS,
            'Previous native both-index comparison has changed')

    basis = read(root / OLD_BASIS)
    points = basis.get('running_couplings')
    require(isinstance(points, list) and len(points) == 7,
            'Missing archived A(-4) 17-parameter boundary and trajectory')
    selected = [points[i] for i in (0, 3, 6)]
    # Ensure nothing outside the quartic beta sector changes at finite background.
    nonquartic_check = numerically_compare_beta_sets(old['report_betas'],
                        fresh['report_betas'], selected, NONQUARTICS)
    # More demanding than 55 gaugeless polynomial probes: all NINE complete
    # (g,Y,lambda)-dependent quartic betas must equal the previously verified
    # in-memory both-fresh-index native experiment at nonzero SM couplings.
    full_quartic_check = numerically_compare_beta_sets(
                        witness['report_betas'], fresh['report_betas'], selected,
                        QUARTICS)
    hessian = beta_monomials_exact()
    scalar = inspect_polynomial(fresh, hessian)
    require(scalar['native_vs_exact_hessian_failing_probes'] == 0 and
            scalar['discrepant_monomial_count'] == 0 and
            scalar['max_reference_hessian_minus_fourth_derivative_abs'] < 2e-8,
            'Corrected native scalar beta fails independent exact 55-probe reference')

    compatible_path = output_root / 'A_m4_fixed_17beta_compatible.json'
    compatible_path.write_text(json.dumps(fresh, indent=2) + '\n', encoding='utf8')
    start = {'model': 'A_-4', 'MSbar_mu_F_gev': float(basis['sampled_scales_gev'][0]),
             'muF_renormalizable_state': points[0]}
    initial = make_state(start)
    mu_s = float(basis['sampled_scales_gev'][-1])
    running = eft.integrate_eft1(root, initial, compatible_path, mu_s,
                                  n_samples=samples)
    out_basis = eft.write_sampled_npz(running, output_root / 'running',
                            source_root=root, source_beta=compatible_path)
    require(out_basis['physical_full_Wilson_running_ready'] is False and
            out_basis['matchete_scalar_and_Wilson_basis_map_certified'] is False,
            'Unsafe physical certification was added to tensor output')
    compare_diag = compare_saved(root, out_basis,
        PRIOR_ODE / 'reference_diagnostic' / 'A_m4_eft1_tensor_basis.json')
    compare_original = compare_saved(root, out_basis,
        PRIOR_ODE / 'native' / 'A_m4_eft1_tensor_basis.json')
    summary = {'format': 't3-a-minus4-permanent-source-fresh-indices-v1',
       'model': 'A_-4',
       'status': 'SOURCE_LEVEL_FRESH_INDEX_REGEN_55_PROBES_PASS_ODE_COMPUTED_NOT_PHYSICAL',
       'source': source,
       'native_rgbeta_17_parameter_inventory': sorted(fresh['report_betas']),
       'native_17_parameter_beta_integrity': True,
       'nonquartic_betas_unchanged_at_nonzero_background': nonquartic_check,
       'all_nine_full_quartic_betas_match_previous_both_fresh_native': full_quartic_check,
       'independent_scalar_only_55_probe_check': {k:v for k,v in scalar.items()
                                                   if k!='all_55_probe_details'},
       'new_beta_sha256': digest(compatible_path),
       'ode_renormalisable_couplings_7_scales_computed': True,
       'fixed_minus_previous_reference_diagnostic': {k:v for k,v in compare_diag.items()
                                                       if k!='all_scales'},
       'fixed_minus_original_native': {k:v for k,v in compare_original.items()
                                       if k!='all_scales'},
       'physical_matchete_to_rgbeta_quartic_map_certified': False,
       'matchete_GammaCC_to_weyl_factor_measured': False,
       'physical_C_ijab_at_muF': None,
       'physical_full_Wilson_ODE_enabled': False,
       'complete_one_loop_F_threshold_matching_certified': False,
       'scalar_threshold_matching_performed': False,
       'EFT_order': 'UV -> integrate out ALL fermions -> EFT1(SM+S1+S2) -> integrate out ALL scalars simultaneously -> SMEFT'}
    with (output_root / 'coupling_comparison.csv').open('w', newline='', encoding='utf8') as f:
        writer = csv.DictWriter(f, fieldnames=list(compare_diag['all_scales'][0]))
        writer.writeheader()
        writer.writerows(compare_diag['all_scales'])
    (output_root / 'summary.json').write_text(
        json.dumps(summary, indent=2, allow_nan=False) + '\n', encoding='utf8')
    return summary


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, default=ROOT)
    p.add_argument('--output-root', type=Path,
                   default=Path('output/a_minus4_next_batch/fresh_index_production'))
    args = p.parse_args(argv)
    try:
        report = evaluate(args.root, args.output_root)
    except (RuntimeError, ValueError, FileNotFoundError) as exc:
        print('A(-4) permanent-source batch BLOCKED:', exc, file=sys.stderr)
        return 1
    print('Corrected production RGBeta: all 17 beta channels validated')
    print('Full nine quartic betas vs previous both-fresh native: PASS')
    print('Independent scalar-only Hessian: 0/55 discrepancies')
    print('Seven-scale A(-4) renormalisable EFT1 ODE: PASS')
    print('Endpoint corrected-minus-previous-diagnostic lambdaH2Adj:',
          report['fixed_minus_previous_reference_diagnostic']['endpoint'][
              'fixed_minus_previous_lambdaH2Adj'])
    print('Physical Matchete/Wilson certification: NO')
    print('Report:', (args.root / args.output_root / 'summary.json').resolve())
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
