"""Read-only diagnostic of corrected_llss16_v4 EFT1 LLSS contributions.

Project path: scripts/AnalyseT3EFT1Contributions.py

Reads ONLY existing cases.json files. Never refits, evolves RGEs, changes
benchmarks, or substitutes archived fermion sources for fitted boundaries.

The saved cases include Weinberg LLSS logarithms and finite scalar bubbles,
resolved by quartic and fermion-tree source. They do NOT contain the raw,
fully indexed C_LLSS(mu_F) tensors, or independently matched hard C5(mu_F).
Those quantities are reported as unavailable, not invented or inferred from C5(mu_S).

PowerShell (repository root):
 python -X utf8 -m scripts.AnalyseT3EFT1Contributions --run-id corrected_llss16_v4
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CAMPAIGNS = ROOT / 'output/studies/general_scalar_extended/corrected_llss_eft1_freeze'
REPORTS = ROOT / 'Reports/output/thesis_qualitative/eft1_llss_contributions'
MODELS = tuple(f'{family}_{a:+d}' for family, alphas in {
    'A': (-4, -2, 0), 'B': (-3, -1, 1), 'C': (-3, -1, 1),
    'D': (-2, 0, 2), 'E': (-4, -2, 0, 2),
}.items() for a in alphas)
OBSERVABLES = ('delta_m21_sq_ev2', 'delta_m3l_sq_ev2',
               'sin2_theta12', 'sin2_theta13', 'sin2_theta23')


def tag(model: str) -> str:
    family, alpha = model.split('_')
    n = int(alpha)
    return f'{family}_{"m" if n < 0 else "p"}{abs(n)}'


def load(path: Path) -> dict:
    with path.open(encoding='utf-8-sig') as handle:
        result = json.load(handle)
    if not isinstance(result, dict):
        raise ValueError(f'Expected JSON object: {path}')
    return result


def matrix(data: dict, what: str) -> np.ndarray:
    if not isinstance(data, dict) or not {'real', 'imag'} <= set(data):
        raise ValueError(f'Missing complex matrix {what}')
    arr = np.asarray(data['real'], float) + 1j * np.asarray(data['imag'], float)
    if arr.shape != (3, 3) or not np.isfinite(arr).all():
        raise ValueError(f'Invalid finite 3x3 matrix {what}: {arr.shape}')
    return arr


def norm(a: np.ndarray) -> float:
    return float(np.linalg.norm(a))


def scalar_projection(a: np.ndarray, total: np.ndarray) -> float:
    den = norm(total)**2
    return float(np.vdot(total, a).real / den) if den else float('nan')


def write_csv(path: Path, rows: list[dict], columns: list[str]) -> None:
    with path.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def rel_residual(actual: np.ndarray, expected: np.ndarray) -> float:
    return norm(actual - expected) / max(norm(expected), 1e-40)


def observe(case: dict, name: str) -> float:
    pred = case['prediction']
    value = float(pred[name])
    if not math.isfinite(value):
        raise ValueError(f'Nonfinite observable {name}')
    return value


def analyse(model: str, path: Path):
    doc = load(path)
    if doc.get('model') != model or doc.get('status') != 'Computed':
        raise ValueError(f'{model}: unavailable or mismatched computed record')
    cases = doc['cases']
    for key in ('1111', '1101', '1011', '1001'):
        if key not in cases:
            raise ValueError(f'{model}: missing {key} case')
    full, log_off, ren_off, both_off = (cases[k] for k in
                                        ('1111', '1101', '1011', '1001'))
    log = matrix(full['c5_induced_wilson_log'], 'LLSS logarithm')
    bubble = matrix(full['c5_joint_scalar_bubble'], 'scalar bubble')
    at_s = matrix(full['c5_at_scalar_threshold'], 'C5(mu_S)')
    if rel_residual(log + bubble, at_s) > 1e-8:
        raise ValueError(f'{model}: C5 at scalar threshold != log + bubble')
    rows = []
    for entry in full['per_general_interaction']:
        category = entry['category']
        if category not in ('quartic_decomposition', 'fermion_tree_decomposition'):
            raise ValueError(f'{model}: unknown decomposition {category}')
        name = str(entry['name'])
        l = matrix(entry['c5_leading_log'], f'{name} leading log')
        b = matrix(entry['c5_joint_scalar_bubble'], f'{name} bubble')
        t = matrix(entry['c5_total'], f'{name} total')
        if rel_residual(l + b, t) > 1e-8:
            raise ValueError(f'{model}: {name} total != log + bubble')
        rows.append({'model': model, 'category': category, 'name': name,
                     'log_norm_GeV_inv': norm(l),
                     'bubble_norm_GeV_inv': norm(b),
                     'combined_norm_GeV_inv': norm(t),
                     'log_projection_fraction': scalar_projection(l, log),
                     'bubble_projection_fraction': scalar_projection(b, bubble),
                     'log_frobenius_fraction': norm(l)/max(norm(log),1e-40),
                     '_log': l, '_bubble': b})
    # Different decomposition categories are ALTERNATIVE sums, not additive to each other.
    for category in ('quartic_decomposition', 'fermion_tree_decomposition'):
        subset = [x for x in rows if x['category'] == category]
        if not subset:
            raise ValueError(f'{model}: missing {category}')
        for key, expected in (('_log', log), ('_bubble', bubble)):
            got = sum((x[key] for x in subset), np.zeros((3,3), complex))
            if rel_residual(got, expected) > 1e-7:
                raise ValueError(f'{model}: {category} {key} sum disagrees with total')
    for row in rows:
        del row['_log']
        del row['_bubble']
    primary = [x for x in rows if x['category'] == 'quartic_decomposition']
    primary.sort(key=lambda x: -abs(x['log_projection_fraction']))
    obs = []
    for name in OBSERVABLES:
        vals = {key: observe(cases[key], name) for key in
                ('1111', '1101', '1011', '1001')}
        obs.append({'model': model, 'observable': name,
                    'full_1111': vals['1111'],
                    'LLSS_log_off_1101': vals['1101'],
                    'renormalisable_EFT1_off_1011': vals['1011'],
                    'both_EFT1_off_1001': vals['1001'],
                    'full_minus_log_off': vals['1111'] - vals['1101'],
                    'full_minus_renormalisable_off': vals['1111'] - vals['1011'],
                    'full_minus_both_off': vals['1111'] - vals['1001']})
    summary = {'model': model, 'status': 'Computed',
               'log_norm_GeV_inv': norm(log),
               'bubble_norm_GeV_inv': norm(bubble),
               'scalar_threshold_c5_norm_GeV_inv': norm(at_s),
               'log_over_bubble_norm': norm(log)/max(norm(bubble),1e-40),
               'log_over_scalar_threshold_c5_norm': norm(log)/max(norm(at_s),1e-40),
               'leading_quartic': primary[0]['name'],
               'leading_quartic_projection_fraction': primary[0]['log_projection_fraction'],
               'quartic_channel_count': len(primary),
               'fermion_source_count': sum(x['category'] == 'fermion_tree_decomposition'
                                           for x in rows),
               'C5_muF_retained_LLSS_sector': 'zero by construction; not hard one-loop matching',
               'C5_muF_full_one_loop': 'not independently computed',
               'C_LLSS_muF_fitted_tensor': 'not stored in corrected cases.json',
               'source_to_Weyl_normalisation_certified': False}
    matrices = {'model': model, 'case': '1111',
                'C5_EFT1_LLSS_leading_log_GeV_inv': full['c5_induced_wilson_log'],
                'C5_scalar_finite_bubble_GeV_inv': full['c5_joint_scalar_bubble'],
                'C5_at_scalar_threshold_GeV_inv': full['c5_at_scalar_threshold'],
                'C5_muF_hard_matched': None, 'C_LLSS_muF_fitted': None,
                'explanation': ('The fitted muF LLSS six-index Wilson tensor is not saved '
                                'in corrected_llss16_v4 cases. It must be reconstructed '
                                'using the fitted UV muF state and original F-tree source. '
                                'No historical reference tensor was substituted.')}
    return summary, rows, obs, matrices


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-id', default='corrected_llss16_v4')
    parser.add_argument('--input', type=Path, default=None,
                        help='Existing corrected campaign directory')
    parser.add_argument('--output', type=Path, default=None,
                        help='New destination, never writes to source campaign')
    parser.add_argument('--models', default=','.join(MODELS))
    args = parser.parse_args(argv)
    if not args.run_id.replace('_','').replace('-','').isalnum():
        parser.error('Invalid run ID')
    requested = [x.strip() for x in args.models.split(',') if x.strip()]
    if not requested or len(requested) != len(set(requested)) or not set(requested) <= set(MODELS):
        parser.error('Invalid model selection')
    source = args.input or CAMPAIGNS / args.run_id
    dest = args.output or REPORTS / args.run_id
    if not source.is_dir():
        raise FileNotFoundError(f'Missing saved run: {source}')
    if source.resolve() == dest.resolve() or source.resolve() in dest.resolve().parents:
        raise ValueError('Output must not be inside source campaign')
    dest.mkdir(parents=True, exist_ok=True)
    summaries, channels, obsrows = [], [], []
    matrices = {}
    failures = []
    for model in requested:
        path = source / tag(model) / 'cases.json'
        try:
            summary, rows, obs, matrix_data = analyse(model, path)
            summaries.append(summary)
            channels.extend(rows)
            obsrows.extend(obs)
            matrices[model] = matrix_data
            print(f'{model}: log={summary["log_norm_GeV_inv"]:.6e} GeV^-1; '
                  f'log/bubble={summary["log_over_bubble_norm"]:.3g}; '
                  f'leading={summary["leading_quartic"]} '
                  f'({100*summary["leading_quartic_projection_fraction"]:.3f}%)')
        except (ValueError, KeyError, FileNotFoundError, TypeError) as error:
            failures.append({'model':model,'error':f'{type(error).__name__}: {error}'})
            print(f'{model}: FAILED: {error}')
    write_csv(dest / 'all16_EFT1_summary.csv', summaries,
              ['model','status','log_norm_GeV_inv','bubble_norm_GeV_inv',
               'scalar_threshold_c5_norm_GeV_inv','log_over_bubble_norm',
               'log_over_scalar_threshold_c5_norm','leading_quartic',
               'leading_quartic_projection_fraction','quartic_channel_count',
               'fermion_source_count','C5_muF_retained_LLSS_sector',
               'C5_muF_full_one_loop','C_LLSS_muF_fitted_tensor',
               'source_to_Weyl_normalisation_certified'])
    write_csv(dest / 'all16_channel_contributions.csv', channels,
              ['model','category','name','log_norm_GeV_inv','bubble_norm_GeV_inv',
               'combined_norm_GeV_inv','log_projection_fraction',
               'bubble_projection_fraction','log_frobenius_fraction'])
    write_csv(dest / 'all16_EFT1_ablation_observables.csv', obsrows,
              ['model','observable','full_1111','LLSS_log_off_1101',
               'renormalisable_EFT1_off_1011','both_EFT1_off_1001',
               'full_minus_log_off','full_minus_renormalisable_off',
               'full_minus_both_off'])
    (dest / 'all16_Weinberg_matrices.json').write_text(
        json.dumps(matrices, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    (dest / 'manifest.json').write_text(json.dumps({
        'source_campaign':str(source), 'requested_models':requested,
        'computed_models':[x['model'] for x in summaries], 'failures':failures,
        'physics_rerun':False, 'fitted_Yukawas_modified':False,
        'scope': 'Existing LLSS leading-log plus scalar-bubble channel decomposition',
        'EFT_sequence': 'UV -> all fermions -> EFT1 -> all scalars -> SMEFT',
        'muF_GeV':1e5, 'muS_GeV':1e3,
        'warnings':[
            'No fitted six-index C_LLSS(muF) tensor is stored in cases.json.',
            'No independently computed hard C5(muF) is available.',
            'The C5(muF)=0 convention applies only to the retained LLSS-sector approximation.',
            'Frobenius percentages are not additive; signed projection fractions are additive within each category.',
            'Quartic and fermion-tree decompositions are separate descriptions of the same amplitudes.',
            'The generated logarithm uses the implemented fixed-order LLSS approximation, not complete EFT1 dimension-five mixing.',
            'Spinor normalization and general scalar basis identification remain uncertified.',
        ]}, indent=2) + '\n',encoding='utf-8')
    print(f'Output: {dest}')
    return 2 if failures else 0


if __name__ == '__main__':
    raise SystemExit(main())
