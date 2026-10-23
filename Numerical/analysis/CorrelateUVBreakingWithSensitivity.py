"""Compare saved UV beta-function scaling breaking with UV neutrino sensitivities.

Read-only input, CSV/JSON output. Does not run RGEs or refit models.
The beta diagnostic is an instantaneous Frobenius-norm rate at the UV boundary;
low-energy sensitivity is an absolute, fitted-baseline local log derivative.
They are not the same observable, and correlations do not demonstrate causality.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OBS = ('delta_m21_sq_ev2', 'delta_m3l_sq_ev2', 'sin2_theta12', 'sin2_theta13', 'sin2_theta23')
MODES = ('optimised', 'shared_a_m2')
DIRECTIONS = ('T1', 'T2', 'T3')
CATEGORIES = ('gauge', 'yukawa', 'scalar', 'mixed', 'other')


def load_csv(path: Path) -> list[dict]:
    if not path.is_file():
        raise FileNotFoundError(f'Missing saved data: {path}')
    with path.open(newline='', encoding='utf-8-sig') as f:
        return list(csv.DictReader(f))


def save_csv(path: Path, rows: list[dict], names: list[str]):
    with path.open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=names, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def corr(x, y, kind='pearson'):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.size < 4 or not np.isfinite(x).all() or not np.isfinite(y).all():
        return None
    if kind == 'spearman':
        x = rank(x)
        y = rank(y)
    xx, yy = x - x.mean(), y - y.mean()
    den = math.sqrt(float(np.dot(xx, xx) * np.dot(yy, yy)))
    return float(np.dot(xx, yy) / den) if den > 1e-25 * max(float(np.dot(x, x)), float(np.dot(y, y)), 1e-30) else None


def rank(a):
    indices = np.argsort(a, kind='mergesort')
    result = np.empty(len(a), dtype=float)
    k = 0
    while k < len(a):
        j = k + 1
        while j < len(a) and a[indices[j]] == a[indices[k]]:
            j += 1
        result[indices[k:j]] = (k + j - 1) / 2 + 1
        k = j
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gamma-minus', type=float, default=0.9)
    parser.add_argument('--gamma-plus', type=float, default=1.1)
    parser.add_argument('--beta-csv', type=Path, default=ROOT / 'output/studies/lambda_5e-9_refit/uv_scaling_term_decomposition/category_breaking.csv')
    parser.add_argument('--sensitivity-csv', type=Path, default=ROOT / 'output/studies/fermion_start_lambda_comparison/uv_sensitivity_comparison/local_slopes_by_model.csv')
    parser.add_argument('--output', type=Path, default=ROOT / 'output/studies/lambda_5e-9_refit/uv_breaking_sensitivity_correlation')
    args = parser.parse_args()
    if not (0 < args.gamma_minus < 1 < args.gamma_plus):
        parser.error('Require 0 < gamma-minus < 1 < gamma-plus.')
    betas = load_csv(args.beta_csv)
    sensitivities = load_csv(args.sensitivity_csv)
    indexed = defaultdict(dict)
    for b in betas:
        key = (b['mode'], b['model'], b['direction'])
        g = float(b['gamma'])
        if math.isclose(g, args.gamma_minus, rel_tol=1e-8, abs_tol=1e-12):
            if 'minus' in indexed[key]:
                raise ValueError(f'Duplicate gamma-minus: {key}')
            indexed[key]['minus'] = b
        if math.isclose(g, args.gamma_plus, rel_tol=1e-8, abs_tol=1e-12):
            if 'plus' in indexed[key]:
                raise ValueError(f'Duplicate gamma-plus: {key}')
            indexed[key]['plus'] = b
    joined = []
    for s in sensitivities:
        key = (s['mode'], s['model'], s['direction'])
        pair = indexed.get(key, {})
        if set(pair) != {'minus', 'plus'}:
            raise ValueError(f'Missing saved beta diagnostics at both requested gamma values for {key}')
        if s['observable'] not in OBS:
            continue
        record = dict(mode=key[0], model=key[1], direction=key[2], observable=s['observable'],
                      uv_abs_log_slope=float(s['uv_abs_slope']),
                      fermion_abs_log_slope=float(s['fermion_abs_slope']),
                      gamma_minus=args.gamma_minus, gamma_plus=args.gamma_plus)
        for col in ('delta_betaP_over_16pi2_P',) + tuple(f'{c}_over_16pi2_P' for c in CATEGORIES):
            values = [float(pair[q][col]) for q in ('minus', 'plus')]
            if not all(math.isfinite(v) and v >= 0 for v in values):
                raise ValueError(f'Invalid breaking norm {key}/{col}: {values}')
            # Finite displacement ratio makes beta-rate diagnostic comparable across selected gamma steps.
            record[col + '_minus'] = values[0]
            record[col + '_plus'] = values[1]
            record[col + '_per_log_gamma_avg'] = sum(values[i] / abs(math.log(g)) for i, g in enumerate((args.gamma_minus, args.gamma_plus))) / 2
        joined.append(record)
    expected = len(MODES) * 16 * len(DIRECTIONS) * len(OBS)
    if len(joined) != expected or len(set((r['mode'],r['model'],r['direction'],r['observable']) for r in joined)) != len(joined):
        raise ValueError(f'Expected {expected} unique joined metric rows, found {len(joined)}')
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    save_csv(out/'joined_by_model.csv', joined, list(joined[0]))
    summary = []
    for mode in MODES:
        for direction in (*DIRECTIONS, 'ALL'):
            for observable in OBS:
                subset = [r for r in joined if r['mode']==mode and r['observable']==observable and (direction=='ALL' or r['direction']==direction)]
                x = [r['delta_betaP_over_16pi2_P_per_log_gamma_avg'] for r in subset]
                y = [r['uv_abs_log_slope'] for r in subset]
                rmax = max(subset, key=lambda r:r['uv_abs_log_slope'])
                bmax = max(subset, key=lambda r:r['delta_betaP_over_16pi2_P_per_log_gamma_avg'])
                summary.append(dict(mode=mode, direction=direction, observable=observable, n=len(subset),
                                    pearson=corr(x,y), spearman=corr(x,y,'spearman'),
                                    max_uv_sensitivity=rmax['uv_abs_log_slope'], max_uv_sensitivity_model=rmax['model'],
                                    max_beta_breaking=bmax['delta_betaP_over_16pi2_P_per_log_gamma_avg'], max_beta_breaking_model=bmax['model']))
    save_csv(out/'correlations.csv', summary, list(summary[0]))
    metadata = {
        'row_count':len(joined), 'gamma_minus':args.gamma_minus,'gamma_plus':args.gamma_plus,
        'beta_input':str(args.beta_csv),'sensitivity_input':str(args.sensitivity_csv),
        'method':'Mean of ||Delta beta_P||/(16pi²||P||)/|ln(gamma)| for two gamma displacements, compared against absolute UV low-energy log slope',
        'important_limits':['The instantaneous beta breaking is at the UV boundary, not integrated to mu_F.',
                            'The response is an absolute low-energy observable slope, and beta breaking is a nonnegative matrix norm; signed directional effects are lost.',
                            'The previous UV refit and fermion-start inputs are separately fitted.',
                            'ALL-direction correlations pool T1/T2/T3 and are descriptive; individual-direction correlations control for direction.',
                            'No formal statistical significance or causal attribution is claimed.']}
    (out/'summary.json').write_text(json.dumps(metadata,indent=2)+'\n',encoding='utf-8')
    print('Joined rows:', len(joined), '| correlations:', len(summary))
    for r in summary:
        if r['direction']=='ALL':
            print(f"{r['mode']:12} {r['observable']:23} Pearson={r['pearson']} Spearman={r['spearman']}")
    print('Saved:',out)


if __name__=='__main__':
    main()
