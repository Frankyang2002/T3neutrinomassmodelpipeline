"""Rank existing fitted Yukawa pairs by number of T3 models they describe.

Reads the authoritative 16x16 cross-evaluation matrix already calculated by
RunCrossCompatibilityClusters; performs no new physics evaluations or fits.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.stats import chi2 as chi2_distribution

PROJECT_ROOT = Path(__file__).resolve().parents[2]
ROOT = PROJECT_ROOT / 'output' / 'studies' / 'small_lambda_1e-7' / 'cross_compatibility'


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--matrix', type=Path, default=ROOT / 'chi2_matrix.csv')
    parser.add_argument('--confidence', type=float, default=0.95,
                        help='Nominal chi-square CDF probability for 5 observables (default 0.95)')
    parser.add_argument('--cutoff', type=float, default=None,
                        help='Override nominal chi-square cutoff with a chosen value')
    args = parser.parse_args()
    if not 0 < args.confidence < 1:
        parser.error('--confidence must be between 0 and 1')
    if not args.matrix.is_file():
        parser.error(f'Missing {args.matrix}. Run RunCrossCompatibilityClusters first.')

    with args.matrix.open('r', newline='', encoding='utf-8-sig') as fh:
        reader = csv.reader(fh)
        names = next(reader)
        matrix = np.asarray([[float(v) for v in row] for row in reader if row], dtype=float)
    n = len(names)
    if matrix.shape != (n, n):
        parser.error(f'Expected NxN chi2 matrix; got {matrix.shape}, N={n}')
    if n != 16 or len(set(names)) != 16:
        parser.error('Expected 16 distinct model names in matrix header')

    cutoff = args.cutoff if args.cutoff is not None else float(chi2_distribution.ppf(args.confidence, df=5))
    if not np.isfinite(cutoff) or cutoff <= 0:
        parser.error('Cutoff must be positive and finite')
    candidates = []
    for j, donor in enumerate(names):
        passing = np.isfinite(matrix[:, j]) & (matrix[:, j] <= cutoff)
        passing_indices = np.flatnonzero(passing)
        candidates.append({
            'donor_model': donor,
            'count': int(len(passing_indices)),
            'fraction': float(len(passing_indices) / n),
            'worst_passing_chi2': float(np.max(matrix[passing, j])) if np.any(passing) else None,
            'total_passing_chi2': float(np.sum(matrix[passing, j])),
            'models_passing': [names[i] for i in passing_indices],
            'chi2_by_model': {names[i]: (float(matrix[i, j]) if np.isfinite(matrix[i, j]) else None)
                              for i in range(n)},
        })
    candidates.sort(key=lambda d: (-d['count'], d['total_passing_chi2'], d['donor_model']))
    output = ROOT / 'max_coverage'
    output.mkdir(parents=True, exist_ok=True)
    result = {
        'method': 'Maximum model coverage among 16 existing individually fitted Yukawa donor pairs',
        'note': ('Nominal 5-observable chi-square cutoff is a screening rule, not a calibrated '
                 '95% confidence statement for fitted parameters or multiple tested donors.'),
        'confidence_input': args.confidence,
        'chi2_degrees_for_screening': 5,
        'chi2_cutoff_inclusive': cutoff,
        'models': names,
        'ranking': candidates,
    }
    (output / 'coverage.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    with (output / 'ranking.csv').open('w', newline='', encoding='utf-8') as fh:
        writer = csv.DictWriter(fh, fieldnames=['rank', 'donor_model', 'count', 'fraction', 'worst_passing_chi2', 'total_passing_chi2'])
        writer.writeheader()
        for rank, c in enumerate(candidates, 1):
            writer.writerow({k: c[k] for k in writer.fieldnames if k != 'rank'} | {'rank': rank})
    best = candidates[0]
    with (output / 'best_model_comparison.csv').open('w', newline='', encoding='utf-8') as fh:
        writer = csv.DictWriter(fh, fieldnames=['model', 'chi2', 'passes_cutoff'])
        writer.writeheader()
        for name in names:
            score = best['chi2_by_model'][name]
            writer.writerow({'model': name, 'chi2': score,
                             'passes_cutoff': score is not None and score <= cutoff})

    print(f'Nominal {100*args.confidence:.1f}% five-observable chi2 cutoff: {cutoff:.6g}')
    print(f'Best single existing Yukawa pair: {best["donor_model"]}')
    print(f'Models passing: {best["count"]}/{n}')
    print('Top donors:')
    for c in candidates[:8]:
        print(f'  {c["count"]:2d}/{n}: {c["donor_model"]} (worst passing chi2={c["worst_passing_chi2"]})')
    print('Models passing with best pair:')
    for model in best['models_passing']:
        print(f'  {model}: chi2={best["chi2_by_model"][model]:.6g}')
    print(f'Files: {output}')
    print('This ranks existing donor matrices only; it does not optimise a new common pair.')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
