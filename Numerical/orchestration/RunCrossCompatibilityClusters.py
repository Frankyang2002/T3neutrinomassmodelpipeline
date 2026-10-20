"""Cluster T3 models by full cross-evaluated compatibility of fitted Yukawas.

No joint optimisation is claimed: each cluster uses an existing fitted Yukawa pair
from a member model. A model belongs to that cluster only when its own full
numerical evaluation passes the specified chi2 threshold.
"""
from __future__ import annotations

import argparse
import copy
import csv
import json
from pathlib import Path

import numpy as np

from Numerical.orchestration import RunSmallLambdaRefit as base

ROOT = base.STUDY / 'cross_compatibility'


def load_models():
    models = []
    for row in base.original_rows():
        name = row['model']
        cfg = base.load(base.STUDY / 'fitted' / (name + '.json'))
        ordinary = cfg['base_state']['ordinary']
        fixed_lambda = 5e-8 if base.cls(row) == 'E' else 1e-7
        actual = float(ordinary['lambdaT3']['real'])
        if not np.isclose(actual, fixed_lambda, rtol=1e-10, atol=0):
            raise ValueError(f'{name}: lambdaT3={actual}, expected {fixed_lambda}')
        y = np.concatenate([base.matr(ordinary[key]).real.ravel() for key in ('y1', 'y2')])
        models.append({'name': name, 'row': row, 'config': cfg, 'y': y})
    if len(models) != 16:
        raise ValueError(f'Expected 16 fitted models; found {len(models)}')
    return models


def evaluate_matrix(models, artifact_map):
    n = len(models)
    chi2 = np.full((n, n), np.inf)
    failures = []
    mapping = base.load(artifact_map) if artifact_map else None
    for target, model in enumerate(models):
        artifacts = base.locate_artifacts(model['name'], mapping)
        parameters, evaluator = base.get_evaluator(model['config'], artifacts)
        for donor, source in enumerate(models):
            try:
                result = base.evaluate_point(evaluator, parameters, source['y'])
                residual = np.asarray(result.whitened_residual, dtype=float)
                if residual.shape != (5,) or not np.all(np.isfinite(residual)):
                    raise ValueError(f'Invalid residual {residual!r}')
                chi2[target, donor] = float(residual @ residual)
            except Exception as exc:
                failures.append({'model': model['name'], 'donor': source['name'],
                                 'error': f'{type(exc).__name__}: {exc}'})
        print(f'[{target+1}/{n}] {model["name"]}: '
              f'{np.count_nonzero(np.isfinite(chi2[target]))}/{n} successful probes', flush=True)
    return chi2, failures


def greedy_cover(chi2, cutoff):
    """Cover models with compatible real fitted donor pairs; no forced cluster sizes."""
    unassigned = set(range(chi2.shape[0]))
    clusters = []
    while unassigned:
        candidates = []
        for donor in range(chi2.shape[1]):
            covered = [i for i in sorted(unassigned) if np.isfinite(chi2[i, donor])
                       and chi2[i, donor] < cutoff]
            if not covered:
                continue
            candidates.append((-len(covered), float(sum(chi2[i, donor] for i in covered)),
                               donor, covered))
        if not candidates:
            raise RuntimeError('Some models fail the chi2 cutoff even with their own Yukawas')
        _, _, donor, covered = min(candidates)
        clusters.append({'donor': donor, 'members': covered})
        unassigned.difference_update(covered)
    return clusters


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--chi2-max', type=float, default=10.0,
                        help='Strict upper bound on each model chi2; default 10')
    parser.add_argument('--artifact-map', type=Path)
    args = parser.parse_args()
    if not np.isfinite(args.chi2_max) or args.chi2_max <= 0:
        parser.error('--chi2-max must be finite and positive')
    ROOT.mkdir(parents=True, exist_ok=True)
    models = load_models()
    chi2, failures = evaluate_matrix(models, args.artifact_map)
    np.savetxt(ROOT / 'chi2_matrix.csv', chi2, delimiter=',',
               header=','.join(m['name'] for m in models), comments='')
    base.save(ROOT / 'probe_failures.json', failures)
    groups = greedy_cover(chi2, args.chi2_max)
    summary = {'method': 'Greedy coverage with existing fitted Yukawa donor matrices',
               'chi2_cutoff_strict': args.chi2_max, 'fixed_lambda': {'A-D': 1e-7, 'E': 5e-8},
               'notes': 'Not a joint fit or globally optimal clustering; raw real Yukawa basis.',
               'clusters': []}
    rows = []
    for idx, group in enumerate(groups, start=1):
        donor = group['donor']
        donor_name = models[donor]['name']
        members = group['members']
        summary['clusters'].append({'cluster': idx, 'yukawa_source': donor_name,
                                    'size': len(members),
                                    'members': [models[i]['name'] for i in members],
                                    'max_chi2': float(max(chi2[i, donor] for i in members))})
        print(f'Cluster {idx}: {len(members)} models, source={donor_name}, '
              f'max chi2={max(chi2[i, donor] for i in members):.6g}', flush=True)
        for i in members:
            m = models[i]
            cfg = copy.deepcopy(m['config'])
            ordinary = cfg['base_state']['ordinary']
            for j, key in enumerate(('y1', 'y2')):
                original_imag = base.matr(ordinary[key]).imag
                base.setmat(ordinary[key], models[donor]['y'][9*j:9*(j+1)].reshape(3, 3)
                            + 1j * original_imag)
            cfg.setdefault('cross_compatibility', {})
            cfg['cross_compatibility'].update({'cluster': idx, 'yukawa_source': donor_name,
                                               'chi2': float(chi2[i, donor])})
            base.save(ROOT / 'shared_configs' / (m['name'] + '.json'), cfg)
            rows.append({'model': m['name'], 'class': base.cls(m['row']),
                         'cluster': idx, 'members': len(members), 'yukawa_source': donor_name,
                         'individual_chi2': m['config'].get('small_lambda_study', {}).get('chi2', ''),
                         'shared_chi2': float(chi2[i, donor])})
    base.save(ROOT / 'clusters.json', summary)
    with (ROOT / 'comparison.csv').open('w', newline='', encoding='utf-8') as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        display = np.log10(np.maximum(chi2, 1e-18))
        display[~np.isfinite(display)] = np.nan
        fig, ax = plt.subplots(figsize=(11, 9))
        im = ax.imshow(display, aspect='auto')
        ax.set_xticks(range(len(models)), [m['name'] for m in models], rotation=90, fontsize=5)
        ax.set_yticks(range(len(models)), [m['name'] for m in models], fontsize=5)
        ax.set_xlabel('Donor Yukawa pair')
        ax.set_ylabel('Evaluated model')
        fig.colorbar(im, ax=ax, label='log10(chi2)')
        fig.tight_layout()
        fig.savefig(ROOT / 'compatibility_heatmap.png', dpi=180)
        plt.close(fig)
    except ImportError:
        print('matplotlib not installed: heatmap skipped', flush=True)
    print(f'Saved {len(groups)} clusters covering all {len(models)} models. Files: {ROOT}', flush=True)
    print('No joint minimisation was performed. Inspect comparison.csv before treating these as common-Yukawa families.', flush=True)


if __name__ == '__main__':
    main()
