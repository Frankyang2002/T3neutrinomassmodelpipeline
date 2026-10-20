"""Budgeted, observable maximum-coverage search for one common real Yukawa pair.

Evaluates interpolations between fitted donor Yukawa pairs. Unlike least_squares,
no finite-difference Jacobian is computed. This is a limited heuristic search,
not a continuous global optimum. Existing physics evaluator is unchanged.

Resume is supported via a per-candidate JSON checkpoint. Interrupt with Ctrl+C.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import time
import numpy as np
from scipy.stats import chi2 as chi2_distribution
from Numerical.orchestration import RunSmallLambdaRefit as base

OUT = base.STUDY / 'maximum_coverage_budgeted'
MATRIX = base.STUDY / 'cross_compatibility' / 'chi2_matrix.csv'


def read_matrix(path, names):
    import csv
    with Path(path).open(newline='', encoding='utf-8-sig') as f:
        rows = list(csv.reader(f))
    header = rows[0]
    a = np.asarray([[float(v) for v in r] for r in rows[1:] if r], dtype=float)
    if len(set(header)) != len(names) or set(header) != set(names) or a.shape != (len(names), len(names)):
        raise ValueError('Invalid 16x16 compatibility matrix or model order')
    ix = [header.index(n) for n in names]
    return a[np.ix_(ix, ix)]


def score_rank(a, cutoff):
    good = a[np.isfinite(a) & (a <= cutoff)]
    return (len(good), -float(good.max()) if len(good) else -float('inf'),
            -float(good.sum()) if len(good) else -float('inf'))


def write_solution(models, vec, scores, cutoff, label):
    import copy
    import csv
    folder = OUT / 'best'
    folder.mkdir(parents=True, exist_ok=True)
    payload = {'origin': label, 'fixed_lambda': {'A-D': 1e-7, 'E': 5e-8},
               'chi2_cutoff': float(cutoff), 'coverage': int(sum(np.isfinite(scores) & (scores <= cutoff))),
               'y1_real': vec[:9].reshape(3, 3).tolist(),
               'y2_real': vec[9:].reshape(3, 3).tolist(),
               'models': [{'model': m['name'], 'chi2': float(v) if np.isfinite(v) else None,
                           'passes': bool(np.isfinite(v) and v <= cutoff)} for m, v in zip(models, scores)],
               'note': 'Nominal five-observable screening, not calibrated parameter confidence'}
    base.save(folder / 'result.json', payload)
    with (folder / 'comparison.csv').open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=['model', 'chi2', 'passes'])
        w.writeheader()
        for r in payload['models']:
            w.writerow(r)
    for m in models:
        cfg = copy.deepcopy(m['config'])
        o = cfg['base_state']['ordinary']
        for j, key in enumerate(('y1', 'y2')):
            base.setmat(o[key], vec[j*9:(j+1)*9].reshape(3,3) + 1j*base.matr(o[key]).imag)
        cfg['maximum_coverage_study'] = {'source': label, 'screening_cutoff': float(cutoff)}
        base.save(folder / 'shared_configs' / (m['name']+'.json'), cfg)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--confidence', type=float, default=.95)
    p.add_argument('--matrix', type=Path, default=MATRIX)
    p.add_argument('--artifact-map', type=Path)
    p.add_argument('--max-proposals', type=int, default=4, help='Maximum fully evaluated new trial Yukawa pairs')
    p.add_argument('--donors', type=int, default=4, help='Top existing donors considered for interpolation')
    p.add_argument('--fractions', default='0.15,0.35', help='Interpolation fractions')
    p.add_argument('--bound', type=float, default=.999)
    args = p.parse_args()
    if not (0 < args.confidence < 1) or args.max_proposals < 0 or args.donors < 2:
        p.error('Invalid confidence or search budget')
    fractions = [float(t) for t in args.fractions.split(',')]
    if not fractions or any(not 0 < t < 1 for t in fractions):
        p.error('fractions must lie between 0 and 1')
    cutoff = float(chi2_distribution.ppf(args.confidence, df=5))
    rows = base.original_rows()
    names = [r['model'] for r in rows]
    configs = [base.load(base.STUDY / 'fitted' / (n+'.json')) for n in names]
    vectors = []
    for row, cfg in zip(rows, configs):
        expected = 5e-8 if base.cls(row) == 'E' else 1e-7
        if not np.isclose(cfg['base_state']['ordinary']['lambdaT3']['real'], expected, atol=0, rtol=1e-10):
            raise ValueError('Incorrect lambda for '+row['model'])
        vectors.append(np.concatenate([base.matr(cfg['base_state']['ordinary'][k]).real.ravel() for k in ('y1','y2')]))
    a = read_matrix(args.matrix, names)
    order = sorted(range(len(names)), key=lambda j: score_rank(a[:, j], cutoff), reverse=True)
    primary = order[0]
    candidates = []
    for j in order[1:min(args.donors, len(order))]:
        for t in fractions:
            candidates.append((j, t))
    candidates = candidates[:args.max_proposals]
    OUT.mkdir(parents=True, exist_ok=True)
    checkpoint_dir = OUT / 'checkpoints'
    checkpoint_dir.mkdir(parents=True, exist_ok=True)
    models = [{'name': n, 'config': cfg} for n, cfg in zip(names, configs)]
    best_vec = vectors[primary].copy()
    best_scores = a[:, primary].copy()
    best_rank = score_rank(best_scores, cutoff)
    best_origin = 'existing donor '+names[primary]
    print(f'Baseline {best_rank[0]}/16: {names[primary]}', flush=True)
    mapping = base.load(args.artifact_map) if args.artifact_map else None
    write_solution(models, best_vec, best_scores, cutoff, best_origin)
    evaluators = [None] * len(models)
    for index, (j, t) in enumerate(candidates, 1):
        name = f'trial_{index:03d}'
        path = checkpoint_dir / (name + '.json')
        vec = (1-t)*vectors[primary] + t*vectors[j]
        if np.max(np.abs(vec)) >= args.bound:
            print(f'{name}: skipped, bounds exceeded', flush=True)
            continue
        if path.is_file():
            saved = base.load(path)
            if saved.get('source') != names[j] or saved.get('fraction') != t:
                raise ValueError('Checkpoint differs from current proposal; delete '+str(path))
            values = np.asarray([float(v) if v is not None else np.inf for v in saved['chi2']], dtype=float)
            if len(values) != len(models):
                raise ValueError('Incomplete checkpoint '+str(path))
            print(f'{name}: restored checkpoint', flush=True)
        else:
            values = np.full(len(models), np.inf)
            print(f'{name}/{len(candidates)}: interpolate {names[primary]} -> {names[j]} t={t}', flush=True)
            start = time.perf_counter()
            for i, m in enumerate(models):
                if evaluators[i] is None:
                    artifacts = base.locate_artifacts(m['name'], mapping)
                    evaluators[i] = base.get_evaluator(m['config'], artifacts)
                params, fn = evaluators[i]
                try:
                    res = base.evaluate_point(fn, params, vec)
                    r = np.asarray(res.whitened_residual, dtype=float)
                    if r.shape != (5,) or not np.all(np.isfinite(r)):
                        raise ValueError('Invalid residual')
                    values[i] = float(r @ r)
                except Exception as exc:
                    print(f'  WARNING {m["name"]}: {type(exc).__name__}: {exc}', flush=True)
                print(f'  [{i+1:02d}/16] {m["name"]}: chi2={values[i]:.5g}; '
                      f'elapsed={time.perf_counter()-start:.1f}s', flush=True)
            base.save(path, {'source': names[j], 'fraction': t,
                             'chi2': [float(v) if np.isfinite(v) else None for v in values]})
        trial_rank = score_rank(values, cutoff)
        print(f'{name}: coverage={trial_rank[0]}/16 (best={best_rank[0]}/16)', flush=True)
        if trial_rank > best_rank:
            best_vec, best_scores, best_rank = vec, values, trial_rank
            best_origin = f'{name} donor={names[j]} fraction={t}'
            write_solution(models, best_vec, best_scores, cutoff, best_origin)
            print(f'NEW BEST {best_rank[0]}/16 (saved)', flush=True)
    print(f'Finished; best={best_rank[0]}/16, source={best_origin}; files={OUT}', flush=True)
    print('This is a finite budget interpolation scan; no global-optimum claim.', flush=True)


if __name__ == '__main__':
    main()
