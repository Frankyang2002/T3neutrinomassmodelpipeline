"""Jointly optimise real Yukawa matrices across balanced T3 model clusters.

Uses already fitted small-lambda configs, existing model-specific RGE/matching
artifacts, and the authoritative evaluator from RunSmallLambdaRefit.
Does not change the physics implementation or overwrite the original fits.
"""
from __future__ import annotations

import argparse
import copy
import csv
import json
from pathlib import Path

import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.optimize import least_squares
from scipy.spatial.distance import pdist

from Numerical.orchestration import RunSmallLambdaRefit as base

OUTPUT = base.STUDY / 'joint_clusters'


def load_models():
    models = []
    for row in base.original_rows():
        name = row['model']
        path = base.STUDY / 'fitted' / (name + '.json')
        if not path.is_file():
            raise FileNotFoundError(f'Missing fitted solution: {path}. Run the previous fit stage first.')
        cfg = base.load(path)
        o = cfg['base_state']['ordinary']
        y = np.r_[base.matr(o['y1']).real.ravel(), base.matr(o['y2']).real.ravel()]
        expected = 5e-8 if base.cls(row) == 'E' else 1e-7
        if not np.isclose(float(o['lambdaT3']['real']), expected, rtol=1e-10, atol=0):
            raise ValueError(f'{name}: incorrect fixed lambda')
        models.append({'name': name, 'row': row, 'config': cfg, 'y': y})
    return models


def initial_groups(points, k, min_members, max_members):
    n = len(points)
    dist = pdist(points)
    labels = np.asarray(fcluster(linkage(dist, method='average'), k, criterion='maxclust'), dtype=int)
    groups = [list(np.flatnonzero(labels == label)) for label in sorted(set(labels))]
    if len(groups) != k:
        raise RuntimeError('Hierarchical clustering yielded fewer clusters than requested')
    # Repair sizes, assigning each transferred model to its nearest eligible centroid.
    for _ in range(10*n):
        sizes = [len(g) for g in groups]
        too_large = [i for i, size in enumerate(sizes) if size > max_members]
        too_small = [i for i, size in enumerate(sizes) if size < min_members]
        if not too_large and not too_small:
            return groups
        if too_small:
            recipient = too_small[0]
            donors = [i for i in range(k) if i != recipient and len(groups[i]) > min_members]
            if not donors:
                raise RuntimeError('Cannot satisfy minimum cluster sizes')
        else:
            donors = too_large
            recipient = -1
        choices = []
        for donor in donors:
            for idx in groups[donor]:
                targets = [recipient] if recipient >= 0 else [j for j in range(k) if j != donor and len(groups[j]) < max_members]
                for dest in targets:
                    if not groups[dest]:
                        cost = 0.0
                    else:
                        centroid = np.mean(points[groups[dest]], axis=0)
                        cost = float(np.linalg.norm(points[idx] - centroid))
                    choices.append((cost, donor, dest, idx))
        if not choices:
            raise RuntimeError('Could not rebalance clusters')
        _, donor, dest, idx = min(choices)
        groups[donor].remove(idx)
        groups[dest].append(idx)
    raise RuntimeError('Cluster balancing did not converge')


def evaluate_safe(ev, params, x):
    try:
        result = base.evaluate_point(ev, params, x)
        resid = np.asarray(result.whitened_residual, dtype=float)
        if resid.shape != (5,) or not np.all(np.isfinite(resid)):
            raise ValueError('Invalid residual')
        return resid
    except Exception:
        return np.full(5, 1e5)


def joint_fit(groups, models, evaluators, max_nfev):
    fitted = []
    for group_index, members in enumerate(groups):
        # Medoid is a better nonlinear starting point than the arithmetic mean.
        # Determine it using the sum of actual full-model residual norms.
        options = [models[i]['y'] for i in members]
        scores = []
        for vector in options:
            scores.append(sum(float(np.dot(r, r)) for i in members
                              for r in [evaluate_safe(*evaluators[i], vector)]))
        x0 = options[int(np.argmin(scores))].copy()
        x0 = np.clip(x0, -0.999, 0.999)
        def residual(x):
            return np.concatenate([evaluate_safe(*evaluators[i], x) for i in members])
        result = least_squares(residual, x0,
                               bounds=(-np.ones(18), np.ones(18)),
                               max_nfev=max_nfev,
                               diff_step=1e-4, xtol=1e-8, ftol=1e-8, gtol=1e-8)
        per_model = {}
        for i in members:
            r = evaluate_safe(*evaluators[i], result.x)
            per_model[models[i]['name']] = float(np.dot(r, r))
        print(f'Cluster {group_index + 1}: {len(members)} models, '
              f'max chi2={max(per_model.values()):.5g}, '
              f'total chi2={sum(per_model.values()):.5g}, '
              f'evaluations={result.nfev}', flush=True)
        fitted.append({'vector': result.x, 'per_model': per_model,
                       'converged': bool(result.success), 'nfev': int(result.nfev),
                       'initial_total_chi2': float(scores[int(np.argmin(scores))])})
    return fitted


def assignment_improvement(groups, fitted, models, evaluators, min_members, max_members):
    # One conservative greedy relocation pass: accept only reductions in total chi2
    # with cluster sizes maintained. This is not guaranteed globally optimal.
    k = len(groups)
    costs = np.zeros((len(models), k))
    for i in range(len(models)):
        for j in range(k):
            r = evaluate_safe(*evaluators[i], fitted[j]['vector'])
            costs[i, j] = float(np.dot(r, r))
    updated = [list(g) for g in groups]
    changes = 0
    while True:
        proposals = []
        for src in range(k):
            if len(updated[src]) <= min_members:
                continue
            for idx in updated[src]:
                for dst in range(k):
                    if dst == src or len(updated[dst]) >= max_members:
                        continue
                    delta = costs[idx, dst] - costs[idx, src]
                    if delta < -1e-7:
                        proposals.append((delta, idx, src, dst))
        if not proposals:
            break
        _, idx, src, dst = min(proposals)
        updated[src].remove(idx)
        updated[dst].append(idx)
        changes += 1
    return updated, changes


def write_results(groups, fits, models, evaluations, k, args):
    folder = OUTPUT / f'k{k}'
    folder.mkdir(parents=True, exist_ok=True)
    manifest = {'method': 'balanced hierarchical initialisation + joint least-squares + constrained relocations',
                'k': k, 'min_members': args.min_members, 'max_members': args.max_members,
                'fixed_lambda': {'A-D': 1e-7, 'E': 5e-8},
                'warning': 'Raw Yukawa basis, no rephasing or basis equivalence; fitted cluster solutions are not guaranteed global minima',
                'clusters': []}
    comparison = []
    for j, members in enumerate(groups):
        fit = fits[j]
        vector = fit['vector']
        entry = {'cluster': j+1, 'models': [models[i]['name'] for i in members],
                 'members': len(members), 'common_yukawa_real': {
                     'y1': vector[:9].reshape(3, 3).tolist(), 'y2': vector[9:].reshape(3, 3).tolist()},
                 'chi2_by_model': fit['per_model'],
                 'total_chi2': sum(fit['per_model'].values()),
                 'max_chi2': max(fit['per_model'].values()),
                 'optimizer_converged': fit['converged'], 'nfev': fit['nfev']}
        manifest['clusters'].append(entry)
        for i in members:
            m = models[i]
            cfg = copy.deepcopy(m['config'])
            o = cfg['base_state']['ordinary']
            for z, key in enumerate(('y1', 'y2')):
                original_imag = base.matr(o[key]).imag
                base.setmat(o[key], vector[z*9:(z+1)*9].reshape(3, 3) + 1j*original_imag)
            cfg['small_lambda_study']['status'] = 'Jointly optimised common Yukawa, evaluated against this model'
            cfg['small_lambda_study']['joint_cluster'] = j+1
            cfg['small_lambda_study']['common_chi2'] = fit['per_model'][m['name']]
            base.save(folder / 'common_configs' / (m['name'] + '.json'), cfg)
            individual = float(m['config']['small_lambda_study']['chi2'])
            comparison.append({'model': m['name'], 'class': base.cls(m['row']),
                               'cluster': j+1, 'cluster_size': len(members),
                               'individual_chi2': individual,
                               'common_chi2': fit['per_model'][m['name']],
                               'common_pass_chi2_10': fit['per_model'][m['name']] < 10,
                               'max_abs_y1': float(np.max(np.abs(vector[:9]))),
                               'max_abs_y2': float(np.max(np.abs(vector[9:])) )})
    base.save(folder / 'joint_clusters.json', manifest)
    with (folder / 'comparison.csv').open('w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=list(comparison[0])); w.writeheader(); w.writerows(comparison)
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        ordered = sorted(comparison, key=lambda r:(r['cluster'], r['model']))
        fig, ax = plt.subplots(figsize=(13, 5))
        ax.bar(np.arange(len(ordered)), [max(1e-15, r['common_chi2']) for r in ordered])
        ax.axhline(10, color='red', ls='--', label='chi2 = 10 illustrative threshold')
        ax.set_yscale('log')
        ax.set_xticks(np.arange(len(ordered)), [r['model'] for r in ordered], rotation=90, fontsize=6)
        ax.set_ylabel('Per-model chi-square using joint cluster Yukawas')
        ax.legend(); fig.tight_layout()
        fig.savefig(folder / 'common_chi2.png', dpi=160); plt.close(fig)
        for entry in manifest['clusters']:
            fig, axs = plt.subplots(1, 2, figsize=(8, 3.8), layout='constrained')
            for ax, y in zip(axs, ('y1', 'y2')):
                matrix = np.asarray(entry['common_yukawa_real'][y])
                im = ax.imshow(np.abs(matrix), vmin=0, vmax=max(.1, np.max(np.abs(matrix))))
                ax.set_title(f'Cluster {entry["cluster"]}: |{y}|')
                for ii in range(3):
                    for jj in range(3):
                        ax.text(jj, ii, f'{abs(matrix[ii,jj]):.2g}', ha='center', va='center', fontsize=8)
                fig.colorbar(im, ax=ax)
            fig.savefig(folder / f'cluster_{entry["cluster"]}_yukawas.png', dpi=150)
            plt.close(fig)
    except ImportError:
        print('matplotlib unavailable; skipped plots', flush=True)
    passing = sum(x['common_pass_chi2_10'] for x in comparison)
    worst = max(x['common_chi2'] for x in comparison)
    print(f'k={k}: {passing}/16 models at chi2<10; worst chi2={worst:.5g}. Results: {folder}', flush=True)
    return passing, worst


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--clusters', type=int, default=3, choices=(3,4))
    p.add_argument('--try-four', action='store_true', help='Also attempt four clusters if 3 cluster result has chi2>=10')
    p.add_argument('--min-members', type=int, default=3)
    p.add_argument('--max-members', type=int, default=7)
    p.add_argument('--max-nfev', type=int, default=100)
    p.add_argument('--rounds', type=int, default=2)
    p.add_argument('--artifact-map', type=Path)
    args = p.parse_args()
    if args.max_nfev < 1 or args.rounds < 1:
        p.error('--max-nfev and --rounds must be positive')
    models = load_models()
    mapping = base.load(args.artifact_map) if args.artifact_map else None
    evaluations = []
    for model in models:
        artifacts = base.locate_artifacts(model['name'], mapping)
        evaluations.append(base.get_evaluator(model['config'], artifacts))
    points = np.asarray([m['y'] for m in models])
    for k in ([args.clusters, 4] if args.try_four and args.clusters == 3 else [args.clusters]):
        if k * args.min_members > len(models) or k * args.max_members < len(models):
            p.error('Infeasible cluster-size constraints')
        groups = initial_groups(points, k, args.min_members, args.max_members)
        fits = None
        for step in range(args.rounds):
            print(f'Optimising k={k}, round={step+1}/{args.rounds}, sizes={[len(g) for g in groups]}', flush=True)
            fits = joint_fit(groups, models, evaluations, args.max_nfev)
            if step != args.rounds-1:
                updated, moved = assignment_improvement(groups, fits, models, evaluations,
                                                        args.min_members, args.max_members)
                print(f'Reassigned {moved} models (subject to size constraints)', flush=True)
                if not moved: break
                groups = updated
        passed, worst = write_results(groups, fits, models, evaluations, k, args)
        if args.try_four and passed == len(models):
            break


if __name__ == '__main__':
    main()
