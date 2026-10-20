"""Independent norm-based study of shared T3 Yukawas (read-only inputs).

Clusters originate from distances between individually fitted Yukawa pairs.
A model counts as supported only if the *same actual donor pair* gives
chi2 <= cutoff in the already evaluated cross-model matrix.
No Yukawa averaging, new RG evaluations, or modification of existing results.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist, squareform
from scipy.stats import chi2 as chi2_distribution

from Numerical.orchestration import RunSmallLambdaRefit as base

DEFAULT_MATRIX = base.STUDY / 'cross_compatibility' / 'chi2_matrix.csv'
OUTPUT = base.STUDY / 'norm_cluster_study'


def read_matrix(path: Path):
    with path.open(newline='', encoding='utf-8-sig') as stream:
        reader = csv.reader(stream)
        names = next(reader)
        data = np.array([[float(item) for item in row] for row in reader if row], float)
    if data.shape != (len(names), len(names)) or len(set(names)) != len(names):
        raise ValueError('Inconsistent cross-evaluation matrix header/shape')
    return names, data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--matrix', type=Path, default=DEFAULT_MATRIX)
    parser.add_argument('--output', type=Path, default=OUTPUT)
    parser.add_argument('--min-models', type=int, default=5)
    parser.add_argument('--confidence', type=float, default=0.95)
    parser.add_argument('--cutoff', type=float)
    parser.add_argument('--radii', type=float, nargs='+', default=[0.05,0.1,0.15,0.2,0.25,0.35,0.5,0.75,1.0])
    args = parser.parse_args()
    if args.min_models < 2 or not 0 < args.confidence < 1 or any(r <= 0 for r in args.radii):
        parser.error('Invalid --min-models, --confidence or --radii')
    cutoff = float(args.cutoff if args.cutoff is not None else chi2_distribution.ppf(args.confidence, 5))
    if not np.isfinite(cutoff) or cutoff <= 0: parser.error('Invalid chi2 cutoff')
    names, chi2 = read_matrix(args.matrix)
    if len(names) != 16: raise ValueError('Expected all 16 models in the cross-evaluation matrix')
    configs = []
    vectors = []
    for name in names:
        path = base.STUDY / 'fitted' / (name + '.json')
        if not path.is_file(): raise FileNotFoundError(f'Missing fitted config: {path}')
        cfg = base.load(path)
        configs.append(cfg)
        ordinary = cfg['base_state']['ordinary']
        vectors.append(np.concatenate([base.matr(ordinary[key]).ravel() for key in ('y1','y2')]))
    vectors = np.asarray(vectors, dtype=complex)
    # Absolute Frobenius metric: sqrt(||Delta Y1||_F^2 + ||Delta Y2||_F^2).
    distances = squareform(pdist(np.concatenate((vectors.real,vectors.imag),axis=1)))
    if not np.allclose(distances, distances.T) or not np.allclose(np.diag(distances), 0):
        raise RuntimeError('Invalid Yukawa distance matrix')
    out = args.output
    out.mkdir(parents=True,exist_ok=True)
    np.savetxt(out/'distance_matrix.csv',distances,delimiter=',',header=','.join(names),comments='')
    rows=[]
    cluster_rows=[]
    for radius in sorted(set(args.radii)):
        # Hierarchical average-linkage clustering, as in the original exploratory study.
        labels=fcluster(linkage(squareform(distances),method='average'),t=radius,criterion='distance')
        for label in sorted(set(labels)):
            members=np.flatnonzero(labels==label)
            # Candidate donors restricted to that norm cluster.
            # Requiring the donor to be in the group keeps interpretation simple.
            options=[]
            for donor in members:
                valid=[int(i) for i in members if np.isfinite(chi2[i,donor]) and chi2[i,donor] <= cutoff]
                valid_scores=[float(chi2[i,donor]) for i in valid]
                options.append((len(valid),sum(valid_scores),int(donor),valid))
            n_good, total, donor, valid=max(options,key=lambda z:(z[0],-z[1],-z[2]))
            cluster_rows.append({'radius':radius,'cluster':int(label),'norm_cluster_size':len(members),
                                 'donor':names[donor],'passing_count':n_good,
                                 'max_passing_chi2':float(max(chi2[i,donor] for i in valid)),
                                 'member_models':';'.join(names[i] for i in members),
                                 'passing_models':';'.join(names[i] for i in valid)})
        # Radius-centred neighbourhoods give an independent, easy-to-explain norm approach.
        for donor in range(len(names)):
            near=np.flatnonzero(distances[:,donor] <= radius+1e-12)
            passed=[int(i) for i in near if np.isfinite(chi2[i,donor]) and chi2[i,donor] <= cutoff]
            passing_chi2=[float(chi2[i,donor]) for i in passed]
            rows.append({'radius':radius,'donor':names[donor], 'norm_neighbours':len(near),
                         'passing_count':len(passed), 'max_passing_chi2':max(passing_chi2),
                         'total_passing_chi2':sum(passing_chi2),
                         'near_models':';'.join(names[i] for i in near),
                         'passing_models':';'.join(names[i] for i in passed)})
    # Pick smallest radius that reaches the target, then most passing, best chi2.
    eligible=[r for r in rows if r['passing_count']>=args.min_models]
    ranked=sorted(eligible if eligible else rows,
                  key=lambda r:(r['radius'] if eligible else -r['passing_count'],
                                -r['passing_count'] if eligible else r['radius'],
                                r['max_passing_chi2'],r['total_passing_chi2'],r['donor']))
    best=ranked[0]
    final={'definition':'combined Frobenius norm between both complex 3x3 Yukawa matrices',
           'selection':'smallest inspected radius giving >= min_models physical passes; donor is an existing fitted pair',
           'chi2_note':'Five-observable nominal chi2 screening cutoff, not a calibrated parameter confidence region',
           'chi2_cutoff':cutoff,'min_models':args.min_models,
           'fixed_lambda':{'A-D':1e-7,'E':5e-8},
           'best':best,'meets_min_models':bool(eligible),
           'backup_unchanged':'cross_compatibility/max_coverage and all fitted solutions',
           'note':'Clusters and neighbourhoods need not be disjoint. Fitted Yukawas use the stored flavour basis.'}
    (out/'summary.json').write_text(json.dumps(final,indent=2)+'\n',encoding='utf-8')
    def write_csv(file,items):
        with (out/file).open('w',newline='',encoding='utf-8') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(items[0]));writer.writeheader();writer.writerows(items)
    write_csv('neighbourhood_scan.csv',rows)
    write_csv('hierarchical_clusters.csv',cluster_rows)
    chosen=names.index(best['donor'])
    comparison=[]
    for i,name in enumerate(names):
        comparison.append({'model':name,'distance_to_donor':float(distances[i,chosen]),
                           'chi2_with_same_yukawas':float(chi2[i,chosen]),
                           'norm_neighbour':bool(distances[i,chosen]<=best['radius']+1e-12),
                           'passes_chi2':bool(np.isfinite(chi2[i,chosen]) and chi2[i,chosen]<=cutoff),
                           'included':bool(distances[i,chosen]<=best['radius']+1e-12 and np.isfinite(chi2[i,chosen]) and chi2[i,chosen]<=cutoff)})
    write_csv('selected_comparison.csv',comparison)
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig,ax=plt.subplots(figsize=(8,5))
        for i,name in enumerate(names):
            ax.scatter(distances[i,chosen],max(chi2[i,chosen],1e-17),marker='o' if comparison[i]['included'] else 'x')
            ax.annotate(str(i+1),(distances[i,chosen],max(chi2[i,chosen],1e-17)),fontsize=7)
        ax.axhline(cutoff,ls='--',color='gray',label=f'chi2 cutoff {cutoff:.2f}')
        ax.axvline(best['radius'],ls=':',color='gray',label=f'norm radius {best["radius"]:g}')
        ax.set(yscale='log',xlabel='Combined Yukawa Frobenius distance',ylabel='Cross-model chi-square',title='Norm similarity versus neutrino compatibility')
        ax.legend();fig.tight_layout();fig.savefig(out/'norm_vs_chi2.png',dpi=180);plt.close(fig)
        fig,ax=plt.subplots(figsize=(8,5))
        radii=sorted(set(row['radius'] for row in rows))
        ax.plot(radii,[max(x['passing_count'] for x in rows if x['radius']==r) for r in radii],marker='o')
        ax.axhline(args.min_models,ls='--',color='gray')
        ax.set(xlabel='Neighbourhood radius',ylabel='Maximum passing models',title='Shared-Yukawa coverage versus norm radius')
        fig.tight_layout();fig.savefig(out/'coverage_vs_radius.png',dpi=180);plt.close(fig)
    except ImportError:
        print('matplotlib unavailable; plot generation skipped')
    print(f'Cutoff: {cutoff:.5g}; target >= {args.min_models} models')
    print(f'Best norm-based candidate: donor={best["donor"]}, radius={best["radius"]:g}')
    print(f'Norm neighbours={best["norm_neighbours"]}, physically passing={best["passing_count"]}/16')
    if not eligible: print('WARNING: None of the inspected radii reaches the required model count')
    for item in comparison:
        if item['included']:print(f' PASS {item["model"]}: d={item["distance_to_donor"]:.4g} chi2={item["chi2_with_same_yukawas"]:.5g}')
    print(f'Files: {out}')

if __name__=='__main__':
    main()
