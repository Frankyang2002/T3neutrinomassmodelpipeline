"""Local 5e-9 T3 refit, donor coverage and Yukawa clustering.

Usage: python -m Numerical.orchestration.RunLambda5e9Refit {prepare,fit,cross,analyse,all}
No symbolic regeneration and no writes to historical 1e-7 study.
"""
from __future__ import annotations
import argparse
import copy
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from scipy.stats import chi2 as chi2_dist
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist, squareform

from Numerical.orchestration import RunSmallLambdaRefit as base

OUT = base.ROOT / 'output' / 'studies' / 'lambda_5e-9_refit'
OLD = base.STUDY
LAMBDA = {'normal': 5e-9, 'E': 2.5e-9}
FIELDS = ('delta_m21_sq_ev2', 'delta_m3l_sq_ev2', 'sin2_theta12', 'sin2_theta13', 'sin2_theta23')


def target_lambda(row):
    return LAMBDA['E' if base.cls(row) == 'E' else 'normal']


def name_path(folder, model):
    return OUT / folder / (model + '.json')


def load_models(artifact_map=None):
    mapping = base.load(artifact_map) if artifact_map else None
    return [(r, base.locate_artifacts(r['model'], mapping)) for r in base.original_rows()]


def prepare():
    """Seed from the previous *individual fits*, rescaled to preserve leading C5.

    Previous fits are starting guesses, NOT accepted new solutions.
    """
    rows = base.original_rows()
    summary = []
    for row in rows:
        name = row['model']
        src = OLD / 'fitted' / (name + '.json')
        if not src.is_file():
            raise FileNotFoundError(f'Missing old fitted config: {src}')
        cfg = base.load(src)
        if cfg.get('representation') != row['representation']:
            raise ValueError(f'{name}: representation mismatch')
        ordinary = cfg['base_state']['ordinary']
        previous = float(ordinary['lambdaT3']['real'])
        target = target_lambda(row)
        if previous <= 0:
            raise ValueError(f'{name}: prior lambda must be positive')
        factor = float(np.sqrt(previous / target))
        for key in ('y1', 'y2'):
            x = base.matr(ordinary[key])
            if np.max(np.abs(x.imag)) > 1e-12:
                raise ValueError(f'{name}: expected real fitted Yukawas')
            base.setmat(ordinary[key], x * factor)
        ordinary['lambdaT3'] = {'real': target, 'imag': 0.0}
        cfg['benchmark_search'] = {'enabled': False}
        cfg['sensitivity'] = {'enabled': False}
        cfg.pop('comparison_point', None)
        cfg.pop('small_lambda_study', None)
        cfg['lambda_5e9_study'] = {
            'status': 'Leading-order rescaled seed only; requires fresh 5e-9 evaluation/refinement',
            'source': str(src.relative_to(base.ROOT)),
            'source_sha256': hashlib.sha256(src.read_bytes()).hexdigest(),
            'lambdaT3': target, 'seed_scale': factor,
        }
        base.save(name_path('prepared', name), cfg)
        summary.append({'model': name, 'class': base.cls(row), 'lambdaT3': target, 'scale': factor,
                        'largest_seed_entry': max(float(np.max(np.abs(base.matr(ordinary[k])))) for k in ('y1', 'y2'))})
    base.save(OUT / 'prepared_manifest.json', summary)
    print('Prepared', len(summary), 'seeds at', OUT)


def get_vector(cfg):
    return np.r_[base.matr(cfg['base_state']['ordinary']['y1']).real.ravel(),
                 base.matr(cfg['base_state']['ordinary']['y2']).real.ravel()]


def put_vector(cfg, v):
    for i, key in enumerate(('y1', 'y2')):
        spec = cfg['base_state']['ordinary'][key]
        base.setmat(spec, v[9*i:9*(i+1)].reshape(3, 3) + 1j * base.matr(spec).imag)


def fit(args):
    """Original Sobol -> best five seeds -> sensitivity-ranked 5/8/12 refinement.

    Always evaluate the exact rescaled seed and preserve previous candidates.
    The old gate/expansion flags are accepted by CLI for compatibility but unused.
    """
    from Numerical.fitting.SobolBenchmarkSearch import (
        BenchmarkParameter, generate_sobol_parameter_points,
        refine_sobol_candidates, _build_evaluator,
    )

    records = []
    for i, (row, artifacts) in enumerate(load_models(args.artifact_map), 1):
        model = row['model']
        if args.model and args.model != model:
            continue
        result_path = name_path('fits', model)
        fitted_path = name_path('fitted', model)
        if args.resume and result_path.is_file() and fitted_path.is_file():
            prior = base.load(result_path)
            if prior.get('status') == 'Success' and float(prior.get('chi2_fit', np.inf)) <= args.accept_chi2:
                print(f'[{i}/16] reuse acceptable {model}: chi2={prior["chi2_fit"]:.7g}', flush=True)
                records.append(prior)
                continue

        prepared_path = name_path('prepared', model)
        cfg = base.load(prepared_path)
        if not np.isclose(float(cfg['base_state']['ordinary']['lambdaT3']['real']), target_lambda(row), rtol=1e-12):
            raise ValueError(f'{model}: wrong lambda in prepared configuration')
        seed = get_vector(cfg)
        if np.any(np.abs(seed) >= args.bound):
            raise ValueError(f'{model}: seed exceeds --bound={args.bound:g}; increase --bound')
        names = [f'{key}_{p}{q}' for key in ('y1', 'y2') for p in range(3) for q in range(3)]
        bindings = [f'ordinary.{key}.real.{p}.{q}' for key in ('y1','y2') for p in range(3) for q in range(3)]
        full_defs = [BenchmarkParameter(n, binding, -args.bound, args.bound)
                     for n, binding in zip(names, bindings)]
        half = np.maximum(args.sobol_fraction*np.abs(seed), args.sobol_floor)
        lo = np.maximum(-args.bound, seed-half)
        hi = np.minimum(args.bound, seed+half)
        local_defs = [BenchmarkParameter(n, binding, float(l), float(h))
                      for n, binding, l, h in zip(names, bindings, lo, hi)]
        ev = _build_evaluator(payload=cfg, parameters=full_defs,
              uv_rgbeta_path=artifacts[0], intermediate_rgbeta_path=artifacts[1],
              final_weinberg_path=artifacts[2], target_path=base.TARGET)

        sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
        signature = {'algorithm':'original_sobol_multistart_no_gate_v1',
                     'model':model, 'prepared_sha256':sha(prepared_path),
                     'artifacts_sha256':[sha(x) for x in artifacts],
                     'target_sha256':sha(base.TARGET),
                     'bound':args.bound, 'sobol_fraction':args.sobol_fraction,
                     'sobol_floor':args.sobol_floor, 'sobol_seed':args.sobol_seed+i}
        checkpoint = name_path('sobol_original', model)
        previous = base.load(checkpoint) if args.resume and checkpoint.is_file() else None
        if previous and previous.get('signature') != signature:
            raise RuntimeError(f'{checkpoint}: incompatible run settings; refusing resume')
        evaluated = int(previous.get('evaluated',0)) if previous else 0
        candidates = list(previous.get('candidates',[])) if previous else []

        def inspect(values, source):
            result = ev(values)
            chi = float(result.chi2)
            if not np.isfinite(chi):
                raise ValueError('Non-finite chi2')
            return {'source':source, 'chi2':chi,
                    'parameters':{k:float(v) for k,v in values.items()},
                    'ordering':str(result.ordering),
                    'prediction':dict(zip(result.observable_names,map(float,result.prediction))),
                    'whitened_residual':np.asarray(result.whitened_residual,dtype=float).tolist()}

        def retain(candidate):
            candidates.append(candidate)
            candidates.sort(key=lambda z:float(z['chi2']))
            del candidates[max(args.keep_best,args.seed_count):]

        def write_checkpoint():
            base.save(checkpoint, {'signature':signature, 'evaluated':evaluated,
                    'candidates':candidates, 'best_chi2':candidates[0]['chi2'] if candidates else None})

        # Re-evaluate exact starting seed, and any earlier fitted Yukawa result.
        if not any(c.get('source')=='exact_rescaled_seed' for c in candidates):
            try: retain(inspect(dict(zip(names,seed.tolist())), 'exact_rescaled_seed'))
            except Exception as exc: print(f'[{i}/16] exact seed failed: {exc}',flush=True)
        if args.resume and fitted_path.is_file() and not any(c.get('source')=='earlier_fit' for c in candidates):
            prior_vec = get_vector(base.load(fitted_path))
            if np.all(np.abs(prior_vec) < args.bound):
                try: retain(inspect(dict(zip(names,prior_vec.tolist())), 'earlier_fit'))
                except Exception as exc: print(f'[{i}/16] earlier fit failed: {exc}',flush=True)
        # Import best saved adaptive-search candidates as well, without assuming
        # a previous checkpoint corresponds to the same physics/parameter ranges.
        prior_checkpoint = name_path('sobol_search', model)
        if args.resume and prior_checkpoint.is_file() and not any(c.get('source')=='adaptive_checkpoint' for c in candidates):
            old = base.load(prior_checkpoint)
            sig = old.get('signature',{})
            compatible = (sig.get('model')==model and
                          sig.get('prepared_sha256')==sha(prepared_path) and
                          sig.get('artifact_sha256')==signature['artifacts_sha256'] and
                          sig.get('target_sha256')==signature['target_sha256'])
            if compatible:
                for old_c in old.get('candidates',[]):
                    vals = old_c.get('parameters',{})
                    if all(n in vals and abs(float(vals[n])) < args.bound for n in names):
                        try: retain(inspect({n:float(vals[n]) for n in names},'adaptive_checkpoint'))
                        except Exception: pass
        write_checkpoint()
        print(f'[{i}/16] {model}: seed/best saved chi2={candidates[0]["chi2"]:.7g}' if candidates
              else f'[{i}/16] {model}: no valid starting point',flush=True)

        # One scrambled Sobol prefix, checkpointed every 64 points; never gate refinement.
        points=generate_sobol_parameter_points(local_defs,args.sobol_points,seed=args.sobol_seed+i)
        for j in range(evaluated,len(points)):
            try: retain(inspect(points[j],f'sobol_local_{j+1}'))
            except Exception as exc:
                if j < 3: print(f'[{i}/16] Sobol {j+1} failed: {exc}',flush=True)
            evaluated=j+1
            if evaluated%min(64,args.sobol_points)==0 or evaluated==len(points):
                write_checkpoint()
                print(f'[{i}/16] {model}: Sobol {evaluated}/{len(points)} best chi2={candidates[0]["chi2"]:.7g}',flush=True)
        if not candidates:
            raise RuntimeError(f'{model}: no valid Sobol or rescaled seed')
        best=candidates[0]
        if best['chi2'] <= args.accept_chi2:
            history={'skipped':True,'reason':'Existing seed reaches target','stages':[]}
        else:
            best,history=refine_sobol_candidates(
                evaluator=ev,observable_names=FIELDS,definitions=full_defs,
                candidates=candidates,target_chi2=args.accept_chi2,
                seed_count=args.seed_count,active_counts=(5,8,12),
                max_nfev=(30,40,60),sensitivity_step=0.02)
        final=ev(best['parameters'])
        fitted=copy.deepcopy(cfg)
        put_vector(fitted,np.asarray([best['parameters'][n] for n in names]))
        score=float(final.chi2)
        fitted['lambda_5e9_study'].update(status='Original Sobol multistart staged refinement',chi2=score)
        record={'model':model,'class':base.cls(row),
                'status':'Success' if score<=args.accept_chi2 else 'PoorFit',
                'lambdaT3':target_lambda(row),'chi2_fit':score,
                'fit_method':'original_sobol_multistart_5_8_12_no_gate',
                'passes_screening':score<=args.accept_chi2,
                'sobol_evaluated':evaluated,'refinement':history,
                'prediction':dict(zip(final.observable_names,map(float,final.prediction)))}
        base.save(fitted_path,fitted)
        base.save(result_path,record)
        records.append(record)
        base.save(OUT/'fit_summary.json',records)
        print(f'[{i}/16] {model}: chi2={score:.7g} pass={record["passes_screening"]}',flush=True)
    print('Fit records:',len(records),'at',OUT)

def cross(args):
    for row in base.original_rows():
        report = name_path('fits', row['model'])
        if not report.is_file() or base.load(report).get('status') != 'Success':
            raise RuntimeError(f"Cannot cross-evaluate: {row['model']} lacks an acceptable fit. Complete fitting or explicitly select a different donor policy.")
    models = load_models(args.artifact_map)
    names = [r['model'] for r, _ in models]
    fitted = [base.load(name_path('fitted', n)) for n in names]
    y = np.vstack([get_vector(d) for d in fitted])
    matrix = np.full((len(names), len(names)), np.nan)
    errors = []
    folder = OUT / 'cross_compatibility'
    folder.mkdir(parents=True, exist_ok=True)
    for i, (row, artifacts) in enumerate(models, 1):
        model = row['model']
        path = folder / 'rows' / (model + '.json')
        if args.resume and path.is_file():
            prior = base.load(path)
            if prior.get('names') == names and len(prior.get('chi2', [])) == len(names):
                matrix[i-1] = [np.nan if v is None else float(v) for v in prior['chi2']]
                print(f'[{i}/16] reuse cross {model}', flush=True)
                continue
        params, evaluator = base.get_evaluator(fitted[i-1], artifacts)
        scores = []
        for j, donor in enumerate(names):
            try:
                fit_result = base.evaluate_point(evaluator, params, y[j])
                val = float(fit_result.chi2)
                if not np.isfinite(val):
                    raise ValueError('Non-finite chi2')
                scores.append(val)
            except Exception as exc:
                scores.append(None)
                errors.append({'target': model, 'donor': donor, 'error': str(exc)})
        matrix[i-1] = [np.nan if x is None else x for x in scores]
        base.save(path, {'model': model, 'names': names, 'chi2': scores})
        print(f'[{i}/16] cross {model}: {sum(v is not None for v in scores)}/16 successful', flush=True)
    with (folder / 'chi2_matrix.csv').open('w', newline='', encoding='utf-8') as f:
        w = csv.writer(f); w.writerow(names)
        w.writerows(matrix.tolist())
    base.save(folder / 'failures.json', errors)
    print('Cross evaluation complete: finite', np.isfinite(matrix).sum(), '/ 256')


def analyse(args):
    names = [r['model'] for r in base.original_rows()]
    fitted = [base.load(name_path('fitted', name)) for name in names]
    y = np.vstack([get_vector(c) for c in fitted])
    folder = OUT / 'cross_compatibility'
    with (folder / 'chi2_matrix.csv').open('r', newline='', encoding='utf-8-sig') as f:
        reader = csv.reader(f)
        header = next(reader)
        if header != names:
            raise ValueError('Cross matrix header differs from manifest order')
        matrix = np.array([[float(v) for v in row] for row in reader], dtype=float)
    if matrix.shape != (16, 16):
        raise ValueError(f'Expected 16x16; got {matrix.shape}')
    cutoff_95 = float(chi2_dist.ppf(.95, 5))
    cutoffs = {'strict_5p86': 5.86, 'nominal_95': cutoff_95}
    ranking = {}
    for label, cutoff in cutoffs.items():
        ranked = []
        for j, donor in enumerate(names):
            passed = np.flatnonzero(np.isfinite(matrix[:, j]) & (matrix[:, j] <= cutoff))
            ranked.append({'donor': donor, 'count': len(passed),
                           'max_passing_chi2': float(np.max(matrix[passed,j])) if len(passed) else None,
                           'total_passing_chi2': float(np.sum(matrix[passed,j])),
                           'passing_models': [names[k] for k in passed]})
        ranked.sort(key=lambda v: (-v['count'], v['total_passing_chi2'], v['donor']))
        ranking[label] = {'cutoff': cutoff, 'ranking': ranked}
        print(label, 'best:', ranked[0]['donor'], ranked[0]['count'], '/16', flush=True)
    distances = squareform(pdist(y, metric='euclidean'))
    hclust = fcluster(linkage(pdist(y), method='average'), t=args.distance, criterion='distance')
    clusters = {}
    for name, label in zip(names, hclust):
        clusters.setdefault(str(int(label)), []).append(name)
    # Donor selection uses full physical chi2 matrix, never cluster averages.
    radius_rows = []
    for j, donor in enumerate(names):
        for radius in (0.02, 0.05, 0.10, 0.20, 0.35, 0.50, 1.0, 2.0):
            near = np.flatnonzero(distances[:,j] <= radius + 1e-12)
            passed = [int(k) for k in near if np.isfinite(matrix[k,j]) and matrix[k,j] <= cutoff_95]
            radius_rows.append({'donor': donor, 'radius': radius, 'near_count': len(near),
                                'physical_pass_count': len(passed),
                                'physical_pass_models': [names[k] for k in passed]})
    radius_rows.sort(key=lambda r: (-r['physical_pass_count'], r['radius'], r['donor']))
    output = {'lambdaT3': LAMBDA, 'names': names, 'cutoffs': ranking,
              'norm_clusters': {'method': 'average linkage raw Yukawa Frobenius distance',
                                'distance_cutoff': args.distance, 'memberships': clusters},
              'norm_neighbourhood_candidates': radius_rows,
              'notes': ['New individually fitted Yukawas at lambda=5e-9 (class E 2.5e-9)',
                        'Donor cross-compatibility tested by full evaluator',
                        'Screening thresholds are not calibrated model confidence regions',
                        'No claim of continuous/global common Yukawa optimum']}
    base.save(OUT / 'analysis.json', output)
    with (OUT / 'donor_ranking.csv').open('w', newline='', encoding='utf-8') as f:
        w = csv.writer(f);w.writerow(['criterion','rank','donor','passing_count','max_passing_chi2','total_passing_chi2'])
        for key,v in ranking.items():
            for rank,r in enumerate(v['ranking'],1):
                w.writerow([key,rank,r['donor'],r['count'],r['max_passing_chi2'],r['total_passing_chi2']])
    print('Results:', OUT)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=('prepare', 'fit', 'cross', 'analyse', 'all'))
    parser.add_argument('--artifact-map', type=Path)
    parser.add_argument('--model', help='Exact model identifier (fit only)')
    parser.add_argument('--bound', type=float, default=4.0, help='Absolute entry bound for Yukawa optimisation')
    parser.add_argument('--max-nfev', type=int, default=300)
    parser.add_argument('--accept-chi2', type=float, default=5.86, help='Early-stop screening cutoff (not a confidence region)')
    parser.add_argument('--distance', type=float, default=0.35)
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--sobol-points', type=int, default=64, help='Power of two: 32,64,128,...')
    parser.add_argument('--sobol-max', type=int, default=1024, help='Total Sobol evaluation cap per model (power of two)')
    parser.add_argument('--refine-gate', type=float, default=20.0, help='Only refine if best seed chi2 <= this heuristic gate')
    parser.add_argument('--sobol-seed', type=int, default=20261021)
    parser.add_argument('--sobol-fraction', type=float, default=0.5, help='Fractional Sobol seed neighbourhood')
    parser.add_argument('--sobol-floor', type=float, default=0.1, help='Minimum per-entry Sobol half-width')
    parser.add_argument('--expand-bounds', type=float, nargs='+', default=[2.,4.], help='Progressively wider exploratory global +- Yukawa search boxes')
    parser.add_argument('--expand-if-above', type=float, default=1000., help='Expand when best chi2 remains above this after a 64-point region')
    parser.add_argument('--seed-count', type=int, default=5)
    parser.add_argument('--keep-best', type=int, default=8)
    parser.add_argument('--original-c5', action='store_true', help='Use original symbolic C5 evaluator for validation')
    args = parser.parse_args()
    if args.sobol_points < 2 or args.sobol_points & (args.sobol_points-1):
        parser.error('--sobol-points must be a power of two >=2')
    if args.sobol_max < args.sobol_points or args.sobol_max & (args.sobol_max-1):
        parser.error('--sobol-max must be a power of two >= --sobol-points')
    if args.sobol_max % args.sobol_points:
        parser.error('--sobol-max must be divisible by --sobol-points')
    if not np.isfinite(args.refine_gate) or args.refine_gate < args.accept_chi2:
        parser.error('--refine-gate must be finite and >= --accept-chi2')
    if args.bound <= 0:
        parser.error('--bound must be positive')
    if any(not np.isfinite(z) or z<=0 for z in args.expand_bounds) or args.expand_bounds != sorted(set(args.expand_bounds)):
        parser.error('--expand-bounds must be positive finite strictly increasing values')
    if not np.isfinite(args.expand_if_above) or args.expand_if_above < args.refine_gate:
        parser.error('--expand-if-above must be >= --refine-gate')
    if not np.isfinite(args.accept_chi2) or args.accept_chi2 <= 0:
        parser.error('--accept-chi2 must be finite and positive')
    if args.stage in ('fit','cross','all') and not args.original_c5:
        from Numerical.orchestration.FastLambda5e9C5 import install
        install()
        print('Fast numeric final-C5 evaluation active (5e-9 runner only).', flush=True)
    if args.stage in ('prepare','all'):
        prepare()
    if args.stage in ('fit','all'):
        fit(args)
    if args.stage in ('cross','all'):
        cross(args)
    if args.stage in ('analyse','all'):
        analyse(args)

if __name__ == '__main__':
    main()
