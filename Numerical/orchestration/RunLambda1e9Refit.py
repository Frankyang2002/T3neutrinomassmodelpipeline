"""Local 1e-9 T3 refit, donor coverage and Yukawa clustering.

Usage: python -m Numerical.orchestration.RunLambda1e9Refit {prepare,fit,cross,analyse,all}
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

OUT = base.ROOT / 'output' / 'studies' / 'lambda_1e-9_refit'
OLD = base.STUDY
LAMBDA = {'normal': 1e-9, 'E': 5e-10}
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
        cfg['lambda_1e9_study'] = {
            'status': 'Seed only; requires independent numerical optimisation',
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
    """Adaptive Sobol boxes followed by ORIGINAL sensitivity-ranked refinement.

    Existing local Sobol checkpoints are migrated in memory without losing
    their best candidates. Wider ranges are diagnostics, not perturbativity claims.
    """
    from Numerical.fitting.SobolBenchmarkSearch import (
        BenchmarkParameter, generate_sobol_parameter_points,
        refine_sobol_candidates, _build_evaluator,
    )
    import warnings

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

        cfg = base.load(name_path('prepared', model))
        if not np.isclose(float(cfg['base_state']['ordinary']['lambdaT3']['real']), target_lambda(row), rtol=1e-12):
            raise ValueError(f'{model}: unexpected lambda')
        seed = get_vector(cfg)
        local_halfwidth = np.maximum(args.sobol_fraction * np.abs(seed), args.sobol_floor)
        stages = [('local', np.maximum(-args.bound, seed-local_halfwidth),
                   np.minimum(args.bound, seed+local_halfwidth))]
        for width in args.expand_bounds:
            stages.append((f'global_+-{width:g}', np.full(18, -width), np.full(18, width)))
        if any(np.any(lo >= hi) for _,lo,hi in stages):
            raise ValueError(f'{model}: invalid local bounds; increase --bound')

        names = [f'{key}_{p}{q}' for key in ('y1','y2') for p in range(3) for q in range(3)]
        checkpoint_path = name_path('sobol_search',model)
        prepared_sha = hashlib.sha256(name_path('prepared', model).read_bytes()).hexdigest()
        artifact_sha = [hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in artifacts]
        target_sha = hashlib.sha256(base.TARGET.read_bytes()).hexdigest()
        stable_sig = {'model':model, 'lambdaT3':target_lambda(row),
                      'artifact_sha256':artifact_sha,'target_sha256':target_sha,
                      'prepared_sha256':prepared_sha,'sobol_seed':args.sobol_seed+i}
        old = base.load(checkpoint_path) if args.resume and checkpoint_path.is_file() else None
        if old:
            sig = old['signature']
            for key,val in stable_sig.items():
                if sig.get(key) != val:
                    raise RuntimeError(f'{checkpoint_path}: changed {key}; refusing incompatible resume')
            if 'stages' in old:
                progress = old['stages']
            else:
                # Migration from local-only Sobol checkpoint, without loss of good candidates.
                original_bounds = sig.get('bounds')
                local_bounds = [stages[0][1].tolist(), stages[0][2].tolist()]
                if not original_bounds or not np.allclose(original_bounds,local_bounds,rtol=0,atol=1e-12):
                    raise RuntimeError(f'{checkpoint_path}: previous local bounds differ; cannot migrate automatically')
                progress = [{'label':'local','evaluated':int(old.get('sobol_evaluated',0)),
                             'failed':int(old.get('sobol_failed',0)), 'completed':False}]
            candidates = list(old.get('candidates',[]))
        else:
            progress=[]
            candidates=[]
        progress_by_label = {p['label']:p for p in progress}
        for label,_,_ in stages:
            progress_by_label.setdefault(label,{'label':label,'evaluated':0,'failed':0,'completed':False})
        progress=[progress_by_label[label] for label,_,_ in stages]
        # Preserve best saved candidate even when it is outside a narrower box.
        if not candidates:
            # Evaluate seed with the SAME numerical physics, at current bounds.
            pass

        def score():
            return min((float(z['chi2']) for z in candidates), default=np.inf)

        def save_checkpoint():
            best=candidates[0]['chi2'] if candidates else None
            base.save(checkpoint_path, {'version':2,'signature':stable_sig,
                     'stages':progress,'candidates':candidates,
                     'sobol_evaluated':sum(p['evaluated'] for p in progress),
                     'sobol_failed':sum(p['failed'] for p in progress),
                     'best_chi2':best,'refine_gate':args.refine_gate})

        chosen_defs=chosen_ev=None
        def make_defs(lo,hi):
            return [BenchmarkParameter(names[k], f'ordinary.{key}.real.{p}.{q}',float(lo[k]),float(hi[k]))
                    for k,(key,p,q) in enumerate((key,p,q) for key in ('y1','y2') for p in range(3) for q in range(3))]
        def make_evaluator(defs):
            return _build_evaluator(payload=cfg, parameters=defs,
                  uv_rgbeta_path=artifacts[0],intermediate_rgbeta_path=artifacts[1],
                  final_weinberg_path=artifacts[2],target_path=base.TARGET)
        def evaluate(ev,values,source,index=None):
            res=ev(values)
            chi=float(res.chi2)
            if not np.isfinite(chi): raise ValueError('nonfinite chi2')
            return {'index':index,'source':source,'chi2':chi,'parameters':dict(values),
                    'ordering':str(res.ordering),
                    'prediction':dict(zip(res.observable_names,map(float,res.prediction))),
                    'whitened_residual':np.asarray(res.whitened_residual,dtype=float).tolist()}
        def retain(candidate):
            candidates.append(candidate)
            candidates.sort(key=lambda z:float(z['chi2']))
            del candidates[max(args.keep_best,args.seed_count):]

        total_budget=args.sobol_max
        stage_cap=args.sobol_points
        for index,(label,lo,hi) in enumerate(stages):
            state=progress[index]
            if score() <= args.refine_gate:
                break
            if state['completed']:
                continue
            completed_before=sum(p['evaluated'] for p in progress)
            if completed_before >= total_budget:
                break
            defs=make_defs(lo,hi)
            ev=make_evaluator(defs)
            if not candidates:
                try:
                    retain(evaluate(ev,dict(zip(names,seed.tolist())),'scaled_seed'))
                except Exception as exc:
                    print(f'[{i}/16] seed evaluation failed: {exc}',flush=True)
            if args.resume and fitted_path.is_file() and not any(z.get('source')=='previous_incomplete_fit' for z in candidates):
                oldvec=get_vector(base.load(fitted_path))
                if np.all(oldvec >= lo) and np.all(oldvec <= hi):
                    try: retain(evaluate(ev,dict(zip(names,oldvec.tolist())),'previous_incomplete_fit'))
                    except Exception as exc: print(f'[{i}/16] saved fit could not be evaluated: {exc}',flush=True)
            if score() <= args.refine_gate:
                chosen_defs,chosen_ev=defs,ev
                save_checkpoint()
                break
            # Each stage gets one 64-point Sobol block. The final stage can
            # use the remaining overall cap; every stage is deterministic.
            cap = stage_cap if index < len(stages)-1 else total_budget
            maximum = min(cap, state['evaluated']+max(0,total_budget-completed_before))
            if maximum < state['evaluated']:
                maximum=state['evaluated']
            # Sobol generator requires power-of-two count; produce a sufficient
            # deterministic prefix and use only the remaining indices.
            count=1 << max(1,(maximum-1).bit_length())
            if maximum > 0 and state['evaluated'] < maximum:
                pts=generate_sobol_parameter_points(defs,count,seed=args.sobol_seed+i+1009*index)
                for j in range(state['evaluated'],maximum):
                    try:
                        retain(evaluate(ev,pts[j],f'sobol_{label}',j+1))
                    except Exception as exc:
                        state['failed']+=1
                        if state['failed']<=3:
                            print(f'[{i}/16] {label} point {j+1} failed: {type(exc).__name__}: {exc}',flush=True)
                    state['evaluated']=j+1
                    if state['evaluated']%args.sobol_points==0 or score()<=args.refine_gate:
                        save_checkpoint()
                        print(f'[{i}/16] {model} {label} Sobol {state["evaluated"]}: best chi2={score():.7g}',flush=True)
                    if score() <= args.refine_gate:
                        break
            if score() > args.expand_if_above and index<len(stages)-1:
                state['completed']=True
                save_checkpoint()
                print(f'[{i}/16] {model}: expand after {state["evaluated"]} {label} points; '
                      f'best chi2={score():.7g} > {args.expand_if_above:g}',flush=True)
                continue
            chosen_defs,chosen_ev=defs,ev
            if index<len(stages)-1:
                # Intermediate region was not horrible; stay in it and stop
                # expanding, but record failure if gate was not achieved.
                state['completed']=True
                save_checkpoint()
                break
            state['completed']=True
            save_checkpoint()

        if not candidates:
            record={'model':model,'class':base.cls(row),'status':'NoValidSeed',
                    'chi2_fit':None,'sobol_evaluated':sum(x['evaluated'] for x in progress)}
            base.save(result_path,record);records.append(record)
            base.save(OUT/'fit_summary.json',records)
            continue
        if score()>args.refine_gate:
            record={'model':model,'class':base.cls(row),'status':'NeedsBetterSobolSeed',
                    'best_seed_chi2':score(),'chi2_fit':None,
                    'sobol_evaluated':sum(x['evaluated'] for x in progress),
                    'sobol_failed':sum(x['failed'] for x in progress),
                    'refine_gate':args.refine_gate,'sobol_max':args.sobol_max,
                    'region_progress':progress}
            base.save(result_path,record);records.append(record)
            base.save(OUT/'fit_summary.json',records)
            print(f'[{i}/16] {model}: best {score():.7g} > gate {args.refine_gate:g}; no refinement',flush=True)
            continue
        # Use the widest reached search region for local optimisation, but only
        # candidates INSIDE that region can be transformed to unit coordinates.
        if chosen_defs is None:
            reached=max((j for j,x in enumerate(progress) if x['evaluated']>0),default=0)
            chosen_defs=make_defs(stages[reached][1],stages[reached][2]);chosen_ev=make_evaluator(chosen_defs)
        def in_bounds(item):
            return all(p.low<=float(item['parameters'][p.name])<=p.high for p in chosen_defs)
        usable=[c for c in candidates if in_bounds(c)]
        if not usable:
            raise RuntimeError(f'{model}: best candidates outside active refinement bounds')
        best=min(usable,key=lambda z:float(z['chi2']))
        if best['chi2']<=args.accept_chi2:
            history={'skipped':True,'reason':'seed already meets target','stages':[]}
        else:
            best,history=refine_sobol_candidates(evaluator=chosen_ev,observable_names=FIELDS,
                definitions=chosen_defs,candidates=sorted(usable,key=lambda x:x['chi2'])[:args.keep_best],
                target_chi2=args.accept_chi2,seed_count=args.seed_count,
                active_counts=(5,8,12),max_nfev=(30,40,60),sensitivity_step=0.02)
        final=chosen_ev(best['parameters'])
        fitted=copy.deepcopy(cfg)
        put_vector(fitted,np.asarray([best['parameters'][name] for name in names],dtype=float))
        fitted['lambda_1e9_study'].update(status='Adaptive Sobol and original staged refinement',
            chi2=float(final.chi2),accept_chi2=args.accept_chi2)
        record={'model':model,'class':base.cls(row),'status':'Success' if float(final.chi2)<=args.accept_chi2 else 'PoorFit',
            'lambdaT3':target_lambda(row),'chi2_fit':float(final.chi2),
            'fit_method':'adaptive_sobol_original_multistart_5_8_12',
            'passes_screening':bool(float(final.chi2)<=args.accept_chi2),
            'sobol_evaluated':sum(x['evaluated'] for x in progress),
            'region_progress':progress,'refinement':history,
            'prediction':dict(zip(final.observable_names,map(float,final.prediction)))}
        base.save(fitted_path,fitted);base.save(result_path,record)
        records.append(record);base.save(OUT/'fit_summary.json',records)
        print(f'[{i}/16] {model}: chi2={final.chi2:.7g}, pass={record["passes_screening"]}',flush=True)
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
              'notes': ['New individually fitted Yukawas at lambda=1e-9 (class E 5e-10)',
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
    parser.add_argument('--bound', type=float, default=2.0, help='Absolute entry bound for Yukawa optimisation')
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
    parser.add_argument('--expand-bounds', type=float, nargs='+', default=[2.,10.,100.], help='Progressively wider global +- Yukawa search boxes')
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
        from Numerical.orchestration.FastLambda1e9C5 import install
        install()
        print('Fast numeric final-C5 evaluation active (10^-9 runner only).', flush=True)
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
