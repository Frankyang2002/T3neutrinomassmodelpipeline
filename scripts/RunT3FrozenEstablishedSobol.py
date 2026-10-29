"""Frozen-UV 16x2 conditional T3 fit using the established Sobol multistart stages.

Does not alter existing fits, RGEs, matching, or screening thresholds.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import re
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from scripts import RunT3GeneralScalar48Campaign as bulk
from Numerical.fitting.SobolBenchmarkSearch import (
    BenchmarkParameter, generate_sobol_parameter_points, refine_sobol_candidates,
)

CASES = ('small_uv_frozen_5e-9', 'large_uv_frozen_0p5')
FORMAT = 't3-frozen-established-sobol-escalating-v2'
ACCEPT = 5.86


def atomic_json(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    temp.replace(path)


def source_fit(source, model, case):
    path = bulk.OUT / source / case / bulk.safe(model) / 'fit.json'
    if not path.is_file():
        return None
    obj = bulk.read(path)
    if obj.get('model') != model or obj.get('case') != case:
        raise ValueError(f'Previous fit identity mismatch: {path}')
    return obj


def candidate_from_row(row, origin):
    if row is None or row.get('chi2') is None or not isinstance(row.get('parameters'), dict):
        return None
    return {'chi2': float(row['chi2']), 'parameters': dict(row['parameters']),
            'whitened_residual': row.get('whitened_residual'),
            'prediction': row.get('prediction'), 'source': origin}


def choose_candidates(rows, keep):
    valid = [r for r in rows if r is not None and math.isfinite(float(r['chi2']))]
    valid.sort(key=lambda r: float(r['chi2']))
    return valid[:keep]


def run_one(args, model, case, registry, output):
    entry, new, reused = bulk.inventory(registry, model)
    if reused:
        prior = bulk.existing_fit(case, model)
        return {'model':model,'case':case,'status':'ExistingSpecialFit' if prior else 'SpecialFitUnavailable',
                'chi2': (prior or {}).get('best_mass_safe',{}).get('chi2'),
                'note':'B/C mixed-quadratic model uses its established special evaluator; not bulk evaluator'}
    prior = source_fit(args.source_campaign, model, case)
    prior_best = (prior or {}).get('best_mass_safe')
    prior_status = (prior or {}).get('status')
    if prior_status in ('CandidateConstrainedAccepted', 'ExistingConditionalAccepted') and not args.refit_accepted:
        return {'model':model,'case':case,'status':'PreservedAccepted',
                'chi2': (prior_best or {}).get('chi2'), 'source_campaign':args.source_campaign}
    if output.joinpath('result.json').exists():
        return bulk.read(output/'result.json')

    history = bulk.historical(case, model, entry)
    if history['uv_on']:
        raise ValueError('Frozen-only search received UV-running case')
    beta_path = bulk.OUT / args.source_campaign / 'betas' / f'{bulk.safe(model)}.json'
    payload = bulk.read(beta_path)
    kernel, match_path, kind = bulk.matching(model)
    from Numerical.orchestration import RunSmallLambdaRefit as refit
    from Numerical.fitting.OscillationFit import OscillationFitTarget
    target = OscillationFitTarget.from_json(refit.TARGET)
    evaluator = bulk.Evaluator(case,model,entry,new,history,payload,kernel,target,args.extra_real)
    lo, hi = bulk.b21.search_box(history['seed'],local_fraction=1.5,
                                  floor_fraction=.05,abs_yukawa_max=math.sqrt(4*math.pi))
    definitions = tuple(BenchmarkParameter(name=n,binding='unused',low=float(a),high=float(b))
                        for n,a,b in zip(bulk.NAMES,lo,hi))
    def physical(values):
        return np.asarray([values[n] for n in bulk.NAMES],dtype=float)
    rows = []
    if candidate_from_row(prior_best,'previous_best'):
        rows.append(candidate_from_row(prior_best,'previous_best'))
    # Import the complete saved top-candidate pool from an earlier run.
    if args.resume_from_run:
        old_dir = bulk.OUT / args.resume_from_run / case / bulk.safe(model)
        if not old_dir.is_dir():
            raise FileNotFoundError(f'Previous run case missing: {old_dir}')
        for old_cp in sorted(old_dir.glob('round_*.json')):
            saved = bulk.read(old_cp)
            if saved.get('model') != model or saved.get('case') != case:
                raise ValueError(f'Checkpoint identity mismatch: {old_cp}')
            rows.extend(saved.get('retained', []))
        old_result = old_dir / 'result.json'
        if old_result.exists():
            saved = bulk.read(old_result)
            if saved.get('model') != model or saved.get('case') != case:
                raise ValueError(f'Previous result identity mismatch: {old_result}')
            if saved.get('best'):
                rows.append(saved['best'])
    # Save previous *fitted* point explicitly as a high-value multistart seed.
    initial = evaluator.evaluate(history['seed'],'historical_seed')
    rows.append(candidate_from_row(initial,'historical_seed'))
    rows = choose_candidates(rows, args.keep_best)
    if args.refine_only:
        if not rows:
            raise RuntimeError('No prior mass-safe Yukawa seed for --refine-only')
        def fit_adapter(values):
            trial = evaluator.evaluate(physical(values),'full_18_parameter_refinement')
            if trial.get('chi2') is None:
                raise ValueError(trial.get('reason','point rejected by physical screen'))
            return SimpleNamespace(chi2=trial['chi2'], ordering='NO',
                observable_names=tuple(target.observable_names),
                prediction=np.asarray([trial['prediction'][n] for n in target.observable_names]),
                whitened_residual=np.asarray(trial['whitened_residual'],dtype=float))
        checkpoint = output / 'full18_refinement.json'
        if checkpoint.is_file():
            saved = bulk.read(checkpoint)
            rows.extend(saved.get('retained', []))
        elif rows[0]['chi2'] > args.target:
            try:
                improved, details = refine_sobol_candidates(
                    evaluator=fit_adapter, observable_names=target.observable_names,
                    definitions=definitions, candidates=rows, target_chi2=args.target,
                    seed_count=args.seed_count, active_counts=(18,),
                    max_nfev=(args.nfev_18,), sensitivity_step=args.sensitivity_step)
                rows.append(improved)
                rows = choose_candidates(rows,args.keep_best)
                atomic_json(checkpoint, {'format':FORMAT,'model':model,'case':case,
                    'retained':rows,'refinement':details,'nfev_18':args.nfev_18,
                    'sensitivity_step':args.sensitivity_step})
            except Exception as exc:
                atomic_json(checkpoint, {'format':FORMAT,'model':model,'case':case,
                    'retained':rows,'refinement':{'status':'NumericalFailure','reason':str(exc)[:500]}})
                raise
        print(f'{model} {case} full18 best_chi2={rows[0]["chi2"]}',flush=True)
    seed_base = int.from_bytes(hashlib.sha256(f'{args.seed}:{model}:{case}'.encode()).digest()[:4],'big')
    for round_number in (() if args.refine_only else range(1,args.rounds+1)):
        round_points = points_for_round(args.sobol_points, round_number, args.max_sobol_points)
        checkpoint = output / f'round_{round_number:03d}.json'
        if checkpoint.exists():
            saved = bulk.read(checkpoint)
            rows.extend(saved.get('retained',[]))
            rows=choose_candidates(rows,args.keep_best)
            if rows and rows[0]['chi2'] <= args.target:
                print(f'{model} {case}: accepted from checkpoint round={round_number}',flush=True)
                break
            continue
        points = generate_sobol_parameter_points(definitions,round_points,
                                                  seed=(seed_base+round_number) % (2**32))
        round_rows = []
        for i, point in enumerate(points,1):
            result = evaluator.evaluate(physical(point),f'round_{round_number}_sobol_{i}')
            candidate = candidate_from_row(result,f'sobol_round_{round_number}_{i}')
            if candidate is not None:
                round_rows.append(candidate)
        seeds = choose_candidates(rows+round_rows, max(args.keep_best,args.seed_count))
        refinement = None
        if seeds and float(seeds[0]['chi2'])>args.target and args.refine:
            # Adapt the existing evaluator to the interface of the *original*
            # Numerical.fitting.SobolBenchmarkSearch.refine_sobol_candidates.
            def fit_adapter(values):
                trial = evaluator.evaluate(physical(values),'established_local_refinement')
                if trial.get('chi2') is None:
                    raise ValueError(trial.get('reason','point rejected by physical screen'))
                return SimpleNamespace(chi2=trial['chi2'], ordering='NO',
                    observable_names=tuple(target.observable_names),
                    prediction=np.asarray([trial['prediction'][n] for n in target.observable_names]),
                    whitened_residual=np.asarray(trial['whitened_residual'],dtype=float))
            try:
                improved,refinement = refine_sobol_candidates(
                    evaluator=fit_adapter,observable_names=target.observable_names,
                    definitions=definitions,candidates=seeds,target_chi2=args.target,
                    seed_count=args.seed_count,active_counts=(5,8,12,18),
                    max_nfev=(args.nfev_5,args.nfev_8,args.nfev_12,args.nfev_18),
                    sensitivity_step=args.sensitivity_step)
                rows.append(improved)
            except Exception as exc:
                refinement={'status':'NumericalFailure','reason':str(exc)[:500]}
        rows.extend(round_rows)
        rows=choose_candidates(rows,args.keep_best)
        best=rows[0] if rows else None
        atomic_json(checkpoint,{'format':FORMAT,'model':model,'case':case,'round':round_number,
                    'retained':rows,'refinement':refinement,'sobol_points':round_points,
                    'evaluations':evaluator.count,'mass_safe':evaluator.mass_safe,
                    'numerical_failures':evaluator.failures})
        print(f'{model} {case} round={round_number}/{args.rounds} sobol_points={round_points} best_chi2={None if best is None else best["chi2"]}',flush=True)
        if best is not None and best['chi2']<=args.target:
            break
    best = choose_candidates(rows,1)
    best=best[0] if best else None
    status = ('CandidateConstrainedAccepted' if best and best['chi2']<=args.target else
              'CandidateConstrainedPoorFit' if best else 'NoMassSafeNeutrinoFit')
    result={'format':FORMAT,'model':model,'case':case,'status':status,
            'chi2':best['chi2'] if best else None,'best':best,
            'historical_lambdaT3_fixed':history['lambda_fixed'],
            'matching_kind':kind,'matching_sha256':bulk.sha(match_path),
            'beta_sha256':bulk.sha(beta_path),
            'original_fitter':'Numerical.fitting.SobolBenchmarkSearch.refine_sobol_candidates',
            'candidate_general_quartic_only':True,'sequential_matching_certified':False,
            'search_budget_exhausted':bool(status!='CandidateConstrainedAccepted'),
            'evaluations':evaluator.count,'mass_safe':evaluator.mass_safe,
            'numerical_failures':evaluator.failures}
    if best:
        cfg=bulk.update_yukawas(history,physical(best['parameters']))
        cfg['conditional_general_scalar_fit']={'model':model,'case':case,
             'chi2':best['chi2'],'sequential_matching_certified':False,
             'source':'established_sobol_multistart'}
        atomic_json(output/'fitted_config.json',cfg)
    atomic_json(output/'result.json',result)
    return result


def points_for_round(initial, round_number, maximum):
    if initial < 2 or initial & (initial-1) or maximum < initial or maximum & (maximum-1):
        raise ValueError('Sobol counts must be powers of two and maximum >= initial')
    return min(maximum, initial * (2 ** (round_number - 1)))


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source-campaign',default='overnight_six_v1_bulk48')
    p.add_argument('--run-id',required=True)
    p.add_argument('--sobol-points',type=int,default=128)
    p.add_argument('--rounds',type=int,default=8,help='Maximum rounds per model/scenario before halting')
    p.add_argument('--max-sobol-points',type=int,default=4096,help='Cap on Sobol points in each round')
    p.add_argument('--continue-on-failure',action='store_true',help='Move to next case despite failing target after all rounds')
    p.add_argument('--seed',type=int,default=20261030)
    p.add_argument('--seed-count',type=int,default=5)
    p.add_argument('--keep-best',type=int,default=12)
    p.add_argument('--nfev-5',type=int,default=30)
    p.add_argument('--nfev-8',type=int,default=40)
    p.add_argument('--nfev-12',type=int,default=60)
    p.add_argument('--nfev-18',type=int,default=400,help='Evaluation budget for final full-18 stage')
    p.add_argument('--sensitivity-step',type=float,default=0.001,help='Unit-coordinate sensitivity step')
    p.add_argument('--resume-from-run',default='',help='Import top seeds from previous run, preserving original outputs')
    p.add_argument('--refine-only',action='store_true',help='Skip Sobol and run full 18-parameter refinement from imported seeds')
    p.add_argument('--extra-real',type=float,default=-0.01)
    p.add_argument('--target',type=float,default=ACCEPT)
    p.add_argument('--models',default='',help='Optional comma-separated subset')
    p.add_argument('--refit-accepted',action='store_true')
    p.add_argument('--no-refine',dest='refine',action='store_false')
    args=p.parse_args(argv)
    if not re.fullmatch('[A-Za-z0-9][A-Za-z0-9_-]{0,63}',args.run_id):p.error('invalid run-id')
    if args.sobol_points<2 or args.sobol_points & (args.sobol_points-1):p.error('Sobol points must be a power of two')
    if args.rounds<1 or args.seed_count<1 or args.keep_best<args.seed_count:p.error('invalid search budget')
    if not (1<=args.nfev_18<=10000 and 0<args.sensitivity_step<=0.25):p.error('invalid full-18 refinement settings')
    if args.refine_only and not args.resume_from_run:p.error('--refine-only requires --resume-from-run')
    if args.resume_from_run and not re.fullmatch('[A-Za-z0-9][A-Za-z0-9_-]{0,63}',args.resume_from_run):p.error('invalid resume-from-run')
    try: points_for_round(args.sobol_points,1,args.max_sobol_points)
    except ValueError as exc:p.error(str(exc))
    registry=bulk.read(bulk.REGISTRY)
    if len(registry['models'])!=16:raise RuntimeError('Expected 16 T3 models')
    selected=set(args.models.split(',')) if args.models else set(registry['models'])
    if not selected.issubset(registry['models']):p.error('unknown model')
    base=bulk.OUT/args.run_id
    base.mkdir(parents=True,exist_ok=True)
    records=[]
    halted=False
    for model in registry['models']:
        if model not in selected:continue
        for case in CASES:
            output=base/case/bulk.safe(model)
            output.mkdir(parents=True,exist_ok=True)
            try:
                rec=run_one(args,model,case,registry,output)
            except Exception as exc:
                rec={'model':model,'case':case,'status':'FailedOrUnavailable','chi2':None,
                     'reason':str(exc)[:700]}
                atomic_json(output/'error.json',rec)
            records.append(rec)
            if rec.get('status') in ('CandidateConstrainedPoorFit','NoMassSafeNeutrinoFit','FailedOrUnavailable','SpecialFitUnavailable') and not args.continue_on_failure:
                halted=True
                print(f'HALTED at {model} {case}: {rec["status"]}; did not proceed to next case. See saved checkpoints. ',flush=True)
            atomic_json(base/'summary.json',{'format':FORMAT,'rows':records})
            print(f'CASE {model} {case} {rec["status"]} chi2={rec.get("chi2")}',flush=True)
            if halted: break
        if halted: break
    if records:
        with (base/'summary.csv').open('w',newline='',encoding='utf-8') as f:
            w=csv.DictWriter(f,fieldnames=('model','case','status','chi2','reason'),extrasaction='ignore')
            w.writeheader();w.writerows(records)
    print('Saved',base/'summary.csv',flush=True)
    return 2 if halted else 0

if __name__=='__main__':
    raise SystemExit(main())
