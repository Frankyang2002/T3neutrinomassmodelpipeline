#!/usr/bin/env python3
"""One optimisation per T3 model, stopping on verified chi2 < 0.01.

The fit uses 18 REAL UV Yukawa entries y1/y2 and the existing five NuFIT
uncorrelated Gaussian residuals. It does NOT fit general complex Yukawas.

For each supported model, either stored A(-4) or D(-2) UV Yukawa benchmark
can be selected as the INITIAL CONDITION. 'best-of-two' evaluates both
benchmarks then starts ONE least-squares optimisation from the better one.
This is 16 model fits (not 32 model/seed fits).

A(-4) and D(-2) use the latest B11/B13/B15/B16 *conditional* source sector
with fresh UV/EFT1 RG at every trial and the existing SMEFT -> PMNS steps.
Other 14 lack a corresponding full-tensor per-trial source backend and are
explicitly BLOCKED; no old rough/analytic-triangle fallback is used.

Use --plan before --run. All outputs are NOT physical C5 certification.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import sys
import traceback

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0,str(ROOT))

from Numerical.parallel_t3.complete.all_model_conditional_optimizer import (
    FORMAT, BENCHMARK, B13, B16, SEEDS, SUPPORTED, DEFAULT_TARGET_CHI2,
    LatestTwoModelEvaluator, SourceBlocked, all_models, jobs, choose_start,
    optimize_job, vector_from_config,
)


def save(path:Path,data:dict)->None:
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_name(path.name+'.partial')
    tmp.write_text(json.dumps(data,indent=2,sort_keys=True,allow_nan=False,
                              default=str)+'\n',encoding='utf-8')
    tmp.replace(path)


def selected_models(raw:list[str])->tuple[str,...]:
    inventory=all_models()
    if raw==['all']:
        return inventory
    if not raw or any(x not in inventory for x in raw) or len(set(raw))!=len(raw):
        raise SourceBlocked('Invalid or repeated --models; use all or explicit model names')
    return tuple(raw)


def dependencies(root:Path,model:str)->list[Path]:
    from Numerical.parallel_t3.complete.two_model_conditional_uvfit_batch12 import (
        B1,B2,B11,LAURENT,SOFT_PATHS)
    from Numerical.parallel_t3.complete.two_model_joint_scalar_batch9 import (
        quartic_source_path)
    from Numerical.parallel_t3.complete.two_model_master_wilson_batch13 import (
        D_EXTRA,REGISTRY)
    required=[root/BENCHMARK,root/B13/model/'result.json',
              root/B16/model/'result.json',root/B11/model/'result.json',
              root/B1/model/'manifest.json',root/B1/model/'trajectories.npz',
              root/B1/model/'fermion_matching_input.json',
              root/B1/model/'scalar_matching_input.json',
              root/B2/model/'one_loop_operator_ledger.json',
              root/LAURENT/model/'manifest.json',root/SOFT_PATHS[model],
              quartic_source_path(root,model)]
    if model=='D_-2': required.extend((root/D_EXTRA,root/REGISTRY))
    return required


def resume_vector(from_root:Path,model:str,engine:LatestTwoModelEvaluator,
                  bound:float)->np.ndarray|None:
    path=from_root/model/'result.json'
    if not path.is_file():
        return None
    doc=json.loads(path.read_text(encoding='utf-8-sig'))
    if doc.get('model')!=model or not isinstance(doc.get('results'),dict):
        raise SourceBlocked(f'{model}: cannot resume a non-successful model result')
    if doc.get('source_files_sha256')!=engine.source_hashes:
        raise SourceBlocked(f'{model}: source inputs changed since previous optimisation')
    if float(doc.get('signed_factor_assumed',float('nan')))!=engine.signed_factor:
        raise SourceBlocked(f'{model}: signed normalization changed since last fit')
    x=np.asarray(doc['results'].get('best_18_real_UV_yukawas'),float)
    if x.shape!=(18,) or not np.isfinite(x).all() or np.any(np.abs(x)>=bound):
        raise SourceBlocked(f'{model}: invalid saved Yukawa state or different bounds')
    return x


def main(argv:list[str]|None=None)->int:
    p=argparse.ArgumentParser(description=__doc__)
    action=p.add_mutually_exclusive_group(required=True)
    action.add_argument('--plan',action='store_true')
    action.add_argument('--preflight',action='store_true')
    action.add_argument('--run',action='store_true')
    p.add_argument('--models',nargs='+',default=['all'],help='all or explicit model labels')
    p.add_argument('--initial-seed',choices=('best-of-two',)+SEEDS,
                   default='best-of-two',help='Choose one starting point per model')
    p.add_argument('--allow-partial',action='store_true',
                   help='Run models with a backend; explicitly report unsupported others BLOCKED')
    p.add_argument('--target-chi2',type=float,default=DEFAULT_TARGET_CHI2,
                   help='Stop when verified five-observable Gaussian chi2 is strictly below this')
    p.add_argument('--max-nfev',type=int,default=100,
                   help='Max SciPy residual calls per MODEL fit (Jacobians cause more pipeline runs)')
    p.add_argument('--yukawa-bound',type=float,default=0.005)
    p.add_argument('--regularization',type=float,default=0.,
                   help='Optional squared distance to the chosen UV seed; chi2 target excludes this')
    p.add_argument('--signed-factor',type=float,default=1.,
                   help='Explicit UNVERIFIED GammaCC/Weyl normalization')
    p.add_argument('--resume-from',type=Path,default=None,
                   help='Previous v2 output; start each model from its best saved Yukawas')
    p.add_argument('--output',type=Path,
                   default=Path('output/parallel_t3/one_fit_per_model_target_001_v1'))
    p.add_argument('--debug',action='store_true')
    args=p.parse_args(argv)
    if not np.isfinite(args.target_chi2) or not 0<args.target_chi2<1e10:
        p.error('--target-chi2 must be positive and finite')
    if not np.isfinite(args.yukawa_bound) or not 0<args.yukawa_bound<.2:
        p.error('Invalid --yukawa-bound')
    if args.max_nfev<1 or args.max_nfev>10000:
        p.error('--max-nfev must be in [1,10000]')
    if not np.isfinite(args.signed_factor) or not 0<abs(args.signed_factor)<1000:
        p.error('--signed-factor must be nonzero finite')
    if not np.isfinite(args.regularization) or args.regularization<0:
        p.error('--regularization must be nonnegative finite')
    try:
        models=selected_models(args.models)
        queue=jobs(models)
    except (ValueError,SourceBlocked) as exc:
        p.error(str(exc))
    if args.plan:
        print(json.dumps({'format':FORMAT,'models':queue,'total_model_fits':len(queue),
            'supported_model_fits':sum(x['supported'] for x in queue),
            'starting_point_policy':args.initial_seed,
            'target_chi2_strict':args.target_chi2,
            'source_matching_only_conditional':True},indent=2))
        return 0
    missing=[]
    for model in models:
        if model not in SUPPORTED:
            continue
        for path in dependencies(ROOT,model):
            if not path.is_file(): missing.append(str(path.relative_to(ROOT)))
            if args.preflight:
                print(('FOUND ' if path.is_file() else 'MISSING ')+str(path),flush=True)
    unsupported=[x for x in queue if not x['supported']]
    if args.preflight:
        if unsupported:
            print(f'{len(unsupported)} model(s) currently lack matching backend.',flush=True)
            for entry in unsupported: print('BLOCKED_BACKEND',entry['model'],flush=True)
        return 2 if missing or (unsupported and not args.allow_partial) else 0
    if unsupported and not args.allow_partial:
        p.error('14 of 16 full-tensor matching backends are not implemented. '
                'Use --models A_-4 D_-2 or explicitly acknowledge with --allow-partial.')
    if missing:
        p.error('Missing required source files: '+', '.join(missing[:8]))
    output=(ROOT/args.output).resolve()
    allowed=(ROOT/'output/parallel_t3').resolve()
    if allowed not in output.parents:
        p.error('Output must be a NEW child of output/parallel_t3')
    if output.exists():
        p.error('Refusing to overwrite existing output directory: '+str(output))
    previous=(ROOT/args.resume_from).resolve() if args.resume_from is not None else None
    if previous is not None and not previous.is_dir():
        p.error('--resume-from directory does not exist')
    cfg=json.loads((ROOT/BENCHMARK).read_text(encoding='utf-8-sig'))
    rows=[]
    for entry in queue:
        model=entry['model']
        if not entry['supported']:
            row={'model':model,'status':'BLOCKED_NoCurrentFullGeneralMatchingBackend',
                 'reason':entry['reason'],'target_met':False,
                 'physical_C5_certified':False}
            print('BLOCKED',model,row['reason'],flush=True)
        else:
            try:
                print('INITIALIZE',model,'B11+B13+B15+B16',flush=True)
                engine=LatestTwoModelEvaluator(ROOT,model,signed_factor=args.signed_factor)
                model_original=vector_from_config(cfg,model)
                baseline_rel=engine.verify_fresh_baseline(model_original)
                original_seed,initial,seed_chi2=choose_start(
                    engine.evaluate,cfg,policy=args.initial_seed)
                start_source=original_seed
                if previous is not None:
                    old=resume_vector(previous,model,engine,args.yukawa_bound)
                    if old is not None:
                        initial=old
                        start_source='resumed_previous_best'
                print('FIT',model,'start',start_source,'target chi2 <',args.target_chi2,
                      'max_nfev',args.max_nfev,flush=True)
                fit=optimize_job(engine.evaluate,initial,max_nfev=args.max_nfev,
                                 bound=args.yukawa_bound,
                                 regularization=args.regularization,
                                 target_chi2=args.target_chi2)
                row={'model':model,'status':fit['status'],
                     'target_met':fit['target_met_after_fresh_benchmark_validation'],
                     'initial_seed':start_source,
                     'original_seed_baseline_chi2':seed_chi2,
                     'reference_B16_relative_C5':engine.baseline_relative,
                     'fresh_RG_baseline_relative_C5':baseline_rel,
                     'source_files_sha256':engine.source_hashes,
                     'signed_factor_assumed':args.signed_factor,
                     'results':fit,'physical_C5_certified':False}
                print('DONE',model,'best chi2',fit['best_chi2'],
                      'target reached',row['target_met'],flush=True)
            except Exception as exc:
                row={'model':model,'status':'BLOCKED_SourceOrNumericalEvaluation',
                     'reason':f'{type(exc).__name__}: {exc}',
                     'target_met':False,'physical_C5_certified':False}
                if args.debug: row['traceback']=traceback.format_exc(limit=20)
                print('BLOCKED',model,row['reason'],file=sys.stderr,flush=True)
        save(output/model/'result.json',row)
        rows.append({'model':model,'status':row['status'],
                     'target_met':row['target_met'],
                     'best_chi2':(row.get('results') or {}).get('best_chi2'),
                     'reason':row.get('reason')})
        summary={'format':FORMAT,
            'created_utc':datetime.now(timezone.utc).isoformat(),
            'source_revision_checked':'61422855d36c9b2a964dbcfc158b7766d3e3b44e',
            'all_fermions_first_all_scalars_joint':True,
            'target_chi2_strict':args.target_chi2,
            'models_requested':len(queue),'model_fits_target_reached':sum(r['target_met'] for r in rows),
            'model_fits_blocked':sum(r['status'].startswith('BLOCKED') for r in rows),
            'model_fits_above_target':sum(r['status'].startswith('Conditional')
                                          and not r['target_met'] for r in rows),
            'physical_C5_certified':False,
            'latest_source_matching_only_A_-4_and_D_-2':True,
            'results':rows}
        save(output/'summary.json',summary)
    print('Summary:',output/'summary.json',flush=True)
    # A partial batch must never be mistaken for 16 successful target fits.
    return 0 if len(rows)==len(queue) and all(x['target_met'] for x in rows) else 2


if __name__=='__main__':
    raise SystemExit(main())
