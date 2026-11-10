#!/usr/bin/env python3
"""Run one independent 18-real-Yukawa conditional rough fit per T3 model.

This IS NOT a full physical finite-MSbar one-loop fit. A model is blocked on
missing source inputs, non-Weinberg projections or failed numerical matching.
Outputs are new; no existing archive/source files are overwritten.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import traceback

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0,str(ROOT))

from Numerical.parallel_t3.complete.all16_rough_yukawa_fit import (
    FORMAT, BENCHMARK, FreshAll16Evaluator, PipelineBlocked,
    all_models, normalize_models, optimize_one, preflight, read_json, save_json)

DEFAULT_OUTPUT = Path('output/parallel_t3/all16_sobol_yukawa_chi2_001_v2')


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    select=p.add_mutually_exclusive_group(required=True)
    select.add_argument('--plan',action='store_true')
    select.add_argument('--preflight',action='store_true')
    select.add_argument('--run',action='store_true')
    p.add_argument('--models',nargs='+',default=None,
                   help='Subset of 16 models; default ALL 16')
    p.add_argument('--project-root',type=Path,default=ROOT)
    p.add_argument('--batch1',type=str,
                   default='output/batch1_general_contract/batch1_native_v1')
    p.add_argument('--output',type=Path,default=DEFAULT_OUTPUT)
    p.add_argument('--target-chi2',type=float,default=0.01)
    p.add_argument('--max-nfev',type=int,default=12,
                   help='SciPy objective iterations PER local least-squares start')
    p.add_argument('--sobol-points',type=int,default=16,
                   help='Scrambled Sobol trials after selecting a benchmark seed')
    p.add_argument('--screen-points',type=int,default=18,
                   help='Direct Yukawa-coordinate sensitivity probes (0 to 18)')
    p.add_argument('--active-schedule',type=int,nargs='+',default=[5,9,18],
                   help='Active parameter counts for successive local fit stages')
    p.add_argument('--starts',type=int,default=2,
                   help='Independent local least-squares starting points per stage')
    p.add_argument('--perturbation',type=float,default=1.0,
                   help='Half-width about benchmark: perturbation * max(abs(y),2e-5)')
    p.add_argument('--seed',type=int,default=17,
                   help='Reproducible base Sobol seed; model index is added')
    p.add_argument('--max-trials',type=int,default=250,
                   help='Hard cap on ALL unique fresh-pipeline trials per model, including seeds')
    p.add_argument('--yukawa-bound',type=float,default=0.005,
                   help='Absolute bound on each actual UV Yukawa matrix entry')
    p.add_argument('--gamma-factor',type=float,default=1.,
                   help='ASSUMED, not certified Matchete GammaCC-to-Weyl conversion')
    p.add_argument('--mass-bilinears',choices=('zero','nonzero-exploratory'),
                   default='zero',help='zero matches thesis shared UV benchmark; nonzero-exploratory uses B mass seed')
    p.add_argument('--wilson-solver',choices=('RK2','DOP853'),default='RK2')
    p.add_argument('--wilson-steps',type=int,default=2)
    p.add_argument('--mu-low',type=float,default=100.)
    p.add_argument('--full-scalar-potential',action='store_true',
                   help='BLOCK when native quartics lack RGBeta-to-tensor conversion')
    p.add_argument('--keep-trial-files',action='store_true',
                   help='Preserve LARGE per-trial intermediate outputs for debugging')
    p.add_argument('--debug',action='store_true')
    args=p.parse_args(argv)
    root=args.project_root.resolve()
    models=normalize_models(args.models)
    print(f'T3 rough one-fit-per-model: {len(models)} model(s); target chi2 < {args.target_chi2:g}',flush=True)
    if args.plan:
        print(json.dumps({'format':FORMAT,'models':list(models),
              'model_count':len(models),'model_optimisations':len(models),
              'initial_seed_policy':'two fixed benchmark evaluations; select better; ONE model optimisation',
              'search_method':'scrambled Sobol + sensitivity screening + bounded least-squares multistart',
              'sobol_points':args.sobol_points,'screen_points':args.screen_points,
              'active_parameter_stages':args.active_schedule,'starts_per_stage':args.starts,
              'fitting_parameters':'18 REAL UV Yukawa entries per model',
              'target_chi2_strict':args.target_chi2,
              'RG_and_Wilson_recomputed_per_trial':True,
              'finite_F_one_loop_matching_included':False,
              'physical_predictions':False},indent=2),flush=True)
        return 0
    entries=[preflight(root,m,batch1=args.batch1) for m in models]
    if args.preflight:
        print(json.dumps({'all_inputs_found':all(e['input_ready'] for e in entries),
                         'models':entries},indent=2),flush=True)
        return 0 if all(e['input_ready'] for e in entries) else 2
    missing_by_model={e['model']:e['missing'] for e in entries if not e['input_ready']}
    for label,missing in missing_by_model.items():
        print('MISSING SOURCE INPUTS',label,missing,flush=True)
    # Do not abort the entire campaign: write a blocked row for each missing
    # model and continue with the other source-backed model evaluations.
    if not (0 < args.target_chi2 < 100 and args.max_nfev >= 1 and
            args.max_trials >= 2 and args.wilson_steps >= 2 and
            args.sobol_points >= 0 and 0 <= args.screen_points <= 18 and
            args.starts >= 1 and args.active_schedule and
            all(1 <= n <= 18 for n in args.active_schedule) and
            0 < args.perturbation <= 20):
        p.error('Invalid optimisation or ODE settings')
    cfg=read_json(root/BENCHMARK)
    out=(root/args.output).resolve()
    allowed=(root/'output/parallel_t3').resolve()
    if out==allowed or allowed not in out.parents:
        p.error('Output must be a new child directory of output/parallel_t3')
    if out.exists():
        p.error('Output exists; choose a new --output (no overwrite)')
    out.mkdir(parents=True)
    rows=[]
    for model in models:
        print('RUN',model,'fresh all-16 rough UV+F+EFT1+S+SMEFT optimisation',flush=True)
        model_file=out/model/'result.json'
        def checkpoint(doc):
            save_json(out/model/'checkpoint.json',doc)
            b=doc.get('best',{})
            print('BEST',model,'chi2',b.get('chi2'),'eval',b.get('eval'),flush=True)
        try:
            if model in missing_by_model:
                raise PipelineBlocked('Missing model-specific source inputs: '+
                                      ', '.join(missing_by_model[model]))
            evaluator=FreshAll16Evaluator(
                root,model,cfg,batch1=args.batch1,
                gamma_factor=args.gamma_factor,
                mu_low=args.mu_low,vev=float(cfg['vev_GeV']),
                reduced_potential=not args.full_scalar_potential,
                mass_bilinears=args.mass_bilinears,
                wilson_solver=args.wilson_solver,
                wilson_steps=args.wilson_steps,
                keep_trial_files=args.keep_trial_files)
            result=optimize_one(evaluator,cfg,model=model,
                target_chi2=args.target_chi2,max_nfev=args.max_nfev,
                max_trials=args.max_trials,yukawa_bound=args.yukawa_bound,
                sobol_points=args.sobol_points,screen_points=args.screen_points,
                active_schedule=tuple(args.active_schedule),starts=args.starts,
                perturbation=args.perturbation,rng_seed=args.seed,
                callback=checkpoint)
        except Exception as exc:
            result={'status':'BLOCKED_SOURCE_OR_NUMERICAL','model':model,
                    'target_met':False,'error':f'{type(exc).__name__}: {exc}',
                    'physical_C5_certified':False}
            if args.debug:
                result['traceback_tail']=traceback.format_exc(limit=15)
        result.update({'format':FORMAT,'model':model,
            'requested_target_chi2':args.target_chi2,
            'calculation_is_conditional_not_physical':True,
            'input_benchmark_file':str(BENCHMARK),
            'wilson_solver':args.wilson_solver,
            'reduced_potential':not args.full_scalar_potential,
            'mass_bilinears_policy':args.mass_bilinears,
            'optimizer_search_method':'Sobol+sensitivity+bounded_TRF',
            'sobol_points_requested':args.sobol_points,
            'active_schedule_requested':args.active_schedule})
        save_json(model_file,result)
        rows.append({'model':model,'status':result['status'],
                     'target_met':bool(result.get('target_met',False)),
                     'best_chi2':result.get('best_chi2'),
                     'error':result.get('error')})
        save_json(out/'summary.json',{'format':FORMAT+'-summary','models':rows,
                 'count_completed':len(rows),'count_requested':len(models),
                 'count_target_met':sum(x['target_met'] for x in rows),
                 'conditional_not_physical':True,
                 'threshold_sequence':['UV','all F integrated together','EFT1 with both S',
                                       'both S jointly at MS','SMEFT'],
                 'GammaCC_factor_assumed':args.gamma_factor})
        print(model,result['status'],'chi2',result.get('best_chi2'),flush=True)
    print('Summary:',out/'summary.json',flush=True)
    return 0 if all(r['target_met'] for r in rows) else 1

if __name__=='__main__':
    raise SystemExit(main())
