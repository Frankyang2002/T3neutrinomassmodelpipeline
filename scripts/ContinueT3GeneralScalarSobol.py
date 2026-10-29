"""Bounded, resumable *new-seed* Sobol continuation of the 16x3 candidate campaign.

Uses RunT3GeneralScalar48Campaign.run_case directly (its original 18-real-Yukawa
Sobol, mass screens, model-specific matching and least-squares refinement).
Previously accepted cases are reused; unresolved eligible cases get fresh seeds.
No historical fits, source campaign, beta exports or original scripts are edited.

This is conditional general-quartic matching, NOT fully matched general-scalar EFT.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import re
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

from scripts import RunT3GeneralScalar48Campaign as bulk

ROOT = bulk.ROOT
OUT = bulk.OUT
CASES = tuple(bulk.STUDIES)
ACCEPTED = frozenset(('CandidateConstrainedAccepted', 'ExistingConditionalAccepted'))
TARGET = 5.86


def require(test, message):
    if not test:
        raise ValueError('Sobol continuation: '+message)


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, document):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name+'.tmp')
    tmp.write_text(json.dumps(document, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    tmp.replace(path)


def acceptable(row):
    try:
        score = float(row.get('chi2'))
        return row.get('status') in ACCEPTED and math.isfinite(score) and score <= TARGET
    except (ValueError, TypeError):
        return False


def choose(current, proposed):
    """Never replace an accepted result; otherwise retain the lower finite chi2."""
    if acceptable(current):
        return current
    if acceptable(proposed):
        return proposed
    def score(row):
        try:
            x=float(row.get('chi2'))
            return x if math.isfinite(x) else math.inf
        except (TypeError, ValueError):
            return math.inf
    return proposed if score(proposed) < score(current) else current


def selected_cells(rows, models, max_rounds, attempts, cases=CASES):
    return [(m, c) for m in models for c in cases
            if not acceptable(rows[(m,c)]) and attempts.get((m,c),0) < max_rounds]


def run(args):
    require(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,63}',args.run_id), 'invalid run ID')
    require(args.source_campaign != args.run_id, 'source and destination must differ')
    require(args.sobol_points>=2 and args.sobol_points & (args.sobol_points-1)==0,
            '--sobol-points must be a power of two')
    require(1<=args.max_rounds<=1000, 'invalid maximum rounds')
    require(0<=args.max_nfev<=1000, 'invalid local optimisation budget')
    require(0 < abs(args.extra_real) <= 0.5, 'invalid extra quartic boundary')
    registry = bulk.read(bulk.REGISTRY)
    models = tuple(registry['models'])
    if args.models:
        requested = tuple(dict.fromkeys(x.strip() for x in args.models.split(',') if x.strip()))
        require(requested and all(m in models for m in requested), 'unknown model selection')
        models = tuple(m for m in models if m in requested)
    cases = tuple(x.strip() for x in args.cases.split(',') if x.strip()) if args.cases else CASES
    require(cases and len(set(cases))==len(cases) and all(c in CASES for c in cases), 'invalid scenario selection')
    source = OUT/args.source_campaign
    require((source/'campaign.json').is_file(), 'missing source campaign')
    original = bulk.read(source/'campaign.json')
    require(original.get('format')==bulk.FORMAT and len(original.get('rows',[]))==48,
            'source is not the original 48-case campaign')
    original_rows = {(r['model'], r['study']):r for r in original['rows']}
    require(len(original_rows)==48, 'duplicate source model/scenario')
    destination = OUT/args.run_id
    signature = {'algorithm':'fresh-seed-run_case-v1', 'source_campaign':args.source_campaign,
                 'source_campaign_sha256':sha(source/'campaign.json'),
                 'sobol_points':args.sobol_points,'max_nfev':args.max_nfev,
                 'seed':args.seed,'extra_real':args.extra_real,'target_chi2':TARGET,
                 'models':list(models),'cases':list(cases)}
    signature_path = destination/'continuation_settings.json'
    if destination.exists():
        require(args.resume or args.report_only, 'destination exists; supply --resume')
        require(signature_path.is_file() and read(signature_path)==signature,
                'resume settings or source campaign changed; use new run-id')
    else:
        require(not args.report_only, 'cannot report before initialization')
        destination.mkdir(parents=True)
        save(signature_path,signature)
    # Existing beta JSON is copied by hash, never exported or regenerated.
    for model in models:
        if model in bulk.REUSED:
            continue
        name = bulk.safe(model)+'.json'
        src = source/'betas'/name
        dst = destination/'betas'/name
        require(src.is_file(), 'missing saved beta: '+str(src))
        if dst.exists():
            require(sha(dst)==sha(src), 'beta changed since continuation '+model)
        else:
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src,dst)
    progress_path = destination/'continuation_progress.json'
    if progress_path.exists():
        progress=read(progress_path)
        require(progress.get('settings_sha256')==sha(signature_path), 'checkpoint signature mismatch')
    else:
        progress={'settings_sha256':sha(signature_path),'attempts':{},'rows':original['rows'],
                  'history':[], 'last_updated_utc':None}
    rows={(r['model'],r['study']):r for r in progress['rows']}
    require(set(rows)==set(original_rows), 'invalid checkpoint cells')
    attempts = {(m,c):int(progress['attempts'].get(m+'|'+c,0)) for m,c in rows}
    # Materialize original per-case fit/config for existing accepted A(-2) donors
    # and all other standalone cases. The B/C historical reuses have their own
    # separately hashed source records, and have no bulk per-cell fit file.
    for model,case in rows:
        if model in bulk.REUSED:
            continue
        dest_case=destination/case/bulk.safe(model)
        old_case=source/case/bulk.safe(model)
        if not (dest_case/'fit.json').exists() and attempts[(model,case)]==0:
            if (old_case/'fit.json').is_file():
                dest_case.mkdir(parents=True,exist_ok=True)
                shutil.copy2(old_case/'fit.json',dest_case/'fit.json')
                if (old_case/'fitted_config.json').is_file():
                    shutil.copy2(old_case/'fitted_config.json',dest_case/'fitted_config.json')

    def checkpoint():
        progress['rows']=[rows[(m,c)] for m in registry['models'] for c in CASES]
        progress['attempts']={m+'|'+c:v for (m,c),v in attempts.items() if v}
        progress['last_updated_utc']=datetime.now(timezone.utc).isoformat()
        save(progress_path,progress)
        aggregate={'format':bulk.FORMAT,'rows':progress['rows'],
                   'model_count':16,'study_count':3,'conditional_only':True,
                   'general_quadratic_extensions_not_in_bulk_runner':True,
                   'no_sequential_wilson_matching':True}
        save(destination/'campaign.json',aggregate)
        with (destination/'campaign.csv').open('w',encoding='utf-8',newline='') as handle:
            writer=csv.DictWriter(handle,fieldnames=list(bulk.row_of({'model':next(iter(registry['models'])),
                'case':CASES[0],'status':'NotRunInThisCampaign'}).keys()))
            writer.writeheader();writer.writerows(progress['rows'])
        accepted=sum(acceptable(r) for r in rows.values())
        missing=len(rows)-accepted
        (destination/'summary.md').write_text(
            '# General-scalar conditional Sobol continuation\n\n'
            f'Accepted: {accepted}/48, unresolved: {missing}/48.\n\n'
            f'Fresh-seed attempts: {sum(attempts.values())}.\n\n'
            'Conditional general-quartic physics only; additional holomorphic quadratic matching '
            'and full sequential Wilson EFT matching are not included.\n',encoding='utf-8')

    checkpoint()
    for model,case in selected_cells(rows,models,args.max_rounds,attempts,cases):
        if args.report_only:
            continue
        if model in bulk.REUSED:
            print('SKIP special matched model',model,case,'(existing B/C candidate workflow only)',flush=True)
            continue
        while not acceptable(rows[(model,case)]) and attempts[(model,case)]<args.max_rounds:
            round_no=attempts[(model,case)]+1
            seed=args.seed + 1009 * list(registry['models']).index(model) + 9176*CASES.index(case) + 100003*round_no
            folder=destination/'attempts'/case/bulk.safe(model)/f'round_{round_no:03d}'
            require(not folder.exists(), 'checkpoint round exists without recorded attempt: '+str(folder))
            try:
                entry,new,reused=bulk.inventory(registry,model)
                require(not reused, 'special reused quartic model')
                beta_path=destination/'betas'/(bulk.safe(model)+'.json')
                payload=bulk.read(beta_path)
                kernel,matching_path,kind=bulk.matching(model)
                result=bulk.run_case(model,case,registry,payload,kernel,matching_path,kind,
                    beta_path,folder,sobol_points=args.sobol_points,
                    max_nfev=args.max_nfev,random_seed=seed,extra_real=args.extra_real)
                candidate=bulk.row_of(result)
                selected=choose(rows[(model,case)],candidate)
                if selected is candidate:
                    rows[(model,case)]=candidate
                    chosen_file=destination/case/bulk.safe(model)/'fit.json'
                    chosen_file.parent.mkdir(parents=True,exist_ok=True)
                    shutil.copy2(folder/'fit.json',chosen_file)
                    cfg=folder/'fitted_config.json'
                    if cfg.is_file():shutil.copy2(cfg,chosen_file.parent/'fitted_config.json')
                status=result['status']
                score=(result.get('best_mass_safe') or {}).get('chi2')
            except Exception as exc:
                status='RoundFailed'
                score=None
                folder.mkdir(parents=True,exist_ok=True)
                save(folder/'failure.json',{'error':str(exc)[:2000],'model':model,'case':case})
            attempts[(model,case)]=round_no
            progress['history'].append({'model':model,'case':case,'round':round_no,
                                        'seed':seed,'status':status,'chi2':score,
                                        'selected_status':rows[(model,case)]['status'],
                                        'selected_chi2':rows[(model,case)]['chi2']})
            checkpoint()
            print(f'{model} {case} round={round_no}/{args.max_rounds} status={status} '
                  f'chi2={score} best={rows[(model,case)]["chi2"]}',flush=True)
    checkpoint()
    accepted=sum(acceptable(r) for r in rows.values())
    print('COMPLETED bounded continuation:',accepted,'/48 accepted; output:',destination/'campaign.csv',flush=True)
    if accepted!=48:
        print('Unresolved cells are retained as poor fit / no mass-safe result. '
              'Additional Sobol rounds are not a guarantee of convergence.',flush=True)
    return 0


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source-campaign',default='overnight_six_v1_bulk48')
    p.add_argument('--run-id',required=True)
    p.add_argument('--sobol-points',type=int,default=128)
    p.add_argument('--max-nfev',type=int,default=100)
    p.add_argument('--max-rounds',type=int,default=8)
    p.add_argument('--seed',type=int,default=20261101)
    p.add_argument('--extra-real',type=float,default=-0.01)
    p.add_argument('--models',help='Optional comma-separated target model names')
    p.add_argument('--cases',help='Comma-separated target scenarios, e.g. small_uv_frozen_5e-9,large_uv_frozen_0p5')
    p.add_argument('--resume',action='store_true')
    p.add_argument('--report-only',action='store_true')
    return run(p.parse_args(argv))

if __name__=='__main__':
    sys.exit(main())
