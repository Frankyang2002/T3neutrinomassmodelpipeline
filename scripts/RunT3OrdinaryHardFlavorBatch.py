"""Numerically lift the saved ordinary scalar hard C5 templates to NF=3.

Read-only inputs: completed Matchete hard templates and the existing full
renormalisable running snapshots. No wolframscript or new matching is run.

    python -X utf8 -m scripts.RunT3OrdinaryHardFlavorBatch --run-id hard_flavor_01
    python -X utf8 -m scripts.RunT3OrdinaryHardFlavorBatch --run-id hard_flavor_01 --resume
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import sys
import traceback

from Numerical.running.T3OrdinaryHardFlavor import (
    MODELS, TAGS, HardFlavorError, compare_archived_v6, evaluate,
    file_sha, read_json, running_parameters, verify_hard_stage,
)

ROOT = Path(__file__).resolve().parents[1]
HARD_ROOT = ROOT/'output/general_sequential_pilot/ordinary_scalar_hard/ordinary_hard_01'
RUNNING_ROOT = ROOT/'output/general_observable/all16_v1'
SOURCE_ROOT = ROOT/'output/general_observable/fermion_first_all16_v1'
REFERENCE_ROOT = (ROOT/'output/general_observable/lo_joint_scalar_matching/'
                  'exploratory_all16_v6/C5_per_unknown_NGamma_exploratory')
OUTPUT_ROOT = ROOT/'output/general_sequential_pilot/ordinary_hard_flavor'


def save(path: Path, doc: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(doc, indent=2, ensure_ascii=False,
                               allow_nan=False)+'\n', encoding='utf-8')


def models_from_csv(raw: str) -> tuple[str, ...]:
    parsed = tuple(x.strip() for x in raw.split(',') if x.strip())
    if not parsed or len(set(parsed))!=len(parsed) or any(x not in MODELS for x in parsed):
        raise HardFlavorError(f'Only unique comma-separated {MODELS} supported')
    return parsed


def input_paths(model: str, args: argparse.Namespace) -> dict[str, Path]:
    tag = TAGS[model]
    hard = args.hard_root/tag
    return {'hard_dir':hard,
            'hard_candidate':hard/'hard_c5_candidate.txt',
            'hard_report':hard/'hard_scalar_report.json',
            'hard_audit':hard/'hard_scalar_input_audit.json',
            'hard_provenance':hard/'provenance.json',
            'running':args.running_root/tag/'running.json',
            'archived':args.reference_root/(tag+'.json')}


def preflight(model: str, args: argparse.Namespace) -> tuple[dict, dict]:
    paths = input_paths(model,args)
    source_root = args.source_root if args.verify_fermion_sources else None
    src = verify_hard_stage(paths['hard_dir'],model,source_root=source_root)
    raw = read_json(paths['running'])
    p = running_parameters(raw,model)
    source_ids = {key:file_sha(path) for key,path in paths.items()
                  if key not in ('hard_dir','archived')}
    meta = {'status':'Ready', 'model':model,
            'native_matching_source':src,
            'frozen_running_input_sha256':source_ids,
            'reference_present':paths['archived'].is_file(),
            'archived_reference_sha256':file_sha(paths['archived']) if paths['archived'].is_file() else None,
            'heavy_masses_gev':p['mf'].tolist(),
            'scalar_matching_scale_gev':p['mu_s'],
            'source_root_verified':args.verify_fermion_sources}
    return meta,raw


def main(argv=None) -> int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--models',default='D_+2,A_-4')
    p.add_argument('--run-id',default='hard_flavor_01')
    p.add_argument('--hard-root',type=Path,default=HARD_ROOT)
    p.add_argument('--running-root',type=Path,default=RUNNING_ROOT)
    p.add_argument('--source-root',type=Path,default=SOURCE_ROOT)
    p.add_argument('--reference-root',type=Path,default=REFERENCE_ROOT)
    p.add_argument('--output-root',type=Path,default=OUTPUT_ROOT)
    p.add_argument('--vev',type=float,default=246.22)
    p.add_argument('--check-only',action='store_true')
    p.add_argument('--resume',action='store_true')
    p.add_argument('--verify-fermion-sources',action='store_true',
                   help='Verify original saved F-only tree/seed source hashes also')
    p.add_argument('--require-reference-parity',action='store_true',
                   help='Block if absent/mismatching archived independent LLSS kernels')
    args=p.parse_args(argv)
    try:
        models=models_from_csv(args.models)
        if not re.fullmatch(r'[A-Za-z0-9_-]+',args.run_id):
            raise HardFlavorError('--run-id must be a simple directory label')
        if args.vev <= 0 or args.vev != args.vev:
            raise HardFlavorError('--vev must be positive and finite')
    except (HardFlavorError,ValueError) as e:
        p.error(str(e))

    prepared=[]
    for model in models:
        try:
            meta,raw=preflight(model,args)
            print(f'{model}: verified Matchete hard template; '
                  f'3 nondegenerate MF eigenvalues; no matching rerun',flush=True)
            prepared.append((model,meta,raw))
        except (OSError,ValueError,KeyError,TypeError) as e:
            print(f'{model}: BLOCKED {type(e).__name__}: {e}',flush=True)
            return 2
    if args.check_only:
        return 0
    out=(args.output_root/args.run_id).resolve()
    settings={'models':list(models),'vev':args.vev,
              'hard_root':str(args.hard_root.resolve()),
              'running_root':str(args.running_root.resolve()),
              'reference_root':str(args.reference_root.resolve()),
              'verify_fermion_sources':args.verify_fermion_sources,
              'source_root':str(args.source_root.resolve()),
              'require_reference_parity':args.require_reference_parity}
    if out.exists():
        if not args.resume or not (out/'settings.json').is_file() or \
                read_json(out/'settings.json')!=settings:
            p.error('Output already exists or resume settings disagree; use new --run-id')
    else:
        if args.resume:
            p.error('Cannot resume an output directory that does not yet exist')
        out.mkdir(parents=True,exist_ok=False)
        save(out/'settings.json',settings)
    rows=[]
    errors=0
    for model,meta,raw in prepared:
        target=out/TAGS[model]
        result_path=target/'result.json'
        if args.resume and result_path.is_file():
            prior=read_json(result_path)
            if prior.get('input_validation')!=meta:
                print(f'{model}: BLOCKED: source fingerprints changed; refusing stale resume')
                return 2
            if prior.get('status')=='HardThreeFlavorAndDirectLogComputed_NotFullC5':
                print(f'{model}: reused validated saved numerical result')
                rows.append(prior)
                continue
        try:
            result=evaluate(model=model,
                            stage_folder=input_paths(model,args)['hard_dir'],
                            running_doc=raw,vev_gev=args.vev,
                            source_root=args.source_root if args.verify_fermion_sources else None)
            result['input_validation']=meta
            ref_path=input_paths(model,args)['archived']
            if ref_path.is_file():
                parity=compare_archived_v6(result,ref_path)
            else:
                parity={'status':'ArchivedReferenceMissing',
                        'reference_path':str(ref_path),
                        'comparison_is_not_required_for_computation':True}
            result['independent_archived_LLSS_crosscheck']=parity
            if args.require_reference_parity and parity['status']!='SourceBackedLLSSKernelAgreement':
                raise HardFlavorError(f'{model}: independent archived LLSS parity failed: {parity}')
            save(result_path,result)
            masses=result['hard_plus_running_subset_metrics']['conditional_takagi_masses_ev_ascending']
            print(f'{model}: three-flavor hard+direct-log LLSS subset '
                  f'{masses} eV; archive={parity["status"]}',flush=True)
            rows.append(result)
        except (OSError,ValueError,KeyError,TypeError) as exc:
            errors+=1
            error={'model':model,'status':'Blocked',
                   'error':f'{type(exc).__name__}: {exc}',
                   'traceback':traceback.format_exc(limit=8)}
            save(target/'error.json',error)
            rows.append(error)
            print(f'{model}: BLOCKED: {error["error"]}',flush=True)
    summary={'format':'t3-ordinary-three-fermion-hard-direct-log-campaign-v1',
             'created_utc':datetime.now(timezone.utc).isoformat(),
             'requested_models':list(models),
             'computed_count':sum(r['status']=='HardThreeFlavorAndDirectLogComputed_NotFullC5'
                                  for r in rows),
             'error_count':errors,
             'full_sequential_c5_certified_count':0,
             'new_wolfram_matching_performed':False,
             'hard_three_eigenvalue_sum_performed':errors==0,
             'full_inherited_fermion_one_loop_included':False,
             'independent_MSbar_counterterm_proved':False,
             'physical_SU2_spinor_normalization_proved':False,
             'archive_matrices_never_double_counted':True,
             'models':[{'model':r['model'],'status':r['status'],
                        'hard_plus_running_subset_metrics':r.get('hard_plus_running_subset_metrics'),
                        'independent_archived_LLSS_crosscheck':r.get('independent_archived_LLSS_crosscheck'),
                        'error':r.get('error')} for r in rows]}
    save(out/'summary.json',summary)
    with (out/'comparison.csv').open('w',newline='',encoding='utf-8') as handle:
        writer=csv.DictWriter(handle,fieldnames=['model','status','m1_eV','m2_eV','m3_eV',
                                                  'hard_pole_norm_GeV_inv',
                                                  'scale_cancellation_residual',
                                                  'reference_parity','error'])
        writer.writeheader()
        for row in rows:
            x=row.get('hard_plus_running_subset_metrics',{})
            masses=x.get('conditional_takagi_masses_ev_ascending',[None]*3)
            checks=row.get('checks',{})
            writer.writerow({'model':row['model'],'status':row['status'],
                'm1_eV':masses[0],'m2_eV':masses[1],'m3_eV':masses[2],
                'hard_pole_norm_GeV_inv':checks.get('pole_residue_norm_gev_inverse'),
                'scale_cancellation_residual':checks.get('hard_plus_direct_log_matching_scale_variation_residual'),
                'reference_parity':row.get('independent_archived_LLSS_crosscheck',{}).get('status'),
                'error':row.get('error')})
    print('Summary:',out/'summary.json',flush=True)
    return 0 if errors==0 else 2


if __name__=='__main__':
    sys.exit(main())
