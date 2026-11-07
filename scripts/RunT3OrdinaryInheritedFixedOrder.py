"""Batch 4: read-only inherited Weinberg + fixed-order pole audit for D(+2), A(-4).

Consumes Batch 2 native scalar hard templates, Batch 3 validated NF=3
numerical outputs, and archived source-generated first-fermion one-loop EFT1.
No new Wolfram Match; no modification to any prior benchmark.

Usage:
  python -X utf8 -m scripts.RunT3OrdinaryInheritedFixedOrder --check-only
  python -X utf8 -m scripts.RunT3OrdinaryInheritedFixedOrder --run-id inherited_pole_01
  python -X utf8 -m scripts.RunT3OrdinaryInheritedFixedOrder --run-id inherited_pole_01 --resume
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import re
import sys
import traceback

from Numerical.running.T3OrdinaryHardFlavor import (
    MODELS, TAGS, HardFlavorError, file_sha, read_json,
    running_parameters, verify_hard_stage,
)
from Numerical.running.T3OrdinaryInheritedFixedOrder import (
    ARCHIVE_NAMES, ARCHIVE_FILE, audit_inherited_source, assemble_fixed_order,
)

ROOT = Path(__file__).resolve().parents[1]
HARD_ROOT = ROOT/'output/general_sequential_pilot/ordinary_scalar_hard/ordinary_hard_01'
BATCH3_ROOT = ROOT/'output/general_sequential_pilot/ordinary_hard_flavor/hard_flavor_01'
RUNNING_ROOT = ROOT/'output/general_observable/all16_v1'
ARCHIVE_ROOT = ROOT/'output/full16_symbolic'
OUTPUT_ROOT = ROOT/'output/general_sequential_pilot/ordinary_inherited_pole'
STATUS = 'InheritedWeinbergAndCommonCouplingPoleCancellationComputed_Conditional'


def write_json(path: Path, obj: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False,
                               allow_nan=False)+'\n', encoding='utf-8')


def chosen_models(raw: str) -> tuple[str, ...]:
    chosen = tuple(s.strip() for s in raw.split(',') if s.strip())
    if not chosen or len(set(chosen)) != len(chosen) or any(s not in MODELS for s in chosen):
        raise HardFlavorError(f'Use unique comma-separated models drawn from {MODELS}')
    return chosen


def inputs(model: str, args: argparse.Namespace) -> dict[str, Path]:
    return {
        'hard_folder':args.hard_root/TAGS[model],
        'running':args.running_root/TAGS[model]/'running.json',
        'batch3':args.batch3_root/TAGS[model]/'result.json',
        'archived_fermion_one_loop':args.archive_root/ARCHIVE_NAMES[model]/ARCHIVE_FILE,
    }


def preflight(model: str, args: argparse.Namespace) -> tuple[dict,dict,dict]:
    paths = inputs(model,args)
    hard = verify_hard_stage(paths['hard_folder'],model)
    archived = audit_inherited_source(paths['archived_fermion_one_loop'],model)
    running = read_json(paths['running'])
    running_parameters(running,model)
    prior = read_json(paths['batch3'])
    if prior.get('status') != 'HardThreeFlavorAndDirectLogComputed_NotFullC5' or prior.get('model') != model:
        raise HardFlavorError(f'{model}: expected completed Batch 3 result')
    reference = prior.get('independent_archived_LLSS_crosscheck',{})
    if args.require_batch3_parity and reference.get('status') != 'SourceBackedLLSSKernelAgreement':
        raise HardFlavorError(f'{model}: Batch 3 archive parity unavailable/failed')
    if prior.get('source_template',{}).get('source_hard_report_sha256') != hard['source_hard_report_sha256']:
        raise HardFlavorError(f'{model}: Batch 3 source report hash changed')
    input_validation = prior.get('input_validation',{})
    hashes = input_validation.get('frozen_running_input_sha256',{})
    if hashes.get('running') != file_sha(paths['running']):
        raise HardFlavorError(f'{model}: Batch 3 running snapshot hash mismatch')
    files = {
        'native_hard_candidate':paths['hard_folder']/'hard_c5_candidate.txt',
        'native_hard_report':paths['hard_folder']/'hard_scalar_report.json',
        'native_hard_audit':paths['hard_folder']/'hard_scalar_input_audit.json',
        'native_hard_provenance':paths['hard_folder']/'provenance.json',
        'running':paths['running'],
        'batch3':paths['batch3'],
        'archived_fermion_one_loop':paths['archived_fermion_one_loop'],
    }
    meta = {'model':model,'status':'Ready',
            'sha256':{name:file_sha(path) for name,path in files.items()},
            'source_inherited':archived,
            'native_hard_expression_sha256':hard['source_expression_sha256_normalized'],
            'batch3_archive_parity':reference.get('status','not_available'),
            'new_matching_performed':False}
    return meta,running,prior


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--models',default='D_+2,A_-4')
    parser.add_argument('--run-id',default='inherited_pole_01')
    parser.add_argument('--hard-root',type=Path,default=HARD_ROOT)
    parser.add_argument('--batch3-root',type=Path,default=BATCH3_ROOT)
    parser.add_argument('--running-root',type=Path,default=RUNNING_ROOT)
    parser.add_argument('--archive-root',type=Path,default=ARCHIVE_ROOT)
    parser.add_argument('--output-root',type=Path,default=OUTPUT_ROOT)
    parser.add_argument('--vev',type=float,default=246.22)
    parser.add_argument('--require-batch3-parity',action='store_true',
                        help='Require that Batch 3 agreed with archived independent LLSS kernel')
    parser.add_argument('--check-only',action='store_true')
    parser.add_argument('--resume',action='store_true')
    args = parser.parse_args(argv)
    try:
        models = chosen_models(args.models)
        if not re.fullmatch(r'[A-Za-z0-9_-]+',args.run_id):
            raise HardFlavorError('--run-id must be a simple directory label')
        if not math.isfinite(args.vev) or args.vev <= 0:
            raise HardFlavorError('--vev must be positive finite')
    except (HardFlavorError,ValueError) as exc:
        parser.error(str(exc))

    ready = []
    for model in models:
        try:
            meta, running, prior = preflight(model,args)
            print(f'{model}: inherited F-threshold Weinberg triplets verified; '
                  'native scalar hard and Batch 3 reused',flush=True)
            ready.append((model,meta,running,prior))
        except (ValueError,TypeError,KeyError,OSError) as exc:
            print(f'{model}: BLOCKED {type(exc).__name__}: {exc}',flush=True)
            return 2
    if args.check_only:
        return 0

    dest = (args.output_root/args.run_id).resolve()
    settings = {'models':list(models),
                'vev':args.vev,
                'hard_root':str(args.hard_root.resolve()),
                'batch3_root':str(args.batch3_root.resolve()),
                'running_root':str(args.running_root.resolve()),
                'archive_root':str(args.archive_root.resolve()),
                'require_batch3_parity':args.require_batch3_parity}
    if dest.exists():
        if (not args.resume or not (dest/'settings.json').is_file()
                or read_json(dest/'settings.json')!=settings):
            parser.error('Output already exists or resume settings differ; choose a new --run-id')
    else:
        if args.resume:
            parser.error('Cannot resume a run directory that does not exist')
        dest.mkdir(parents=True, exist_ok=False)
        write_json(dest/'settings.json',settings)
    rows = []
    for model,meta,running,prior in ready:
        target = dest/TAGS[model]
        file = target/'result.json'
        if args.resume and file.is_file():
            saved = read_json(file)
            if saved.get('input_validation') != meta:
                print(f'{model}: BLOCKED changed source fingerprint on resume',flush=True)
                return 2
            if saved.get('status') == STATUS:
                print(f'{model}: reused validated inherited/pole result',flush=True)
                rows.append(saved)
                continue
        try:
            result = assemble_fixed_order(
                model=model,
                inherited_path=inputs(model,args)['archived_fermion_one_loop'],
                hard_stage_folder=inputs(model,args)['hard_folder'],
                running_doc=running,batch3_result=prior,vev_gev=args.vev)
            result['input_validation']=meta
            write_json(file,result)
            masses=result['common_fixed_order_masses']['conditional_takagi_masses_ev_ascending']
            print(f'{model}: inherited one-loop included; common-coupling pole cancelled; '
                  f'conditional masses {masses} eV',flush=True)
            rows.append(result)
        except (ValueError,TypeError,KeyError,OSError,OverflowError) as exc:
            err={'model':model,'status':'Blocked','error':f'{type(exc).__name__}: {exc}',
                 'traceback':traceback.format_exc(limit=10)}
            write_json(target/'error.json',err)
            rows.append(err)
            print(f'{model}: BLOCKED {err["error"]}',flush=True)
    errors = sum(r['status']!='InheritedWeinbergAndCommonCouplingPoleCancellationComputed_Conditional'
                 for r in rows)
    summary={'format':'t3-ordinary-inherited-fermion-pole-audit-campaign-v1',
             'created_utc':datetime.now(timezone.utc).isoformat(),
             'model_count':len(rows),
             'conditional_fixed_order_computed_count':len(rows)-errors,
             'error_count':errors,
             'new_native_wolfram_matching_performed':False,
             'source_legacy_F_one_loop_matchete_weinberg_terms_audited':errors==0,
             'source_native_scalar_hard_matchete_template_audited':errors==0,
             'independent_native_NF3_F_one_loop_re_matching_performed':False,
             'physical_full_sequential_c5_certified_count':0,
             'absolute_SU2_spinor_normalization_certified':False,
             'models':[{'model':r['model'],'status':r['status'],
                        'masses':r.get('common_fixed_order_masses',{}).get('conditional_takagi_masses_ev_ascending'),
                        'pole_residual':r.get('checks',{}).get('common_pole_cancellation_relative_residual'),
                        'mixed_endpoint_pole_residual':r.get('checks',{}).get('mixed_endpoint_uncancelled_pole_relative_residual'),
                        'error':r.get('error')} for r in rows]}
    write_json(dest/'summary.json',summary)
    with (dest/'comparison.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=['model','status','m1_eV','m2_eV','m3_eV',
                                            'common_pole_residual','evolved_pole_residual',
                                            'error'])
        writer.writeheader()
        for result in rows:
            masses=result.get('common_fixed_order_masses',{}).get(
                'conditional_takagi_masses_ev_ascending',[None,None,None])
            checks=result.get('checks',{})
            writer.writerow({'model':result['model'],'status':result['status'],
                'm1_eV':masses[0],'m2_eV':masses[1],'m3_eV':masses[2],
                'common_pole_residual':checks.get('common_pole_cancellation_relative_residual'),
                'evolved_pole_residual':checks.get('mixed_endpoint_uncancelled_pole_relative_residual'),
                'error':result.get('error')})
    print('Summary:',dest/'summary.json',flush=True)
    return 0 if errors==0 else 2

if __name__=='__main__':
    sys.exit(main())
