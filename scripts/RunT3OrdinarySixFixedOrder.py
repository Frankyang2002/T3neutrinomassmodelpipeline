"""Batch 5 guarded fixed-order Weinberg numerical assembly for six ordinary T3.

Run AFTER RunT3GeneralOrdinaryScalarHardBatch5 has produced source-backed
one-eigenvalue scalar templates. This runner does NO new Wolfram matching.
Outputs are conditional, restricted to the holomorphic Weinberg sector.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import re
import traceback

from Numerical.running.T3OrdinaryHardFlavor import (
    HardFlavorError, file_sha, read_json, running_parameters,
)
from Numerical.running.T3OrdinarySixFixedOrder import (
    MODELS, TAGS, ARCHIVE, audit_inherited, evaluate,
    inherited_prefactor, verify_native_hard,
)

ROOT=Path(__file__).resolve().parents[1]
DEFAULT_HARD=ROOT/'output/general_sequential_pilot/ordinary_scalar_hard_batch5/ordinary6_hard_01'
DEFAULT_RUNNING=ROOT/'output/general_observable/all16_v1'
DEFAULT_ARCHIVE=ROOT/'output/full16_symbolic'
DEFAULT_FFIRST=ROOT/'output/general_observable/fermion_first_all16_v1'
DEFAULT_REFERENCE=(ROOT/'output/general_observable/lo_joint_scalar_matching/'
                   'exploratory_all16_v6/C5_per_unknown_NGamma_exploratory')
DEFAULT_OUTPUT=ROOT/'output/general_sequential_pilot/ordinary_six_fixed_order'
STATUS='InheritedHardAndCommonCouplingPoleCancelled_Conditional'


def save(path:Path,data:dict)->None:
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(data,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')


def select_models(raw:str)->tuple[str,...]:
    selection=tuple(x.strip() for x in raw.split(',') if x.strip())
    if not selection or len(set(selection))!=len(selection) or any(x not in MODELS for x in selection):
        raise HardFlavorError(f'Choose unique model names from {MODELS}')
    return selection


def paths(model:str,args:argparse.Namespace)->dict[str,Path]:
    root=args.archive_root/ARCHIVE[model]/'data'
    tag=TAGS[model]
    return {
        'hard':args.hard_root/tag,
        'candidate':args.hard_root/tag/'hard_c5_candidate.txt',
        'hard_report':args.hard_root/tag/'hard_scalar_report.json',
        'hard_audit':args.hard_root/tag/'hard_scalar_input_audit.json',
        'provenance':args.hard_root/tag/'provenance.json',
        'running':args.running_root/tag/'running.json',
        'inherited':root/'eft1_after_F_one_loop_inputform.txt',
        'pole':root/'c5_pole_rge_consistency.json',
        'reference':args.reference_root/(tag+'.json'),
    }


def preflight(model:str,args:argparse.Namespace)->tuple[dict,dict,dict]:
    ip=paths(model,args)
    pole=read_json(ip['pole'])
    factor=inherited_prefactor(pole,model)
    hard=verify_native_hard(ip['hard'],model,pole,
                            args.source_root if args.verify_fermion_sources else None)
    source=audit_inherited(ip['inherited'],model,factor)
    running=read_json(ip['running'])
    running_parameters(running,model)
    files={key:val for key,val in ip.items() if key not in ('hard','reference')}
    meta={'model':model,'status':'Ready','source_sha256':{k:file_sha(v)
          for k,v in files.items()},'verified_hard':hard,
          'inherited_source':source,'archive_reference_present':ip['reference'].is_file(),
          'archive_reference_sha256':file_sha(ip['reference']) if ip['reference'].is_file() else None,
          'fermion_sources_hashed':args.verify_fermion_sources}
    return meta,pole,running


def main(argv=None)->int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--models',default=','.join(MODELS))
    parser.add_argument('--run-id',default='ordinary6_fixed_01')
    parser.add_argument('--hard-root',type=Path,default=DEFAULT_HARD)
    parser.add_argument('--running-root',type=Path,default=DEFAULT_RUNNING)
    parser.add_argument('--archive-root',type=Path,default=DEFAULT_ARCHIVE)
    parser.add_argument('--source-root',type=Path,default=DEFAULT_FFIRST)
    parser.add_argument('--reference-root',type=Path,default=DEFAULT_REFERENCE)
    parser.add_argument('--output-root',type=Path,default=DEFAULT_OUTPUT)
    parser.add_argument('--vev',type=float,default=246.22)
    parser.add_argument('--verify-fermion-sources',action='store_true',
                        help='Check every native F-tree source fingerprint')
    parser.add_argument('--require-reference-parity',action='store_true',
                        help='Fail on archived conditional LLSS kernel disagreement')
    parser.add_argument('--check-only',action='store_true')
    parser.add_argument('--resume',action='store_true')
    args=parser.parse_args(argv)
    try:
        selected=select_models(args.models)
        if not re.fullmatch(r'[A-Za-z0-9_-]+',args.run_id):
            raise HardFlavorError('Invalid --run-id')
        if not math.isfinite(args.vev) or args.vev<=0:
            raise HardFlavorError('Positive finite --vev required')
    except (HardFlavorError,ValueError) as exc:
        parser.error(str(exc))
    settings={'models':list(selected),'vev':args.vev,
              'hard_root':str(args.hard_root.resolve()),
              'archive_root':str(args.archive_root.resolve()),
              'running_root':str(args.running_root.resolve()),
              'source_root':str(args.source_root.resolve()),
              'reference_root':str(args.reference_root.resolve()),
              'require_reference_parity':args.require_reference_parity,
              'verify_fermion_sources':args.verify_fermion_sources}
    output=(args.output_root/args.run_id).resolve()
    if not args.check_only:
        if output.exists():
            if (not args.resume or not (output/'settings.json').is_file() or
                read_json(output/'settings.json')!=settings):
                parser.error('Output exists or resume settings changed; choose new --run-id')
        else:
            if args.resume: parser.error('Cannot resume missing output directory')
            output.mkdir(parents=True,exist_ok=False)
            save(output/'settings.json',settings)
    results=[]
    for model in selected:
        ip=paths(model,args)
        try:
            meta,pole,running=preflight(model,args)
            print(f'{model}: source-specific inherited and hard template verified',flush=True)
        except (ValueError,TypeError,KeyError,OSError,RuntimeError) as exc:
            print(f'{model}: BLOCKED source preflight: {type(exc).__name__}: {exc}',flush=True)
            if args.check_only: continue
            err={'model':model,'status':'Blocked','error':f'{type(exc).__name__}: {exc}',
                 'traceback':traceback.format_exc(limit=6)}
            save(output/TAGS[model]/'error.json',err)
            results.append(err)
            continue
        if args.check_only:
            results.append({'model':model,'status':'Ready'})
            continue
        dest=output/TAGS[model]/'result.json'
        if args.resume and dest.is_file():
            saved=read_json(dest)
            if saved.get('input_validation')!=meta:
                print(f'{model}: BLOCKED saved source fingerprints changed',flush=True)
                return 2
            if saved.get('status')==STATUS:
                print(f'{model}: reused validated fixed-order result',flush=True)
                results.append(saved)
                continue
        try:
            result=evaluate(model,ip['hard'],running,pole,ip['inherited'],vev=args.vev,
                   reference_path=ip['reference'],
                   source_root=args.source_root if args.verify_fermion_sources else None)
            if (args.require_reference_parity and
                result['archive_comparison']['status']!='SourceBackedLLSSKernelAgreement'):
                raise HardFlavorError(f'{model}: archived exploratory kernel parity failed: '
                        f"{result['archive_comparison']}")
            result['input_validation']=meta
            save(dest,result)
            results.append(result)
            masses=result['conditional_masses']['conditional_takagi_masses_ev_ascending']
            print(f'{model}: pole cancelled, masses {masses} eV; '
                  f"archive={result['archive_comparison']['status']}",flush=True)
        except (ValueError,TypeError,KeyError,OSError,RuntimeError) as exc:
            err={'model':model,'status':'Blocked',
                 'error':f'{type(exc).__name__}: {exc}',
                 'traceback':traceback.format_exc(limit=10)}
            save(output/TAGS[model]/'error.json',err)
            results.append(err)
            print(f'{model}: BLOCKED numerical assembly: {err["error"]}',flush=True)
    count=sum(x['status']==('Ready' if args.check_only else STATUS) for x in results)
    if args.check_only:
        print(f'Preflight {count}/{len(selected)}')
        return 0 if count==len(selected) else 2
    summary={'format':'t3-ordinary-six-hard-inherited-source-validated-v1',
        'created_utc':datetime.now(timezone.utc).isoformat(),
        'requested_models':list(selected),'computed_count':count,
        'blocked_count':len(selected)-count,'physical_full_C5_certified_count':0,
        'new_native_fermion_matching_performed':False,
        'new_native_scalar_matching_performed':False,
        'source_hard_one_eigenvalue_templates_required':True,
        'models':[{'model':x['model'],'status':x['status'],
            'conditional_masses':x.get('conditional_masses',{}).get('conditional_takagi_masses_ev_ascending'),
            'archive_comparison':x.get('archive_comparison',{}).get('status'),
            'error':x.get('error')} for x in results]}
    save(output/'summary.json',summary)
    with (output/'comparison.csv').open('w',newline='',encoding='utf-8') as stream:
        writer=csv.DictWriter(stream,fieldnames=['model','status','m1_eV','m2_eV','m3_eV',
            'archive_parity','pole_cancellation_residual','scale_variation_residual','error'])
        writer.writeheader()
        for row in results:
            masses=row.get('conditional_masses',{}).get('conditional_takagi_masses_ev_ascending',[None]*3)
            checks=row.get('checks',{})
            writer.writerow({'model':row['model'],'status':row['status'],
                'm1_eV':masses[0],'m2_eV':masses[1],'m3_eV':masses[2],
                'archive_parity':row.get('archive_comparison',{}).get('status'),
                'pole_cancellation_residual':checks.get('common_coupling_pole_cancellation_residual'),
                'scale_variation_residual':checks.get('matching_scale_variation_relative_residual'),
                'error':row.get('error')})
    print('Summary:',output/'summary.json',flush=True)
    return 0 if count==len(selected) else 2

if __name__=='__main__':
    raise SystemExit(main())
