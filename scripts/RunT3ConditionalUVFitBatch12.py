#!/usr/bin/env python3
"""A(-4)/D(-2) pole-closed source-sector -> SMEFT -> optional full-UV-RG fit.

Read-only on existing source and prior results. All output goes to a NEW
output/parallel_t3 child. Does NOT represent a physical complete T3 C5.

python -X utf8 -m scripts.RunT3ConditionalUVFitBatch12 --all --preflight
python -X utf8 -m scripts.RunT3ConditionalUVFitBatch12 --all --signed-factor 1 --debug
python -X utf8 -m scripts.RunT3ConditionalUVFitBatch12 --all --signed-factor 1 --fit --max-nfev 18 --output output/parallel_t3/thesis_conditional_uvfit_batch12_fit_v1 --debug
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import traceback

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))

from Numerical.parallel_t3.complete.two_model_conditional_uvfit_batch12 import (
    FORMAT, MODELS, B11, B1, B2, LAURENT, SOFT_PATHS, BENCHMARK, run_one,
)

DEST=Path('output/parallel_t3/thesis_conditional_uvfit_batch12_v1')


def write(path:Path,data:dict)->None:
    path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_name(path.name+'.partial')
    temp.write_text(json.dumps(data,indent=2,sort_keys=True,
                    ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')
    temp.replace(path)


def check_inputs(root:Path,models:tuple[str,...],b11:Path,*,fit:bool)->bool:
    required=[root/BENCHMARK]
    if fit:
        required.append(root/'scripts/RunT3ThesisGeneralBatch1.py')
    for model in models:
        required.extend([
            root/b11/model/'result.json', root/SOFT_PATHS[model],
            root/B1/model/'manifest.json',
            root/B1/model/'trajectories.npz',
            root/B1/model/'fermion_matching_input.json',
            root/B1/model/'scalar_matching_input.json',
            root/B2/model/'one_loop_operator_ledger.json',
            root/B2/model/'manifest.json',
            root/LAURENT/model/'laurent_audit.json',
            root/LAURENT/model/'manifest.json',
        ])
    ok=True
    for path in required:
        status=path.is_file()
        print(('FOUND ' if status else 'MISSING ')+str(path),flush=True)
        ok &= status
    return ok


def main(argv:list[str]|None=None)->int:
    p=argparse.ArgumentParser(description=__doc__)
    sel=p.add_mutually_exclusive_group(required=True)
    sel.add_argument('--all',action='store_true')
    sel.add_argument('--model',choices=MODELS)
    p.add_argument('--project-root',type=Path,default=ROOT)
    p.add_argument('--batch11',type=Path,default=B11)
    p.add_argument('--output',type=Path,default=DEST)
    p.add_argument('--signed-factor',type=float,required=False,default=1.,
                   help='EXPLICIT conditional GammaCC-to-Weyl multiplier; 1 is NOT certified')
    p.add_argument('--fit',action='store_true',
                   help='Re-run full renormalizable UV and EFT1 RG at each of 18-real UV Yukawa trials')
    p.add_argument('--max-nfev',type=int,default=25)
    p.add_argument('--yukawa-bound',type=float,default=.005)
    p.add_argument('--preflight',action='store_true')
    p.add_argument('--debug',action='store_true')
    args=p.parse_args(argv)
    root=args.project_root.resolve()
    out=(root/args.output).resolve()
    parent=(root/'output/parallel_t3').resolve()
    if parent not in out.parents or out==parent:
        p.error('Output must be a new subdirectory of output/parallel_t3')
    if not 0<abs(args.signed_factor)<1000 or args.max_nfev<1:
        p.error('Bad signed-factor or max-nfev')
    selected=MODELS if args.all else (args.model,)
    if args.preflight:
        return 0 if check_inputs(root,selected,args.batch11,fit=args.fit) else 2
    if out.exists():p.error('Refusing to overwrite existing '+str(out))
    rows=[]
    for model in selected:
        try:
            print('RUN',model, 'reference and UV fit' if args.fit else 'reference SMEFT',flush=True)
            report=run_one(root,model,signed_factor=args.signed_factor,
                           fit=args.fit,max_nfev=args.max_nfev,
                           yukawa_bound=args.yukawa_bound,batch11_dir=args.batch11)
            ok=report['status'].startswith('Conditional')
        except (Exception, ) as exc:
            ok=False
            report={'format':FORMAT,'model':model,
                    'status':'Blocked_SourceOrNumericalInterface',
                    'reason':f'{type(exc).__name__}: {exc}',
                    'physical_neutrino_prediction_claimed':False}
            if args.debug:report['traceback_tail']=traceback.format_exc(limit=20)
        write(out/model/'result.json',report)
        row={'model':model,'status':report['status'],
             'reason':report.get('reason'),
             'nufit_chi2':(report.get('baseline') or {}).get('nufit_chi2_gaussian')}
        rows.append(row)
        print(model,report['status'],flush=True)
        if not ok:print(report['reason'],file=sys.stderr,flush=True)
    count=sum(x['status'].startswith('Conditional') for x in rows)
    write(out/'summary.json',{'format':FORMAT+'-summary','models':rows,
         'source_matching_only_conditional':True,
         'signed_factor_assumed':args.signed_factor,
         'UV_18_real_fit_requested':args.fit,
         'full_finite_MSbar_C5_certified':False,
         'success_count':count,'expected_count':len(selected),
         'prior_outputs_modified':False})
    print('Summary:',out/'summary.json')
    return 0 if count==len(selected) else 2


if __name__=='__main__':raise SystemExit(main())
