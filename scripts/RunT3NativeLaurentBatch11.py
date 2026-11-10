"""Source-native A(-4)/D(-2) fermion-hard + 2-scalar Laurent closure.

python -X utf8 -m scripts.RunT3NativeLaurentBatch11 --all --preflight
python -X utf8 -m scripts.RunT3NativeLaurentBatch11 --all --debug
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import traceback

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))

from Numerical.parallel_t3.complete.two_model_native_laurent_batch11 import (
    B10, FORMAT, MODELS, SOFT_PATHS, NATIVE_STATUS, B1_SUBDIR,run_one,
)

DEST=Path('output/parallel_t3/thesis_native_laurent_batch11_v1')


def write(path:Path,doc:dict):
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_name(path.name+'.partial')
    tmp.write_text(json.dumps(doc,indent=2,sort_keys=True,
                  ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')
    tmp.replace(path)


def main(argv=None)->int:
    p=argparse.ArgumentParser(description=__doc__)
    sel=p.add_mutually_exclusive_group(required=True)
    sel.add_argument('--all',action='store_true')
    sel.add_argument('--model',choices=MODELS)
    p.add_argument('--project-root',type=Path,default=ROOT)
    p.add_argument('--batch10',type=Path,default=B10,
                   help='Prior successful Batch10 fixed output, relative to project root')
    p.add_argument('--output',type=Path,default=DEST,
                   help='New output subfolder under output/parallel_t3/')
    p.add_argument('--preflight',action='store_true')
    p.add_argument('--debug',action='store_true')
    p.add_argument('--illustrative-smeft',action='store_true',
                   help='Optional, unverified source-normalization example; NOT physics')
    args=p.parse_args(argv)
    root=args.project_root.resolve()
    out=(root/args.output).resolve()
    safe=(root/'output/parallel_t3').resolve()
    if safe not in out.parents or out==safe:
        p.error('Output must be a child of output/parallel_t3/')
    selected=MODELS if args.all else (args.model,)
    if args.preflight:
        missing=False
        for model in selected:
            for f in [root/args.batch10/model/'result.json',
                      root/SOFT_PATHS[model],root/NATIVE_STATUS[model],
                      root/B1_SUBDIR/model/'manifest.json']:
                exists=f.is_file();missing |= not exists
                print(('FOUND ' if exists else 'MISSING ')+str(f),flush=True)
        return 2 if missing else 0
    if out.exists():p.error('Refusing to overwrite '+str(out))
    rows=[]
    for model in selected:
        try:
            result=run_one(root,model,batch10=args.batch10,
                       include_illustrative=args.illustrative_smeft)
            success=result['status'].startswith('NativeHardSoftLaurentAndScaleLogClosed')
        except (ValueError,OSError,RuntimeError,ImportError,TypeError,
                KeyError,ArithmeticError) as exc:
            success=False
            result={'format':FORMAT,'model':model,
                    'status':'Blocked_NativeLaurentOrConventionGap',
                    'reason':type(exc).__name__+': '+str(exc),
                    'physical_MSbar_C5_claimed':False}
            if args.debug:result['traceback_tail']=traceback.format_exc(limit=18)
        write(out/model/'result.json',result)
        rows.append({'model':model,'status':result['status'],
                     'pole_closure':result.get('source_pole_closure_relative'),
                     'scale_invariance':result.get('source_scale_partition_invariance_relative'),
                     'reason':result.get('reason')})
        print(model,result['status'],flush=True)
        if not success:print(result['reason'],file=sys.stderr,flush=True)
    completed=sum(x['status'].startswith('NativeHardSoftLaurentAndScaleLogClosed')
                  for x in rows)
    write(out/'summary.json',{'format':FORMAT+'-summary',
          'models':rows,'success_count':completed,'expected_models':len(rows),
          'physical_MSbar_C5_claimed':False,'previous_outputs_changed':False})
    print('Summary:',out/'summary.json')
    return 0 if completed==len(rows) else 2

if __name__=='__main__':raise SystemExit(main())
