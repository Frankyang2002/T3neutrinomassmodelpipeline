#!/usr/bin/env python3
"""Batch 16 native/ RGBeta scalar tensor-coordinate span and conditional rerun.

No source or archived output is modified. B16 outputs are in a separate tree.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import json
import sys
import traceback

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0,str(ROOT))

from Numerical.parallel_t3.complete.two_model_scalar_coordinate_batch16 import (
    MODELS, FORMAT, ROOT_B13, run_one)
from Numerical.parallel_t3.complete.two_model_joint_scalar_batch9 import (
    quartic_source_path)

DEFAULT=Path('output/parallel_t3/thesis_scalar_coordinate_batch16_v1')


def main(argv=None)->int:
    parser=argparse.ArgumentParser(description=__doc__)
    which=parser.add_mutually_exclusive_group(required=True)
    which.add_argument('--all',action='store_true')
    which.add_argument('--model',choices=MODELS)
    parser.add_argument('--preflight',action='store_true')
    parser.add_argument('--skip-smeft',action='store_true',
                        help='Skip conditional 100 GeV neutrino post-processing')
    parser.add_argument('--output',type=Path,default=DEFAULT)
    parser.add_argument('--debug',action='store_true')
    args=parser.parse_args(argv)
    models=MODELS if args.all else (args.model,)
    required=[]
    for model in models:
        required.extend([
          ROOT/ROOT_B13/model/'result.json',
          ROOT/ROOT_B13/model/'full_tensor_candidate.npz',
          quartic_source_path(ROOT,model),
        ])
    missing=[f for f in required if not f.is_file()]
    for f in required:
        if args.preflight or f in missing:
            print('MISSING' if f in missing else 'FOUND',f,flush=True)
    if missing:return 1
    if args.preflight:
        print('Batch16 preflight passed (not a physics calculation).',flush=True)
        return 0
    dest=(args.output if args.output.is_absolute() else ROOT/args.output)
    records=[]
    for model in models:
        try:
            print('RUN',model,'native/RGBeta quartic coordinate span',flush=True)
            out=run_one(ROOT,model,output=dest/model,
                        include_smeft=not args.skip_smeft)
            records.append({
                'model':model,'status':out['status'],
                'missing_reference_directions':out['unresolved_reference_directions'],
                'batch13_reproduction_relative':out['batch13_reproduction_relative'],
                'complete_physical_C5':False,
            })
            print(model,out['status'],flush=True)
        except Exception as exc:
            records.append({'model':model,'status':'Blocked_SourceOrBasisGap',
                            'error':str(exc),'error_type':type(exc).__name__})
            print(model,'Blocked_SourceOrBasisGap',type(exc).__name__+':',exc,flush=True)
            if args.debug:traceback.print_exc()
    dest.mkdir(parents=True,exist_ok=True)
    summary={
        'format':FORMAT+'-summary','models':records,
        'completed_count':sum(x['status'].startswith('ScalarCoordinate') for x in records),
        'expected_count':len(records),'physical_C5_certified':False,
        'prior_archived_outputs_modified':False,
    }
    (dest/'summary.json').write_text(json.dumps(summary,indent=2,sort_keys=True,
        allow_nan=False)+'\n',encoding='utf-8')
    print('Summary:',dest/'summary.json',flush=True)
    return 0 if summary['completed_count']==len(records) else 1

if __name__=='__main__':
    raise SystemExit(main())
