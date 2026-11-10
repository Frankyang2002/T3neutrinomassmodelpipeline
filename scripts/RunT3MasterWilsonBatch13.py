#!/usr/bin/env python3
"""Batch 13: all-F-first A(-4)/D(-2) formal full-tensor EFT1 Wilson mixing.

Use --preflight to verify source files without a lengthy tensor calculation.
Results are conditional, not physical C5. Original Batches 9/11/12 untouched.
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

from Numerical.parallel_t3.complete.two_model_master_wilson_batch13 import (
    B9,B11,D_EXTRA,REGISTRY,MODELS,FORMAT,run_one)

DEFAULT=Path('output/parallel_t3/thesis_master_wilson_batch13_v1')


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    s=p.add_mutually_exclusive_group(required=True)
    s.add_argument('--all',action='store_true')
    s.add_argument('--model',choices=MODELS)
    p.add_argument('--preflight',action='store_true')
    p.add_argument('--output',type=Path,default=DEFAULT)
    p.add_argument('--master-to-native-factor',type=float,default=1.,
                   help='Explicit formal master-to-native factor; NOT certified')
    p.add_argument('--signed-factor',type=float,default=1.,
                   help='Explicit formal GammaCC to SMEFT factor; NOT certified')
    p.add_argument('--beta-mode',choices=('fonseca_eq22_six_term',
                   'ward_covariant_candidate'),default='fonseca_eq22_six_term')
    p.add_argument('--skip-smeft',action='store_true')
    p.add_argument('--debug',action='store_true')
    a=p.parse_args(argv)
    models=MODELS if a.all else (a.model,)
    needed=[]
    for model in models:
        slug=model.replace('_-','_m')
        needed.extend((ROOT/B9/model/'result.json',ROOT/B9/model/(slug+'.npz'),
                       ROOT/B9/model/(slug+'.json'),ROOT/B11/model/'result.json'))
    if 'D_-2' in models:
        needed.extend((ROOT/D_EXTRA,ROOT/REGISTRY))
    for f in needed:
        if not f.is_file():
            print('MISSING',f,flush=True)
        elif a.preflight:
            print('FOUND',f,flush=True)
    if not all(f.is_file() for f in needed): return 1
    if a.preflight:
        print('Batch13 preflight passed (not a physics calculation).',flush=True)
        return 0
    dest=a.output if a.output.is_absolute() else ROOT/a.output
    records=[]
    for model in models:
        try:
            print('RUN',model,'full formal EFT1 master-tensor beta',flush=True)
            o=run_one(ROOT,model,master_to_native=a.master_to_native_factor,
                      signed_factor=a.signed_factor,beta_mode=a.beta_mode,
                      include_smeft=not a.skip_smeft,output=dest/model)
            records.append({'model':model,'status':o['status'],
                            'master_native_log_discrepancy':o['master_vs_native_log_relative_difference'],
                            'source_full_ADM_certified':False})
            print(model,o['status'],flush=True)
        except Exception as e:
            records.append({'model':model,'status':'Blocked_SourceOrConventionGap',
                            'error':str(e),'error_type':type(e).__name__})
            print(model,'Blocked_SourceOrConventionGap',type(e).__name__+':',e,flush=True)
            if a.debug: traceback.print_exc()
    dest.mkdir(parents=True,exist_ok=True)
    summary={'format':FORMAT+'-summary','models':records,
             'successful_count':sum(r['status'].startswith('FullTensor') for r in records),
             'expected_count':len(records),'physical_C5_certified':False,
             'signed_factor_assumed':a.signed_factor,
             'master_to_native_factor_assumed':a.master_to_native_factor,
             'B9_B11_B12_outputs_modified':False}
    target=dest/'summary.json'
    target.write_text(json.dumps(summary,indent=2,sort_keys=True,allow_nan=False)+'\n',
                      encoding='utf-8')
    print('Summary:',target,flush=True)
    return 0 if summary['successful_count']==len(records) else 1


if __name__=='__main__':
    raise SystemExit(main())
