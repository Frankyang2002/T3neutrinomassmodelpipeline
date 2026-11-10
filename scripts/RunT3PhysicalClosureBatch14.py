#!/usr/bin/env python3
"""Batch 14: source-driven combined audit of the last three physical EFT gaps.

This audits, but does NOT assume or certify, GammaCC/Weyl normalization,
complete one-loop finite matching and complete gauge-covariant EFT1 ADM.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import traceback

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
from Numerical.parallel_t3.complete.two_model_physical_closure_batch14 import (
    MODELS,FORMAT,preflight_paths,run_one)

DEFAULT=Path('output/parallel_t3/thesis_physical_closure_batch14_v1')

def main(argv=None)->int:
    p=argparse.ArgumentParser(description=__doc__)
    s=p.add_mutually_exclusive_group(required=True)
    s.add_argument('--all',action='store_true')
    s.add_argument('--model',choices=MODELS)
    p.add_argument('--preflight',action='store_true')
    p.add_argument('--output',type=Path,default=DEFAULT)
    p.add_argument('--skip-native-quartic-map',action='store_true',
                   help='Only run native matching and Ward checks; leave coupling-map probe undone')
    p.add_argument('--generator-details',action='store_true',
                   help='Save residual for each SU2/U1/SU3 generator')
    p.add_argument('--debug',action='store_true')
    a=p.parse_args(argv)
    models=MODELS if a.all else (a.model,)
    needed={model:preflight_paths(ROOT,model) for model in models}
    missing=[]
    for model,files in needed.items():
        for path in files:
            if not path.is_file():
                print('MISSING',path,flush=True)
                missing.append(path)
            elif a.preflight:
                print('FOUND',path,flush=True)
    if missing: return 1
    if a.preflight:
        print('Batch14 preflight passed (not a physics calculation).',flush=True)
        return 0
    output=a.output if a.output.is_absolute() else ROOT/a.output
    results=[]
    for model in models:
        print('RUN',model,'three physical closure gates',flush=True)
        try:
            record=run_one(ROOT,model,output=output/model,
                compare_native_quartics=not a.skip_native_quartic_map,
                include_generator_details=a.generator_details)
            r={'model':model,'status':record['status'],
               'ward_maxima':record['gauge_Wilson_gate']['max_relative_Ward_residuals'],
               'native_F_sector_pole_relative':record['finite_matching_gate'][
                   'source_channel_hard_soft_pole_residual_relative'],
               'physical_C5_certified':False}
            results.append(r)
            print(model,record['status'],flush=True)
        except Exception as e:
            results.append({'model':model,'status':'Blocked_AuditExecution',
                            'error':str(e),'exception':type(e).__name__})
            print(model,'Blocked_AuditExecution',str(e),file=sys.stderr,flush=True)
            if a.debug:traceback.print_exc()
    output.mkdir(parents=True,exist_ok=True)
    summary={'format':FORMAT+'-summary','models':results,
             'audit_completed':sum(x['status'].startswith('ThreePhysics') for x in results),
             'total':len(results),'physical_C5_certified':False,
             'prior_outputs_modified':False}
    target=output/'summary.json'
    target.write_text(json.dumps(summary,sort_keys=True,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('Summary:',target,flush=True)
    return 0 if summary['audit_completed']==len(results) else 1

if __name__=='__main__':
    raise SystemExit(main())
