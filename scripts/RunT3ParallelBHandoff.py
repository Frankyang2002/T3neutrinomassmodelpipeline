"""Audit saved B all-16 UV/EFT1 plus F-tree outputs for A/C consumption.

Does NOT run RGBeta/Wolfram, Wilson evolution, or scalar matching.
Prints per-model verification statuses and writes only B's handoff_index.json.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))

from Numerical.parallel_t3.background.contracts import MODELS,tag
from Numerical.parallel_t3.background.handoff import audit_all


def main(argv=None) -> int:
    ap=argparse.ArgumentParser(description=__doc__)
    group=ap.add_mutually_exclusive_group(required=True)
    group.add_argument('--all',action='store_true')
    group.add_argument('--model',choices=MODELS)
    ap.add_argument('--run-id',default='B_all16_qualitative_v1')
    ap.add_argument('--no-write',action='store_true',
                    help='Validate only; do not create a handoff index')
    args=ap.parse_args(argv)
    models=MODELS if args.all else (args.model,)
    output=None if args.no_write else (ROOT/'output'/'parallel_t3'/'background'/
                                       args.run_id/'handoff_index.json')
    result=audit_all(ROOT,args.run_id,models,output=output)
    for model,row in result['model_results'].items():
        print(f"{model}: {row['status']}")
        if row['status'].startswith('Blocked_'):
            print(f"  {row['error']}")
        elif row.get('source_drift_warnings'):
            for note in row['source_drift_warnings']:
                print(f'  WARNING: {note}')
    print(f"Handoff verified: {result['passed_count']} / {result['requested_count']}")
    if output:print('Handoff index:',output)
    return 0 if result['passed_count']==result['requested_count'] else 2


if __name__=='__main__':
    raise SystemExit(main())
