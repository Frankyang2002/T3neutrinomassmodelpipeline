"""Batch 6 native HHSS EWSB Hessian candidate-alias UV/EFT diagnostic.

Run from repository root:
 python -X utf8 scripts/RunLONativePortalBatch6.py --check-only
 python -X utf8 scripts/RunLONativePortalBatch6.py --run-id source01
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0,str(ROOT))
from Numerical.lo_certification.benchmark.native_portal_batch6 import preflight,calculate


def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--project-root',type=Path,default=ROOT)
    ap.add_argument('--batch4-run-id',default='source01')
    ap.add_argument('--batch5-run-id',default='source01')
    ap.add_argument('--run-id',default='source01')
    ap.add_argument('--check-only',action='store_true')
    args=ap.parse_args(argv)
    try:
        result=(preflight(args.project_root,args.batch4_run_id,args.batch5_run_id)
                if args.check_only else calculate(args.project_root,args.batch4_run_id,
                                                  args.batch5_run_id,args.run_id))
    except (ValueError,OSError,RuntimeError,KeyError,ImportError,TypeError) as exc:
        print('BATCH6 BLOCKED: '+str(exc),file=sys.stderr)
        return 2
    print(json.dumps(result,indent=2,allow_nan=False))
    return 0


if __name__=='__main__':
    raise SystemExit(main())
