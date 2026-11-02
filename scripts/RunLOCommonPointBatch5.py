"""Compare fitted A_-4 direct UV loop against Batch4 conditional EFT, source controlled.

Usage from project root:
 python -X utf8 scripts/RunLOCommonPointBatch5.py --check-only
 python -X utf8 scripts/RunLOCommonPointBatch5.py --run-id source01
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0,str(ROOT))
from Numerical.lo_certification.benchmark.common_point_batch5 import (preflight,calculate)

def main(argv=None)->int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project-root',type=Path,default=ROOT)
    parser.add_argument('--batch4-run-id',default='source01')
    parser.add_argument('--check-only',action='store_true')
    parser.add_argument('--run-id',default='source01')
    a=parser.parse_args(argv)
    try:
        obj=(preflight(a.project_root,a.batch4_run_id) if a.check_only else
             calculate(a.project_root,a.batch4_run_id,a.run_id))
    except (OSError,ValueError,RuntimeError,KeyError,ImportError,TypeError) as e:
        print('BATCH5 BLOCKED: '+str(e),file=sys.stderr)
        return 2
    print(json.dumps(obj,indent=2,allow_nan=False))
    return 0

if __name__=='__main__':
    raise SystemExit(main())
