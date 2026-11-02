"""Run A_-4 common fitted Yukawas through canonical RGE -> F-tree -> joint scalars.

PowerShell: python -X utf8 scripts/RunLOCommonInputBatch4.py --check-only
            python -X utf8 scripts/RunLOCommonInputBatch4.py --run-id source01
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from Numerical.lo_certification.benchmark.common_input_batch4 import (
    Batch4Error, preflight, calculate)


def main(argv=None) -> int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root',type=Path,default=ROOT)
    p.add_argument('--canonical-run-id',default='full16_general_v1',
                   help='Existing canonical full-general symbolic source (read only)')
    p.add_argument('--check-only',action='store_true')
    p.add_argument('--run-id',default='source01',help='New unique output ID')
    args=p.parse_args(argv)
    try:
        if args.check_only:
            result=preflight(args.project_root,symbolic_run_id=args.canonical_run_id)
        else:
            result=calculate(args.project_root,args.run_id,
                             symbolic_run_id=args.canonical_run_id)
        print(json.dumps(result,indent=2,allow_nan=False))
        return 0
    except (ValueError,RuntimeError,KeyError,FileNotFoundError,ImportError) as exc:
        print('BATCH4 BLOCKED: '+str(exc),file=sys.stderr)
        return 2


if __name__=='__main__':
    raise SystemExit(main())
