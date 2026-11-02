"""Run Batch7 native Matchete scalar-potential Hessian / SU2 weight-basis audit."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0,str(ROOT))

from Numerical.lo_certification.benchmark import native_scalar_basis_batch7 as batch


def main(argv=None) -> int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root',type=Path,default=ROOT)
    p.add_argument('--batch4-run-id',default='source01')
    p.add_argument('--batch5-run-id',default='source01')
    p.add_argument('--batch6-run-id',default='source01')
    p.add_argument('--run-id',default='source01')
    p.add_argument('--check-only',action='store_true')
    a=p.parse_args(argv)
    try:
        if a.check_only:
            info=batch.preflight(a.project_root,a.batch4_run_id,a.batch5_run_id,a.batch6_run_id)
        else:
            info=batch.calculate(a.project_root,a.batch4_run_id,a.batch5_run_id,
                                 a.batch6_run_id,a.run_id)
    except (ValueError,RuntimeError,KeyError,TypeError,OSError) as exc:
        print('BATCH7 BLOCKED:',str(exc),file=sys.stderr)
        return 2
    print(json.dumps(info,indent=2,allow_nan=False))
    return 0


if __name__=='__main__':
    raise SystemExit(main())
