"""Run an all-16 exploratory six-term full-tensor EFT1 Wilson trial.

Use --model A_-4 for the first smoke test, then omit --model for all 16.
All inputs are pre-existing, no Wolfram or Matchete recomputation.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))

from Numerical.running.All16ConditionalFullWilsonTrial import (
    DEFAULT_RELEASE,DEFAULT_V6,DEFAULT_RUNNING,run_all,
)


def main(argv=None)->int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model',help='optional single model A_-4, B_-1, E_+0, ...')
    parser.add_argument('--steps',type=int,default=4,help='midpoint integration steps (1..64)')
    parser.add_argument('--release-dir',type=Path,default=DEFAULT_RELEASE)
    parser.add_argument('--v6-dir',type=Path,default=DEFAULT_V6)
    parser.add_argument('--running-dir',type=Path,default=DEFAULT_RUNNING)
    parser.add_argument('--output-root',type=Path,default=Path('output/general_observable/all16_full_wilson_reduced_potential_v1'))
    parser.add_argument('--vev',type=float,default=246.22)
    parser.add_argument('--save-tensors',action='store_true',
                        help='write complete C(MF),C(MS) NPZ for each model (~MB per file)')
    args=parser.parse_args(argv)
    try:
        result=run_all(ROOT,args.output_root,models=[args.model] if args.model else None,
                       release_dir=args.release_dir,v6_dir=args.v6_dir,
                       running_dir=args.running_dir,steps=args.steps,
                       vev=args.vev,save_tensors=args.save_tensors)
    except (OSError,ValueError,TypeError,KeyError,RuntimeError) as exc:
        print('BLOCKED:',type(exc).__name__,str(exc),file=sys.stderr)
        return 2
    return 0 if result['models_blocked']==0 else 2

if __name__=='__main__':
    raise SystemExit(main())
