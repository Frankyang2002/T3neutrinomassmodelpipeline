"""Command-line rough Stream-1 saved-B input handoff for A/C.

This script intentionally does NOT run RGBeta, Matchete, the Wilson ODE or C matching.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

# Allow direct script invocation from the project root on Windows.
ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))

from Numerical.parallel_t3.rough_inputs.builder import Options,run_batch,MODELS,InputBlocked


def main(argv=None) -> int:
    p=argparse.ArgumentParser(description=__doc__)
    group=p.add_mutually_exclusive_group(required=True)
    group.add_argument('--model',help='A_-4 etc; A_-4 is executed first when --all')
    group.add_argument('--all',action='store_true')
    p.add_argument('--root',type=Path,default=ROOT)
    p.add_argument('--b-run-id',default='B_all16_qualitative_v1')
    p.add_argument('--run-id',default='rough_fermion_first_v1')
    p.add_argument('--gamma-factor',type=float,default=1.,
                   help='Explicit unverified signed GammaCC -> Weyl benchmark factor')
    p.add_argument('--assume-formal-identity',action='store_true',
                   help='Provisional Matchete real scalar source -> canonical identity')
    p.add_argument('--assume-sm-weyl-yukawa',action='store_true',
                   help='Provisional standard SU(2) SM Yukawa embedding')
    p.add_argument('--reduced-potential',action='store_true',
                   help='Project unmapped nonzero running quartics to explicit documented zero')
    p.add_argument('--hard-c5-zero-benchmark',action='store_true',
                   help='Input C5_hard(MF)=0 by benchmark choice, NEVER physical-vanishing claim')
    p.add_argument('--hard-c5-root',type=Path,default=None,
                   help='Relative directory containing <tag>/hard_C5.npz and companion manifest')
    p.add_argument('--overwrite',action='store_true')
    args=p.parse_args(argv)
    if args.model and args.model not in MODELS:p.error('Unknown production model '+args.model)
    if args.run_id in ('.','..') or not args.run_id.replace('_','').replace('-','').isalnum():
        p.error('Unsafe run ID')
    opt=Options(root=args.root.resolve(),B_run_id=args.b_run_id,out_run_id=args.run_id,
        gamma_factor=args.gamma_factor,
        assume_formal_identity=args.assume_formal_identity,
        assume_sm_weyl_yukawa=args.assume_sm_weyl_yukawa,
        reduced_potential=args.reduced_potential,
        hard_c5_zero_benchmark=args.hard_c5_zero_benchmark,
        hard_c5_root=args.hard_c5_root,overwrite=args.overwrite)
    models=list(MODELS) if args.all else [args.model]
    try:result=run_batch(models,opt)
    except InputBlocked as e:
        print('FATAL:',e,file=sys.stderr)
        return 2
    print(f"Written: {result['count_written']}; blocked: {result['count_blocked']}")
    return 0 if result['count_blocked']==0 else 1

if __name__=='__main__':raise SystemExit(main())
