"""Run the full 45-Weyl A(-4) conditional Wilson EFT1 ODE after fixed RGBeta.

PowerShell from repository root:
  python scripts/RunAminus4ConditionalFullWilson.py --output-root output/a_minus4_next_batch/conditional_full_wilson_v1

No Wolfram rerun. Does NOT integrate out scalars or report a physical mass.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path: sys.path.insert(0,str(ROOT))

from Numerical.running.Aminus4ConditionalFullWilson import ConditionalBlocked,run


def main(argv=None)->int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=ROOT)
    p.add_argument('--output-root',type=Path,
                   default=Path('output/a_minus4_next_batch/conditional_full_wilson_v1'))
    p.add_argument('--points',type=int,default=3)
    p.add_argument('--rtol',type=float,default=5e-5)
    p.add_argument('--atol',type=float,default=1e-20)
    p.add_argument('--method',choices=('RK23','RK45','DOP853'),default='RK23')
    args=p.parse_args(argv)
    try:
        report=run(args.root,args.output_root,
                   points=args.points,rtol=args.rtol,atol=args.atol,method=args.method)
    except (ConditionalBlocked, ValueError,RuntimeError,KeyError,OSError) as exc:
        print('A(-4) conditional Wilson blocked: '+str(exc),file=sys.stderr)
        return 2
    print('Full conditional Wilson basis: 45 Weyl x 12 real; 80,730 complex independent coordinates')
    print('Boundary/master tensor dictionary: PASS')
    print('Coupling background: corrected native RGBeta, seven scales')
    print('Six-term reference-candidate Wilson ODE:',report['ODE_solver']['success'])
    print('Ward residual (MF):',report['gauge_ward_initial_max_relative'])
    print('Ward residual (MS):',max(report['gauge_ward_endpoint_max_relative'].values()))
    print('Full C relative change:',report['full_tensor_relative_change'])
    print('Generated LLHH (not SMEFT C5) max |GeV^-1|:',
          report['max_abs_induced_LLHH_not_SMEFT_C5_gev_inv'])
    print('Matchete signed CC, scalar matching, neutrino masses: NOT CERTIFIED')
    print('Output:',(Path(args.root)/args.output_root/'summary.json').resolve())
    return 0


if __name__=='__main__':raise SystemExit(main())
