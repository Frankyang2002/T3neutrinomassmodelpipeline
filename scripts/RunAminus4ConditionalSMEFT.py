"""One common scalar threshold, conditional A(-4) C5 and neutrino mass proxies.

No Wilson ODE or Wolfram rerun: use corrected production RGBeta plus the
completed full 45-Weyl conditional Wilson EFT1 trajectory.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
from Numerical.joint_scalar_matching.Aminus4ConditionalSMEFT import (
    run, DEFAULT_WILSON,DEFAULT_OUTPUT,ConditionalScalarBlocked)


def main(argv=None)->int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=ROOT)
    p.add_argument('--wilson-root',type=Path,default=DEFAULT_WILSON)
    p.add_argument('--output-root',type=Path,default=DEFAULT_OUTPUT)
    p.add_argument('--vev-gev',type=float,default=246.22)
    args=p.parse_args(argv)
    try:
        report=run(args.root,args.output_root,wilson_root=args.wilson_root,
                   vev_gev=args.vev_gev)
    except (ConditionalScalarBlocked,ValueError,RuntimeError,KeyError,OSError,
            FloatingPointError) as exc:
        print('A(-4) conditional joint-scalar matching BLOCKED: '+str(exc),file=sys.stderr)
        return 2
    m=report['matching']; n=m['neutrino_mass_proxy']
    print('F-only threshold -> corrected EFT1 Wilson -> one COMMON scalar threshold: PASS')
    print('Integrated scalars: 2 real S1 and 6 real S2 components; 36 eigenstate pairs')
    print('MSbar heavy scalar masses (GeV):',m['heavy_scalar_masses_gev'])
    print('EFT1 induced LLHH max |C5| (GeV^-1):',m['max_abs_direct_C5_gev_inv'])
    print('Finite 2-scalar bubble max |delta C5| (GeV^-1):',m['max_abs_scalar_bubble_C5_gev_inv'])
    print('Combined conditional max |C5| (GeV^-1):',m['max_abs_conditional_C5_gev_inv'])
    print('Scalar-bubble SU2 residual:',m['scalar_bubble_weinberg_tensor_residual'])
    print('MSbar threshold log-derivative residual:',m['matching_scale_log_derivative_relative_residual'])
    print('Conditional mass proxies at v=%.5f GeV (eV):'%n['vev_gev'],
          n['takagi_singular_values_eV'])
    print('Physical Matchete/SMEFT/neutrino-mass certification: NO')
    print('Report:',(Path(args.root)/args.output_root/'summary.json').resolve())
    return 0

if __name__=='__main__':raise SystemExit(main())
