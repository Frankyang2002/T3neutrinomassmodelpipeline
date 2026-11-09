"""Read-only checkpoint for v20 + v21; never authorizes Wilson evolution."""
from __future__ import annotations
import argparse,sys
from pathlib import Path
from Numerical.parallel_t3.wilson.bc_mass_bridge_gate_v22 import audit
from Numerical.parallel_t3.wilson.bc_mass_eigen_native_bridge_v21 import BCBridgeBlocked
from scripts.RunT3ParallelABCMassEigenNativeV21 import resolved,write_or_verify
BASE=Path('output/parallel_t3/wilson')
def main(argv=None)->int:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--project-root',type=Path,default=Path('.'))
    ap.add_argument('--v20',type=Path,default=BASE/'full_scalar_census_conditional_v20.json')
    ap.add_argument('--v21',type=Path,default=BASE/'bc_native_mass_eigen_complex_conditional_v21.json')
    ap.add_argument('--output',type=Path,default=BASE/'bc_quartic_mass_readiness_conditional_v22.json')
    ap.add_argument('--resume',action='store_true')
    a=ap.parse_args(argv);root=a.project_root.resolve()
    try:
        report=audit(root,resolved(root,a.v20),resolved(root,a.v21))
        out=resolved(root,a.output)
        print('v22:',write_or_verify(out,report,resume=a.resume),out)
        print('v22:',report['status'])
        print('v22: physical hard matching, signed Weyl, one-loop ADM, scalar threshold and SMEFT remain BLOCKED')
        return 0
    except (BCBridgeBlocked,ValueError,KeyError,TypeError,OSError) as exc:
        print('v22 BLOCKED:',exc,file=sys.stderr);return 7
if __name__=='__main__':raise SystemExit(main())
