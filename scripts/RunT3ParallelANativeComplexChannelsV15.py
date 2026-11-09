"""Save immutable v15 complex/conjugate evidence derived from native v14 JSON."""
from __future__ import annotations
import argparse
from pathlib import Path
import sys
from Numerical.parallel_t3.wilson.native_complex_channel_v15 import audit_from_native
from Numerical.parallel_t3.wilson.native_context_quartic_v14 import V14Blocked
from scripts.RunT3ParallelAQuarticNativeContextAuditV14 import write_immutable


def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--project-root',type=Path,default=Path('.'))
    ap.add_argument('--native',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args(argv)
    try:
        root=a.project_root.resolve()
        report=audit_from_native(root,a.native)
        out=a.output if a.output.is_absolute() else root/a.output
        write_immutable(out,report)
        print('v15:',report['status'],'native direct/conjugate unit responses',
              report['counts'].get('unit_direct_zero_conjugate',0),'/',report['counts']['attempted'])
        print('v15 complex native evidence only; physical normalization unverified.')
        return 0
    except (V14Blocked,ValueError,KeyError,OSError,TypeError) as exc:
        print('v15 BLOCKED:',exc,file=sys.stderr)
        return 7

if __name__=='__main__':raise SystemExit(main())
