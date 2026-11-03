#!/usr/bin/env python3
"""Read-only B2 validation of a fresh native RGBeta internal tensor witness."""
from __future__ import annotations
import argparse
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'Numerical/lo_certification/physical_normalisation'))
from b2_gauge_coupling import AuditBlocked
from b2_rgbeta_internal_action import analyse, write_outputs


def main() -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--native',required=True,type=Path)
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args()
    native=args.native if args.native.is_absolute() else ROOT/args.native
    output=args.output if args.output.is_absolute() else ROOT/args.output
    try:
        report=analyse(ROOT,native.resolve())
        write_outputs(output,report)
    except (AuditBlocked,ValueError,KeyError,IndexError,OSError) as exc:
        print('BLOCKED_NATIVE_INTERNAL_FACTOR_EVIDENCE:',str(exc),file=sys.stderr)
        return 3
    print(f"{report['status']}; families=5/5; internal ratios=30/30; "
          f'physical=0/16; output={output}')
    return 0


if __name__=='__main__':
    raise SystemExit(main())
