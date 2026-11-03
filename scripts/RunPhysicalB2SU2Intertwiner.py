#!/usr/bin/env python3
"""Read-only G2 SU(2) intertwiner audit; output in a fresh directory.

G1 native and RGBeta formal reports must already exist. Native Mathematica
is not re-run. This proves only a conditional orientation equivalence, not
a physical Matchete/RGBeta source dictionary.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'Numerical/lo_certification/physical_normalisation'))
from b2_gauge_coupling import AuditBlocked
from b2_su2_intertwiner import (FORMAL_DEFAULT,NATIVE_DEFAULT,OUTPUT_DEFAULT,
                                run,write_report)


def main(argv=None)->int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--native',type=Path,default=NATIVE_DEFAULT)
    parser.add_argument('--formal',type=Path,default=FORMAL_DEFAULT)
    parser.add_argument('--output',type=Path,default=OUTPUT_DEFAULT)
    args=parser.parse_args(argv)
    as_abs=lambda p: p if p.is_absolute() else ROOT/p
    try:
        result=run(ROOT,as_abs(args.native),as_abs(args.formal))
        write_report(as_abs(args.output),result)
    except (AuditBlocked,OSError,ValueError,KeyError,TypeError,IndexError) as exc:
        print('BLOCKED_B2_G2_INTERTWINER:',exc,file=sys.stderr)
        return 3
    print(f"{result['status']}; G1=30/30; Matchete=16/16; "
          f"B/C mixing=6/6; native RGBeta numeric axes=NO; "
          f"physical=0/16; output={as_abs(args.output)}")
    return 0


if __name__=='__main__':
    raise SystemExit(main())
