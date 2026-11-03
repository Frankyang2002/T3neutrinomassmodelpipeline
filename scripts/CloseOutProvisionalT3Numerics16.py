#!/usr/bin/env python3
"""Close the source-conditional T3 numerical phase from existing frozen files."""
from __future__ import annotations
import argparse
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from Numerical.lo_patterns_16.final_numerical_closeout import (  # noqa: E402
    CloseoutBlocked, DEST, write,
)


def main() -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=ROOT,
                        help='Repository root containing original saved outputs')
    parser.add_argument('--output',type=Path,default=DEST,
                        help='NEW output directory beneath lo_pattern_study16')
    args=parser.parse_args()
    try:
        report,files=write(args.root,args.output)
    except (CloseoutBlocked,OSError,KeyError,IndexError,TypeError,ValueError) as exc:
        print('BLOCKED_FINAL_NUMERICAL_CLOSEOUT:',exc,file=sys.stderr)
        return 2
    target=args.output if args.output.is_absolute() else args.root/args.output
    print('FINAL_NUMERICAL_CLOSEOUT_SOURCE_CONDITIONAL_COMPLETE; '
          f'models={len(report["results"])}/16; scans={report["scan_count"]}/288; '
          f'texture_controls={len(report["controls"])}/64; '
          'coefficient_scenarios=5/5; physical=0/16; '
          f'files={len(files)}; output={target.resolve()}')
    return 0


if __name__=='__main__':
    sys.exit(main())
