"""One-command source-conditional all-16 T3 direct-UV pattern study. No EFT."""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from Numerical.lo_patterns_16.pattern_study import run_study, OUTPUT_REL


def main() -> int:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--project-root',type=Path,default=ROOT)
    ap.add_argument('--input-dir',type=Path,help='Existing full frozen candidate Weyl JSON inputs (16 signed names)')
    ap.add_argument('--output-dir',type=Path,help='NEW output directory inside owned lo_pattern_study16')
    ap.add_argument('--physical-dictionary',type=Path,help='Reserved for coordinator certificate; fail-closed until B1/B2 trust contract exists')
    ap.add_argument('--no-scans',action='store_true')
    ap.add_argument('--no-plots',action='store_true')
    args=ap.parse_args()
    root=args.project_root.resolve()
    try:
        result=run_study(root,args.output_dir or root/OUTPUT_REL,args.input_dir,args.physical_dictionary,
                         do_scans=not args.no_scans,do_plots=not args.no_plots)
    except (OSError,ValueError,RuntimeError) as exc:
        print(f'PATTERN STUDY BLOCKED: {type(exc).__name__}: {exc}',file=sys.stderr)
        return 2
    counts=result['counts']
    print('Conditional direct UV:',counts['conditional_evaluated'],'/16;',
          'blocked:',counts['blocked'],'/16;',
          'physically certified: 0/16')
    return 0 if counts['blocked']==0 and counts['scan_blocked']==0 else 2


if __name__=='__main__':raise SystemExit(main())
