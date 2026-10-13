"""One-command presentation-only T3 display build from saved numerical outputs.

Uses the existing BuildDisplayResults implementation and then adds the
interactive correlation report. Does not execute a scan, fit, RGE or matching.
"""
from __future__ import annotations
import argparse
from pathlib import Path
from Numerical.plotting.BuildDisplayResults import (
    PROJECT_ROOT, DEFAULT_DATA_INPUT, DEFAULT_REPORT_INPUT, DEFAULT_OUTPUT, build_display,
)
from Numerical.plotting.BuildInteractiveCorrelations import build, gather


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data-input',type=Path,default=DEFAULT_DATA_INPUT)
    p.add_argument('--report-input',type=Path,default=DEFAULT_REPORT_INPUT)
    p.add_argument('--output',type=Path,default=DEFAULT_OUTPUT)
    p.add_argument('--correlation-source',type=Path,
                   default=PROJECT_ROOT/'Reports/output/full/fitted_sensitivity')
    args=p.parse_args(argv)
    points,_=gather(args.correlation_source)
    # Prevalidate so the existing display folder isn't deleted if no data exists.
    if not points:
        p.error('No saved successful correlation scan points in '+str(args.correlation_source)+
                '; no output replaced. Supply an existing sensitivity scan directory.')
    target=Path(args.output).resolve()
    build_display(args.data_input,args.report_input,target)
    build(args.correlation_source,target/'interactive_correlations.html')
    print('Combined display:',target)
    return 0
if __name__=='__main__':raise SystemExit(main())
