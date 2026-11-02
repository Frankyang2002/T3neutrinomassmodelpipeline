"""Run V4 source-explicit named HHSS coupling projection audit (read-only sources)."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from Numerical.lo_certification.couplings.named_hhss_batch import run_named_hhss_batch


def main()->int:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--repo-root',default='.')
    ap.add_argument('--run-name',default='ws2_named_hhss_v1')
    ap.add_argument('--v3-run-name',default='ws2_quartic_completion_v1')
    ap.add_argument('--model',action='append',help='Optional model selection, e.g. B_-1')
    args=ap.parse_args()
    if not args.run_name or args.run_name in ('.','..') or '/' in args.run_name or '\\' in args.run_name:
        ap.error('--run-name must be one bare directory name')
    root=Path(args.repo_root).resolve()
    out=root/'output/general_observable/lo_certification_couplings'/args.run_name
    result=run_named_hhss_batch(root,out,tuple(args.model) if args.model else None,args.v3_run_name)
    print(json.dumps({'output_directory':str(out),'models_processed':result['models_processed'],
      'extra_HHSS_source_symbols_examined':result['extra_HHSS_source_symbols_examined'],
      'all16_source_projections_passed':result['all16_source_projections_passed'],
      'physical_named_maps_certified':result['physical_named_maps_certified'],
      'physical_C5_calculated':False},indent=2))
    return 0 if result['all16_source_projections_passed'] else 1


if __name__=='__main__':raise SystemExit(main())
