"""Run Workstream 2 V5 all-16 exact Matchete CG/RGBeta source crosswalk."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from Numerical.lo_certification.couplings.named_dictionary_batch_v5 import run_v5_batch


def main()->int:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--repo-root',default='.')
    ap.add_argument('--run-name',default='ws2_matchete_rgbeta_dictionary_v1')
    ap.add_argument('--v4-run-name',default='ws2_named_hhss_v1')
    ap.add_argument('--model',action='append',help='Optional model slug, e.g. B_-1')
    args=ap.parse_args()
    if not args.run_name or args.run_name in ('.','..') or '/' in args.run_name or '\\' in args.run_name:
        ap.error('--run-name must be one bare directory name')
    root=Path(args.repo_root).resolve()
    out=root/'output/general_observable/lo_certification_couplings'/args.run_name
    result=run_v5_batch(root,out,tuple(args.model) if args.model else None,args.v4_run_name)
    print(json.dumps({'output_directory':str(out),
                      'models_processed':result['models_processed'],
                      'all16_source_and_CG_crosswalk_passed':result['all16_source_and_CG_crosswalk_passed'],
                      'extra_RGBeta_quartic_declarations_checked':result['extra_RGBeta_quartic_declarations_checked'],
                      'extra_HHSS_source_symbols_crosswalked':result['extra_HHSS_source_symbols_crosswalked'],
                      'native_Matchete_CG_tensors_parsed':result['native_Matchete_CG_tensors_parsed'],
                      'matching_native_extra_HHSS_channels':result['matching_native_extra_HHSS_channels'],
                      'physically_certified_HHSS_coupling_maps':0,
                      'physical_C5_calculated':False},indent=2))
    return 0 if result['all16_source_and_CG_crosswalk_passed'] else 1

if __name__=='__main__':raise SystemExit(main())
