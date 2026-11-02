"""Workstream-2 all-16 baseline HHSS portal completion and extra SU2 reference audit.

Usage: python -m scripts.RunLOQuarticCompletion --repo-root . --run-name ws2_quartic_completion_v1
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path


def main(argv=None)->int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root',type=Path,default=Path('.'))
    parser.add_argument('--run-name',default='ws2_quartic_completion_v1')
    parser.add_argument('--model',action='append',default=None)
    args=parser.parse_args(argv)
    if not args.run_name or args.run_name in ('.','..') or '/' in args.run_name or '\\' in args.run_name:
        parser.error('--run-name must be one new directory basename')
    root=args.repo_root.resolve()
    from Numerical.lo_certification.couplings.quartic_completion_batch import run_completion_batch
    output=root/'output/general_observable/lo_certification_couplings'/args.run_name
    result=run_completion_batch(root,output,tuple(args.model) if args.model else None)
    print(json.dumps({
        'output_directory':str(output),
        'models_processed':result['models_processed'],
        'all16_model_audits_passed':result['all16_model_audits_passed'],
        'baseline_portal_directions_completed':result['baseline_portal_directions_completed'],
        'extra_HHSS_registry_entries':result['extra_HHSS_registry_entries'],
        'physical_C5_calculated':False,
    },indent=2))
    return 0 if result['all16_model_audits_passed'] else 2

if __name__=='__main__':
    raise SystemExit(main())
