"""General all-16 Workstream-2 coupling/basis dictionary and existing-solver gates.

From repository root:
    python -m scripts.RunLOCouplingCertification --repo-root .
Writes ONLY output/general_observable/lo_certification_couplings/<new-run>.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path


def main(argv=None) -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root',type=Path,default=Path('.'))
    parser.add_argument('--run-name',default='ws2_all16_v1')
    parser.add_argument('--model',action='append',default=None,
                        help='Repeat for a subset; default: all 16')
    parser.add_argument('--allow-local-version-drift',action='store_true',
                        help='Explicitly bypass inspected GitHub blob-SHA check; blocked certification')
    parser.add_argument('--existing-master-coupling-root',type=Path,
                        help='Optional already-complete existing seven-term solver coupling JSON directory')
    parser.add_argument('--existing-joint-background-root',type=Path,
                        help='Optional already-complete existing joint scalar background JSON directory')
    args=parser.parse_args(argv)
    from Numerical.lo_certification.couplings.batch import run_batch, model_paths, write_json
    from Numerical.joint_scalar_matching.fermion_boundary_bridge import slug
    from Numerical.joint_scalar_matching.adapter import ALL_MODELS
    from Numerical.lo_certification.couplings.interfaces import check_existing_solver_inputs
    root=args.repo_root.resolve()
    if not args.run_name or args.run_name in ('.','..') or '/' in args.run_name or '\\' in args.run_name:
        parser.error('--run-name must be a new single path component')
    if bool(args.existing_master_coupling_root) != bool(args.existing_joint_background_root):
        parser.error('Both existing solver input roots must be supplied together')
    output=root/'output/general_observable/lo_certification_couplings'/args.run_name
    result=run_batch(root,output,tuple(args.model) if args.model else None,
                     reject_drift=not args.allow_local_version_drift)
    if args.existing_master_coupling_root:
        release=root/'output/general_observable/lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f'
        for row in result['models']:
            model=row['model']; tag=slug(model)
            boundary=release/tag/'boundary.json'
            master=(root/args.existing_master_coupling_root/tag/'couplings.json').resolve()
            background=(root/args.existing_joint_background_root/f'{tag}.json').resolve()
            check_existing_solver_inputs(model,boundary,master,background,
                                         output/tag/'existing_solver_interface.json')
    print(json.dumps({'output_directory':str(output),
                      'models_processed':result['model_count'],
                      'physical_C5_certified':False,
                      'verified_source_models':sum(x['native_CG_status']=='verified' for x in result['models']),
                      'positive_illustrative_mass_models':sum(x['positive_illustrative_mass'] for x in result['models']),
                      'physical_blocked_models':sum(x['physical_status']=='blocked' for x in result['models'])},indent=2))
    return 0 if result['model_count']==16 and all(x['native_CG_status']=='verified'
                                               for x in result['models']) else 2


if __name__=='__main__':
    raise SystemExit(main())
