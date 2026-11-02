"""All-16 conditional EFT1 seven-term source + one two-scalar threshold.

Without --coupling-root and --background-root: write 16 source-checked input
contracts, no physical numbers.  With BOTH roots: validate all external tensors
before exclusive writes, then output one conditional O(hbar) result per model.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path: sys.path.insert(0,str(ROOT))

from RGE.running.intermediate.ProvisionalSourceKernel import MODELS, tag, SourceKernelError
from RGE.running.intermediate.SourceBackedEFT1Batch import read_json, digest
from RGE.running.intermediate.IntegratedEFT1JointScalarBridge import (
    FORMAT, CONTRACT_FORMAT, contracts, evaluate,
)

RELEASE=Path('output/general_observable/lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f')
OUTPUT=Path('output/general_observable/lo_eft1_mixing/integrated_eft1_joint_v1')


def main(argv=None) -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--release-root',type=Path,default=RELEASE)
    parser.add_argument('--coupling-root',type=Path,default=None,
                        help='Read-only <tag>/couplings.json, exact source SHA in each file')
    parser.add_argument('--background-root',type=Path,default=None,
                        help='Read-only <tag>.json + <tag>.npz, same muS / quartics')
    parser.add_argument('--output-root',type=Path,default=OUTPUT)
    parser.add_argument('--model',choices=MODELS)
    parser.add_argument('--trajectory-points',type=int,default=21)
    parser.add_argument('--check-only',action='store_true')
    a=parser.parse_args(argv)
    if (a.coupling_root is None)!=(a.background_root is None):
        parser.error('--coupling-root and --background-root must be supplied together')
    out=a.output_root.resolve()
    allowed=(ROOT/'output/general_observable/lo_eft1_mixing').resolve()
    if out==allowed or allowed not in out.parents:
        parser.error(f'Outputs restricted to a NEW subdirectory under {allowed}')
    release=a.release_root.resolve()
    source_summary=read_json(release/'summary.json')
    if source_summary.get('format') != 't3-lo-fermion-boundary-source-resolved-provisional-v1' or \
       source_summary.get('n_models_verified') != 16 or \
       source_summary.get('absolute_spinor_normalization_certified') is not False or \
       source_summary.get('physical_eft1_workers_must_reject_this_as_certified_input') is not True:
        raise SourceKernelError('Source release certification/format changed; must re-audit')
    records=source_summary.get('models')
    if not isinstance(records,list) or len(records)!=16 or \
       {r.get('model') for r in records}!=set(MODELS):
        raise SourceKernelError('Expected exactly the sixteen known fermion-first boundaries')
    by_model={r['model']:r for r in records}
    selected=[a.model] if a.model else list(MODELS)
    numerical=a.coupling_root is not None
    prepared=[]
    rows=[]
    for model in selected:
        s=by_model[model]
        path=Path(s['path'])
        if path!=Path(tag(model))/'boundary.json':
            raise SourceKernelError(f'{model}: unsafe or unexpected release path')
        full=release/path
        boundary=read_json(full)
        if boundary.get('model')!=model or \
           len(boundary.get('sources',[]))!=s['PL_sources'] or \
           boundary.get('mu_F_gev')!=s['mu_F_gev']:
            raise SourceKernelError(f'{model}: source release manifest does not match boundary')
        sha=digest(full)
        contract=contracts(boundary,sha)
        prepared.append((out/tag(model)/'integrated_input_contract.json',contract))
        row={'model':model,'source_boundary_sha256':sha,
             'source_count':len(boundary['sources']),
             'status':'CrossWorkstreamContractOnly_MissingActualCouplingsAndScalarBackgrounds'}
        if numerical:
            coupling_path=a.coupling_root/tag(model)/'couplings.json'
            background_path=a.background_root/f'{tag(model)}.json'
            result=evaluate(full,coupling_path,background_path,
                            trajectory_points=a.trajectory_points)
            prepared.append((out/tag(model)/'conditional_eft1_joint_scalar.json',result))
            row['status']=result['status']
            row['quartic_scale_cancellation_residual']=result['quartic_muS_log_cancellation_relative_residual']
            row['full_SU2_quartic_residual']=result['full_SU2_quartic_residual']
            row['full_SU2_bubble_residual']=result['full_SU2_joint_bubble_residual']
        rows.append(row)
    summary={'format':FORMAT,'status':('ConditionalAllModelsProcessed_NotPhysicalC5' if numerical else
            'AllSourceBackedContractsPrepared_NoPhysicalInputs'),
            'source_release_summary_sha256':digest(release/'summary.json'),
            'model_count':len(rows),'models_evaluated_conditionally':len(rows) if numerical else 0,
            'physical_C5_computed':False,'GammaCC_factor':None,
            'seven_term_eft1_running_included':numerical,
            'one_joint_two_scalar_threshold_included':numerical,
            'source_boundary_contains_tree_LLSS_only':True,
            'hard_fermion_one_loop_matching_included':False,
            'complete_dim5_operator_inventory_certified':False,
            'renormalization_convention':'(1/2) C5 (L.H)(L.H); MSbar; real scalar quartic V=lambda*phi4/4!',
            'models':rows}
    prepared.append((out/'summary.json',summary))
    if a.check_only:
        print(json.dumps(summary,indent=2,allow_nan=False));return 0
    if out.exists():
        raise FileExistsError('Refusing existing output directory; choose a new --output-root')
    # Every model has been validated before writing anything.
    out.mkdir(parents=True,exist_ok=False)
    for path,obj in prepared:
        path.parent.mkdir(parents=True,exist_ok=True)
        with path.open('x',encoding='utf8') as stream:
            json.dump(obj,stream,indent=2,allow_nan=False)
            stream.write('\n')
    print(f'Wrote {len(rows)} '+('conditional joint EFT1 results' if numerical else
                                'cross-workstream input contracts')+f' to {out}')
    print('Physical C5 computed: false; signed GammaCC, hard F and full EFT1 physical inputs unverified')
    return 0


if __name__=='__main__':
    try:raise SystemExit(main())
    except (SourceKernelError,ValueError,KeyError,FileNotFoundError,
            FileExistsError,RuntimeError,ImportError) as error:
        print('BLOCKED:',error,file=sys.stderr)
        raise SystemExit(2)
