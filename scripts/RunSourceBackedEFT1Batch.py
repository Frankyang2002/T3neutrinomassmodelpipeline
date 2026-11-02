"""All-16 source-backed EFT1 audit, optional seven-term conditional diagnostics.

No quartic/Yukawa/gauge parameters are inferred from source or Stage22 files.
--coupling-root is OPTIONAL; without it only missing-input contracts are written.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0,str(ROOT))

from RGE.running.intermediate.SourceBackedEFT1Batch import (
    FORMAT, COUPLING_FORMAT, SOURCE_TO_MASTER, SourceKernelError, MODELS,
    tag, digest, read_json, required_couplings, evaluate, stage22_inventory,
)

DEFAULT_RELEASE='output/general_observable/lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f'
DEFAULT_OUTPUT='output/general_observable/lo_eft1_mixing/source_backed_batch_v1'


def main(argv=None) -> int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--release-root',type=Path,default=Path(DEFAULT_RELEASE))
    p.add_argument('--output-root',type=Path,default=Path(DEFAULT_OUTPUT))
    p.add_argument('--coupling-root',type=Path,default=None,
                   help='External explicit, diagnostic-only model JSON files <tag>/couplings.json')
    p.add_argument('--stage22-diagnostics',type=Path,default=None,
                   help='Optional read-only Stage22 symbolic geometry; never treated as couplings')
    p.add_argument('--model',choices=MODELS,default=None)
    p.add_argument('--trajectory-points',type=int,default=21)
    p.add_argument('--check-only',action='store_true')
    args=p.parse_args(argv)
    out=args.output_root.resolve()
    allowed=(ROOT/'output/general_observable/lo_eft1_mixing').resolve()
    if out==allowed or allowed not in out.parents:
        raise SourceKernelError(f'Outputs restricted to new subdirectory beneath {allowed}')
    source_root=args.release_root.resolve()
    _summary=source_root/'summary.json'
    src_summary=read_json(_summary)
    if src_summary.get('format')!='t3-lo-fermion-boundary-source-resolved-provisional-v1' or \
        src_summary.get('n_models_verified')!=16 or \
        src_summary.get('absolute_spinor_normalization_certified') is not False or \
        src_summary.get('physical_eft1_workers_must_reject_this_as_certified_input') is not True:
        raise SourceKernelError('Unexpected release source certification/format, refusing processing')
    model_rows=src_summary.get('models')
    if not isinstance(model_rows,list) or {r.get('model') for r in model_rows}!=set(MODELS):
        raise SourceKernelError('Source release must document the 16 expected model names')
    records={r['model']:r for r in model_rows}
    selected=(args.model,) if args.model else MODELS
    stage22=None
    if args.stage22_diagnostics is not None:
        stage22=stage22_inventory(read_json(args.stage22_diagnostics))
    outputs=[]
    rows=[]
    for model in selected:
        row=records[model]; rel=Path(row['path'])
        if rel.is_absolute() or '..' in rel.parts or rel!=Path(tag(model))/'boundary.json':
            raise SourceKernelError(f'{model}: unexpected source path')
        path=source_root/rel
        boundary=read_json(path); sha=digest(path)
        if boundary.get('model')!=model or boundary.get('status')!='NumericalLLSSBySourceVerified_SpinorNormalizationOpen':
            raise SourceKernelError(f'{model}: unexpected source boundary status')
        if len(boundary.get('sources',[]))!=row['PL_sources'] or \
            boundary.get('mu_F_gev')!=row['mu_F_gev']:
            raise SourceKernelError(f'{model}: source release summary mismatch')
        contract=required_couplings(boundary,sha)
        if stage22 is not None:
            contract['stage22_extra_quartic_geometry']=stage22[model]
        files=[(out/tag(model)/'required_couplings_contract.json',contract)]
        result_status='ContractOnly_NoPhysicalCouplings'
        if args.coupling_root is not None:
            input_path=args.coupling_root/tag(model)/'couplings.json'
            coupling=read_json(input_path)
            result=evaluate(boundary,sha,coupling,trajectory_points=args.trajectory_points)
            files.append((out/tag(model)/'seven_term_diagnostic.json',result))
            result_status=result['status']
        outputs.extend(files)
        rows.append({'model':model,'source_count':len(boundary['sources']),
                     'source_boundary_sha256':sha,
                     'status':result_status,
                     'n_real_scalars':boundary['real_scalar_dimension']})
    summary={
        'format':FORMAT,'status':'SourceBackedNonPhysicalEFT1Batch',
        'source_release_summary_sha256':digest(_summary),'model_count':len(rows),
        'coupling_input_format':COUPLING_FORMAT if args.coupling_root else None,
        'stage22_geometry_audited':stage22 is not None,
        'stage22_is_physical_quartic_input':False,
        'per_NGamma_conditional_diagnostics_computed':args.coupling_root is not None,
        'physical_C5_computed':False,
        'source_to_master_status':SOURCE_TO_MASTER,
        'hard_fermion_one_loop_included':False,'joint_scalar_matching_included':False,
        'models':rows,
    }
    outputs.append((out/'summary.json',summary))
    if args.check_only:
        print(json.dumps(summary,indent=2));return 0
    if any(path.exists() for path,_ in outputs):
        raise FileExistsError('Output exists; refusing to overwrite a previous batch')
    for path,obj in outputs:
        path.parent.mkdir(parents=True,exist_ok=True)
        with path.open('x',encoding='utf8') as f:
            json.dump(obj,f,indent=2,allow_nan=False)
            f.write('\n')
    print(f'Wrote {len(rows)} source-backed model contracts to {out}')
    print('Conditional seven-term diagnostics:',str(args.coupling_root is not None).lower())
    print('Physical C5 computed: false (GammaCC and physical quartics unresolved)')
    return 0


if __name__=='__main__':
    try:
        raise SystemExit(main())
    except (SourceKernelError,ValueError,KeyError,FileNotFoundError,FileExistsError) as exc:
        print(f'BLOCKED: {exc}',file=sys.stderr)
        raise SystemExit(2)
