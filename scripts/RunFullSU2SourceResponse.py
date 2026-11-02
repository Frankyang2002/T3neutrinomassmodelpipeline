"""CLI for conditional all-16 source-resolved EFT1 SU(2) kernel and Ward audit."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

PROJECT=Path(__file__).resolve().parents[1]
if str(PROJECT) not in sys.path:
    sys.path.insert(0,str(PROJECT))

from RGE.running.intermediate.FullSU2SourceResponse import (
    digest, load_json, run_diagnostic,
)
from RGE.running.intermediate.ProvisionalSourceKernel import (
    MODELS, SourceKernelError, tag,
)

DEFAULT_RELEASE='output/general_observable/lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f'
DEFAULT_OUTPUT='output/general_observable/lo_eft1_mixing/source_su2_response_v1'


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--release-root',type=Path,default=Path(DEFAULT_RELEASE))
    parser.add_argument('--output-root',type=Path,default=Path(DEFAULT_OUTPUT))
    parser.add_argument('--model',choices=MODELS)
    parser.add_argument('--quartic-root',type=Path,
                        help='OPTIONAL diagnostic quartic inputs, <tag>/quartic.json')
    parser.add_argument('--mu-S-gev',type=float,default=None)
    parser.add_argument('--check-only',action='store_true',help='Validate, print summary, write nothing')
    args=parser.parse_args(argv)
    release=args.release_root.resolve()
    out=args.output_root.resolve()
    allowed=(PROJECT/'output/general_observable/lo_eft1_mixing').resolve()
    if out!=allowed and allowed not in out.parents:
        raise SourceKernelError(f'Output root must be inside {allowed}')
    if not release.is_dir():
        raise FileNotFoundError(f'Fermion source release missing: {release}')
    summary_path=release/'summary.json'
    summary=load_json(summary_path)
    if summary.get('absolute_spinor_normalization_certified') is not False or \
       summary.get('physical_eft1_workers_must_reject_this_as_certified_input') is not True:
        raise SourceKernelError('Unexpected release certificate; inspect before processing')
    rows=summary.get('models')
    if not isinstance(rows,list):
        raise SourceKernelError('Missing release model inventory')
    catalog={r['model']:r for r in rows}
    if args.model is None and set(catalog)!=set(MODELS):
        raise SourceKernelError('Release must contain exactly 16 models')
    selected=(args.model,) if args.model else MODELS
    results=[]
    report_rows=[]
    for model in selected:
        row=catalog.get(model)
        if row is None or row.get('status')!='NumericalLLSSBySourceVerified_SpinorNormalizationOpen':
            raise SourceKernelError(f'{model}: missing/changed release row')
        relative=Path(row['path'])
        if relative.is_absolute() or '..' in relative.parts or relative != Path(tag(model))/'boundary.json':
            raise SourceKernelError(f'{model}: unsafe release path')
        source_path=release/relative
        boundary=load_json(source_path)
        if boundary.get('model')!=model or len(boundary.get('sources',[]))!=row.get('PL_sources'):
            raise SourceKernelError(f'{model}: source count/model mismatch')
        quartic=None
        if args.quartic_root is not None:
            qpath=args.quartic_root/tag(model)/'quartic.json'
            quartic=load_json(qpath)
        value=run_diagnostic(boundary, source_sha256=digest(source_path),
                             quartic=quartic,mu_s_gev=args.mu_S_gev)
        results.append((model,value))
        rec={'model':model, 'source_count':value['source_count'],
             'source_boundary_sha256':value['source_boundary_sha256'],
             'status':value['status']}
        if 'total' in value:
            rec['conditional_SU2_Ward_pass']=value['total']['SU2_Ward_tensor_algebra_pass']
            rec['conditional_SU2_Ward_relative_residual']=value['total']['weinberg_residual_relative']
        report_rows.append(rec)
    batch={'format':'t3-eft1-full-su2-source-response-batch-v1',
           'status':'FullSU2SourceResponsesOnly_NotPhysicalC5',
           'model_count':len(results),
           'fermion_release_summary_sha256':digest(summary_path),
           'quartic_input_kind':'none' if args.quartic_root is None else 'diagnostic_only',
           'physical_C5_computed':False,
           'source_spinor_normalization_certified':False,
           'models':report_rows}
    if args.check_only:
        print(json.dumps(batch,indent=2));return 0
    files=[out/tag(m)/'full_su2_response.json' for m,_ in results]+[out/'summary.json']
    if any(p.exists() for p in files):
        raise FileExistsError('Refusing overwrite of source SU2 response output')
    for model,value in results:
        path=out/tag(model)/'full_su2_response.json'
        path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    (out/'summary.json').write_text(json.dumps(batch,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(f'Wrote {len(results)} source-resolved full-SU2 responses to {out}')
    print('Physical C5 computed: false; spinor/quartic normalization remains uncertified')
    return 0


if __name__=='__main__':
    try:
        sys.exit(main())
    except (SourceKernelError,FileExistsError,FileNotFoundError,ValueError) as e:
        print(f'BLOCKED: {e}',file=sys.stderr)
        sys.exit(2)
