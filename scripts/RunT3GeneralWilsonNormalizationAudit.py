"""Read-only 16-model LLSS chirality/normalization audit.

Requires corrected B/C PL tensors, not the earlier Proj[1] selection.
Does not invent the undetermined GammaCC -> Weyl bilinear sign/factor.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

from Numerical.running.GeneralT3WilsonEvolution import (
    MODELS, tag, verify_batch, read_json, save_json, digest,
)
from Numerical.running.GeneralT3WilsonNormalizationAudit import (
    FORMAT, STATUS, audit_model, master_evidence,
)
from scripts.RunT3GeneralWilsonPolynomialBatch import verified_source_path

ROOT=Path(__file__).resolve().parents[1]
DEFAULT_BOUNDARY=Path('output/general_observable/ordered_wilson_tensor_v4')
DEFAULT_POLYNOMIAL=Path('output/general_observable/general_wilson_polynomial_batch_v2')
DEFAULT_OUTPUT=Path('output/general_observable/general_wilson_normalization_audit_v1')


def absolute(path: Path) -> Path:
    return path.resolve() if path.is_absolute() else (ROOT/path).resolve()


def main(argv: list[str] | None = None) -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    sel=parser.add_mutually_exclusive_group(required=True)
    sel.add_argument('--all',action='store_true')
    sel.add_argument('--models',help='Comma-separated subset, still verifies all-16 boundary')
    parser.add_argument('--boundary-root',type=Path,default=DEFAULT_BOUNDARY)
    parser.add_argument('--polynomial-root',type=Path,default=DEFAULT_POLYNOMIAL)
    parser.add_argument('--output-root',type=Path,default=DEFAULT_OUTPUT)
    args=parser.parse_args(argv)
    models=MODELS if args.all else tuple(x.strip() for x in args.models.split(',') if x.strip())
    if not models or len(set(models))!=len(models) or not set(models)<=set(MODELS):
        parser.error('Unknown/duplicate model selection')
    boundary=absolute(args.boundary_root)
    polynomials=absolute(args.polynomial_root)
    output=absolute(args.output_root)
    if len({boundary,polynomials,output})!=3 or any(
        a in b.parents for a in (boundary,polynomials,output)
        for b in (boundary,polynomials,output) if a!=b):
        parser.error('Input/output directories must be distinct and non-nested')
    try:
        checked=verify_batch(boundary)
        conventions=master_evidence(ROOT)
    except (OSError,ValueError,KeyError,RuntimeError,ImportError) as exc:
        print(f'BLOCKED: {type(exc).__name__}: {exc}',file=sys.stderr)
        return 2
    reports=[]
    for model in models:
        try:
            ordered_path=boundary/tag(model)/'ordered_tensor.json'
            polynomial_path=polynomials/tag(model)/'polynomial_bridge.json'
            ordered=read_json(ordered_path)
            bridge=read_json(polynomial_path)
            source_path=verified_source_path(model,ordered,boundary)
            source=read_json(source_path)
            if bridge.get('ordered_boundary_sha256')!=digest(ordered_path):
                raise ValueError('Polynomial bridge was built from another tensor version')
            if bridge.get('operator_source_sha256')!=digest(source_path):
                raise ValueError('Polynomial bridge source checksum mismatch')
            report=audit_model(model,ordered,source,bridge)
            report.update({
                'ordered_tensor_sha256':digest(ordered_path),
                'polynomial_bridge_sha256':digest(polynomial_path),
                'matchete_operator_source_sha256':digest(source_path),
            })
            dest=output/tag(model)/'normalization_audit.json'
            save_json(dest,report)
            row={'model':model,'status':report['status'],
                 'PL_terms':report['native_matchete_chirality']['native_PL_matched_term_count'],
                 'PR_terms_ignored':report['native_matchete_chirality']['native_PR_conjugate_term_count'],
                 'independent_coordinates':report['master_polynomial_coordinate_check']['independent_symmetric_coordinate_count'],
                 'report':str(dest)}
        except (OSError,ValueError,TypeError,KeyError,RuntimeError,IndexError,ImportError) as exc:
            row={'model':model,'status':'Blocked','reason':f'{type(exc).__name__}: {exc}'}
            save_json(output/tag(model)/'error.json',row)
        reports.append(row)
        print(f'{model}: {row["status"]}' +
              ('; '+row['reason'] if 'reason' in row else ''),flush=True)
        save_json(output/'summary.json',{
            'format':FORMAT,
            'status':'AllRequestedPLAndPolynomialStructureVerified_SpinorNormalizationOpen'
                      if len(reports)==len(models) and all(x['status']==STATUS for x in reports)
                      else 'PartialOrBlocked',
            'models_requested':len(models),
            'models_structure_verified':sum(x['status']==STATUS for x in reports),
            'models_blocked':sum(x['status']=='Blocked' for x in reports),
            'full_16_model_batch':args.all,
            'master_and_matchete_source_conventions':conventions,
            'absolute_spinor_normalization_certified':False,
            'physical_eft1_wilson_running_computed':False,
            'results':reports,
        })
    good=sum(x['status']==STATUS for x in reports)
    print(f'PL/polynomial normalization structure: {good}/{len(models)}; '
          'absolute GammaCC/Weyl normalization: NOT CERTIFIED; EFT1 running: NOT CERTIFIED')
    return 0 if good==len(models) else 2

if __name__=='__main__':
    raise SystemExit(main())
