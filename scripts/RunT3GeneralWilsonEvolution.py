"""Inspect all 16 ordered T3 Wilson boundaries; conditionally evolve exact gamma.

By default, only preflight is performed: it never manufactures a gamma matrix.
Use --mode evolve --gamma-root <directory> only with independently supplied
full gamma samples in the exact ordered-source basis.

Example:
 python -X utf8 -m scripts.RunT3GeneralWilsonEvolution --all --mode preflight
"""
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
import sys

from Numerical.running.GeneralT3WilsonEvolution import (
    MODELS, boundary_coordinates, aggregate_source_coordinates, digest, evolve_given_gamma,
    read_json, save_json, tag, verify_batch,
)

ROOT=Path(__file__).resolve().parents[1]
DEFAULT_INPUT=Path('output/general_observable/ordered_wilson_tensor_v3')
DEFAULT_OUTPUT=Path('output/general_observable/wilson_eft1_evolution_preflight_v1')


def resolve(path:Path) -> Path:
    return path.resolve() if path.is_absolute() else (ROOT/path).resolve()


def run(args:argparse.Namespace) -> int:
    inp=resolve(args.boundary_root)
    out=resolve(args.output_root)
    inspected=verify_batch(inp)
    selection=MODELS if args.all else (args.model,)
    rows=[]
    for model in selection:
        info=inspected[model]
        row={k:v for k,v in info.items() if k not in ('coordinates','initial')}
        raw=read_json(Path(info['tensor_file']))
        source_labels,source_values,_=boundary_coordinates(raw,model)
        combined_labels,combined_values=aggregate_source_coordinates(source_labels,source_values)
        combined_path=out/tag(model)/'summed_ordered_boundary.json'
        save_json(combined_path, {
            'format':'t3-summed-ordered-source-wilson-boundary-v1',
            'model':model,'ordered_tensor_sha256':info['tensor_sha256'],
            'mu_F_gev':info['mu_F_gev'],
            'coordinate_basis':'ordered_Matchete_real_components_with_flavour',
            'coordinates':combined_labels,
            'real':combined_values.real.tolist(),
            'imag':combined_values.imag.tolist(),
            'uncombined_source_coordinate_count':len(source_labels),
            'combined_coordinate_count':len(combined_labels),
            'canonical_C_ijab_normalization_certified':False,
            'physical_eft1_wilson_running_certified':False,
        })
        row['summed_ordered_boundary_file']=str(combined_path)
        row['summed_ordered_coordinate_count']=len(combined_labels)
        if args.mode=='preflight':
            row['status']='BoundaryVerifiedGammaNotSupplied'
            row['full_gamma_available']=False
            if args.emit_coordinates:
                save_json(out/tag(model)/'gamma_basis_requirements.json',{
                    'format':'t3-gamma-coordinate-requirements-v1',
                    'model':model,
                    'ordered_tensor_sha256':info['tensor_sha256'],
                    'mu_F_gev':info['mu_F_gev'],
                    'source_coordinates':combined_labels,
                    'source_term_coordinates_combined_by_exact_order':True,
                    'closure_coordinates_must_be_supplied_by_gamma_source':True,
                    'gamma_matrix_derived':False,
                })
        else:
            source=resolve(args.gamma_root)/tag(model)/'gamma.json'
            try:
                if not source.is_file():
                    raise FileNotFoundError(f'Gamma source not found: {source}')
                g=read_json(source)
                raw=read_json(Path(info['tensor_file']))
                result=evolve_given_gamma(raw,model,info['tensor_sha256'],g,
                                         npoints=args.points,max_dim=args.max_gamma_dim)
                save_json(out/tag(model)/'conditional_gamma_running.json',result)
                row.update({'status':result['status'],
                            'gamma_sha256':digest(source),
                            'closure_zero_initial_coordinates':result['closure_zero_initial_coordinates'],
                            'relative_change':result['relative_change'],
                            'trajectory_file':str(out/tag(model)/'conditional_gamma_running.json')})
            except (FileNotFoundError,KeyError,ValueError,RuntimeError,TypeError) as exc:
                row.update({'status':'GammaBlocked','reason':f'{type(exc).__name__}: {exc}'})
        rows.append(row)
        print(f"{model}: {row['status']}"+(f"; {row['reason']}" if 'reason' in row else ''),flush=True)
    summary={'format':'t3-all16-ordered-wilson-gamma-handoff-v1',
             'mode':args.mode,
             'status':('VerifiedBoundaryPreflightNoGamma' if args.mode=='preflight'
                       else 'ConditionalGammaTransport' if all(x['status']=='ConditionalExplicitGammaNumericalTransport' for x in rows)
                       else 'ConditionalGammaIncomplete'),
             'selected_model_count':len(selection),
             'verified_complete_16_model_boundary':True,
             'conditional_gamma_evolved_count':sum(r['status']=='ConditionalExplicitGammaNumericalTransport' for r in rows),
             'blocked_gamma_count':sum(r['status']=='GammaBlocked' for r in rows),
             'physical_eft1_wilson_running_certified':False,
             'canonical_wilson_normalization_certified':False,
             'scalar_threshold_c5_computed':False,
             'results':rows}
    save_json(out/'summary.json',summary)
    print(f"Verified matching tensors: 16/16; conditional gamma evolved: {summary['conditional_gamma_evolved_count']}/{len(selection)}; gamma blocked: {summary['blocked_gamma_count']}; physical EFT1 running: NOT CERTIFIED")
    return 1 if summary['blocked_gamma_count'] else 0


def main(argv:list[str]|None=None) -> int:
    ap=argparse.ArgumentParser(description=__doc__)
    choice=ap.add_mutually_exclusive_group(required=True)
    choice.add_argument('--all',action='store_true')
    choice.add_argument('--model',choices=MODELS)
    ap.add_argument('--mode',choices=('preflight','evolve'),default='preflight')
    ap.add_argument('--boundary-root',type=Path,default=DEFAULT_INPUT)
    ap.add_argument('--output-root',type=Path,default=DEFAULT_OUTPUT)
    ap.add_argument('--gamma-root',type=Path)
    ap.add_argument('--points',type=int,default=41)
    ap.add_argument('--max-gamma-dim',type=int,default=2500)
    ap.add_argument('--emit-coordinates',action='store_true')
    args=ap.parse_args(argv)
    if args.mode=='evolve' and args.gamma_root is None:
        ap.error('--mode evolve requires --gamma-root; no gamma is inferred')
    try:return run(args)
    except (OSError,ValueError,KeyError,TypeError) as exc:
        print('BLOCKED: '+f'{type(exc).__name__}: {exc}',file=sys.stderr)
        return 2

if __name__=='__main__':raise SystemExit(main())
