"""Workstream A standalone smoke and strict B->A Wilson transport CLI.

Run from repository root: python -m scripts.RunT3ParallelA --synthetic A_-4
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
from Numerical.parallel_t3.wilson.basis import ALL_MODELS
from Numerical.parallel_t3.wilson.workflow import (
    BlockedInput, synthetic_case, run_wilson, write_results, run_from_files)


def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__)
    options=ap.add_mutually_exclusive_group(required=True)
    options.add_argument('--synthetic',choices=ALL_MODELS)
    options.add_argument('--synthetic-all16',action='store_true')
    options.add_argument('--model',choices=ALL_MODELS)
    ap.add_argument('--root',type=Path,default=Path('.'))
    ap.add_argument('--background-npz',type=Path)
    ap.add_argument('--background-json',type=Path)
    ap.add_argument('--boundary-npz',type=Path)
    ap.add_argument('--boundary-json',type=Path)
    ap.add_argument('--out',type=Path,default=Path('output/parallel_t3/wilson'))
    ap.add_argument('--adm-variant',choices=('fonseca_eq22_six_term',
                     'ward_covariant_candidate','source_v1_unstarred'),
                     default='fonseca_eq22_six_term')
    ap.add_argument('--solver',choices=('RK2','DOP853'),default='RK2')
    ap.add_argument('--steps',type=int,default=3)
    ap.add_argument('--assume-formal-identity',action='store_true',
                    help='FOR QUALITATIVE TESTS ONLY: map unknown scalar/spinor normalization to identity/1')
    args=ap.parse_args(argv)
    models=ALL_MODELS if args.synthetic_all16 else (args.synthetic or args.model,)
    status={}
    for model in models:
        folder=args.out / model.replace('+','p').replace('-','m')
        try:
            if args.synthetic or args.synthetic_all16:
                basis,initial,bg,seed=synthetic_case(model)
                result=run_wilson(basis,initial,bg,seed['MF_gev'],seed['MS_gev'],
                                  adm_variant=args.adm_variant,solver=args.solver,steps=args.steps)
                report=write_results(result,basis,folder,
                     {'synthetic_fixture':True, 'synthetic_parameters':seed,
                      'basis_not_Matchete_RGBeta_certified':True,
                      'corresponding_native_F_formal_boundary_used':False})
            else:
                required=(args.background_npz,args.background_json,args.boundary_npz,args.boundary_json)
                if not all(required):
                    raise BlockedInput('Blocked_RequireFourBInputPaths')
                report=run_from_files(args.root,model,*required,folder,
                         adm_variant=args.adm_variant,solver=args.solver,steps=args.steps,
                         assume_formal_identity=args.assume_formal_identity)
            print(f"{model}: {report['status']}, C5_EFT1(MS) norm = {report['C5_at_MS_norm']:.6e}; {folder}")
            status[model]={'status':report['status'],
                           'C5_EFT1_MS_frobenius':report['C5_at_MS_norm']}
        except (ValueError,RuntimeError,FileExistsError) as exc:
            msg=str(exc);status[model]={'status':'Blocked','reason':msg}
            print(f'{model}: Blocked: {msg}')
            if not args.synthetic_all16:raise
    if args.synthetic_all16:
        args.out.mkdir(parents=True,exist_ok=True)
        p=args.out/'all16_qualitative_summary.json'
        if p.exists():raise FileExistsError('Refusing to overwrite '+str(p))
        p.write_text(json.dumps(status,indent=2)+'\n',encoding='utf8')
    if any(k['status']=='Blocked' for k in status.values()):return 2
    return 0


if __name__=='__main__':raise SystemExit(main())
