"""Run the full coupled EFT1 dimension-five Wilson tensor, MF -> MS.

Usage from repository root:
  python -m scripts.RunT3RoughStream2 --model A_-4 --c5-boundary zero
  python -m scripts.RunT3RoughStream2 --all --c5-boundary zero
  python -m scripts.RunT3RoughStream2 --all --probe-b-only

A numerical output is written only when canonical input fields and IDs exist.
Saved native Workstream B data can be inspected, but are never fake-converted.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

from Numerical.parallel_t3.wilson.basis import ALL_MODELS
from Numerical.parallel_t3.wilson.rough_stream2 import (
    INPUT_ROOT, OUTPUT_ROOT, ALLOWED_ADMS, run_one, probe_saved_b,
    slug,
)


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    choice=p.add_mutually_exclusive_group(required=True)
    choice.add_argument('--model',choices=ALL_MODELS)
    choice.add_argument('--all',action='store_true')
    p.add_argument('--project-root',type=Path,default=Path('.'))
    p.add_argument('--input-root',type=Path,default=INPUT_ROOT)
    p.add_argument('--output-root',type=Path,default=OUTPUT_ROOT)
    p.add_argument('--b-root',type=Path,default=Path(
        'output/parallel_t3/background/B_all16_qualitative_v1'))
    p.add_argument('--probe-b-only',action='store_true',
                   help='Read saved native B tensors and hashes; never claim full Wilson running')
    p.add_argument('--c5-boundary',choices=('require','zero'),default='require',
                   help='"zero" explicitly sets absent hard C5=0; otherwise require supplied C5')
    p.add_argument('--adm',choices=ALLOWED_ADMS,default='fonseca_eq22_six_term')
    p.add_argument('--rtol',type=float,default=1e-6)
    p.add_argument('--atol',type=float,default=1e-21)
    p.add_argument('--steps',type=int,default=6,
                   help='number of recorded log-scale intervals (DOP853 adaptive internally)')
    p.add_argument('--assume-formal-identity',action='store_true',
                   help='Explicit unverified Matchete scalar/spinor identity assumption; rough only')
    p.add_argument('--overwrite',action='store_true',
                   help='Replace only Stream2-owned output/model directories')
    args=p.parse_args(argv)
    if args.steps<1 or not 0<args.rtol<1 or not 0<args.atol<1:
        p.error('Need steps>=1 and 0<rtol,atol<1')
    if args.all:
        order=('A_-4',)+tuple(m for m in ALL_MODELS if m!='A_-4')
    else:order=(args.model,)
    root=args.project_root.resolve()
    out=root/args.output_root
    out.mkdir(parents=True,exist_ok=True)
    summary_path=out/('b_probe_summary.json' if args.probe_b_only else 'batch_summary.json')
    if summary_path.exists() and not args.overwrite:
        p.error(f'{summary_path} exists; use --overwrite for a new run')
    summary={
        'format':'t3-rough-stream2-batch-summary-v1',
        'selected_adm':args.adm,'solver':'DOP853',
        'model_order':list(order),
        'boundary_choice':args.c5_boundary,
        'background_couplings':'Stream1 canonical running scale samples only',
        'is_saved_B_probe_only':bool(args.probe_b_only),
        'results':{},'nonzero_C5_from_pure_LLSS_count':0,
        'successful_model_count':0,
    }
    for model in order:
        try:
            if args.probe_b_only:
                result=probe_saved_b(root,args.b_root,model)
                print(f'{model}: B read-only probe: {result["status"]}',flush=True)
            else:
                result=run_one(root,model,input_root=args.input_root,
                               output_root=args.output_root,
                               c5_boundary=args.c5_boundary,
                               adm_variant=args.adm,rtol=args.rtol,atol=args.atol,
                               steps=args.steps,overwrite=args.overwrite,
                               assume_formal_identity=args.assume_formal_identity)
                result={
                    'status':'computed_candidate',
                    'input_context_ids':{key:result[key] for key in (
                        'matching_context_id','matched_seed_id','lagrangian_source_id')},
                    'C5_mf_norm':result['C5_EFT1_mf_norm_gev_inverse'],
                    'C5_ms_norm':result['C5_EFT1_ms_norm_gev_inverse'],
                    'LLSS_mf_norm':result['LLSS_mf_norm_gev_inverse'],
                    'LLSS_ms_norm':result['LLSS_ms_norm_gev_inverse'],
                    'generated_from_zero_C5':result['C5_generated_from_zero_boundary'],
                    'generated_from_pure_LLSS':result['C5_generated_from_LLSS_only_boundary'],
                    'backend':result['solver']['backends'],
                    'output':str(root/args.output_root/slug(model))}
                if result['generated_from_pure_LLSS']:
                    summary['nonzero_C5_from_pure_LLSS_count']+=1
                summary['successful_model_count']+=1
                print(f'{model}: C5(MF)={result["C5_mf_norm"]:.6e} '
                      f'C5(MS)={result["C5_ms_norm"]:.6e} '
                      f'LLSS(MF)={result["LLSS_mf_norm"]:.6e}',flush=True)
        except Exception as e:
            # Batch must attempt all16; no synthetic replacement for failures.
            result={'status':'blocked','error_type':type(e).__name__,
                    'reason':str(e)}
            try:
                b=probe_saved_b(root,args.b_root,model)
                result['saved_Workstream_B_probe']=b
            except Exception as be:
                result['saved_Workstream_B_probe_error']=str(be)
            print(f'{model}: BLOCKED: {type(e).__name__}: {e}',file=sys.stderr,flush=True)
        summary['results'][model]=result
        # Persist attempts incrementally; a long batch crash still leaves evidence.
        summary_path.write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    if args.probe_b_only:
        return 0 if all(x['status']!='B_saved_input_missing'
                        for x in summary['results'].values()) else 2
    return 0 if summary['successful_model_count']==len(order) else 2


if __name__=='__main__':
    raise SystemExit(main())
