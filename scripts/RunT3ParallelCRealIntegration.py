"""Read-only A/B preflight; optional C-only simultaneous matching and SMEFT.

Default behaviour is PREVIEW, not synthetic generation. Point --wilson-root and
--background-root to the exact A/B endpoint run directories, then pass --run.
Only C-owned result paths are written. No A/B file is ever rewritten.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from Numerical.parallel_t3.threshold.models import ALL_MODELS, normalize_model
from Numerical.parallel_t3.threshold.preflight import PREFLIGHT_FORMAT, preflight_pair
from Numerical.parallel_t3.threshold.pipeline import run_case


def _slug(model: str) -> str:
    return model.replace('_-', '_m').replace('_+', '_p')


def _manifest(root: Path, model: str, name: str) -> Path:
    """Find a single unambiguous manifest in the explicit run root only."""
    stem = 'wilson' if name == 'A' else 'background'
    tag = _slug(model)
    names = (f'{stem}.json', f'{stem}_at_MS.json') if name == 'A' else (
        'background.json', 'background_at_MS.json')
    candidates = [root/parent/name for parent in (tag, f'models/{tag}') for name in names]
    existing = [p for p in candidates if p.is_file()]
    if len(existing) > 1:
        raise ValueError(f'Ambiguous {name} manifests for {model}; specify a single run root: {existing}')
    return existing[0] if existing else candidates[0]


def _json_write(path: Path, document: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(document, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def _norms(path: Path) -> tuple[np.ndarray,np.ndarray]:
    with np.load(path, allow_pickle=False) as z:
        return np.array(z['C5_SMEFT_MS']),np.array(z['C5_scalar_bubble_MS'])


def _comparison(output_root: Path) -> dict:
    x,x_b = _norms(output_root/'strict_one_loop'/'C5_smeft_at_MS.npz')
    y,y_b = _norms(output_root/'rge_improved'/'C5_smeft_at_MS.npz')
    def metrics(a, b):
        return {'difference_frobenius_gev_inverse': float(np.linalg.norm(b-a)),
                'relative_to_larger_norm': float(np.linalg.norm(b-a)/max(
                    np.linalg.norm(a),np.linalg.norm(b),1e-100))}
    return {'format': 't3-parallel-c-loop-counting-comparison-v2',
            'C5_SMEFT_MS_rge_improved_minus_strict': metrics(x,y),
            'scalar_bubble_rge_improved_minus_strict': metrics(x_b,y_b),
            'fixed_order_interpretation': ('The difference in LLSS bubble contractions uses '
                'one-loop-evolved rather than order-zero Wilson coefficients and is formally '
                'higher order in the fixed-order expansion. Other effects are not certified.'),
            'no_direct_LLSS_to_C5_log_added_twice': True}


def parser() -> argparse.ArgumentParser:
    p=argparse.ArgumentParser(description='C real A/B endpoint preflight and simultaneous threshold')
    p.add_argument('--wilson-root',type=Path,default=Path('output/parallel_t3/wilson'),
                   help='Exact A output run directory, containing <model_slug>/wilson.json')
    p.add_argument('--background-root',type=Path,default=Path('output/parallel_t3/background'),
                   help='Exact B output run directory, containing <model_slug>/background.json')
    p.add_argument('--wilson',type=Path,help='One explicit A wilson manifest')
    p.add_argument('--background',type=Path,help='One explicit B background manifest')
    p.add_argument('--models',nargs='*',default=None,
                   help='Model tags (default all 16); e.g. A_-4 B_-1')
    p.add_argument('--output',type=Path,default=Path('output/parallel_t3/threshold/real_ab_v2'))
    p.add_argument('--mode',choices=('strict_one_loop','rge_improved','both'),default='both')
    p.add_argument('--run',action='store_true',help='Actually match and evolve after preflight')
    p.add_argument('--mu-low',type=float,default=173.)
    p.add_argument('--vev',type=float,default=246.22)
    p.add_argument('--match-only',action='store_true',help='Skip SMEFT RGE explicitly')
    p.add_argument('--allow-projection-residual',action='store_true',
                   help='Diagnostic only; nonzero SU2 residual is not pure Weinberg')
    p.add_argument('--permit-synthetic',action='store_true',
                   help='Test this real-integration CLI against explicit synthetic fixture inputs')
    p.add_argument('--overwrite',action='store_true',help='Allow replacing existing C result files')
    return p


def main(argv=None) -> int:
    args=parser().parse_args(argv)
    if bool(args.wilson) != bool(args.background):
        parser().error('--wilson and --background must be provided together')
    if args.wilson and args.models and len(args.models)!=1:
        parser().error('Explicit A/B manifests require exactly one model (or omit --models)')
    if args.wilson:
        if args.models:
            models=[normalize_model(args.models[0])]
        else:
            try:
                models=[normalize_model(json.loads(args.wilson.read_text(encoding='utf-8-sig'))['model'])]
            except (OSError,ValueError,KeyError) as exc:
                parser().error(f'Explicit A manifest missing/invalid model: {exc}')
    else:
        names=args.models or list(ALL_MODELS)
        models=[normalize_model(m) for m in names]
    if len(models) != len(set(models)):
        parser().error('Duplicate model requested')
    modes=('strict_one_loop','rge_improved') if args.mode=='both' else (args.mode,)
    report={'format':'t3-parallel-c-real-integration-summary-v2',
            'mode':args.mode,'run_requested':bool(args.run),
            'models_requested':models, 'A_root':str(args.wilson_root),
            'B_root':str(args.background_root),
            'synthetic_allowed_only_for_explicit_test':bool(args.permit_synthetic),
            'model_results':{},'counts':{}}
    any_blocked=False
    for model in models:
        entry={'modes':{}, 'status':'blocked'}
        report['model_results'][model]=entry
        try:
            w,b=(args.wilson,args.background) if args.wilson else (
                _manifest(args.wilson_root,model,'A'),
                _manifest(args.background_root,model,'B'))
        except (OSError,ValueError) as exc:
            entry['reason']=f'Manifest discovery: {exc}'
            any_blocked=True
            print(model, 'BLOCKED', entry['reason'],flush=True)
            continue
        for mode in modes:
            pf=preflight_pair(w,b,mode=mode,mu_low=args.mu_low,
                              allow_projection_residual=args.allow_projection_residual,
                              expected_model=model,permit_synthetic=args.permit_synthetic)
            out=args.output/_slug(model)
            if args.mode=='both':out=out/mode
            entry['modes'][mode]={'preflight': pf,
                                  'preflight_status': pf['status']}
            if pf['status']!='ready_conditional':
                any_blocked=True
                print(model,mode,'BLOCKED:',pf['blockers'][0] if pf['blockers'] else 'unknown',flush=True)
                continue
            if not args.run:
                print(model,mode,'READY_CONDITIONAL (preflight only)',flush=True)
                continue
            try:
                result=run_case(wilson_manifest=w,background_manifest=b,
                                output_dir=out,mu_low=args.mu_low,vev=args.vev,
                                mode=mode,allow_projection_residual=args.allow_projection_residual,
                                match_only=args.match_only,overwrite=args.overwrite)
                entry['modes'][mode]['result_manifest']=str(out/'manifest.json')
                entry['modes'][mode]['result_status']=result['status']
                print(model,mode,'MATCHED_CONDITIONAL',out/'manifest.json',flush=True)
            except Exception as exc:
                any_blocked=True
                entry['modes'][mode]['result_status']='blocked'
                entry['modes'][mode]['run_error']=f'{type(exc).__name__}: {exc}'
                print(model,mode,'BLOCKED running:',exc,flush=True)
        if (args.mode=='both' and args.run and
                all(entry['modes'].get(mode,{}).get('result_status')=='NumericalCandidateNotCertified'
                    for mode in modes)):
            try:
                entry['mode_comparison']=_comparison(args.output/_slug(model))
                _json_write(args.output/_slug(model)/'counting_comparison.json',
                            entry['mode_comparison'])
            except (OSError,ValueError,KeyError) as exc:
                any_blocked=True
                entry['mode_comparison_error']=f'{type(exc).__name__}: {exc}'
        passed=all(entry['modes'].get(mode,{}).get(
            'result_status' if args.run else 'preflight_status')==(
            'NumericalCandidateNotCertified' if args.run else 'ready_conditional')
            for mode in modes)
        entry['status']='success_conditional' if passed else 'blocked'
    report['counts']={'model_count':len(models),
                      'all_requested_modes_ready_or_completed':sum(
                          data['status']=='success_conditional' for data in report['model_results'].values()),
                      'blocked':sum(data['status']=='blocked' for data in report['model_results'].values())}
    _json_write(args.output/'real_integration_summary.json',report)
    print('SUMMARY',report['counts'],'->',args.output/'real_integration_summary.json',flush=True)
    return 1 if any_blocked else 0


if __name__=='__main__':
    sys.exit(main())
