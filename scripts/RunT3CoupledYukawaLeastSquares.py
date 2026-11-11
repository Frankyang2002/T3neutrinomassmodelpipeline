"""Bounded, resumable local least-squares refinement of saved coupled-EFT1 Yukawas.

Project path: scripts/RunT3CoupledYukawaLeastSquares.py
Requires scripts/RunT3CoupledYukawaOvernight.py (fingerprint-compatible version).
Uses the exact existing case-1111 EFT1 evaluator; no production source changes.
The five whitened oscillation residuals, not sqrt(chi2), are minimised.
This is exploratory and uses the cached/provisional matching approximations.
"""
from __future__ import annotations

import argparse
import contextlib
import hashlib
import json
import math
import shutil
import tempfile
import time
import traceback
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from scripts import RunT3CoupledYukawaOvernight as baseline

ROOT = baseline.ROOT
BASE = baseline.BASE
DEFAULT_OUT = BASE / 'coupled_yukawa_least_squares'
MIXED = baseline.MIXED
SPECIAL = baseline.SPECIAL


class BudgetReached(Exception):
    pass


def read(path):
    return baseline.read_json(path)


def save(path, data):
    return baseline.write_json(path, data)


def signature_file(path):
    return baseline.sha(path)


def run_one(model, args, source_root, output_root, registry, release, channels, canonical, target):
    from scripts import RunT3GeneralAll16FromReference as general

    tag = general.tag(model)
    source = source_root / tag
    previous = BASE / 'coupled_yukawa_overnight' / args.from_run / tag / 'best.json'
    if not previous.is_file():
        raise FileNotFoundError(f'{model}: missing overnight checkpoint {previous}')
    old = read(previous)
    if old.get('model') != model or old.get('status') != 'Computed':
        raise ValueError(f'{model}: incompatible overnight checkpoint')
    x_seed = np.asarray(old['yukawas_flat'], float)
    if x_seed.shape != (18,) or not np.isfinite(x_seed).all():
        raise ValueError(f'{model}: invalid saved Yukawa vector')
    if np.max(np.abs(x_seed)) > args.bound + 1e-10:
        raise ValueError(f'{model}: saved Yukawa exceeds --bound={args.bound}')
    config = read(source / 'fitted_config.json')
    overrides = read(source / 'general_uv_overrides.json') if model not in MIXED else None
    source_hash = {
        'source_config': signature_file(source/'fitted_config.json'),
        'source_overrides': signature_file(source/'general_uv_overrides.json') if model not in MIXED else None,
        'overnight_best': signature_file(previous),
        'target': signature_file(args.target),
    }
    run = output_root / tag
    checkpoint = run/'best.json'
    config_path = run/'best_fitted_config.json'
    overrides_path = run/'best_general_uv_overrides.json'
    history_path = run/'evaluations.jsonl'
    manifest = {
        'format':'t3-coupled-yukawa-least-squares-v1',
        'model': model,
        'from_run': args.from_run,
        'source_run':args.source_run,
        'canonical_run_id':canonical,
        'bound':args.bound,
        'active_counts':list(args.active_counts),
        'max_evals':args.max_evals,
        'fd_fraction':args.fd_fraction,
        'spinor_factor':args.spinor_factor,
        'vev':args.vev,
        'input_sha256':source_hash,
        'matched_scope':'coupled_EFT1_1111_tree_F_seven_term_plus_scalar_bubble',
    }
    run.mkdir(parents=True, exist_ok=True)
    settings_path=run/'settings.json'
    if settings_path.exists():
        if not args.resume:
            raise FileExistsError(f'{model}: prior optimisation exists; pass --resume or fresh --run-id')
        if read(settings_path) != manifest:
            raise ValueError(f'{model}: saved settings differ; choose another run ID')
    else:
        save(settings_path,manifest)

    fn, restore, _ = baseline.setup_physics(model, args.response_root)
    chol=np.linalg.cholesky(target.covariance)
    names=tuple(target.observable_names)
    compat=baseline.symbolic_export_compatibility(model,canonical,args.allow_stale_symbolic_source)
    evaluations=0
    failures=0
    best=None
    x_best=x_seed.copy()
    stop_message=''

    def archive_best(x,chi2,record,origin):
        nonlocal best,x_best
        x_best=np.asarray(x,float).copy()
        best={
            'format':'t3-coupled-yukawa-least-squares-v1',
            'model':model,'status':'Computed','chi2':float(chi2),
            'seed_chi2':float(old['chi2']),
            'yukawas_flat':x_best.tolist(),
            'prediction_1111':record['prediction'],
            'masses_ev':record.get('masses_ev_ascending'),
            'from_run':args.from_run,
            'improved_from_overnight':float(chi2) < float(old['chi2']),
            'eval_count':evaluations,'failed_evals':failures,
            'last_improvement':origin,
            'source_fingerprint_compatibility':args.allow_stale_symbolic_source,
        }
        save(checkpoint,best)
        p=baseline.params_from_vector(x_best)
        save(config_path,general.config_with_yukawas(config,p))
        if model not in MIXED:
            save(overrides_path,general.override_yukawas(overrides,p))
        print(f'{model}: BEST chi2={chi2:.8g} ({origin})',flush=True)

    try:
        with compat as compat_report, tempfile.TemporaryDirectory(prefix='t3_ls_'+tag+'_',dir=str(run)) as tmp:
            compat_path=run/'symbolic_export_compatibility.json'
            if compat_path.is_file() and read(compat_path)!=compat_report:
                raise ValueError(f'{model}: cached symbolic inputs changed since last refinement; use fresh --run-id')
            if not compat_path.exists():save(compat_path,compat_report)
            saved=Path(tmp)
            if model in MIXED:
                shutil.copy2(source/'result.json',saved/'result.json')

            if checkpoint.is_file():
                best=read(checkpoint)
                x_best=np.asarray(best['yukawas_flat'],float)
            else:
                # Initial checkpoint is only seeded after a fresh EFT1 evaluation.
                pass

            def evaluate(x,origin):
                nonlocal evaluations,failures
                if args.deadline is not None and time.monotonic() >= args.deadline:
                    raise BudgetReached('Global time budget reached')
                if evaluations >= args.max_evals:
                    raise BudgetReached('Per-model evaluation budget reached')
                x=np.asarray(x,float)
                if x.shape!=(18,) or not np.isfinite(x).all() or np.any(abs(x)>args.bound+1e-9):
                    raise ValueError('Yukawa candidate outside numerical bounds')
                evaluations+=1
                p=baseline.params_from_vector(x)
                save(saved/'fitted_config.json',general.config_with_yukawas(config,p))
                if model not in MIXED:
                    save(saved/'general_uv_overrides.json',general.override_yukawas(overrides,p))
                try:
                    if model in MIXED:
                        cases,_ =fn(model,saved,(release,channels),args.reference_running,registry,
                                    spinor_factor=args.spinor_factor,vev=args.vev)
                    else:
                        cases,_ =fn(model,saved,canonical,(release,channels),registry,
                                    spinor_factor=args.spinor_factor,
                                    wolframscript=args.wolframscript,vev=args.vev)
                    if set(cases)!={'1111'}:
                        raise ValueError(f'{model}: objective returned other cases: {list(cases)}')
                    record=cases['1111']
                    obs=record['prediction']
                    vector=np.asarray([obs[key] for key in names],float)
                    residual=np.linalg.solve(chol,vector-target.central_values)
                    chi2=float(residual@residual)
                    if not (np.isfinite(residual).all() and math.isfinite(chi2)):
                        raise ValueError('Nonfinite objective')
                    native=float(record['chi2'])
                    if not math.isclose(chi2,native,rel_tol=1e-6,abs_tol=1e-6):
                        raise ValueError(f'whitened objective differs from pipeline chi2: {chi2} vs {native}')
                    baseline.append_jsonl(history_path,{'eval':evaluations,'origin':origin,'status':'Computed','chi2':chi2})
                    if best is None or chi2 < best['chi2']-1e-10:
                        archive_best(x,chi2,record,origin)
                    return residual
                except (KeyboardInterrupt,SystemExit,BudgetReached):
                    raise
                except Exception as exc:
                    failures+=1
                    baseline.append_jsonl(history_path,{'eval':evaluations,'origin':origin,'status':'Failed',
                                                        'error':f'{type(exc).__name__}: {exc}'})
                    print(f'{model} evaluation {evaluations}: {type(exc).__name__}: {exc}',flush=True)
                    # Propagate: nonfinite or incorrect forward calculations are NOT valid fit residuals.
                    raise

            # Always reproduce the last best before further fitting; no implicit trust in stale results.
            try:
                initial=evaluate(x_best,'recheck_best')
            except BudgetReached as exc:
                stop_message=str(exc)
                return {'model':model,'chi2':float(best['chi2']) if best else None,'stopped':stop_message}

            for active in args.active_counts:
                if args.deadline is not None and time.monotonic() >= args.deadline:break
                if evaluations>=args.max_evals or float(best['chi2'])<=args.target_chi2:break
                # Approximate sensitivity of each Yukawa coordinate in whitened residual space.
                # This is a ranking heuristic only, not a physical gradient claim.
                anchor=x_best.copy()
                try:
                    base_r=evaluate(anchor,f'stage_{active}_anchor')
                except BudgetReached as exc:
                    stop_message=str(exc)
                    break
                scores=np.zeros(18)
                for i in range(18):
                    if evaluations>=args.max_evals: break
                    step=max(abs(anchor[i])*args.fd_fraction,args.fd_floor)
                    direction=1 if anchor[i]+step<=args.bound else -1
                    changed=anchor.copy()
                    changed[i]=np.clip(anchor[i]+direction*step,-args.bound,args.bound)
                    if changed[i]==anchor[i]:continue
                    try:
                        r=evaluate(changed,f'jac_rank_{active}_{i}')
                        scores[i]=np.linalg.norm((r-base_r)/(changed[i]-anchor[i]))
                    except BudgetReached:
                        break
                if evaluations>=args.max_evals:
                    stop_message='Per-model evaluation budget reached'
                    break
                ranking=np.argsort(-scores)
                selected=np.sort(ranking[:min(active,18)])
                if np.max(scores)<=0:
                    print(f'{model}: no reliable coordinate sensitivity; ending stage',flush=True)
                    break
                # Keep coordinates not selected at current best values.
                anchor=x_best.copy()
                lb=np.full(selected.size,-args.bound)
                ub=np.full(selected.size,args.bound)
                x0=np.clip(anchor[selected],lb+1e-12,ub-1e-12)
                print(f'{model}: refine active={active}, selected={selected.tolist()}, best={best["chi2"]:.8g}',flush=True)

                def residual_sub(z):
                    x=anchor.copy();x[selected]=z
                    return evaluate(x,f'ls_active_{active}')

                try:
                    least_squares(residual_sub,x0,bounds=(lb,ub),
                                  method='trf',max_nfev=args.stage_nfev,
                                  x_scale='jac',diff_step=args.ls_diff_step,
                                  ftol=1e-8,xtol=1e-8,gtol=1e-8,verbose=0)
                except BudgetReached as exc:
                    stop_message=str(exc)
                    print(f'{model}: {stop_message}',flush=True)
                    break
                # Failure of one branch should not corrupt saved best.
                save(run/f'stage_{active}_summary.json',{
                    'stage':active,'eval_count':evaluations,'best_chi2':float(best['chi2']),
                    'selected':selected.tolist(),
                })
            final={'model':model,'best_chi2':float(best['chi2']),
                   'overnight_chi2':float(old['chi2']),
                   'eval_count_this_invocation':evaluations,
                   'failures':failures,'stopped':stop_message or 'stages_completed'}
            save(run/'summary.json',final)
            return final
    finally:
        restore()


def main(argv=None):
    from scripts import RunT3GeneralAll16FromReference as general
    from scripts import RunT3GeneralSequentialFitAblation as original
    from Numerical.fitting.OscillationFit import OscillationFitTarget
    from Numerical.orchestration.RunSmallLambdaRefit import TARGET
    from Numerical.joint_scalar_matching.fermion_boundary_bridge import load_release_manifest,load_source

    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-id',required=True)
    p.add_argument('--from-run',default='yukawa_overnight_compat_v2')
    p.add_argument('--source-run',default='thesis16_merged_v1')
    p.add_argument('--models',default='D_+0,C_+1,A_+0,D_-2')
    p.add_argument('--active-counts',default='5,8,12')
    p.add_argument('--max-evals',type=int,default=120)
    p.add_argument('--stage-nfev',type=int,default=20)
    p.add_argument('--target-chi2',type=float,default=5.86)
    p.add_argument('--fd-fraction',type=float,default=.01)
    p.add_argument('--fd-floor',type=float,default=1e-5)
    p.add_argument('--ls-diff-step',type=float,default=1e-4)
    p.add_argument('--bound',type=float,default=3.)
    p.add_argument('--max-hours',type=float,default=0.)
    p.add_argument('--response-root',type=Path,default=ROOT/'output/audits/batch3b_quartic_fixed_20261111_234059')
    p.add_argument('--release-dir',type=Path,default=ROOT/'output/general_observable/lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f')
    p.add_argument('--v4-channel-dir',type=Path,default=ROOT/'output/general_observable/lo_joint_scalar_matching/source_quartic_v4/native_quartic_channels')
    p.add_argument('--quartic-seed-dir',type=Path,default=ROOT/'output/full16_symbolic')
    p.add_argument('--registry',type=Path,default=ROOT/'RGE/general/generated/general_rgbeta_quartic_registry.json')
    p.add_argument('--reference-running',type=Path,default=ROOT/'output/general_observable/all16_v1')
    p.add_argument('--target',type=Path,default=TARGET)
    p.add_argument('--vev',type=float,default=246.22)
    p.add_argument('--spinor-factor',type=float,default=1.)
    p.add_argument('--wolframscript',default='wolframscript')
    p.add_argument('--allow-stale-symbolic-source',action='store_true')
    p.add_argument('--resume',action='store_true')
    p.add_argument('--check-only',action='store_true')
    args=p.parse_args(argv)
    if not args.run_id.replace('_','').replace('-','').isalnum():p.error('Invalid run-id')
    models=tuple(x.strip() for x in args.models.split(',') if x.strip())
    if not models or len(set(models))!=len(models) or set(models)-set(general.MODELS):p.error('Invalid models')
    try:args.active_counts=tuple(int(x) for x in args.active_counts.split(','))
    except Exception:p.error('Invalid --active-counts')
    if not args.active_counts or any(not 1<=x<=18 for x in args.active_counts):p.error('active-counts must lie in 1..18')
    if (args.max_evals<1 or args.stage_nfev<1 or args.bound<=0 or args.target_chi2<=0
            or args.fd_fraction<=0 or args.fd_floor<=0 or args.ls_diff_step<=0 or args.max_hours<0):p.error('Invalid optimisation options')
    args.deadline=time.monotonic()+args.max_hours*3600 if args.max_hours else None
    fit_root=baseline.SOURCE/args.source_run
    if not (fit_root/'settings.json').is_file():raise FileNotFoundError('Missing fitted source run')
    canonical=read(fit_root/'settings.json').get('canonical_run_id','full16_general_v1')
    root=DEFAULT_OUT/args.run_id
    if root.exists() and not (args.resume or args.check_only):raise FileExistsError(f'Output exists: {root}')
    target=OscillationFitTarget.from_json(args.target)
    if args.check_only:
        ok=True
        for model in models:
            tag=general.tag(model)
            folder=fit_root/tag
            prior=BASE/'coupled_yukawa_overnight'/args.from_run/tag/'best.json'
            required=[folder/'fitted_config.json',prior,args.registry,args.release_dir/tag/'boundary.json',args.v4_channel_dir/(tag+'.json')]
            if model not in MIXED:required.append(folder/'general_uv_overrides.json')
            if model in SPECIAL:required.append(args.response_root)
            missing=[str(x) for x in required if not x.exists()]
            print(model,'READY' if not missing else 'MISSING',missing,flush=True)
            ok=ok and not missing
        return 0 if ok else 2
    root.mkdir(parents=True,exist_ok=True)
    registry=read(args.registry)
    releases=load_release_manifest(args.release_dir)
    rows=[]
    for model in models:
        tag=general.tag(model)
        try:
            ref=load_source(releases[model],expected_model=model)
            ch=original._load_channels(model,ref,v4=args.v4_channel_dir,seeds=args.quartic_seed_dir)
            row=run_one(model,args,fit_root,root,registry,ref,ch,canonical,target)
        except (KeyboardInterrupt,SystemExit):raise
        except Exception as exc:
            row={'model':model,'status':'Failed','error':f'{type(exc).__name__}: {exc}',
                 'traceback':traceback.format_exc(limit=7)}
            save(root/tag/'error.json',row)
            print(f'{model}: ERROR {row["error"]}',flush=True)
        rows.append(row)
        save(root/'summary.json',{'format':'t3-coupled-yukawa-least-squares-v1','models':rows})
    print('Saved:',root/'summary.json',flush=True)
    return 2 if any(x.get('status')=='Failed' for x in rows) else 0


if __name__=='__main__':
    raise SystemExit(main())
