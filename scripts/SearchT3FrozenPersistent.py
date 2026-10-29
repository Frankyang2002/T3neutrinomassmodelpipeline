"""Persistent Sobol + local-18 frozen-UV Yukawa search.

Runs indefinitely until target chi2, Ctrl-C, or optional --max-cycles.
No UV running and no modification of physical screening/matching.
"""
from __future__ import annotations
import argparse
import json
import math
import re
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from scipy.stats import qmc
from scripts import RunT3GeneralScalar48Campaign as bulk
from scripts.RefineT3Frozen18Stable import eligible_candidate, load_seeds, write_json

FORMAT = 't3-frozen-persistent-sobol-v1'
CASES = ('small_uv_frozen_5e-9', 'large_uv_frozen_0p5')


def sobol_batch(lo, hi, points, seed):
    if points < 2 or points & (points-1):
        raise ValueError('Sobol points must be a power of two >=2')
    draws = qmc.Sobol(d=len(lo), scramble=True, seed=seed).random_base2(int(math.log2(points)))
    return lo + draws * (hi-lo)


def choose_seeds(rows, count, lo=None, hi=None, min_distance=0.08):
    """Select the best point plus fitness-ranked, geometrically distinct starts.

    Distance is RMS displacement in normalized [0,1] Yukawa coordinates.
    If no sufficiently separated candidates exist, retain the best available
    distinct points rather than fabricating starts.
    """
    ranked = sorted((r for r in rows if eligible_candidate(r)),
                    key=lambda r: float(r['chi2']))
    if not ranked or count <= 0:
        return []
    names = bulk.NAMES
    lo = np.zeros(len(names)) if lo is None else np.asarray(lo, dtype=float)
    hi = np.ones(len(names)) if hi is None else np.asarray(hi, dtype=float)
    width = hi - lo
    if np.any(width <= 0):
        raise ValueError('Invalid search bounds')

    def coords(row):
        return (np.array([float(row['parameters'][n]) for n in names])-lo)/width

    selected = [ranked[0]]
    vectors = [coords(ranked[0])]
    used = {tuple(vectors[0])}
    # Prefer competitive candidates, while requiring genuinely different starts.
    pool = ranked[1:max(32, 20*count)]
    for row in pool:
        vec = coords(row)
        key = tuple(vec)
        if key in used:
            continue
        distance = min(np.linalg.norm(vec-other)/np.sqrt(len(vec)) for other in vectors)
        if distance >= min_distance:
            selected.append(row)
            vectors.append(vec)
            used.add(key)
            if len(selected) >= count:
                return selected
    # Fallback: choose the most distant distinct candidate from the available pool.
    while len(selected) < count:
        choices=[]
        for row in ranked[1:]:
            vec=coords(row)
            if tuple(vec) in used:
                continue
            distance=min(np.linalg.norm(vec-other)/np.sqrt(len(vec)) for other in vectors)
            choices.append((distance, -float(row['chi2']), row, vec))
        if not choices:
            break
        _, _, row, vec = max(choices, key=lambda item:(item[0],item[1]))
        selected.append(row)
        vectors.append(vec)
        used.add(tuple(vec))
    return selected


def mass_scale_factor(prediction, target_dm31, *, max_factor=1e5):
    """Common Yukawa multiplier assuming fixed-scalar bilinear C5.

    Return None for undefined, wrong-sign, or unrepresentable atmospheric gaps.
    Does NOT certify scaling: evaluate the transformed candidate authoritatively.
    """
    try:
        dm = float(prediction['delta_m3l_sq_ev2'])
        target = float(target_dm31)
        if not (np.isfinite(dm) and np.isfinite(target) and dm * target > 0):
            return None
        factor = (abs(target) / abs(dm)) ** .25
        if not np.isfinite(factor) or not (1 / max_factor <= factor <= max_factor):
            return None
        return factor
    except (KeyError, ValueError, TypeError, OverflowError, ZeroDivisionError):
        return None


def normalize_mass_candidate(row, lo, hi, target_dm31, evaluate, *, rtol=.03):
    """Mass-normalize a valid point and verify quartic scaling numerically.

    Physical screening is always performed by evaluate.  Never replace a valid
    original point with an invalid transformed point.  The scale ratio check
    excludes candidates for which the assumed homogeneity is not supported.
    """
    if not eligible_candidate(row) or not isinstance(row.get('prediction'), dict):
        return row, 'invalid_input'
    factor = mass_scale_factor(row['prediction'], target_dm31)
    if factor is None:
        return row, 'unscalable'
    x = np.array([float(row['parameters'][n]) for n in bulk.NAMES])
    scaled = factor * x
    if not np.all(np.isfinite(scaled)) or np.any(scaled < lo) or np.any(scaled > hi):
        return row, 'outside_search_box'
    if abs(factor - 1) < 1e-10:
        return row, 'already_normalized'
    transformed = evaluate(scaled)
    if not eligible_candidate(transformed):
        return row, 'screen_rejected'
    try:
        observed = float(transformed['prediction']['delta_m3l_sq_ev2'])
        predicted = float(row['prediction']['delta_m3l_sq_ev2']) * factor**4
        if not np.isfinite(observed) or abs(observed - predicted) > rtol * max(abs(predicted), 1e-30):
            return row, 'scaling_mismatch'
    except (KeyError, TypeError, ValueError, OverflowError):
        return row, 'missing_normalized_prediction'
    return transformed, 'normalized'


def parse_radii(value):
    radii = tuple(float(x.strip()) for x in value.split(','))
    if not radii or not all(np.isfinite(x) and 0 < x <= .2 for x in radii):
        raise ValueError('--radii requires comma-separated radii in (0,0.2]')
    return radii


def run(args):
    if args.case not in CASES:
        raise ValueError('Frozen scenarios only')
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,63}', args.run_id):
        raise ValueError('Invalid run id')
    if args.points < 2 or args.points & (args.points-1):
        raise ValueError('--points must be a power of 2')
    radii = parse_radii(args.radii)
    if not (0 < args.mass_scaling_rtol <= 0.25):
        raise ValueError('Invalid --mass-scaling-rtol')
    if not (0 <= args.min_seed_distance <= 1):
        raise ValueError('Invalid --min-seed-distance')
    if not (1 <= args.seeds <= 20 and args.nfev >= 20 and 0 < args.radius <= 0.2 and 0 < args.diff_step <= .01):
        raise ValueError('Invalid optimisation options')
    reg = bulk.read(bulk.REGISTRY)
    entry, new, reused = bulk.inventory(reg, args.model)
    if reused:
        raise ValueError('Mixed-scalar B(-1)/C(-1) require their special evaluator')
    hist = bulk.historical(args.case, args.model, entry)
    if hist['uv_on']:
        raise ValueError('UV-running case not allowed')
    beta_path = bulk.OUT / args.source_campaign / 'betas' / (bulk.safe(args.model) + '.json')
    payload = bulk.read(beta_path)
    kernel, matching_path, matching_kind = bulk.matching(args.model)
    from Numerical.orchestration import RunSmallLambdaRefit as refit
    from Numerical.fitting.OscillationFit import OscillationFitTarget
    target = OscillationFitTarget.from_json(refit.TARGET)
    evaluator = bulk.Evaluator(args.case, args.model, entry, new, hist, payload, kernel, target, args.extra_real)
    lo, hi = bulk.b21.search_box(hist['seed'], local_fraction=1.5,
                                 floor_fraction=.05, abs_yukawa_max=math.sqrt(4*math.pi))
    width = hi-lo
    target_dm31 = float(target.central_values[target.observable_names.index('delta_m3l_sq_ev2')])
    dest = bulk.OUT / args.run_id / args.case / bulk.safe(args.model)
    state_file = dest / 'state.json'
    history_file = dest / 'cycles.jsonl'
    dest.mkdir(parents=True, exist_ok=True)
    provenance = {'source_campaign':args.source_campaign, 'beta_sha256':bulk.sha(beta_path),
                  'matching_sha256':bulk.sha(matching_path), 'historical_sha256':bulk.sha(hist['cfg_path']),
                  'points':args.points,'seed':args.seed,'radius':args.radius,'diff_step':args.diff_step,
                  'local_seeds':args.seeds,'nfev':args.nfev,'target':args.target,
                  'radii':list(radii),'min_seed_distance':args.min_seed_distance,
                  'extra_real':args.extra_real,'mass_normalization':True,
                  'mass_normalization_rtol':args.mass_scaling_rtol}
    if state_file.exists():
        if not args.resume:
            raise FileExistsError('Existing run: use --resume or a new --run-id')
        old = bulk.read(state_file)
        if old.get('inputs') != provenance or old.get('model') != args.model or old.get('case') != args.case:
            raise ValueError('Checkpoint inputs differ; use a new --run-id')
        cycle = int(old['completed_cycles'])
        best = old.get('best')
    else:
        if args.resume:
            raise FileNotFoundError('No checkpoint exists to resume')
        cycle = 0
        best = None
        for item, _ in load_seeds(args.model,args.case,args.source_campaign,args.from_run):
            if best is None or float(item['chi2']) < float(best['chi2']):
                best = item
        # Check imported best on the current evaluator before using its chi2 as trusted.
        if best:
            x = np.array([best['parameters'][n] for n in bulk.NAMES])
            checked = evaluator.evaluate(x,'imported_best')
            best = checked if eligible_candidate(checked) else None
        write_json(state_file, {'format':FORMAT,'model':args.model,'case':args.case,
                                'inputs':provenance,'completed_cycles':cycle,'best':best,'status':'Running'})
    print(f'RESUME {args.model} {args.case} cycles={cycle} best={None if best is None else best["chi2"]}',flush=True)

    def consider(row):
        nonlocal best
        if eligible_candidate(row) and (best is None or float(row['chi2']) < float(best['chi2'])):
            best = row
            print(f'NEW BEST chi2={best["chi2"]:.9g}',flush=True)

    while best is None or float(best['chi2']) > args.target:
        if args.max_cycles and cycle >= args.max_cycles:
            print('BUDGET REACHED; resume with larger --max-cycles or omit it',flush=True)
            break
        number = cycle + 1
        # Deterministic unique scrambled Sobol batches across resumed cycles.
        draw_seed = int((args.seed + 104729 * number) % (2**32))
        candidates = []
        normalization_counts = {}
        if best is not None:
            candidates.append(best)
        print(f'CYCLE {number} Sobol={args.points} seed={draw_seed}',flush=True)
        for i,x in enumerate(sobol_batch(lo,hi,args.points,draw_seed)):
            row = evaluator.evaluate(x,f'cycle_{number}_sobol_{i+1}')
            if eligible_candidate(row):
                # Keep the original only if no verified mass-normalized counterpart
                # is available. Both are evaluated with the same physics code.
                consider(row)
                normalized, norm_status = normalize_mass_candidate(
                    row, lo, hi, target_dm31,
                    lambda xx, i=i: evaluator.evaluate(xx, f'cycle_{number}_massnorm_{i+1}'),
                    rtol=args.mass_scaling_rtol)
                normalization_counts[norm_status] = normalization_counts.get(norm_status, 0) + 1
                candidates.append(normalized)
                consider(normalized)
                if best and float(best['chi2']) <= args.target:
                    break
        seeds = choose_seeds(candidates,args.seeds,lo,hi,args.min_seed_distance)
        stages=[]
        if best is None or float(best['chi2'])>args.target:
            for idx, seed_row in enumerate(seeds,1):
                # Fresh local window around selected point. xtol ends this attempt, not campaign.
                seed_x=np.array([seed_row['parameters'][n] for n in bulk.NAMES])
                z=np.clip((seed_x-lo)/width,0,1)
                radius=radii[(number-1+idx-1)%len(radii)]
                low=np.maximum(0,z-radius)
                high=np.minimum(1,z+radius)
                low=np.maximum(0,np.minimum(low,z-1e-12))
                high=np.minimum(1,np.maximum(high,z+1e-12))
                if np.any(high<=low):
                    continue
                errors=[0]
                start_chi2=float(seed_row["chi2"])
                last_valid=[None]
                def residual(q):
                    row=evaluator.evaluate(lo+width*q, f'cycle_{number}_local_{idx}')
                    consider(row)
                    if eligible_candidate(row):
                        last_valid[0] = row
                    v=np.asarray(row.get('whitened_residual',[]),dtype=float)
                    if v.shape!=(5,) or not np.all(np.isfinite(v)):
                        errors[0]+=1
                        return np.full(5,1e5)
                    return v
                sol=least_squares(residual,z,bounds=(low,high),method='trf',jac='3-point',
                                  diff_step=args.diff_step,x_scale='jac',max_nfev=args.nfev,
                                  ftol=1e-12,xtol=1e-12,gtol=1e-12)
                final_x=lo+width*sol.x
                final_row=evaluator.evaluate(final_x,f'cycle_{number}_local_{idx}_final')
                consider(final_row)
                final_good=eligible_candidate(final_row)
                final_residual=final_row.get('whitened_residual') if final_good else None
                final_chi2=float(final_row['chi2']) if final_good else None
                displacement=float(np.linalg.norm(sol.x-z)/np.sqrt(len(z)))
                stages.append({'seed':idx,'radius':radius,'start_chi2':start_chi2,
                               'start_parameters':seed_row['parameters'],
                               'start_whitened_residual':seed_row.get('whitened_residual'),
                               'final_chi2':final_chi2,
                               'final_parameters':final_row.get('parameters') if final_good else None,
                               'final_whitened_residual':final_residual,
                               'normalized_rms_displacement':displacement,
                               'nfev':int(sol.nfev),'message':str(sol.message),
                               'status':int(sol.status),'rejected_evaluations':errors[0]})
                print(f'  LOCAL {idx}/{len(seeds)} radius={radius} start={start_chi2:.6g} '
                      f'final={final_chi2} displacement={displacement:.4g} '
                      f'nfev={sol.nfev} {sol.message}',flush=True)
                if best and float(best['chi2'])<=args.target:
                    break
        cycle=number
        record={'cycle':cycle,'best_chi2':None if best is None else float(best['chi2']),
                'sobol_points':args.points,'normalization_counts':normalization_counts,
                'local_stages':stages,'total_evaluations':evaluator.count}
        with history_file.open('a',encoding='utf-8') as f:
            f.write(json.dumps(record,allow_nan=False)+'\n')
        status='Accepted' if best and float(best['chi2'])<=args.target else 'Running'
        write_json(state_file,{'format':FORMAT,'model':args.model,'case':args.case,
                               'inputs':provenance,'completed_cycles':cycle,'best':best,
                               'status':status,'total_evaluations_this_invocation':evaluator.count})
        if best:
            cfg=bulk.update_yukawas(hist,np.array([best['parameters'][n] for n in bulk.NAMES]))
            cfg['conditional_general_scalar_fit']={'model':args.model,'case':args.case,
                'chi2':best['chi2'],'source':'persistent_sobol_full18',
                'sequential_matching_certified':False}
            write_json(dest/'fitted_config.json',cfg)
        print(f'CHECKPOINT cycle={cycle} best_chi2={record["best_chi2"]} massnorm={normalization_counts}',flush=True)
    if best and float(best['chi2'])<=args.target:
        print(f'ACCEPTED {args.model} {args.case} chi2={best["chi2"]}',flush=True)
    return 0


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-id',required=True)
    p.add_argument('--model',default='A_+0')
    p.add_argument('--case',default='small_uv_frozen_5e-9',choices=CASES)
    p.add_argument('--source-campaign',default='overnight_six_v1_bulk48')
    p.add_argument('--from-run',action='append',default=[])
    p.add_argument('--points',type=int,default=256)
    p.add_argument('--seeds',type=int,default=3)
    p.add_argument('--nfev',type=int,default=300)
    p.add_argument('--radius',type=float,default=.02,help='Legacy option; radii schedule overrides it')
    p.add_argument('--radii',default='0.02,0.05,0.10')
    p.add_argument('--min-seed-distance',type=float,default=.08)
    p.add_argument('--diff-step',type=float,default=1e-6)
    p.add_argument('--mass-scaling-rtol',type=float,default=.03)
    p.add_argument('--seed',type=int,default=20261030)
    p.add_argument('--target',type=float,default=5.86)
    p.add_argument('--extra-real',type=float,default=-.01)
    p.add_argument('--max-cycles',type=int,default=0,help='0 means unlimited cycles')
    p.add_argument('--resume',action='store_true')
    return run(p.parse_args(argv))

if __name__=='__main__':
    raise SystemExit(main())
