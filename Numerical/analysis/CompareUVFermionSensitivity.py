"""Compare saved UV-refit (111) and fermion-start (011) T1/T2/T3 slopes.

Read-only analysis: no fits, RGE solves or matching calculations. UV scans are
loaded from output/studies/lambda_5e-9_refit/comparison/{mode}/sensitivity.
Both scan families are evaluated using the SAME T1/T2/T3 definitions, but may
have different baseline fitted inputs and different gamma grids.

A local log slope is estimated by quadratic interpolation through gamma=1 and
the nearest successful point on either side:
    d ln|O| / d ln(gamma) at gamma=1.
The fermion-start gamma spacing is much wider than the UV scan spacing. This
is a finite-difference estimate, not a controlled 111-vs-011 ablation.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OBS = ('delta_m21_sq_ev2','delta_m3l_sq_ev2','sin2_theta12','sin2_theta13','sin2_theta23')
DIRECTIONS = ('T1','T2','T3')
MODES = ('optimised','shared_a_m2')
UV_ROOT = ROOT / 'output/studies/lambda_5e-9_refit/comparison'
FER_ROOT = ROOT / 'output/studies/fermion_start_lambda_comparison/5e-9/reports_existing'
OUT = ROOT / 'output/studies/fermion_start_lambda_comparison/uv_sensitivity_comparison'


def read(path: Path) -> dict:
    return json.loads(path.read_text(encoding='utf-8-sig'))


def samples(scan: dict, direction: str, observable: str) -> list[tuple[float,float]]:
    result = []
    for entry in scan.get('directions',{}).get(direction,{}).get('points',[]):
        if entry.get('status','Success') != 'Success':
            continue
        try:
            g = float(entry['gamma'])
            o = float(entry['prediction'][observable])
        except (KeyError, TypeError, ValueError):
            continue
        if math.isfinite(g) and g > 0 and math.isfinite(o):
            result.append((g,o))
    return sorted(result)


def local_slope(entries: list[tuple[float,float]]) -> dict:
    """Three-point log-space derivative with unequal left/right spacings."""
    centres = [(g,v) for g,v in entries if abs(g-1.0)<1e-10]
    left = [(g,v) for g,v in entries if g < 1.0-1e-10]
    right = [(g,v) for g,v in entries if g > 1.0+1e-10]
    if len(centres)!=1 or not left or not right:
        raise ValueError('Need exactly one gamma=1 point and successful points on both sides')
    g0,v0=centres[0]
    gm,vm=max(left,key=lambda z:z[0])
    gp,vp=min(right,key=lambda z:z[0])
    if min(abs(vm),abs(v0),abs(vp))<=0:
        raise ValueError('Logarithmic observable derivative undefined at zero')
    xm,xp=math.log(gm),math.log(gp)
    ym=math.log(abs(vm)/abs(v0))
    yp=math.log(abs(vp)/abs(v0))
    slope=(ym*xp*xp-yp*xm*xm)/(xm*xp*(xp-xm))
    return {'slope':slope,'gamma_minus':gm,'gamma_plus':gp,'baseline':v0,
            'n_valid_points':len(entries), 'max_fractional_triplet':max(abs(vm/v0-1),abs(vp/v0-1))}


def definition(scan: dict, direction: str) -> str:
    return str(scan.get('directions',{}).get(direction,{}).get('definition','')).replace(' ','').lower()


def save_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main(argv=None) -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--uv-root',type=Path,default=UV_ROOT)
    parser.add_argument('--fermion-root',type=Path,default=FER_ROOT)
    parser.add_argument('--output',type=Path,default=OUT)
    args=parser.parse_args(argv)
    rows=[];warnings=[];counts={}
    for mode in MODES:
        uv_dir=args.uv_root/mode/'sensitivity'
        ferm_dir=args.fermion_root/mode/'sensitivity'
        ufiles={p.parent.name:p for p in uv_dir.glob('*/sensitivity_diagnostics.json')}
        ffiles={p.parent.name:p for p in ferm_dir.glob('*/sensitivity_diagnostics.json')}
        paired=sorted(set(ufiles)&set(ffiles))
        counts[mode]={'uv_scans':len(ufiles),'fermion_scans':len(ffiles),'paired_models':len(paired)}
        for model in sorted(set(ufiles)^set(ffiles)):
            warnings.append(f'{mode}/{model}: missing in '+('UV' if model not in ufiles else 'fermion-start'))
        for model in paired:
            try:
                u=read(ufiles[model]);f=read(ffiles[model])
            except (ValueError,OSError) as ex:
                warnings.append(f'{mode}/{model}: cannot read JSON: {ex}');continue
            if u.get('status')!='Success' or f.get('status')!='Success':
                warnings.append(f'{mode}/{model}: scan status not Success');continue
            for d in DIRECTIONS:
                ud=definition(u,d);fd=definition(f,d)
                if not ud or not fd or ud != fd:
                    warnings.append(f'{mode}/{model}/{d}: nonmatching or absent direction definitions');continue
                for obs in OBS:
                    try:
                        um=local_slope(samples(u,d,obs))
                        fm=local_slope(samples(f,d,obs))
                    except ValueError as ex:
                        warnings.append(f'{mode}/{model}/{d}/{obs}: {ex}');continue
                    rows.append({'mode':mode,'model':model,'direction':d,'observable':obs,
                        'uv_log_slope':um['slope'],'fermion_log_slope':fm['slope'],
                        'uv_abs_slope':abs(um['slope']),'fermion_abs_slope':abs(fm['slope']),
                        'uv_minus_fermion_slope':um['slope']-fm['slope'],
                        'uv_baseline':um['baseline'],'fermion_baseline':fm['baseline'],
                        'uv_gamma_minus':um['gamma_minus'],'uv_gamma_plus':um['gamma_plus'],
                        'fermion_gamma_minus':fm['gamma_minus'],'fermion_gamma_plus':fm['gamma_plus'],
                        'uv_n_valid_points':um['n_valid_points'],'fermion_n_valid_points':fm['n_valid_points'],
                        'uv_triplet_max_fraction':um['max_fractional_triplet'],
                        'fermion_triplet_max_fraction':fm['max_fractional_triplet']})
    fields=['mode','model','direction','observable','uv_log_slope','fermion_log_slope',
            'uv_abs_slope','fermion_abs_slope','uv_minus_fermion_slope',
            'uv_baseline','fermion_baseline','uv_gamma_minus','uv_gamma_plus',
            'fermion_gamma_minus','fermion_gamma_plus','uv_n_valid_points',
            'fermion_n_valid_points','uv_triplet_max_fraction','fermion_triplet_max_fraction']
    output=args.output
    output.mkdir(parents=True,exist_ok=True)
    save_csv(output/'local_slopes_by_model.csv',rows,fields)
    summary=[]
    for mode in MODES:
        for obs in OBS:
            sub=[r for r in rows if r['mode']==mode and r['observable']==obs]
            if not sub:continue
            umax=max(sub,key=lambda r:r['uv_abs_slope'])
            fmax=max(sub,key=lambda r:r['fermion_abs_slope'])
            summary.append({'mode':mode,'observable':obs,'paired_model_directions':len(sub),
                'max_uv_abs_log_slope':umax['uv_abs_slope'],'uv_peak_model':umax['model'],
                'uv_peak_direction':umax['direction'],
                'max_fermion_abs_log_slope':fmax['fermion_abs_slope'],
                'fermion_peak_model':fmax['model'],'fermion_peak_direction':fmax['direction']})
    sfields=['mode','observable','paired_model_directions','max_uv_abs_log_slope',
             'uv_peak_model','uv_peak_direction','max_fermion_abs_log_slope',
             'fermion_peak_model','fermion_peak_direction']
    save_csv(output/'summary.csv',summary,sfields)
    meta={'uv_root':str(args.uv_root),'fermion_root':str(args.fermion_root),
          'counts':counts,'metric_rows':len(rows),'warnings':warnings,
          'method':'Three-point quadratic derivative of ln|observable| vs ln(gamma) at gamma=1 using nearest successful neighbors.',
          'limitations':['Different fitted baselines: not controlled 111 vs 011 ablation.',
          'Different gamma step sizes: local derivative estimates have different truncation errors.',
          'Near-zero fermion slopes (~1e-14) are numerical precision noise.']}
    (output/'summary.json').write_text(json.dumps(meta,indent=2)+'\n',encoding='utf-8')
    print('Counts:',counts,'| metric rows:',len(rows),'| warnings:',len(warnings))
    for s in summary:
        print(f"{s['mode']:12} {s['observable']:23} UV={s['max_uv_abs_log_slope']:.6g}  fermion={s['max_fermion_abs_log_slope']:.6g}")
    print('Saved:',output)
    return 0 if rows else 2


if __name__=='__main__':
    raise SystemExit(main())
