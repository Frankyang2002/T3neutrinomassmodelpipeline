"""Study B: T3 heavy-fermion mass/threshold-scale scan with frozen controls.

Run from repository root:
 python -m Numerical.orchestration.RunHeavyFermionScaleComparison
 python -m Numerical.orchestration.RunHeavyFermionScaleComparison --run --limit 2
 python -m Numerical.orchestration.RunHeavyFermionScaleComparison --run

The same fitted model and common-lambda point is used at each mass factor.
Both MF and mu_fermion_threshold_gev are scaled; scalar threshold, UV scale,
Yukawas and lambdaT3 remain unchanged. All observables are evaluated from the
existing RGE pipeline and the authoritative C5 matching expression.
The frozen control changes MF and its matching scale identically but turns off
all RG evolution. It is a sensitivity diagnostic, not a full fixed-order EFT.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
from copy import deepcopy
from pathlib import Path

import numpy as np

from Numerical.orchestration.ResumePipelineNumericalResults import PROJECT_ROOT, OUTPUT_ROOT, REQUIRED, discover_cases
from Numerical.orchestration.RunFittedSensitivityComparison import best_candidate, OBS
from Numerical.orchestration.RunCommonScalarComparison import deform_common_scalar, evaluate
from Numerical.orchestration.RunCommonScalarRGControls import fixed_point
from Numerical.fitting.OscillationFit import OscillationFitTarget

TARGET = PROJECT_ROOT / 'configs' / 'nufit_6_1_ic24_no.json'


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def scale_heavy_mass(payload, factor):
    """Multiply every MF matrix element and mu_F by the same positive factor."""
    result = deepcopy(payload)
    spec = result['base_state']['ordinary']['MF']
    spec['real'] = (np.asarray(spec['real'], dtype=float) * factor).tolist()
    if 'imag' in spec:
        spec['imag'] = (np.asarray(spec['imag'], dtype=float) * factor).tolist()
    result['scales']['mu_fermion_threshold_gev'] = float(result['scales']['mu_fermion_threshold_gev']) * factor
    return result


def check_scales(payload):
    s = payload['scales']
    u = float(payload['base_state']['mu_gev'])
    f = float(s['mu_fermion_threshold_gev'])
    ss = float(s['mu_scalar_threshold_gev'])
    low = float(s['mu_low_gev'])
    if not u > f > ss > low > 0:
        raise ValueError(f'Threshold hierarchy violated: mu_UV={u:g}, mu_F={f:g}, mu_S={ss:g}, mu_low={low:g}')
    return u, f, ss, low


def save_plots(rows, folder):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print('matplotlib unavailable: CSV/JSON are complete')
        return
    good = [r for r in rows if r['status'] == 'Success']
    if not good:
        print('No successful scan points; no figures created')
        return
    folder.mkdir(parents=True, exist_ok=True)
    palette = plt.get_cmap('tab20')
    models = sorted(set(r['model'] for r in good))
    ys = sorted(set(r['reference_yukawa_target'] for r in good))
    labels = {'delta_m21_sq_ev2':r'$\Delta m^2_{21}$ (eV$^2$)',
              'delta_m3l_sq_ev2':r'$\Delta m^2_{3\ell}$ (eV$^2$)',
              'sin2_theta12':r'$\sin^2\theta_{12}$',
              'sin2_theta13':r'$\sin^2\theta_{13}$',
              'sin2_theta23':r'$\sin^2\theta_{23}$'}
    for y in ys:
        for ob in OBS:
            fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharex=True)
            for j, model in enumerate(models):
                pts = sorted((r for r in good if r['model']==model and r['reference_yukawa_target']==y),
                             key=lambda r:r['mass_factor'])
                if not pts: continue
                c = palette(j % 20)
                x = [p['mu_fermion_threshold_gev'] for p in pts]
                axes[0].plot(x, [p['running_'+ob] for p in pts], '-o', lw=1.2, ms=3, color=c, label=model)
                axes[0].plot(x, [p['fixed_'+ob] for p in pts], '--', lw=1.1, color=c)
                axes[1].plot(x, [p['relative_percent_'+ob] for p in pts], '-o', lw=1.2, ms=3, color=c)
            axes[0].set(ylabel=labels[ob], title='RGE (solid) and frozen parameters (dashed)')
            axes[1].set(ylabel=r'$100\,(X_{\rm RG}-X_{\rm fixed})/X_{\rm fixed}$ (%)',title='Isolated RGE correction')
            for ax in axes:
                ax.set_xscale('log')
                ax.set_xlabel(r'$\mu_F$ (GeV); $M_F$ scaled by the same factor')
                ax.grid(alpha=.25)
            axes[0].legend(fontsize=5, ncol=2)
            fig.suptitle(f'Study B: reference Yukawa target = {y:g}, {ob}')
            fig.tight_layout()
            fig.savefig(folder/f'yref_{y:g}_{ob}.png', dpi=180)
            plt.close(fig)
    print('Saved figures:', folder)


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--scenario',default='smallY_smallL')
    p.add_argument('--run',action='store_true')
    p.add_argument('--limit',type=int,default=0)
    p.add_argument('--mass-factors',nargs='+',type=float,default=[0.5,1.,2.])
    p.add_argument('--reference-yukawas',nargs='+',type=float,default=[0.01,1.])
    p.add_argument('--force',action='store_true')
    args=p.parse_args(argv)
    if args.limit<0 or not args.mass_factors or not args.reference_yukawas or any(not math.isfinite(x) or x<=0 for x in args.mass_factors+args.reference_yukawas):
        p.error('Require nonnegative --limit, and positive finite mass factors/reference Yukawas')
    if len(set(args.mass_factors))!=len(args.mass_factors) or len(set(args.reference_yukawas))!=len(args.reference_yukawas):
        p.error('Scan settings must not contain duplicates')
    scan_root=PROJECT_ROOT/'Reports'/'output'/'full'/'common_scalar_comparison'/args.scenario
    out=PROJECT_ROOT/'Reports'/'output'/'full'/'heavy_fermion_scale_comparison'/args.scenario
    cases={key:case for sc,key,case in discover_cases(OUTPUT_ROOT) if sc==args.scenario}
    keys=sorted(x.name for x in scan_root.iterdir() if x.is_dir() and (x/'common_scalar_scan.json').is_file()) if scan_root.is_dir() else []
    keys=keys[:args.limit or None]
    selected={}
    for key in keys:
        best=best_candidate(args.scenario,key)
        if best is None or key not in cases: continue
        paths=[cases[key]/'data'/name for name in REQUIRED]
        if any(not x.is_file() for x in paths): continue
        scan=read_json(scan_root/key/'common_scalar_scan.json')
        points={float(pt['reference_yukawa_target']):pt for pt in scan['points'] if pt['status']=='Success'}
        if not all(any(np.isclose(y,v,rtol=1e-10,atol=0) for v in points) for y in args.reference_yukawas): continue
        selected[key]=(read_json(best[2]),paths,points)
    print(f'Eligible: {len(selected)} models; mass factors={args.mass_factors}; reference Yukawas={args.reference_yukawas}')
    if not selected: return 1
    if not args.run:
        print('Dry run only; use --run to compute the heavy-mass scan')
        return 0
    if not TARGET.is_file(): raise FileNotFoundError(TARGET)
    target=OscillationFitTarget.from_json(TARGET)
    out.mkdir(parents=True,exist_ok=True)
    output_file=out/'running_vs_fixed_mass_scan.json'
    settings={'scenario':args.scenario,'mass_factors':args.mass_factors,'reference_yukawas':args.reference_yukawas,'models':sorted(selected)}
    if output_file.is_file() and not args.force:
        old=read_json(output_file)
        if old.get('settings')!=settings: raise RuntimeError('Existing scan has different settings. Use --force')
        rows=old['points']
        print('Reused saved results; use --force to recompute')
    else:
        rows=[]
        for key,(base,paths,points) in selected.items():
            for y in args.reference_yukawas:
                pt=next(v for k,v in points.items() if np.isclose(y,k,rtol=1e-10,atol=0))
                lam=complex(pt['lambdaT3']['real'],pt['lambdaT3']['imag'])
                common=deform_common_scalar(base,lam)
                for factor in args.mass_factors:
                    row={'model':key,'reference_yukawa_target':y,'mass_factor':factor,'status':'Failed','error':'',
                         'lambdaT3_real':float(lam.real),'lambdaT3_imag':float(lam.imag)}
                    try:
                        modified=scale_heavy_mass(common,factor)
                        mu_uv,mu_f,mu_s,mu_low=check_scales(modified)
                        row.update(mu_uv_gev=mu_uv,mu_fermion_threshold_gev=mu_f,mu_scalar_threshold_gev=mu_s,
                                   mu_low_gev=mu_low)
                        running=evaluate(modified,paths,target)['prediction']
                        fixed=fixed_point(modified,paths[2],target)
                        for ob in OBS:
                            rg=float(running[ob]);fr=float(fixed[ob])
                            row['running_'+ob]=rg;row['fixed_'+ob]=fr
                            row['delta_'+ob]=rg-fr
                            row['relative_percent_'+ob]=100*(rg-fr)/fr if fr else float('nan')
                        if factor==1.:
                            row['max_abs_baseline_difference']=max(abs(float(pt['prediction'][ob])-float(running[ob])) for ob in OBS)
                        row['status']='Success'
                    except Exception as exc:
                        row['error']=f'{type(exc).__name__}: {exc}'
                    rows.append(row)
                    print(f'{key} y_ref={y:g} mass_factor={factor:g} {row["status"]} {row["error"]}',flush=True)
        output_file.write_text(json.dumps({'settings':settings,'points':rows},indent=2)+'\n',encoding='utf-8')
    with (out/'running_vs_fixed_mass_scan.csv').open('w',newline='',encoding='utf-8') as f:
        fields=list(dict.fromkeys(k for r in rows for k in r))
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader();writer.writerows(rows)
    save_plots(rows,out/'plots')
    good=sum(r['status']=='Success' for r in rows)
    baseline=[r['max_abs_baseline_difference'] for r in rows if r['status']=='Success' and 'max_abs_baseline_difference' in r]
    print(f'Successful {good}/{len(rows)}; max absolute baseline prediction difference={max(baseline) if baseline else "N/A"}')
    print('Output:',out)
    return 0 if good==len(rows) else 1

if __name__=='__main__':
    raise SystemExit(main())
