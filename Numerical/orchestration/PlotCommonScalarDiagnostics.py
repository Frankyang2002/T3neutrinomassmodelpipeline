"""Plot existing common-lambda T3 scans; never rerun or modify physics outputs.

Usage from repository root:
 python -m Numerical.orchestration.PlotCommonScalarDiagnostics \
   --input Reports/output/full/common_scalar_comparison/smallY_smallL \
   --output Reports/output/full/common_scalar_diagnostics/smallY_smallL

Input should contain per-model folders with common_scalar_scan.json.
Plots use *reference-model* Yukawa targets as labels, not the actual Yukawa of each model.
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OBS = {
    'delta_m21_sq_ev2': r'$\Delta m^2_{21}\ [\mathrm{eV}^2]$',
    'delta_m3l_sq_ev2': r'$\Delta m^2_{3\ell}\ [\mathrm{eV}^2]$',
    'sin2_theta12': r'$\sin^2\theta_{12}$',
    'sin2_theta13': r'$\sin^2\theta_{13}$',
    'sin2_theta23': r'$\sin^2\theta_{23}$',
}

def load_scans(root: Path):
    scans=[]
    for f in sorted(root.glob('*/common_scalar_scan.json')):
        d=json.loads(f.read_text(encoding='utf-8-sig'))
        if not isinstance(d.get('points'), list):
            raise ValueError(f'Missing points: {f}')
        scans.append(d)
    if not scans:
        raise FileNotFoundError(f'No */common_scalar_scan.json found under {root}')
    target_ref=scans[0]['reference_yukawa_targets']
    common_ref=scans[0]['common_lambdaT3']
    for d in scans[1:]:
        if d['reference_yukawa_targets']!=target_ref or d['common_lambdaT3']!=common_ref:
            raise ValueError('Scans do not share identical scalar-coupling grid')
    return scans

def generate(root: Path, out: Path):
    scans=load_scans(root)
    out.mkdir(parents=True, exist_ok=True)
    targets=np.asarray(scans[0]['reference_yukawa_targets'],dtype=float)
    palette=plt.get_cmap('tab20')
    models=[d['model'] for d in scans]
    rows=[]
    for d in scans:
        for p in d['points']:
            r={'model':d['model'], 'reference_yukawa_target':p['reference_yukawa_target'],
               'common_lambda_real':p['lambdaT3']['real'], 'common_lambda_imag':p['lambdaT3']['imag'],
               'actual_yukawa_max_singular_value':p['actual_yukawa_max_singular_value'],
               'status':p['status'], 'chi2':p.get('chi2')}
            r.update(p.get('prediction') or {})
            rows.append(r)
    with (out/'evaluated_points.csv').open('w',newline='',encoding='utf-8') as f:
        keys=['model','reference_yukawa_target','common_lambda_real','common_lambda_imag',
              'actual_yukawa_max_singular_value','status','chi2',*OBS]
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)

    fig, axs=plt.subplots(2,3,figsize=(15,8.3))
    for j,(key,label) in enumerate(OBS.items()):
        ax=axs.flat[j]
        for i,d in enumerate(scans):
            x=[];y=[]
            for p in d['points']:
                if p['status']=='Success' and p['prediction'] is not None:
                    x.append(p['reference_yukawa_target']); y.append(p['prediction'][key])
            ax.plot(x,y,marker='o',markersize=3,lw=1.2,color=palette(i%20),label=models[i])
        ax.set_xscale('log');ax.set_xlabel('Reference Yukawa target (not model-specific $y$)')
        ax.set_ylabel(label);ax.grid(alpha=.25)
    ax=axs.flat[5]
    for i,d in enumerate(scans):
        x=np.asarray([p['reference_yukawa_target'] for p in d['points']],float)
        y=[p['actual_yukawa_max_singular_value'] for p in d['points']]
        ax.plot(x,y,marker='o',markersize=3,lw=1.2,color=palette(i%20))
    ax.axhline(1,ls='--',color='black',lw=1,label='Imposed input cutoff')
    ax.set_xscale('log');ax.set_yscale('log');ax.grid(alpha=.25)
    ax.set_xlabel('Reference Yukawa target');ax.set_ylabel('Actual max singular value of $Y_1,Y_2$')
    handles,labels=axs.flat[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='lower center',ncol=4,fontsize=7,frameon=False)
    fig.suptitle('Common scalar coupling: low-energy neutrino observables (16 T3 models)',fontsize=14)
    fig.tight_layout(rect=[0,0.16,1,.96]);fig.savefig(out/'common_scalar_observables.png',dpi=180);plt.close(fig)

    fig,axs=plt.subplots(1,3,figsize=(15,4.9))
    keys=('delta_m21_sq_ev2','delta_m3l_sq_ev2','sin2_theta13')
    for i,d in enumerate(scans):
        points={round(float(p['reference_yukawa_target']),10):p for p in d['points']}
        base=points.get(.01)
        if not base or base['status']!='Success':continue
        for j,key in enumerate(keys):
            xx=[];yy=[]
            b=float(base['prediction'][key]);
            for p in d['points']:
                if p['status']=='Success' and p['prediction'] is not None:
                    xx.append(p['reference_yukawa_target'])
                    yy.append(100*(p['prediction'][key]/b-1) if b!=0 else np.nan)
            axs[j].plot(xx,yy,marker='o',markersize=3,lw=1.2,color=palette(i%20),label=models[i])
    for j,key in enumerate(keys):
        axs[j].set_xscale('log');axs[j].set_xlabel('Reference Yukawa target')
        axs[j].set_ylabel('Change from $y_{ref}=10^{-2}$ [%]')
        axs[j].set_title(OBS[key]);axs[j].axhline(0,color='black',lw=.7);axs[j].grid(alpha=.25)
    fig.suptitle('Parameter-scan sensitivity, NOT isolated RGE corrections',fontsize=13)
    fig.tight_layout(rect=[0,.03,1,.92]);fig.savefig(out/'changes_vs_reference_point.png',dpi=180);plt.close(fig)

    n_success=sum(p['status']=='Success' for d in scans for p in d['points'])
    report={'model_count':len(scans),'scan_points_total':sum(len(d['points']) for d in scans),
            'successful_points':n_success,'skipped_or_failed_points':sum(len(d['points']) for d in scans)-n_success,
            'reference_targets':targets.tolist(), 'common_lambda':scans[0]['common_lambdaT3'],
            'caveat':'Point-to-point changes are not isolated RGE corrections. No gauge-vs-Yukawa beta decomposition or fixed-parameter control data are present.'}
    (out/'diagnostics_summary.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args();generate(a.input,a.output)
if __name__=='__main__':main()
