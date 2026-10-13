"""Three-panel T3 flavour-correlation scans; read-only reference benchmarks.

Run: python -m Numerical.orchestration.PlotCommonYukawaThreePanels --run --points 5
Resume: add --resume. A unique output folder protects prior results.
"""
from __future__ import annotations
import argparse
import csv
import json
import os
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from Numerical.orchestration.RunCommonParameterComparisonV2 import (
    PROJECT_ROOT, SCENARIO, OBS, candidate_for_target, models_and_cases, load,
)
from Numerical.core.RGBetaEvaluator import load_rgbeta_payload
from Numerical.running.FinalC5TrajectoryAdapter import FinalC5TrajectoryEvaluator
from Numerical.running.T3Trajectory import run_renormalisable_t3_trajectory, continue_with_weinberg_running
from Numerical.fitting.ScanCLI import build_uv_state_from_config
from Numerical.fitting.OscillationFit import OscillationFitTarget, evaluate_mass_matrix_fit
from physics.NeutrinoTrajectory import neutrino_mass_from_c5, charged_lepton_mass_basis_matrix
from physics.NeutrinoObservables import takagi_factorization

INPUT = PROJECT_ROOT / 'Reports/output/full/common_yukawa_least_squares' / SCENARIO / 'optimized_common_input.json'
OUTPUT = PROJECT_ROOT / 'Reports/output/full/common_yukawa_three_panel_scans' / SCENARIO
DIRECTIONS = ('S1', 'S2', 'S3')
COLORS = {'S1': '#187ac5', 'S2': '#e46b36', 'S3': '#268f53'}
FIELDS = ['model', 'direction', 't', 'status', 'error', 'chi2', *OBS, 'delta_cp_deg']


def deform(base, direction, t):
    p = deepcopy(base)
    y1 = p['base_state']['ordinary']['y1']
    y2 = p['base_state']['ordinary']['y2']
    a = np.asarray(y1['real'], dtype=float) + 1j*np.asarray(y1.get('imag', np.zeros((3,3))), dtype=float)
    b = np.asarray(y2['real'], dtype=float) + 1j*np.asarray(y2.get('imag', np.zeros((3,3))), dtype=float)
    # Flavour-sensitive, common-to-all-representations deformations.
    # S1: first lepton-flavour row; S2: second row; S3: complex phase on first row.
    if direction == 'S1':
        a[0, :] *= 1.0 + t
        b[0, :] *= 1.0 + t
    elif direction == 'S2':
        a[1, :] *= 1.0 + t
        b[1, :] *= 1.0 + t
    elif direction == 'S3':
        a[0, :] *= np.exp(1j*t)
        b[0, :] *= np.exp(1j*t)
    else:
        raise ValueError(direction)
    for spec, value in ((y1, a), (y2, b)):
        spec['real'] = value.real.tolist()
        spec['imag'] = value.imag.tolist()
    return p


def delta_cp_deg(mass_matrix):
    """PDG Dirac phase from rephasing-invariant PMNS quartets, in [0,360)."""
    u, _, residual = takagi_factorization(mass_matrix)
    if residual > 1e-7:
        raise ValueError(f'Takagi residual too large: {residual}')
    s13 = abs(u[0,2]); c13 = np.sqrt(max(0., 1-s13*s13))
    if c13 < 1e-10:
        return float('nan')
    s12 = abs(u[0,1])/c13; c12 = abs(u[0,0])/c13
    s23 = abs(u[1,2])/c13; c23 = abs(u[2,2])/c13
    denom = 2*s12*c12*s23*c23*s13
    jdenom = s12*c12*s23*c23*s13*c13*c13
    if min(denom, jdenom) < 1e-12:
        return float('nan')
    cosd = (abs(u[1,0])**2 - s12*s12*c23*c23 - c12*c12*s23*s23*s13*s13)/denom
    J = np.imag(u[0,0]*u[1,1]*np.conj(u[0,1])*np.conj(u[1,0]))
    sind = J/jdenom
    return float(np.degrees(np.arctan2(sind, np.clip(cosd, -1, 1))) % 360)


def evaluator(case, payload, target):
    data = case/'data'
    from Numerical.orchestration.ResumePipelineNumericalResults import REQUIRED
    paths = [data/name for name in REQUIRED]
    uv = load_rgbeta_payload(paths[0]); eft = load_rgbeta_payload(paths[1])
    c5_builder = FinalC5TrajectoryEvaluator(paths[2])
    scales = payload['scales']; numerical = payload.get('numerical', {})
    view = SimpleNamespace(raw={**payload, 'scan': {'bindings': {}}})
    state = build_uv_state_from_config(view, {}).validated()
    rtol = float(numerical.get('rtol', 1e-8)); atol = float(numerical.get('atol', 1e-11))
    trajectory = run_renormalisable_t3_trajectory(
        state, uv, eft,
        mu_fermion_threshold_gev=float(scales['mu_fermion_threshold_gev']),
        mu_scalar_threshold_gev=float(scales['mu_scalar_threshold_gev']),
        rtol=rtol, atol=atol)
    c5 = np.asarray(c5_builder({}, trajectory), dtype=complex)
    full = continue_with_weinberg_running(trajectory, c5, float(scales['mu_low_gev']), rtol=rtol, atol=atol)
    mass = neutrino_mass_from_c5(full.weinberg.K, vev_gev=float(numerical.get('vev_gev',246.22)))
    mass = charged_lepton_mass_basis_matrix(mass, full.weinberg.ye)
    fit = evaluate_mass_matrix_fit(mass, target)
    return float(fit.chi2), dict(zip(fit.observable_names, map(float, fit.prediction))), delta_cp_deg(mass)


def plot_model(model, rows, target, directory):
    center = dict(zip(target.observable_names, map(float,target.central_values)))
    sigma = dict(zip(target.observable_names, np.sqrt(np.diag(target.covariance))))
    panels = [
        ('sin2_theta13','delta_m3l_sq_ev2',r'$\sin^2\theta_{13}$',r'$\Delta m_{31}^2\;[\mathrm{eV}^2]$'),
        ('sin2_theta23','delta_m3l_sq_ev2',r'$\sin^2\theta_{23}$',r'$\Delta m_{31}^2\;[\mathrm{eV}^2]$'),
        ('sin2_theta13','delta_cp_deg',r'$\sin^2\theta_{13}$',r'$\delta_{\rm CP}\;[\mathrm{deg}]$'),
    ]
    fig, axes = plt.subplots(1,3,figsize=(15,4.6),constrained_layout=True)
    for ax,(x,y,xlabel,ylabel) in zip(axes,panels):
        if x in center:
            for n,alpha in ((3,.12),(1,.26)):
                ax.axvspan(center[x]-n*sigma[x],center[x]+n*sigma[x],color='#178dd1',alpha=alpha)
        if y in center:
            for n,alpha in ((3,.12),(1,.26)):
                ax.axhspan(center[y]-n*sigma[y],center[y]+n*sigma[y],color='#178dd1',alpha=alpha)
        for direction in DIRECTIONS:
            group=sorted((r for r in rows if r['direction']==direction and r['status']=='Success'),key=lambda r:float(r['t']))
            if not group: continue
            xx=np.asarray([float(r[x]) for r in group]); yy=np.asarray([float(r[y]) for r in group])
            if y=='delta_cp_deg':
                yy=np.rad2deg(np.unwrap(np.deg2rad(yy)))
                # Show physically equivalent phase branches without connecting across a 360-degree wrap.
                yy=yy-360*np.floor(np.median(yy)/360)
            ax.plot(xx,yy,'-o',lw=1.4,ms=3,color=COLORS[direction],label=direction)
            mid=min(range(len(group)),key=lambda i:abs(float(group[i]['t'])))
            ax.scatter([xx[mid]],[yy[mid]],marker='x',s=60,color=COLORS[direction],zorder=5)
        if y=='delta_cp_deg':
            ax.set_ylim(-10,370)
            ax.set_yticks([0,90,180,270,360])
            ax.text(.02,.02,'No experimental CP band supplied',transform=ax.transAxes,fontsize=7)
        ax.set_xlabel(xlabel); ax.set_ylabel(ylabel); ax.grid(alpha=.2)
    axes[0].legend(fontsize=8)
    fig.suptitle(model)
    fig.savefig(directory/f'{model}_three_panels.png',dpi=160)
    plt.close(fig)


def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--run',action='store_true')
    ap.add_argument('--resume',action='store_true')
    ap.add_argument('--points',type=int,default=5)
    ap.add_argument('--range',type=float,default=.15,help='S1/S2 fractional row deformation; S3 phase in radians')
    args=ap.parse_args(argv)
    if args.points<3 or args.points%2!=1 or not 0<args.range<1:
        ap.error('Require odd points >=3 and 0 < range < 1')
    configs,cases=models_and_cases()
    common=load(INPUT); ordinary=common['shared_ordinary']; source=load(configs[common['source_model']])
    source=deepcopy(source)
    for key in ('y1','y2','lambdaT3','MF','mS1Sq','mS2Sq'):
        source['base_state']['ordinary'][key]=deepcopy(ordinary[key])
    payloads={name:candidate_for_target(source,load(path)) for name,path in configs.items()}
    target=OscillationFitTarget.from_json(PROJECT_ROOT/'configs/nufit_6_1_ic24_no.json')
    if target.ordering!='NO' or tuple(target.observable_names)!=OBS:
        raise ValueError('Unexpected NuFIT target')
    print('Input:',INPUT,'\nOutput:',OUTPUT,'\nEvaluations:',len(payloads)*3*args.points,flush=True)
    if not args.run: return 0
    manifest={'input':str(INPUT),'points':args.points,'range':args.range,'directions':list(DIRECTIONS),'models':sorted(payloads)}
    if OUTPUT.exists():
        if not args.resume: raise FileExistsError(f'Existing output: {OUTPUT}; use --resume')
        if not (OUTPUT/'manifest.json').exists() or load(OUTPUT/'manifest.json')!=manifest:
            raise ValueError('Resume manifest mismatch')
    else:
        OUTPUT.mkdir(parents=True,exist_ok=False)
        (OUTPUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    csv_path=OUTPUT/'scan_points.csv'
    existing={}
    if csv_path.exists():
        with csv_path.open(newline='',encoding='utf-8') as f:
            for r in csv.DictReader(f):
                existing[(r['model'],r['direction'],round(float(r['t']),12))]=r
    rows=[]
    def checkpoint():
        temp=csv_path.with_suffix('.tmp')
        with temp.open('w',newline='',encoding='utf-8') as f:
            w=csv.DictWriter(f,fieldnames=FIELDS);w.writeheader();w.writerows(rows)
        os.replace(temp,csv_path)
    for name in sorted(payloads):
        subset=[]
        for direction in DIRECTIONS:
            for t in np.linspace(-args.range,args.range,args.points):
                key=(name,direction,round(float(t),12))
                if key in existing and existing[key]['status']=='Success':
                    r=existing[key]; print(name,direction,round(t,4),'Reused',flush=True)
                else:
                    r={'model':name,'direction':direction,'t':float(t),'status':'Success','error':'','chi2':'',**{k:'' for k in OBS},'delta_cp_deg':''}
                    try:
                        chi2,pred,cp=evaluator(cases[name],deform(payloads[name],direction,float(t)),target)
                        r.update(chi2=chi2,**pred,delta_cp_deg=cp)
                    except Exception as e:
                        r.update(status='Failed',error=f'{type(e).__name__}: {e}')
                    print(name,direction,round(t,4),r['status'],r['chi2'],flush=True)
                rows.append(r);subset.append(r);checkpoint()
        plot_model(name,subset,target,OUTPUT)
    summary={**manifest,'evaluations':len(rows),'successes':sum(r['status']=='Success' for r in rows),
             'regions':'Gaussian marginal 1sigma/3sigma only; no official 2D NuFIT contours',
             'cp_experimental':'Not plotted: CP experimental data not provided',
             'scan_definitions':{'S1':'multiply first lepton row of Y1,Y2 by (1+t)',
                                 'S2':'multiply second lepton row of Y1,Y2 by (1+t)',
                                 'S3':'multiply first lepton row of Y1,Y2 by exp(i*t)'}}
    (OUTPUT/'scan_summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    print('Finished:',OUTPUT,flush=True)
    return 0 if summary['successes']==len(rows) else 2

if __name__=='__main__':
    raise SystemExit(main())
