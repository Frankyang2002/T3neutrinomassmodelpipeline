"""Compute fixed-parameter controls for an existing common-lambda T3 scan.

Uses fitted UV configs, existing final-C5 matching expression and oscillation
observable extraction. The no-running control freezes ALL renormalised
parameters at their UV values through matching and to the observation scale,
while retaining the one-loop radiative C5 matching itself. This is a
well-defined numerical diagnostic, not a full fixed-order calculation.

PowerShell from the repository root:
  python -m Numerical.orchestration.studies.scalar_rge.RunCommonScalarRGControls --dry-run
  python -m Numerical.orchestration.studies.scalar_rge.RunCommonScalarRGControls --run --limit 2
  python -m Numerical.orchestration.studies.scalar_rge.RunCommonScalarRGControls --run

Results are separate from the original scan and do not overwrite it.
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
from dataclasses import replace
from types import SimpleNamespace

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from Numerical.orchestration.support.ResumePipelineNumericalResults import PROJECT_ROOT, OUTPUT_ROOT, REQUIRED, discover_cases
from Numerical.orchestration.studies.fitted_rge.RunFittedSensitivityComparison import best_candidate
from Numerical.fitting.ScanCLI import build_uv_state_from_config
from Numerical.fitting.OscillationFit import OscillationFitTarget, evaluate_mass_matrix_fit
from Numerical.running.FermionThresholdBoundary import build_intermediate_scalar_boundary
from Numerical.running.FinalC5TrajectoryAdapter import evaluate_final_c5_on_trajectory
from physics.NeutrinoTrajectory import neutrino_mass_from_c5, charged_lepton_mass_basis_matrix

OBS = ('delta_m21_sq_ev2', 'delta_m3l_sq_ev2', 'sin2_theta12', 'sin2_theta13', 'sin2_theta23')
TARGET = PROJECT_ROOT / 'configs' / 'nufit_6_1_ic24_no.json'


def load(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def uv_view(payload):
    data = dict(payload)
    data.setdefault('scan', {'bindings': {}})
    return SimpleNamespace(raw=data)


def fixed_point(payload, c5_path, target):
    """Evaluate matching with UV parameters frozen at all thresholds.

    Retain the same explicit matching scales, but do not solve any RGEs.
    The heavy mass basis and scalar masses are handled by the production
    final-C5 adapter rather than by a hand-written loop formula.
    """
    uv = build_uv_state_from_config(uv_view(payload), {}).validated()
    scales = payload['scales']
    mu_f = float(scales['mu_fermion_threshold_gev'])
    mu_s = float(scales['mu_scalar_threshold_gev'])
    # The production boundary validates that its input carries the matching-scale
    # label. For a frozen control, only relabel the scale; do NOT evolve couplings.
    # Dataclass replacement retains all physical parameters unchanged.
    uv_at_f = replace(uv, mu_gev=mu_f).validated()
    intermediate_at_f = build_intermediate_scalar_boundary(
        uv_at_f, matching_scale_gev=mu_f,
    )
    intermediate = replace(intermediate_at_f, mu_gev=mu_s).validated()
    # Adapter accesses these exact public trajectory attributes.
    frozen = SimpleNamespace(uv_threshold_state=uv_at_f,
                             eft1_threshold_state=intermediate,
                             mu_scalar_threshold_gev=mu_s)
    c5 = evaluate_final_c5_on_trajectory(c5_path, frozen)
    vev = float(payload.get('numerical', {}).get('vev_gev', 246.22))
    mnu = neutrino_mass_from_c5(c5, vev_gev=vev)
    mnu_flavour = charged_lepton_mass_basis_matrix(mnu, uv.sm.ye)
    fit = evaluate_mass_matrix_fit(mnu_flavour, target)
    return dict(zip(fit.observable_names, map(float, fit.prediction)))


def create_plots(rows, folder):
    folder.mkdir(parents=True, exist_ok=True)
    models = sorted({r['model'] for r in rows if r['status']=='Success'})
    if not models:
        print('No successful fixed controls: skipping plots.')
        return
    colors = plt.get_cmap('tab20')
    labels = {'delta_m21_sq_ev2':r'$\Delta m_{21}^2$',
              'delta_m3l_sq_ev2':r'$\Delta m_{3\ell}^2$',
              'sin2_theta12':r'$\sin^2\theta_{12}$',
              'sin2_theta13':r'$\sin^2\theta_{13}$',
              'sin2_theta23':r'$\sin^2\theta_{23}$'}
    for obs in OBS:
        fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(13, 5))
        for j, model in enumerate(models):
            pts = sorted([r for r in rows if r['model']==model and r['status']=='Success'],
                         key=lambda r:r['reference_yukawa_target'])
            if not pts: continue
            x = [p['reference_yukawa_target'] for p in pts]
            ax0.plot(x, [p['running_'+obs] for p in pts], '-', color=colors(j % 20), lw=1.4, label=model)
            ax0.plot(x, [p['fixed_'+obs] for p in pts], '--', color=colors(j % 20), lw=1.1)
            ax1.plot(x, [p['delta_'+obs] for p in pts], '-o', color=colors(j % 20), lw=1.2, ms=3)
        ax0.set(title='Solid: RGE-improved; dashed: frozen parameters', ylabel=labels[obs])
        ax1.set(title='RGE effect (running minus frozen)', ylabel=r'$\Delta X_{\rm RG}$')
        for ax in (ax0, ax1):
            ax.set_xscale('log'); ax.set_xlabel('Reference-model Yukawa target')
            ax.grid(alpha=.25)
        ax0.legend(fontsize=5, ncol=2, loc='best')
        fig.suptitle('Common $\\lambda_{T3}$, compensating Yukawas: '+obs)
        fig.tight_layout()
        fig.savefig(folder/(obs+'_running_vs_fixed.png'), dpi=180)
        plt.close(fig)
    # Actual Yukawa strength, not calibrated reference strength, varies by model.
    for obs in OBS:
        fig, ax = plt.subplots(figsize=(8,5))
        for j, model in enumerate(models):
            pts = sorted([r for r in rows if r['model']==model and r['status']=='Success'],
                         key=lambda r:r['actual_yukawa_max_singular_value'])
            if pts:
                ax.plot([p['actual_yukawa_max_singular_value'] for p in pts],
                        [p['delta_'+obs] for p in pts], '-o', lw=1.2, ms=3,
                        color=colors(j % 20), label=model)
        ax.set(xscale='log', xlabel='Actual maximal UV Yukawa singular value',
               ylabel=r'$\Delta X_{\rm RG}$', title='RGE correction vs actual Yukawa: '+obs)
        ax.grid(alpha=.25); ax.legend(fontsize=5, ncol=2)
        fig.tight_layout(); fig.savefig(folder/(obs+'_correction_vs_actual_yukawa.png'), dpi=180)
        plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--scenario', default='smallY_smallL')
    ap.add_argument('--run', action='store_true', help='Compute fixed controls')
    ap.add_argument('--dry-run', action='store_true', help='Check input availability only')
    ap.add_argument('--limit', type=int, default=0)
    ap.add_argument('--force', action='store_true')
    args = ap.parse_args(argv)
    if args.limit < 0: ap.error('--limit must be nonnegative')
    scan_root = PROJECT_ROOT/'Reports'/'output'/'full'/'common_scalar_comparison'/args.scenario
    out = PROJECT_ROOT/'Reports'/'output'/'full'/'common_scalar_rg_controls'/args.scenario
    cases = {key:path for scenario,key,path in discover_cases(OUTPUT_ROOT) if scenario==args.scenario}
    keys = sorted(p.name for p in scan_root.iterdir() if p.is_dir() and (p/'common_scalar_scan.json').is_file()) if scan_root.is_dir() else []
    keys = keys[:args.limit or None]
    if not keys: raise RuntimeError('No common_scalar_scan.json files found under '+str(scan_root))
    if not TARGET.is_file(): raise FileNotFoundError(TARGET)
    configs = {}
    for key in keys:
        fit = best_candidate(args.scenario,key)
        if key not in cases or fit is None: raise RuntimeError(f'Missing fitted configuration or model: {key}')
        matching = cases[key]/'data'/REQUIRED[2]
        if not matching.is_file(): raise FileNotFoundError(matching)
        configs[key] = (fit[2], matching)
    print('Ready:',len(keys),'model(s); uploaded scans are used as running predictions.')
    if args.dry_run or not args.run: return 0
    target = OscillationFitTarget.from_json(TARGET)
    out.mkdir(parents=True, exist_ok=True)
    rows=[]
    for key in keys:
        config, matching = configs[key]
        base=load(config)
        scan=load(scan_root/key/'common_scalar_scan.json')
        for point in scan['points']:
            y=float(point['reference_yukawa_target'])
            row={'model':key,'reference_yukawa_target':y,
                 'actual_yukawa_max_singular_value':point.get('actual_yukawa_max_singular_value'),
                 'lambdaT3_real':point['lambdaT3']['real'],
                 'lambdaT3_imag':point['lambdaT3']['imag'],
                 'status':'Failed','error':''}
            try:
                if point['status']!='Success': raise RuntimeError('Original running point not successful')
                # Use the same deformation function as the original scan.
                from Numerical.orchestration.studies.scalar_rge.RunCommonScalarComparison import deform_common_scalar
                new_lambda=complex(row['lambdaT3_real'], row['lambdaT3_imag'])
                modified=deform_common_scalar(base,new_lambda)
                fixed=fixed_point(modified, matching, target)
                running=point['prediction']
                for ob in OBS:
                    row['running_'+ob]=float(running[ob]); row['fixed_'+ob]=float(fixed[ob])
                    row['delta_'+ob]=float(running[ob]-fixed[ob])
                    row['relative_percent_'+ob]=(100*(running[ob]-fixed[ob])/fixed[ob]
                                                  if fixed[ob]!=0 else float('nan'))
                row['status']='Success'
            except Exception as exc:
                row['error']=f'{type(exc).__name__}: {exc}'
            rows.append(row)
            print(f"{key} y_ref={y:g} {row['status']} {row['error']}", flush=True)
    fields=list(dict.fromkeys(k for row in rows for k in row))
    with (out/'running_vs_fixed.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader();writer.writerows(rows)
    (out/'running_vs_fixed.json').write_text(json.dumps(rows,indent=2)+'\n',encoding='utf-8')
    create_plots(rows, out/'plots')
    good=sum(r['status']=='Success' for r in rows)
    print(f'Saved {good}/{len(rows)} fixed-reference comparisons to {out}')
    return 0 if good==len(rows) else 1

if __name__=='__main__': raise SystemExit(main())
