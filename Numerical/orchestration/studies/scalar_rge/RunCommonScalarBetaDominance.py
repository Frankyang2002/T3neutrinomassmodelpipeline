"""Gauge/Yukawa decomposition of T3 UV Yukawa beta functions.

Uses the existing fitted candidates and successful common-scalar scan points.
No modification of RGE equations, model fits, or GitHub files.

From the repository root:
  python -m Numerical.orchestration.studies.scalar_rge.RunCommonScalarBetaDominance --dry-run
  python -m Numerical.orchestration.studies.scalar_rge.RunCommonScalarBetaDominance --run --limit 2
  python -m Numerical.orchestration.studies.scalar_rge.RunCommonScalarBetaDominance --run

For each saved UV state, evaluate exported RGBeta beta_y1,beta_y2 with:
  B_full = B(g, lambda, Yukawa)
  B_g0   = B(0, lambda, Yukawa)
  B_l0   = B(0, 0, Yukawa)
  B_y0   = B(0, 0, 0)
Then gauge=B_full-B_g0, scalar=B_g0-B_l0, Yukawa=B_l0-B_y0,
residual=B_y0. This is an algebraic ablation decomposition of the
exported one-loop expressions, not a separate beta-function derivation.
Only the UV (above fermion threshold) is assessed.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from Numerical.orchestration.support.ResumePipelineNumericalResults import (
    PROJECT_ROOT, OUTPUT_ROOT, REQUIRED, discover_cases,
)
from Numerical.orchestration.studies.fitted_rge.RunFittedSensitivityComparison import best_candidate
from Numerical.orchestration.studies.scalar_rge.RunCommonScalarComparison import deform_common_scalar, _view
from Numerical.fitting.ScanCLI import build_uv_state_from_config
from Numerical.core.RGBetaEvaluator import (
    evaluate_inputform_expression, load_rgbeta_payload, state_environment,
)
from Numerical.running.UVRunner import run_uv_segment

DEFAULT_SOURCE = PROJECT_ROOT / 'Reports' / 'output' / 'full' / 'common_scalar_comparison'
DEFAULT_OUTPUT = PROJECT_ROOT / 'Reports' / 'output' / 'full' / 'common_scalar_beta_dominance'
BETAS = ('y1', 'y2')
GAUGE = ('gY', 'g2', 'g3')
YUKAWA = ('y1', 'y2', 'yu', 'yd', 'ye', 'h')


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def find_expression(payload, coupling):
    report = payload['report_betas']
    if coupling in report:
        return str(report[coupling]), coupling
    # Do not guess a key if the export uses different naming.
    raise KeyError(f'Missing exact report_betas[{coupling!r}] (available: {list(report)})')


def beta_parts(expr, environment):
    base = dict(environment)
    nog = dict(base)
    for key in GAUGE:
        if key in nog:
            nog[key] = 0.0
    nol = dict(nog)
    for key in list(nol):
        if key.startswith('lambda'):
            nol[key] = 0.0
    noy = dict(nol)
    for key in YUKAWA:
        if key in noy:
            noy[key] = np.zeros_like(noy[key], dtype=complex)
    full = np.asarray(evaluate_inputform_expression(expr, base), dtype=complex)
    without_g = np.asarray(evaluate_inputform_expression(expr, nog), dtype=complex)
    without_g_lambda = np.asarray(evaluate_inputform_expression(expr, nol), dtype=complex)
    without_all = np.asarray(evaluate_inputform_expression(expr, noy), dtype=complex)
    if full.shape != (3, 3):
        raise ValueError(f'Expected Yukawa beta 3x3 matrix; got {full.shape}')
    gauge = full - without_g
    scalar = without_g - without_g_lambda
    yukawa = without_g_lambda - without_all
    residual = without_all
    reconstruction = gauge + scalar + yukawa + residual
    error = float(np.linalg.norm(full-reconstruction, ord='fro'))
    if error > 1e-9 * max(1.0, float(np.linalg.norm(full, ord='fro'))):
        raise RuntimeError(f'Beta decomposition did not reconstruct full beta: {error}')
    return full, gauge, yukawa, scalar, residual, error


def norm_pair(values):
    return float(np.sqrt(sum(np.linalg.norm(v, ord='fro')**2 for v in values)))


def measure(payload, state):
    environment = state_environment(state)
    parts = [beta_parts(find_expression(payload, k)[0], environment) for k in BETAS]
    full, gauge, yukawa, scalar, residual = [norm_pair([p[i] for p in parts]) for i in range(5)]
    gauge_ratio = float(yukawa / gauge) if gauge > 0 else None
    max_y = max(float(np.linalg.svd(getattr(state, k), compute_uv=False)[0]) for k in BETAS)
    return {'mu_gev': float(state.mu_gev), 'actual_yukawa_max_singular_value': max_y,
            'beta_full_norm': full, 'beta_gauge_norm': gauge,
            'beta_yukawa_norm': yukawa, 'beta_scalar_norm': scalar,
            'beta_residual_norm': residual, 'ratio_yukawa_to_gauge': gauge_ratio,
            'max_reconstruction_error': max(p[5] for p in parts)}


def plot_results(rows, plotdir):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print('Matplotlib unavailable; CSV and JSON written.')
        return
    valid = [r for r in rows if r['status'] == 'Success']
    if not valid:
        return
    plotdir.mkdir(parents=True, exist_ok=True)
    targets = sorted({r['reference_yukawa_target'] for r in valid})
    cmap = plt.get_cmap('viridis')
    models = sorted({r['model'] for r in valid})
    for name in ('ratio_vs_scale', 'beta_norms_vs_scale'):
        fig, axes = plt.subplots(1, len(targets), figsize=(4.0 * len(targets), 4.8), sharey=(name == 'ratio_vs_scale'))
        for j, target in enumerate(targets):
            ax = axes[j]
            group = [r for r in valid if r['reference_yukawa_target'] == target]
            for i, model in enumerate(models):
                rs = sorted([r for r in group if r['model'] == model], key=lambda v: v['mu_gev'])
                if not rs: continue
                x = [r['mu_gev'] for r in rs]
                color = cmap(i/max(1,len(models)-1))
                if name == 'ratio_vs_scale':
                    ax.plot(x, [r['ratio_yukawa_to_gauge'] for r in rs], color=color, lw=1.2,
                            label=model if j == len(targets)-1 else None)
                else:
                    ax.plot(x, [r['beta_gauge_norm'] for r in rs], color=color, alpha=.7, lw=1.1)
                    ax.plot(x, [r['beta_yukawa_norm'] for r in rs], color=color, alpha=.7, ls='--', lw=1.1)
            ax.set_xscale('log'); ax.set_yscale('log'); ax.grid(alpha=.25)
            ax.set_xlabel(r'$\mu\ (\mathrm{GeV})$'); ax.set_title(rf'$y_{{\rm ref}}={target:g}$')
            if name == 'ratio_vs_scale': ax.axhline(1, color='black', ls=':', lw=1)
        axes[0].set_ylabel(r'$\|\beta_Y^{\rm Yukawa}\|_F/\|\beta_Y^{\rm gauge}\|_F$' if name == 'ratio_vs_scale' else r'$\|\beta_{Y}\|_F$ (one-loop numerator)')
        if name == 'ratio_vs_scale':
            handles, labels = axes[-1].get_legend_handles_labels()
            fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(.5,0.0), ncol=4, fontsize=6)
        else:
            fig.text(.5,.015,'Solid: gauge-dependent; dashed: pure Yukawa',ha='center',fontsize=9)
        fig.tight_layout(rect=(0,.10 if name == 'ratio_vs_scale' else .05,1,1))
        fig.savefig(plotdir/(name+'.png'), dpi=170)
        plt.close(fig)
    # Highest-scale UV input diagnostic: each model contributes one point per target.
    fig, ax = plt.subplots(figsize=(9,5.5))
    for i, model in enumerate(models):
        entries = []
        for target in targets:
            group = [r for r in valid if r['model']==model and r['reference_yukawa_target']==target]
            if group: entries.append(max(group,key=lambda r:r['mu_gev']))
        ax.plot([r['actual_yukawa_max_singular_value'] for r in entries],
                [r['ratio_yukawa_to_gauge'] for r in entries],'-o',lw=1,ms=3,label=model)
    ax.set_xscale('log');ax.set_yscale('log');ax.axhline(1,ls=':',c='black')
    ax.set(xlabel='Actual largest Yukawa singular value at UV input',
           ylabel=r'$\|\beta_Y^{\rm Yukawa}\|_F/\|\beta_Y^{\rm gauge}\|_F$')
    ax.grid(alpha=.25);ax.legend(fontsize=5,ncol=2,loc='best')
    fig.tight_layout();fig.savefig(plotdir/'uv_input_ratio_vs_yukawa.png',dpi=180);plt.close(fig)
    print('Wrote 3 plots to',plotdir)


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scenario',default='smallY_smallL')
    parser.add_argument('--source-root',type=Path,default=DEFAULT_SOURCE)
    parser.add_argument('--output-root',type=Path,default=DEFAULT_OUTPUT)
    parser.add_argument('--points',type=int,default=9,help='Number of logarithmic UV scales including endpoints')
    parser.add_argument('--limit',type=int,default=0,help='Number of models; zero means all')
    parser.add_argument('--run',action='store_true')
    parser.add_argument('--dry-run',action='store_true')
    args=parser.parse_args(argv)
    if args.points<2 or args.limit<0: parser.error('Require --points >= 2 and --limit >= 0')
    available={k:case for scenario,k,case in discover_cases(OUTPUT_ROOT) if scenario==args.scenario}
    scans=args.source_root/args.scenario
    selected=[]
    for key in sorted(available):
        scan_path=scans/key/'common_scalar_scan.json'
        fit=best_candidate(args.scenario,key)
        beta_path=available[key]/'data'/'uv_rgbeta_rge.json'
        if scan_path.is_file() and fit is not None and beta_path.is_file():
            selected.append((key,scan_path,fit[2],beta_path))
    selected=selected[:args.limit or None]
    print(f'Found {len(selected)} models with scans, fitted candidates and UV RGBeta artifacts.')
    if not selected:
        print('No eligible data; check source-root and local fitted/artifact directories.')
        return 1
    for key,scan_path,fit_path,beta_path in selected:
        scan=read_json(scan_path)
        n=sum(p.get('status')=='Success' for p in scan['points'])
        print(f'  {key}: {n} successful points')
    if not args.run:
        print('Dry run: use --run to evaluate beta decomposition.')
        return 0
    output=args.output_root/args.scenario
    output.mkdir(parents=True,exist_ok=True)
    rows=[];failures=[]
    for key,scan_path,fit_path,beta_path in selected:
        scan=read_json(scan_path);base=read_json(fit_path);payload=load_rgbeta_payload(beta_path)
        for point in scan['points']:
            if point.get('status')!='Success': continue
            target=float(point['reference_yukawa_target'])
            complex_lambda=complex(point['lambdaT3']['real'],point['lambdaT3']['imag'])
            try:
                deformed=deform_common_scalar(base,complex_lambda)
                state=build_uv_state_from_config(_view(deformed),{}).validated()
                mu_f=float(deformed['scales']['mu_fermion_threshold_gev'])
                scales=np.geomspace(float(state.mu_gev),mu_f,args.points)
                trajectory=run_uv_segment(state,payload,mu_f,save_scales_gev=scales,
                    rtol=float(deformed.get('numerical',{}).get('rtol',1e-8)),
                    atol=float(deformed.get('numerical',{}).get('atol',1e-11)))
                for i in range(trajectory.n_points):
                    sample=measure(payload,trajectory.state_at_index(i))
                    sample.update(model=key,reference_yukawa_target=target,status='Success')
                    rows.append(sample)
                print(f'{key} y_ref={target:g} Success ({trajectory.n_points} UV points)',flush=True)
            except Exception as exc:
                msg=f'{type(exc).__name__}: {exc}'
                failures.append(dict(model=key,reference_yukawa_target=target,status='Failed',error=msg))
                print(f'{key} y_ref={target:g} Failed {msg}',flush=True)
    if rows:
        with (output/'beta_decomposition.csv').open('w',newline='',encoding='utf-8') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    (output/'beta_decomposition.json').write_text(json.dumps({'rows':rows,'failures':failures,'decomposition':'gauge=full-g0; scalar=g0-g0lambda0; Yukawa=g0lambda0-g0lambda0Y0; residual=g0lambda0Y0','scope':'UV above fermion threshold only'},indent=2)+'\n',encoding='utf-8')
    plot_results(rows,output/'plots')
    print(f'Completed {len(rows)} saved UV states; failed scan points: {len(failures)}. Output: {output}')
    return 1 if failures else 0

if __name__=='__main__':
    raise SystemExit(main())
