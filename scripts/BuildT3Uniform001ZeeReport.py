"""Reuse the existing Zee Gaussian-oval report with NEW uniform-0.01 T3 results.

Modes:
 - reference: CURRENTLY AVAILABLE, draws experimental ovals with all 16
   independently refitted new benchmark reference centres. No T1/T2/T3 paths
   or model-response covariance ovals are invented from a single point.
 - sensitivity: invokes UNMODIFIED scripts.PlotT3GeneralZeeOvals.main on genuine
   new sensitivity points.csv/summary.json calculated by
   scripts.RunT3Uniform001ZeeSensitivity.

Outputs are under Reports/output/thesis_general_5_3/zee_uniform001_ovals_v1/.
No existing report or source is modified.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
CAMPAIGNS=ROOT/'output/studies/thesis_qualitative/uniform_scalar_001/campaigns'
SENSITIVITY=ROOT/'output/studies/thesis_qualitative/uniform_scalar_001/zee_sensitivity'
REPORTS=ROOT/'Reports/output/thesis_general_5_3/zee_uniform001_ovals_v1'
QUARTIC=0.01


def load_reference_table(path: Path, summary_path: Path, expected_models: tuple[str,...]) -> dict:
    """Only accepted NEW reference fits; no old trajectories or new scan points."""
    meta=json.loads(summary_path.read_text(encoding='utf-8-sig'))
    if not (meta.get('stage')=='reference' and
            meta.get('schema')=='t3-thesis-uniform-001-reference-and-scan-v1' and
            meta.get('uniform_scalar_quartic_uv')==QUARTIC and
            meta.get('conjugate_yukawas_uv_zero') is True):
        raise ValueError('Not the saved uniform-0.01 reference campaign')
    summaries=meta.get('model_summaries',[])
    accepted={m['model'] for m in summaries if m.get('status')=='Accepted'}
    if accepted!=set(expected_models) or len(summaries)!=len(expected_models):
        raise ValueError('The 16 new reference fits are not all accepted')
    out={}
    with path.open(newline='',encoding='utf-8-sig') as file:
        for row in csv.DictReader(file):
            m=row['model']
            if m not in expected_models or m in out or row['stage']!='reference' or \
               row['status'] not in ('AcceptedApproxFit','AcceptedRescaledSeed'):
                raise ValueError('Missing/duplicate/nonaccepted new reference: '+m)
            for k in ('chi2','lambdaT3','sin2_theta12','delta_m21_sq_ev2'):
                if not math.isfinite(float(row[k])):
                    raise ValueError(f'Invalid {m} {k}')
            if not 0<=float(row['chi2'])<=5.86:
                raise ValueError(f'{m}: reference is not an accepted fit')
            if not math.isclose(float(row['lambdaT3']),0.05 if m.startswith('E_') else 0.1,
                                rel_tol=0,abs_tol=1e-13):
                raise ValueError(f'{m}: wrong reference lambdaT3')
            out[m]=row
    if set(out)!=set(expected_models):
        raise ValueError('Reference comparison.csv is incomplete')
    return out


def preview_grid(table: dict, nufit: dict, image_path: Path) -> None:
    """Reuse existing plot module experimental contours, palette, and model titles."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    from scripts import PlotT3GeneralZeeOvals as original

    x='sin2_theta12'; y='delta_m21_sq_ev2'
    r1=original.JOINT_RADII['68% joint']
    r3=original.JOINT_RADII['99.73% joint']
    juno=original.JUNO_2026
    # Exactly the same experiment-frame construction as the existing report:
    # original NuFIT 3sigma marginal spans + JUNO 1sigma + actual reference.
    values={x:[],y:[]}
    for exp,radius in ((nufit,r3),(juno,r1)):
        for var in (x,y):
            mu=float(exp['central_values'][var])*original.UNIT_SCALE[var]
            sigma=float(exp['one_sigma_errors'][var])*original.UNIT_SCALE[var]
            values[var].extend((mu-radius*sigma,mu+radius*sigma))
    for row in table.values():
        for var in (x,y):
            values[var].append(float(row[var])*original.UNIT_SCALE[var])
    bounds={}
    for var,arr in values.items():
        lo,hi=min(arr),max(arr)
        mid=.5*(lo+hi)
        width=hi-lo
        bounds[var]=(mid-.58*width,mid+.58*width)
    fig,axes=plt.subplots(4,4,figsize=(16.8,15.4))
    for ax,model in zip(axes.flat,original.MODEL_ORDER):
        ax.set_facecolor('white')
        ax.set_xlim(*bounds[x]);ax.set_ylim(*bounds[y])
        original._exp_ellipse(ax,nufit,x,y,r3,facecolor=original.PALETTE['nufit_3'],
                             edgecolor='#94B9D8',alpha=.38,lw=1.0,zorder=.2)
        original._exp_ellipse(ax,nufit,x,y,r1,facecolor=original.PALETTE['nufit_1'],
                             edgecolor='#28639B',alpha=.32,lw=1.3,zorder=.3)
        original._exp_ellipse(ax,juno,x,y,r1,facecolor=original.PALETTE['juno_1'],
                             edgecolor='#A14B22',alpha=.28,lw=1.35,zorder=.4)
        val=table[model]
        ax.plot(float(val[x])*original.UNIT_SCALE[x],
                float(val[y])*original.UNIT_SCALE[y],marker='o',linestyle='none',
                markersize=4.5,markerfacecolor='#289344',markeredgecolor='white',
                markeredgewidth=.6,zorder=5)
        ax.set(xlabel=original.LABELS[x],ylabel=original.LABELS[y])
        ax.xaxis.label.set_size(9.3);ax.yaxis.label.set_size(9.3)
        ax.set_title(original._model_title(model),fontsize=11.1,pad=7)
        ax.grid(alpha=.16,color='#5C6471',linewidth=.55,zorder=.1)
        ax.tick_params(labelsize=8.1,direction='in',length=3)
        for spine in ax.spines.values():
            spine.set_linewidth(.75);spine.set_color('#666666')
    fig.subplots_adjust(left=.075,right=.975,top=.855,bottom=.055,wspace=.42,hspace=.45)
    fig.text(.5,.984,'Overview: experimental Gaussian ovals + uniform-0.01 reference fits',
             ha='center',va='top',fontsize=16)
    fig.text(.5,.955,
             '16 newly refitted general-T3 models | BSM quartics = 0.01; Yconj(UV) = 0',
             ha='center',va='top',fontsize=10)
    handles=[
        Patch(facecolor=original.PALETTE['nufit_3'],edgecolor='#94B9D8',
              label='NuFIT 6.1 ~3σ (joint Gaussian)'),
        Patch(facecolor=original.PALETTE['nufit_1'],edgecolor='#28639B',
              label='NuFIT 6.1 ~1σ (joint Gaussian)'),
        Patch(facecolor=original.PALETTE['juno_1'],edgecolor='#A14B22',
              label='JUNO 2026 1σ (measured, joint Gaussian)'),
        Line2D([],[],color='#289344',marker='o',markersize=4.5,lw=0,
               label='New refitted UV=0.01 baseline'),
    ]
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,.93),
               ncol=4,fontsize=8.2,frameon=True,edgecolor='#CCCCCC',fancybox=False)
    image_path.parent.mkdir(parents=True,exist_ok=True)
    original._save_figure(fig,image_path,dpi=170)
    plt.close(fig)


def main(argv=None) -> int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mode',choices=('reference','sensitivity'),default='reference')
    p.add_argument('--reference-run-id',default='uniform001_all16_reference_v1')
    p.add_argument('--sensitivity-run-id',default='uniform001_zee_all16_v1')
    p.add_argument('--report-dir',type=Path,default=REPORTS)
    p.add_argument('--nufit',type=Path,default=ROOT/'configs/nufit_6_1_ic24_no.json')
    args=p.parse_args(argv)
    from scripts import PlotT3GeneralZeeOvals as original
    if args.mode=='reference':
        folder=CAMPAIGNS/args.reference_run_id
        tab=load_reference_table(folder/'reference_comparison.csv',
                                 folder/'reference_summary.json',original.MODEL_ORDER)
        nufit=original.load_nufit(args.nufit)
        out=args.report_dir/'model_grid_uniform001_reference_only.png'
        preview_grid(tab,nufit,out)
        args.report_dir.mkdir(parents=True,exist_ok=True)
        (args.report_dir/'reference_preview_provenance.json').write_text(
            json.dumps({'format':'t3-uniform001-reference-only-preview-v1',
                        'reference_campaign':args.reference_run_id,
                        'reference_accepted':len(tab),'extra_bsm_quartics':QUARTIC,
                        'uv_yukawa_conjugate_zero':True,
                        'reference_figure_only':True,
                        'new_T1_T2_T3_sensitivity_computed':False,
                        'no_model_response_ovals_without_gamma_trajectory':True,
                        'contour_code':'scripts.PlotT3GeneralZeeOvals._exp_ellipse',
                        'physics_matching_scope':'conditional one-scale, not sequential SMEFT'},
                        indent=2)+'\n',encoding='utf-8')
        print(f'Saved reference-only preview: {out}')
        return 0
    source=SENSITIVITY/args.sensitivity_run_id
    if not (source/'points.csv').is_file() or not (source/'summary.json').is_file():
        p.error('No new uniform-0.01 sensitivity results; run RunT3Uniform001ZeeSensitivity first')
    summ=json.loads((source/'summary.json').read_text(encoding='utf-8-sig'))
    if (summ.get('schema')!='t3-general-zee-compensated-sensitivity-uniform001-v1'
        or summ.get('uniform_quartic')!=QUARTIC
        or summ.get('conjugate_yukawas_uv_zero') is not True
        or len(summ.get('models',[]))!=16
        or any(s.get('status')!='Complete' for s in summ['models'])):
        p.error('Uniform-0.01 sensitivity campaign is incomplete or wrong provenance')
    return original.main(['--source',str(source),'--nufit',str(args.nufit),
                          '--output-dir',str(args.report_dir/'full_sensitivity')])


if __name__=='__main__':
    raise SystemExit(main())
