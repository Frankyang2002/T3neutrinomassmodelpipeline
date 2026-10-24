"""Recreate thesis figures from *saved* UV/fermion-start T3 comparisons.

No RGE integrations, fits, or source modifications. Run at repository root:
    python -m Numerical.analysis.BuildExistingRGComparisonThesisFigures --source local

Use --source snapshot to regenerate the archived October 2026 GitHub extraction.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'Reports' / 'output' / 'thesis_assets'
SNAPSHOT = OUT / 'data' / 'verified_existing_comparison_snapshot.json'
CATEGORIES = [('A', 0, 3), ('B', 3, 6), ('C', 6, 9), ('D', 9, 12), ('E', 12, 16)]
OBS = ('delta_m21_sq_ev2', 'delta_m3l_sq_ev2', 'sin2_theta12', 'sin2_theta13', 'sin2_theta23')
OBSLABELS = [r'$\Delta m^2_{21}$',r'$\Delta m^2_{3\ell}$',r'$\sin^2\theta_{12}$',r'$\sin^2\theta_{13}$',r'$\sin^2\theta_{23}$']
LIMIT = 5.86  # Fit-screen criterion used in study, NOT a calibrated confidence interval.


def load(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def source_local():
    base = ROOT / 'output' / 'studies' / 'fermion_start_lambda_comparison'
    casenames = ('5e-9', '0p5')
    summaries = {}
    for case in casenames:
        summaries[case] = {}
        for mode in ('optimised', 'shared_a_m2'):
            path = base / case / 'comparison' / mode / 'summary.json'
            data = load(path)
            assert data['case'] == case and data['mode'] == mode and data['uv_running'] is False, path
            items = data['models']
            if len(items) != 16:
                raise ValueError(f'{path} must contain all 16 models')
            summaries[case][mode] = {row['model']: row for row in items}
    models = list(summaries['5e-9']['optimised'])
    for case in casenames:
        for mode in ('optimised','shared_a_m2'):
            if set(summaries[case][mode]) != set(models):
                raise ValueError(f'Missing models in {case}/{mode}')
    fitdata = {}
    for case in casenames:
        fitdata[case] = {
            mode: [float(summaries[case][mode][m]['chi2']) for m in models]
            for mode in ('optimised','shared_a_m2')}
        fitdata[case]['dm3_shared'] = [
            float(summaries[case]['shared_a_m2'][m]['prediction']['delta_m3l_sq_ev2']) for m in models]

    slpath = base / 'uv_sensitivity_comparison' / 'summary.csv'
    with slpath.open(newline='',encoding='utf-8-sig') as f:
        lines = list(csv.DictReader(f))
    uv = {}
    for mode in ('optimised','shared_a_m2'):
        byobs = {r['observable']:r for r in lines if r['mode']==mode}
        uv[mode] = {
            'uv':[float(byobs[x]['max_uv_abs_log_slope']) for x in OBS],
            'fermion':[float(byobs[x]['max_fermion_abs_log_slope']) for x in OBS]}
    abpath = ROOT / 'output/studies/lambda_5e-10_refit/diagnostics/a_m2_exact_seed_uv_ablation/ablation.json'
    ab = load(abpath)
    if set(ab['increment_percent']) != {'UV','EFT1','Weinberg'}:
        raise ValueError('Staged ablation does not match expected stages')
    staged = {s: [float(ab['increment_percent'][s][o]) for o in OBS]
              for s in ('UV','EFT1','Weinberg')}
    return {
        'provenance': {'source':'local_saved_outputs', 'note':'Read saved production studies; no fitting or running performed'},
        'models':models, 'fits':fitdata,
        'uv_sensitivity':{'observables':list(OBS),'modes':uv},
        'staged_ablation':{'model':ab['model'],'lambdaT3':ab['lambdaT3'],
                           'full_chi2':ab['stages']['111']['chi2'],
                           'observables':list(OBS),'stage_increments_percent':staged},
    }


def _style():
    plt.rcParams.update({
        'font.family':'DejaVu Sans', 'font.size':10,
        'axes.titlesize':12, 'axes.labelsize':10,
        'savefig.facecolor':'white', 'figure.facecolor':'white',
        'axes.spines.top':False, 'axes.spines.right':False,
        'pdf.fonttype':42, 'ps.fonttype':42,
    })


def _save(fig, name, directory):
    directory.mkdir(parents=True, exist_ok=True)
    for ext in ('png','pdf'):
        fig.savefig(directory / f'{name}.{ext}', dpi=210, bbox_inches='tight')
    plt.close(fig)


def _groups(ax):
    for index, (cat, a, b) in enumerate(CATEGORIES):
        if index % 2 == 0:
            ax.axvspan(a-.5,b-.5,color='0.97',zorder=-2)
        if b !=16:
            ax.axvline(b-.5,color='0.8',lw=.8,zorder=-1)
        ax.text((a+b-1)/2,1.018,cat,transform=ax.get_xaxis_transform(),
                ha='center',va='bottom',color='0.35',size=9)


def _labels(models):
    result=[]
    for m in models:
        a=m.split('_alpha_')[-1]
        result.append(a.replace('m','−').replace('p','+'))
    return result


def fig11(data, directory):
    modes=data['uv_sensitivity']['modes']
    fig, axs=plt.subplots(1,2,figsize=(12.3,4.0),sharey=True)
    x=np.arange(5)
    for ax,(mode,title) in zip(axs,[('optimised','Individually fitted Yukawas'),('shared_a_m2','Shared A(−2) Yukawas')]):
        u=np.asarray(modes[mode]['uv'],dtype=float)
        f=np.asarray(modes[mode]['fermion'],dtype=float)
        ax.semilogy(x-.08,u,'o-',markersize=6,lw=2,label='UV on (111)')
        ax.semilogy(x+.08,f,'s--',markersize=5,lw=1.7,label='UV frozen (011)')
        ax.axhline(1e-12,color='0.5',ls=':',lw=1)
        ax.set_xticks(x,OBSLABELS)
        ax.set_ylim(1e-16,20)
        ax.set_title(title)
        ax.grid(axis='y',alpha=.23)
        ax.legend(loc='upper right',fontsize=9)
    axs[0].set_ylabel(r'Maximum $\left|d\ln|O|/d\ln\gamma\right|$')
    fig.suptitle('Saved UV-on versus fermion-start rescaling sensitivities',y=1.06,fontweight='semibold')
    fig.text(.5,-.07,'Different fitted baselines and gamma spacings: descriptive comparison, not a same-seed ablation. '
             'Fermion-start values near $10^{-14}$ are numerical zero.',ha='center',fontsize=9)
    _save(fig,'fig11_saved_uv_vs_fermion_sensitivity',directory)


def fig12(data,directory):
    d=data['staged_ablation']
    vals=d['stage_increments_percent']
    fig,axs=plt.subplots(1,2,figsize=(12.5,4.1))
    for ax,ii,title in zip(axs,[[0,1],[2,3,4]],['Mass-squared splittings','Mixing observables']):
        x=np.arange(len(ii));w=.25
        for n,mode in enumerate(('UV','EFT1','Weinberg')):
            y=[vals[mode][i] for i in ii]
            ax.bar(x+(n-1)*w,y,w,label=mode)
        ax.axhline(0,color='0.3',lw=.8)
        ax.set_xticks(x,[OBSLABELS[i] for i in ii])
        ax.set_ylabel('Change from preceding stage (%)')
        ax.set_title(title)
        ax.grid(axis='y',alpha=.2)
        ax.legend(fontsize=9)
    fig.suptitle(r'Exact-seed staged RGE ablation: T3-A(−2), $\lambda_{T3}=5\times10^{-10}$',y=1.06,fontweight='semibold')
    fig.text(.5,-.08,r'Same input for 000→100→110→111; full-stage $\chi^2\approx1.57\times10^4$ (poor fit). '
             'Stage changes are successive, not additive.',ha='center',fontsize=9)
    _save(fig,'fig12_same_seed_staged_ablation',directory)


def fig13(data,directory):
    fits=data['fits'];models=data['models'];x=np.arange(len(models))
    fig,axs=plt.subplots(2,1,figsize=(13.3,7.5),sharex=True)
    for ax,mode,title in zip(axs,['optimised','shared_a_m2'],['Optimised separately for each representation','Shared A(−2) Yukawas within each coupling case']):
        for case,name,marker in [('5e-9',r'$\lambda_{T3}=5\times10^{-9}$ (E: $2.5\times10^{-9}$)','o'),('0p5',r'$\lambda_{T3}=0.5$ (E: 0.25)','s')]:
            y=np.log10(1+np.asarray(fits[case][mode]))
            ax.plot(x,y,marker+'-',label=name,markersize=5,lw=1.3)
        ax.axhline(np.log10(1+LIMIT),color='0.5',ls=':',lw=1.2,label=r'$\chi^2=5.86$ fit screen')
        ax.set_title(title)
        ax.set_ylabel(r'$\log_{10}(1+\chi^2)$')
        ax.grid(axis='y',alpha=.22)
        ax.legend(fontsize=9,loc='upper left')
        _groups(ax)
    axs[-1].set_xticks(x,_labels(models),rotation=0)
    axs[-1].set_xlabel('Representation hypercharge label, grouped by T3 class (A–E)')
    fig.suptitle('Saved low- and high-scalar benchmark fit quality (UV frozen)',y=1.02,fontweight='semibold')
    fig.text(.5,-.02,'These are fitted targets, not independent neutrino predictions. '
             'The threshold is a study fit screen, not a calibrated confidence interval.',ha='center',fontsize=9)
    _save(fig,'fig13_low_high_scalar_fit_quality',directory)


def fig14(data,directory):
    fits=data['fits'];models=data['models'];x=np.arange(len(models))
    fig,ax=plt.subplots(figsize=(13.1,4.4))
    for case,name,marker in [('5e-9',r'$\lambda_{T3}=5\times10^{-9}$ (E: $2.5\times10^{-9}$)','o'),('0p5',r'$\lambda_{T3}=0.5$ (E: 0.25)','s')]:
        v=np.asarray(fits[case]['dm3_shared'])
        ax.plot(x,v/v[0],marker+'-',lw=1.7,markersize=5,label=name)
    ax.axhline(1,color='0.35',ls='--',lw=1)
    _groups(ax)
    ax.set_xticks(x,_labels(models))
    ax.set_ylabel(r'$\Delta m^2_{3\ell}(\mathrm{model})/\Delta m^2_{3\ell}(\mathrm{A(−2)})$')
    ax.set_xlabel('Representation hypercharge label, grouped by T3 class (A–E)')
    ax.set_ylim(0.2,1.18)
    ax.grid(axis='y',alpha=.24)
    ax.legend(loc='lower left',fontsize=9)
    fig.suptitle('Representation-dependent mass splitting at shared Yukawas (UV frozen)',y=1.03,fontweight='semibold')
    fig.text(.5,-.055,'Yukawas shared across representations WITHIN each coupling case; donor Yukawas differ BETWEEN cases. '
             'Not a pure scalar-coupling ablation.',ha='center',fontsize=9)
    _save(fig,'fig14_shared_yukawa_representation_effect',directory)


def table(data,folder):
    folder.mkdir(parents=True,exist_ok=True)
    models=data['models'];fits=data['fits']
    with (folder/'table_existing_scalar_comparison.csv').open('w',encoding='utf-8',newline='') as f:
        w=csv.writer(f)
        w.writerow(['model','chi2_opt_5e-9','chi2_opt_0p5','chi2_shared_5e-9','chi2_shared_0p5',
                    'dm3_shared_5e-9_ev2','dm3_shared_0p5_ev2','dm3_ratio_5e-9_to_A','dm3_ratio_0p5_to_A'])
        for i,m in enumerate(models):
            w.writerow([m,fits['5e-9']['optimised'][i],fits['0p5']['optimised'][i],
                        fits['5e-9']['shared_a_m2'][i],fits['0p5']['shared_a_m2'][i],
                        fits['5e-9']['dm3_shared'][i],fits['0p5']['dm3_shared'][i],
                        fits['5e-9']['dm3_shared'][i]/fits['5e-9']['dm3_shared'][0],
                        fits['0p5']['dm3_shared'][i]/fits['0p5']['dm3_shared'][0]])
    with (folder/'table_uv_sensitivity_maxima.csv').open('w',encoding='utf-8',newline='') as f:
        w=csv.writer(f);w.writerow(['mode','observable','max_abs_log_slope_uv','max_abs_log_slope_fermion_start'])
        for mode,d in data['uv_sensitivity']['modes'].items():
            for i,o in enumerate(OBS):w.writerow([mode,o,d['uv'][i],d['fermion'][i]])
    with (folder/'table_staged_ablation.csv').open('w',encoding='utf-8',newline='') as f:
        w=csv.writer(f);w.writerow(['observable','UV_increment_pct','EFT1_increment_pct','Weinberg_increment_pct'])
        for i,o in enumerate(OBS):w.writerow([o]+[data['staged_ablation']['stage_increments_percent'][s][i]
                                              for s in ('UV','EFT1','Weinberg')])


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',choices=('local','snapshot'),default='local')
    a=p.parse_args(argv)
    d=load(SNAPSHOT) if a.source=='snapshot' else source_local()
    if len(d['models'])!=16 or d['uv_sensitivity']['observables']!=list(OBS):
        raise ValueError('Unexpected source data or observable order')
    for case in ('5e-9','0p5'):
        for column in ('optimised','shared_a_m2','dm3_shared'):
            v=d['fits'][case][column]
            if len(v)!=16 or not np.isfinite(v).all():raise ValueError('Invalid fit data')
    _style()
    dest=OUT/'figures'
    fig11(d,dest);fig12(d,dest);fig13(d,dest);fig14(d,dest)
    table(d,OUT/'tables')
    for case in ('5e-9','0p5'):
        n=sum(float(z)<=LIMIT for z in d['fits'][case]['optimised'])
        print(f'{case}: {n}/16 optimised benchmarks pass study chi2 <= {LIMIT}')
    print('Wrote 4 figures in PNG and PDF, plus 3 CSV tables to',OUT)
    return 0

if __name__=='__main__':
    raise SystemExit(main())
