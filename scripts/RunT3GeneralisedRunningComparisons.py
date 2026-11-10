#!/usr/bin/env python3
"""Recompute Figure 5.2/5.3 UV, EFT1, SMEFT stage-freeze effects.

The FULL outputs are read from the 32-case fixed-Yukawa thesis calculation.
For each model and reference, three NEW whole-pipeline counterfactuals disable:
  UV:       one-loop RGBeta UV evolution (constant UV state at MF)
  EFT1:     scalar RGBeta evolution AND full dimension-five Wilson evolution
  SMEFT:    one-loop SM+Weinberg running between MS and mu_low
All non-disabled stages and both threshold matchings still execute afresh.
No coefficient is fitted, manufactured or copied between representations.

This is an explicitly exploratory, not gauge-certified, comparison.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
from contextlib import ExitStack
from pathlib import Path
import shutil
import sys
import traceback
from unittest.mock import patch

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path: sys.path.insert(0,str(ROOT))

from Numerical.parallel_t3.complete.all16_rough_yukawa_fit import (
    BENCHMARK, FreshAll16Evaluator, all_models, read_json, save_json)
from Numerical.parallel_t3.rough_inputs import builder as builder_module
from Numerical.parallel_t3.threshold import preflight as preflight_module
from Numerical.parallel_t3.threshold import rough_numerics as numerics_module
from Numerical.parallel_t3.threshold import pipeline as pipeline_module
from Numerical.parallel_t3.wilson import workflow, rough_stream2
from Numerical.running import UVRunner, IntermediateScalarRunner
from Numerical.core.StateVector import pack_uv_state
from Numerical.running.IntermediateScalarStateVector import pack_intermediate_scalar_state
from scripts.RunT3All16ExploratoryFit import (
    corrected_quartic_unit_tensors,consistent_project_c5,
    exploratory_preflight,exploratory_rough_case)
from scripts.BuildT3GeneralisedThesisTables import load_reference

STAGES=('UV','EFT1','SMEFT')
KEYS=('sin2_theta12','sin2_theta13','sin2_theta23','dm21_sq_eV2','dm31_sq_eV2')


def _frozen_uv(state,payload,mu_final_gev,*,save_scales_gev=None,**kwargs):
    """Return actual UVRunningResult geometry with constant UV couplings."""
    y0,layout=pack_uv_state(state)
    scales=np.asarray(save_scales_gev if save_scales_gev is not None
        else [state.mu_gev,mu_final_gev],float)
    if not (len(scales)>=2 and np.isclose(scales[0],state.mu_gev)
            and np.isclose(scales[-1],mu_final_gev)):
        raise ValueError('Frozen UV output scales mismatch')
    return UVRunner.UVRunningResult(
        mu_initial_gev=float(state.mu_gev),mu_final_gev=float(mu_final_gev),
        layout=layout,mu_gev=scales,t=np.log(scales),
        y=np.repeat(y0[:,None],len(scales),axis=1),
        solver_success=True,solver_message='UV couplings intentionally held constant',
        nfev=0,njev=0,nlu=0)


def _frozen_eft1(state,payload,mu_final_gev,*,save_scales_gev=None,**kwargs):
    """Return actual IntermediateScalarRunningResult geometry with frozen couplings."""
    y0,layout=pack_intermediate_scalar_state(state)
    scales=np.asarray(save_scales_gev if save_scales_gev is not None
        else [state.mu_gev,mu_final_gev],float)
    if not (len(scales)>=2 and np.isclose(scales[0],state.mu_gev)
            and np.isclose(scales[-1],mu_final_gev)):
        raise ValueError('Frozen EFT1 output scales mismatch')
    return IntermediateScalarRunner.IntermediateScalarRunningResult(
        mu_initial_gev=float(state.mu_gev),mu_final_gev=float(mu_final_gev),
        layout=layout,mu_gev=scales,t=np.log(scales),
        y=np.repeat(y0[:,None],len(scales),axis=1),
        solver_success=True,solver_message='EFT1 scalar couplings intentionally held constant',
        nfev=0,njev=0,nlu=0)


def _frozen_wilson_beta(c,y,lam,gauge=(),*,variant='fonseca_eq22_six_term',**kwargs):
    """Seven exact zero beta channels, keeping tensor geometry and threshold code."""
    from Numerical.parallel_t3.wilson.adm import SYMMETRIC_TERMS
    return ({name:np.zeros_like(c,dtype=np.complex128) for name in SYMMETRIC_TERMS},
            'exactly_frozen_EFT1_Wilson_beta')


def _frozen_smeft(c5,b_npz,bmeta,*,mu_s,mu_low,**kwargs):
    """Identity propagation of the *matched* C5, not a new matching contribution."""
    from Numerical.parallel_t3.threshold.matching import symmetric3
    return {'c5_low':symmetric3(c5,'frozen SMEFT C5'),
            'nfev':0,'solver_success':True,'running_sm_parameters':{}}


def _cleanup(root,token):
    if not (isinstance(token,str) and token.startswith('fit16_')
            and len(token)==16 and token[6:].isalnum()):
        raise ValueError('Unsafe scratch token')
    for part in ('background','rough_inputs','rough_wilson','rough_results','fit16_adapted'):
        p=root/'output/parallel_t3'/part/token
        if p.is_dir() and p.name==token: shutil.rmtree(p)


def _relative(full,frozen):
    try:
        f=float(full); v=float(frozen)
    except (TypeError,ValueError):
        return None
    if not (math.isfinite(f) and math.isfinite(v)) or v==0:
        return None
    return 100*(f-v)/v


def _plot(rows,out,reference):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.colors import TwoSlopeNorm
    except ImportError:
        print('matplotlib not installed: CSV/JSON produced, figures skipped',flush=True)
        return
    models=list(all_models())
    stage_order=list(STAGES)
    figs=[('Figure_5_2_running_angles',KEYS[:3],r'Change in mixing angles (%)'),
          ('Figure_5_3_running_splittings',KEYS[3:],r'Change in mass-squared differences (%)')]
    for name,keys,title in figs:
        mat=np.full((len(models),len(keys)*len(stage_order)),np.nan)
        for i,model in enumerate(models):
            for j,key in enumerate(keys):
                for k,stage in enumerate(stage_order):
                    z=next((r for r in rows if r['reference']==reference and
                        r['model']==model and r['stage_frozen']==stage),None)
                    if z is not None and z.get(key+'_pct') is not None:
                        mat[i,j*3+k]=z[key+'_pct']
        fig,ax=plt.subplots(figsize=(max(9,len(keys)*3.5),8))
        finite=mat[np.isfinite(mat)]
        lim=max(float(np.quantile(abs(finite),0.95)) if finite.size else 1.,1e-5)
        norm=TwoSlopeNorm(vmin=-lim,vcenter=0,vmax=lim)
        im=ax.imshow(np.ma.masked_invalid(mat),cmap='coolwarm',norm=norm,aspect='auto')
        ax.set_yticks(np.arange(len(models)),[m.replace('_','(')+')' for m in models])
        labels=[f'{key.replace("sin2_theta","sin²θ").replace("dm21_sq_eV2","Δm²21").replace("dm31_sq_eV2","Δm²31")}\n{stage}'
                for key in keys for stage in STAGES]
        ax.set_xticks(np.arange(len(labels)),labels,rotation=45,ha='right',fontsize=8)
        ax.set_title(f'{title} — {reference} fixed Yukawas (exploratory)')
        ax.set_xticks(np.arange(-.5,mat.shape[1],1),minor=True)
        ax.set_yticks(np.arange(-.5,mat.shape[0],1),minor=True)
        ax.grid(which='minor',color='grey',alpha=.2)
        ax.tick_params(which='minor',bottom=False,left=False)
        fig.colorbar(im,ax=ax,label='100 × (full − stage-frozen) / stage-frozen')
        fig.tight_layout()
        for ext in ('png','pdf'):
            fig.savefig(out/f'{name}_{reference}.{ext}',dpi=220)
        plt.close(fig)


def _write_analysis(rows,dest,reference):
    lines=['# Section 5.2: individual running-stage effects', '',
       f'Reference Yukawa matrices: {reference}; no per-model reoptimisation.',
       'All cases rerun both threshold matchings and all non-frozen RGE stages.',
       'EFT1-off means both scalar RGBeta evolution AND full dim-5 Wilson evolution are frozen.',
       'Definition: delta_X(O) = 100 [O_full - O_X-off]/O_X-off.',
       'Undefined percentages (zero frozen denominator or failed calculation) are N/A.',
       'The results remain conditional exploratory EFT benchmarks.', '',
       '| Model | Stage off | Δm²21 (%) | Δm²31 (%) | sin²θ12 (%) | sin²θ13 (%) | sin²θ23 (%) |',
       '|---|---|---:|---:|---:|---:|---:|']
    for r in rows:
        if r['reference']!=reference:continue
        values=[r.get(k+'_pct') for k in ('dm21_sq_eV2','dm31_sq_eV2',*KEYS[:3])]
        lines.append('| '+' | '.join([r['model'],r['stage_frozen']]+[
            ('N/A' if v is None else f'{v:.5g}') for v in values])+' |')
    lines += ['', '## Numerical stage-effect summaries', '']
    for stage in STAGES:
        good=[r for r in rows if r['reference']==reference and r['stage_frozen']==stage and
              r['status']=='VALID_EXPLORATORY']
        lines.append(f'- {stage} frozen: {len(good)}/16 valid counterfactuals.')
        for key in KEYS:
            magnitudes=[abs(r[key+'_pct']) for r in good if r.get(key+'_pct') is not None]
            if magnitudes:
                winner=max((r for r in good if r.get(key+'_pct') is not None),
                    key=lambda r:abs(r[key+'_pct']))
                lines.append(f'  - {key}: median absolute change '
                             f'{float(np.median(magnitudes)):.5g}%, '
                             f'maximum {abs(winner[key+"_pct"]):.5g}% in {winner["model"]}.')
    (dest/f'Section_5_2_results_{reference}.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    tex=[r'\subsection*{Running-stage dependence: updated numerical comparison}',
        f'For the fixed {reference.replace("_","(")+ ")"} Yukawa benchmark, the',
        'UV, intermediate-EFT, and SMEFT running segments are disabled one at a time.',
        r'We define $\delta_X(O)=100(O_{\mathrm{full}}-O_{X\mathrm{-off}})/O_{X\mathrm{-off}}$.',
        'The individual counterfactuals recompute both threshold matching steps.',
        'The EFT1-off calculation freezes both the renormalisable scalar background and',
        'the coupled dimension-five Wilson coefficients.',
        r'\begin{table}[htbp]',r'\centering',r'\scriptsize',
        r'\begin{tabular}{llrrrrr}',r'\hline',
        r'Model & Disabled & $\Delta m^2_{21}$ & $\Delta m^2_{31}$ & $\sin^2\theta_{12}$ & $\sin^2\theta_{13}$ & $\sin^2\theta_{23}$ \\',r'\hline']
    for r in rows:
        if r['reference']!=reference:continue
        val=[r.get(k+'_pct') for k in ('dm21_sq_eV2','dm31_sq_eV2',*KEYS[:3])]
        tex.append(r['model'].replace('_','(')+')'+' & '+r['stage_frozen']+' & '+' & '.join(
            '--' if v is None else f'{v:.4g}' for v in val)+r' \\')
    tex += [r'\hline',r'\end{tabular}',
            r'\caption{Percentage changes from disabling one running stage while keeping the other stages active. Fixed-Yukawa, exploratory calculation.}',
            f'\\label{{tab:generalised_running_{reference.lower().replace("_","")}}}',
            r'\end{table}']
    (dest/f'Section_5_2_results_{reference}.tex').write_text('\n'.join(tex)+'\n',encoding='utf-8')


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--tables',type=Path,required=True,help='Completed 32-case thesis table directory')
    p.add_argument('--reference-results',nargs='+',type=Path,required=True)
    p.add_argument('--reference',choices=('D_-2','A_-4','both'),default='D_-2')
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--batch1',default='output/batch1_general_contract/batch1_native_v1')
    p.add_argument('--wilson-solver',choices=('RK2','DOP853'),default='RK2')
    p.add_argument('--wilson-steps',type=int,default=2)
    p.add_argument('--keep-trial-files',action='store_true')
    p.add_argument('--debug',action='store_true')
    args=p.parse_args(argv)
    root=ROOT.resolve(); tables=(root/args.tables).resolve(); output=(root/args.output).resolve()
    if not tables.is_dir() or not (tables/'results.json').is_file():
        p.error('Run the 32-case thesis tables generator first (--tables)')
    if output.exists() or (root/'output/parallel_t3').resolve() not in output.parents:
        p.error('--output must be a NEW directory below output/parallel_t3')
    baseline=read_json(tables/'results.json')
    provenance=read_json(tables/'provenance.json')
    if provenance.get('wilson_solver')!=args.wilson_solver or provenance.get('wilson_steps')!=args.wilson_steps:
        p.error('Wilson solver/steps differ from the baseline 32-case run')
    refs=('D_-2','A_-4') if args.reference=='both' else (args.reference,)
    refdirs=[(root/path).resolve() for path in args.reference_results]
    references={r:load_reference(r,refdirs) for r in refs}
    cfg=read_json(root/BENCHMARK)
    for r in refs:
        if str(references[r]['source'])!=provenance['reference_benchmarks'][r]['result_json']:
            p.error(f'{r}: selected source file does not match the 32-case baseline')
    output.mkdir(parents=True)
    summary={'format':'t3-generalised-stage-freeze-v1',
             'definition':'100*(full-stage_frozen)/stage_frozen',
             'reference':list(refs),'original_32case_tables':str(tables),
             'conditional_not_physical':True,'stage_frozen_options':list(STAGES),
             'eft1_freeze_includes_full_Wilson_running':True,
             'rows':[]}
    save_json(output/'provenance.json',summary)
    rows=[]
    for reference in refs:
        for model in all_models():
            original=baseline[reference].get(model,{})
            baseline_observables=original.get('observables') or {}
            if original.get('status')!='VALID_EXPLORATORY':
                print(f'WARNING {reference} {model}: no valid full-running baseline',flush=True)
            for stage in STAGES:
                print(f'RUN {reference} -> {model}, freeze {stage}',flush=True)
                trial_token=None
                record={'reference':reference,'model':model,'stage_frozen':stage,
                        'status':'BLOCKED','full_chi2':original.get('chi2')}
                try:
                    if original.get('status')!='VALID_EXPLORATORY':
                        raise RuntimeError('Full-running baseline for this model is missing')
                    evaluator=FreshAll16Evaluator(root,model,cfg,batch1=args.batch1,
                        gamma_factor=float(provenance['gamma_factor_assumed']),
                        mu_low=float(provenance['mu_low_gev']),vev=float(cfg['vev_GeV']),
                        reduced_potential=True,mass_bilinears=provenance['mass_bilinears'],
                        wilson_solver=args.wilson_solver,wilson_steps=args.wilson_steps,
                        keep_trial_files=True)
                    with ExitStack() as stack:
                        stack.enter_context(patch.object(builder_module,'quartic_unit_tensors',corrected_quartic_unit_tensors))
                        stack.enter_context(patch.object(workflow,'project_c5',consistent_project_c5))
                        stack.enter_context(patch.object(rough_stream2,'project_c5',consistent_project_c5))
                        stack.enter_context(patch.object(preflight_module,'preflight_pair',exploratory_preflight))
                        stack.enter_context(patch.object(numerics_module,'run_rough_case',exploratory_rough_case))
                        if stage=='UV':
                            stack.enter_context(patch.object(UVRunner,'run_uv_segment',_frozen_uv))
                        elif stage=='EFT1':
                            stack.enter_context(patch.object(IntermediateScalarRunner,'run_intermediate_scalar_segment',_frozen_eft1))
                            stack.enter_context(patch.object(workflow,'beta',_frozen_wilson_beta))
                        elif stage=='SMEFT':
                            stack.enter_context(patch.object(pipeline_module,'evolve_existing_smeft',_frozen_smeft))
                        trial=evaluator(references[reference]['x'].copy())
                        trial_token=trial.meta.get('scratch_token')
                    obs=trial.observables
                    record.update(status='VALID_EXPLORATORY',chi2=trial.chi2,
                        residuals=trial.residuals,observables=obs)
                    for key in KEYS:
                        record[key+'_full']=baseline_observables.get(key)
                        record[key+'_frozen']=obs.get(key)
                        record[key+'_pct']=_relative(baseline_observables.get(key),obs.get(key))
                    if not args.keep_trial_files: _cleanup(root,trial_token)
                except Exception as exc:
                    record['error']=f'{type(exc).__name__}: {exc}'
                    if args.debug: record['traceback']=traceback.format_exc(limit=16)
                    print('BLOCKED',reference,model,stage,record['error'],flush=True)
                rows.append(record)
                save_json(output/reference/model/f'{stage}.json',record)
                save_json(output/'results.json',{'rows':rows,'definition':summary['definition']})
                print('DONE',reference,model,stage,record['status'],flush=True)
    columns=['reference','model','stage_frozen','status','chi2','full_chi2']+[
        f'{key}_{part}' for key in KEYS for part in ('full','frozen','pct')]+['error']
    with (output/'running_comparisons.csv').open('w',encoding='utf-8',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=columns,extrasaction='ignore')
        writer.writeheader();writer.writerows(rows)
    for reference in refs:
        _write_analysis(rows,output,reference)
        _plot(rows,output,reference)
    completed=sum(r['status']=='VALID_EXPLORATORY' for r in rows)
    save_json(output/'summary.json',{'requested':len(rows),'completed':completed,
        'failed':len(rows)-completed,'reference':list(refs),
        'conditional_not_physical':True,'definition':summary['definition']})
    print(f'Running comparisons: {completed}/{len(rows)} cases, output={output}',flush=True)
    return 0 if completed==len(rows) else 1

if __name__=='__main__':
    raise SystemExit(main())
