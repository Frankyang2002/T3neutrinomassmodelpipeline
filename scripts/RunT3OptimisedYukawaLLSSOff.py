#!/usr/bin/env python3
"""Benchmark-consistent LLSS-off ablation for saved optimised A(-4), D(-2) Yukawas.

Project path: scripts/RunT3OptimisedYukawaLLSSOff.py
Depends on the existing RefreshOptimisedYukawaThesisTables.py and
RunT3OriginalSimultaneousLLSSHeatmap.py; neither is modified.
UV -> integrate all fermions -> EFT1 -> integrate both scalars -> SMEFT.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from unittest.mock import patch

import numpy as np

from scripts import RefreshOptimisedYukawaThesisTables as refresh
from scripts import RunT3OriginalSimultaneousLLSSHeatmap as ablate
from scripts import RunT3SimultaneousEFT1All16 as nonmixed
from scripts import RunT3SimultaneousEFT1MixedBC as mixed
from RGE.running.intermediate import GeneralEFT1C5Mixing as beta

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'output/studies/general_scalar_extended/optimal_yukawa_thesis_refresh'
ORDER = refresh.ORDER
OBS = refresh.OBS
REFS = ('A_-4','D_-2')
TAGS = {'A_-4':'Am4','D_-2':'Dm2'}


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_source(source):
    manifest = read(source/'manifest.json')
    results = read(source/'results.json')
    if manifest.get('full_cases') is not True:
        raise ValueError('Saved source must contain all 16 U/E/W/S cases')
    if set(manifest.get('references',[])) != set(REFS):
        raise ValueError('Saved source must contain both donors')
    if set(manifest.get('models',[])) != set(ORDER):
        raise ValueError('Saved source must contain all 16 models')
    restriction=manifest.get('extra_yukawas')
    if restriction not in ('zero_and_uv_rhs_frozen','target-model-saved'):
        raise ValueError(f'Unrecognised conjugate-Yukawa policy: {restriction}')
    donors=refresh.donor_inputs()
    for donor in REFS:
        for filename, digest in manifest['donor_hashes'][donor].items():
            if donors[donor]['hashes'][filename] != digest:
                raise ValueError(f'{donor}: donor fit source has changed: {filename}')
        for model in ORDER:
            item=results.get(donor,{}).get(model,{})
            if item.get('status')!='Computed' or not all(k in item.get('cases',{}) for k in ('1111','1001')):
                raise ValueError(f'Missing saved full/frozen cases: {donor}/{model}')
    return manifest,results


def predict(case, observable):
    v=float(case['prediction'][observable])
    if not math.isfinite(v):
        raise ValueError('Nonfinite observable '+observable)
    return v


def compare(saved,off,isolator):
    rows=[]; audits=[]
    for ref in REFS:
        for model in ORDER:
            old=saved[ref][model]['cases']
            new=off.get(ref,{}).get(model,{})
            if new.get('status')!='Computed':
                raise ValueError(f'LLSS-off run failed: {ref}/{model}: {new.get("error", "missing")}')
            trial=new['cases']['1111']
            if not trial.get('eft1_on') or not trial.get('wilson_log_on'):
                raise ValueError(f'{ref}/{model}: counterfactual is not the full-running trajectory')
            old_b=ablate._decode_complex(old['1111']['c5_joint_scalar_bubble'])
            new_b=ablate._decode_complex(trial['c5_joint_scalar_bubble'])
            bubble=float(np.linalg.norm(old_b-new_b)/max(np.linalg.norm(old_b),1e-300))
            if bubble > 5e-5:
                raise ValueError(f'{ref}/{model}: finite scalar matching changed: {bubble:.3e}')
            # Strictly check the source was removed from a nontrivial EFT1 RHS.
            diag=isolator.per_model.get(model,{})
            if diag.get('rhs_calls',0)==0 or diag.get('max_rhs_source_norm',0)<=0:
                raise ValueError(f'{ref}/{model}: missing LLSS source in ablated RHS')
            row={'benchmark':ref,'model':model,'status':'VALID_EXPLORATORY',
                 'scalar_bubble_relative_residual':bubble}
            for obs in OBS:
                full=predict(old['1111'],obs)
                frozen=predict(old['1001'],obs)
                trial_val=predict(trial,obs)
                if frozen==0:raise ValueError(f'{ref}/{model}: zero 1001 observable {obs}')
                row.update({obs+'_FULL':full,obs+'_FROZEN':frozen,
                            obs+'_LLSS_OFF':trial_val,
                            obs+'_pct_running':100*(full-frozen)/frozen,
                            obs+'_pct_mixing_off':100*(trial_val-frozen)/frozen})
            rows.append(row)
            audits.append({'benchmark':ref,'model':model,'source_rhs_calls':diag['rhs_calls'],
                           'source_max_norm':diag['max_rhs_source_norm'],
                           'bubble_relative_residual':bubble})
    return rows,audits


def write_csv(path,rows):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def render(rows,out,args):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import SymLogNorm
    from matplotlib.transforms import blended_transform_factory
    for ref in REFS:
        subset={r['model']:r for r in rows if r['benchmark']==ref}
        mat=np.full((len(ORDER),4),np.nan)
        for i,model in enumerate(ORDER):
            r=subset[model]
            for j,obs in enumerate(OBS[:2]):
                mat[i,2*j]=r[obs+'_pct_running']
                mat[i,2*j+1]=r[obs+'_pct_mixing_off']
        tag=TAGS[ref]
        dest=out/tag
        dest.mkdir(parents=True,exist_ok=True)
        write_csv(dest/'four_column_masssq.csv',[
            {'model':model,
             'delta_m21_running_pct':mat[i,0], 'delta_m21_mixing_off_pct':mat[i,1],
             'delta_m3l_running_pct':mat[i,2], 'delta_m3l_mixing_off_pct':mat[i,3]}
            for i,model in enumerate(ORDER)])
        vmax=max(1.,float(np.nanmax(np.abs(mat))))
        norm=SymLogNorm(linthresh=1.,vmin=-vmax,vmax=vmax,base=10)
        cmap=plt.get_cmap('RdBu_r').copy()
        fig,ax=plt.subplots(figsize=(args.figure_width,args.row_height*len(ORDER)+2))
        im=ax.imshow(mat,aspect='auto',norm=norm,cmap=cmap)
        ax.set_yticks(range(len(ORDER)),[m.replace('_','(')+')' for m in ORDER],fontsize=args.tick_font_size)
        ax.set_xticks(range(4),['Running','Mixing off']*2,fontsize=args.tick_font_size)
        ax.set_ylabel('T3 model',fontsize=args.axis_font_size)
        ax.tick_params(axis='both',length=0)
        trans=blended_transform_factory(ax.transData,ax.transAxes)
        for x,label in ((.5,r'$\Delta m^2_{21}$'),(2.5,r'$\Delta m^2_{3\ell}$')):
            ax.text(x,1.025,label,transform=trans,ha='center',va='bottom',fontsize=args.top_font_size)
        ax.set_xticks(np.arange(-.5,4,1),minor=True)
        ax.set_yticks(np.arange(-.5,len(ORDER),1),minor=True)
        ax.grid(which='minor',color='white',linewidth=.7)
        ax.tick_params(which='minor',bottom=False,left=False)
        ax.axvline(1.5,color='black',linewidth=1.2)
        for y in (3,6,9,12):ax.axhline(y-.5,color='black',linewidth=1.2)
        for i in range(len(ORDER)):
            for j in range(4):
                v=mat[i,j]
                label=(f'{v:+.1e}' if abs(v)>=1e4 or abs(v)<.01 else
                       f'{v:+.1f}' if abs(v)>=10 else f'{v:+.2f}' if abs(v)>=1 else f'{v:+.3f}')
                rgba=cmap(norm(v))
                luminance=.2126*rgba[0]+.7152*rgba[1]+.0722*rgba[2]
                ax.text(j,i,label,ha='center',va='center',fontsize=args.cell_font_size,
                        color='white' if luminance<.48 else 'black')
        cb=fig.colorbar(im,ax=ax,fraction=.035,pad=.025)
        cb.ax.tick_params(labelsize=args.tick_font_size)
        cb.set_label('Signed change (%) — symmetric log colours',fontsize=args.axis_font_size)
        fig.subplots_adjust(left=.14,right=.88,bottom=.07,top=.94)
        for ext in ('png','pdf'):
            fig.savefig(dest/f'four_column_masssq.{ext}',dpi=300,bbox_inches='tight')
        plt.close(fig)
        print(ref,'plot:',dest/'four_column_masssq.png',flush=True)


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source-run-id',default='optimal_yukawa_thesis_heatmaps_3stage_v1')
    p.add_argument('--counterfactual-run-id',default='llss_off_optimised_both_v1')
    p.add_argument('--output',type=Path,default=BASE/'llss_off_comparison_both_v1')
    p.add_argument('--allow-stale-symbolic-source',action='store_true')
    p.add_argument('--resume',action='store_true')
    p.add_argument('--plot-only',action='store_true')
    p.add_argument('--figure-width',type=float,default=12.)
    p.add_argument('--row-height',type=float,default=.54)
    p.add_argument('--cell-font-size',type=float,default=12.)
    p.add_argument('--tick-font-size',type=float,default=12.)
    p.add_argument('--top-font-size',type=float,default=15.)
    p.add_argument('--axis-font-size',type=float,default=13.)
    args=p.parse_args(argv)
    for key in ('figure_width','row_height','cell_font_size','tick_font_size','top_font_size','axis_font_size'):
        if not 0<float(getattr(args,key))<100:p.error('Invalid positive '+key)
    if not all(c.isalnum() or c in '_-' for c in args.source_run_id+args.counterfactual_run_id):
        p.error('Unsafe run ID')
    source=BASE/args.source_run_id
    manifest,saved=verify_source(source)
    out=args.output.resolve()
    comparison=out/'all16_both_benchmarks.csv'
    if args.plot_only:
        with comparison.open(newline='',encoding='utf-8') as f:rows=list(csv.DictReader(f))
        for row in rows:
            for k,v in row.items():
                if k.endswith(('_pct_running','_pct_mixing_off')):row[k]=float(v)
        render(rows,out,args)
        return 0
    if out.exists() and any(out.iterdir()):
        p.error('Choose a fresh --output directory; use --plot-only for existing results')
    if args.counterfactual_run_id==args.source_run_id:
        p.error('Counterfactual run ID must differ from saved baseline')
    run_dir=BASE/args.counterfactual_run_id
    if run_dir.exists() and not args.resume:
        p.error('Counterfactual run already exists; use --resume or a new ID')
    if manifest['extra_yukawas']=='zero_and_uv_rhs_frozen':
        restriction=['--zero-extra-yukawas']
    else:restriction=[]
    options=['--run-id',args.counterfactual_run_id,'--source-run',manifest['source_run'],
             '--canonical-run-id',manifest['canonical_run'],'--full-cases',
             '--response-root',manifest['cached_quartic_context'],*restriction]
    if args.allow_stale_symbolic_source:options+=['--allow-stale-symbolic-source']
    if args.resume:options+=['--resume']
    class ExactRHSIsolator(ablate.LLSSIsolator):
        def beta(self, ctx, c):
            terms=self.original_beta(ctx,c)
            if not self.in_integration:
                return terms
            # The simultaneous integrator evaluates calculate_master_beta once
            # for the mu_F source-parity check before its first ODE RHS call.
            # Preserve that evaluation exactly; ablate only the ODE RHS.
            if self.skip_initial_parity:
                self.skip_initial_parity = False
                return terms
            src=ablate.source_in_original_basis(ctx,c)
            norm=float(np.linalg.norm(src))
            diag=self.per_model.setdefault(ctx.model,{
                'initial_source_norm':norm,'max_rhs_source_norm':0.,'rhs_calls':0})
            diag['rhs_calls']+=1
            diag['max_rhs_source_norm']=max(diag['max_rhs_source_norm'],norm)
            self.rhs_calls+=1
            self.maximum_source_norm=max(self.maximum_source_norm,norm)
            terms['scalar_pair']=terms['scalar_pair']-src
            return terms
    isolator=ExactRHSIsolator(beta.calculate_master_beta)
    orig_nonmixed=nonmixed.integrate_simultaneously
    orig_mixed=mixed.integrate_mixed
    def wrap_nonmixed(*a,**kw):
        with isolator.integration(eft_on=kw.get('eft_on',True)):
            return orig_nonmixed(*a,**kw)
    def wrap_mixed(*a,**kw):
        with isolator.integration(eft_on=kw.get('eft_on',True)):
            return orig_mixed(*a,**kw)
    # The donor fit check is deliberately inapplicable to an LLSS-off counterfactual.
    # The source donor inputs remain checked by verify_source() and refresh.donor_inputs().
    with patch.object(beta,'calculate_master_beta',isolator.beta), \
         patch.object(nonmixed,'integrate_simultaneously',wrap_nonmixed), \
         patch.object(mixed,'integrate_mixed',wrap_mixed), \
         patch.object(refresh,'validate_donor',lambda *a,**kw:None):
        rc=refresh.main(options)
    if rc!=0:
        raise RuntimeError(f'LLSS-off computation returned {rc}; inspect {run_dir}/results.json')
    new=read(run_dir/'results.json')
    # Source must not be altered during the counterfactual.
    verify_source(source)
    rows,audits=compare(saved,new,isolator)
    out.mkdir(parents=True,exist_ok=False)
    write_csv(comparison,rows)
    write_csv(out/'source_audit.csv',audits)
    (out/'provenance.json').write_text(json.dumps({
        'source_run':str(source),'source_manifest_sha256':sha(source/'manifest.json'),
        'source_results_sha256':sha(source/'results.json'),
        'counterfactual_run':str(run_dir),
        'donor_input_hashes':manifest['donor_hashes'],
        'extra_yukawas_policy':manifest['extra_yukawas'],
        'comparison_baseline':'same benchmark 1001 from saved source',
        'running':'100*(saved_1111-saved_1001)/saved_1001',
        'mixing_off':'100*(counterfactual_1111-saved_1001)/saved_1001',
        'ablation':'only HHSS x LLSS contribution to scalar_pair LLHH beta',
        'physics_status':'conditional exploratory, not certified hard one-loop matching',
        'donor_validation_in_counterfactual':'not applicable; source fit hashes validated',
    },indent=2)+'\n',encoding='utf-8')
    render(rows,out,args)
    return 0

if __name__=='__main__':
    raise SystemExit(main())
