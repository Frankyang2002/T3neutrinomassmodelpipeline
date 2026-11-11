from pathlib import Path
import csv, math, argparse, json, hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

MODELS=['A_-4','A_-2','A_+0','B_-3','B_-1','B_+1','C_-3','C_-1','C_+1','D_-2','D_+0','D_+2','E_-4','E_-2','E_+0','E_+2']
METRICS=[('delta_m21_sq_ev2',r'$\Delta m^2_{21}$'),('delta_m3l_sq_ev2',r'$\Delta m^2_{3\ell}$'),('sin2_theta12',r'$\sin^2\theta_{12}$'),('sin2_theta13',r'$\sin^2\theta_{13}$'),('sin2_theta23',r'$\sin^2\theta_{23}$')]
STAGES=['UV','EFT1','SMEFT']

def fail(msg): raise SystemExit('ERROR: '+msg)
def pct(a,b):
    if not math.isfinite(a) or not math.isfinite(b) or b==0:fail('Invalid full/frozen observables')
    return 100*(a-b)/b

def old_records(path):
    if not path.is_file():fail('Missing locked source: '+str(path))
    out={}
    with path.open(newline='',encoding='utf-8-sig') as f:
        for r in csv.DictReader(f):
            if r['switch'] not in ('UV','EFT1','SMEFT') or r['observable'] not in dict(METRICS): continue
            key=(r['model'],r['switch'],r['observable'])
            if key in out:fail('Duplicate old result: '+str(key))
            try:a=float(r['full_value']);b=float(r['frozen_value']);d=float(r['relative_change_percent'])
            except ValueError:fail('Missing original value: '+str(key))
            if not math.isclose(pct(a,b),d,rel_tol=1e-8,abs_tol=1e-8):fail('Old percentage mismatch: '+str(key))
            out[key]=(a,b,d)
    if len(out)!=16*3*5:fail(f'Old running data incomplete: {len(out)}/240')
    return out

def combined_records(path,old):
    if path is None:return {}
    if not path.is_file():fail('Missing combined-EFT1 file: '+str(path))
    out={}
    with path.open(newline='',encoding='utf-8-sig') as f:
        reader=csv.DictReader(f)
        required={'model','observable','full_value','eft1_all_frozen_value'}
        if not required.issubset(set(reader.fieldnames or [])):fail('Combined EFT1 CSV must contain '+str(required))
        for r in reader:
            k=(r['model'],r['observable'])
            if k in out:fail('Duplicate EFT1 result: '+str(k))
            if k[0] not in MODELS or k[1] not in dict(METRICS):fail('Unknown model/observable: '+str(k))
            a=float(r['full_value']); b=float(r['eft1_all_frozen_value'])
            olda=old[(k[0],'EFT1',k[1])][0]
            if not math.isclose(a,olda,rel_tol=1e-8,abs_tol=1e-12):fail('FULL benchmark changed: '+str(k))
            out[k]=(a,b,pct(a,b))
    if len(out)!=80:fail(f'Combined EFT1 results incomplete: {len(out)}/80')
    return out

def build(args):
    root=Path(args.project_root).resolve()
    source=root/'Reports/output/thesis_qualitative/representation_running/csv/frozen_stage_observable_changes.csv'
    old=old_records(source)
    combined=combined_records(Path(args.combined_eft1).resolve() if args.combined_eft1 else None,old)
    if not combined and not args.allow_incomplete:fail('Combined EFT1 counterfactual missing. Provide --combined-eft1 from an actual 16-model rerun; never substitute the old EFT1 background-only or newer generalised outputs. Add --allow-incomplete to generate a diagnostic graph with missing EFT1 cells.')
    out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=True)
    data=[]
    for model in MODELS:
        for stage in STAGES:
            for metric,label in METRICS:
                a,b,d=(combined[(model,metric)] if stage=='EFT1' and combined else old[(model,stage,metric)] if stage!='EFT1' else (old[(model,stage,metric)][0],float('nan'),float('nan')))
                data.append({'model':model,'stage_off':stage,'observable':metric,'full_value':a,'stage_frozen_value':b,'relative_change_percent':d,'status':'RECOMPUTED_ALL_EFT1_FROZEN' if stage=='EFT1' and combined else 'NOT_COMPUTED' if stage=='EFT1' else 'LOCKED_OLD_DATA'})
    with (out/'running_all_stages_comparison.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
    for keyset,filename,title in [(METRICS[:2],'mass_squared_running_all_frozen.png','Signed percentage changes in mass-squared differences'),(METRICS[2:],'mixing_running_all_frozen.png','Signed percentage changes in mixing-angle observables')]:
        ncols=len(keyset)*3
        grid=np.array([[next(r['relative_change_percent'] for r in data if r['model']==model and r['stage_off']==stage and r['observable']==metric) for stage in STAGES for metric,_ in keyset] for model in MODELS])
        finite=grid[np.isfinite(grid)];lim=max(1.,np.max(np.abs(finite)))
        fig,ax=plt.subplots(figsize=(max(9,ncols*1.25),9.0))
        cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('lightgray')
        im=ax.imshow(np.ma.masked_invalid(grid),cmap=cmap,norm=TwoSlopeNorm(vmin=-lim,vcenter=0,vmax=lim),aspect='auto')
        ax.set_xticks(np.arange(ncols),[f'{stage}\n{label}' for stage in STAGES for _,label in keyset]);ax.set_yticks(np.arange(16),MODELS)
        ax.set_title(title+(' [INCOMPLETE: EFT1 NOT COMPUTED]' if not combined else ' [ALL PARAMETERS FROZEN BY STAGE]'),pad=17)
        ax.set_xlabel('Disabled running stage and oscillation observable');ax.set_ylabel('T3 model')
        for i in range(16):
            for j in range(ncols):
                v=grid[i,j];label='N/A' if not np.isfinite(v) else f'{v:.3g}'
                ax.text(j,i,label,ha='center',va='center',fontsize=8,color=('white' if np.isfinite(v) and abs(v)>.55*lim else 'black'))
        ax.set_xticks(np.arange(-.5,ncols,1),minor=True);ax.set_yticks(np.arange(-.5,16,1),minor=True);ax.grid(which='minor',color='black',linewidth=.5);ax.tick_params(which='minor',bottom=False,left=False)
        for z in [3,6,9,12]:ax.axhline(z-.5,color='black',linewidth=1.5)
        for z in range(len(keyset),ncols,len(keyset)):ax.axvline(z-.5,color='black',linewidth=1.5)
        cb=fig.colorbar(im,ax=ax);cb.set_label(r'$100\,(O_{\mathrm{full}}-O_{X\mathrm{-off}})/O_{X\mathrm{-off}}$ (%)')
        fig.tight_layout();fig.savefig(out/filename,dpi=230,bbox_inches='tight');plt.close(fig)
    manifest={'source':str(source.relative_to(root)),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'combined_EFT1_provided':bool(combined),'EFT1_definition':'freeze ALL EFT1 beta functions, renormalisable AND every dimension-five Wilson component (including Weinberg where present)','SMEFT_definition':'freeze ALL SMEFT running including Weinberg C5','UV_and_SMEFT_reused_from_locked_dataset':True,'numerical_physics_recomputed_by_this_script':False,'graph_status':'COMPLETE_WITH_RECOMPUTED_EFT1' if combined else 'INCOMPLETE_FOR_DIAGNOSTIC_ONLY'}
    (out/'provenance.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    print('Generated:',out);print('STATUS:',manifest['graph_status']);print('Fitted inputs and matching must be identical to original for new EFT1 rerun.')

if __name__=='__main__':
    p=argparse.ArgumentParser(description='Regenerate thesis stage-freeze heatmaps without mixing pipelines or substituting background-only for combined EFT1 freeze.')
    p.add_argument('--project-root',default='.')
    p.add_argument('--combined-eft1',help='CSV from freshly recomputed original 16-model EFT1-all-frozen counterfactuals; columns model,observable,full_value,eft1_all_frozen_value')
    p.add_argument('--output',default='Reports/output/thesis_qualitative/representation_running_v2_all_frozen')
    p.add_argument('--allow-incomplete',action='store_true',help='Render placeholder N/A for missing combined EFT1 results; not for thesis submission')
    build(p.parse_args())
