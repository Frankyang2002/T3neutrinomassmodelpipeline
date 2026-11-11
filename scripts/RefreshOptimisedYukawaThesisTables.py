"""Refresh 16-model thesis tables with A(-4)/D(-2) least-squares Yukawas.

Uses the CURRENT coupled EFT1 runners, not the older one-scale report paths.
UV -> all fermions -> EFT1 -> all scalars -> SMEFT.
Conditional/exploratory: provisional quartic basis, finite hard matching omitted.
Never changes fit inputs, symbolic exports, historical results or GitHub.

From repo root:
 python -X utf8 -m scripts.RefreshOptimisedYukawaThesisTables --check-only
 python -X utf8 -m scripts.RefreshOptimisedYukawaThesisTables --run-id optimal_yukawa_thesis_refresh_v1
 python -X utf8 -m scripts.RefreshOptimisedYukawaThesisTables --run-id optimal_yukawa_thesis_refresh_v1 --full-cases
 python -X utf8 -m scripts.RefreshOptimisedYukawaThesisTables --run-id optimal_yukawa_thesis_refresh_v1 --render-only

By default evaluates only case 1111, sufficient for masses/mixing tables.
--full-cases evaluates all 16 U/E/W/S switches for stage-freeze heatmaps.
Each model retains its own fitted scalar/SM and extra-Yukawa coordinates;
ONLY the real 3x3 y1 and y2 matrices are transferred from the donor.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import shutil
import tempfile
import traceback

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
STUDIES = ROOT / 'output/studies/general_scalar_extended'
FIT_ROOT = STUDIES / 'coupled_yukawa_least_squares'
DEFAULT_SOURCE = STUDIES / 'thesis_general_all16/thesis16_merged_v1'
DEFAULT_OUT = STUDIES / 'optimal_yukawa_thesis_refresh'
DONOR_DIR = {
    'A_-4': FIT_ROOT / 'ls_A4_target001_v1/A_m4',
    'D_-2': FIT_ROOT / 'ls_D2_target001_v1/D_m2',
}
ORDER = ('A_-4','A_-2','A_+0','B_-3','B_-1','B_+1',
         'C_-3','C_-1','C_+1','D_-2','D_+0','D_+2',
         'E_-4','E_-2','E_+0','E_+2')
OBS = ('delta_m21_sq_ev2','delta_m3l_sq_ev2',
       'sin2_theta12','sin2_theta13','sin2_theta23')
SCHEMA = 't3-optimised-coupled-eft1-thesis-refresh-v1'


def read(path):
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(path)
    return json.loads(path.read_text(encoding='utf-8-sig'))


def write(path, content):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name+'.tmp')
    temp.write_text(json.dumps(content,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    temp.replace(path)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def matrix(obj):
    real = np.asarray(obj['real'],dtype=float)
    imag = np.asarray(obj['imag'],dtype=float)
    z = real+1j*imag
    if z.shape!=(3,3) or not np.isfinite(z).all():
        raise ValueError('Expected finite complex 3x3 matrix')
    return z


def label(name):
    fam,a=name.split('_')
    n=int(a)
    return f'{fam}({n:+d})' if n != 0 else f'{fam}(0)'


def texmodel(name):
    return label(name).replace('(',r'($').replace(')',r'$)')


def fmt(x):
    if x is None:return '--'
    v=float(x)
    if not math.isfinite(v):return '--'
    if v==0:return '0'
    if abs(v)<1e-4 or abs(v)>=1e4:
        power=int(math.floor(math.log10(abs(v))))
        mantissa=v/(10.0**power)
        return rf'${mantissa:.4g}\\times10^{{{power}}}$'.replace('\\\\','\\')
    return f'{v:.7g}'


def table(path,caption,ident,cols,units,rows):
    n=len(cols)
    out=[r'\begin{table}[p]',r'\centering',rf'\caption{{{caption}}}',
         rf'\label{{{ident}}}',r'\begingroup',r'\footnotesize',
         r'\setlength{\tabcolsep}{3.2pt}',
         r'\begin{tabular}{l'+'r'*(n-1)+'}',r'\hline',
         ' & '.join(cols)+r' \\',
         ' & '.join(units)+r' \\',r'\hline']
    for row in rows:out.append(' & '.join(row)+r' \\')
    out += [r'\hline',r'\end{tabular}',r'\endgroup',r'\end{table}','']
    path.write_text('\n'.join(out),encoding='utf-8')


def donor_inputs():
    result={}
    for model,folder in DONOR_DIR.items():
        best=read(folder/'best.json')
        cfg=read(folder/'best_fitted_config.json')
        extra=read(folder/'best_general_uv_overrides.json')
        if best.get('model')!=model or best.get('status')!='Computed':
            raise ValueError(f'{model}: not a valid completed donor fit')
        vector=np.asarray(best['yukawas_flat'],dtype=float)
        if vector.shape!=(18,) or not np.isfinite(vector).all():
            raise ValueError(f'{model}: invalid 18 real Yukawa values')
        y1,y2=vector[:9].reshape(3,3),vector[9:].reshape(3,3)
        for name,arr in (('y1',y1),('y2',y2)):
            if not np.allclose(matrix(cfg['base_state']['ordinary'][name]),arr,rtol=0,atol=1e-18):
                raise ValueError(f'{model}: config does not reproduce best.json {name}')
            if not np.allclose(np.asarray([[complex(v['real'],v['imag']) for v in row]
                                        for row in extra[name]]),arr,rtol=0,atol=1e-18):
                raise ValueError(f'{model}: overrides do not reproduce best.json {name}')
        result[model]={'best':best,'config':cfg,'overrides':extra,
                       'y1':y1,'y2':y2,'paths':{n:str(folder/n) for n in
                         ('best.json','best_fitted_config.json','best_general_uv_overrides.json')},
                       'hashes':{n:sha(folder/n) for n in
                         ('best.json','best_fitted_config.json','best_general_uv_overrides.json')}}
    return result


def extract_row(case,ref,model):
    pred=case['prediction']
    if any(key not in pred for key in OBS):
        raise ValueError(f'{ref}->{model}: missing observables')
    masses=case['masses_ev_ascending']
    if len(masses)!=3 or any(float(x)<0 or not math.isfinite(float(x)) for x in masses):
        raise ValueError(f'{ref}->{model}: invalid masses')
    return {'reference':ref,'model':model,'status':'Computed',
            'chi2':float(case['chi2']),
            'm1_ev':float(masses[0]),'m2_ev':float(masses[1]),'m3_ev':float(masses[2]),
            **{name:float(pred[name]) for name in OBS},
            'c5_low':case['c5_at_low_scale']}


def validate_donor(row,donor,*,chi2_rtol=0.03):
    """Reject a reference that is not actually reproduced by the selected runner."""
    reference=donor['best']
    predicted=reference['prediction_1111']
    for name in OBS:
        x=row[name];y=float(predicted[name])
        if not math.isclose(x,y,rel_tol=2e-4,abs_tol=2e-9):
            raise RuntimeError(f'{row["model"]}: donor baseline not reproduced: {name} {x} vs {y}')
    for key,expected in zip(('m1_ev','m2_ev','m3_ev'),reference['masses_ev']):
        if not math.isclose(row[key],float(expected),rel_tol=2e-4,abs_tol=1e-9):
            raise RuntimeError(f'{row["model"]}: donor mass baseline not reproduced: {key}')
    # Chi2 close to zero is extremely sensitive to integrator tolerances; do not
    # use relative fit-chi2 comparison when both are effectively zero.
    if row['chi2']>1e-5:
        raise RuntimeError(f'{row["model"]}: donor chi2 unexpectedly large {row["chi2"]}')


def recompute_epsilon(rows):
    indexed={(r['reference'],r['model']):r for r in rows if r.get('status')=='Computed'}
    for ref in DONOR_DIR:
        donor=indexed.get((ref,ref))
        if donor is None:
            continue
        reference=matrix(donor['c5_low'])
        d=np.vdot(reference,reference)
        if d.real<=0:raise ValueError('Zero donor C5 matrix')
        for target in ORDER:
            r=indexed.get((ref,target))
            if r is None:continue
            z=matrix(r['c5_low'])
            k=np.vdot(reference,z)/d
            norm=np.linalg.norm(z)
            r['k_re']=float(k.real);r['k_im']=float(k.imag)
            r['epsilon_C5']=float(np.linalg.norm(z-k*reference)/norm) if norm else None


def render(rows,out):
    from math import degrees,asin,sqrt
    report=out/'tables'
    report.mkdir(parents=True,exist_ok=True)
    for ref in ('D_-2','A_-4'):
        selected={r['model']:r for r in rows if r['reference']==ref}
        mass=[];mix=[]
        for model in ORDER:
            r=selected.get(model)
            m=texmodel(model)
            if not r or r['status']!='Computed':
                mass.append([m]+['--']*7)
                mix.append([m]+['--']*7)
                continue
            mass.append([m,*[fmt(r[k]) for k in ('m1_ev','m2_ev','m3_ev')],
                fmt(r['delta_m21_sq_ev2']/1e-5),fmt(r['delta_m3l_sq_ev2']/1e-3),
                fmt(r['chi2']),fmt(r.get('epsilon_C5'))])
            angles=[degrees(asin(sqrt(min(1,max(0,r[k]))))) for k in
                   ('sin2_theta12','sin2_theta13','sin2_theta23')]
            mix.append([m,*[fmt(r[k]) for k in ('sin2_theta12','sin2_theta13','sin2_theta23')],
                        *[fmt(x) for x in angles],fmt(r['chi2'])])
        tag='Dm2' if ref=='D_-2' else 'Am4'
        table(report/f'{tag}_masses.tex',
              rf'T3 neutrino masses with the optimised {label(ref)} Yukawa matrices. The matrices are held fixed across all models; $\epsilon_{{C_5}}$ compares the low-scale Wilson matrix with the {label(ref)} reference.',
              f'tab:appC_{tag}_masses',
              ['Model','$m_1$','$m_2$','$m_3$',r'$\Delta m^2_{21}$',
               r'$\Delta m^2_{3\ell}$',r'$\chi^2$',r'$\epsilon_{C_5}$'],
              ['',r'\multicolumn{3}{c}{(eV)}',r'$(10^{-5}\,\mathrm{eV}^2)$',
               r'$(10^{-3}\,\mathrm{eV}^2)$','',''],mass)
        table(report/f'{tag}_mixing.tex',
              rf'T3 neutrino mixing observables with the optimised {label(ref)} Yukawa matrices, held fixed across the models.',
              f'tab:appC_{tag}_mixing',
              ['Model',r'$\sin^2\theta_{12}$',r'$\sin^2\theta_{13}$',
               r'$\sin^2\theta_{23}$',r'$\theta_{12}$',r'$\theta_{13}$',
               r'$\theta_{23}$',r'$\chi^2$'],
              ['','','','',r'$(^\circ)$',r'$(^\circ)$',r'$(^\circ)$',''],mix)
    csv_path=out/'observables.csv'
    columns=('reference','model','status','m1_ev','m2_ev','m3_ev',*OBS,'chi2','epsilon_C5','k_re','k_im','error')
    with csv_path.open('w',encoding='utf-8',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=columns,extrasaction='ignore')
        writer.writeheader();writer.writerows(rows)
    print('Tables:',report,flush=True)


def plot_heatmaps(results,out,target):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    sigmas=np.sqrt(np.diag(np.asarray(target.covariance,dtype=float)))
    if len(sigmas)!=5 or np.any(sigmas<=0):
        raise ValueError('Invalid fitting target covariance')
    cols=('0111','1011','1101','1110')
    names=('UV off','EFT1 off','Wilson off','SMEFT off')
    folder=out/'heatmaps';folder.mkdir(parents=True,exist_ok=True)
    for ref in ('D_-2','A_-4'):
        for k,name in enumerate(OBS):
            arr=np.full((len(ORDER),4),np.nan)
            for i,model in enumerate(ORDER):
                case=results.get(ref,{}).get(model,{})
                if case.get('status')!='Computed':continue
                records=case['cases']
                if not all(c in records for c in (*cols,'1111')):continue
                on=float(records['1111']['prediction'][name])
                for j,c in enumerate(cols):
                    arr[i,j]=(float(records[c]['prediction'][name])-on)/sigmas[k]
            fig,ax=plt.subplots(figsize=(7.7,8.2))
            # Robust symmetric visual scale; NaNs remain visibly masked.
            finite=np.abs(arr[np.isfinite(arr)])
            lim=float(np.percentile(finite,95)) if finite.size else 1.0
            lim=max(lim,1e-8)
            image=ax.imshow(np.ma.masked_invalid(arr),aspect='auto',cmap='coolwarm',
                            vmin=-lim,vmax=lim)
            ax.set_xticks(range(4),names)
            ax.set_yticks(range(len(ORDER)),[label(x) for x in ORDER])
            ax.set_title(f'{label(ref)} fixed Yukawas: {name}\n'
                         '$(O_{\\mathrm{off}}-O_{1111})/\\sigma_O$')
            fig.colorbar(image,ax=ax,label='Signed change / target 1σ')
            fig.tight_layout()
            image_path=folder/f'{"Dm2" if ref=="D_-2" else "Am4"}_{name}.png'
            fig.savefig(image_path,dpi=200)
            plt.close(fig)
    print('Heatmaps:',folder,flush=True)


def args_parser():
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--run-id',default='optimal_yukawa_thesis_refresh_v1')
    p.add_argument('--source-run',default='thesis16_merged_v1')
    p.add_argument('--canonical-run-id',default=None)
    p.add_argument('--response-root',type=Path,default=ROOT/'output/audits/batch3b_quartic_fixed_20261111_234059')
    p.add_argument('--release-dir',type=Path,default=ROOT/'output/general_observable/lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f')
    p.add_argument('--v4-channel-dir',type=Path,default=ROOT/'output/general_observable/lo_joint_scalar_matching/source_quartic_v4/native_quartic_channels')
    p.add_argument('--quartic-seed-dir',type=Path,default=ROOT/'output/full16_symbolic')
    p.add_argument('--registry',type=Path,default=ROOT/'RGE/general/generated/general_rgbeta_quartic_registry.json')
    p.add_argument('--reference-running',type=Path,default=ROOT/'output/general_observable/all16_v1')
    p.add_argument('--allow-stale-symbolic-source',action='store_true')
    p.add_argument('--run-only',choices=('all','A_-4','D_-2'),default='all')
    p.add_argument('--models',default='all',help='all or comma-separated models for pilots')
    p.add_argument('--full-cases',action='store_true',help='Compute 16 running combinations per model for heatmaps')
    p.add_argument('--check-only',action='store_true')
    p.add_argument('--render-only',action='store_true')
    p.add_argument('--resume',action='store_true',help='Reuse model outputs with matching manifest')
    return p


def main(argv=None):
    args=args_parser().parse_args(argv)
    from scripts import RunT3CoupledYukawaOvernight as baseline
    from scripts import RunT3GeneralAll16FromReference as general
    from scripts import RunT3GeneralSequentialFitAblation as original
    from Numerical.joint_scalar_matching.fermion_boundary_bridge import load_release_manifest,load_source
    from Numerical.fitting.OscillationFit import OscillationFitTarget
    from Numerical.orchestration.RunSmallLambdaRefit import TARGET
    from scripts import RunT3SimultaneousEFT1All16 as simultaneous
    from scripts import RunT3SimultaneousEFT1MixedBC as mixed

    if not args.run_id.replace('-','').replace('_','').isalnum():
        raise ValueError('Unsafe --run-id')
    selected=ORDER if args.models=='all' else tuple(x.strip() for x in args.models.split(',') if x.strip())
    if not selected or len(selected)!=len(set(selected)) or set(selected)-set(ORDER):
        raise ValueError('Invalid --models')
    refs=('D_-2','A_-4') if args.run_only=='all' else (args.run_only,)
    inputs=donor_inputs()
    src_root=STUDIES/'thesis_general_all16'/args.source_run
    settings=read(src_root/'settings.json')
    canonical=args.canonical_run_id or settings.get('canonical_run_id','full16_general_v1')
    out=DEFAULT_OUT/args.run_id
    target=OscillationFitTarget.from_json(TARGET)
    sig={'format':SCHEMA,'references':list(refs),'models':list(selected),
         'full_cases':args.full_cases,'source_run':args.source_run,'canonical_run':canonical,
         'donor_hashes':{m:v['hashes'] for m,v in inputs.items()},
         'target_sha256':sha(TARGET),
         'cached_quartic_context':str(args.response_root.resolve()),
         'extra_yukawas':'target-model-saved',
         'conditional_eft1':True}
    manifest=out/'manifest.json'
    if args.check_only:
        print('Input donors validated against best.json:',list(inputs))
        missing=[]
        for model in selected:
            folder=src_root/general.tag(model)
            for p in (folder/'fitted_config.json',folder/'result.json',args.release_dir/general.tag(model)/'boundary.json',args.v4_channel_dir/(general.tag(model)+'.json')):
                if not p.is_file():missing.append(str(p))
            if model not in baseline.MIXED and not (folder/'general_uv_overrides.json').is_file():
                missing.append(str(folder/'general_uv_overrides.json'))
        for p in (args.registry,args.response_root,TARGET):
            if not p.exists():missing.append(str(p))
        print('Missing prerequisites:',*missing,sep='\n  ')
        return 2 if missing else 0
    if args.render_only:
        existing=read(manifest)
        if existing!=sig:raise ValueError('Manifest differs; do not reformat unrelated dataset')
        data=read(out/'results.json')
        rows=[v['row'] for r in refs for m,v in data[r].items() if 'row' in v]
        recompute_epsilon(rows);render(rows,out)
        if args.full_cases:plot_heatmaps(data,out,target)
        return 0
    if manifest.exists():
        if read(manifest)!=sig:raise ValueError('Existing output manifest differs; choose a new --run-id')
        if not args.resume:raise FileExistsError('Run exists; specify --resume')
        results=read(out/'results.json')
    else:
        if out.exists() and any(out.iterdir()):raise FileExistsError('Nonempty output directory without manifest')
        out.mkdir(parents=True,exist_ok=True)
        write(manifest,sig)
        results={ref:{} for ref in refs}
    registry=read(args.registry)
    release=load_release_manifest(args.release_dir)
    for ref in refs:
        for model in selected:
            existing=results.setdefault(ref,{}).get(model)
            if existing and existing.get('status')=='Computed':
                print('REUSED',ref,model,flush=True);continue
            print('RUN',ref,'->',model,flush=True)
            folder=src_root/general.tag(model)
            if model==ref:
                config=inputs[ref]['config']
                overrides=inputs[ref]['overrides']
            else:
                config=read(folder/'fitted_config.json')
                overrides=read(folder/'general_uv_overrides.json') if model not in baseline.MIXED else None
                values={f'{name}_{i}{j}':float(inputs[ref][name][i,j])
                        for name in ('y1','y2') for i in range(3) for j in range(3)}
                config=general.config_with_yukawas(config,values)
                if overrides is not None:overrides=general.override_yukawas(overrides,values)
            try:
                with tempfile.TemporaryDirectory(prefix='t3_refresh_') as temp:
                    temp=Path(temp)
                    write(temp/'fitted_config.json',config)
                    if overrides is not None:write(temp/'general_uv_overrides.json',overrides)
                    if model in baseline.MIXED:
                        shutil.copy2(folder/'result.json',temp/'result.json')
                    static=load_source(release[model],expected_model=model)
                    channels=original._load_channels(model,static,v4=args.v4_channel_dir,seeds=args.quartic_seed_dir)
                    runner,restore,_=baseline.setup_physics(model,args.response_root)
                    try:
                        with baseline.symbolic_export_compatibility(
                            model,canonical,args.allow_stale_symbolic_source) as compat:
                            if args.full_cases:
                                fn=(mixed._mixed_simultaneous if model in baseline.MIXED
                                    else simultaneous._nonmixed_general)
                            else:fn=runner
                            if model in baseline.MIXED:
                                cases,provenance=fn(model,temp,(static,channels),
                                    args.reference_running,registry,spinor_factor=1.,vev=246.22)
                            else:
                                cases,provenance=fn(model,temp,canonical,(static,channels),registry,
                                    spinor_factor=1.,wolframscript='wolframscript',vev=246.22)
                    finally:restore()
                    expected={f'{u}{e}{w}{s}' for u in '01' for e in '01' for w in '01' for s in '01'} if args.full_cases else {'1111'}
                    if set(cases)!=expected:raise ValueError('Wrong set of running switches')
                    row=extract_row(cases['1111'],ref,model)
                    if model==ref:validate_donor(row,inputs[ref])
                    result={'status':'Computed','row':row,'cases':cases,
                            'provenance':provenance,'symbolic_export_compatibility':compat,
                            'target_model_inputs':str(folder),'donor':inputs[ref]['paths']}
            except (KeyboardInterrupt,SystemExit):raise
            except Exception as exc:
                result={'status':'Blocked','error':f'{type(exc).__name__}: {exc}',
                        'traceback':traceback.format_exc(limit=10)}
                print('BLOCKED',ref,model,result['error'],flush=True)
            results[ref][model]=result
            write(out/'results.json',results)
            print(ref,model,result['status'],flush=True)
    rows=[]
    for ref in refs:
        for model in selected:
            entry=results[ref].get(model,{})
            if entry.get('status')=='Computed':rows.append(entry['row'])
            else:rows.append({'reference':ref,'model':model,'status':'Blocked',
                              'error':entry.get('error','not evaluated')})
    recompute_epsilon(rows)
    render(rows,out)
    if args.full_cases:plot_heatmaps(results,out,target)
    complete=sum(r['status']=='Computed' for r in rows)
    write(out/'summary.json',{'format':SCHEMA,'completed':complete,'requested':len(rows),
                              'conditional_not_certified':True,
                              'donor_baselines_verified':{ref:results[ref].get(ref,{}).get('status')=='Computed' for ref in refs},
                              'errors':[r for r in rows if r['status']!='Computed']})
    return 0 if complete==len(rows) and all(results[ref].get(ref,{}).get('status')=='Computed' for ref in refs) else 2

if __name__=='__main__':
    raise SystemExit(main())
