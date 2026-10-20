"""Local-only T3 small-lambda study. Requires existing model-specific RGE/C5 JSONs.

Stages: prepare -> fit -> cluster. No symbolic generation and no GitHub writes.
"""
from __future__ import annotations
import argparse
import copy
import csv
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
from scipy.optimize import least_squares
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist

ROOT = Path(__file__).resolve().parents[2]
REF = ROOT / 'Numerical/benchmarks/fitted_reference'
STUDY = ROOT / 'output/studies/small_lambda_1e-7'
TARGET = ROOT / 'configs/nufit_6_1_ic24_no.json'


def load(p):
    return json.loads(Path(p).read_text(encoding='utf-8-sig'))


def save(p, d):
    p = Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(d, indent=2) + '\n', encoding='utf-8')


def matr(d):
    return np.asarray(d['real'], dtype=float) + 1j * np.asarray(d.get('imag', np.zeros((3, 3))), dtype=float)


def setmat(d, x):
    d['real'] = np.asarray(x.real, dtype=float).tolist()
    d['imag'] = np.asarray(x.imag, dtype=float).tolist()


def cls(row):
    d = row['representation']
    dims = (d['d_s1'], d['d_s2'], d['d_f'])
    return {(1,3,2):'A', (2,2,1):'B', (2,2,3):'C', (3,1,2):'D', (3,3,2):'E'}[dims]


def original_rows():
    m = load(REF / 'benchmark_manifest.json')
    rows = m['models']
    if len(rows) != 16 or len({r['model'] for r in rows}) != 16:
        raise ValueError('Expected 16 unique fitted reference models')
    return rows


def prepare():
    records = []
    for row in original_rows():
        name = row['model']
        src = REF / (name + '.json')
        payload = load(src)
        kind = payload['kind']
        if kind != 't3_pipeline_numerical_results_v1' or payload['representation'] != row['representation']:
            raise ValueError(f'Unexpected config schema/representation: {src}')
        if row['representation'].get('shared_scalar'):
            raise ValueError(f'Unsupported shared-scalar model {name}')
        old_l = float(payload['base_state']['ordinary']['lambdaT3']['real'])
        new_l = 5e-8 if cls(row) == 'E' else 1e-7
        if old_l <= 0:
            raise ValueError(f'Expected positive fitted lambda in {name}')
        factor = np.sqrt(old_l / new_l)
        modified = copy.deepcopy(payload)
        o = modified['base_state']['ordinary']
        for key in ('y1','y2'):
            setmat(o[key], matr(o[key]) * factor)
        o['lambdaT3'] = {'real':new_l, 'imag':0.0}
        modified['benchmark_search'] = {'enabled':False}
        modified['sensitivity'] = {'enabled':False}
        # Historical fixed-grid display metadata is not the new fitted point.
        modified.pop('comparison_point', None)
        modified['small_lambda_study'] = {
            'source':str(src.relative_to(ROOT)), 'source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),
            'original_lambda':old_l, 'lambdaT3':new_l,
            'yukawa_initial_scale_factor':float(factor),
            'status':'Analytical rescaling seed only; not yet a fit',
        }
        out = STUDY / 'prepared' / (name+'.json')
        save(out, modified)
        records.append({'model':name,'class':cls(row),'lambda':new_l,'yukawa_factor':factor,'config':str(out.relative_to(ROOT))})
    save(STUDY / 'prepared_manifest.json',records)
    print(f'Prepared {len(records)} UV seeds in {STUDY / "prepared"}. These are NOT new fitted solutions.')


def locate_artifacts(model, artifact_map):
    if artifact_map:
        entry = artifact_map.get(model)
        if entry is None:
            raise FileNotFoundError(f'Missing artifact-map entry for {model}')
        result = []
        for key in ('uv_rgbeta','eft1_rgbeta','final_weinberg'):
            p = Path(entry[key]); p = p if p.is_absolute() else ROOT / p
            if not p.is_file(): raise FileNotFoundError(p)
            result.append(p)
        return result
    # Exact model-directory match: refuse ambiguity rather than guessing artifact identity.
    base = ROOT / 'output/full/comparison/smallY_smallL' / model
    names = ('uv_rgbeta_rge.json','eft1_rgbeta_rge.json','final_weinberg_coefficient.json')
    result=[]
    for name in names:
        hits = list(base.rglob(name)) if base.is_dir() else []
        if len(hits) != 1:
            raise FileNotFoundError(f'{model}: expected exactly one {name} under {base}; found {len(hits)}. Supply --artifact-map.')
        result.append(hits[0])
    return result


def get_evaluator(cfg, artifacts):
    # Reuse exactly the repository's numerical evaluator, RGEs, and C5 matching.
    from Numerical.fitting.SobolBenchmarkSearch import BenchmarkParameter, _build_evaluator
    params = [BenchmarkParameter(f'{y}_{i}{j}', f'ordinary.{y}.real.{i}.{j}', -1., 1.)
              for y in ('y1','y2') for i in range(3) for j in range(3)]
    evalfn = _build_evaluator(payload=cfg, parameters=params, uv_rgbeta_path=artifacts[0],
                             intermediate_rgbeta_path=artifacts[1], final_weinberg_path=artifacts[2],
                             target_path=TARGET)
    return params,evalfn


def evaluate_point(evaluator, params, x):
    return evaluator({p.name:float(v) for p,v in zip(params,x)})


def fit(args):
    rows=original_rows()
    mapping=load(args.artifact_map) if args.artifact_map else None
    report=[]
    for row in rows:
        model=row['model']
        if args.model and model != args.model: continue
        cfg=load(STUDY/'prepared'/(model+'.json'))
        if abs(cfg['base_state']['ordinary']['lambdaT3']['real']-(5e-8 if cls(row)=='E' else 1e-7))>1e-20:
            raise ValueError('Incorrect new lambda in '+model)
        artifacts=locate_artifacts(model,mapping)
        params,evaluator=get_evaluator(cfg,artifacts)
        o=cfg['base_state']['ordinary']
        x0=np.concatenate([matr(o[y]).real.ravel() for y in ('y1','y2')])
        if np.any(np.abs(x0)>=1): raise ValueError(f'{model}: rescaled seed exceeds fit bound +/-1')
        # Optimise 18 real Yukawa entries. Imaginary parts are kept at their stored values.
        def residual(x):
            try:
                r=np.asarray(evaluate_point(evaluator,params,x).whitened_residual,dtype=float)
                if r.shape!=(5,) or not np.all(np.isfinite(r)): raise ValueError('Invalid residual')
                return r
            except Exception:
                return np.full(5,1e6)
        seed_resid=residual(x0)
        sol=least_squares(residual,x0,bounds=(-np.ones(18),np.ones(18)),max_nfev=args.max_nfev,
                          diff_step=1e-4,xtol=1e-9,ftol=1e-9,gtol=1e-9)
        valid=False
        try:
            result=evaluate_point(evaluator,params,sol.x)
            valid=np.isfinite(result.chi2)
        except Exception as exc:
            print(f'FIT ERROR {model}: {exc}',file=sys.stderr);result=None
        if not valid:
            report.append({'model':model,'status':'Failed'});continue
        fitted=copy.deepcopy(cfg)
        for index,y in enumerate(('y1','y2')):
            setmat(fitted['base_state']['ordinary'][y],sol.x[index*9:(index+1)*9].reshape(3,3)+1j*matr(o[y]).imag)
        fitted['small_lambda_study']['status']='Numerically optimised' if sol.success else 'Best effort: optimizer did not converge'
        fitted['small_lambda_study']['chi2']=float(result.chi2)
        fitted['small_lambda_study']['max_nfev']=args.max_nfev
        save(STUDY/'fitted'/(model+'.json'),fitted)
        summary={'model':model,'class':cls(row),'lambda':float(o['lambdaT3']['real']),
                 'chi2_start':float(np.dot(seed_resid,seed_resid)), 'chi2_fit':float(result.chi2),
                 'converged':bool(sol.success),'nfev':int(sol.nfev),
                 'prediction':dict(zip(result.observable_names,map(float,result.prediction))),
                 'artifacts':list(map(str,artifacts))}
        save(STUDY/'fits'/(model+'.json'),summary);report.append(summary)
        print(model,'chi2:',summary['chi2_fit'],'success:',sol.success,flush=True)
    save(STUDY/'fit_summary.json',report)
    print('Saved',len(report),'model fit records')


def cluster(args):
    rows=original_rows(); fitted=[]
    for row in rows:
        p=STUDY/'fitted'/(row['model']+'.json')
        if p.is_file(): fitted.append((row,load(p)))
    if len(fitted)<2: raise RuntimeError('Need at least two successful fitted model configs')
    # Exploratory *raw-entry* distances; no flavour-basis or rephasing alignment.
    vectors=np.array([np.concatenate([matr(c['base_state']['ordinary'][y]).real.ravel()
                                      for y in ('y1','y2')]) for _,c in fitted]);
    distances=pdist(vectors,metric='euclidean')
    labels=fcluster(linkage(distances,method='average'),t=args.distance,criterion='distance')
    memberships={}
    for (row,cfg),lab in zip(fitted,labels):memberships.setdefault(str(int(lab)),[]).append(row['model'])
    save(STUDY/'clusters.json',{'method':'average-linkage Euclidean raw real Yukawa entries',
        'distance_threshold':args.distance,'warning':'No rephasing, basis alignment, or physical texture equivalence established',
        'memberships':memberships})
    mapping=load(args.artifact_map) if args.artifact_map else None
    comparisons=[]
    # Mean matrices define an exploratory common-Yukawa candidate; never relabel it a refit.
    for label,models in memberships.items():
        indices=[k for k,(row,_) in enumerate(fitted) if row['model'] in models]
        mean=np.mean(vectors[indices],axis=0)
        for k in indices:
            row,cfg=fitted[k]
            common=copy.deepcopy(cfg)
            for i,y in enumerate(('y1','y2')):
                setmat(common['base_state']['ordinary'][y],mean[i*9:(i+1)*9].reshape(3,3)+1j*matr(cfg['base_state']['ordinary'][y]).imag)
            common['small_lambda_study']['status']='Unoptimised cluster mean Yukawa';common['small_lambda_study']['cluster']=label
            save(STUDY/'common_yukawa'/ (row['model']+'.json'),common)
            record={'model':row['model'],'cluster':label,'members':len(models),'fitted_chi2':cfg['small_lambda_study']['chi2']}
            try:
                artifacts=locate_artifacts(row['model'],mapping)
                params,ev=get_evaluator(common,artifacts)
                value=evaluate_point(ev,params,mean)
                record['common_chi2']=float(value.chi2)
                record['common_prediction']=dict(zip(value.observable_names,map(float,value.prediction)))
                record['status']='Evaluated'
            except Exception as exc:
                record['status']='Unavailable';record['error']=str(exc)
            comparisons.append(record)
    save(STUDY/'common_comparison.json',comparisons)
    with (STUDY/'common_comparison.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=['model','cluster','members','fitted_chi2','common_chi2','status','error']);w.writeheader();w.writerows([{k:r.get(k,'') for k in w.fieldnames} for r in comparisons])
    try:
        import matplotlib;matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig,ax=plt.subplots(figsize=(8,5));ax.scatter(range(len(fitted)),labels)
        ax.set_xticks(range(len(fitted)),[r['model'] for r,_ in fitted],rotation=90,fontsize=6)
        ax.set_ylabel('Exploratory cluster ID');fig.tight_layout()
        fig.savefig(STUDY/'cluster_assignments.png',dpi=160);plt.close(fig)
        fig,ax=plt.subplots(figsize=(9,4))
        ax.plot([r.get('fitted_chi2',np.nan) for r in comparisons],label='Individually fitted chi2',marker='o')
        ax.plot([r.get('common_chi2',np.nan) for r in comparisons],label='Cluster mean chi2',marker='x')
        ax.set_yscale('symlog',linthresh=1e-3);ax.set_ylabel('Approx. chi2');ax.legend();fig.tight_layout()
        fig.savefig(STUDY/'chi2_comparison.png',dpi=160);plt.close(fig)
        plots_dir=STUDY/'yukawa_matrix_plots'
        plots_dir.mkdir(parents=True,exist_ok=True)
        for row,cfg in fitted:
            model=row['model']
            common=load(STUDY/'common_yukawa'/(model+'.json'))
            fig,axes=plt.subplots(2,2,figsize=(8,7),layout='constrained')
            for j,name in enumerate(('y1','y2')):
                for i,(tag,payload) in enumerate((('Individual',cfg),('Cluster mean',common))):
                    arr=np.abs(matr(payload['base_state']['ordinary'][name]))
                    im=axes[i,j].imshow(arr,vmin=0,vmax=max(0.1,float(np.max(arr))),cmap='viridis')
                    axes[i,j].set_title(f'{tag} {name.upper()}')
                    axes[i,j].set_xticks((0,1,2));axes[i,j].set_yticks((0,1,2))
                    for a in range(3):
                        for b in range(3):axes[i,j].text(b,a,f'{arr[a,b]:.2g}',ha='center',va='center',fontsize=7,color='white' if arr[a,b]>0.05 else 'black')
            fig.suptitle(model)
            fig.savefig(plots_dir/(model+'.png'),dpi=150);plt.close(fig)
    except ImportError:print('matplotlib unavailable: skipped figures')
    print('Clusters',memberships,'\nResults:',STUDY/'common_comparison.csv')


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('stage',choices=('prepare','fit','cluster'))
    ap.add_argument('--artifact-map',type=Path,help='Optional JSON: {model: {uv_rgbeta,eft1_rgbeta,final_weinberg}}')
    ap.add_argument('--max-nfev',type=int,default=300)
    ap.add_argument('--model',help='Only fit this exact manifest model ID')
    ap.add_argument('--distance',type=float,default=0.35,help='Exploratory Euclidean clustering threshold')
    a=ap.parse_args()
    {'prepare':prepare,'fit':lambda:fit(a),'cluster':lambda:cluster(a)}[a.stage]()
if __name__=='__main__':main()
