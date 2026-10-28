"""Batch 17: conditional B/C UV-running Sobol + constrained Yukawa/lambda fit.

Reuses Batch 16 candidate UV/EFT1 -> one-scale C5 -> low-energy evaluator.
No C5 is evaluated for points rejected by the mass screen. Every successful
candidate has the *actual* five-observable chi2 recomputed. The 18 real Yukawa
entries and log10(lambdaT3) are freely optimised inside local Sobol windows.
The matched C5 convention and grouped thresholds remain UNCERTIFIED.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
from hashlib import sha256
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import least_squares

from scripts import RunT3GeneralScalarBCExploratoryFit as b11
from scripts import RunT3GeneralScalarBCPilot as b9
from scripts import RunT3ZeroExtraBoundaryBenchmark as b7
from scripts import DiagnoseT3BCUVScalarHierarchy as b14
from scripts import ProbeT3BCUVNeutrinoNormalization as b16
from scripts import AnalyzeT3BCMatchingConventions as b10

ROOT = b11.ROOT
CASE = 'small_uv_running_5e-9'
FORMAT = 't3-bc-uv-hierarchy-constrained-sobol-refit-v1'
MODELS = b16.MODELS
NAMES = b11.NAMES + ('log10_lambdaT3',)
DEFAULT_FACTORS = {'B_-1': (0.005,0.01,0.015,0.02),
                   'C_-1': (0.005,0.01,0.015)}
PASS_CHI2 = 5.86


def require(test: bool, message: str) -> None:
    if not test: raise ValueError('B/C constrained refit: '+message)


def validate_options(*, factors:tuple[float,...], sources:tuple[str,...],
                     sobol_points:int, seed_count:int, max_nfev:int,
                     local_fraction:float, floor_fraction:float,
                     lambda_low:float, lambda_high:float,
                     yukawa_abs_max:float, ratio_max:float,
                     scalar_min:float, scalar_max:float) -> None:
    require(bool(factors) and len(factors)<=12 and len(set(factors))==len(factors)
            and all(math.isfinite(x) and 0<x<=0.2 for x in factors),
            'invalid seed factors')
    require(bool(sources) and set(sources)<=set(b16.SOURCES)
            and len(sources)==len(set(sources)), 'invalid seed sources')
    require(isinstance(sobol_points,int) and sobol_points>=2
            and (sobol_points & (sobol_points-1))==0,
            'Sobol points must be a power of two >= 2')
    require(1<=seed_count<=8 and 1<=max_nfev<=2000,'invalid local fit budgets')
    require(0<local_fraction<=3 and 0<floor_fraction<=1,
            'invalid local perturbation fractions')
    require(1e-14<=lambda_low<lambda_high<=1.,'invalid topology-coupling interval')
    require(0<yukawa_abs_max<=3.,'invalid Yukawa absolute cap')
    require(0<ratio_max<1 and 0<scalar_min<scalar_max,
            'invalid scalar/fermion mass screen')


def destination(root:Path, run_id:str) -> Path:
    rid=b9.parse_run_id(run_id)
    base=(Path(root)/'output/studies/general_scalar_extended'/CASE).resolve()
    dst=(base/rid).resolve()
    require(dst.is_relative_to(base) and not dst.exists(),
            'invalid/colliding output ID: '+str(dst))
    return dst


def candidate_record(historical:dict, params:np.ndarray,
                     *,lambda_low:float,lambda_high:float) -> tuple[dict,float]:
    x=np.asarray(params,dtype=float)
    require(x.shape==(19,) and np.isfinite(x).all(),'expected 19 finite parameters')
    lam=float(10.**x[-1])
    require(lambda_low*(1-1e-12)<=lam<=lambda_high*(1+1e-12) and lam<=1.,
            'lambdaT3 outside parameter bounds')
    lam=float(np.clip(lam,lambda_low,lambda_high))
    out=deepcopy(historical)
    out['lambdaT3']=complex(lam)
    out['cfg']['base_state']['ordinary']['lambdaT3']=b7.complex_json(complex(lam))
    return b11.trial_record(out,x[:18]),lam


def whiten(target:Any, prediction:dict[str,float]) -> np.ndarray:
    names=tuple(target.observable_names)
    vals=np.asarray([float(prediction[n]) for n in names],dtype=float)
    require(vals.shape==(5,) and np.isfinite(vals).all(),'nonfinite oscillation observables')
    central=np.asarray(target.central_values,dtype=float)
    cov=np.asarray(target.covariance,dtype=float)
    residual=np.linalg.solve(np.linalg.cholesky(cov),vals-central)
    require(residual.shape==(5,) and np.isfinite(residual).all(),'nonfinite residual')
    return residual


def screen_penalty(row:dict,*,ratio_max:float,scalar_min:float,
                   scalar_max:float)->np.ndarray:
    """Nonzero continuous-scale failure penalty; never a fake chi2 candidate.

    This guides least-squares back toward a *screened* point. We never claim
    its squared norm to be a NuFIT chi2, or keep it as a best fit.
    """
    rho=row.get('rhoF')
    hi=row.get('scalar_heavy_muS_tev')
    lo=row.get('scalar_light_muS_tev')
    def number(x, default):
        try:
            x=float(x)
            return x if math.isfinite(x) else default
        except (TypeError,ValueError): return default
    r=number(rho,2*ratio_max)
    upper=number(hi,2*scalar_max)
    lower=number(lo,scalar_min/2)
    return np.array([1e5*(1+max(0,r/ratio_max-1)),
                     1e5*(1+max(0,upper/scalar_max-1)),
                     1e5*(1+max(0,scalar_min/max(lower,1e-100)-1)),
                     1e5,1e5],dtype=float)


class ScreenedEvaluator:
    """All calls use live Batch16 RG, matching and low-energy calculations."""
    def __init__(self,model:str,source:dict,registry:dict,bundle:dict,
                 compiled:Any,target:Any,*,q:float,mu:float,
                 ratio_max:float,scalar_min:float,scalar_max:float,
                 lambda_low:float,lambda_high:float):
        self.model=model;self.source=source;self.registry=registry
        self.bundle=bundle;self.compiled=compiled;self.target=target
        self.q=q;self.mu=mu;self.ratio_max=ratio_max
        self.scalar_min=scalar_min;self.scalar_max=scalar_max
        self.lambda_low=lambda_low;self.lambda_high=lambda_high
        self.calls=0;self.screened=0;self.failed=0;self.match_calls=0
        self.best=None;self.cache={}

    def evaluate(self,x:np.ndarray,*,label:str)->dict:
        arr=np.asarray(x,dtype=float)
        require(arr.shape==(19,) and np.isfinite(arr).all(),'malformed trial vector')
        key=arr.tobytes()
        if key in self.cache:return self.cache[key]
        self.calls+=1
        row={'label':label,'parameters':dict(zip(NAMES,map(float,arr))),
             'chi2':None,'screen_pass':False,'prediction':None}
        try:
            record,lam=candidate_record(self.source,arr,
                              lambda_low=self.lambda_low,lambda_high=self.lambda_high)
            row['lambdaT3_uv']=lam
            result=b16.run_one(self.model,record,self.registry,self.bundle,
                        self.compiled,self.target,quartic_fraction=self.q,
                        mu2_fraction=self.mu,ratio_max=self.ratio_max,
                        scalar_min_tev=self.scalar_min,scalar_max_tev=self.scalar_max)
            row.update(result)
            if result.get('screen_pass') and result.get('chi2') is not None:
                residual=whiten(self.target,result['prediction'])
                chi2=float(np.dot(residual,residual))
                require(np.isclose(chi2,float(result['chi2']),rtol=1e-7,atol=1e-7),
                        'target chi2 does not agree with NuFIT evaluator')
                row['chi2']=chi2
                row['whitened_residual']=residual.tolist()
                self.screened+=1;self.match_calls+=1
                if self.best is None or chi2 < self.best['chi2']:
                    self.best=deepcopy(row)
            else:
                row['whitened_residual']=screen_penalty(result,
                      ratio_max=self.ratio_max,scalar_min=self.scalar_min,
                      scalar_max=self.scalar_max).tolist()
        except Exception as exc:
            self.failed+=1
            row.update(status='NumericalTrialFailed',reason=str(exc)[:500],
                       whitened_residual=[1e6]*5)
        self.cache[key]=row
        return row


def seed_anchors(model:str,historical:dict,saved:dict,*,
                 factors:tuple[float,...],sources:tuple[str,...],
                 lambda_low:float,lambda_high:float,
                 yukawa_abs_max:float,local_fraction:float,
                 floor_fraction:float) -> list[dict]:
    original=float(complex(historical['lambdaT3']).real)
    require(original>0,'nonpositive historical lambdaT3')
    anchors=[]
    for source_name in sources:
        for f in factors:
            seed=b16.seeded_yukawa(saved,historical,source_name=source_name,factor=f)
            require(np.max(np.abs(seed))<yukawa_abs_max,
                    'seed exceeds Yukawa cap; increase --yukawa-abs-max')
            centre_lam=float(np.clip(original/f**2,lambda_low,lambda_high))
            lo_half=np.maximum(np.abs(seed)*local_fraction,
                               np.max(np.abs(seed))*floor_fraction)
            low=np.maximum(-yukawa_abs_max,seed-lo_half)
            high=np.minimum(yukawa_abs_max,seed+lo_half)
            require(np.all(low<high),'collapsed Yukawa Sobol bound')
            llow=max(math.log10(lambda_low),math.log10(centre_lam)-2.)
            lhigh=min(math.log10(lambda_high),math.log10(centre_lam)+2.)
            require(llow<lhigh,'collapsed lambda Sobol bound')
            bounds=np.vstack((np.r_[low,llow],np.r_[high,lhigh]))
            x=np.r_[seed,math.log10(centre_lam)]
            anchors.append({'model':model,'source':source_name,'factor':f,
                            'seed':x,'bounds':bounds})
    return anchors


def sobol_points(anchor:dict,n:int,*,random_seed:int)->list[np.ndarray]:
    from Numerical.fitting.SobolBenchmarkSearch import (BenchmarkParameter,
          generate_sobol_parameter_points)
    lo,hi=anchor['bounds']
    defs=tuple(BenchmarkParameter(name,'ordinary.'+name,float(a),float(b))
            for name,a,b in zip(NAMES,lo,hi))
    generated=generate_sobol_parameter_points(defs,n,seed=random_seed)
    return [np.asarray([float(row[k]) for k in NAMES],dtype=float)
            for row in generated]


def local_refinement(evaluator:ScreenedEvaluator,entry:dict,anchor:dict,
                     *,max_nfev:int)->dict:
    """All 19 variables in bounded unit coordinates; no frozen Yukawa subset."""
    lo,hi=anchor['bounds'];width=hi-lo
    start=np.array([entry['parameters'][k] for k in NAMES],dtype=float)
    z0=np.clip((start-lo)/width,1e-8,1-1e-8)
    records=[]
    def residual(z):
        x=lo+np.asarray(z,dtype=float)*width
        row=evaluator.evaluate(x,label='least_squares')
        records.append((row['status'],row['chi2']))
        return np.asarray(row['whitened_residual'],dtype=float)
    sol=least_squares(residual,z0,bounds=(np.zeros(19),np.ones(19)),
             method='trf',jac='2-point',max_nfev=max_nfev,
             xtol=1e-6,ftol=1e-6,gtol=1e-6)
    return {'anchor_source':anchor['source'],'anchor_factor':anchor['factor'],
            'start_chi2':float(entry['chi2']),
            'max_nfev':max_nfev,'optimiser_nfev':int(sol.nfev),
            'optimiser_success':bool(sol.success),
            'optimiser_message':str(sol.message),
            'function_calls_recorded':len(records),
            'screened_function_calls':sum(v is not None for _,v in records),
            'best_screened_chi2_so_far':(evaluator.best['chi2'] if evaluator.best else None)}


def run_model(*,model:str,run_id:str,root:Path=ROOT,
              factors:tuple[float,...]|None=None,
              sources:tuple[str,...]=b16.SOURCES,
              sobol_n:int=8,seed_count:int=2,max_nfev:int=40,
              random_seed:int=20261028,
              local_fraction:float=0.8,floor_fraction:float=0.03,
              lambda_low:float=1e-9,lambda_high:float=0.1,
              yukawa_abs_max:float=0.15,
              ratio_max:float=0.5,scalar_min:float=0.5,
              scalar_max:float=5.0,
              prior_run_id:str|None=None) -> dict:
    require(model in MODELS,'only B_-1/C_-1 supported')
    fs=DEFAULT_FACTORS[model] if factors is None else tuple(factors)
    validate_options(factors=fs,sources=tuple(sources),sobol_points=sobol_n,
       seed_count=seed_count,max_nfev=max_nfev,local_fraction=local_fraction,
       floor_fraction=floor_fraction,lambda_low=lambda_low,lambda_high=lambda_high,
       yukawa_abs_max=yukawa_abs_max,ratio_max=ratio_max,
       scalar_min=scalar_min,scalar_max=scalar_max)
    dest=destination(root,run_id)
    from Numerical.experimental.GeneralT3ScalarRGState import load_registry
    from Numerical.orchestration import RunSmallLambdaRefit as base
    from Numerical.fitting.OscillationFit import OscillationFitTarget
    from scripts.RunGeneralT3NumericalC5 import DEFAULT_SOURCE
    registry=load_registry(root/'RGE/general/generated/general_rgbeta_quartic_registry.json')
    saved,historical,settings=b14.fit_source(root,model,
        prior_run_id or b14.DEFAULT_RUNS[model],registry)
    require(settings['quartic_fraction']==0.1 and settings['mu2_fraction']==0,
            'current constrained runs require the tested q=0.1, mu2=0 prescription')
    bundle=b11.beta_bundle(model,historical,registry)
    coeff=Path(DEFAULT_SOURCE)/({'B_-1':'B_m1_final','C_-1':'C_m1_final'}[model]
                   +'_one_loop_probe_weinberg_ordered_coefficient.inputform.txt')
    compiled=b11.compile_source(coeff)
    b10.check_convention_evidence(b7.read_json(b10.EVIDENCE23),
          b7.read_json(b10.EVIDENCE26),models=(model,))
    target=OscillationFitTarget.from_json(base.TARGET)
    evaluator=ScreenedEvaluator(model,historical,registry,bundle,compiled,target,
           q=settings['quartic_fraction'],mu=settings['mu2_fraction'],
           ratio_max=ratio_max,scalar_min=scalar_min,scalar_max=scalar_max,
           lambda_low=lambda_low,lambda_high=lambda_high)
    anchors=seed_anchors(model,historical,saved,factors=fs,sources=tuple(sources),
           lambda_low=lambda_low,lambda_high=lambda_high,
           yukawa_abs_max=yukawa_abs_max,local_fraction=local_fraction,
           floor_fraction=floor_fraction)
    candidates=[];counts={};errors=[]
    for i,anchor in enumerate(anchors):
        points=[anchor['seed'],*sobol_points(anchor,sobol_n,
                  random_seed=random_seed+1009*i)]
        for j,x in enumerate(points):
            r=evaluator.evaluate(x,label=f'anchor_{i+1}_point_{j}')
            counts[r['status']]=counts.get(r['status'],0)+1
            if r['chi2'] is not None:
                candidates.append((r,i))
            elif len(errors)<15:
                errors.append({'anchor':i+1,'point':j,'reason':r.get('reason','screen failed'),
                               'status':r['status']})
            if j==0 or (r['chi2'] is not None and
                    r['chi2']==evaluator.best['chi2']):
                print(f"{model}: {anchor['source']} f={anchor['factor']:g} "
                      f"point {j}/{sobol_n}: {r['status']} "
                      f"chi2={r['chi2']} rhoF={r.get('rhoF')}",flush=True)
    candidates.sort(key=lambda row:row[0]['chi2'])
    history=[]
    if candidates and candidates[0][0]['chi2']>PASS_CHI2:
        for row,i in candidates[:seed_count]:
            history.append(local_refinement(evaluator,row,anchors[i],max_nfev=max_nfev))
            print(f"{model}: local {len(history)} best screened "
                  f"chi2={evaluator.best['chi2'] if evaluator.best else None}",flush=True)
            if evaluator.best is not None and evaluator.best['chi2']<=PASS_CHI2:break
    best=evaluator.best
    status=('CandidateConstrainedAccepted' if best and best['chi2']<=PASS_CHI2
            else 'CandidateConstrainedPoorFit' if best else 'NoMassSafeNeutrinoFit')
    paths=[settings['saved_fit_path'],settings['saved_config_path'],
           historical['cfg_path'],historical['fit_path'],coeff,base.TARGET,
           root/'RGE/general/generated/general_rgbeta_quartic_registry.json',
           *bundle['files'].values(),b10.EVIDENCE23,b10.EVIDENCE26]
    provenance={str(p.relative_to(root)) if p.is_relative_to(root) else str(p):
               sha256(p.read_bytes()).hexdigest() for p in paths}
    report={'format':FORMAT,'case':CASE,'model':model,'run_id':run_id,
       'status':status,'accepted_threshold_chi2':PASS_CHI2,
       'model_count':1,'quartic_fraction':settings['quartic_fraction'],
       'mu2_fraction':settings['mu2_fraction'],
       'mass_screen':{'rhoF_max':ratio_max,'scalar_min_tev':scalar_min,
                      'scalar_max_tev':scalar_max},
       'search':{'sobol_points_per_anchor':sobol_n,'anchors':len(anchors),
           'seed_sources':list(sources),'seed_factors':list(fs),
           'local_seed_count':seed_count,'local_max_nfev':max_nfev,
           'local_fraction':local_fraction,'floor_fraction':floor_fraction,
           'lambda_low':lambda_low,'lambda_high':lambda_high,
           'yukawa_abs_max':yukawa_abs_max,'random_seed':random_seed},
       'run_counts':{'evaluator_calls':evaluator.calls,
             'mass_safe_c5_evaluations':evaluator.screened,
             'failed_numerical_evaluations':evaluator.failed,
             'status_counts':counts},
       'best_mass_safe':best,'optimiser':history,'trial_errors':errors,
       'source_sha256':provenance,
       'physical_matchete_rgbeta_map_certified':False,
       'general_scalar_sequential_matching_performed':False,
       'neutrino_chi2_is_candidate_only':True,
       'historical_outputs_overwritten':False}
    dest.mkdir(parents=True,exist_ok=False)
    b7.save_json(dest/'fit.json',report)
    if best is not None:
        params=np.asarray([best['parameters'][k] for k in NAMES],dtype=float)
        fitted,lam=candidate_record(historical,params,
                    lambda_low=lambda_low,lambda_high=lambda_high)
        cfg=fitted['cfg']
        cfg['benchmark_search']={'enabled':False}
        cfg['sensitivity']={'enabled':False}
        cfg['constrained_general_scalar_uv_fit']={
            'status':status,'case':CASE,'model':model,'run_id':run_id,
            'candidate_chi2':best['chi2'],'rhoF':best['rhoF'],
            'lambdaT3_uv':lam,'physical_matching_certified':False,
            'sequential_matching_performed':False}
        b7.save_json(dest/'fitted_config.json',cfg)
    (dest/'summary.md').write_text(summary_markdown(report),encoding='utf-8')
    print('Saved:',dest/'summary.md',flush=True)
    return report


def summary_markdown(report:dict)->str:
    best=report['best_mass_safe']
    lines=['# Batch 17 — constrained UV-running B/C Sobol + optimisation','',
      f"Model: `{report['model']}`; run: `{report['run_id']}`; status: `{report['status']}`.",
      f"Candidate accepted if chi2 <= {PASS_CHI2}, rhoF <= {report['mass_screen']['rhoF_max']}, "
      f"and mS(1 TeV) in [{report['mass_screen']['scalar_min_tev']}, "
      f"{report['mass_screen']['scalar_max_tev']}] TeV.",'',
      f"Sobol: {report['search']['anchors']} anchors × "
      f"{report['search']['sobol_points_per_anchor']} points each.",
      f"Mass-safe matched candidate evaluations: {report['run_counts']['mass_safe_c5_evaluations']}; "
      f"failed numerical evaluations: {report['run_counts']['failed_numerical_evaluations']}.",'']
    if best is None:
        lines+=['No mass-safe point reached the conditional neutrino evaluator.','']
    else:
        lines += [f"Best screened chi2: **{best['chi2']:.8g}**",
          f"rhoF: {best.get('rhoF'):.6g}",
          f"scalar masses at 1 TeV: {best.get('scalar_light_muS_tev'):.6g}–"
          f"{best.get('scalar_heavy_muS_tev'):.6g} TeV",
          f"UV lambdaT3: {best['lambdaT3_uv']:.8g}",'',
          '| Observable | Value |','|---|---:|',
          *[f'| {k} | {float(v):.10g} |' for k,v in best['prediction'].items()],'']
    lines += ['## Interpretation','',
       '- All five NuFIT residuals and mass screens are recomputed at each successful trial.',
       '- An invalid mass screen is never included as a valid chi2 point.',
       '- The sampled UV mass ordering and 1 TeV scalar window are numerical conditions, not full threshold certification.',
       '- The Matchete–RGBeta coupling normalization, one-scale C5 and sequential matching remain conditional.',
       '- No historical or earlier exploratory outputs are overwritten.','']
    return '\n'.join(lines)+'\n'


def main() -> None:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--model',required=True,choices=MODELS)
    p.add_argument('--run-id',required=True)
    p.add_argument('--seed-factors',help='CSV seed factors <= 0.2; default model-specific Batch15 screen passes')
    p.add_argument('--seed-source',choices=('both',)+b16.SOURCES,default='both')
    p.add_argument('--sobol-points',type=int,default=8)
    p.add_argument('--seed-count',type=int,default=2)
    p.add_argument('--max-nfev',type=int,default=40)
    p.add_argument('--seed',type=int,default=20261028)
    p.add_argument('--local-fraction',type=float,default=0.8)
    p.add_argument('--floor-fraction',type=float,default=0.03)
    p.add_argument('--lambda-low',type=float,default=1e-9)
    p.add_argument('--lambda-high',type=float,default=0.1)
    p.add_argument('--yukawa-abs-max',type=float,default=0.15)
    p.add_argument('--rho-max',type=float,default=0.5)
    p.add_argument('--scalar-min-tev',type=float,default=0.5)
    p.add_argument('--scalar-max-tev',type=float,default=5.)
    p.add_argument('--previous-run-id')
    a=p.parse_args()
    factors=(b16.parse_positive_csv(a.seed_factors,'Yukawa factors',maximum=0.2)
             if a.seed_factors else None)
    sources=b16.SOURCES if a.seed_source=='both' else (a.seed_source,)
    run_model(model=a.model,run_id=a.run_id,factors=factors,sources=sources,
       sobol_n=a.sobol_points,seed_count=a.seed_count,max_nfev=a.max_nfev,
       random_seed=a.seed,local_fraction=a.local_fraction,
       floor_fraction=a.floor_fraction,lambda_low=a.lambda_low,
       lambda_high=a.lambda_high,yukawa_abs_max=a.yukawa_abs_max,
       ratio_max=a.rho_max,scalar_min=a.scalar_min_tev,scalar_max=a.scalar_max_tev,
       prior_run_id=a.previous_run_id)

if __name__=='__main__': main()
