"""Exploratory B(-1)/C(-1) general-scalar Yukawa refits; NEW outputs only.

Candidate native RGBeta UV/EFT1 including seven additional quartics and the
holomorphic scalar quadratic. Evaluate saved Matchete ONE-SCALE ordered C5 in
the actual mixed-scalar eigenbasis. The two direct quartics use the Stage23/26
CONDITIONAL SIGNED map lambda_M = -lambda_RGBeta (NOT PHYSICALLY CERTIFIED).

For each trial point, rerun the appropriate UV/EFT1 RGEs, diagonalise MF in the
existing heavy-fermion mass-basis convention at muF, evolve SM+Weinberg from
muS to 100 GeV, and score NuFIT five-observable chi^2. Fits are not physical
sequential general-potential predictions. No previous outputs are modified.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
from hashlib import sha256
import json
import math
from pathlib import Path
import re
from typing import Any

import numpy as np

from scripts import RunT3ZeroExtraBoundaryBenchmark as b7
from scripts import RunT3GeneralScalarBCPilot as b9
from scripts import AnalyzeT3BCMatchingConventions as b10

ROOT = Path(__file__).resolve().parents[1]
OUT_ROOT = ROOT / 'output/studies/general_scalar_extended'
FORMAT = 't3-bc-general-scalar-exploratory-yukawa-refit-v1'
NAMES = tuple(f'{field}_{i}{j}' for field in ('y1','y2') for i in range(3) for j in range(3))
FIELDS = ('delta_m21_sq_ev2','delta_m3l_sq_ev2','sin2_theta12','sin2_theta13','sin2_theta23')


def require(condition:bool, msg:str)->None:
    if not condition:
        raise ValueError('Exploratory B/C refit: '+msg)


def real_finite(x:Any, name:str)->float:
    z=float(x)
    require(math.isfinite(z),name+' must be finite')
    return z


def fit_destination(case:str, run_id:str, *,root:Path=ROOT)->Path:
    require(case in b9.CASES,'unsupported historical case')
    rid=b9.parse_run_id(run_id)
    base=(Path(root)/'output/studies/general_scalar_extended').resolve()
    dest=(base/case/rid).resolve()
    require(dest.is_relative_to(base),'output path escapes general-scalar tree')
    require(not dest.exists(),'run ID already exists: '+str(dest))
    return dest


def seed_vector(source:dict)->np.ndarray:
    ordinary=source['cfg']['base_state']['ordinary']
    from Numerical.orchestration import RunSmallLambdaRefit as base
    matrices=[base.matr(ordinary[k]) for k in ('y1','y2')]
    require(all(z.shape==(3,3) and np.max(abs(z.imag))<1e-13 for z in matrices),
            'current pilot optimises the historical 18 real Yukawa components only')
    return np.concatenate([z.real.ravel() for z in matrices])


def trial_record(source:dict, values:np.ndarray)->dict:
    from Numerical.orchestration import RunSmallLambdaRefit as base
    v=np.asarray(values,dtype=float)
    require(v.shape==(18,) and np.isfinite(v).all(),'expected 18 finite real Yukawa entries')
    changed=dict(source)
    changed['cfg']=deepcopy(source['cfg'])
    for idx,k in enumerate(('y1','y2')):
        base.setmat(changed['cfg']['base_state']['ordinary'][k],
                    v[idx*9:(idx+1)*9].reshape(3,3).astype(complex))
    return changed


def beta_bundle(model:str, source:dict, registry:dict)->dict:
    """Validate and cache input beta expressions; never drop extra beta keys."""
    prepared=b9.prepare_initial(source,registry,'general_quartics',0.1,0.0)
    _,layout,_,_=prepared
    new,reused=source['new'],source['reused']
    files={stage:b7.beta_path(model,stage) for stage in ('UV','EFT1')}
    payloads={stage:b7.read_json(path) for stage,path in files.items()}
    # The EFT1 layout is derived by actual fermion-threshold projection below.
    uv=b7.check_payload(payloads['UV'],model=model,stage='UV',rep=layout.representation,
           base_names=tuple(layout.names[:-len(new)]),new_names=new,reused_names=reused)
    return {'files':files,'payloads':payloads,'UV':uv}


def integrate_point(model:str,record:dict,registry:dict, bundle:dict,
                    *,quartic_fraction:float,mu2_fraction:float)->tuple[dict,Any,Any]:
    """Every trial integrates full native candidate state, not historical RG curves."""
    from Numerical.experimental.GeneralT3ScalarRGState import (
       unpack_general_scalar_state,project_general_after_fermion,pack_general_scalar_state)
    y0,ly0,extra,mu0=b9.prepare_initial(record,registry,
         'general_quartics_and_mixing' if mu2_fraction>0 else 'general_quartics',
         quartic_fraction,mu2_fraction)
    new,reused=record['new'],record['reused']
    if record['uv_on']:
        uv=b7.integrate_segment(y0,
             lambda t,y:b7.rhs_at(t,y,layout=ly0,registry=registry,
                                   betas=bundle['UV'],stage='UV'),1e7,1e5,n_points=3)
        yf=uv.y[:,-1]; uv_nfev=int(uv.nfev)
    else:
        yf=y0.copy(); uv_nfev=0
    heavy_state=unpack_general_scalar_state(yf[:ly0.size],ly0,registry,
                                            uv=True,mu_gev=1e5).baseline
    # Heavy-mass basis must be evaluated AT muF, especially if UV running is ON.
    intermediate=project_general_after_fermion(
         unpack_general_scalar_state(yf[:ly0.size],ly0,registry,uv=True,mu_gev=1e5),registry)
    ye,lye=pack_general_scalar_state(intermediate,registry)
    ye=np.concatenate((ye,yf[ly0.size:]))
    betas=b7.check_payload(bundle['payloads']['EFT1'],model=model,stage='EFT1',
         rep=lye.representation,base_names=tuple(lye.names[:-len(new)]),
         new_names=new,reused_names=reused)
    eft=b7.integrate_segment(ye,
         lambda t,y:b7.rhs_at(t,y,layout=lye,registry=registry,
                               betas=betas,stage='EFT1'),1e5,1e3,n_points=3)
    end=b7.snapshot(eft.y[:,-1],lye,registry,stage='EFT1',mu_gev=1e3,
                    new=new,reused=reused)
    eft_state=unpack_general_scalar_state(eft.y[:-2,-1],lye,registry,
                                          uv=False,mu_gev=1e3).baseline
    return {'eft1_end':end,'uv_nfev':uv_nfev,'eft1_nfev':int(eft.nfev),
            'boundary_extra_quartic_fraction':quartic_fraction,
            'boundary_mu2_fraction':mu2_fraction},heavy_state,eft_state


def heavy_flavour_at_muF(state:Any)->tuple[np.ndarray,np.ndarray,np.ndarray,float]:
    from Numerical.running.FinalC5TrajectoryAdapter import _heavy_mass_basis_data
    x=_heavy_mass_basis_data(state,offdiagonal_atol=1e-7,imaginary_atol=1e-9,
                             takagi_residual_rtol=1e-8,takagi_residual_atol=1e-7)
    require(np.all(np.asarray(x.masses_gev)>0),'invalid heavy masses after mass-basis rotation')
    return np.asarray(x.masses_gev),np.asarray(x.y1),np.asarray(x.y2),float(max(x.takagi_residual,x.basis_residual))


def compile_source(source:Path):
    from scripts.RunGeneralT3NumericalC5 import compile_coefficient
    require(source.is_file(),'missing saved Matchete coefficient: '+str(source))
    return compile_coefficient(source.read_text(encoding='utf-8'))


def c5_at_endpoint(model:str,endpoint:dict,heavy:tuple,compiled)->np.ndarray:
    """Saved Matchete ordered coefficient with *conditional* signed direct maps."""
    import mpmath as mp
    from Numerical.experimental.GeneralT3MixedScalarMassBasis import MixedDoubletMassBasis
    info=b10.endpoint({'eft1_end':endpoint})
    basis=MixedDoubletMassBasis.build(model,info['m1'],info['m2'],info['mu'])
    require(not basis.tachyonic and not basis.degenerate and
            float(min(basis.eigenvalues_sq))>0,'mixed scalar mass eigenvalues invalid')
    mf,y1,y2,_=heavy
    require(np.asarray(mf).shape==(3,) and
            np.asarray(y1).shape==(3,3) and np.asarray(y2).shape==(3,3),
            'heavy basis flavour shape mismatch')
    mp.mp.dps=45
    c=lambda z:mp.mpc(str(complex(z).real),str(complex(z).imag))
    env={f'T3MEU{i+1}{j+1}':c(basis.U[i,j]) for i in range(2) for j in range(2)}
    env.update({'MSX':mp.mpf(str(math.sqrt(basis.eigenvalues_sq[0]))),
                'MSY':mp.mpf(str(math.sqrt(basis.eigenvalues_sq[1]))),
                'lambdaT3':c(info['lambdaT3']),
                'hbar':1/(16*mp.pi**2),'conj':mp.conj,'sqrt':mp.sqrt,'log':mp.log})
    for name,symbol in b7.DIRECT.items():
        env[symbol]=c(-info['quartics'][name])  # Stage23/26 CONDITIONAL map.
    A=[[mp.mpc(0) for _ in range(3)] for _ in range(3)]
    for i in range(3):
        for j in range(3):
            for k in range(3):
                local=dict(env,MF=mp.mpf(str(mf[k])),
                          y1a=c(y1[i,k]),y1b=c(y1[j,k]),
                          y2a=c(y2[i,k]),y2b=c(y2[j,k]))
                value=eval(compiled,{'__builtins__':{}},local)
                require(bool(mp.isfinite(value.real) and mp.isfinite(value.imag)),
                        'nonfinite ordered one-scale coefficient')
                A[i][j]+=value
    C=np.array([[complex(A[i][j]+A[j][i]) for j in range(3)]
                for i in range(3)],dtype=complex)
    require(np.isfinite(C.real).all() and np.isfinite(C.imag).all()
            and np.allclose(C,C.T,rtol=1e-10,atol=1e-17),
            'invalid symmetric C5')
    return C


def low_energy_fit(c5:np.ndarray,eft_state:Any,target:Any,*,vev:float=246.22):
    from Numerical.running.ScalarThresholdBoundary import (
         build_final_sm_boundary,build_sm_weinberg_initial_conditions)
    from Numerical.running.SMWeinbergEvolution import evolve_sm_weinberg
    from Numerical.fitting.ParameterScan import neutrino_mass_from_c5,charged_lepton_mass_basis_matrix
    from Numerical.fitting.OscillationFit import evaluate_mass_matrix_fit
    boundary=build_final_sm_boundary(eft_state,matching_scale_gev=1e3)
    init=build_sm_weinberg_initial_conditions(boundary,c5)
    low=evolve_sm_weinberg(init,1e3,1e2,rtol=1e-8,atol=1e-11)
    mass=neutrino_mass_from_c5(low.K,vev_gev=vev)
    flavour=charged_lepton_mass_basis_matrix(mass,low.ye)
    fit=evaluate_mass_matrix_fit(flavour,target)
    require(math.isfinite(float(fit.chi2)),'nonfinite five-observable chi2')
    return fit,low


class CandidateEvaluator:
    def __init__(self,case:str,model:str,source:dict,registry:dict,bundle:dict,
                 compiled:Any,target:Any,*,quartic_fraction:float,mu2_fraction:float):
        self.case,self.model,self.source=case,model,source
        self.registry,self.bundle,self.compiled,self.target=registry,bundle,compiled,target
        self.quartic_fraction=quartic_fraction;self.mu2_fraction=mu2_fraction
        self.calls=0;self.failed=0

    def evaluate(self,params:dict,*,detail:bool=False):
        vector=np.array([real_finite(params[k],k) for k in NAMES])
        trial=trial_record(self.source,vector)
        trajectory,uv_state,eft_state=integrate_point(self.model,trial,self.registry,self.bundle,
               quartic_fraction=self.quartic_fraction,mu2_fraction=self.mu2_fraction)
        heavy=heavy_flavour_at_muF(uv_state)
        c5=c5_at_endpoint(self.model,trajectory['eft1_end'],heavy,self.compiled)
        fit,low=low_energy_fit(c5,eft_state,self.target,
                         vev=float(self.source['cfg'].get('numerical',{}).get('vev_gev',246.22)))
        self.calls+=1
        if detail:
            return fit,{'trajectory':trajectory,'MF_muF_gev':heavy[0].tolist(),
                'heavy_mass_basis_residual':heavy[3],
                'C5_muS_gev_inv':b10.matrix_json(c5),
                'C5_low_gev_inv':b10.matrix_json(low.K)}
        return fit

    def __call__(self,point):
        return self.evaluate(point)


def result_record(fit:Any,point:dict,source:str)->dict:
    return {'source':source,'chi2':float(fit.chi2),
            'parameters':{k:float(point[k]) for k in NAMES},
            'prediction':dict(zip(fit.observable_names,map(float,fit.prediction))),
            'ordering':str(fit.ordering),
            'whitened_residual':np.asarray(fit.whitened_residual,dtype=float).tolist()}


def local_parameters(seed:np.ndarray,bound:float,fraction:float,floor:float):
    from Numerical.fitting.SobolBenchmarkSearch import BenchmarkParameter
    defs=[];local=[]
    for i,name in enumerate(NAMES):
        binding=f'ordinary.{name.split("_")[0]}.real.{name[-2]}.{name[-1]}'
        defs.append(BenchmarkParameter(name,binding,-bound,bound))
        half=max(abs(seed[i])*fraction,floor)
        lo=max(-bound,float(seed[i]-half));hi=min(bound,float(seed[i]+half))
        require(lo<hi,'empty local Sobol region for '+name)
        local.append(BenchmarkParameter(name,binding,lo,hi))
    return defs,local


def inputs_manifest(source:dict,model:str,case:str,bundle:dict,coefficient:Path,*,
                    quartic_fraction:float,mu2_fraction:float,seed:int,
                    sobol_points:int,bound:float,fraction:float,floor:float,seed_count:int)->dict:
    files=[source['cfg_path'],source['fit_path'],b7.REGISTRY,
           b10.EVIDENCE23,b10.EVIDENCE26,coefficient,
           bundle['files']['UV'],bundle['files']['EFT1']]
    from Numerical.orchestration import RunSmallLambdaRefit as base
    files.append(base.TARGET)
    return {'format':FORMAT,'case':case,'model':model,
           'files_sha256':{Path(f).relative_to(ROOT).as_posix():sha256(Path(f).read_bytes()).hexdigest()
                            for f in files},
           'fit_options':{'quartic_fraction':quartic_fraction,'mu2_fraction':mu2_fraction,
                'sobol_seed':seed,'sobol_points':sobol_points,'bound':bound,
                'fraction':fraction,'floor':floor,'seed_count':seed_count},
           'conditional_map':'lambda_M = -lambda_RGBeta for two direct quartics',
           'physical_coupling_map_certified':False,
           'sequential_general_scalar_matching':False,
           'historical_input_unchanged':True}


def new_run(case:str,model:str,run_id:str,*,quartic_fraction:float,mu2_fraction:float,
            sobol_points:int,seed:int,bound:float,sobol_fraction:float,sobol_floor:float,
            seed_count:int,fit:bool,root:Path=ROOT)->dict:
    from Numerical.experimental.GeneralT3ScalarRGState import load_registry
    from Numerical.orchestration import RunSmallLambdaRefit as base
    from Numerical.fitting.OscillationFit import OscillationFitTarget
    from Numerical.fitting.SobolBenchmarkSearch import (
         generate_sobol_parameter_points,refine_sobol_candidates)
    from scripts.RunGeneralT3NumericalC5 import DEFAULT_SOURCE
    dest=fit_destination(case,run_id,root=root)
    require(model in b9.MODELS,'only B_-1/C_-1 supported')
    require(sobol_points>=2 and sobol_points&(sobol_points-1)==0,
            '--sobol-points must be power of two >= 2')
    require(1<=seed_count<=8,'--seed-count must be 1..8')
    bound=real_finite(bound,'bound')
    require(bound>0,'bound must be positive')
    sf=b9.finite_fraction(sobol_fraction,'sobol_fraction',maximum=1.0)
    floor=real_finite(sobol_floor,'floor');require(0<floor<bound,'invalid Sobol floor')
    q=b9.finite_fraction(quartic_fraction,'quartic_fraction',maximum=1.0)
    u=b9.finite_fraction(mu2_fraction,'mu2_fraction',maximum=0.2)
    registry=load_registry(b7.REGISTRY)
    source=b9.study_source(case,model,registry,root=root)
    mapping=b10.check_convention_evidence(b7.read_json(b10.EVIDENCE23),
                     b7.read_json(b10.EVIDENCE26),models=(model,))
    # Same saved coefficient location as the verified Batch10 evaluator.
    coefficient=Path(DEFAULT_SOURCE)/({'B_-1':'B_m1_final',
              'C_-1':'C_m1_final'}[model] +
              '_one_loop_probe_weinberg_ordered_coefficient.inputform.txt')
    compiled=compile_source(coefficient)
    beta=beta_bundle(model,source,registry)
    target=OscillationFitTarget.from_json(base.TARGET)
    evaluator=CandidateEvaluator(case,model,source,registry,beta,compiled,target,
                  quartic_fraction=q,mu2_fraction=u)
    seed_vector0=seed_vector(source)
    require(np.max(np.abs(seed_vector0))<bound,
            'historical Yukawa seed exceeds bound; raise --bound')
    defs,local=local_parameters(seed_vector0,bound,sf,floor)
    seed_point=dict(zip(NAMES,map(float,seed_vector0)))
    manifest=inputs_manifest(source,model,case,beta,coefficient,
          quartic_fraction=q,mu2_fraction=u,seed=seed,sobol_points=sobol_points,
          bound=bound,fraction=sf,floor=floor,seed_count=seed_count)
    first=evaluator(seed_point)
    candidates=[result_record(first,seed_point,'historical_yukawa_seed_new_general_scalar_RG')]
    failures=[]
    if fit:
        points=generate_sobol_parameter_points(local,sobol_points,seed=seed)
        for idx,point in enumerate(points):
            try:
                res=evaluator(point)
                candidates.append(result_record(res,point,f'sobol_{idx+1}'))
                candidates.sort(key=lambda x:x['chi2']);candidates=candidates[:max(8,seed_count)]
            except Exception as exc:
                failures.append({'point':idx+1,'reason':str(exc)[:300]})
        best=candidates[0]
        if best['chi2']>b9.FIT_CUTOFF:
            best,history=refine_sobol_candidates(
               evaluator=evaluator,observable_names=FIELDS,definitions=defs,
               candidates=candidates,target_chi2=b9.FIT_CUTOFF,seed_count=seed_count,
               active_counts=(5,8,12),max_nfev=(20,30,40),sensitivity_step=.02)
        else:
            history={'skipped':True,'reason':'accepted existing candidate'}
        final,details=evaluator.evaluate(best['parameters'],detail=True)
        picked=best['parameters']
    else:
        final,details=evaluator.evaluate(seed_point,detail=True)
        picked=seed_point
        history={'skipped':True,'reason':'probe mode, no optimisation'}
    accepted=bool(float(final.chi2)<=b9.FIT_CUTOFF)
    result={'format':FORMAT,'model':model,'case':case,'run_id':run_id,
        'status':'ExploratoryAcceptedFit' if accepted and fit else ('ExploratoryProbe' if not fit else 'ExploratoryPoorFit'),
        'chi2_initial':float(first.chi2),'chi2_final':float(final.chi2),
        'passes_approx_chi2_5p86':accepted,
        'historical_fit_chi2':float(source['chi2']),
        'historical_fitted_config':source['cfg_path'].relative_to(ROOT).as_posix(),
        'initial_parameters':seed_point,'best_parameters':picked,
        'prediction':dict(zip(final.observable_names,map(float,final.prediction))),
        'whitened_residual':np.asarray(final.whitened_residual,dtype=float).tolist(),
        'optimiser':history,'sobol_success_count':len(candidates)-1,
        'sobol_failure_count':len(failures),'sobol_failures':failures,
        'candidate_evaluation_count':evaluator.calls,
        'final_diagnostics':details,
        'source_manifest':manifest,
        'UV_running_enabled':source['uv_on'],
        'input_extra_quartic_fraction_of_abs_lambdaT3':q,
        'input_mu2_fraction_of_min_scalar_m2':u,
        'physical_matchete_rgbeta_map_certified':False,
        'one_scale_matchete_expression_reused':True,
        'full_sequential_general_scalar_matching_performed':False,
        'fits_not_independent_predictions':True,
        'historical_outputs_overwritten':False}
    # Persist a full, distinguishable new fitted configuration: historical file untouched.
    fitted_config=trial_record(source,np.array([picked[k] for k in NAMES],dtype=float))['cfg']
    fitted_config['benchmark_search']={'enabled':False}
    fitted_config['sensitivity']={'enabled':False}
    fitted_config['exploratory_general_scalar_fit']={
        'status':result['status'], 'case':case, 'model':model, 'run_id':run_id,
        'source_config_sha256':b7.digest(source['cfg_path']),
        'source_historical_fit_sha256':b7.digest(source['fit_path']),
        'chi2_general_scalar':float(final.chi2),
        'boundary_extra_quartic_fraction':q,
        'boundary_quadratic_mixing_fraction':u,
        'conditional_signed_direct_c5_map_only':True,
        'full_sequential_general_scalar_matching_performed':False,
        'physical_general_scalar_fit_certified':False,
        'original_historical_fit_metadata_retained_for_provenance_only':True,
    }
    dest.mkdir(parents=True,exist_ok=False)
    b7.save_json(dest/f'{b9.safe_tag(model)}_exploratory_fitted_config.json',fitted_config)
    b7.save_json(dest/f'{b9.safe_tag(model)}_exploratory_fit.json',result)
    b7.save_json(dest/'summary.json',{'format':FORMAT,'case':case,'model':model,
                 'run_id':run_id,'chi2':float(final.chi2),
                 'passes_approx_chi2_5p86':accepted,
                 'physical_c5_certified':False,
                 'historical_outputs_overwritten':False})
    (dest/'summary.md').write_text(summary_markdown(result),encoding='utf-8')
    return result


def summary_markdown(r:dict)->str:
    return '\n'.join([
       '# Exploratory general-scalar Yukawa refit (B/C)', '',
       f"Historical scenario: `{r['case']}`; model: `{r['model']}`; run: `{r['run_id']}`",'',
       f"Mode: `{r['status']}`. Candidate map lambda_M=-lambda_RGBeta, **uncertified**.",
       'RG evolution and heavy-fermion mass-basis conversion are recalculated at each trial point.',
       'One-scale Matchete C5 is evolved from 1 TeV to 100 GeV by SM+Weinberg RGEs.', '',
       f"Historical chi2 = {r['historical_fit_chi2']:.9g}",
       f"General-scalar seed chi2 = {r['chi2_initial']:.9g}",
       f"General-scalar final chi2 = {r['chi2_final']:.9g}",
       f"Approximate chi2 <= 5.86: {'YES' if r['passes_approx_chi2_5p86'] else 'NO'}",'',
       '| Observable | Fit result |','|---|---:|',
       *[f'| {k} | {v:.10g} |' for k,v in r['prediction'].items()], '',
       f"Candidate evaluations: {r['candidate_evaluation_count']}; Sobol failures: {r['sobol_failure_count']}", '',
       '**Not a certified physical general-scalar fit.** No new Matchete matching, common field-action certification, or sequential mixed-scalar threshold matching has been performed.',
       'Historical saved fits and all prior general-scalar results remain unchanged.', ''])


def main()->None:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--case',choices=tuple(b9.CASES),required=True)
    ap.add_argument('--model',choices=b9.MODELS,required=True)
    ap.add_argument('--run-id',required=True)
    ap.add_argument('--fit',action='store_true',help='Run local Sobol and existing staged refinement (default is seed probe only).')
    ap.add_argument('--quartic-fraction',type=float,default=0.1)
    ap.add_argument('--mu2-fraction',type=float,default=0.)
    ap.add_argument('--sobol-points',type=int,default=8)
    ap.add_argument('--seed',type=int,default=20261028)
    ap.add_argument('--bound',type=float,default=3.0)
    ap.add_argument('--sobol-fraction',type=float,default=0.5)
    ap.add_argument('--sobol-floor',type=float,default=1e-6)
    ap.add_argument('--seed-count',type=int,default=3)
    a=ap.parse_args()
    r=new_run(a.case,a.model,a.run_id,quartic_fraction=a.quartic_fraction,
         mu2_fraction=a.mu2_fraction,sobol_points=a.sobol_points,seed=a.seed,
         bound=a.bound,sobol_fraction=a.sobol_fraction,sobol_floor=a.sobol_floor,
         seed_count=a.seed_count,fit=a.fit)
    print(f"{r['case']}: {r['model']} {r['status']} chi2={r['chi2_final']:.8g}")
    print(f"Saved: {OUT_ROOT/r['case']/r['run_id']/'summary.md'}")

if __name__=='__main__':
    main()
