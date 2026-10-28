"""Reduced-scalar, zero-additional-boundary, 16-model x 3-study comparison.

Input: immutable original fitted scenario configs and existing native RGBeta/
model-specific C5 sources. New calculations only under
output/studies/reduced_scalar_sobol/<run_id>/.

For each selected model and scenario, perform Sobol exploration, then ALWAYS
run local bounded least_squares from the best viable Sobol/seed starts.
Independently evaluate scenario-specific fixed A(-2) Yukawa matrices. No extra
quartic or extra quadratic has a nonzero input boundary. Extra quartics are
nevertheless evolved in the candidate one-loop UV/EFT1 beta equations.

This is a restricted leading-one-loop-C5 approximation, not general matching:
RG-generated B/C extra quartics are not inserted into fixed-order C5; those
insertions have a second formal loop factor under the zero boundary.

Usage: python -m scripts.RunT3ReducedScalarSobolComparison --run-id pilot_v1 \\
           --models A_-2,B_-1 --sobol-points 8 --max-nfev 20 --starts 1
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
import traceback
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.optimize import least_squares
from scipy.stats import qmc

from scripts import RunT3GeneralScalar48Campaign as bulk
from scripts import Stage31CoupledUVRunning as s31
from scripts import Stage30CoupledA0UVRunning as s30

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'output/studies/reduced_scalar_sobol'
CASES = tuple(bulk.STUDIES)
FIT_CUTOFF = 5.86
FORMAT = 't3-reduced-zero-extra-sobol-shared-a2-v1'
OBS = ('delta_m21_sq_ev2','delta_m3l_sq_ev2','sin2_theta12',
       'sin2_theta13','sin2_theta23')


def require(condition, reason):
    if not condition:
        raise ValueError('Reduced scalar Sobol: '+reason)


def valid_run_id(value):
    require(bool(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,63}',value)),
            'run-id must contain letters, digits, dash or underscore')
    return value


def choose_csv_option(value, available):
    if not value:
        return tuple(available)
    chosen=tuple(x.strip() for x in value.split(',') if x.strip())
    require(bool(chosen) and len(set(chosen))==len(chosen) and
            set(chosen)<=set(available),'invalid or duplicate selection: '+str(value))
    return chosen


def fit_status(best):
    if best is None:
        return 'NoMassSafeNeutrinoFit'
    return ('Accepted' if best['chi2'] <= FIT_CUTOFF else 'PoorFit')


def local_bounds(seed, fraction):
    from scripts import FitT3BPlus1GeneralThreeStudies as b21
    return b21.search_box(np.asarray(seed,dtype=float),
         local_fraction=fraction,floor_fraction=.05,
         abs_yukawa_max=math.sqrt(4*math.pi))


def initial_trials(seed,lo,hi,sobol_count,random_seed):
    """Exactly one old seed plus 2^k scrambled Sobol starts; no optimisation."""
    require(sobol_count>=2 and (sobol_count & (sobol_count-1))==0,
            'Sobol points must be a power of two >= 2')
    require(np.shape(seed)==(18,) and np.shape(lo)==(18,) and np.shape(hi)==(18,),
            'Yukawa search requires 18 real components')
    require(np.all(np.isfinite(seed)) and np.all(hi>lo),
            'nonfinite or reversed Yukawa search bounds')
    draws=qmc.Sobol(d=18,scramble=True,seed=random_seed).random_base2(
           int(math.log2(sobol_count)))
    return [('historical_seed',np.asarray(seed,dtype=float))]+[
        (f'sobol_{i+1}',lo+(hi-lo)*z) for i,z in enumerate(draws)]


def rank_valid_trials(rows,starts):
    """Only genuinely evaluated neutrino candidates are eligible for fitting."""
    good=[r for r in rows if r.get('status')=='NeutrinoEvaluated'
          and r.get('chi2') is not None and math.isfinite(float(r['chi2']))
          and all(k in r.get('parameters',{}) for k in bulk.NAMES)]
    return sorted(good,key=lambda x:x['chi2'])[:starts]


def write_json(path,data):
    path=Path(path)
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('x',encoding='utf-8') as f:
        json.dump(data,f,indent=2,allow_nan=False)
        f.write('\n')


def source_meta(paths):
    return {str(Path(p).relative_to(ROOT)).replace('\\','/'):bulk.sha(p)
            for p in paths}


def verify_existing_beta(model,entry,new,reused,base,*,wolframscript,timeout):
    """Read already-produced native EFT1 beta exports; export B/C if absent.

    Existing bulk export was deliberately never produced for B/C(-1); their
    single new vertex and six reused structural quartics require native export.
    """
    old=ROOT/'output/studies/general_scalar_extended/bulk_48/bulk48_v1/betas'/(
        bulk.safe(model)+'.json')
    if old.is_file():
        data=bulk.read(old)
        require(data.get('model')==model and data.get('status')=='Success' and
            tuple(data.get('new_quartic_symbols',[]))==new and
            set(data.get('reused_quartic_symbols',[]))==set(reused),
            'prior native EFT1 beta export metadata changed: '+model)
        return data,old
    output=base/'betas'/f'{bulk.safe(model)}.json'
    data,path=bulk.export_betas(model,entry,new,reused,base,wolframscript,
                                 timeout,resume=output.exists())
    return data,path


class OfficialLegacyKernel:
    """Use the authoritative legacy evaluator; do not use failing fast parity.

    Legacy one-scale and hierarchical hard+direct results remain incomparable
    as matching *schemes*. This does not magically repair legacy matching.
    """
    def __init__(self,path):
        self.path=Path(path)

    def parity_check(self,uv,eft):
        return None  # no surrogate kernel to cross-compare

    def evaluate(self,uv,eft):
        from Numerical.running.FinalC5TrajectoryAdapter import extract_final_c5_inputs
        from Numerical.running.SMWeinbergStage import evaluate_symbolic_c5
        view=SimpleNamespace(uv_threshold_state=uv,eft1_threshold_state=eft,
                             mu_scalar_threshold_gev=1e3)
        c5=np.asarray(evaluate_symbolic_c5(
            self.path,extract_final_c5_inputs(view).as_config()),dtype=complex)
        require(c5.shape==(3,3) and np.isfinite(c5.real).all() and
            np.isfinite(c5.imag).all() and
            np.allclose(c5,c5.T,rtol=1e-8,atol=1e-17),
            'legacy C5 nonfinite, wrong shape or nonsymmetric')
        return c5


def matching(model):
    kernel,path,kind=bulk.matching(model)
    if kind=='legacy-common-scale':
        kernel=OfficialLegacyKernel(path)
    return kernel,path,kind


class ReducedEvaluator:
    """Reuses checked bulk scalar ODEs with an EXACT zero-extra input boundary.

    B/C reused quartics are explicitly cleared in the ordinary state. Native
    new quartics are zero complex sidecars. Nonzero running values are evolved,
    not artificially projected to zero at every integration step.
    """
    def __init__(self,case,model,entry,new,reused,hist,payload,kernel,target):
        self.case,self.model=case,model
        self.entry,self.new,self.reused=entry,new,reused
        self.hist,self.payload,self.kernel,self.target=hist,payload,kernel,target
        self.q0=tuple(0j for _ in new)
        self.reused_boundary={key:0j for key in reused}
        self.count=self.mass_safe=self.screen_pass=self.failures=0
        self.best=self.parity=None
        self.problems=[]
        self.uvspec=None
        cfg=bulk.update_yukawas(hist,hist['seed'])
        vec,layout=s31.prepare_state(cfg,entry,reused,self.reused_boundary)
        if hist['uv_on']:
            if model=='A_+0':
                p=bulk.BETA_DIR/'stage30_A_p0_full_uv_betas.json'
                beta=s30.validate_export(bulk.read(p),new,layout.names)
            else:
                p=bulk.BETA_DIR/f'stage31_{bulk.safe(model)}_full_uv_betas.json'
                beta=s31.validate_export(bulk.read(p),entry,new,reused,layout.names)
            self.uvspec=(beta,p)
        else:
            from Numerical.core.StateVector import unpack_uv_state
            uv=unpack_uv_state(vec,layout,mu_gev=1e5)
            self.frozen=bulk.run_eft(uv,self.q0,payload,model,new,reused)

    def evaluate(self,x,label):
        self.count+=1
        row={'label':label,'status':'Rejected','chi2':None,'reason':'',
             'prediction':None}
        try:
            cfg=bulk.update_yukawas(self.hist,x)
            vec,layout=s31.prepare_state(cfg,self.entry,self.reused,self.reused_boundary)
            from Numerical.core.StateVector import unpack_uv_state
            if self.hist['uv_on']:
                uv,qf,nuv=bulk.run_uv(vec,layout,self.uvspec[0],self.new,self.q0)
                eft,masses,qs,neft=bulk.run_eft(
                    uv,qf,self.payload,self.model,self.new,self.reused)
            else:
                uv=unpack_uv_state(vec,layout,mu_gev=1e5)
                eft,masses,qs,neft=self.frozen
                qf=self.q0;nuv=0
            chk=bulk.testscreen(uv,eft,masses,qf,qs)
            row.update(chk)
            if chk['screen']:
                self.screen_pass+=1
                if self.parity is None and not isinstance(self.kernel,OfficialLegacyKernel):
                    self.parity=self.kernel.parity_check(uv,eft)
                c5=self.kernel.evaluate(uv,eft)
                from scripts.RunT3GeneralScalarBCExploratoryFit import low_energy_fit
                from scripts import FitT3BPlus1GeneralThreeStudies as b21
                fit,_=low_energy_fit(c5,eft,self.target,
                    vev=float(cfg.get('numerical',{}).get('vev_gev',246.22)))
                residual=b21.residual_from_fit(fit)
                row.update({'status':'NeutrinoEvaluated','chi2':float(fit.chi2),
                    'prediction':dict(zip(fit.observable_names,map(float,fit.prediction))),
                    'whitened_residual':residual.tolist(),
                    'parameters':dict(zip(bulk.NAMES,map(float,x))),
                    'uv_nfev':nuv,'eft1_nfev':neft})
                self.mass_safe+=1
                if self.best is None or row['chi2']<self.best['chi2']:
                    self.best=deepcopy(row)
            if row['chi2'] is None and len(self.problems)<8:
                self.problems.append({'label':label,'reason':row['reason']})
        except (ValueError,TypeError,RuntimeError,ArithmeticError,OverflowError,
                FloatingPointError,np.linalg.LinAlgError) as exc:
            self.failures+=1
            row['status']='NumericalFailure';row['reason']=str(exc)[:500]
            if len(self.problems)<8:
                self.problems.append({'label':label,'reason':row['reason']})
        return row


def a2_reference(case,registry,source,base):
    """Return strictly scenario-matched A(-2) Yukawas; no silent fallback."""
    entry=registry['models']['A_-2']
    hist=bulk.historical(case,'A_-2',entry)
    if source=='historical':
        return np.asarray(hist['seed'],dtype=float),{
            'source':'historical_A_-2_scenario_fit',
            'cfg_path':str(hist['cfg_path'].relative_to(ROOT)),
            'cfg_sha256':bulk.sha(hist['cfg_path']),
            'source_fit_sha256':bulk.sha(hist['fit_path'])}
    saved=base/case/bulk.safe('A_-2')/'fit.json'
    if not saved.is_file():
        return None,{'source':'new_reduced_A_-2_fit','reason':'No A_-2 fit.json in this run'}
    obj=bulk.read(saved)
    best=obj.get('best_mass_safe') or {}
    if obj.get('status')!='Accepted' or not all(k in best.get('parameters',{})
                                             for k in bulk.NAMES):
        return None,{'source':'new_reduced_A_-2_fit',
                     'reason':'New A_-2 result not accepted'}
    return np.array([best['parameters'][k] for k in bulk.NAMES],dtype=float),{
        'source':'new_reduced_A_-2_fit','fit_path':str(saved.relative_to(ROOT)),
        'fit_sha256':bulk.sha(saved)}


def previous_best_yukawas(case, model, source_run):
    """Require a complete accepted warm-start result from an immutable prior run.

    Only its Yukawa vector is reused. The target model's scalar inputs,
    representations, beta functions and matching are always rebuilt normally.
    """
    valid_run_id(source_run)
    path=OUT/source_run/case/bulk.safe(model)/'fit.json'
    require(path.is_file(), 'missing warm-start fit: '+str(path))
    saved=bulk.read(path)
    best=saved.get('best_mass_safe') or {}
    require(saved.get('model')==model and saved.get('case')==case and
            saved.get('status')=='Accepted' and
            all(k in best.get('parameters',{}) for k in bulk.NAMES),
            'warm-start source is not an accepted matching model/scenario: '+str(path))
    y=np.asarray([best['parameters'][k] for k in bulk.NAMES],dtype=float)
    require(y.shape==(18,) and np.isfinite(y).all() and
            np.max(np.abs(y))<=math.sqrt(4*math.pi),
            'warm-start Yukawa values are nonfinite or outside perturbative search box')
    return y,{'source':'accepted_previous_reduced_fit',
              'fit_path':str(path.relative_to(ROOT)).replace('\\','/'),
              'fit_sha256':bulk.sha(path),
              'previous_chi2':float(best['chi2'])},path


def extended_bounds(seed, fraction, extra_trials):
    """Include accepted warm/A2 starts without silently clipping their Yukawas.

    Bounds remain inside the original absolute Yukawa cap, but encompass all
    supplied valid starts, including transferred A(-2) matrices.
    """
    lo,hi=local_bounds(seed,fraction)
    cap=math.sqrt(4*math.pi)
    for label,vec in extra_trials:
        x=np.asarray(vec,dtype=float)
        require(x.shape==(18,) and np.isfinite(x).all(),
                'invalid extra candidate '+str(label))
        require(np.max(np.abs(x))<cap,
                'extra candidate exceeds absolute Yukawa search cap: '+str(label))
        margin=np.maximum(1e-4,0.01*np.abs(x))
        lo=np.maximum(-cap,np.minimum(lo,x-margin))
        hi=np.minimum(cap,np.maximum(hi,x+margin))
    require(np.all(lo<hi),'extended candidate bounds invalid')
    return lo,hi


def optimise(evaluator,hist,*,sobol_points,starts,max_nfev,
             seed,local_fraction,extra_trials=()):
    lo,hi=extended_bounds(hist['seed'],local_fraction,extra_trials)
    initial=initial_trials(hist['seed'],lo,hi,sobol_points,seed)
    # Warm and shared-reference points are explicit, separately labelled trial
    # inputs. Deduplicate so two identical points do not consume local starts.
    unique=[]
    for label,x in initial[:1]+list(extra_trials)+initial[1:]:
        x=np.asarray(x,dtype=float)
        if not any(np.array_equal(x,prior) for _,prior in unique):
            unique.append((label,x))
    first=[evaluator.evaluate(x,label) for label,x in unique]
    chosen=rank_valid_trials(first,starts)
    local=[]
    for i,start in enumerate(chosen):
        x=np.array([start['parameters'][k] for k in bulk.NAMES],dtype=float)
        z0=np.clip((x-lo)/(hi-lo),1e-8,1.-1e-8)
        def residual(z):
            point=lo+(hi-lo)*z
            result=evaluator.evaluate(point,f'local_{i+1}')
            return np.array(result.get('whitened_residual') or [1e4]*5,dtype=float)
        result=least_squares(residual,z0,bounds=(np.zeros(18),np.ones(18)),
                 max_nfev=max_nfev,xtol=1e-6,ftol=1e-6,gtol=1e-6)
        local.append({'start':start['label'],'start_chi2':start['chi2'],
                      'success':bool(result.success),'message':str(result.message),
                      'nfev':int(result.nfev),'cost':float(result.cost)})
    best=deepcopy(evaluator.best)
    if best is not None:
        status=fit_status(best)
    elif getattr(evaluator,'screen_pass',0)>0:
        status='MatchingOrEvaluationFailure'
    elif evaluator.failures==evaluator.count:
        status='NumericalFailure'
    else:
        status='NoMassSafeNeutrinoFit'
    return {'status':status,'best_mass_safe':best,
        'sobol_trials':first,'optimiser_starts':local,
        'counts':{'evaluations':evaluator.count,'mass_safe':evaluator.mass_safe,
                  'screen_pass':getattr(evaluator,'screen_pass',0),
                  'numerical_failures':evaluator.failures},
        'first_problems':deepcopy(evaluator.problems),
        'extra_trial_labels':[name for name,_ in extra_trials],
        'local_optimisation_attempted':bool(local),
        'local_optimisation_skipped_reason':(
            None if local else 'No mass-safe Sobol/seed candidate')}


def row(model,case,record):
    best=record.get('best_mass_safe') or {}
    shared=record.get('shared_a2') or {}
    prediction=best.get('prediction') or {}
    sharedprediction=shared.get('prediction') or {}
    return {'model':model,'case':case,'status':record.get('status'),
        'chi2_optimised':best.get('chi2'),'chi2_shared_a2':shared.get('chi2'),
        'shared_status':shared.get('status'),
        'best_source':best.get('label'),
        'sobol_points':record.get('settings',{}).get('sobol_points'),
        'local_starts':len(record.get('optimiser_starts') or []),
        'evaluations':record.get('counts',{}).get('evaluations'),
        'mass_safe_evaluations':record.get('counts',{}).get('mass_safe'),
        'boundary_lambdaT3':record.get('lambdaT3_input'),
        'boundary_scale_gev':record.get('extra_zero_boundary_scale_gev'),
        'c5_source_kind':record.get('c5_source_kind'),
        **{f'fit_{k}':prediction.get(k) for k in OBS},
        **{f'shared_{k}':sharedprediction.get(k) for k in OBS},
        'error':record.get('error')}


def settings(args):
    return {'runner_sha256':bulk.sha(Path(__file__)),
            'bulk_runner_sha256':bulk.sha(Path(bulk.__file__)),
            'scalar_registry_sha256':bulk.sha(bulk.REGISTRY),
            'target_sha256':bulk.sha(ROOT/'configs/nufit_6_1_ic24_no.json'),
            'sobol_points':args.sobol_points,'starts':args.starts,
            'max_nfev':args.max_nfev,'seed':args.seed,
            'local_fraction':args.local_fraction,'a2_source':args.a2_source,
            'a2_as_start':args.a2_as_start,
            'warm_start_run':args.warm_start_run,
            'fit_cutoff':FIT_CUTOFF,
            'extra_quartic_boundary_real_imag':[0.,0.],
            'extra_quadratic_boundary':0.}


def process_case(model,case,registry,base,args):
    from Numerical.fitting.OscillationFit import OscillationFitTarget
    from Numerical.orchestration import RunSmallLambdaRefit as refit
    entry,new,reused=bulk.inventory(registry,model)
    hist=bulk.historical(case,model,entry)
    target=OscillationFitTarget.from_json(refit.TARGET)
    require(target.ordering=='NO' and len(target.observable_names)==5,
            'historical oscillation target changed')
    payload,beta_path=verify_existing_beta(model,entry,new,reused,base,
            wolframscript=args.wolframscript,timeout=args.timeout_sec)
    kernel,match_path,match_kind=matching(model)
    evaluator=ReducedEvaluator(case,model,entry,new,reused,hist,payload,kernel,target)
    extra_trials=[]
    extra_provenance=[]
    extra_paths=[]
    if args.warm_start_run:
        warm,warm_meta,warm_path=previous_best_yukawas(
            case,model,args.warm_start_run)
        extra_trials.append(('previous_individual_fit',warm))
        extra_provenance.append(warm_meta)
        extra_paths.append(warm_path)
    if args.a2_as_start and model!='A_-2':
        reference,ref_meta=a2_reference(case,registry,args.a2_source,base)
        require(reference is not None,
                'A(-2) candidate unavailable before optimising '+model+
                ': '+str(ref_meta.get('reason','unknown reason')))
        extra_trials.append(('A_-2_reference_fit',reference))
        extra_provenance.append(ref_meta)
        if 'fit_path' in ref_meta:
            extra_paths.append(ROOT/ref_meta['fit_path'])
        elif 'cfg_path' in ref_meta:
            extra_paths.append(ROOT/ref_meta['cfg_path'])
    seed_offset=sum((i+1)*ord(c) for i,c in enumerate(model+'|'+case))
    outcome=optimise(evaluator,hist,sobol_points=args.sobol_points,
          starts=args.starts,max_nfev=args.max_nfev,
          seed=args.seed+seed_offset,local_fraction=args.local_fraction,
          extra_trials=tuple(extra_trials))
    evaluation_counts=deepcopy(outcome['counts'])
    if args.a2_source=='new-fit' and model=='A_-2':
        own=outcome.get('best_mass_safe') or {}
        if outcome['status']=='Accepted' and all(k in own.get('parameters',{})
                                             for k in bulk.NAMES):
            a2_vec=np.array([own['parameters'][k] for k in bulk.NAMES],dtype=float)
            a2_meta={'source':'this_case_new_reduced_A_-2_fit'}
        else:
            a2_vec=None
            a2_meta={'source':'this_case_new_reduced_A_-2_fit',
                     'reason':'New A_-2 result not accepted'}
    else:
        a2_vec,a2_meta=a2_reference(case,registry,args.a2_source,base)
    if a2_vec is None:
        shared={'status':'ReferenceUnavailable','chi2':None,
                'reason':a2_meta['reason']}
    else:
        # Record separately; shared point must never change optimised best.
        test=evaluator.evaluate(a2_vec,'shared_A_-2')
        shared={'status':test['status'],'chi2':test['chi2'],
                'reason':test.get('reason'),'prediction':test.get('prediction'),
                'rhoF':test.get('rhoF'),
                'scalar_masses_muS_tev':test.get('scalar_masses_muS_tev')}
    sources=[hist['cfg_path'],hist['fit_path'],bulk.REGISTRY,beta_path,
             match_path,refit.TARGET]
    if evaluator.uvspec is not None:sources.append(evaluator.uvspec[1])
    sources.extend(extra_paths)
    outcome.update({'format':FORMAT,'model':model,'case':case,
        'lambdaT3_input':hist['lambda_fixed'],
        'extra_zero_boundary_scale_gev':1e7 if hist['uv_on'] else 1e5,
        'new_additional_quartic_symbols':list(new),
        'reused_additional_quartic_symbols_zeroed':list(reused),
        'extra_quadratic_boundary_exact_zero':True,
        'quartic_one_loop_running_included':True,
        'additional_quadratic_evolution_included':False,
        'general_scalar_one_loop_matching_included':False,
        'zero_boundary_relevance_argument':'docs/T3_REDUCED_SCALAR_POTENTIAL_CONCLUSION.md',
        'c5_source_kind':match_kind,'c5_source':str(match_path.relative_to(ROOT)),
        'shared_a2':shared,'shared_a2_source':a2_meta,
        'extra_optimisation_seed_sources':extra_provenance,
        'optimisation_counts_before_shared_evaluation':evaluation_counts,
        'matching_parity_error':evaluator.parity,
        'settings':settings(args),'source_sha256':source_meta(sources)})
    return outcome


def save_case(base,model,case,record):
    folder=base/case/bulk.safe(model)
    write_json(folder/'fit.json',record)
    best=record.get('best_mass_safe') or {}
    if best and all(k in best.get('parameters',{}) for k in bulk.NAMES):
        entry=bulk.read(bulk.REGISTRY)['models'][model]
        hist=bulk.historical(case,model,entry)
        x=np.array([best['parameters'][k] for k in bulk.NAMES],dtype=float)
        cfg=bulk.update_yukawas(hist,x)
        cfg['reduced_scalar_zero_extra_boundary']={
            'all_additional_quartics_exact_zero':True,
            'all_additional_quadratics_exact_zero':True,
            'boundary_scale_gev':record['extra_zero_boundary_scale_gev'],
            'reused_rgbeta_quartics_zeroed_in_runner':record[
                'reused_additional_quartic_symbols_zeroed'],
            'new_rgbeta_quartics_zero_sidecars_in_runner':record[
                'new_additional_quartic_symbols'],
            'NOTE':'This config alone does not encode sidecar states. Use the runner.'}
        # Keep config's reused ordinary couplings consistent with the zero UV input.
        for k in record['reused_additional_quartic_symbols_zeroed']:
            cfg['base_state']['ordinary'][k]={'real':0.,'imag':0.}
        write_json(folder/'fitted_config.json',cfg)


def summarise(base,selected_cases,selected_models,registry,args):
    rows=[]
    for case in selected_cases:
        for model in selected_models:
            path=base/case/bulk.safe(model)/'fit.json'
            if path.is_file():
                obj=bulk.read(path)
                rows.append(row(model,case,obj))
            else:
                rows.append(row(model,case,{'status':'NotRun'}))
    csvpath=base/'comparison.csv'
    with csvpath.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    meta={'format':FORMAT,'rows':rows,'cases':list(selected_cases),
          'models':list(selected_models),'settings':settings(args),
          'note':'Fixed-order leading one-loop C5 with zero-extra boundary; '
                 'not certified general matching; legacy source is nonuniform.'}
    tmp=base/'comparison.json.tmp'
    tmp.write_text(json.dumps(meta,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    tmp.replace(base/'comparison.json')
    accepted=sum(r['status']=='Accepted' for r in rows)
    shared=sum(r['chi2_shared_a2'] is not None for r in rows)
    lines=['# Reduced scalar zero-extra boundary: Sobol + A(-2) comparison','',
        f'Completed: {sum(r["status"]!="NotRun" for r in rows)}/{len(rows)}; '
        f'accepted: {accepted}; shared A(-2) evaluated: {shared}.','',
        'Additional scalar quartics and quadratic coefficients set zero at '
        'the scenario input boundary; additional quartics are evolved under '
        'one-loop RGEs but not inserted into strict fixed-order leading C5.',
        'No matching-convention certification; 4 legacy C5 source schemes differ.',
        '', '| Model | Scenario | Individual chi2 | Shared A(-2) chi2 | Status |',
        '|---|---|---:|---:|---|']
    for r in rows:
        val=lambda x:'—' if x is None else f'{x:.7g}'
        lines.append(f'| {r["model"]} | {r["case"]} | '
                     f'{val(r["chi2_optimised"])} | '
                     f'{val(r["chi2_shared_a2"])} | {r["status"]} |')
    (base/'summary.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print(f'SUMMARY {accepted}/{len(rows)} accepted, {shared} shared evaluations '
          f'-> {csvpath}',flush=True)


def run(args):
    valid_run_id(args.run_id)
    require(args.sobol_points>=2 and (args.sobol_points & (args.sobol_points-1))==0,
            '--sobol-points must be a power of two >= 2')
    require(1<=args.starts<=8 and 1<=args.max_nfev<=1000,
            'starts must be 1..8 and max-nfev must be 1..1000')
    require(0<args.local_fraction<=3.,'invalid local search fraction')
    if args.warm_start_run:
        valid_run_id(args.warm_start_run)
        require(args.warm_start_run!=args.run_id,
                'warm-start source must be a different run ID')
    registry=bulk.read(bulk.REGISTRY)
    require(len(registry.get('models',{}))==16,'registry must contain 16 models')
    selected_models=choose_csv_option(args.models,registry['models'])
    selected_cases=choose_csv_option(args.cases,CASES)
    require(bool(selected_models and selected_cases),'empty selection')
    if args.a2_as_start:
        require(args.a2_source=='new-fit',
                '--a2-as-start requires --a2-source new-fit for this campaign')
    if args.a2_source=='new-fit':
        require('A_-2' in selected_models,
             'new-fit A(-2) source requires A_-2 in --models')
        selected_models=('A_-2',)+tuple(x for x in selected_models if x!='A_-2')
    base=OUT/args.run_id
    cfg_path=base/'settings.json'
    if base.exists():
        require(args.resume,'output directory exists; use --resume or new run-id')
        require(cfg_path.is_file() and bulk.read(cfg_path)=={
            'settings':settings(args),'models':list(selected_models),
            'cases':list(selected_cases)},
            'existing run settings/selection differ; use new run-id')
    else:
        base.mkdir(parents=True,exist_ok=False)
        write_json(cfg_path,{'settings':settings(args),
                            'models':list(selected_models),
                            'cases':list(selected_cases)})
    for case in selected_cases:
        for model in selected_models:
            dest=base/case/bulk.safe(model)/'fit.json'
            if dest.is_file():
                print('RESUME',model,case,flush=True)
                continue
            try:
                obj=process_case(model,case,registry,base,args)
            except Exception as exc:
                obj={'format':FORMAT,'model':model,'case':case,
                    'status':'FailedOrUnavailable','error':str(exc)[:1200],
                    'settings':settings(args),
                    'traceback':traceback.format_exc(limit=8)}
            # An exception never fabricates a neutrino chi2.
            save_case(base,model,case,obj)
            print(model,case,obj['status'],
                  (obj.get('best_mass_safe') or {}).get('chi2'),flush=True)
            summarise(base,selected_cases,selected_models,registry,args)
    summarise(base,selected_cases,selected_models,registry,args)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--run-id',required=True)
    ap.add_argument('--resume',action='store_true')
    ap.add_argument('--models',help='Comma-separated model keys; default all 16')
    ap.add_argument('--cases',help='Comma-separated scenario keys; default all 3')
    ap.add_argument('--sobol-points',type=int,default=32)
    ap.add_argument('--starts',type=int,default=2,
                    help='Best valid Sobol/seed candidates for local optimisation')
    ap.add_argument('--max-nfev',type=int,default=60)
    ap.add_argument('--local-fraction',type=float,default=1.5)
    ap.add_argument('--seed',type=int,default=20261028)
    ap.add_argument('--a2-source',choices=('historical','new-fit'),
                    default='historical',help='Scenario-matched shared Yukawa source')
    ap.add_argument('--warm-start-run',
                    help='Previous accepted reduced-scalar run ID; use its fitted Yukawas as extra starts')
    ap.add_argument('--a2-as-start',action='store_true',
                    help='Evaluate fitted A(-2) Yukawas as an individual optimisation start in each other model')
    ap.add_argument('--wolframscript',default='wolframscript')
    ap.add_argument('--timeout-sec',type=int,default=3600)
    args=ap.parse_args()
    run(args)

if __name__=='__main__':main()
