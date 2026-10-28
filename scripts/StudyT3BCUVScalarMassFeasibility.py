"""Batch 15: conditional B/C UV scalar-mass feasibility screen; NO neutrino fit.

Reuses the fitted, source-hash-checked Batch 14 inputs, native UV and EFT1
RGBeta and their established DOP853 solver. Scans positive UV scalar masses
and a common scale of the 18 saved real Yukawas; MF and scalar quartic
boundary prescription are unchanged. Skips EFT1 if the formal UV endpoint
already has inverted running scalar/fermion ordering. A small beta-sector
finite-substitution identity identifies the local heavy Y/MF dependence.

Every result is candidate RG physics. Neither a physical threshold nor a
new neutrino fit is calculated. No historical files are modified.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
import csv
from hashlib import sha256
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from scripts import RunT3GeneralScalarBCPilot as b9
from scripts import RunT3GeneralScalarBCExploratoryFit as b11
from scripts import DiagnoseT3BCUVScalarHierarchy as b14
from scripts import RunT3ZeroExtraBoundaryBenchmark as b7

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'output/studies/general_scalar_extended/diagnostics'
FORMAT = 't3-bc-uv-positive-mass-feasibility-screen-v1'
CASES = ('B_-1','C_-1')
MASS_KEYS = ('mS1Sq','mS2Sq')
CSV_FIELDS = ('model','yukawa_factor','uv_scalar_mass_factor','status',
              'uv_scalar1_tev','uv_scalar2_tev','scalar_light_muF_tev',
              'scalar_heavy_muF_tev','min_fermion_muF_tev','ratio_muF',
              'scalar_light_muS_tev','scalar_heavy_muS_tev',
              'ordering_at_muF','strict_hierarchy_at_muF','near_1tev_at_muS',
              'screen_pass','uv_nfev','eft1_nfev','reason')


def require(ok: bool, message: str) -> None:
    if not ok:
        raise ValueError('T3 B/C UV feasibility: '+message)


def parse_factors(raw: str, label: str, *, maximum: float=10.) -> tuple[float,...]:
    try:
        entries=raw.split(',')
        result=tuple(float(z.strip()) for z in entries)
    except (ValueError, AttributeError) as exc:
        raise ValueError(label+' must be comma-separated positive numbers') from exc
    require(bool(result) and len(result)<=12 and all(math.isfinite(z) and 0<z<=maximum for z in result),
            label+' must contain 1–12 finite positive factors')
    require(len(set(result))==len(result),label+' contains duplicate factors')
    require(1. in result,label+' must include 1.0 for the historical-input replay')
    return tuple(sorted(result))


def candidate_record(source:dict, saved:dict, yukawa_factor:float,
                     scalar_mass_factor:float) -> dict:
    """Only mutate a deep copy, preserving the historical 1.14:1 scalar ratio."""
    require(math.isfinite(yukawa_factor) and 0<yukawa_factor<=10.,'invalid Yukawa factor')
    require(math.isfinite(scalar_mass_factor) and 0<scalar_mass_factor<=10.,
            'invalid positive scalar mass factor')
    values=saved['best_parameters']
    original=np.asarray([float(values[k]) for k in b11.NAMES],dtype=float)
    require(original.shape==(18,) and np.isfinite(original).all(),'invalid saved Yukawas')
    record=b11.trial_record(source,original*yukawa_factor)
    ordinary=record['cfg']['base_state']['ordinary']
    for name in MASS_KEYS:
        historic=float(source['cfg']['base_state']['ordinary'][name])
        require(math.isfinite(historic) and historic>0,'nonpositive historical '+name)
        ordinary[name]=historic*scalar_mass_factor**2
    return record


def scalar_spectrum(vector:np.ndarray, layout:Any, registry:dict,
                    *,scale:float,stage:str, heavy_at_F:tuple[float,...]|None=None)->dict:
    from Numerical.experimental.GeneralT3ScalarRGState import unpack_general_scalar_state
    state=unpack_general_scalar_state(np.asarray(vector[:layout.size],float),layout,registry,
                                      uv=(stage=='UV'),mu_gev=scale)
    m1,m2=float(state.baseline.mS1Sq),float(state.baseline.mS2Sq)
    mu2=complex(float(vector[layout.size]),float(vector[layout.size+1]))
    eig=np.linalg.eigvalsh(np.array([[m1,mu2],[mu2.conjugate(),m2]],dtype=complex))
    # Preserve numerical negatives in the report, rather than taking abs(sqrt).
    masses=[math.sqrt(float(x)) if x>0 else None for x in eig]
    if stage=='UV':
        f=np.linalg.svd(np.asarray(state.baseline.MF,complex),compute_uv=False)
    else:
        require(heavy_at_F is not None,'missing threshold heavy masses')
        f=np.asarray(heavy_at_F,dtype=float)
    require(np.isfinite(f).all() and np.min(f)>0,'invalid positive F spectrum')
    return {'mS1Sq_gev2':m1,'mS2Sq_gev2':m2,
            'scalar_eigen_mass_squared_gev2':list(map(float,eig)),
            'scalar_masses_gev':masses,'min_fermion_mass_gev':float(np.min(f)),
            'max_fermion_mass_gev':float(np.max(f)),
            'positive_scalar_masses':bool(np.min(eig)>0),
            'mu2_abs_gev2':float(abs(mu2))}


def beta_sector_at_uv(vector:np.ndarray,layout:Any,registry:dict,betas:dict,
                      *,scale:float)->dict:
    """Exact inclusion/exclusion in Y and MF at one frozen RG point.

    This measures *algebraic dependence* of the mass beta expression on these
    variables, not a diagrammatic or gauge-invariant operator decomposition.
    """
    from Numerical.experimental.GeneralT3ScalarRGState import unpack_general_scalar_state
    from Numerical.core.RGBetaEvaluator import state_environment,evaluate_inputform_expression
    state=unpack_general_scalar_state(np.asarray(vector[:layout.size],float),layout,registry,
                                      uv=True,mu_gev=scale)
    env=state_environment(state.baseline)
    env.update(state.extra)
    env[b7.MASS_SYMBOL]=complex(float(vector[layout.size]),float(vector[layout.size+1]))
    all_y={'y1':np.zeros_like(env['y1']),'y2':np.zeros_like(env['y2'])}
    all_m={'MF':np.zeros_like(env['MF'])}

    def eval_pair(y_on:bool,m_on:bool)->dict:
        e=dict(env)
        if not y_on:e.update(all_y)
        if not m_on:e.update(all_m)
        result={}
        for name in MASS_KEYS:
            value=evaluate_inputform_expression(betas[name],e)
            require(not isinstance(value,np.ndarray) or value.size==1,
                    'scalar mass beta unexpectedly matrix valued')
            z=complex(np.asarray(value).item())
            require(math.isfinite(z.real) and math.isfinite(z.imag) and
                    abs(z.imag)<=1e-7*max(1.,abs(z.real)),
                    'non-real/non-finite scalar mass beta '+name)
            result[name]=float(z.real/(16*math.pi**2))
        return result

    whole=eval_pair(True,True)
    neither=eval_pair(False,False)
    y_only=eval_pair(True,False)
    m_only=eval_pair(False,True)
    sectors={}
    for name in MASS_KEYS:
        v=whole[name]; b=neither[name]
        y=y_only[name]-b
        m=m_only[name]-b
        mixed=v-y_only[name]-m_only[name]+b
        require(math.isclose(b+y+m+mixed,v,rel_tol=1e-12,abs_tol=1e-5),
                'beta inclusion/exclusion closure failed')
        sectors[name]={
            'full_d_m2_d_lnmu_gev2':v,
            'no_y_no_MF':b,
            'y_dependence_at_MF0':y,
            'MF_dependence_at_y0':m,
            'joint_y_MF_dependence':mixed,
            'total_heavy_related':v-b,
            'identity_verified':True}
    return {'mu_gev':scale,'terms':sectors,
            'interpretation':'Frozen-state algebraic Y/MF dependence, NOT unique Feynman-diagram beta classes'}


def evaluate_uv(model:str,record:dict,registry:dict,bundle:dict,
                *,q:float,mix:float,grid:int)->dict:
    from Numerical.experimental.GeneralT3ScalarRGState import unpack_general_scalar_state
    y0,layout,_,_=b9.prepare_initial(record,registry,
        'general_quartics_and_mixing' if mix>0 else 'general_quartics',q,mix)
    rhs=lambda t,y:b7.rhs_at(t,y,layout=layout,registry=registry,
                              betas=bundle['UV'],stage='UV')
    uv=b7.integrate_segment(y0,rhs,1e7,1e5,n_points=grid)
    initial=scalar_spectrum(y0,layout,registry,scale=1e7,stage='UV')
    yf=np.asarray(uv.y[:,-1],dtype=float)
    final=scalar_spectrum(yf,layout,registry,scale=1e5,stage='UV')
    return {'initial':initial,'final':final,'first_vector':y0,
            'final_vector':yf,'layout':layout,'uv_nfev':int(uv.nfev)}


def evaluate_eft1(model:str,uv:dict,record:dict,registry:dict,bundle:dict,
                 *,grid:int)->dict:
    from Numerical.experimental.GeneralT3ScalarRGState import (
        unpack_general_scalar_state,project_general_after_fermion,pack_general_scalar_state)
    layout=uv['layout'];yf=uv['final_vector']
    state=unpack_general_scalar_state(yf[:layout.size],layout,registry,
                                      uv=True,mu_gev=1e5)
    f=tuple(sorted(float(x) for x in np.linalg.svd(state.baseline.MF,compute_uv=False)))
    after_f=project_general_after_fermion(state,registry)
    vec,eft_layout=pack_general_scalar_state(after_f,registry)
    vec=np.concatenate((vec,yf[layout.size:]))
    new,reused=record['new'],record['reused']
    betas=b7.check_payload(bundle['payloads']['EFT1'],model=model,stage='EFT1',
            rep=eft_layout.representation,
            base_names=tuple(eft_layout.names[:-len(new)]),
            new_names=new,reused_names=reused)
    rhs=lambda t,y:b7.rhs_at(t,y,layout=eft_layout,registry=registry,
                              betas=betas,stage='EFT1')
    sol=b7.integrate_segment(vec,rhs,1e5,1e3,n_points=grid)
    final=scalar_spectrum(sol.y[:,-1],eft_layout,registry,
                          scale=1e3,stage='EFT1',heavy_at_F=f)
    return {'final':final,'eft1_nfev':int(sol.nfev)}


def classify(uv:dict,eft:dict|None,*,ratio_max:float,
             scalar_min_tev:float,scalar_max_tev:float)->dict:
    end=uv['final']; masses=end['scalar_masses_gev'];
    good=bool(end['positive_scalar_masses'] and all(x is not None for x in masses))
    ratio=(float(max(masses)/end['min_fermion_mass_gev']) if good else None)
    ordering=bool(good and ratio<1)
    strict=bool(ordering and ratio<=ratio_max)
    low=None;close=False
    if eft is not None:
        low=eft['final']['scalar_masses_gev']
        close=bool(eft['final']['positive_scalar_masses'] and
                   all(x is not None and scalar_min_tev*1000<=x<=scalar_max_tev*1000
                       for x in low))
    passed=bool(strict and close)
    return {'ratio_muF':ratio,'ordering_at_muF':ordering,
            'strict_hierarchy_at_muF':strict,'near_1tev_at_muS':close,
            'screen_pass':passed,
            'reason':('CandidatePassedScreenOnly' if passed else
                      'NonPositiveScalarAtMuF' if not good else
                      'ScalarAtOrAboveFAtMuF' if not ordering else
                      'HierarchyMarginNotReached' if not strict else
                      'ScalarNotInTargetWindowAtMuS')}


def model_scan(*,model:str,source:dict,saved:dict,settings:dict,registry:dict,
               yukawa_factors:tuple[float,...],mass_factors:tuple[float,...],
               ratio_max:float,scalar_min_tev:float,scalar_max_tev:float,
               grid:int)->dict:
    baseline=candidate_record(source,saved,1.,1.)
    bundle=b11.beta_bundle(model,baseline,registry)
    trials=[];beta={};passes=[]
    # Deterministic baseline first, so a failed historical replay stops before
    # a misleading scan is written. Other coordinates are numerical probes.
    combos=[(1.,1.)]+[(y,s) for y in yukawa_factors for s in mass_factors
                       if (y,s)!=(1.,1.)]
    for y,s in combos:
        original=(y,s)==(1.,1.)
        row={'model':model,'yukawa_factor':y,'uv_scalar_mass_factor':s,
             'saved_chi2_not_reevaluated':float(saved['chi2_final']),
             'new_oscillation_fit_performed':False}
        try:
            trial=candidate_record(source,saved,y,s)
            uv=evaluate_uv(model,trial,registry,bundle,
                q=settings['quartic_fraction'],mix=settings['mu2_fraction'],grid=grid)
            start,end=uv['initial'],uv['final']
            if original:
                # This is an independent UV segment replay, not just using
                # stored endpoint numbers as the new result.
                f=np.sort(np.asarray(saved['final_diagnostics']['MF_muF_gev'],float))
                from Numerical.experimental.GeneralT3ScalarRGState import unpack_general_scalar_state
                state=unpack_general_scalar_state(uv['final_vector'][:uv['layout'].size],
                         uv['layout'],registry,uv=True,mu_gev=1e5)
                got=np.sort(np.linalg.svd(state.baseline.MF,compute_uv=False))
                require(np.allclose(f,got,rtol=1e-5,atol=1e-4),
                        'baseline heavy masses failed saved-fit replay')
                beta={'uv_start':beta_sector_at_uv(uv['first_vector'],uv['layout'],
                      registry,bundle['UV'],scale=1e7),
                      'uv_at_muF':beta_sector_at_uv(uv['final_vector'],uv['layout'],
                      registry,bundle['UV'],scale=1e5)}
            row.update({'uv_nfev':uv['uv_nfev'],
                        'uv_scalar1_tev':start['scalar_masses_gev'][0]/1000,
                        'uv_scalar2_tev':start['scalar_masses_gev'][1]/1000,
                        'scalar_light_muF_tev':(end['scalar_masses_gev'][0]/1000 if end['scalar_masses_gev'][0] is not None else None),
                        'scalar_heavy_muF_tev':(end['scalar_masses_gev'][1]/1000 if end['scalar_masses_gev'][1] is not None else None),
                        'min_fermion_muF_tev':end['min_fermion_mass_gev']/1000,
                        'muF_scalar_eigen_mass_squared_gev2':end['scalar_eigen_mass_squared_gev2']})
            # Inverted-ordering cases do not get an EFT1 evolution. Their old
            # 1 TeV endpoint is explicitly a *formal extrapolation* only.
            check=classify(uv,None,ratio_max=ratio_max,
                           scalar_min_tev=scalar_min_tev,scalar_max_tev=scalar_max_tev)
            eft=None
            if check['ordering_at_muF']:
                eft=evaluate_eft1(model,uv,trial,registry,bundle,grid=grid)
                check=classify(uv,eft,ratio_max=ratio_max,
                               scalar_min_tev=scalar_min_tev,scalar_max_tev=scalar_max_tev)
            row.update(check)
            row['scalar_light_muS_tev']=(eft['final']['scalar_masses_gev'][0]/1000
                    if eft and eft['final']['scalar_masses_gev'][0] is not None else None)
            row['scalar_heavy_muS_tev']=(eft['final']['scalar_masses_gev'][1]/1000
                    if eft and eft['final']['scalar_masses_gev'][1] is not None else None)
            row['eft1_nfev']=eft['eft1_nfev'] if eft else 0
            row['status']='CandidateScreenPass' if check['screen_pass'] else 'CandidateScreenFail'
            if check['screen_pass']:passes.append((y,s))
        except Exception as exc:
            if original:
                raise RuntimeError(model+': baseline UV replay/sector diagnostic failed') from exc
            row.update({'status':'NumericalTrialFailed','reason':str(exc)[:400],
                        'screen_pass':False,'ordering_at_muF':False,
                        'strict_hierarchy_at_muF':False,'near_1tev_at_muS':False})
        trials.append(row)
        print(f'{model}: Y×{y:g}, mUV×{s:g}: {row["status"]}; '
              f'rhoF={row.get("ratio_muF")}; {row.get("reason")}',flush=True)
    return {'model':model,'source_fit_chi2':float(saved['chi2_final']),
            'beta_decomposition':beta,'n_trials':len(trials),
            'n_candidate_screen_passes':len(passes),'candidate_pass_coordinates':passes,
            'trials':trials}


def summary_markdown(data:dict)->str:
    lines=['# Batch 15 — UV scalar mass feasibility (B/C)','',
       'Candidate UV/EFT1 RGBeta only; this is **not a refit** or certified matching.',
       'Yukawa factors multiply all 18 saved fitted real Yukawas, and scalar mass factors multiply',
       'the original positive UV scalar masses (so m² is multiplied by the factor squared).',
       'MF and general-quartic boundary convention are held fixed.','',
       f"Screen: mSmax(muF)/MFmin(muF) <= {data['screen']['ratio_max']:.3g}, "
       f"and all scalar running masses at muS in [{data['screen']['scalar_min_tev']:.3g}, "
       f"{data['screen']['scalar_max_tev']:.3g}] TeV.",
       'A screen pass is NOT an accepted oscillation fit; chi2 is NOT recomputed.','',
       '| Model | Y factor | UV scalar mass factor | rho(muF) | Scalar masses at muS (TeV) | Screen |',
       '|---|---:|---:|---:|---|---|']
    for model,report in data['models'].items():
        for r in report['trials']:
            masses=(f"{r['scalar_light_muS_tev']:.3g}–{r['scalar_heavy_muS_tev']:.3g}"
                    if r.get('scalar_light_muS_tev') is not None and
                       r.get('scalar_heavy_muS_tev') is not None else 'not run / failed')
            rho=f"{r['ratio_muF']:.4g}" if r.get('ratio_muF') is not None else '—'
            lines.append(f"| {model} | {r['yukawa_factor']:.3g} | "
                         f"{r['uv_scalar_mass_factor']:.3g} | {rho} | {masses} | "
                         f"{r['status']} |")
    lines.extend(['','## UV mass-squared beta dependence at fixed state','','Values are d(mS_i²)/dln(mu), in GeV²; signed terms may cancel.',
       'The inclusion/exclusion varies y1,y2 and MF *only during expression evaluation*,',
       'without changing any solved trajectory. This is an algebraic diagnostic, not a',
       'unique diagrammatic decomposition.',''])
    for model,report in data['models'].items():
        for where,block in report['beta_decomposition'].items():
            lines.append(f'### {model}: {where}')
            lines.append('| Mass | Full | No Y/MF | Y only (MF=0) | MF only (Y=0) | Joint Y–MF |')
            lines.append('|---|---:|---:|---:|---:|---:|')
            for name,terms in block['terms'].items():
                lines.append('| '+name+' | '+' | '.join(f'{terms[key]:.5g}' for key in (
                 'full_d_m2_d_lnmu_gev2','no_y_no_MF','y_dependence_at_MF0',
                 'MF_dependence_at_y0','joint_y_MF_dependence'))+' |')
            lines.append('')
    lines.extend(['## Interpretation','',
       '- UV masses, heavy-fermion masses and threshold ordering are RUNNING parameters, not pole masses.',
       '- EFT1 is deliberately skipped for inverted-ordering UV endpoints. No scalar-first matching is performed.',
       '- Scanning Yukawa magnitudes changes neutrino observables; candidate points MUST be refitted with the five-observable chi2.',
       '- Baseline UV evolution is replayed with the source-fingerprint-validated fitted Yukawas.',
       '- If the prior Batch 14 diagnostic exists, its 100 TeV scalar endpoint is independently reproduced.',
       '- The additional quartic mapping and sequential general-potential C5 matching remain uncertified.',
       '- If no candidate passes, the result is limited to this finite positive-UV-mass/Yukawa grid; no no-go theorem follows.',
       '- No historical or previous general-scalar output is overwritten.',''])
    return '\n'.join(lines)


def execute(*,run_id:str,models:tuple[str,...]=CASES,
            yukawa_factors:tuple[float,...]=(0.05,0.25,0.5,1.),
            mass_factors:tuple[float,...]=(0.5,1.,2.),
            ratio_max:float=0.5,scalar_min_tev:float=0.5,
            scalar_max_tev:float=5.,grid:int=3,root:Path=ROOT,
            run_ids:dict|None=None)->dict:
    b9.parse_run_id(run_id)
    dest=root/'output/studies/general_scalar_extended/diagnostics'/run_id
    require(not dest.exists(),'output already exists: use a NEW run ID')
    require(models and len(set(models))==len(models) and set(models)<=set(CASES),
            'model list must contain unique supported B/C models')
    for label,v in [('Yukawa factors',yukawa_factors),('scalar mass factors',mass_factors)]:
        require(v and len(v)<=12 and len(set(v))==len(v) and 1. in v and
                all(math.isfinite(x) and 0<x<=10. for x in v),label+' invalid')
    require(0<ratio_max<1 and math.isfinite(ratio_max),'invalid hierarchy margin')
    require(0<scalar_min_tev<scalar_max_tev and math.isfinite(scalar_max_tev),
            'invalid scalar-mass window')
    require(3<=grid<=20,'grid must be 3–20 points')
    registry=b7.read_json(root/'RGE/general/generated/general_rgbeta_quartic_registry.json')
    results={'format':FORMAT,'status':'CandidateRGFeasibilityScreenOnly',
             'no_chi2_optimization_performed':True,'physical_matching_certified':False,
             'physical_thresholds_computed':False,'historical_outputs_overwritten':False,
             'screen':{'ratio_max':ratio_max,'scalar_min_tev':scalar_min_tev,
                       'scalar_max_tev':scalar_max_tev},
             'input_grid':{'yukawa_factors':yukawa_factors,
                           'uv_scalar_mass_factors':mass_factors},'models':{}}
    for model in models:
        saved,source,settings=b14.fit_source(root,model,
              (run_ids or b14.DEFAULT_RUNS)[model],registry)
        info=model_scan(model=model,source=source,saved=saved,settings=settings,
              registry=registry,yukawa_factors=yukawa_factors,
              mass_factors=mass_factors,ratio_max=ratio_max,
              scalar_min_tev=scalar_min_tev,scalar_max_tev=scalar_max_tev,
              grid=grid)
        info['input_fit_sha256']=sha256(settings['saved_fit_path'].read_bytes()).hexdigest()
        info['input_fitted_config_sha256']=sha256(settings['saved_config_path'].read_bytes()).hexdigest()
        info['source_run_id']=(run_ids or b14.DEFAULT_RUNS)[model]
        info['verified_input_hash_count']=settings['verified_source_files']
        # If the completed Batch 14 diagnostic is present, independently
        # reproduce its 100 TeV scalar/fermion ordering on the original point.
        reference=root/'output/studies/general_scalar_extended/diagnostics/uv_mass_hierarchy_v1/diagnostic.json'
        info['batch14_muF_replay']={'checked':False,'reason':'prior diagnostic not present'}
        if reference.is_file():
            old=b7.read_json(reference)
            prior=old.get('models',{}).get(model)
            require(isinstance(prior,dict) and prior.get('source_run_id')==info['source_run_id'],
                    model+': Batch 14 model/source-run mismatch')
            expected=float(prior['diagnostic']['general_at_muF']['m_scalar_heavy_gev'])/1000.
            actual=float(info['trials'][0]['scalar_heavy_muF_tev'])
            require(math.isclose(expected,actual,rel_tol=5e-6,abs_tol=1e-6),
                    model+': Batch 14 scalar mass at muF not reproduced')
            info['batch14_muF_replay']={'checked':True,'prior_heavy_scalar_tev':expected,
                    'replayed_heavy_scalar_tev':actual,
                    'prior_diagnostic_sha256':sha256(reference.read_bytes()).hexdigest()}
        results['models'][model]=info
    dest.mkdir(parents=True,exist_ok=False)
    with (dest/'scan.csv').open('w',newline='',encoding='utf-8') as handle:
        writer=csv.DictWriter(handle,fieldnames=CSV_FIELDS)
        writer.writeheader()
        for rep in results['models'].values():
            for row in rep['trials']:
                writer.writerow({k:row.get(k) for k in CSV_FIELDS})
    (dest/'screen.json').write_text(json.dumps(results,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    (dest/'summary.md').write_text(summary_markdown(results)+'\n',encoding='utf-8')
    print('Saved:',dest/'summary.md',flush=True)
    return results


def main()->None:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-id',required=True,help='New unique run ID')
    p.add_argument('--model',choices=CASES)
    p.add_argument('--yukawa-factors',default='0.05,0.25,0.5,1')
    p.add_argument('--uv-scalar-factors',default='0.5,1,2')
    p.add_argument('--ratio-max',type=float,default=0.5)
    p.add_argument('--scalar-min-tev',type=float,default=0.5)
    p.add_argument('--scalar-max-tev',type=float,default=5.)
    p.add_argument('--grid',type=int,default=3)
    p.add_argument('--b-run-id',default=b14.DEFAULT_RUNS['B_-1'])
    p.add_argument('--c-run-id',default=b14.DEFAULT_RUNS['C_-1'])
    a=p.parse_args()
    execute(run_id=a.run_id,models=(a.model,) if a.model else CASES,
            yukawa_factors=parse_factors(a.yukawa_factors,'Yukawa factors'),
            mass_factors=parse_factors(a.uv_scalar_factors,'scalar mass factors'),
            ratio_max=a.ratio_max,scalar_min_tev=a.scalar_min_tev,
            scalar_max_tev=a.scalar_max_tev,grid=a.grid,
            run_ids={'B_-1':a.b_run_id,'C_-1':a.c_run_id})


if __name__=='__main__':
    main()
