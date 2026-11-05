"""Uniform-0.01 T3 UV/EFT1/SMEFT running ablation for thesis Figures 5.2/5.3.

Uses the 16 *newly refitted* reference Yukawa configurations (no Yukawa refits),
all allowed BSM quartics except lambdaT3 = +0.01 at muUV, and zero
conjugate Yukawas at muUV. Each stage-off comparison freezes the *full
renormalisable state* over that stage while retaining the current conditional
one-scale full-general C5 kernel at the scalar endpoint. SMEFT off freezes
SM+Weinberg parameters between 1 TeV and 100 GeV.

THIS IS NOT the full sequential fermion -> scalar -> SMEFT Wilson calculation;
no intermediate-EFT LLSS Wilson transport is claimed. No old running
percentages, matrices or plots are used. The existing heatmap plotter owns
all plotting. Source symbolic exports are pinned and verified via the
previously completed fixed-Yukawa table driver.

Run from local repository root; does not edit GitHub or existing physics code.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import replace
import json
import math
from pathlib import Path
import re
import subprocess
import sys
import traceback
from typing import Any

import numpy as np

from scripts import RunT3QualitativeScalarRefits as base
from scripts import RunT3UniformScalarBenchmark as uniform
from scripts import RunT3UniformScalarCampaign as campaign
from scripts import BuildT3Uniform001FixedYukawaTables as fixed_tables

ROOT=Path(__file__).resolve().parents[1]
REFERENCE_ROOT=ROOT/'output/studies/thesis_qualitative/uniform_scalar_001/campaigns/uniform001_all16_reference_v1'
OUTPUT_ROOT=ROOT/'output/studies/thesis_qualitative/uniform_scalar_001/stage_running'
REPORT_ROOT=ROOT/'Reports/output/thesis_general_5_2/uniform001_stage_running'
SCHEMA='t3-uniform001-3stage-conditional-one-scale-running-v1'
MODELS=tuple(base.MODELS)
STAGES=('UV','EFT1','SMEFT')
CASES={'UV':'011','EFT1':'101','SMEFT':'110'}
OBS=('delta_m21_sq_ev2','delta_m3l_sq_ev2','sin2_theta12','sin2_theta13','sin2_theta23')
MU_UV,MU_F,MU_S,MU_LOW=1e7,1e5,1e3,1e2


def select_models(value: str)->tuple[str,...]:
    if value=='all':return MODELS
    if value=='pilot':return ('A_-4','B_-1')
    out=tuple(v.strip() for v in value.split(',') if v.strip())
    if not out or len(out)!=len(set(out)) or set(out)-set(MODELS):
        raise ValueError('Unknown, repeated, or empty model names')
    return out


def tag(model: str)->str:
    return base.tag(model)


def save(path:Path, value:dict)->None:
    base.json_save(path,value)


def read(path:Path)->dict:
    return base.json_read(path)


def signature(ref_result:Path, ref_config:Path, paths:dict[str,Path])->dict:
    return {'reference_result_sha256':base.sha(ref_result),
            'reference_config_sha256':base.sha(ref_config),
            'accepted_source_sha256':{k:base.sha(v) for k,v in sorted(paths.items())}}


def frozen_state(state:Any,mu:float)->Any:
    """Freeze all couplings while updating the bookkeeping scale."""
    return replace(state,mu_gev=mu).validated()


def require_c5(c5:Any)->np.ndarray:
    a=np.asarray(c5,dtype=complex)
    if a.shape!=(3,3) or not np.isfinite(a).all():
        raise ValueError('Invalid 3-by-3 conditional matched C5')
    if np.linalg.norm(a-a.T)>1e-8*max(np.linalg.norm(a),1e-35):
        raise ValueError('Non-symmetric conditional matched C5')
    return (a+a.T)/2.


def predict_at_low(c5:np.ndarray, eft_state:Any, target:Any,
                   *,smeft_on:bool,vev:float)->dict:
    """Identical SM+Weinberg evolution/fit to the reference, with an off switch."""
    from scripts.RunT3GeneralScalarBCExploratoryFit import low_energy_fit
    from Numerical.running.ScalarThresholdBoundary import (
        build_final_sm_boundary,build_sm_weinberg_initial_conditions)
    from Numerical.fitting.ParameterScan import (
        neutrino_mass_from_c5,charged_lepton_mass_basis_matrix)
    from Numerical.fitting.OscillationFit import evaluate_mass_matrix_fit

    if smeft_on:
        fit,low=low_energy_fit(c5,eft_state,target,vev=vev)
        K,ye=low.K,low.ye
    else:
        boundary=build_final_sm_boundary(eft_state,matching_scale_gev=MU_S)
        init=build_sm_weinberg_initial_conditions(boundary,c5)
        K,ye=init.K,init.ye
        mass=neutrino_mass_from_c5(K,vev_gev=vev)
        flavour=charged_lepton_mass_basis_matrix(mass,ye)
        fit=evaluate_mass_matrix_fit(flavour,target)
    mass=neutrino_mass_from_c5(K,vev_gev=vev)
    flavour=charged_lepton_mass_basis_matrix(mass,ye)
    masses=sorted((np.linalg.svd(flavour,compute_uv=False)*1e9).tolist())
    if len(masses)!=3 or not np.isfinite(masses).all():
        raise ValueError('Nonfinite mass eigenvalues')
    p=base.fit_prediction(fit)
    for key in OBS:
        if not math.isfinite(float(p[key])):
            raise ValueError(f'Nonfinite observable {key}')
    return {'prediction':p,'chi2':float(fit.chi2),'neutrino_masses_ev':masses,
            'C5_low_gev_inv':{'real':np.asarray(K).real.tolist(),
                               'imag':np.asarray(K).imag.tolist()}}


def _nonmixed(model:str, evaluator:Any, x:np.ndarray,target:Any,vev:float)->dict:
    from scripts import RunT3GeneralAll16FromReference as general
    from Numerical.orchestration.GeneralRunningBenchmark import initial_state
    from Numerical.running.UVRunner import run_uv_segment
    from Numerical.running.FermionThresholdBoundary import build_intermediate_scalar_boundary
    from Numerical.running.IntermediateScalarRunner import run_intermediate_scalar_segment
    from Numerical.running.IntermediateScalarStateVector import pack_intermediate_scalar_state

    y=base.parameter_dict(x)
    _,initial=initial_state(evaluator.uv,mu_uv=MU_UV,mu_f=MU_F,mu_s=MU_S,
                       overrides=general.override_yukawas(evaluator.fixed,y))
    out={}
    for uv_on in (True,False):
        if uv_on:
            sol=run_uv_segment(initial,evaluator.uv,MU_F,rtol=1e-7,atol=1e-10)
            at_f=sol.final_state
        else:at_f=frozen_state(initial,MU_F)
        at_boundary=build_intermediate_scalar_boundary(at_f)
        _,layout=pack_intermediate_scalar_state(at_boundary)
        if set(layout.names)!=set(evaluator.eft.get('report_betas',{})):
            raise ValueError(f'{model}: EFT1 beta/state inventory mismatch')
        for eft_on in ((True,False) if uv_on else (True,)):
            if eft_on:
                sol=run_intermediate_scalar_segment(at_boundary,evaluator.eft,MU_S,
                                                     rtol=1e-7,atol=1e-10)
                at_s=sol.final_state
            else:at_s=frozen_state(at_boundary,MU_S)
            if min(at_s.mS1Sq,at_s.mS2Sq)<=0:
                raise ValueError(f'{model}: nonpositive scalar mass-squared')
            c5=require_c5(evaluator.kernel.evaluate(at_f,at_s))
            sm_modes=(True,False) if (uv_on and eft_on) else (True,)
            for sm_on in sm_modes:
                key=f'{int(uv_on)}{int(eft_on)}{int(sm_on)}'
                out[key]=predict_at_low(c5,at_s,target,smeft_on=sm_on,vev=vev)
    return out


def _mixed(model:str,evaluator:Any,x:np.ndarray,target:Any,vev:float)->dict:
    from scripts import RunT3GeneralScalarBCExploratoryFit as bc
    from Numerical.experimental.GeneralT3ScalarRGState import (
        unpack_general_scalar_state,project_general_after_fermion,
        pack_general_scalar_state)

    trial=bc.trial_record(evaluator.source,np.asarray(x,dtype=float))
    initial,ly0,extra,mu2=uniform.bc_prepare_initial_absolute(
        trial,evaluator.registry,evaluator.mu2_fraction)
    new,reused=trial['new'],trial['reused']
    out={}
    for uv_on in (True,False):
        if uv_on:
            sol=bc.b7.integrate_segment(
                initial,lambda t,y:bc.b7.rhs_at(t,y,layout=ly0,
                  registry=evaluator.registry,betas=evaluator.bundle['UV'],stage='UV'),
                MU_UV,MU_F,n_points=3)
            uv_vec=sol.y[:,-1].copy()
        else:uv_vec=initial.copy()
        at_f=unpack_general_scalar_state(uv_vec[:ly0.size],ly0,
                  evaluator.registry,uv=True,mu_gev=MU_F)
        heavy=at_f.baseline
        eft_initial=project_general_after_fermion(at_f,evaluator.registry)
        eft_vec,lye=pack_general_scalar_state(eft_initial,evaluator.registry)
        eft_vec=np.concatenate((eft_vec,uv_vec[ly0.size:]))
        betas=bc.b7.check_payload(
            evaluator.bundle['payloads']['EFT1'],model=model,stage='EFT1',
            rep=lye.representation,base_names=tuple(lye.names[:-len(new)]),
            new_names=new,reused_names=reused)
        for eft_on in ((True,False) if uv_on else (True,)):
            if eft_on:
                sol=bc.b7.integrate_segment(
                    eft_vec,lambda t,y:bc.b7.rhs_at(t,y,layout=lye,
                      registry=evaluator.registry,betas=betas,stage='EFT1'),
                    MU_F,MU_S,n_points=3)
                end_vec=sol.y[:,-1].copy()
            else:end_vec=eft_vec.copy()
            snap=bc.b7.snapshot(end_vec,lye,evaluator.registry,
                                stage='EFT1',mu_gev=MU_S,new=new,reused=reused)
            at_s=unpack_general_scalar_state(end_vec[:lye.size],lye,
                            evaluator.registry,uv=False,mu_gev=MU_S).baseline
            heavy_info=bc.heavy_flavour_at_muF(heavy)
            c5=require_c5(bc.c5_at_endpoint(model,snap,heavy_info,evaluator.coefficient))
            sm_modes=(True,False) if (uv_on and eft_on) else (True,)
            for sm_on in sm_modes:
                key=f'{int(uv_on)}{int(eft_on)}{int(sm_on)}'
                out[key]=predict_at_low(c5,at_s,target,smeft_on=sm_on,vev=vev)
    return out


def verify_baseline(model:str,cases:dict, reference:dict,
                    *, tolerance:float=1e-5)->dict:
    """Fail closed if any newly computed full-running observable or C5 differs."""
    if set(cases)!=set(('111','011','101','110')):
        raise ValueError(f'{model}: expected four complete running cases')
    current=cases['111']
    for key in OBS:
        actual=float(current['prediction'][key]);saved=float(reference['prediction'][key])
        if not math.isfinite(actual) or abs(actual-saved)>max(1e-10,tolerance*abs(saved)):
            raise ValueError(f'{model}: full-running baseline differs in {key}: '
                             f'{actual:g} versus saved {saved:g}')
    from scripts.BuildT3Uniform001FixedYukawaTables import (
        c5_from_json,low_c5)
    got=c5_from_json(current['C5_low_gev_inv'])
    expect=c5_from_json(low_c5(reference['numerical_details']))
    norm=float(np.linalg.norm(expect))
    if norm<=0:raise ValueError('Saved baseline C5 is zero')
    rel=float(np.linalg.norm(got-expect)/norm)
    if not math.isfinite(rel) or rel>tolerance:
        raise ValueError(f'{model}: saved C5 baseline not reproduced (rel {rel:.3g})')
    return {'C5_low_relative_difference':rel,
            'max_observable_scaled_residual':max(
                abs(float(current['prediction'][k])-float(reference['prediction'][k]))/
                max(abs(float(reference['prediction'][k])),1e-15) for k in OBS)}


def percent_change(full:float, frozen:float)->float:
    if not (math.isfinite(full) and math.isfinite(frozen)) or frozen==0:
        raise ValueError('Cannot calculate relative running change')
    return 100.*(full-frozen)/frozen


def stage_rows(model:str,cases:dict)->list[dict]:
    if set(cases)!=set(('111','011','101','110')):
        raise ValueError('Only complete three-stage running set can produce heatmap rows')
    out=[];full=cases['111']['prediction']
    for stage in STAGES:
        case=CASES[stage];frozen=cases[case]['prediction']
        for observable in OBS:
            a=float(full[observable]);b=float(frozen[observable])
            out.append({'model':model,'switch':stage,'on_case':'111',
                        'off_case':case,'observable':observable,
                        'full_value':a,'frozen_value':b,
                        'absolute_change':a-b,
                        'relative_change_percent':percent_change(a,b),
                        'comparison_kind':'independently_evolved_running_stage',
                        'benchmark':'uniform-0.01-zero-conjugate',
                        'matching_scope':'conditional-one-scale-not-sequential'})
    return out


def _run_model(model:str, selection:dict, target:Any, args:Any)->dict:
    from scripts import RunT3GeneralZeeSensitivity as zee
    ns=argparse.Namespace(fit_run_id='merged',canonical_run_id=args.canonical_run_id)
    original,oldfit,paths=zee.sources(model,ns,selection[model])
    ref_dir=args.reference_root/'reference'/tag(model)
    reference,cfg=campaign.validated_reference(
        ref_dir,model=model,source_paths=paths,manifest_path=zee.MERGED_FILE,
        original_config=original,chi2_cutoff=args.chi2_cutoff)
    x=base.fitted_vector(cfg);lam=base.real_lambda(cfg)
    vev=float(cfg.get('numerical',{}).get('vev_gev',246.22))
    if model in base.MIXED:
        factory=uniform.ModelFactory(model,cfg,oldfit,paths,target,args.canonical_run_id)
        evaluator,_,_,boundary=factory.build(lam)
        cases=_mixed(model,evaluator,x,target,vev)
        mode='saved-general-BC-conditional-mass-basis'
    else:
        factory=fixed_tables.PinnedCanonicalNonMixedFactory(
            model,original,paths,target,args.canonical_run_id)
        fixed_tables.verify_pinned_baseline(factory,reference,cfg)
        evaluator,_,_,boundary=factory.build(lam)
        cases=_nonmixed(model,evaluator,x,target,vev)
        mode='opt-in-frozen-canonical-source-with-reproduction-check'
    checked=verify_baseline(model,cases,reference,tolerance=args.baseline_tolerance)
    return {'schema':SCHEMA,'status':'Computed','model':model,
            'reference_run_id':args.reference_root.name,
            'lambdaT3':lam,'extra_BSM_quartics_uv':uniform.QUARTIC,
            'conjugate_yukawas_uv_zero':True,'scheme':mode,
            'matching_scope':'conditional-one-scale-not-certified-sequential',
            'reference_chi2':float(reference['chi2']),
            'source_fingerprints':signature(ref_dir/'result.json',
                ref_dir/'fitted_config.json',paths),
            'baseline_reproduction':checked,'boundary':boundary,
            'cases':cases}


def _raw_rows(output_root:Path,models:tuple[str,...])->list[dict]:
    values=[]
    for model in models:
        path=output_root/'models'/tag(model)/'result.json'
        d=read(path)
        if d.get('schema')!=SCHEMA or d.get('model')!=model or d.get('status')!='Computed':
            raise ValueError(f'{model}: result lacks matching schema/model/status')
        if d.get('extra_BSM_quartics_uv')!=0.01 or d.get('conjugate_yukawas_uv_zero') is not True:
            raise ValueError(f'{model}: not the requested 0.01/zero-conjugate result')
        if not isinstance(d.get('baseline_reproduction'),dict):
            raise ValueError(f'{model}: full-running baseline not verified')
        values.extend(stage_rows(model,d['cases']))
    return values


def write_heatmap_csv(run_root:Path,models:tuple[str,...])->Path:
    if set(models)!=set(MODELS):
        raise ValueError('Need all 16 verified new-benchmark models before plotting thesis heatmaps')
    rows=_raw_rows(run_root,MODELS)
    if len(rows)!=len(MODELS)*len(STAGES)*len(OBS):
        raise ValueError('Incomplete stage comparison rows')
    dest=run_root/'frozen_stage_observable_changes.csv'
    dest.parent.mkdir(parents=True,exist_ok=True)
    cols=tuple(rows[0])
    with dest.open('w',newline='',encoding='utf-8') as fh:
        writer=csv.DictWriter(fh,fieldnames=cols);writer.writeheader();writer.writerows(rows)
    return dest


def make_plots(run_root:Path,report_root:Path)->list[Path]:
    source=write_heatmap_csv(run_root,MODELS)
    report_root.mkdir(parents=True,exist_ok=True)
    plots=report_root/'plots';plots.mkdir(parents=True,exist_ok=True)
    cmd=[sys.executable,'-X','utf8','-m','scripts.PlotT3OscillationHeatmaps',
         '--input',str(source),'--plot-dir',str(plots),'--csv-dir',str(report_root/'csv')]
    subprocess.run(cmd,cwd=ROOT,check=True)
    provenance={'schema':SCHEMA,'reference_run_id':REFERENCE_ROOT.name,
       'data_run_id':run_root.name,'new_benchmark_extra_scalar_quartic_uv':0.01,
       'conjugate_yukawas_uv_zero':True,
       'matching':'conditional-one-scale-C5-not-sequential-EFT',
       'stage_off_definition':'freeze whole renormalisable state of indicated stage',
       'SMEFT_off_definition':'freeze SM and C5 from 1 TeV to 100 GeV',
       'frozen_csv_sha256':base.sha(source),'model_count':16,
       'model_reference_check':'saved new uniform-0.01 fitted baseline reproduced',
       'old_2026_running_csv_used':False,'original_plotter_unmodified':True}
    save(report_root/'provenance.json',provenance)
    return [plots/'oscillation_mixing_heatmap_percent.png',
            plots/'oscillation_masssq_heatmap_percent.png']


def main(argv=None)->int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--stage',choices=('preflight','compute','plot'),default='preflight')
    p.add_argument('--models',default='pilot',help='pilot, all, or comma-separated model IDs')
    p.add_argument('--reference-root',type=Path,default=REFERENCE_ROOT)
    p.add_argument('--canonical-run-id',default='full16_general_v1')
    p.add_argument('--run-id',default='uniform001_running_v1')
    p.add_argument('--chi2-cutoff',type=float,default=5.86)
    p.add_argument('--baseline-tolerance',type=float,default=1e-5)
    p.add_argument('--resume',action='store_true')
    args=p.parse_args(argv)
    models=select_models(args.models)
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,63}',args.run_id):
        p.error('run-id must use only letters, digits, _ and -')
    if args.baseline_tolerance<=0 or args.baseline_tolerance>1e-3:
        p.error('--baseline-tolerance must be in (0,1e-3]')
    run_root=OUTPUT_ROOT/args.run_id
    report_root=REPORT_ROOT/args.run_id
    if args.stage=='plot':
        files=make_plots(run_root,report_root)
        print('Updated thesis heatmaps:',*(str(x) for x in files),sep='\n  ',flush=True)
        return 0
    from scripts import RunT3GeneralZeeSensitivity as zee
    selection=zee.fit_selections(list(models),'merged')
    ns=argparse.Namespace(fit_run_id='merged',canonical_run_id=args.canonical_run_id)
    for model in models:
        original,oldfit,paths=zee.sources(model,ns,selection[model])
        ref=args.reference_root/'reference'/tag(model)
        campaign.validated_reference(ref,model=model,source_paths=paths,
          manifest_path=zee.MERGED_FILE,original_config=original,
          chi2_cutoff=args.chi2_cutoff)
        print(f'{model}: accepted uniform-0.01 new reference verified',flush=True)
    if args.stage=='preflight':
        print(f'Preflight complete: {len(models)}/{len(models)} validated, no physics run',flush=True)
        return 0
    if run_root.exists() and not args.resume:
        raise FileExistsError(f'{run_root}: existing run, use --resume')
    target=zee.target_object()
    errors=[]
    for model in models:
        out=run_root/'models'/tag(model)/'result.json'
        if args.resume and out.is_file():
            previous=read(out)
            if previous.get('schema')==SCHEMA and previous.get('status')=='Computed':
                _,_,paths=zee.sources(model,ns,selection[model])
                ref=args.reference_root/'reference'/tag(model)
                if previous.get('source_fingerprints')==signature(
                          ref/'result.json',ref/'fitted_config.json',paths):
                    print(f'{model}: already saved (verified input fingerprints)',flush=True)
                    continue
        print(f'{model}: UV/EFT1/SMEFT stage-off evaluations',flush=True)
        try:
            result=_run_model(model,selection,target,args)
            save(out,result)
            print(f'  baseline reproduced; stages UV, EFT1, SMEFT complete',flush=True)
        except Exception as e:
            errors.append(model)
            save(out,{'schema':SCHEMA,'status':'Failed','model':model,
               'reason':f'{type(e).__name__}: {e}',
               'traceback':traceback.format_exc(limit=7)[-3000:]})
            print(f'  FAILED: {type(e).__name__}: {e}',flush=True)
    print(f'Finished {len(models)-len(errors)}/{len(models)} selected models',flush=True)
    if errors:
        print('Blocked: '+', '.join(errors)+'; fix and repeat with --resume',flush=True)
        return 1
    if set(models)==set(MODELS):
        source=write_heatmap_csv(run_root,MODELS)
        print(f'Complete 16/16: {source}',flush=True)
        print(f'Run: python -X utf8 -m scripts.RunT3Uniform001StageRunning '
              f'--stage plot --run-id {args.run_id}',flush=True)
    return 0

if __name__=='__main__':
    raise SystemExit(main())
