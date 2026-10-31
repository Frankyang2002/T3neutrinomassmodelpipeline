"""Checkpointed, additive T3 integration campaign using *existing source physics*.

Stages:
  UV: all-16 full matrix RGBeta candidate UV running with 6 optional Yukawas
  BC: existing B/C full-mass UV->EFT1 runner plus scalar action/eigenbasis
  matching: archived restricted B/C Matchete v23 finite triangle reference
  legacy: optional existing conditional B/C Y/X sequential numerical pipeline
  validation: provenance, continuity, convergence, and certification gates

No existing production files are changed. No physical C5 certificate is created.
The matched B/C reference and coupled mass RGE study do not yet share a
certified physical action dictionary, and are NEVER combined as if they did.

  python -u -X utf8 -m scripts.RunT3FullStackOvernight --preflight
  python -u -X utf8 -m scripts.RunT3FullStackOvernight --all --resume
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from hashlib import sha256
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import time
import traceback
from typing import Any, Mapping

import numpy as np
from scipy.integrate import solve_ivp

from Numerical.experimental.T3GeneralMatrixBeta import (
    EXTRA_MODELS, ExtendedUVEquationSystem, trial_extra_seed, finite_matrix,
)
from Numerical.experimental.T3BCThresholdInteractionBridge import bridge_trajectory

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT/'output/studies/t3_full_stack_overnight_v1'
REGISTRY = ROOT/'RGE/general/generated/general_rgbeta_quartic_registry.json'
STAGE29 = ROOT/'output/general_t3_stage22/stage29_convention_uv_handoff.json'
UV_EXPORT = ROOT/'Reports/output/workstream_a/integrated_general_yukawa_uv/native'
UV_EXPORT_SUMMARY = UV_EXPORT.parent/'summary.json'
RGBETA_MODEL = ROOT/'RGE/running/rgbeta/T3RGBetaModel.wl'
BC_MODELS = ('B_-1', 'C_-1')
CONVENTION = '16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g'
FORMAT = 't3-additive-full-stack-overnight-candidate-v1'


def require(ok: bool, msg: str) -> None:
    if not ok:
        raise ValueError('Overnight T3: '+msg)


def read_json(path: Path) -> dict:
    require(path.is_file(), 'missing source '+str(path))
    obj = json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(obj, dict), 'expected JSON object '+str(path))
    return obj


def digest(path: Path) -> str:
    require(path.is_file(), 'cannot hash missing source '+str(path))
    return sha256(path.read_bytes()).hexdigest()


def write_atomic(path: Path, value: Mapping[str, Any], *, overwrite: bool = False) -> None:
    require(overwrite or not path.exists(), 'refusing overwrite '+str(path))
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name+'.tmp')
    require(not temp.exists(), 'prior incomplete temporary output '+str(temp))
    try:
        with temp.open('w', encoding='utf-8') as f:
            json.dump(value, f, indent=2, allow_nan=False)
            f.write('\n')
        os.replace(temp, path)
    finally:
        if temp.exists():temp.unlink()


def model_specifications(registry: Mapping[str, Any]) -> list[str]:
    models = registry.get('models', {})
    require(len(models) == 16 and set(EXTRA_MODELS) <= set(models),
            'expected 16 T3 models and six extra Yukawa models')
    total = sum(int(x['total_additional_quartics']) for x in models.values())
    extra_quadratics = sum(int(x['quadratics_not_implemented']) for x in models.values())
    require(total == 50 and extra_quadratics == 8,
            'general-scalar inventory changed; new source review required')
    return sorted(models)


def safe_model(model: str) -> str:
    require(bool(re.fullmatch(r'[A-E]_[+-]\d+', model)), 'invalid model name')
    return model.replace('-', 'm').replace('+', 'p')


def uv_dependencies(model: str, registry: dict) -> dict[str, str]:
    from scripts import CoupledUVRunning as stage31
    row = registry['models'][model]
    files = [REGISTRY, STAGE29, UV_EXPORT_SUMMARY, RGBETA_MODEL,
             UV_EXPORT/(safe_model(model)+'.json'),
             stage31.config_path(row),
             ROOT/'Numerical/core/RGBetaEvaluator.py',
             ROOT/'Numerical/core/BetaVector.py',
             ROOT/'Numerical/core/StateVector.py',
             ROOT/'scripts/CoupledUVRunning.py',
             ROOT/'scripts/CoupledA0UVRunning.py',
             ROOT/'Numerical/experimental/T3GeneralMatrixBeta.py']
    return {str(p.relative_to(ROOT)).replace('\\','/'): digest(p) for p in files}


def load_uv_inputs(model: str, registry: dict):
    """Rebuild the existing Stage31 *candidate* fitted boundary without mutation."""
    from scripts import CoupledUVRunning as stage31
    from scripts import CoupledA0UVRunning as stage30
    meta = registry['models'][model]
    require(meta['model'] == model, 'model registry entry mismatch')
    new = tuple(sorted(x['name'] for x in meta['new_rgbeta_symbols']))
    reused = tuple(x['existing_name'] for x in meta['preexisting_rgbeta_symbols'])
    require(len(new)+len(reused) == int(meta['total_additional_quartics']),
            'quartic registry incomplete')
    handoff = read_json(STAGE29)
    stage30.verify_stage29_registry_provenance(handoff, REGISTRY)
    seed = stage31.stage29_values(handoff, meta, new, reused)
    config = read_json(stage31.config_path(meta))
    vec, layout = stage31.prepare_state(config, meta, reused, seed)
    return vec, layout, {key: seed[key] for key in new}


def source_equations(model: str, registry: dict, layout: Any):
    summary = read_json(UV_EXPORT_SUMMARY)
    require(summary.get('all_requested_models_succeeded') is True and
            summary.get('local_model_sha256') == digest(RGBETA_MODEL) and
            len(summary.get('rows', [])) == 16,
            'native UV beta exports are stale relative to current RGBeta model; regenerate Workstream A exports')
    row = registry['models'][model]
    path = UV_EXPORT/(safe_model(model)+'.json')
    report = read_json(path)
    meta = report.get('metadata', {})
    expect = {'dS1': row['dimensions'][0], 'dS2': row['dimensions'][1],
              'dF': row['dimensions'][2], 'alpha': row['alpha'],
              'SharedScalar': False, 'MajoranaFermion':
                  bool(row['alpha']==-1 and row['dimensions'][2]%2==1)}
    # For Dirac representations with alpha != -1, the source must say non-Majorana.
    require(report.get('status') == 'Success' and isinstance(meta, dict),
            'native UV export status invalid for '+model)
    require(all(meta.get(k)==v for k,v in expect.items()),
            'native UV model representation/Majorana mismatch for '+model)
    return ExtendedUVEquationSystem(
        model=model, base_layout=layout,
        extra_quartics=tuple(sorted(x['name'] for x in row['new_rgbeta_symbols'])),
        extra_yukawa=EXTRA_MODELS.get(model),
        expressions=report['report_betas'], source_metadata=meta)


def solve_uv(model: str, *, registry: dict, seed_fraction: float,
             points: int, rtol: float, atol: float, convergence: bool = False) -> dict:
    from Numerical.core.StateVector import unpack_uv_state
    from Numerical.experimental.GeneralT3RunningFermionThresholds import (
        locate_fermion_crossings, mass_matrix_from_vector)
    require(3 <= points <= 1001 and 0 < rtol < .1 and 0 < atol < .1,
            'bad integration sampling/tolerances')
    baseline, layout, quartics = load_uv_inputs(model, registry)
    equations = source_equations(model, registry, layout)
    init_state = unpack_uv_state(baseline,layout,mu_gev=1e7)
    e0 = (trial_extra_seed(init_state.y1, init_state.y2, seed_fraction)
          if equations.extra_yukawa else None)
    y0 = equations.initial_vector(baseline, quartics, e0)
    t0,t1=math.log(1e7),math.log(1e5)
    first_beta = equations.derivative(t0,y0)
    tgrid=np.linspace(t0,t1,points)
    def integrate(tolerance):
        s=solve_ivp(equations.derivative,(t0,t1),y0,
                    method='DOP853',rtol=tolerance,atol=atol,
                    max_step=.35,t_eval=tgrid,dense_output=True)
        require(bool(s.success) and s.y.shape==(equations.size,points) and
                np.isfinite(s.y).all(),model+': coupled UV solve failed: '+s.message)
        return s
    sol=integrate(rtol)
    final=sol.y[:,-1]
    out = {
        'format': FORMAT,'status': 'CandidateFullMatrixRGBetaUVIntegrated',
        'model': model, 'extra_yukawa':equations.extra_yukawa,
        'extra_yukawa_seed':'independent artificial RGBeta matrix, nonphysical' if e0 is not None else None,
        'extra_yukawa_seed_fraction':seed_fraction if e0 is not None else 0.,
        'uv_start_gev':1e7,'uv_stop_gev':1e5,
        'one_loop_convention':CONVENTION,'loop_factor_applied_once':True,
        'ode_method':'DOP853','rtol':rtol,'atol':atol,'nfev':int(sol.nfev),
        'beta_parameter_names':list(equations.all_names),
        'beta_coordinate_count':equations.size,
        'initial_beta_euclidean_norm':float(np.linalg.norm(first_beta)),
        'initial':equations.summary(y0),'final':equations.summary(final),
        'trajectory':[{'mu_gev':float(math.exp(tgrid[i])),
                       'extra_yukawa_frobenius_norm':
                           equations.summary(sol.y[:,i])['extra_yukawa_frobenius_norm'],
                       'new_quartic_max_magnitude':max((abs(v) for v in equations.split(sol.y[:,i])[1].values()),default=0.)}
                       for i in range(points)],
        'physical_matchete_rgbeta_dictionary_certified':False,
        'general_scalar_holomorphic_quadratic_beta_in_this_UV_export':False,
        'fermion_removed_sequentially':False,
        'UV_to_EFT1_matching_performed':False,
        'physical_neutrino_masses_computed':False,
    }
    base_end,_,extrafinal=equations.split(final)
    start=unpack_uv_state(baseline,layout,mu_gev=1e7)
    end=unpack_uv_state(base_end,layout,mu_gev=1e5)
    mf0=np.asarray(start.MF,complex);mf1=np.asarray(end.MF,complex)
    out['MF_initial_singular_masses_gev']=np.sort(np.linalg.svd(mf0,compute_uv=False)).tolist()
    out['MF_final_singular_masses_gev']=np.sort(np.linalg.svd(mf1,compute_uv=False)).tolist()
    out['y1_frobenius_initial_final']=[float(np.linalg.norm(start.y1)),float(np.linalg.norm(end.y1))]
    out['y2_frobenius_initial_final']=[float(np.linalg.norm(start.y2)),float(np.linalg.norm(end.y2))]
    # This is a full-UV diagnostic only: roots must be recomputed after physical decoupling.
    majorana=bool(registry['models'][model]['alpha']==-1 and registry['models'][model]['dimensions'][2]%2==1)
    try:
        crosses=locate_fermion_crossings(
            lambda t: mass_matrix_from_vector(equations.split(sol.sol(t))[0], layout),
            1e7,1e5,majorana=majorana)
        out['candidate_full_UV_fermion_crossings']=crosses
    except (ValueError, RuntimeError) as exc:
        out['candidate_full_UV_fermion_crossings']={'status':'DiagnosticBlocked','reason':str(exc)}
    if convergence:
        improved=integrate(rtol/5)
        final_ref=improved.y[:,-1]
        relative=float(np.linalg.norm(final-final_ref)/max(np.linalg.norm(final_ref),1e-25))
        out['convergence']={'refined_rtol':rtol/5,'relative_full_vector_endpoint_difference':relative,
                            'refined_nfev':int(improved.nfev),
                            'passed':bool(relative<max(10*rtol,1e-8))}
        require(out['convergence']['passed'],model+': UV convergence tolerance failed')
    return out


def stage_bc(model:str, *, output:Path, points:int, mu2_fraction:float,
             mu2_phase:float, higgs_mass_squared:float, resume:bool):
    """Use the existing physically limited coupled B/C five-mass runner."""
    from scripts.RunT3BCCoupledUVToEFT1 import run_model
    folder=output/'bc_full_mass'
    tag=safe_model(model)
    target=folder/f'{tag}_uv_to_eft1_coupled_five_mass_v1.json'
    if target.exists():
        require(resume, 'already exists, use --resume or a new --output-id')
        report=read_json(target)
        require(report.get('model')==model and report.get('solver_success') is True and
                math.isclose(report.get('initial_mu2_fraction_at_uv',-1),mu2_fraction,rel_tol=0,abs_tol=1e-13) and
                math.isclose(report.get('initial_mu2_phase_at_uv',1e99),mu2_phase,rel_tol=0,abs_tol=1e-13) and
                math.isclose(report.get('initial_higgs_mass_squared_gev2_at_uv',1e99),higgs_mass_squared,rel_tol=0,abs_tol=1e-13) and
                len(report.get('uv_trajectory',[]))==points,
                'resume B/C input/length mismatch')
        prior_sources=report.get('input_sha256')
        require(isinstance(prior_sources,dict) and len(prior_sources)>=6,
                'B/C resume lacks source provenance')
        for source, oldhash in prior_sources.items():
            p=Path(source)
            if not p.is_absolute():p=ROOT/p
            require(p.is_file() and digest(p)==oldhash,
                    'B/C five-mass RGE source/benchmark changed: '+str(source))
    else:
        report=run_model(model,output_dir=folder,points=points,
                         mu2_fraction=mu2_fraction,mu2_phase=mu2_phase,
                         higgs_mass_squared=higgs_mass_squared)
    bridge=bridge_trajectory(model,report)
    bridge['source_five_mass_report_sha256']=digest(target)
    return report,bridge


def run_subprocess(root:Path,args:list[str],log_path:Path)->dict:
    """Run an existing reference calculation in an isolated process."""
    require(all(isinstance(a,str) for a in args), 'invalid subprocess argv')
    log_path.parent.mkdir(parents=True,exist_ok=True)
    start=time.monotonic()
    with log_path.open('w', encoding='utf-8') as f:
        f.write('COMMAND: '+ ' '.join(args)+'\n\n');f.flush()
        rc=subprocess.run(args,cwd=root,stdout=f,stderr=subprocess.STDOUT,check=False).returncode
    return {'argv':args,'exit_code':int(rc),'elapsed_s':float(time.monotonic()-start),
            'log':str(log_path.relative_to(root)) if log_path.is_relative_to(root) else str(log_path)}


def reference_matching(output:Path)->dict:
    """Recompute archived B/C restricted matching; never pass B/C RGE masses here."""
    script=ROOT/'scripts/FinalizeT3BCLeadingBenchmark.experimental.py'
    require(script.is_file(), 'missing existing C v23 matching reference')
    path=output/'matching/bc_restricted_v23.json'
    if not path.exists():
        log=run_subprocess(ROOT,[sys.executable,'-u','-X','utf8',str(script),
            '--root',str(ROOT),'--output',str(path),'--low-scale','100'],
            output/'logs/matching_v23.log')
        require(log['exit_code']==0, 'reference matching subprocess failed: '+str(log))
    result=read_json(path)
    # Persisted v23 reports are only reusable while all four source-audit
    # hashes agree; a pass flag from a previous run is not enough.
    source_hashes=result.get('source_audits_sha256')
    require(isinstance(source_hashes,dict) and len(source_hashes)==4,
            'matching source provenance missing')
    for source,recorded in source_hashes.items():
        relative=Path(str(source).replace('\\','/'))
        require(not relative.is_absolute() and '..' not in relative.parts and
                digest(ROOT/relative)==recorded,
                'matching source audit changed: '+str(source))
    require(result.get('source_checks_passed') is True and
            result.get('direct_Matchete_all_heavy_equals_sequential_full_mass') is True and
            result.get('physical_C5_certified') is False,
            'restricted reference matching source checks failed')
    return {'status':result.get('status'), 'file':str(path.relative_to(ROOT)),
            'models':list(result.get('models',{})),
            'direct_equals_sequential':True,
            'physical_C5_certified':False,
            'running_input_from_BC_five_mass_eft1_inserted':False,
            'reason':'reference uses frozen 1e7 GeV inputs and SM-only background; cannot be relabelled as a physical UV/EFT1/EFT2 running result'}


def run_legacy_sequential(model:str,output:Path,resume:bool):
    """Optional historical conditional Y->X calculation, separate from Stage B."""
    dest_y='full_stack_legacy_y_v1'
    dest_x='full_stack_legacy_x_v1'
    path=ROOT/'output/studies/general_scalar_extended/continuous_full_scalar/continuous_v1'
    yd=path/dest_y
    xd=path/dest_x
    tag=safe_model(model)
    yfile=yd/f'{tag}_y_threshold_wilson.json'
    xfile=xd/f'{tag}_eft2_sequential.json'
    logs=[]
    if yfile.exists():
        require(resume and read_json(yfile).get('model') == model and
                read_json(yfile).get('finite_Y_threshold_delta_C_computed') is False,
                'existing Y output requires --resume and expected certification')
    if not yfile.exists():
        logs.append(run_subprocess(ROOT,[sys.executable,'-u','-X','utf8','-m',
            'scripts.RunT3BCYThresholdProjection','--model',model,
            '--run-id','continuous_v1','--output-id',dest_y],
            output/'logs'/f'{tag}_legacy_y.log'))
        require(logs[-1]['exit_code']==0,'legacy Y stage failed')
    if xfile.exists():
        require(resume, 'existing X output requires --resume')
    if not xfile.exists():
        logs.append(run_subprocess(ROOT,[sys.executable,'-u','-X','utf8','-m',
            'scripts.RunT3BCEFT2Sequential','--model',model,
            '--run-id','continuous_v1','--y-output-id',dest_y,'--output-id',dest_x],
            output/'logs'/f'{tag}_legacy_x.log'))
        require(logs[-1]['exit_code']==0,'legacy X stage failed')
    report=read_json(xfile)
    require(report.get('model')==model and
            report.get('physical_neutrino_mass_prediction_certified') is False and
            report.get('leading_log_sequential_candidate_completed') is True,
            'legacy sequential report certification changed')
    return {'status':'HistoricalConditionalY_XSequentialReplayed',
            'source':'continuous_v1 (NOT the new coupled five-mass trajectory)',
            'output':str(xfile.relative_to(ROOT)),
            'mnu_ev_descending':report['takagi_masses_ev_descending'],
            'physical_neutrino_masses_certified':False,'subprocesses':logs}


def preflight(models:list[str],registry:dict)->dict:
    out={}
    for m in models:
        try:
            deps=uv_dependencies(m,registry)
            vector,layout,quartics=load_uv_inputs(m,registry)
            eq=source_equations(m,registry,layout)
            require(set(eq.all_names)==set(read_json(UV_EXPORT/(safe_model(m)+'.json'))['report_betas']),
                    'incomplete native beta coverage')
            out[m]={'status':'ReadyCandidateUV','real_coordinates':eq.size,
                    'extra_yukawa':eq.extra_yukawa,'source_sha256':deps}
        except Exception as exc:
            out[m]={'status':'Blocked','reason':f'{type(exc).__name__}: {exc}'}
    return out


def parse_args():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument('--all',action='store_true',help='run selected stages')
    a.add_argument('--preflight',action='store_true',help='validate 16 model inputs without ODE solving')
    a.add_argument('--models',default='',help='comma-separated subset; default all 16')
    a.add_argument('--phase',choices=('all','uv','bc','matching','legacy'),default='all')
    a.add_argument('--output-id',default='t3_full_stack_overnight_v1')
    a.add_argument('--resume',action='store_true')
    a.add_argument('--extra-seed-fraction',type=float,default=0.3)
    a.add_argument('--uv-points',type=int,default=49)
    a.add_argument('--uv-rtol',type=float,default=2e-7)
    a.add_argument('--uv-atol',type=float,default=1e-10)
    a.add_argument('--convergence',action='store_true')
    a.add_argument('--bc-points',type=int,default=49)
    a.add_argument('--bc-higgs-mass-squared',type=float,default=-7800.)
    a.add_argument('--bc-mu2-fraction',type=float,default=.01)
    a.add_argument('--bc-mu2-phase',type=float,default=.6)
    a.add_argument('--include-legacy-sequential',action='store_true',
                   help='also replay original conditional Y/X benchmark on its old trajectory')
    a.add_argument('--fail-fast',action='store_true')
    return a.parse_args()


def run():
    args=parse_args()
    require(args.all or args.preflight,'select --all or --preflight')
    require(bool(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,48}',args.output_id)),
            'unsafe output-id')
    rootout=ROOT/'output/studies'/args.output_id
    registry=read_json(REGISTRY)
    allmodels=model_specifications(registry)
    models=allmodels if not args.models else [x.strip() for x in args.models.split(',')]
    require(len(set(models))==len(models) and all(x in allmodels for x in models),
            'invalid selected models')
    if args.preflight:
        checks=preflight(models,registry)
        for m in models:print(m,checks[m]['status'],flush=True)
        print('Preflight ready:',sum(v['status']=='ReadyCandidateUV' for v in checks.values()),'/',len(models))
        return 0 if all(v['status']=='ReadyCandidateUV' for v in checks.values()) else 2
    phases=({'uv','bc','matching'} if args.phase=='all' else {args.phase})
    if args.include_legacy_sequential:phases.add('legacy')
    rootout.mkdir(parents=True,exist_ok=True)
    report={'format':FORMAT,'status':'Running','source_commit':'current local files hashed individually',
            'models':models,'phases':sorted(phases),
            'parameters':{key:getattr(args,key) for key in ('extra_seed_fraction','uv_points','uv_rtol','uv_atol','convergence','bc_points','bc_higgs_mass_squared','bc_mu2_fraction','bc_mu2_phase')},
            'results':{},'errors':{},
            'physical_C5_certified':False,
            'full_16_model_sequential_matching_performed':False}
    def fail(name,e):
        report['errors'][name]={'type':type(e).__name__,'message':str(e),
                                'traceback':traceback.format_exc()[-3200:]}
        print('BLOCKED',name,str(e),file=sys.stderr,flush=True)
        if args.fail_fast:raise e
    for model in models:
        row={}
        report['results'][model]=row
        tag=safe_model(model)
        if 'uv' in phases:
            try:
                deps=uv_dependencies(model,registry)
                path=rootout/'uv'/f'{tag}_full_matrix_uv.json'
                if path.exists():
                    require(args.resume,'prior UV output exists; use --resume or new --output-id')
                    result=read_json(path)
                    require(result.get('source_sha256')==deps and result.get('status')=='CandidateFullMatrixRGBetaUVIntegrated' and
                            result.get('rtol')==args.uv_rtol and result.get('atol')==args.uv_atol and
                            len(result.get('trajectory',[]))==args.uv_points and
                            result.get('extra_yukawa_seed_fraction') == (args.extra_seed_fraction if model in EXTRA_MODELS else 0.),
                            'UV resume input/source mismatch')
                    require(not args.convergence or result.get('convergence',{}).get('passed') is True,
                            'UV resume missing convergence check')
                else:
                    result=solve_uv(model,registry=registry,seed_fraction=args.extra_seed_fraction,
                                    points=args.uv_points,rtol=args.uv_rtol,atol=args.uv_atol,
                                    convergence=args.convergence)
                    result['source_sha256']=deps
                    write_atomic(path,result)
                row['uv']={'status':result['status'],'path':str(path.relative_to(ROOT)),
                            'nfev':result['nfev'],
                            'extra_yukawa':result['extra_yukawa'],
                            'fermion_candidate_thresholds':result.get('candidate_full_UV_fermion_crossings')}
                print(model,'UV PASS nfev',result['nfev'],flush=True)
            except Exception as e:fail(model+'/UV',e)
        if model in BC_MODELS and 'bc' in phases:
            try:
                bc,bridge=stage_bc(model,output=rootout,points=args.bc_points,
                    mu2_fraction=args.bc_mu2_fraction,mu2_phase=args.bc_mu2_phase,
                    higgs_mass_squared=args.bc_higgs_mass_squared,resume=args.resume)
                bpath=rootout/'mass_basis'/f'{tag}_mixed_mass_basis.json'
                if bpath.exists():
                    oldbridge=read_json(bpath)
                    require(args.resume and oldbridge.get('status')==bridge['status'] and
                            oldbridge.get('source_five_mass_report_sha256') == bridge['source_five_mass_report_sha256'],
                            'existing mass-basis bridge source/status mismatch')
                else:write_atomic(bpath,bridge)
                row['bc']={'status':bc['status'],'path':str((rootout/'bc_full_mass'/f'{tag}_uv_to_eft1_coupled_five_mass_v1.json').relative_to(ROOT)),
                           'nfev_UV_EFT1':[bc['uv_solver_nfev'],bc['eft1_solver_nfev']],
                           'bridge':str(bpath.relative_to(ROOT)),
                           'physical_matching_certified':False}
                print(model,'B/C full mass and eigenbasis PASS',flush=True)
            except Exception as e:fail(model+'/BC',e)
        if model in BC_MODELS and 'legacy' in phases:
            try:
                row['legacy']=run_legacy_sequential(model,rootout,args.resume)
                print(model,'historical candidate Y/X PASS',flush=True)
            except Exception as e:fail(model+'/legacy',e)
        # Save progress after each model, including failed stages, without touching prior outputs.
        write_atomic(rootout/'progress.json',report,overwrite=True)
    if 'matching' in phases and any(m in BC_MODELS for m in models):
        try:
            report['restricted_bc_matching_reference']=reference_matching(rootout)
            print('B/C reference restricted C v23 PASS (frozen inputs)',flush=True)
        except Exception as e:fail('B/C/restricted_matching',e)
    report['status']='CompletedCandidateStagesPhysicalMatchingNotCertified' if not report['errors'] else 'PartialWithBlockedStages'
    report['remaining_physics']=[
        'Workstream A physical Matchete/RGBeta action and Yukawa basis certificate',
        'All 16 general Yukawa/bilinear physical matching coefficients beyond archived restricted B/C triangle',
        'Running UV/EFT1/EFT2 parameter and Wilson matching with finite field/EOM threshold shifts',
        'Physical matching callbacks after individual heavy-fermion and scalar eigenstate decoupling',
        'True full three-generation EFT Wilson anomalous dimension and all-model C5'
    ]
    write_atomic(rootout/'summary.json',report,overwrite=args.resume)
    write_atomic(rootout/'progress.json',report,overwrite=True)
    csv_path=rootout/'status.csv'
    require(args.resume or not csv_path.exists(),'status.csv exists; use --resume')
    with csv_path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=('model','uv','bc','legacy','any_error','physical_C5_certified'))
        w.writeheader()
        for model in models:
            row=report['results'][model]
            w.writerow({'model':model,'uv':row.get('uv',{}).get('status','NotRun'),
                        'bc':row.get('bc',{}).get('status','NotRun'),
                        'legacy':row.get('legacy',{}).get('status','NotRun'),
                        'any_error':any(k.startswith(model+'/') for k in report['errors']),
                        'physical_C5_certified':False})
    print('SUMMARY',report['status'],'at',rootout/'summary.json',flush=True)
    return 0 if not report['errors'] else 3


if __name__=='__main__':
    try:raise SystemExit(run())
    except (ValueError, OSError, KeyError, TypeError, json.JSONDecodeError) as e:
        print('BLOCKED:',str(e),file=sys.stderr,flush=True)
        raise SystemExit(2)
