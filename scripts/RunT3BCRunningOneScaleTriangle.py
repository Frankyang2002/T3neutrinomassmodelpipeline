"""Evaluate source-verified restricted B/C Matchete C5 at the RUNNING UV endpoint.

This is NOT the historical frozen v23 benchmark: MF, Yukawas, scalar masses,
mixing and quartics are all taken from the SAME coupled B/C 1e5 GeV UV state.
No y or MF is borrowed from 1e3 GeV scalar-only EFT1. Direct v22 ordered
matching is used only as an exploratory common-scale one-loop coefficient.
A complete renormalized UV->EFT1->EFT2 physical C5 is not provided.

Usage:
  python -u -X utf8 -m scripts.RunT3BCRunningOneScaleTriangle --all
  python -u -X utf8 -m scripts.RunT3BCRunningOneScaleTriangle --all --scan --resume
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import sys
import traceback

import numpy as np

from Numerical.experimental.T3BCRunningOneScaleTriangle import (
    MODELS,require,uv_boundary,assemble_ordered_c5,evaluate_one_scale,packed,
)

ROOT=Path(__file__).resolve().parents[1]
FORMAT='t3-bc-common-scale-running-uv-restricted-triangle-campaign-v1'
BC_PARENT=ROOT/'output/studies/t3_full_stack_overnight_v1/bc_full_mass'
OUTPUT=ROOT/'output/studies/t3_bc_running_one_scale_triangle_v1'
DIRECT=ROOT/'scripts/DirectVsSequentialOneLoop.experimental.py'
V18=ROOT/'scripts/EvaluateT3BCThreeFermionTriangle.experimental.py'
V22=ROOT/'output/audits/t3_bc_direct_vs_sequential_full_mass_v22/summary.json'
V23=ROOT/'output/audits/t3_bc_consolidated_leading_reference_v23/summary.json'


def digest(path:Path)->str:
    require(path.is_file(),'missing source '+str(path))
    return sha256(path.read_bytes()).hexdigest()


def read(path:Path)->dict:
    require(path.is_file(),'missing input '+str(path))
    obj=json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(obj,dict),'expected JSON object '+str(path))
    return obj


def write(path:Path,obj:dict,*,resume:bool)->None:
    require(resume or not path.exists(),'output exists; use --resume or another output-id: '+str(path))
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_name(path.name+'.tmp')
    require(not tmp.exists(),'stale temporary report '+str(tmp))
    try:
        with tmp.open('w',encoding='utf8') as f:
            json.dump(obj,f,indent=2,allow_nan=False)
            f.write('\n')
        os.replace(tmp,path)
    finally:
        if tmp.exists():tmp.unlink()


def import_source(path:Path,name:str):
    require(path.is_file(),'missing source script '+str(path))
    spec=importlib.util.spec_from_file_location(name,path)
    require(spec is not None and spec.loader is not None,
            'cannot import '+str(path))
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def model_tag(model:str)->str:
    require(model in MODELS,'only B_-1 and C_-1 supported')
    return model.replace('-','m')


def validate_bc_source(report:dict)->None:
    # Same source-digest policy as the overnight runner, including saved
    # native RGBeta UV/EFT1 betas and original fitted benchmark config.
    hashes=report.get('input_sha256')
    require(isinstance(hashes,dict) and len(hashes)>=6,
            'B/C source hash inventory absent')
    for key,value in hashes.items():
        path=Path(key)
        if not path.is_absolute():path=ROOT/path
        require(digest(path)==value,'B/C RGE input changed: '+str(path))


def load_bc(model:str, file:Path):
    report=read(file)
    validate_bc_source(report)
    boundary=uv_boundary(model,report)
    return report,boundary


def archived_source(model:str,v22_module,v18_module,v22_report:dict):
    stem=model_tag(model)
    require(v22_report.get('all_direct_vs_sequential_full_mass_agree') is True and
            v22_report.get('physical_C5') is None and
            v22_report.get('full_running_performed') is False,
            'v22 physical/provenance gate changed')
    row=v22_report['models'][stem]
    require(row.get('direct_vs_sequential_full_mass_agreement') is True,
            'v22 model consistency not passed')
    source,report,_,_=v22_module.source_gates(ROOT,stem)
    require(digest(source)==row['one_loop_direct_coefficient_sha256'] and
            digest(report)==row['all_heavy_direct_json_sha256'],
            'v22 source changed since archived direct vs sequential agreement')
    compiled,denom,counts,qcounts=v22_module.compile_source(
        source.read_text(encoding='utf-8-sig'),model)
    require(len(compiled)==480 and set(qcounts)==set(v22_module.Q_NAMES),
            'native LLHH source inventory changed')
    def coefficients(mf,mx,my,U,quartics):
        return v22_module.direct_ordered(compiled,denom,mf,mx,my,U,quartics)
    params=read(ROOT/'output/general_t3_mass_eigen_matchete'/f'{stem}_final_parameters.json')
    require(params['model']==model,'archived Matchete model mismatch')
    frozen=v18_module.load_frozen_reference(ROOT,stem)
    mx,my=map(float,params['heavy_masses_gev'])
    U=v18_module.as_matrix(params['unitary_U'],(2,2),'unitary_U')
    reference_quartics=v22_module.scenarios(frozen['quartics'][v22_module.Q_NAMES[0]])['baseline_only']
    reference_c5,_=assemble_ordered_c5(
       np.asarray(frozen['mf'],float),frozen['y1'],frozen['y2'],mx,my,U,
       reference_quartics,coefficients)
    stored=np.array([[complex(*z) for z in line]
                     for line in row['direct_full_3x3_C5_GeV_inverse_NOT_PHYSICAL']],complex)
    reference_rel=float(np.linalg.norm(reference_c5-stored)/max(np.linalg.norm(stored),1e-250))
    require(reference_rel<2e-8,'recompiled direct expression fails archived v22 frozen C5 parity')
    return coefficients,{
      'archived_frozen_C5_parity_relative':reference_rel,
      'native_direct_source_sha256':digest(source),
      'native_direct_matching_report_sha256':digest(report),
      'v22_comparison_report_sha256':digest(V22),
      'v18_evaluator_source_sha256':digest(V18),
      'restricted_direct_source_term_count':len(compiled),
      'restricted_quartic_inventory':qcounts,
      'electroweak_vev_reference_gev':float(frozen['vev_gev']),
      'reference_frozen_C5_3x3_gev_inverse_NOT_PHYSICAL':packed(stored)
    }


def running_direct_sequential_parity(model:str,boundary:dict,coefficients,
                                     v22_module,v18_module)->dict:
    """Re-evaluate independent F/Y/X triangle weights for RUNNING input U,
    quartics, MF, MX, MY. This is still within archived restricted theory.
    """
    stem=model_tag(model)
    v16=v22_module.import_module(
        ROOT/'scripts/AssembleT3BCFormalFiniteScaleEnvelope.experimental.py',
        't3_running_v16_'+stem)
    v17=v22_module.import_module(
        ROOT/'scripts/CompareT3BCFiniteEnvelopeToTriangle.experimental.py',
        't3_running_v17_'+stem)
    weights_sym=v22_module.sequential_weights(ROOT,stem,v16,v17,v18_module)
    scalar=boundary['scalar_basis']
    mx,my=map(float,np.sqrt(scalar.eigenvalues_sq))
    sub=v22_module.substitution(scalar.U,boundary['quartics'])
    numeric_weights=tuple(v18_module.polynomial_values(z,sub) for z in weights_sym)
    rows=[]
    for index,mf in enumerate(boundary['mf']):
        direct=coefficients(float(mf),mx,my,scalar.U,boundary['quartics'])
        sequential=v22_module.expected_ordered(numeric_weights,float(mf),mx,my,v17)
        residual,_,_=v22_module.consistency(direct,sequential)
        rows.append({'heavy_column':index,'mass_gev':float(mf),
                     'direct_sequential_relative_residual':residual,
                     'pass':residual<2e-8})
    maxres=max(row['direct_sequential_relative_residual'] for row in rows)
    require(maxres<2e-8,'running-input direct vs F/Y/X source triangle mismatch: '+str(maxres))
    return {'status':'RunningInputRestrictedFullMassTriangleParityPassed',
            'per_heavy_column':rows,
            'max_relative_residual':maxres,
            'both_sides_use_same_running_U_quartics_and_three_MF':True,
            'physical_EFT_counterterms_and_running_certified':False}


def matrix_json_to_numpy(obj):
    return np.asarray(obj['real'],float)+1j*np.asarray(obj['imag'],float)


def scan_cases(fraction:float,phase:float)->dict[str,tuple[float,float]]:
    return {'zero_independent_mu':(0.,0.),
            'phase_reversed_mu':(fraction,-phase),
            'larger_independent_mu':(min(2.5*fraction,.45),phase)}


def run_scenario(model:str,label:str,fraction:float,phase:float,output:Path,
                 points:int,*,resume:bool):
    from scripts.RunT3BCCoupledUVToEFT1 import run_model
    dirname=output/'bc_scenarios'/model_tag(model)/label
    target=dirname/f'{model_tag(model)}_uv_to_eft1_coupled_five_mass_v1.json'
    if target.exists():
        require(resume,'scenario exists; use --resume '+str(target))
        report=read(target)
        require(report.get('solver_success') is True and report.get('model')==model
                and report.get('initial_mu2_fraction_at_uv')==fraction
                and report.get('initial_mu2_phase_at_uv')==phase
                and len(report.get('uv_trajectory',[]))==points,
                'existing scenario does not match its requested parameters')
    else:
        report=run_model(model,output_dir=dirname,points=points,
                         mu2_fraction=fraction,mu2_phase=phase)
    validate_bc_source(report)
    return target,uv_boundary(model,report)


def run_model(model:str,*,source:Path,output:Path,scan:bool,
              points:int,resume:bool,v22_module,v18_module,v22_report:dict):
    report,boundary=load_bc(model,source)
    coefficients,reference=archived_source(model,v22_module,v18_module,v22_report)
    baseline=evaluate_one_scale(boundary,coefficients,
                               vev_gev=reference['electroweak_vev_reference_gev'])
    baseline['running_bc_source_sha256']=digest(source)
    baseline['frozen_source_parity']=reference
    baseline['running_input_direct_vs_sequential_parity']=running_direct_sequential_parity(
        model,boundary,coefficients,v22_module,v18_module)
    baseline['common_running_parameters_verified_same_UV_scale']=True
    baseline['uv_source_has_nonzero_extra_mixed_quadratic']=abs(boundary['scalar_basis'].mu_sq)>0
    baseline['matching_source_includes_extra_mixed_quadratic_as_independent_UV_vertex']=False
    dest=output/f'{model_tag(model)}_running_one_scale.json'
    # Explicit no-relabel policy: running direct mass-eigenstate extrapolation
    # cannot become physical matching by changing source filenames.
    if dest.exists():
        require(resume,'output exists '+str(dest))
        old=read(dest)
        require(old.get('running_bc_source_sha256')==digest(source) and
                old.get('frozen_source_parity',{}).get('native_direct_source_sha256')==
                    reference['native_direct_source_sha256'] and
                old.get('physical_C5_certified') is False,
                'candidate output source changed; choose a new output-id')
    else:
        write(dest,baseline,resume=False)
    cases=[]
    if scan:
        frac=float(report['initial_mu2_fraction_at_uv'])
        phase=float(report['initial_mu2_phase_at_uv'])
        a=np.array(baseline['reference_mnu_singular_ev_descending_NOT_PHYSICAL'],float)
        for label,(f,p) in scan_cases(frac,phase).items():
            source2,b2=run_scenario(model,label,f,p,output,points,resume=resume)
            val=evaluate_one_scale(b2,coefficients,
                        vev_gev=reference['electroweak_vev_reference_gev'])
            val['scenario_label']=label
            val['running_input_direct_vs_sequential_parity']=running_direct_sequential_parity(
                model,b2,coefficients,v22_module,v18_module)
            val['source_BC_report_sha256']=digest(source2)
            val['mixed_mass_source_input']={'mu2_fraction_at_UV':f,'phase':p}
            new=np.array(val['reference_mnu_singular_ev_descending_NOT_PHYSICAL'])
            c0=np.asarray(matrix_json_to_numpy(baseline['candidate_ordered_source_C5_3x3_gev_inverse_NOT_PHYSICAL']),complex)
            c1=np.asarray(matrix_json_to_numpy(val['candidate_ordered_source_C5_3x3_gev_inverse_NOT_PHYSICAL']),complex)
            val['relative_candidate_C5_change_vs_baseline']=float(
                np.linalg.norm(c1-c0)/max(np.linalg.norm(c0),1e-250))
            val['relative_reference_mass_singular_spectrum_change_vs_baseline']=float(
                np.linalg.norm(new-a)/max(np.linalg.norm(a),1e-250))
            file=output/'scenarios'/model_tag(model)/(label+'.json')
            if file.exists():
                require(resume and read(file).get('source_BC_report_sha256')==digest(source2),
                        'scenario output source changed')
            else:write(file,val,resume=False)
            cases.append({'label':label,'path':str(file.relative_to(ROOT)),
                          'relative_change':val['relative_reference_mass_singular_spectrum_change_vs_baseline'],
                          'masses_ev_NOT_PHYSICAL':val['reference_mnu_singular_ev_descending_NOT_PHYSICAL']})
    return {'model':model,'status':baseline['status'],
            'running_common_scale_gev':baseline['matching_input_scale_gev'],
            'source_five_mass_report':str(source.relative_to(ROOT)) if source.is_relative_to(ROOT) else str(source),
            'source_five_mass_sha256':digest(source),
            'result_path':str(dest.relative_to(ROOT)),
            'frozen_parity_residual':reference['archived_frozen_C5_parity_relative'],
            'reference_masses_ev_NOT_PHYSICAL':baseline['reference_mnu_singular_ev_descending_NOT_PHYSICAL'],
            'muX_gev':baseline['muX_gev'],'muY_gev':baseline['muY_gev'],
            'scenarios':cases,
            'physical_C5_certified':False}


def main() -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    choice=parser.add_mutually_exclusive_group(required=True)
    choice.add_argument('--all',action='store_true')
    choice.add_argument('--model',choices=MODELS)
    parser.add_argument('--bc-source-dir',type=Path,default=BC_PARENT)
    parser.add_argument('--output-id',default='t3_bc_running_one_scale_triangle_v1')
    parser.add_argument('--scan',action='store_true',help='rerun B/C native RGBeta for three independent mixed-mass settings')
    parser.add_argument('--bc-points',type=int,default=49)
    parser.add_argument('--resume',action='store_true')
    parser.add_argument('--fail-fast',action='store_true')
    args=parser.parse_args()
    require(re.fullmatch('[a-zA-Z0-9][a-zA-Z0-9_-]{0,54}',args.output_id) is not None,
            'invalid output-id')
    require(3<=args.bc_points<=200,'BC points must be in [3,200]')
    os.chdir(ROOT)  # existing v22 verification uses project-relative paths
    output=ROOT/'output/studies'/args.output_id
    output.mkdir(parents=True,exist_ok=True)
    v22=import_source(DIRECT,'t3_v22_direct_for_running_candidate')
    v18=import_source(V18,'t3_v18_reference_for_running_candidate')
    evidence=read(V22)
    if V23.is_file():
        v23=read(V23)
        require(v23.get('physical_C5_certified') is False and
                v23.get('direct_Matchete_all_heavy_equals_sequential_full_mass') is True,
                'v23 evidence status changed unexpectedly')
    models=MODELS if args.all else (args.model,)
    result={'format':FORMAT,'status':'Running','models':{},'errors':{},
            'source_mode':'Matchete native 480-term restricted full-mass triangle with running UV endpoint parameters',
            'full_sequential_EFT_matching_performed':False,
            'physical_C5_certified':False,
            'full_general_T3_Lagrangian_certified':False,
            'notes':'No EFT1/EFT2 Wilson evolution or finite field/EOM matching is inferred from these one-scale numbers.'}
    for model in models:
        tag=model_tag(model)
        src=args.bc_source_dir/f'{tag}_uv_to_eft1_coupled_five_mass_v1.json'
        try:
            row=run_model(model,source=src,output=output,scan=args.scan,
                points=args.bc_points,resume=args.resume,v22_module=v22,
                v18_module=v18,v22_report=evidence)
            result['models'][model]=row
            print(f'{model} PASS: common-scale={row["running_common_scale_gev"]:.7g} GeV; '
                  f'X,Y={row["muX_gev"]:.7g},{row["muY_gev"]:.7g} GeV; '
                  'restricted candidate triangle; physical C5=NOT CERTIFIED',flush=True)
        except Exception as exc:
            result['errors'][model]={'type':type(exc).__name__,'reason':str(exc),
                                      'traceback':traceback.format_exc()[-2500:]}
            print(f'{model} BLOCKED: {exc}',file=sys.stderr,flush=True)
            if args.fail_fast:raise
        write(output/'progress.json',result,resume=True)
    result['status']='CompleteRunningUVOneScaleCandidateNoPhysicalC5' if not result['errors'] else 'PartialBlocked'
    summary=output/'summary.json'
    write(summary,result,resume=args.resume)
    print('SUMMARY',result['status'],summary,flush=True)
    return 0 if not result['errors'] else 3

if __name__=='__main__':
    try:raise SystemExit(main())
    except (OSError,ValueError,KeyError,RuntimeError,TypeError,OverflowError) as exc:
        print('BLOCKED:',str(exc),file=sys.stderr,flush=True)
        raise SystemExit(2)
