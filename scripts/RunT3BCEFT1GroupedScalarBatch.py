"""B/C EFT1 coupled Wilson evolution through ONE grouped scalar threshold.

Reuses historical continuous_v1 coupled EFT1 trajectory, not the new five-mass
UV trajectory. No second intermediate EFT, no Y->X running, no invented finite
matching. Grouped scalar decoupling at existing EFT1 endpoint is conditional.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
import traceback
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
MODELS=('B_-1','C_-1')
FORMAT='t3-bc-single-eft1-grouped-scalar-conditional-v1'
STATUS='HistoricalEFT1CoupledWilsonGroupedScalarCandidateNoFiniteMatching'

def require(ok,msg):
    if not ok:raise ValueError('B/C grouped scalar: '+msg)

def read(p):
    require(p.is_file(),'missing '+str(p))
    v=json.loads(p.read_text(encoding='utf-8-sig'))
    require(isinstance(v,dict),'invalid JSON object '+str(p))
    return v

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def unpack(obj):
    require(isinstance(obj,dict) and set(obj)=={'real','imag'},'packed complex matrix schema')
    a=np.asarray(obj['real'],float);b=np.asarray(obj['imag'],float)
    require(a.shape==b.shape==(3,3) and np.isfinite(a).all() and np.isfinite(b).all(),'invalid 3x3')
    return a+1j*b

def check_boundary(batch1:dict):
    require(batch1.get('status')=='ConditionalUVToEFT1SourceBoundaryAudited',
            'Batch 1 audit not successful')
    # Batch1 is a provenance-only input. It is NOT injected numerically into
    # historical EFT1 running, which uses its own saved UV boundary.
    require(batch1.get('physical_C5_certified') is not True,'unexpected physical certificate')

def validate_output(model, result, *,expected_samples, expected_rtol):
    from scripts import RunT3BCCoupledWilsonSequential as coupled
    require(result.get('model')==model and result.get('format')==coupled.FORMAT
            and result.get('status')==coupled.STATUS,'wrong coupled model/status')
    require(result.get('one_generation_gauge_index_wilson_evolution_completed') is True,
            'EFT1 Wilson evolution missing')
    require(result.get('eft1_samples')==expected_samples
            and float(result.get('wilson_rtol'))==expected_rtol,'solver configuration mismatch')
    require(result.get('finite_scalar_threshold_matching_included') is False
            and result.get('full_general_scalar_matching_certified') is False
            and result.get('full_sequential_matching_performed') is False,
            'upstream unexpectedly claims full matching')
    require(result.get('additional_quartic_mapping_physically_certified') is False
            and result.get('full_three_generation_wilson_anomalous_dimension_included') is False,
            'upstream caveats changed')
    matrix=unpack(result['coupled_C5_at_scalar_threshold_gev_inverse'])
    denom=max(np.linalg.norm(matrix),1.e-250)
    require(np.linalg.norm(matrix-matrix.T)/denom<1e-10,
            'candidate C5 not symmetric')
    muF=float(result['hard_threshold_grouped_scale_gev']);muS=float(result['scalar_threshold_gev'])
    require(math.isfinite(muF) and math.isfinite(muS) and muF>muS>0,
            'invalid grouped threshold hierarchy')
    mass=np.asarray(result['takagi_masses_ev_descending'],float)
    require(mass.shape==(3,) and np.isfinite(mass).all()
            and np.all(mass>=0) and np.all(np.diff(mass)<=0),'invalid candidate mass spectrum')
    require(float(result['scalar_threshold_weinberg_residual'])<1e-4,
            'large Weinberg projector residual')
    for key in ('direct_source_parity','stage17_extra_source_parity'):
        require(float(result['source_normalization_checks'][key])<3e-4,
                'source parity failed '+key)
    return {'model':model,'mu_F_gev':muF,'mu_S_grouped_gev':muS,
            'candidate_C5_at_grouped_scalar_threshold_gev_inverse':result['coupled_C5_at_scalar_threshold_gev_inverse'],
            'candidate_masses_ev_descending':mass.tolist(),
            'weinberg_projection_residual':float(result['scalar_threshold_weinberg_residual']),
            'direct_source_parity':float(result['source_normalization_checks']['direct_source_parity']),
            'stage17_extra_source_parity':float(result['source_normalization_checks']['stage17_extra_source_parity']),
            'eft1_scalar_nfev':int(result['scalar_ode_nfev']),
            'eft1_wilson_nfev':int(result['wilson_ode_nfev']),
            'scalar_matching_type':'grouped scalar threshold at EFT1 endpoint, no Y->X EFT',
            'finite_scalar_threshold_correction':'UNKNOWN_NOT_ZERO',
            'physical_C5_certified':False,'full_matching_certified':False}

def run(*,folder,out,batch1_path,samples=65,rtol=5e-7,resume=False,run_coupled=False):
    from scripts import RunT3BCCoupledWilsonSequential as coupled
    audit=read(batch1_path);check_boundary(audit)
    require(folder.is_dir(),'historical continuous EFT1 folder missing '+str(folder))
    require(17<=samples<=129 and samples%2==1,'samples must be odd 17..129')
    require(1e-10<rtol<1e-3,'rtol out of range')
    require(not out.exists() or resume,'output exists; use --resume or new output ID')
    src=folder/coupled.DEFAULT_OUTPUT_ID
    results={};errors={}
    for model in MODELS:
        try:
            tag=model.replace('-','m')
            p=src/(tag+'_coupled_wilson.json')
            if p.is_file():
                saved=read(p)
                # Require matching source sample count, no silent reuse of other solver settings.
                row=validate_output(model,saved,expected_samples=samples,expected_rtol=rtol)
            else:
                require(run_coupled,'missing coupled EFT1 result '+str(p)+
                        '; rerun with --run-coupled (requires prior direct/extra reports)')
                result=coupled.calculate(model,folder=folder,out=src,samples=samples,
                                         rtol=rtol,resume=resume)
                row=validate_output(model,result,expected_samples=samples,expected_rtol=rtol)
            row['historical_eft1_coupled_source_sha256']=digest(p)
            results[model]=row
            print(model,'grouped muS=',row['mu_S_grouped_gev'],
                  'masses eV=',row['candidate_masses_ev_descending'],flush=True)
        except Exception as e:
            errors[model]={'type':type(e).__name__,'reason':str(e),
                           'traceback':traceback.format_exc(limit=5)}
            print(model,'BLOCKED',str(e),flush=True)
    output={'format':FORMAT,'status':STATUS if not errors else 'PartialBlocked',
            'models':results,'errors':errors,
            'batch1_audit_sha256':digest(batch1_path),
            'input_historical_continuous_folder':str(folder),
            'same_trajectory_as_batch1_five_mass_running_verified':False,
            'batch1_boundary_numerically_inserted':False,
            'one_intermediate_EFT_only':True,'separate_X_only_EFT_created':False,
            'finite_scalar_threshold_matching_computed':False,
            'independent_F_EFT1_counterterms_computed':False,
            'complete_three_family_EFT1_Wilson_running':False,
            'physical_C5_certified':False,'physical_neutrino_mass_prediction':False,
            'next_gate':'Recompute the same EFT1 Wilson system using the latest five-mass UV boundary and close finite scalar matching.'}
    out.mkdir(parents=True,exist_ok=True)
    target=out/'summary.json'
    require(not target.exists() or resume,'summary exists; choose new output ID')
    tmp=target.with_suffix('.json.tmp')
    tmp.write_text(json.dumps(output,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    tmp.replace(target)
    print('SUMMARY',output['status'],'physical C5=NOT CERTIFIED','output:',target)
    return 2 if errors else 0

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--continuous-run-id',default='continuous_v1')
    parser.add_argument('--output-id',default='t3_bc_eft1_grouped_scalar_batch_v1')
    parser.add_argument('--batch1-path',type=Path,default=ROOT/'output/studies/t3_bc_fermion_threshold_batch_v1/summary.json')
    parser.add_argument('--samples',type=int,default=65)
    parser.add_argument('--rtol',type=float,default=5e-7)
    parser.add_argument('--resume',action='store_true')
    parser.add_argument('--run-coupled',action='store_true')
    a=parser.parse_args()
    import re
    require(bool(re.fullmatch('[A-Za-z0-9][A-Za-z0-9_-]{0,65}',a.continuous_run_id)),'invalid continuous id')
    require(bool(re.fullmatch('[A-Za-z0-9][A-Za-z0-9_-]{0,65}',a.output_id)),'invalid output id')
    from scripts import RunT3FullScalarContinuous as continuous
    return run(folder=continuous.OUTPUT/a.continuous_run_id,
               out=ROOT/'output/studies'/a.output_id,batch1_path=a.batch1_path,
               samples=a.samples,rtol=a.rtol,resume=a.resume,run_coupled=a.run_coupled)
if __name__=='__main__':
    try:raise SystemExit(main())
    except (ValueError,OSError,ImportError) as e:
        print('BLOCKED:',e,file=sys.stderr);raise SystemExit(2)
