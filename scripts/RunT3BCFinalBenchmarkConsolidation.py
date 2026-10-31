"""B/C batch 3: compare five-mass SM-only proxy to historical EFT1 result.

Run: python -m scripts.RunT3BCFinalBenchmarkConsolidation
No historical EFT1 coefficients are injected into the five-mass trajectory.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import re
import sys
import traceback
from Numerical.experimental import T3BCFinalBenchmarkConsolidation as check
from Numerical.experimental import T3BCMegaPhysicalBridge as bridge

ROOT=Path(__file__).resolve().parents[1]
FIVE=ROOT/'output/studies/t3_full_stack_overnight_v1/bc_full_mass'
ONE=ROOT/'output/studies/t3_bc_running_one_scale_triangle_v1'
HIST=ROOT/'output/studies/t3_bc_eft1_grouped_scalar_batch_v1/summary.json'
B1=ROOT/'output/studies/t3_bc_fermion_threshold_batch_v1/summary.json'
OUTPUT=ROOT/'output/studies/t3_bc_final_benchmark_consolidation_v1'


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def require(ok,msg):check.require(ok,msg)

def put(path,obj,resume):
    require(resume or not path.exists(),'refusing to overwrite '+str(path))
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix('.json.tmp')
    tmp.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    tmp.replace(path)


def run(out=OUTPUT,*,rtol=1e-8,resume=False,history=HIST,five=FIVE,one=ONE,batch1=B1):
    require(1e-11<rtol<1e-5,'invalid numerical tolerance')
    audit=bridge.load(batch1)
    require(audit.get('status')=='ConditionalUVToEFT1SourceBoundaryAudited' and
            audit.get('physical_C5_certified') is not True,
            'missing successful Batch 1 source audit')
    historical=bridge.load(history)
    results={};errors={}
    for model in check.MODELS:
        stem=model.replace('-','m')
        src=five/(stem+'_uv_to_eft1_coupled_five_mass_v1.json')
        runpath=one/(stem+'_running_one_scale.json')
        try:
            historical_row=check.validate_grouped(model,historical)
            source=bridge.load(src);running=bridge.load(runpath)
            _evidence,_sm=bridge.check_running_pair(model,source,running,
                                             source_path=src,running_path=runpath)
            # Reuse source-native SM evolution; do not substitute EFT1 fields.
            a=bridge.run_sm_only_proxy(model,source,running,rtol=rtol)
            b=bridge.run_sm_only_proxy(model,source,running,rtol=rtol/10)
            row=check.compare(model,historical_row,running,a,b)
            row['source_sha256']={
                'five_mass_UV_EFT1':sha(src),'one_scale_triangle':sha(runpath),
                'historical_EFT1_grouped':sha(history),'fermion_threshold_audit':sha(batch1)}
            row['same_uv_fermion_endpoint_verified']=True
            row['running_source_parity_verified']=True
            row['sm_only_evolution_solver_calls']=[int(a['solver_nfev']),int(b['solver_nfev'])]
            results[model]=row
            print(model,'SM-only 100 GeV proxy (eV)',row['same_five_mass_UV_input_SM_only_100gev_proxy_masses_ev_descending'],
                  'relative convergence',row['SM_only_proxy_convergence_relative_C5'],flush=True)
        except Exception as exc:
            errors[model]={'type':type(exc).__name__,'reason':str(exc),
                           'traceback':traceback.format_exc(limit=6)}
            print(model,'BLOCKED',str(exc),flush=True)
    report={'format':'t3-bc-five-mass-vs-historical-eft1-comparison-v1',
        'status':check.STATUS if not errors else 'PartialBlocked',
        'models':results,'errors':errors,
        'one_intermediate_EFT_only':True,
        'historical_and_five_mass_trajectories_identical':False,
        'historical_eft1_wilson_boundary_used_in_new_five_mass_trajectory':False,
        'same_trajectory_EFT1_wilson_running_completed':False,
        'finite_scalar_threshold_matching_computed':False,
        'physical_C5_certified':False,'physical_mass_prediction_certified':False,
        'interpretation':'Five-mass UV SM-only evolution is a reference, not EFT1 transport; historical grouped EFT1 values belong to different UV inputs.'}
    put(out/'summary.json',report,resume)
    print('SUMMARY',report['status'],'physical C5=NOT CERTIFIED',out/'summary.json',flush=True)
    return 0 if not errors else 2


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output-id',default='t3_bc_final_benchmark_consolidation_v1')
    p.add_argument('--rtol',type=float,default=1e-8)
    p.add_argument('--resume',action='store_true')
    p.add_argument('--history',type=Path,default=HIST)
    p.add_argument('--batch1',type=Path,default=B1)
    a=p.parse_args()
    require(re.fullmatch('[A-Za-z0-9][A-Za-z0-9_-]{0,65}',a.output_id) is not None,
            'unsafe output ID')
    return run(out=ROOT/'output/studies'/a.output_id,rtol=a.rtol,
               resume=a.resume,history=a.history,batch1=a.batch1)

if __name__=='__main__':
    try:raise SystemExit(main())
    except (ValueError,OSError,KeyError,TypeError) as ex:
        print('BLOCKED',ex,file=sys.stderr);raise SystemExit(2)
