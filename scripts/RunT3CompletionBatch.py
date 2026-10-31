#!/usr/bin/env python3
"""One-command batch for verified 16-model RGBeta UV + physical readiness report.

No production source modification, no overwrite of existing UV/EFT/Weinberg
outputs, no assumption of a physical Matchete/RGBeta coupling dictionary.

Usage:
  python scripts/RunT3CompletionBatch.py --self-test
  python scripts/RunT3CompletionBatch.py --preflight
  python scripts/RunT3CompletionBatch.py --six
  python scripts/RunT3CompletionBatch.py --all
  python scripts/RunT3CompletionBatch.py --model A_p0
  python scripts/RunT3CompletionBatch.py --all --no-phase-check
"""
from __future__ import annotations

import argparse
from datetime import datetime,timezone
from hashlib import sha256
import json
import os
from pathlib import Path
import sys
import traceback

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))

from Numerical.experimental.T3GeneralYukawaCoupledUV import (
    SPECS,EXTRAS,CONFIGS,REPORTS,OUTPUT,DynamicLayout,kind,load_case,
    prepare_state,endpoint_for,integrate,run_model,require,self_test as uv_selftest,
)
from Numerical.experimental.T3GeneralYukawaPhaseAudit import (
    test_rephase_covariance,self_test as phase_selftest,
)

PRODUCTION=ROOT/'RGE/running/rgbeta/T3RGBetaModel.wl'
EXPORT_SUMMARY=ROOT/'Reports/output/workstream_a/integrated_general_yukawa_uv/summary.json'
PHYSICAL_PROOFS={
    'formal_A0':ROOT/'Reports/output/workstream_a/yukawa_a0_action_bridge/a0_action_bridge.json',
    'native_spinor':ROOT/'Reports/output/workstream_a/yukawa_spinor_phase/a0_spinor_phase.json',
    'mass_Ward':ROOT/'Reports/output/workstream_a/a0_dirac_mass_ward/a0_dirac_mass_ward.json',
    'general_RGBeta_16':EXPORT_SUMMARY,
    'one_step_6':ROOT/'Reports/output/workstream_a/general_yukawa_leading_log/summary.json',
}


def read_json(path:Path)->dict:
    x=json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(x,dict),f'{path} must contain a JSON object')
    return x


def source_fingerprint()->str:
    require(PRODUCTION.is_file(),'Missing production RGBeta source')
    content=PRODUCTION.read_text(encoding='utf-8-sig')
    require('T3RGBetaAddGeneralFermionYukawas[' in content and
            'T3_EXPERIMENTAL_GENERAL_FERMIONS' in content,
            'General-Yukawa source not installed')
    require('eps[T3RGBetaRep[dF], #1, #2]' in content,
            'Heavy SU(2) doublet mass epsilon correction not installed')
    return sha256(PRODUCTION.read_bytes()).hexdigest()


def verify_export_freshness(*, strict:bool=True)->dict:
    actual=source_fingerprint()
    require(EXPORT_SUMMARY.is_file(),'Missing 16-model production export summary. Rerun scripts/RunT3ExtraYukawaStack.py --all')
    summary=read_json(EXPORT_SUMMARY)
    require(summary.get('mode')=='all' and summary.get('all_requested_models_succeeded') is True,
            'RGBeta export summary does not certify all 16 models')
    recorded=summary.get('local_model_sha256')
    if recorded is not None:
        require(recorded==actual,'Production model changed since RGBeta exports; regenerate all 16')
    elif strict:
        raise ValueError('RGBeta export summary has no model source SHA')
    names={r.get('model') for r in summary.get('rows',[]) if r.get('status')=='VerifiedRGBetaUVReport'}
    require(names==set(SPECS),'16-model RGBeta export summary incomplete or duplicate')
    return {'status':'16RGExportsSourceAndMetadataValidated',
            'production_rgbeta_sha256':actual,
            'export_source_sha_match':recorded==actual,
            'all16_report_count':len(names)}


def physical_gate()->dict:
    found={}
    for key,path in PHYSICAL_PROOFS.items():
        try:
            obj=read_json(path)
            found[key]={'path':str(path.relative_to(ROOT)),
                        'status':obj.get('status','NoStatus'),'exists':True}
        except (OSError,ValueError,KeyError) as exc:
            found[key]={'path':str(path.relative_to(ROOT)),
                        'exists':False,'error':str(exc)}
    return {
        'status':'PhysicalCrossPackageDictionaryNotCertified',
        'source_reports':found,
        'fixed_convention_in_previous_tests':'C=+i gamma2 gamma0, epsilon12=+1, conditional eta=+1',
        'remaining_physical_dependencies':[
            'Native absolute sign of charge conjugation and RGBeta SU2 tensor axes in common field basis',
            'One simultaneous y1/y2/extra-Yukawa/MF physical action dictionary across doublet and triplet components',
            'Matched extended heavy-sector finite threshold coefficients and light scalar threshold operators',
        ],
        'legacy_neutrino_predictions_not_equivalent_to_extended_neutrino_predictions':True,
        'extended_neutrino_mass_matching_executed':False,
        'extended_neutrino_oscillation_observables_executed':False,
    }


def atomic_json(path:Path,data:dict):
    path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_name(path.name+'.tmp')
    temp.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
    temp.replace(path)


def preflight_one(model:str, *, seed_factor:float, phase_check:bool=True)->dict:
    from Numerical.experimental.T3GeneralYukawaCoupledUV import (
        evaluate_betas,pack_derivative
    )
    data,cfg=load_case(model)
    st=prepare_state(model,data,cfg,seed_factor)
    layout=DynamicLayout(data['report_betas'].keys())
    y=layout.pack(st)
    beta=evaluate_betas(data['report_betas'],st)
    packed=pack_derivative(layout,beta)
    require(y.shape==packed.shape,'RHS shape mismatch')
    a,b=endpoint_for(cfg,st)
    phase=None
    if phase_check and SPECS[model][3]==-1:
        phase={'status':'NotApplicableMajoranaFermion',
               'reason':'Neutral self-conjugate heavy fermion has a Majorana mass; independent F/Fc vectorlike U(1) rephase is not a symmetry'}
    elif phase_check:
        phase=test_rephase_covariance(data['report_betas'],st)
    return {'status':'NativeExportFullRHSNumericallyEvaluated',
            'model':model,'couplings':len(layout.names),'real_ode_dimension':layout.size,
            'uv_gev':a,'endpoint_gev':b,'rhs_euclidean':float(__import__('numpy').linalg.norm(packed)),
            'phase_covariance':phase,
            'physical_matchete_dictionary_certified':False}


def self_test():
    uv_selftest();phase_selftest()
    require(len(SPECS)==16 and len(EXTRAS)==6,'Model inventory invalid')
    for model in SPECS:
        d1,d2,df,alpha=SPECS[model]
        require(all(x in (1,2,3) for x in (d1,d2,df)),'Unexpected SU2 representation')
        expected=(model[0] in ('A','D','E') and alpha in (-2,0))
        require((model in EXTRAS)==expected,'Wrong extra Yukawa assignment')
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        p=Path(td)/'summary.json'
        atomic_json(p,{'status':'OK','values':[1,2,3]})
        require(read_json(p)['values']==[1,2,3],'Atomic JSON roundtrip failed')
    print('PASS: complete catalogue, source-gate structure, artifact serialization')


def run(mode:str,*,selected:str|None=None,seed_factor:float=.3,
        n_save:int=9,phase_check:bool=True,continue_on_error:bool=True,plots:bool=True)->int:
    fresh=verify_export_freshness()
    models=([selected] if mode=='model' else
            list(EXTRAS) if mode=='six' else list(SPECS))
    require(bool(models),'Empty model selection')
    OUTPUT.mkdir(parents=True,exist_ok=True)
    records=[]
    for i,model in enumerate(models,1):
        record={'model':model}
        try:
            diagnostic=preflight_one(model,seed_factor=seed_factor,phase_check=phase_check)
            record['preflight']=diagnostic
            if mode!='preflight':
                detail_path=OUTPUT/'models'/(model+'.json')
                # Always recalculate, never reuse stale trajectory under new input.
                full=run_model(model,seed_factor=seed_factor,n_save=n_save,output=detail_path)
                record['status']='CoupledGeneralYukawaRGBetaOnlyUVPassed'
                record['nfev']=full['nfev'];record['real_ode_dimension']=full['real_ode_size']
                record['trajectory_report']=str(detail_path.relative_to(ROOT))
                record['extra_yukawa']=full['extra_yukawa']
                record['control']=full['seeded_vs_zero_control']
                print(f'[{i}/{len(models)}] PASS {model}: {full["real_ode_size"]} real ODE variables; '
                      f'{full["nfev"]} RHS evaluations; extra={full["extra_yukawa"]}',flush=True)
            else:
                record['status']='FullRGBetaRHSPreflightPassed'
                print(f'[{i}/{len(models)}] RHS PASS {model}: '
                      f'{diagnostic["real_ode_dimension"]} real variables',flush=True)
        except Exception as exc:
            record['status']='Blocked'
            record['error']=f'{type(exc).__name__}: {exc}'
            record['exception_type']=type(exc).__name__
            print(f'[{i}/{len(models)}] BLOCKED {model}: {exc}',file=sys.stderr,flush=True)
        records.append(record)
        atomic_json(OUTPUT/'summary.json',{
            'format':'t3-all-remaining-large-stack-v1',
            'status':'RunningOrPartial',
            'mode':mode,'generated_at_utc':datetime.now(timezone.utc).isoformat(),
            'source_fingerprint':fresh,'rows':records,
            'physical_gate':physical_gate(),
        })
        if record['status']=='Blocked' and not continue_on_error:break
    passed=all(row['status']!='Blocked' for row in records) and len(records)==len(models)
    phase_warnings=[r['model'] for r in records if r.get('preflight',{}).get('phase_covariance',{}).get('status')=='RGBetaVectorlikeFieldRephaseDiagnosticFailed']
    report={
        'phase_covariance_warnings':phase_warnings,
        'format':'t3-all-remaining-large-stack-v1',
        'status':('AllRequestedFullCoupledUVRGBetaOnlyPassed' if passed and mode!='preflight'
                  else 'AllRequestedFullRGBetaRHSPassed' if passed else
                  'PartialFullCoupledUVRGBetaOnlyBlocked'),
        'mode':mode,'generated_at_utc':datetime.now(timezone.utc).isoformat(),
        'source_fingerprint':fresh,
        'n_requested_models':len(models),'n_passed_models':sum(r['status']!='Blocked' for r in records),
        'rows':records,'physical_gate':physical_gate(),
        'full_ode_running_performed':mode!='preflight',
        'neutrino_mass_prediction_performed':False,
        'warning':'RGBeta-only coordinates from legacy fitted reference numerical values; no Matchete map, threshold EFT matching, or extended neutrino-mass calculation.',
    }
    atomic_json(OUTPUT/'summary.json',report)
    if mode!='preflight':
        try:
            from scripts.SummarizeT3CompletionBatch import summarize
            diag=summarize(OUTPUT,plot=plots)
            report['aggregate_csv']=str(Path(diag['csv']).relative_to(ROOT))
            atomic_json(OUTPUT/'summary.json',report)
        except (ValueError,OSError,KeyError,ImportError) as exc:
            report['aggregation_warning']=f'{type(exc).__name__}: {exc}'
            atomic_json(OUTPUT/'summary.json',report)
            print('WARNING: CSV/figure aggregation failed: '+str(exc),file=sys.stderr)
    print(f'\nSUMMARY {report["status"]}: {report["n_passed_models"]}/{len(models)} models')
    print('Report: '+str(OUTPUT/'summary.json'))
    print('Physical extended neutrino-mass matching: BLOCKED on absolute cross-package dictionary.')
    return 0 if passed else 4


def main()->int:
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--self-test',action='store_true')
    group.add_argument('--preflight',action='store_true')
    group.add_argument('--six',action='store_true')
    group.add_argument('--all',action='store_true')
    group.add_argument('--model',choices=list(SPECS))
    parser.add_argument('--seed-factor',type=float,default=0.3)
    parser.add_argument('--n-save',type=int,default=9)
    parser.add_argument('--no-phase-check',action='store_true')
    parser.add_argument('--no-plots',action='store_true',help='Skip optional figures; still emit CSV')
    parser.add_argument('--stop-on-error',action='store_true')
    a=parser.parse_args()
    if a.self_test:self_test();return 0
    mode='model' if a.model else 'preflight' if a.preflight else 'six' if a.six else 'all'
    return run(mode,selected=a.model,seed_factor=a.seed_factor,n_save=a.n_save,
               phase_check=not a.no_phase_check,continue_on_error=not a.stop_on_error,
               plots=not a.no_plots)

if __name__=='__main__':
    try:sys.exit(main())
    except (ValueError,OSError,KeyError,TypeError,json.JSONDecodeError) as exc:
        print('BLOCKED '+type(exc).__name__+': '+str(exc),file=sys.stderr)
        sys.exit(3)
