"""Run native experimental Matchete EFT1(H,X,Y) -> EFT2(H,X), Y heavy only.

The upstream canonical F-first EOMEFT, CG reconstruction and strict
CheckLagrangian machinery are reused from the installed project. The original
Wolfram program is never modified. A complete derived standalone Wolfram
program is written beside it, after strict source-structure checks.

The saved mass-eigen parameter point is the project's symbolic probe point,
NOT the continuous_v1 running mu_Y benchmark. Raw operators MUST NOT be
inserted into numerical Wilson coefficients without matching their coupling
conventions, scale, subtraction scheme, and complete EFT basis.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
from hashlib import sha256
from pathlib import Path
from typing import Any

from Numerical.experimental.T3BCYOnlyMatcheteBuilder import (
    BASE_NAME, GENERATED_NAME, build,
)
from scripts.RunGeneralT3EFT1ScalarThresholdProbe import (
    _verified_tree_input,
)
from scripts.RunGeneralT3MassEigenMatchete import ROOT, INPUT, OUTPUT, prepare_model
from scripts.RunGeneralT3FermionFirstEFT1Tree import OUTDIR as TREE_OUTDIR, _prior_certificate
from scripts.RunGeneralT3FermionFirstEFT1OneLoop import OUTDIR as ONELOOP_OUTDIR

MODELS=('B_-1','C_-1')
NATIVE=ROOT/'Lagrangian'/BASE_NAME
GENERATED=ROOT/'Lagrangian'/GENERATED_NAME
OUT=OUTPUT/'eft1_scalar_Y_only_probe'
STATUS='YOnlyOneLoopRawEFT2SectorsExtractedNotPhysicalMatching'
FORMAT='t3-bc-y-only-matchete-eigenfield-probe-v1'


def require(test: bool, message: str) -> None:
    if not test:
        raise RuntimeError('Y-only Matchete runner: '+message)


def sha(path: Path) -> str:
    return sha256(Path(path).read_bytes()).hexdigest()


def json_load(path:Path) -> dict:
    data=json.loads(Path(path).read_text(encoding='utf-8-sig'))
    require(isinstance(data,dict), 'expected JSON dictionary '+str(path))
    return data


def run_model(model:str, *, out:Path, scope:str, preview:bool,
              wolframscript:str, timeout_sec:int)->dict[str,Any]:
    require(model in MODELS,'unexpected model')
    require(scope in ('full','leading-c5'),'unexpected matching scope')
    require(timeout_sec>0,'invalid timeout')
    source_info=build(NATIVE, GENERATED)
    prepared=prepare_model(model,which='final',folder=INPUT,outdir=OUTPUT)
    mass_basis_cert=_prior_certificate(model,'final',OUTPUT,prepared['payload']['unitary_U'])
    canonical,tree_cert,one_loop_cert=_verified_tree_input(
        model,'final',TREE_OUTDIR,ONELOOP_OUTDIR)
    require(canonical.is_file() and tree_cert.is_file() and one_loop_cert.is_file(),
            'missing certified F-first inputs')
    relative=lambda p: str(Path(p).relative_to(ROOT))
    stem=model.replace('+','p').replace('-','m')+'_final_y_only_'+scope.replace('-','_')
    dest=out/(stem+'.json')
    log=out/(stem+'.log')
    state={
        'format':'t3-bc-y-only-matchete-driver-v1',
        'model':model,'scope':scope,'status':'PreflightCompletedMatcheteNotRun',
        'mass_basis_certificate':relative(mass_basis_cert),
        'saved_F_first_eft1_source':relative(canonical),
        'saved_F_first_tree_certificate':relative(tree_cert),
        'saved_F_first_one_loop_certificate':relative(one_loop_cert),
        'parameter_source':relative(prepared['parameters']),
        'matchete_native_source':source_info,
        'matchete_parameter_point_equals_running_muY':False,
        'renormalized_sequential_Y_matching_certified':False,
        'physical_matching_certified':False,
        'running_numeric_coupling_substitution_performed':False,
        'full_EFT1_operator_matching_combined':False,
    }
    out.mkdir(parents=True,exist_ok=True)
    if preview:
        preflight=out/(stem+'_preflight.json')
        require(not preflight.exists(), 'preflight output exists: '+str(preflight))
        preflight.write_text(json.dumps(state,indent=2)+'\n',encoding='utf8')
        print(f'{model}: Y-only Matchete preflight passed; generated {relative(GENERATED)}',flush=True)
        return state
    require(not dest.exists() and not log.exists(),
            'output already exists; use a new --output-dir')
    require(shutil.which(wolframscript) or Path(wolframscript).is_file(),
            'wolframscript executable missing')
    command=[wolframscript,'-file',str(GENERATED),
             str(prepared['parameters']),str(canonical),str(dest)]
    env=dict(os.environ)
    env.pop('T3_EXPERIMENTAL_GENERAL_QUARTICS',None)
    env.pop('T3_EFT1_LEADING_C5_ONLY',None)
    if scope=='leading-c5':
        env['T3_EFT1_LEADING_C5_ONLY']='1'
    try:
        run=subprocess.run(command,cwd=ROOT,env=env,text=True,capture_output=True,
                           timeout=timeout_sec,check=False,errors='replace')
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError(f'{model}: Wolfram execution timeout; progress file may be at {dest.with_suffix(".progress.json")}') from exc
    log.write_text('COMMAND: '+repr(command)+'\nEXIT: '+str(run.returncode)+
                   '\nSTDOUT:\n'+run.stdout+'\nSTDERR:\n'+run.stderr,encoding='utf8')
    require(dest.is_file(), f'{model}: no Matchete report; see {log}')
    report=json_load(dest)
    require(run.returncode==0 and report.get('status')==STATUS and
            report.get('format')==FORMAT and report.get('model')==model,
            f'{model}: Matchete returned incomplete or failed stage; see {log}; status={report.get("status")}')
    require(report.get('heavy_scalar_fields')==['NewScalarY'] and
            report.get('light_scalar_fields')==['H','NewScalarX'] and
            report.get('heavy_Y_absent_and_light_X_retained_in_matched_eft') is True and
            report.get('tree_matching_completed') is True and
            report.get('one_loop_matching_completed') is True and
            report.get('physical_matching_performed') is False,
            f'{model}: incomplete heavy/light field-inventory checks')
    counts=report.get('raw_EFT2_operator_sector_counts')
    require(isinstance(counts,dict) and set(counts)=={'LLHH','LLHX','LLXX'} and
            all(isinstance(v,int) and v>=0 for v in counts.values()),
            f'{model}: invalid EFT2 raw operator inventory')
    for suffix in ('_sm_tree.inputform.txt', '_sm_tree_plus_loop.inputform.txt',
                   '_scalar_one_loop_correction.inputform.txt',
                   '_raw_llhh.inputform.txt','_raw_llhx.inputform.txt',
                   '_raw_llxx.inputform.txt'):
        p=out/(stem+suffix)
        require(p.is_file() and p.stat().st_size>0,
                f'{model}: missing/empty Matchete EFT2 output '+str(p))
    matched=(out/(stem+'_sm_tree_plus_loop.inputform.txt')).read_text(encoding='utf8')
    require('Field[NewScalarX' in matched and 'Field[NewScalarY' not in matched and
            'Field[NewFermion' not in matched,
            f'{model}: generated matching did not retain X while removing Y/F')
    state.update({
        'status':'RawYOnlyMatcheteOutputsExtractedNotPhysicalMatching',
        'result':relative(dest),
        'report_sha256':sha(dest),
        'raw_operator_sector_counts':counts,
        'raw_operator_sector_supported_by_scope':{
            'LLHH':True,'LLHX':scope=='full','LLXX':scope=='full'},
        'outdir':relative(out),
        'full_eft_checked':scope=='full',
        'symbolic_operator_analysis_performed':True,
    })
    certificate=out/(stem+'_driver_certificate.json')
    certificate.write_text(json.dumps(state,indent=2)+'\n',encoding='utf8')
    print(f'{model}: native Y-only Matchete raw sectors {counts}; physical matching NOT CERTIFIED',flush=True)
    return state


def main()->int:
    ap=argparse.ArgumentParser(description=__doc__)
    g=ap.add_mutually_exclusive_group(required=True)
    g.add_argument('--all',action='store_true')
    g.add_argument('--model',choices=MODELS)
    ap.add_argument('--scope',choices=('full','leading-c5'),default='full')
    ap.add_argument('--preflight-only',action='store_true')
    ap.add_argument('--output-dir',type=Path,default=OUT)
    ap.add_argument('--wolframscript',default='wolframscript')
    ap.add_argument('--timeout-sec',type=int,default=7200)
    args=ap.parse_args()
    models=MODELS if args.all else (args.model,)
    for model in models:
        run_model(model,out=args.output_dir.resolve(),scope=args.scope,
                  preview=args.preflight_only,wolframscript=args.wolframscript,
                  timeout_sec=args.timeout_sec)
    return 0

if __name__=='__main__':raise SystemExit(main())
