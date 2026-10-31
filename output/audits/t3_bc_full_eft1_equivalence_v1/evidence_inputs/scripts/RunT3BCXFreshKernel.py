"""Attempt native X-only matching in a fresh Wolfram kernel from a saved Y-tree EOMEFT.

No Y rematching in this invocation, no inherited Y finite matching in X
source, no RG-coupling substitution, and no physical-matching certification.
Both models are processed independently even if one fails. All prior source
hashes and diagnostics must match the saved v2 preflight.

  python -u -X utf8 -m scripts.RunT3BCXFreshKernel --all --preflight-only
  python -u -X utf8 -m scripts.RunT3BCXFreshKernel --all --check-x
  python -u -X utf8 -m scripts.RunT3BCXFreshKernel --all --match-x
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import subprocess

ROOT=Path(__file__).resolve().parents[1]
ROOT_MAT=ROOT/'output/general_t3_mass_eigen_matchete'
ORIGINAL_WL=ROOT/'Lagrangian/RunT3BCXYNativeFollowup.experimental.wl'
FRESH_WL=ROOT/'Lagrangian/RunT3BCXFromSavedY.experimental.wl'
V2_DIR=ROOT_MAT/'eft1_scalar_xy_native_followup_v2'
OUT=ROOT_MAT/'eft1_scalar_xy_fresh_kernel_v3'
F_TREE=ROOT_MAT/'eft1_fermion_first_tree'
MODELS=('B_-1','C_-1')
V2_CHECK='XFollowupEFT2XHeavyCheckPassedOnly'
Y_STATUS='YOnlyIndependentCanonicalRawEFT2NotCertified'
MATCH_STATUS='XFollowupYTreeOnlyXRawMatchingNotCertified'
CHECK_STATUS='XFollowupEFT2XHeavyCheckPassedOnly'
FRESH_Y_STATUS='SavedYTreeEOMConsumedForFreshXNoYRematch'
FORMAT='t3-bc-fresh-x-only-native-attempt-v3'
FLAGS={'ordinary_check','hermiticity','canonical','heavy_mass_basis','gauge_invariant','all_objects_defined'}


def sha(path:Path)->str:
    return sha256(path.read_bytes()).hexdigest()


def require(value:bool,why:str)->None:
    if not value:raise RuntimeError('Fresh X-only native: '+why)


def parse(path:Path)->dict:
    require(path.is_file(),f'missing {path}')
    obj=json.loads(path.read_text(encoding='utf-8'))
    require(isinstance(obj,dict),f'invalid JSON object: {path}')
    return obj


def passed(check:object)->bool:
    return isinstance(check,dict) and set(check)==FLAGS and all(v is True for v in check.values())


def model_paths(model:str)->dict[str,Path]:
    require(model in MODELS,'invalid model')
    tag=model.replace('+','p').replace('-','m')
    stem=f'{tag}_final_xy_native_diagnose_x'
    xbase=f'{stem}_x_followup'
    return {
        'original_wl':ORIGINAL_WL,
        'fresh_wl':FRESH_WL,
        'parameters':ROOT_MAT/f'{tag}_final_parameters.json',
        'fermion_eom':F_TREE/f'{tag}_final_eft1_fermion_first_tree_eft1_transition_eom.inputform.txt',
        'prior_provenance':V2_DIR/(stem+'.provenance.json'),
        'prior_y_report':V2_DIR/(stem+'.json'),
        'prior_x_report':V2_DIR/(xbase+'.json'),
        'prior_y_tree':V2_DIR/(xbase+'_y_tree_eomeft.inputform.txt'),
        'prior_y_full':V2_DIR/(xbase+'_y_full_eomeft.inputform.txt'),
    }


def preflight(model:str)->dict:
    p=model_paths(model)
    for name,path in p.items():require(path.is_file(),f'{model}: missing {name}: {path}')
    v2=parse(p['prior_provenance']); x=parse(p['prior_x_report']); y=parse(p['prior_y_report'])
    require(v2.get('model')==model and v2.get('status')==V2_CHECK and
            v2.get('all_x_checks_passed') is True and
            v2.get('physical_c5_certified') is False,
            f'{model}: v2 native structural X validation not certified')
    require(y.get('status')==Y_STATUS and y.get('model')==model and
            y.get('physical_matching_performed') is False,
            f'{model}: v2 Y-only Matchete evidence invalid')
    require(x.get('status')==V2_CHECK and passed(x.get('rebased_x_heavy_eft2_detailed_check')) and
            x.get('selected_x_input_full_matchete_check_passed') is True and
            x.get('matching_input_selection')=='heavy_free_rebased_y_tree_eomeft' and
            x.get('physical_matching_certified') is False,
            f'{model}: v2 rebase check not independently passed')
    old=v2.get('source_sha256')
    actual={'script':sha(p['original_wl']), 'parameters':sha(p['parameters']),
            'first_fermion_eomeft':sha(p['fermion_eom'])}
    require(old==actual,f'{model}: v2 source hashes changed; cannot reuse EOMEFT ({old} vs {actual})')
    tree=p['prior_y_tree'].read_text(encoding='utf-8')
    full=p['prior_y_full'].read_text(encoding='utf-8')
    for name,data in (('tree',tree),('full',full)):
        require(len(data)>500 and 'Field[NewScalarX' in data and
                'Field[NewScalarY' not in data and 'Field[NewFermion' not in data,
                f'{model}: saved Y {name} EOMEFT field inventory invalid')
    source=p['fresh_wl'].read_text(encoding='utf-8')
    for token in ('T3_BC_X_FRESH_EOM_TREE','XFollowupTreeMatchingInProgress',
                  'XFollowupOneLoopMatchingInProgress',
                  'SavedYTreeEOMConsumedForFreshXNoYRematch'):
        require(token in source,f'{model}: changed fresh Wolfram source {token}')
    require(parse(p['parameters']).get('model')==model,
            f'{model}: wrong prepared mass-basis model')
    return {'model':model, 'source_sha256':actual,
            'fresh_wl_sha256':sha(p['fresh_wl']),
            'saved_y_tree_eomeft_sha256':sha(p['prior_y_tree']),
            'saved_y_full_eomeft_sha256':sha(p['prior_y_full']),
            'v2_validation_status':V2_CHECK,
            'physical_matching_certified':False}


def classify(y:dict,x:dict,rc:int)->tuple[str,bool]:
    status=x.get('status','X_FOLLOWUP_REPORT_MISSING')
    good=rc==0 and y.get('status')==FRESH_Y_STATUS and \
         status==MATCH_STATUS and passed(x.get('rebased_x_heavy_eft2_detailed_check')) and \
         x.get('matching_input_selection')=='heavy_free_rebased_y_tree_eomeft' and \
         x.get('physical_matching_certified') is False
    return status,good


def run_one(model:str,*,mode:str,out:Path=OUT,wolframscript:str='wolframscript',
            timeout:int=7200,overwrite:bool=False)->dict:
    require(mode in ('check','match'),f'bad mode {mode}')
    prov=preflight(model)
    paths=model_paths(model)
    tag=model.replace('+','p').replace('-','m')
    stem=f'{tag}_final_x_from_saved_y_{mode}'
    out.mkdir(parents=True,exist_ok=True)
    target=out/(stem+'.json')
    xreport=out/(stem+'_x_followup.json')
    progress=out/(stem+'.progress.json')
    report_path=out/(stem+'.provenance.json')
    log=out/(stem+'.log')
    existing=[target,xreport,progress,report_path,log]
    require(overwrite or not any(p.exists() for p in existing),
            f'{model}: previous output exists, use --overwrite or new --output-dir')
    if overwrite:
        for p in out.glob(stem+'*'):
            if p.is_file():p.unlink()
    env=dict(os.environ)
    env.pop('T3_EFT1_LEADING_C5_ONLY',None)
    env.pop('T3_EXPERIMENTAL_GENERAL_QUARTICS',None)
    env['T3_BC_Y_INDEPENDENT_CANONICAL']='1'
    env['T3_BC_Y_INDEPENDENT_CHECK_ONLY']='0'
    env['T3_BC_X_FOLLOWUP_MATCH']='1' if mode=='match' else '0'
    env['T3_BC_X_FRESH_EOM_TREE']=str(paths['prior_y_tree'].resolve())
    env['T3_BC_X_FRESH_EOM_FULL']=str(paths['prior_y_full'].resolve())
    cmd=[wolframscript,'-file',str(paths['fresh_wl']),str(paths['parameters']),
         str(paths['fermion_eom']),str(target)]
    try:
        proc=subprocess.run(cmd,cwd=ROOT,env=env,timeout=timeout,
                            capture_output=True,text=True,encoding='utf-8',
                            errors='replace',check=False)
        rc=proc.returncode
        stdout=proc.stdout;stderr=proc.stderr
    except subprocess.TimeoutExpired as exc:
        rc=124;stdout=(exc.stdout or b'');stderr=(exc.stderr or b'')
        if isinstance(stdout,bytes):stdout=stdout.decode('utf-8',errors='replace')
        if isinstance(stderr,bytes):stderr=stderr.decode('utf-8',errors='replace')
        stderr+='\nPython subprocess timeout; Wolfram stage checkpoint preserved.'
    log.write_text('command: '+' '.join(cmd)+'\nexit='+str(rc)+'\nstdout:\n'+stdout+
                   '\nstderr:\n'+stderr,encoding='utf8')
    y=parse(target) if target.is_file() else {}
    x=parse(xreport) if xreport.is_file() else {}
    progress_info=parse(progress) if progress.is_file() else {}
    status,matching_completed=classify(y,x,rc)
    if mode=='check':
        structurally_passed=rc==0 and y.get('status')==FRESH_Y_STATUS and \
            status==CHECK_STATUS and passed(x.get('rebased_x_heavy_eft2_detailed_check'))
    else:structurally_passed=matching_completed
    result={**prov,'format':FORMAT,'mode':mode,
        'status':status,'exit_code':rc,
        'fresh_x_matchete_matching_completed':matching_completed,
        'structural_x_check_passed':structurally_passed,
        'checkpoint':progress_info.get('last_stage','NOT_SAVED'),
        'last_stage_history':progress_info.get('stage_history',[])[-10:],
        'y_stage_status':y.get('status','NOT_SAVED'),
        'x_matchete_tree_status':x.get('matchete_tree_stage_status'),
        'x_matchete_one_loop_status':x.get('matchete_one_loop_stage_status'),
        'x_rebased_checks':x.get('rebased_x_heavy_eft2_detailed_check'),
        'run_log':str(log),
        'x_report':str(xreport) if xreport.exists() else None,
        'reused_verified_y_tree_eomeft':True,
        'y_one_loop_inherited_in_x_source':False,
        'f_first_eomeft_equivalence_certified':False,
        'running_boundary_substituted':False,
        'renormalized_sequential_matching_certified':False,
        'physical_matching_certified':False}
    report_path.write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print(f'{model}: X-only fresh kernel {mode}: status={status}; '
          f'exit={rc}; stage={result["checkpoint"]}; '
          f'matching={"RAW PASS" if matching_completed else "NOT COMPLETED"}; '
          'physical=NOT CERTIFIED',flush=True)
    return result


def main()->int:
    p=argparse.ArgumentParser(description=__doc__)
    grp=p.add_mutually_exclusive_group(required=True)
    grp.add_argument('--all',action='store_true')
    grp.add_argument('--model',choices=MODELS)
    work=p.add_mutually_exclusive_group(required=True)
    work.add_argument('--preflight-only',action='store_true')
    work.add_argument('--check-x',action='store_true')
    work.add_argument('--match-x',action='store_true')
    p.add_argument('--output-dir',type=Path,default=OUT)
    p.add_argument('--wolframscript',default='wolframscript')
    p.add_argument('--timeout-sec',type=int,default=7200)
    p.add_argument('--overwrite',action='store_true')
    args=p.parse_args()
    failed=False
    for model in (MODELS if args.all else (args.model,)):
        try:
            if args.preflight_only:
                print(f'{model}: preflight PASS '+json.dumps(preflight(model)),flush=True)
            else:
                info=run_one(model,mode='match' if args.match_x else 'check',
                    out=args.output_dir.resolve(),wolframscript=args.wolframscript,
                    timeout=args.timeout_sec,overwrite=args.overwrite)
                if not info['structural_x_check_passed']:failed=True
        except Exception as exc:
            failed=True
            print(f'{model}: ERROR {exc}',flush=True)
    return 1 if failed else 0

if __name__=='__main__':raise SystemExit(main())
