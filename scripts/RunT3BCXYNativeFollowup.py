"""Fail-closed experimental native Y->X Matchete structural matching.

This is NOT a physical two-threshold matching result. The native Y-only source
is independently reconstructed rather than EOMEFT-equivalent to F-first
matching. The X matching starts with Y's TREE-ONLY EOMEFT, deliberately
excluding inherited one-loop Y matching. It uses original non-running
prepared input values, not the continuous_v1 running threshold boundaries.

python -m scripts.RunT3BCXYNativeFollowup --all --check-x
python -m scripts.RunT3BCXYNativeFollowup --all --match-x
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import subprocess

ROOT=Path(__file__).resolve().parents[1]
WL=ROOT/'Lagrangian/RunT3BCXYNativeFollowup.experimental.wl'
PARAM=ROOT/'output/general_t3_mass_eigen_matchete'
TREE=PARAM/'eft1_fermion_first_tree'
OUT=PARAM/'eft1_scalar_xy_native_followup_v2'
MODELS=('B_-1','C_-1')
Y_STATUS='YOnlyIndependentCanonicalRawEFT2NotCertified'
CHECK_STATUS='XFollowupEFT2XHeavyCheckPassedOnly'
MATCH_STATUS='XFollowupYTreeOnlyXRawMatchingNotCertified'
FORMAT='t3-bc-xy-native-followup-driver-v2'


X_FLAGS=frozenset(('ordinary_check','hermiticity','canonical',
                     'heavy_mass_basis','gauge_invariant','all_objects_defined'))


def full_native_check(check:object)->bool:
    """Require exactly the six full Matchete validation booleans."""
    return isinstance(check,dict) and set(check)==X_FLAGS and \
        all(value is True for value in check.values())


def certified_x_input_selection(report:dict)->bool:
    """Strict *structural* X-only input gate; not a physical certificate.

    Unchanged EOMEFT is eligible only with direct X-heavy Matchete checks
    and an exactly zero light-versus-heavy free-Lagrangian difference.
    The rebased expression must independently pass all six Matchete checks.
    Neither branch silently uses MatchedEFT.
    """
    selection=report.get('matching_input_selection')
    if report.get('selected_x_input_full_matchete_check_passed') is not True:
        return False
    if selection=='heavy_free_rebased_y_tree_eomeft':
        return full_native_check(report.get('rebased_x_heavy_eft2_detailed_check'))
    if selection=='unchanged_y_tree_eomeft_proven_equal_free_lagrangians':
        return full_native_check(report.get('direct_x_heavy_eft2_detailed_check')) and \
            report.get('free_lagrangians_algebraically_identical') is True
    return False


def require(condition:bool,message:str)->None:
    if not condition:
        raise RuntimeError('Native XY sequential diagnostic: '+message)


def sha(path:Path)->str:
    return sha256(path.read_bytes()).hexdigest()


def input_paths(model:str)->dict[str,Path]:
    require(model in MODELS,'unsupported model')
    tag=model.replace('+','p').replace('-','m')
    return {'script':WL,'parameters':PARAM/(tag+'_final_parameters.json'),
            'first_fermion_eomeft':TREE/(
              tag+'_final_eft1_fermion_first_tree_eft1_transition_eom.inputform.txt')}


def preflight(model:str)->dict[str,str]:
    paths=input_paths(model)
    for name,path in paths.items():
        require(path.is_file(),f'{model}: missing {name}: {path}')
    source=paths['script'].read_text(encoding='utf-8')
    for token in ('T3BCXFollowup[treeResult,loopResult,outPath]',
                  'T3_BC_X_FOLLOWUP_MATCH',
                  '"saved_fermion_EOMEFT_equivalence_certified"->False',
                  '"y_loop_inherited_in_x_matching"->False',
                  '"running_threshold_couplings_applied"->False'):
        require(token in source,f'Missing X-followup source invariant: {token}')
    tree=paths['first_fermion_eomeft'].read_text(encoding='utf-8')
    require(len(tree)>1000 and all(x in tree for x in
         ('Field[NewScalarX','Field[NewScalarY','GammaCC')) and
         'Field[NewFermion' not in tree,f'{model}: changed F-first EOMEFT inventory')
    params=json.loads(paths['parameters'].read_text(encoding='utf-8'))
    require(params.get('model')==model,f'{model}: parameter model mismatch')
    return {key:sha(path) for key,path in paths.items()}


def run_model(model:str, *, mode:str, output_dir:Path=OUT,
              wolframscript:str='wolframscript',timeout:int=7200,
              overwrite:bool=False)->dict:
    require(mode in ('check','match','diagnose'),'unexpected mode')
    hashes=preflight(model)
    tag=model.replace('+','p').replace('-','m')
    stem=tag+'_final_xy_native_'+mode+'_x'
    output_dir.mkdir(parents=True,exist_ok=True)
    target=output_dir/(stem+'.json')
    detail=output_dir/(stem+'_x_followup.json')
    log=output_dir/(stem+'.log')
    provenance=output_dir/(stem+'.provenance.json')
    if not overwrite:
        require(not any(p.exists() for p in (target,detail,log,provenance)),
                f'{model}: prior output exists; choose new --output-dir or --overwrite')
    if mode=='match':
        # A successful diagnostic has the SAME strict native X-heavy checks,
        # so rerunning expensive Y-only matching just to create a certificate
        # is unnecessary. Source hashes must still be identical.
        candidates=[output_dir/(tag+f'_final_xy_native_{m}_x.provenance.json')
                    for m in ('check','diagnose')]
        certs=[json.loads(p.read_text(encoding='utf-8')) for p in candidates if p.is_file()]
        require(any(c.get('status')==CHECK_STATUS and
                    c.get('source_sha256')==hashes and
                    c.get('all_x_checks_passed') is True and
                    c.get('physical_c5_certified') is False for c in certs),
                f'{model}: X-heavy check certificate absent, stale or failed; '
                f'run --diagnose-x (or --check-x) first')
    paths=input_paths(model)
    env=dict(os.environ)
    env.pop('T3_EFT1_LEADING_C5_ONLY',None)
    env.pop('T3_EXPERIMENTAL_GENERAL_QUARTICS',None)
    env['T3_BC_Y_INDEPENDENT_CANONICAL']='1'
    env['T3_BC_Y_INDEPENDENT_CHECK_ONLY']='0'
    env['T3_BC_X_FOLLOWUP_MATCH']='1' if mode=='match' else '0'
    command=[wolframscript,'-file',str(WL),str(paths['parameters']),
             str(paths['first_fermion_eomeft']),str(target)]
    try:
        proc=subprocess.run(command,env=env,cwd=ROOT,encoding='utf-8',
                            errors='replace',text=True,capture_output=True,
                            check=False,timeout=timeout)
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError(f'{model}: Matchete timed out at X followup') from exc
    log.write_text('command: '+' '.join(command)+'\nexit='+str(proc.returncode)+
                   '\nstdout:\n'+proc.stdout+'\nstderr:\n'+proc.stderr,
                   encoding='utf-8')
    yreport=json.loads(target.read_text(encoding='utf-8')) if target.is_file() else {}
    xreport=json.loads(detail.read_text(encoding='utf-8')) if detail.is_file() else {}
    expected=CHECK_STATUS if mode!='match' else MATCH_STATUS
    checks=xreport.get('rebased_x_heavy_eft2_detailed_check',{})
    direct=xreport.get('direct_x_heavy_eft2_detailed_check',{})
    checks_ok=certified_x_input_selection(xreport)
    audit_file=output_dir/(stem+'_x_followup_rebase_crosschecks.json')
    audit=json.loads(audit_file.read_text(encoding='utf-8')) if audit_file.is_file() else {}
    result={'format':FORMAT,'model':model,'status':xreport.get('status','NO_X_REPORT'),
      'mode':mode,'source_sha256':hashes,'returncode':proc.returncode,
      'Y_only_native_matching_status':yreport.get('status','NO_Y_REPORT'),
      'all_x_checks_passed':checks_ok,
      'rebased_matchete_checks':checks,'direct_heavy_matchete_checks':direct,
      'free_lagrangians_algebraically_identical':xreport.get('free_lagrangians_algebraically_identical'),
      'matching_input_selection':xreport.get('matching_input_selection'),
      'rebase_audit_file':audit_file.name if audit else None,
      'raw_x_one_loop_matching_completed':mode=='match' and
         xreport.get('status')==MATCH_STATUS,
      'saved_fermion_eomeft_equivalence_proved':False,
      'inherited_y_one_loop_in_x_source':False,
      'running_boundary_substituted':False,
      'renormalized_sequential_matching_certified':False,
      'physical_c5_certified':False,
      'report_file':target.name,'x_followup_report_file':detail.name,
      'log_file':log.name}
    provenance.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    if mode=='diagnose':
        require(yreport.get('status')==Y_STATUS and xreport and audit,
                f'{model}: Y native stage or X diagnostic missing; see {log}')
        print(f'{model}: X heavy check comparison: '
              f'direct={audit.get("direct_heavy_check_passed")}, '
              f'rebased={audit.get("rebased_heavy_check_passed")}, '
              f'free_equal={audit.get("free_lagrangians_algebraically_identical")}; '
              f'status={xreport.get("status")}; physical matching=NOT CERTIFIED',flush=True)
        return result
    require(proc.returncode==0 and yreport.get('status')==Y_STATUS and
        xreport.get('status')==expected and checks_ok and
        yreport.get('physical_matching_performed') is False and
        xreport.get('physical_matching_certified') is False,
        f'{model}: X continuation did not pass; status={xreport.get("status")}, '
        f'Y status={yreport.get("status")}; '
        f'direct heavy checks={direct}; rebased checks={checks}; '
        f'free_equal={xreport.get("free_lagrangians_algebraically_identical")}; '
        f'see {log} and {detail}')
    print(f'{model}: X followup {mode} PASS; status={expected}; '
          'physical matching=NOT CERTIFIED',flush=True)
    return result


def main()->int:
    p=argparse.ArgumentParser(description=__doc__)
    models=p.add_mutually_exclusive_group(required=True)
    models.add_argument('--all',action='store_true')
    models.add_argument('--model',choices=MODELS)
    actions=p.add_mutually_exclusive_group(required=True)
    actions.add_argument('--preflight-only',action='store_true')
    actions.add_argument('--check-x',action='store_true')
    actions.add_argument('--diagnose-x',action='store_true',help='run both models and export actual heavy-registry/rebase checks even when X check fails')
    actions.add_argument('--match-x',action='store_true')
    p.add_argument('--output-dir',type=Path,default=OUT)
    p.add_argument('--wolframscript',default='wolframscript')
    p.add_argument('--timeout-sec',type=int,default=7200)
    p.add_argument('--overwrite',action='store_true')
    args=p.parse_args()
    targets=MODELS if args.all else (args.model,)
    for model in targets:
        if args.preflight_only:
            print(f'{model}: preflight PASS, source hashes={preflight(model)}',flush=True)
        else:
            run_model(model,mode='match' if args.match_x else ('diagnose' if args.diagnose_x else 'check'),
                      output_dir=args.output_dir.resolve(),
                      wolframscript=args.wolframscript,
                      timeout=args.timeout_sec,overwrite=args.overwrite)
    return 0

if __name__=='__main__':
    raise SystemExit(main())
