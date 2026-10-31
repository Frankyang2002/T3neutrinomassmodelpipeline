#!/usr/bin/env python3
"""Guarded six-model native Matchete extended-C5 extraction batch.

The existing 16-model UV RGBeta ODEs are NOT rerun or overwritten.
The Matchete experiment matches all heavy fields at one scale in a restricted
basis (general quartics, extra Yukawas, but NOT extra scalar bilinears).
This is not yet hierarchical neutrino-mass matching or a RGBeta coordinate map.

Usage (repository root):
  python scripts/RunT3ExtendedMatchingBatch.py --self-test
  python scripts/RunT3ExtendedMatchingBatch.py --preflight
  python scripts/RunT3ExtendedMatchingBatch.py --pilot
  python scripts/RunT3ExtendedMatchingBatch.py --all
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from hashlib import sha256
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / 'Reports/output/workstream_a/extended_c5_native'
NATIVE = ROOT / 'Lagrangian/RunT3ExtendedC5Native.wl'
COUPLED_SUMMARY = ROOT / 'Reports/output/workstream_a/full_coupled_general_yukawa_uv/summary.json'
EXPORT_SUMMARY = ROOT / 'Reports/output/workstream_a/integrated_general_yukawa_uv/summary.json'
SPECS = {
    'A_p0': ('A','p0','y1Conj'), 'A_m2': ('A','m2','y2Conj'),
    'D_p0': ('D','p0','y1Conj'), 'D_m2': ('D','m2','y2Conj'),
    'E_p0': ('E','p0','y1Conj'), 'E_m2': ('E','m2','y2Conj'),
}
SOURCES = (
    'Lagrangian/model/LagrangianBuilder.wl',
    'Lagrangian/interactions/T3Topology.wl',
    'Lagrangian/RunMatching.wl',
    'Lagrangian/RunT3ExtendedC5Native.wl',
    'RGE/running/rgbeta/T3RGBetaModel.wl',
)
SUCCESS = 'NativeOneStepExtendedWeinbergExtracted'


def need(ok: bool, msg: str):
    if not ok: raise ValueError(msg)


def read(path: Path) -> dict:
    obj = json.loads(path.read_text(encoding='utf-8-sig'))
    need(isinstance(obj,dict),str(path)+': expected JSON object')
    return obj


def atomic(path:Path, content:dict):
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_name(path.name+'.tmp')
    tmp.write_text(json.dumps(content,indent=2)+'\n',encoding='utf-8')
    tmp.replace(path)


def fingerprints(root:Path=ROOT)->dict:
    result={}
    for p in SOURCES:
        path=root/p
        need(path.is_file(),'Missing source '+str(path))
        result[p]=sha256(path.read_bytes()).hexdigest()
    model=(root/'Lagrangian/model/LagrangianBuilder.wl').read_text(encoding='utf-8-sig')
    native=(root/'Lagrangian/RunT3ExtendedC5Native.wl').read_text(encoding='utf-8-sig')
    extra=(root/'Lagrangian/interactions/T3Topology.wl').read_text(encoding='utf-8-sig')
    need('BuildT3Lagrangian[model_Association]' in model,
         'Matchete builder source contract changed')
    need('T3GeneralExtraYukawaCount[model_Association]' in extra and
         'BuildT3ExtraYukawaCandidates' in extra,
         'Matchete extra-Yukawa source contract changed')
    for frag in ('RunT3Matching[built["LUV"],5,1]', 'ExtractWeinbergCoefficient',
                 'T3_EXPERIMENTAL_GENERAL_QUADRATICS','native_coefficient_extra_difference_path'):
        need(frag in native,'Native matching script source contract missing '+frag)
    return result


def previous_runs(root:Path=ROOT)->dict:
    coupled=read(root/'Reports/output/workstream_a/full_coupled_general_yukawa_uv/summary.json')
    need(coupled.get('status')=='AllRequestedFullCoupledUVRGBetaOnlyPassed',
         '16-model coupled UV summary not successful')
    rows=coupled.get('rows',[])
    need(len(rows)==16 and all(r.get('status')=='CoupledGeneralYukawaRGBetaOnlyUVPassed'
                               for r in rows),
         'Coupled UV row status not verified for 16 models')
    exported=read(root/'Reports/output/workstream_a/integrated_general_yukawa_uv/summary.json')
    need(exported.get('all_requested_models_succeeded') is True,
         '16-model general RGBeta export summary not successful')
    for model in SPECS:
        record=root/'Reports/output/workstream_a/full_coupled_general_yukawa_uv/models'/f'{model}.json'
        obj=read(record)
        need(obj.get('status')=='CoupledGeneralYukawaRGBetaOnlyUVPassed' and
             obj.get('extra_yukawa')==SPECS[model][2] and
             obj.get('endpoint_above_evolved_heavy_mass') is True,
             model+': missing successful prior coupled RG trajectory')
    return {'status':'Previous16RGESucceededSixNativeTargetsAvailable',
            'coupled_status':coupled['status'],
            'full_coupled_summary_sha256':sha256((root/'Reports/output/workstream_a/full_coupled_general_yukawa_uv/summary.json').read_bytes()).hexdigest(),
            'general_export_sha256':sha256((root/'Reports/output/workstream_a/integrated_general_yukawa_uv/summary.json').read_bytes()).hexdigest()}


def validate_native(model:str, path:Path, *, heavy_flavours:int) -> dict:
    need(model in SPECS, 'Unknown target model')
    obj=read(path)
    need(obj.get('format')=='t3-extended-c5-native-matchete-one-step-v1' and
         obj.get('model')==model and obj.get('heavy_flavours')==heavy_flavours,
         model+': native report format/model/flavors mismatch')
    need(obj.get('status')==SUCCESS,model+': native status '+str(obj.get('status'))+
         ': '+str(obj.get('reason','unavailable')))
    need(obj.get('additional_yukawa')==SPECS[model][2],
         model+': incorrect extra Yukawa recorded')
    for key in ('one_loop_matching_performed', 'native_matchete_lagrangian_check_passed',
                'matched_c5_in_matchete_native_coordinates'):
        need(obj.get(key) is True, model+': native stage did not establish '+key)
    for key in ('absolute_matchete_rgbeta_dictionary_certified',
                'hierarchical_matching_certified','cross_package_numeric_substitution_authorized',
                'physical_neutrino_mass_prediction'):
        need(obj.get(key) is False, model+': native report overclaims '+key)
    evidence={}
    for key in ('native_coefficient_full_path','native_coefficient_zero_extra_path',
                'native_coefficient_extra_difference_path'):
        val=obj.get(key)
        need(isinstance(val,str) and val.strip(),model+': missing native coefficient path')
        filename=Path(val)
        need(filename.is_file(),model+': absent coefficient file '+str(filename))
        txt=filename.read_text(encoding='utf-8-sig').strip()
        need(bool(txt),model+': empty '+key)
        evidence[key]={'path':str(filename),'sha256':sha256(filename.read_bytes()).hexdigest(),
                       'characters':len(txt)}
    return {'status':SUCCESS, 'model':model,'heavy_flavours':heavy_flavours,
            'extra_sector':obj.get('extra_dependent_sector_diagnostic'),
            'coefficient_evidence':evidence,
            'status_scope':'native Matchete one-step restricted-basis coefficient only',
            'absolute_matchete_rgbeta_dictionary_certified':False,
            'hierarchical_matching_certified':False,
            'numeric_3x3_flavor_lift_certified':False}


def tail(log:str,lines:int=16)->str:
    return '\n'.join(log.splitlines()[-lines:])


def run_one(model:str,*,heavy_flavours:int,wl:str,timeout:int,root:Path=ROOT,
            output:Path=OUTPUT)->dict:
    cl,tag,extra=SPECS[model]
    native_output=output/'native'/f'{model}_h{heavy_flavours}.json'
    log=output/'logs'/f'{model}_h{heavy_flavours}.log'
    native_output.parent.mkdir(parents=True,exist_ok=True)
    log.parent.mkdir(parents=True,exist_ok=True)
    if native_output.exists():
        native_output.unlink()  # prevent stale success from a prior failed run
    if log.exists(): log.unlink()
    cmd=[wl,'-file',str(root/'Lagrangian/RunT3ExtendedC5Native.wl'),
         str(native_output),cl,tag,str(heavy_flavours)]
    env=dict(os.environ)
    env.update({'T3_EXPERIMENTAL_GENERAL_FERMIONS':'1',
                'T3_EXPERIMENTAL_GENERAL_QUARTICS':'1',
                'T3_EXPERIMENTAL_FLAVOR_MASSES':'1',
                'T3_EXPERIMENTAL_GENERAL_QUADRATICS':'0'})
    try:
        process=subprocess.run(cmd,cwd=root,env=env,stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT,text=True,errors='replace',
                               timeout=timeout,check=False)
        stdout=process.stdout or ''
        return_code=process.returncode
    except subprocess.TimeoutExpired as exc:
        stdout='TIMEOUT\n'+str(exc)+'\n'+(exc.stdout.decode('utf-8',errors='replace')
                        if isinstance(exc.stdout,bytes) else (exc.stdout or ''))
        return_code=124
    except OSError as exc:
        stdout=type(exc).__name__+': '+str(exc)
        return_code=127
    log.write_text(stdout,encoding='utf-8')
    if return_code!=0:
        result={'status':'Blocked','model':model,'exit_code':return_code,
                'reason':'Native Wolfram stage failed; inspect final log lines',
                'log_path':str(log),'log_tail':tail(stdout)}
        if native_output.exists():
            try:
                native=read(native_output)
                result.update(native_status=native.get('status'),
                              failed_step=native.get('failed_step'),
                              native_reason=native.get('reason'))
            except Exception as exc:result['native_parse_error']=str(exc)
        return result
    try:
        result=validate_native(model,native_output,heavy_flavours=heavy_flavours)
        result.update(exit_code=0, native_path=str(native_output),log_path=str(log))
        return result
    except Exception as exc:
        return {'status':'Blocked','model':model,'exit_code':0,
                'reason':'Invalid successful native report: '+str(exc),
                'log_path':str(log),'log_tail':tail(stdout)}


def self_test():
    need(len(SPECS)==6 and sum(v[2]=='y1Conj' for v in SPECS.values())==3 and
         sum(v[2]=='y2Conj' for v in SPECS.values())==3,
         'Six-model registry mismatch')
    with tempfile.TemporaryDirectory() as td:
        t=Path(td)
        nat=t/'native.json'
        content={'format':'t3-extended-c5-native-matchete-one-step-v1',
                 'model':'A_p0','heavy_flavours':1,'status':SUCCESS,
                 'additional_yukawa':'y1Conj',
                 'one_loop_matching_performed':True,
                 'native_matchete_lagrangian_check_passed':True,
                 'matched_c5_in_matchete_native_coordinates':True,
                 'absolute_matchete_rgbeta_dictionary_certified':False,
                 'hierarchical_matching_certified':False,
                 'cross_package_numeric_substitution_authorized':False,
                 'physical_neutrino_mass_prediction':False}
        for key in ('native_coefficient_full_path','native_coefficient_zero_extra_path',
                    'native_coefficient_extra_difference_path'):
            path=t/(key+'.txt');path.write_text('1 + 2 y1Conj',encoding='utf-8')
            content[key]=str(path)
        atomic(nat,content)
        valid=validate_native('A_p0',nat,heavy_flavours=1)
        need(valid['status']==SUCCESS,'Synthetic native witness rejected')
        content['cross_package_numeric_substitution_authorized']=True
        atomic(nat,content)
        try:validate_native('A_p0',nat,heavy_flavours=1)
        except ValueError:pass
        else:raise AssertionError('Overclaimed physical mapping accepted')
        content['cross_package_numeric_substitution_authorized']=False
        content['model']='B_m1';atomic(nat,content)
        try:validate_native('A_p0',nat,heavy_flavours=1)
        except ValueError:pass
        else:raise AssertionError('Wrong model accepted')
        need(tail('1\n2\n3',2)=='2\n3','Log tail error')
    print('PASS: six-model source registry, native evidence schema, no-physical-overclaim controls')


def main()->int:
    ap=argparse.ArgumentParser(description=__doc__)
    g=ap.add_mutually_exclusive_group(required=True)
    g.add_argument('--self-test',action='store_true')
    g.add_argument('--preflight',action='store_true')
    g.add_argument('--pilot',action='store_true')
    g.add_argument('--all',action='store_true')
    g.add_argument('--proxy-only',action='store_true')
    g.add_argument('--model',choices=list(SPECS))
    ap.add_argument('--heavy-flavours',type=int,choices=(1,3),default=1)
    ap.add_argument('--wolframscript',default='wolframscript')
    ap.add_argument('--timeout-seconds',type=int,default=3600)
    args=ap.parse_args()
    if args.self_test:
        self_test();return 0
    try:
        src=fingerprints()
        prev=previous_runs()
        preflight={'format':'t3-extended-c5-preflight-v1',
                   'status':'NativeMatchingReadyNotYetExecuted',
                   'previous_runs':prev,'source_sha256':src,
                   'native_matchete_rgbeta_dictionary_certified':False,
                   'physical_extended_neutrino_matching_complete':False}
        atomic(OUTPUT/'preflight.json',preflight)
        print('PASS: 16-model UV source provenance and six matched-model inputs verified',flush=True)
        if args.preflight:return 0
        models=([args.model] if args.model else ['A_p0'] if args.pilot else list(SPECS))
        rows=[]
        from Numerical.experimental.T3ExtendedC5KernelProxy import (
            run_one as run_proxy, atomic as save_proxy, OUTPUT as proxy_output)
        run_native = not args.proxy_only
        for n,model in enumerate(models,1):
            print(f'[{n}/{len(models)}] Matching {model} with native Matchete ({args.heavy_flavours} heavy flavours)...',flush=True)
            try:
                proxy=run_proxy(model)
                save_proxy(proxy_output/'models'/f'{model}.json', proxy)
                proxy_status={'status':proxy['status'],
                    'delta_relative':proxy['proxy_delta_relative_to_control'],
                    'report_path':str(proxy_output/'models'/f'{model}.json')}
                print(f'[{n}/{len(models)}] Numerical proxy PASS {model}: '
                      f'{proxy["proxy_delta_relative_to_control"]:.4e}',flush=True)
            except Exception as exc:
                proxy_status={'status':'Blocked','reason':type(exc).__name__+': '+str(exc)}
                print(f'[{n}/{len(models)}] Numerical proxy BLOCKED {model}: {exc}',flush=True)
            if run_native:
                record=run_one(model,heavy_flavours=args.heavy_flavours,wl=args.wolframscript,
                               timeout=args.timeout_seconds)
            else:
                record={'status':'NotRequested','model':model}
            record['numerical_proxy']=proxy_status
            rows.append(record)
            print(f'[{n}/{len(models)}] {record["status"]} {model}',flush=True)
            if record['status']=='Blocked':
                print(record.get('reason',''))
                print(record.get('log_tail',''))
            report={'format':'t3-extended-c5-native-one-step-batch-v1',
                    'status':'RunningOrPartial', 'models_requested':models,
                    'heavy_flavours':args.heavy_flavours,
                    'source_sha256':src,'rows':rows,
                    'physical_matchete_rgbeta_dictionary_certified':False,
                    'hierarchical_extended_neutrino_matching_complete':False,
                    'extended_neutrino_observables_computed':False}
            atomic(OUTPUT/'summary.json',report)
        proxy_ok=len(rows)==len(models) and all(
            r.get('numerical_proxy',{}).get('status')=='OriginalT3YukawaFlavorTensorSensitivityOnly'
            for r in rows)
        complete=proxy_ok and (not run_native or all(r['status']==SUCCESS for r in rows))
        report['status']=(('AllRequestedOneStepNativeC5Extracted' if run_native else
                           'AllRequestedNumericalProxiesCalculated') if complete else
                          'OneStepNativeC5PartialOrBlocked')
        report['generated_at_utc']=datetime.now(timezone.utc).isoformat()
        atomic(OUTPUT/'summary.json',report)
        print('SUMMARY:',report['status'],len([x for x in rows if x['status']==SUCCESS]),'/',len(models))
        print('Report:',OUTPUT/'summary.json')
        print('Original-flavor-tensor numerical proxy: ',sum(
            r.get('numerical_proxy',{}).get('status')=='OriginalT3YukawaFlavorTensorSensitivityOnly'
            for r in rows),'/',len(models))
        print('Physical cross-package coupling dictionary: NOT CERTIFIED')
        print('Extended hierarchical matching / physical neutrino masses: NOT CERTIFIED')
        return 0 if complete else 4
    except Exception as exc:
        print('BLOCKED:',type(exc).__name__+':',exc,file=sys.stderr)
        return 4


if __name__=='__main__':
    raise SystemExit(main())
