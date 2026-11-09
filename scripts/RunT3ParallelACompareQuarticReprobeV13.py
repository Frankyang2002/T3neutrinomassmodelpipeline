"""Read-only comparison of independent native v13 quartic-response probe routes.

Diagnostic only; never certifies the blocked v10 smoke or any physics output.
"""
from __future__ import annotations
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys
from typing import Any

from scripts.RunT3ParallelAQuarticForensicsV13 import ForensicBlocked,require

FORMAT='t3-parallel-a-native-rgbeta-quartic-response-reprobe-v13'
STATUS='NativeResponseForensicsCaptured_NoPhysicsAuthorization'
ROUTES=('isolated_registry','full_registry_then_other_couplings_zero')


def classify(method:dict[str,Any], target:str, source:str)->str:
    if not isinstance(method,dict):return 'malformed_method'
    if method.get('native_execution_returned_failure') is True:
        return 'native_execution_failed'
    if method.get('refinement_failed_or_timeout') is True:
        return 'native_refinement_failed_or_timed_out'
    if method.get('native_coefficient_numeric') is not True:
        return 'native_symbolic_or_non_linear'
    val=method.get('native_coefficient_value')
    if not isinstance(val,(int,float)) or isinstance(val,bool) or abs(val)>1e100:
        return 'nonfinite_or_non_numeric'
    if method.get('linear_in_source') is not True:
        return 'native_non_linear'
    if target==source and abs(float(val)-1)>2e-8:
        return 'native_non_unit_diagonal'
    return 'native_numeric_observation_not_full_certification'


def compare(forensic:dict,reprobe:dict,forensic_sha:str)->dict:
    require(forensic.get('status')=='NativeQuarticSmokeFailureClassified_NotCertified',
            'invalid forensic status')
    require(forensic.get('physical_Wilson_evolution_authorized') is False,
            'forensic input falsely certifies physics')
    require(reprobe.get('format')==FORMAT and reprobe.get('status')==STATUS,
            'reprobe was not native forensic output')
    require(reprobe.get('forensic_input_sha256')==forensic_sha,
            'reprobe does not match selected forensic input')
    require(reprobe.get('physical_Wilson_evolution_authorized') is False,
            'reprobe falsely certifies physics')
    require(reprobe.get('installed_RGBeta_main_file_sha256')==
            forensic.get('installed_RGBeta_main_file_sha256'),
            'native RGBeta installation changed')
    expected=forensic.get('minimal_reprobe_targets')
    require(isinstance(expected,list) and 1<=len(expected)<=3,
            'invalid forensic selection')
    models=reprobe.get('models')
    require(isinstance(models,dict) and set(models)=={x['model'] for x in expected},
            'reprobe model coverage changed')
    rows=[]
    for m,record in models.items():
        probes=record.get('probes')
        targets=[x for x in expected if x['model']==m]
        require(isinstance(probes,list) and len(probes)==len(targets),
                f'{m}: missing live native re-probe')
        for row,target in zip(probes,targets):
            require(row.get('target')==target['target'] and
                    row.get('source')==target['source'],
                    f'{m}: re-probe target/source mismatch')
            modes={k:classify(row.get(k),target['target'],target['source']) for k in ROUTES}
            if modes[ROUTES[0]]!=modes[ROUTES[1]]:
                comparison='isolation_and_full_registry_differ'
            elif modes[ROUTES[0]]=='native_numeric_observation_not_full_certification':
                comparison='both_native_routes_numeric_check_old_v10_failure'
            else:
                comparison='failure_persists_or_both_routes_unresolved'
            rows.append({'model':m,'target':target['target'],'source':target['source'],
                         'old_failure_reason':target['reason'],
                         'route_classifications':modes,'comparison':comparison})
    return {'format':'t3-parallel-a-native-rgbeta-quartic-reprobe-comparison-v13',
        'status':'NativeDifferentialProbeClassified_NotCertified',
        'forensic_input_sha256':forensic_sha,
        'results':rows,
        'native_150_diagonals_certified':False,
        'native_complex_Lam_certified':False,
        'physical_Wilson_evolution_authorized':False}


def main(argv=None)->int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--forensic',type=Path,required=True)
    p.add_argument('--reprobe',type=Path,required=True)
    args=p.parse_args(argv)
    try:
        raw=args.forensic.read_bytes();fs=sha256(raw).hexdigest()
        a=json.loads(raw.decode('utf-8-sig'))
        b=json.loads(args.reprobe.read_text(encoding='utf-8-sig'))
        report=compare(a,b,fs)
        for x in report['results']:
            print(f"v13 {x['model']} {x['target']} <- {x['source']}: "
                  f"old={x['old_failure_reason']}; ")
            for route,value in x['route_classifications'].items():
                print(f'    {route}: {value}')
            print('    comparison: '+x['comparison'])
        print('Native v10 still BLOCKED. Differential response is diagnostic only.')
        return 0
    except (ValueError,TypeError,KeyError,OSError,IndexError,ForensicBlocked) as exc:
        print('v13 comparison BLOCKED: '+str(exc),file=sys.stderr)
        return 2


if __name__=='__main__':raise SystemExit(main())
