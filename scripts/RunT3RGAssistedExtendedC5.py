#!/usr/bin/env python3
"""Recompute six full coupled UV endpoints, then conditionally evaluate Matchete C5.

Requires the user's already verified locally generated native Matchete exports,
RGBeta native UV reports, and fitted reference config. Refuses to claim
physical normalization, EFT1 running, or sequential threshold matching.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path
import json
import sys

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
from Numerical.experimental import T3RGAssistedExtendedC5 as BR

DEFAULT_SCENARIO=ROOT/'Numerical/configs/t3_rg_assisted_c5_convention_hypothesis.json'


def atomic_json(path:Path,obj:dict)->None:
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_name(path.name+'.tmp')
    tmp.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    tmp.replace(path)


def main()->int:
    parser=argparse.ArgumentParser(description=__doc__)
    act=parser.add_mutually_exclusive_group(required=True)
    act.add_argument('--self-test',action='store_true')
    act.add_argument('--pilot',action='store_true',help='Only A_p0')
    act.add_argument('--all',action='store_true',help='Six extended-Yukawa models')
    act.add_argument('--model',choices=list(BR.EXTRA_QUARTIC))
    parser.add_argument('--scenario',type=Path,default=DEFAULT_SCENARIO)
    parser.add_argument('--native-source',type=Path,help='Optional native 24-file ZIP instead of locally exported files')
    parser.add_argument('--output-dir',type=Path,default=BR.OUT)
    args=parser.parse_args()
    if args.self_test:
        BR.self_test()
        from scripts.TestT3RGAssistedExtendedC5 import self_test as integrated_test
        integrated_test()
        return 0
    config=json.loads(args.scenario.read_text(encoding='utf-8-sig'))
    scenario=BR.read_scenario(config)
    models=list(BR.EXTRA_QUARTIC) if args.all else ['A_p0'] if args.pilot else [args.model]
    out=args.output_dir.resolve()
    results=[]
    for j,model in enumerate(models,1):
        try:
            result=BR.run_model(model,scenario,root=ROOT,native_source=args.native_source)
            dest=out/'models'/f'{model}.json'
            atomic_json(dest,result)
            first=result['phase_scan']['scenarios'][0]
            row={'model':model,'status':result['status'],
                 'report_path':str(dest),
                 'endpoint_mu_gev':result['uv_endpoint_gev'],
                 'loop_C5_frobenius_gev_inverse':first['frobenius_C5_total_gev_inverse'],
                 'extra_C5_frobenius_gev_inverse':first['frobenius_C5_extra_gev_inverse'],
                 'hypothesis_masses_ev_phase0':first['conditional_observables']['masses_ev']}
            print(f'[{j}/{len(models)}] CONDITIONAL PASS {model}: '
                  f'||C5||={row["loop_C5_frobenius_gev_inverse"]:.5e} GeV^-1; '
                  f'direct_extra={row["extra_C5_frobenius_gev_inverse"]:.5e}',flush=True)
        except Exception as exc:
            row={'model':model,'status':'BLOCKED',
                 'error':type(exc).__name__+': '+str(exc)}
            print(f'[{j}/{len(models)}] BLOCKED {model}: {row["error"]}',file=sys.stderr,flush=True)
        results.append(row)
        summary={'format':'t3-conditional-rg-assisted-extended-c5-batch-v1',
                 'status':'ConditionalBatchPartialOrBlocked',
                 'generated_at_utc':datetime.now(timezone.utc).isoformat(),
                 'scenario_file':str(args.scenario.resolve()),
                 'scenario_values':config,'models_requested':models,'models_completed':sum(r['status']!='BLOCKED' for r in results),
                 'rows':results,'physical_matchete_rgbeta_map_certified':False,
                 'full_hierarchical_threshold_matching_certified':False,
                 'low_scale_neutrino_predictions_certified':False}
        if len(results)==len(models) and all(r['status']!='BLOCKED' for r in results):
            summary['status']='SixModelConditionalUVToOneStepC5Passed' if len(results)==6 else 'ConditionalUVToOneStepC5Passed'
        atomic_json(out/'summary.json',summary)
    print('Summary:',out/'summary.json')
    print('NOTE: all C5 and neutrino-mass values are CONDITIONAL on unverified basis/mass and matching assumptions.')
    print('The runner STOPS above the first heavy threshold; no EFT1/final Weinberg RG evolution was performed.')
    return 0 if all(r['status']!='BLOCKED' for r in results) else 4


if __name__=='__main__':
    try:raise SystemExit(main())
    except (ValueError,OSError,KeyError,TypeError,json.JSONDecodeError) as exc:
        print('BLOCKED:',type(exc).__name__+':',exc,file=sys.stderr)
        raise SystemExit(3)
