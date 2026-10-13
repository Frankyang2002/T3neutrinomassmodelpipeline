"""Two lambdaT3 benchmark sensitivity scans with validated per-point caching.

Default is read-only inventory. With --run-missing, compute ONLY missing points
using the existing numerical evaluator and precomputed RGBeta/matching artifacts.
No fitting, RGBeta, Matchete or full pipeline is invoked.  The previously
optimised common y1/y2 are preserved; only lambdaT3 and controlled sensitivity
perturbations are changed. Outputs are kept separate from previous results.

From the repository root:
 python -m Numerical.orchestration.RunScalarCouplingCorrelations
 python -m Numerical.orchestration.RunScalarCouplingCorrelations --run-missing
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from copy import deepcopy
from pathlib import Path

from Numerical.orchestration.RunCommonParameterComparisonV2 import (
    PROJECT_ROOT, SCENARIO, OBS, candidate_for_target, models_and_cases,
    load, evaluate,
)
from Numerical.fitting.BenchmarkSensitivity import scaled_config_payload

COUPLINGS = (0.0005, 0.05)
DIRECTIONS = ('T1', 'T2', 'T3')
COMMON_INPUT = (PROJECT_ROOT / 'Reports/output/full/common_yukawa_least_squares'
                / SCENARIO / 'optimized_common_input.json')
OUTPUT = PROJECT_ROOT / 'Reports/output/full/scalar_coupling_correlations' / SCENARIO
FIELDS = ('model','base_model','lambdaT3_reference','direction','gamma','status',
          'error','chi2',*OBS,'input_hash')


def number_tag(value: float) -> str:
    return format(value, '.8g').replace('.', 'p').replace('-', 'm')


def fingerprint(payload, model, direction, gamma):
    """Bind cached predictions to complete numerical input, not only a model label."""
    encoded = json.dumps({'model':model,'direction':direction,'gamma':gamma,
        'payload':payload}, sort_keys=True, separators=(',', ':'), allow_nan=False)
    return hashlib.sha256(encoded.encode('utf-8')).hexdigest()


def scalar(value):
    if isinstance(value, dict):
        if 'real' not in value: raise ValueError('lambdaT3 complex scalar missing real')
        return complex(float(value['real']), float(value.get('imag',0)))
    return complex(float(value),0)


def set_scalar(payload, value):
    spec=payload['base_state']['ordinary']['lambdaT3']
    if not isinstance(spec,dict):
        raise ValueError('Expected complex scalar dictionary for lambdaT3')
    spec['real']=float(value)
    spec['imag']=0.0


def prepare_models():
    configs,cases=models_and_cases()
    if not COMMON_INPUT.is_file():
        raise FileNotFoundError('Missing optimised Yukawa input: '+str(COMMON_INPUT))
    common=load(COMMON_INPUT)
    ordinary=common['shared_ordinary']
    source_name=common['source_model']
    source=deepcopy(load(configs[source_name]))
    for name in ('y1','y2','lambdaT3','MF','mS1Sq','mS2Sq'):
        if name not in ordinary:
            raise KeyError(f'Missing optimised shared input {name}')
        source['base_state']['ordinary'][name]=deepcopy(ordinary[name])
    prepared={name:candidate_for_target(source,load(path)) for name,path in sorted(configs.items())}
    return prepared,cases,common


def read_cache(path):
    found={}
    if not path.exists():return found
    with path.open(newline='',encoding='utf-8-sig') as stream:
        for row in csv.DictReader(stream):
            if row.get('status')!='Success':continue
            key=(row.get('base_model'),row.get('direction'),row.get('gamma'))
            try:
                if not row.get('input_hash') or not math.isfinite(float(row['chi2'])):
                    continue
                if not all(math.isfinite(float(row[k])) for k in OBS):continue
            except (KeyError,ValueError,TypeError):continue
            found[key]=row
    return found


def checkpoint(path,rows):
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix('.csv.tmp')
    with tmp.open('w',newline='',encoding='utf-8') as stream:
        writer=csv.DictWriter(stream,fieldnames=FIELDS)
        writer.writeheader();writer.writerows(rows)
    os.replace(tmp,path)


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-missing',action='store_true',help='Evaluate only uncached points')
    p.add_argument('--points',type=int,default=7)
    p.add_argument('--gamma-min',type=float,default=.9)
    p.add_argument('--gamma-max',type=float,default=1.1)
    p.add_argument('--limit-models',type=int,default=0)
    args=p.parse_args(argv)
    if args.points<3 or args.points%2!=1 or not 0<args.gamma_min<1<args.gamma_max:
        p.error('Require odd --points >=3 and 0<gamma-min<1<gamma-max')
    if args.limit_models<0:p.error('--limit-models must be nonnegative')
    import numpy as np
    target=None
    if args.run_missing:
        from Numerical.fitting.OscillationFit import OscillationFitTarget
        target=OscillationFitTarget.from_json(PROJECT_ROOT/'configs/nufit_6_1_ic24_no.json')
        if tuple(target.observable_names)!=OBS or target.ordering!='NO':
            raise ValueError('Unexpected oscillation-fit target')
    configs,cases,common=prepare_models()
    names=list(configs)[:args.limit_models or None]
    gammas=[float(g) for g in np.linspace(args.gamma_min,args.gamma_max,args.points)]
    tally={'cached':0,'missing':0,'computed':0,'failed':0}
    for lam in COUPLINGS:
        folder=OUTPUT/('lambda_'+number_tag(lam))
        path=folder/'scan_points.csv'
        cache=read_cache(path)
        rows=[]
        for name in names:
            base=deepcopy(configs[name]);set_scalar(base,lam)
            if scalar(base['base_state']['ordinary']['lambdaT3'])!=complex(lam,0):
                raise RuntimeError('Scalar coupling change failed')
            for direction in DIRECTIONS:
                for gamma in gammas:
                    modified=scaled_config_payload(base,direction,gamma)
                    key=(name,direction,format(gamma,'.17g'))
                    h=fingerprint(modified,name,direction,gamma)
                    cached=cache.get(key)
                    # Existing results are reusable only for exactly the same input.
                    if cached is not None and cached['input_hash']==h:
                        rows.append(cached);tally['cached']+=1;continue
                    row=dict(model=f'{name}__lambda_{number_tag(lam)}',base_model=name,
                             lambdaT3_reference=lam,direction=direction,
                             gamma=format(gamma,'.17g'),status='Missing',error='',
                             chi2='',input_hash=h)
                    row.update({o:'' for o in OBS})
                    if args.run_missing:
                        try:
                            chi2,pred=evaluate(modified,cases[name],target)
                            if not math.isfinite(chi2) or not all(math.isfinite(float(pred[o])) for o in OBS):
                                raise ValueError('Non-finite prediction')
                            row.update(status='Success',chi2=chi2,**pred)
                            tally['computed']+=1
                        except Exception as exc:
                            row.update(status='Failed',error=f'{type(exc).__name__}: {exc}')
                            tally['failed']+=1
                    else:tally['missing']+=1
                    rows.append(row)
                    if args.run_missing:checkpoint(path,rows)
        if args.run_missing:
            checkpoint(path,rows)
            manifest={'coupling_reference':lam,'optimised_common_input':str(COMMON_INPUT),
                      'source_model':common['source_model'],'points':args.points,
                      'gamma_min':args.gamma_min,'gamma_max':args.gamma_max,
                      'meaning':'At gamma=1 lambdaT3 equals coupling_reference. T1 and T2 vary it by gamma.',
                      'note':'Same optimised y1/y2 matrices for both scalar benchmarks; no refit.',
                      'cache':'Each successful row carries sha256 of full deformed numerical payload.'}
            (folder/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
        print(f'lambdaT3={lam:g}: {len(rows)} sensitivity points; stored in {path}',flush=True)
    print('Summary:',tally)
    if not args.run_missing:print('Dry run: no files written. Use --run-missing to evaluate missing points.')
    return 2 if tally['failed'] else 0

if __name__=='__main__':
    raise SystemExit(main())
