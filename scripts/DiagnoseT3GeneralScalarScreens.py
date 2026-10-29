"""Read-only, targeted diagnostics for the existing 16x3 general-scalar campaign.

Runs real candidate UV/EFT1 RGEs at historical fitted Yukawas (and optionally
fresh Sobol points). Reports each mass/perturbativity screen independently;
for legacy C5 models compares fast and authoritative matrix elements.
Does NOT alter matching, physics bounds, saved campaigns, or fitted configs.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
import re
from collections import Counter
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from scripts import RunT3GeneralScalar48Campaign as bulk

DEFAULT_CASES = (
    ('A_+0', 'small_uv_running_5e-9'),
    ('A_-2', 'small_uv_running_5e-9'),
    ('B_+1', 'small_uv_running_5e-9'),
    ('A_-4', 'large_uv_frozen_0p5'),
)
FIELDS = ('rhoF','scalar_min_tev','scalar_max_tev','yukawa_sigma','quartic_max')


def finite(x):
    try:
        v=float(x)
        return v if math.isfinite(v) else None
    except (TypeError, ValueError, OverflowError):
        return None


def screen_details(uv, eft, masses, qf, qs):
    """Same bounds/quantities as bulk.testscreen, but separately labelled."""
    from Numerical.running.FermionThresholdBoundary import fermion_threshold_masses
    from Numerical.running.IntermediateScalarRGBetaEvaluator import intermediate_scalar_state_environment
    fermion=np.asarray(fermion_threshold_masses(uv).masses_gev,dtype=float)
    masses=np.asarray(masses,dtype=float)
    scalar_sq=np.asarray((uv.mS1Sq,uv.mS2Sq),dtype=float)
    positive=bool(np.isfinite(fermion).all() and np.isfinite(masses).all() and
                  np.isfinite(scalar_sq).all() and np.min(fermion)>0 and
                  np.min(masses)>0 and np.min(scalar_sq)>0)
    y=max(float(np.linalg.svd(np.asarray(z,dtype=complex),compute_uv=False)[0])
          for z in (uv.y1,uv.y2))
    env=intermediate_scalar_state_environment(eft)
    quartics=[abs(complex(v)) for k,v in env.items() if k.startswith('lambda')]
    quartics.extend(abs(complex(v)) for v in (*qf,*qs))
    qm=float(max(quartics))
    rho=float(np.sqrt(scalar_sq.max())/fermion.min()) if positive else None
    min_mass=finite(masses.min()/1000)
    max_mass=finite(masses.max()/1000)
    checks={
       'positive_masses': positive,
       'rhoF_le_0p5': bool(positive and rho<=0.5),
       'scalar_min_ge_0p5TeV':bool(positive and masses.min()>=500),
       'scalar_max_le_5TeV':bool(positive and masses.max()<=5000),
       'yukawa_sigma_le_sqrt4pi':bool(math.isfinite(y) and y<=math.sqrt(4*math.pi)),
       'quartic_max_le_4pi':bool(math.isfinite(qm) and qm<=4*math.pi),
    }
    return {'passed':all(checks.values()), 'checks':checks,
            'failed_checks':[k for k,passed in checks.items() if not passed],
            'rhoF':finite(rho),'scalar_min_tev':min_mass,'scalar_max_tev':max_mass,
            'yukawa_sigma':finite(y),'quartic_max':finite(qm),
            'limits':{'rhoF':0.5,'scalar_min_tev':0.5,'scalar_max_tev':5.0,
                      'yukawa_sigma':math.sqrt(4*math.pi),'quartic_max':4*math.pi}}


def matrix_json(matrix):
    a=np.asarray(matrix,dtype=complex)
    return {'real':a.real.tolist(),'imag':a.imag.tolist()}


def legacy_parity(kernel, uv, eft):
    from Numerical.running.FinalC5TrajectoryAdapter import extract_final_c5_inputs
    from Numerical.running.SMWeinbergStage import evaluate_symbolic_c5
    view=SimpleNamespace(uv_threshold_state=uv,eft1_threshold_state=eft,
                         mu_scalar_threshold_gev=1e3)
    cfg=extract_final_c5_inputs(view).as_config()
    exact=np.asarray(evaluate_symbolic_c5(kernel.path,cfg),dtype=complex)
    fast=np.asarray(kernel.evaluate(uv,eft),dtype=complex)
    difference=fast-exact
    denom=max(float(np.max(abs(exact))),1e-30)
    err=float(np.max(abs(difference))/denom)
    # Diagnostic includes all matrix entries; never changes the 3e-8 tolerance.
    return {'tolerance':3e-8,'relative_max_error':err,'passes':bool(err<=3e-8),
            'reference_max_abs':float(np.max(abs(exact))),
            'fast_max_abs':float(np.max(abs(fast))),
            'max_abs_difference':float(np.max(abs(difference))),
            'reference_c5':matrix_json(exact),'fast_c5':matrix_json(fast),
            'difference_c5':matrix_json(difference)}


def investigate(model, case, registry, *, extra_real=-0.01, sobol_points=0, seed=20261030):
    from Numerical.fitting.SobolBenchmarkSearch import generate_sobol_parameter_points, BenchmarkParameter
    entry,new,reused=bulk.inventory(registry,model)
    if reused:
        raise ValueError('Model uses separately fitted mixed-quadratic B/C implementation')
    history=bulk.historical(case,model,entry)
    beta_path=bulk.OUT / '_not_used_'
    # Use an existing saved native EFT beta payload. Avoid regenerating Wolfram beta files.
    candidates=[bulk.OUT/'overnight_six_v1_bulk48'/'betas'/f'{bulk.safe(model)}.json',
                bulk.OUT/'bulk48_v1'/'betas'/f'{bulk.safe(model)}.json']
    paths=[path for path in candidates if path.is_file()]
    if not paths:
        raise FileNotFoundError('No saved candidate EFT beta: '+', '.join(map(str,candidates)))
    beta_path=paths[0]
    payload=bulk.read(beta_path)
    kernel,kernel_path,kind=bulk.matching(model)
    q0=tuple(complex(extra_real) for _ in new)
    frozen=None
    uv_betas=None
    if history['uv_on']:
        vec,layout,_=bulk.uv_base(history['cfg'],entry)
        uv_betas,_=bulk.uv_betas(model,entry,new,layout)
    else:
        vec,layout,unpack=bulk.uv_base(history['cfg'],entry)
        frozen=bulk.run_eft(unpack(vec,layout,mu_gev=1e5),q0,payload,model,new,())
    trials=[('historical_seed',np.array(history['seed'],dtype=float))]
    if sobol_points:
        if sobol_points<2 or sobol_points&(sobol_points-1):
            raise ValueError('--sobol-points must be a power of two >= 2')
        lo,hi=bulk.b21.search_box(history['seed'],local_fraction=1.5,
                      floor_fraction=.05,abs_yukawa_max=math.sqrt(4*math.pi))
        definitions=[BenchmarkParameter(name,f'ordinary.{name.split("_")[0]}.real.{name[-2]}.{name[-1]}',float(l),float(h))
                     for name,l,h in zip(bulk.NAMES,lo,hi)]
        points=generate_sobol_parameter_points(definitions,sobol_points,seed=seed)
        trials.extend((f'sobol_{i+1}',np.array([p[k] for k in bulk.NAMES],dtype=float))
                      for i,p in enumerate(points))
    result={'model':model,'case':case,'source_beta':str(beta_path.relative_to(bulk.ROOT)),
            'matching_source':str(kernel_path.relative_to(bulk.ROOT)),
            'matching_kind':kind,'sobol_points':sobol_points,'trials':[]}
    for name,vector in trials:
        row={'label':name}
        try:
            cfg=bulk.update_yukawas(history,vector)
            vec,layout,unpack=bulk.uv_base(cfg,entry)
            if history['uv_on']:
                uv,qf,nuv=bulk.run_uv(vec,layout,uv_betas,new,q0)
                eft,masses,qs,neft=bulk.run_eft(uv,qf,payload,model,new,())
            else:
                uv=unpack(vec,layout,mu_gev=1e5)
                eft,masses,qs,neft=frozen
                qf=q0;nuv=0
            s=screen_details(uv,eft,masses,qf,qs)
            row.update({'status':'ScreenPassed' if s['passed'] else 'ScreenRejected',
                        'screen':s,'uv_nfev':nuv,'eft1_nfev':neft})
            old=bulk.testscreen(uv,eft,masses,qf,qs)
            if bool(old['screen'])!=s['passed']:
                raise RuntimeError('Diagnostic screen differs from existing bulk.testscreen')
            if s['passed']:
                if kind=='legacy-common-scale':
                    row['legacy_parity']=legacy_parity(kernel,uv,eft)
                    row['status']='ParityPassed' if row['legacy_parity']['passes'] else 'ParityFailed'
                else:
                    row['parity_relative_error']=float(kernel.parity_check(uv,eft))
                    row['status']='ParityPassed'
        except (ValueError,RuntimeError,ArithmeticError,OverflowError,TypeError,
                FloatingPointError,np.linalg.LinAlgError) as exc:
            row.update(status='EvaluationError',error=repr(exc))
        result['trials'].append(row)
        print(model,case,name,row['status'],
              row.get('screen',{}).get('failed_checks',[]),
              row.get('legacy_parity',{}).get('relative_max_error'),flush=True)
    result['counts']=dict(Counter(x['status'] for x in result['trials']))
    result['failed_checks']=dict(Counter(k for x in result['trials']
                                        for k in x.get('screen',{}).get('failed_checks',[])))
    return result


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-id',required=True)
    p.add_argument('--cases',nargs='*',help='MODEL:SCENARIO. Defaults to three full-UV cases and A_-4 high/frozen.')
    p.add_argument('--sobol-points',type=int,default=0,
                   help='Optional fresh diagnostic Sobol points (power of two). Default: historical seed only.')
    p.add_argument('--seed',type=int,default=20261030)
    p.add_argument('--extra-real',type=float,default=-0.01)
    args=p.parse_args(argv)
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,63}',args.run_id):
        p.error('Invalid run ID')
    if args.sobol_points and (args.sobol_points<2 or args.sobol_points&(args.sobol_points-1)):
        p.error('--sobol-points must be a power of two >= 2')
    selected=DEFAULT_CASES if not args.cases else tuple(tuple(x.split(':',1)) for x in args.cases)
    if any(len(x)!=2 or x[1] not in bulk.STUDIES for x in selected):
        p.error('Expected MODEL:SCENARIO entries')
    dest=bulk.ROOT/'output/studies/general_scalar_extended/diagnostics'/args.run_id
    if dest.exists():
        p.error('Diagnostic output already exists; choose a new --run-id')
    registry=bulk.read(bulk.REGISTRY)
    results=[]
    for model,case in selected:
        try:
            results.append(investigate(model,case,registry,extra_real=args.extra_real,
                 sobol_points=args.sobol_points,seed=args.seed))
        except Exception as exc:
            results.append({'model':model,'case':case,'status':'Blocked','error':repr(exc)})
            print(model,case,'BLOCKED',repr(exc),flush=True)
    dest.mkdir(parents=True,exist_ok=False)
    (dest/'diagnostics.json').write_text(json.dumps({'format':'t3-screen-and-parity-diagnostics-v1',
       'results':results,'no_physics_modified':True},indent=2,allow_nan=False)+'\n',encoding='utf-8')
    rows=[]
    for item in results:
        for trial in item.get('trials',[]):
            s=trial.get('screen',{})
            rows.append({'model':item['model'],'case':item['case'],'trial':trial['label'],
             'status':trial['status'],'failed_checks':';'.join(s.get('failed_checks',[])),
             **{k:s.get(k) for k in FIELDS},
             'parity_error':trial.get('legacy_parity',{}).get('relative_max_error',
                             trial.get('parity_relative_error')),
             'error':trial.get('error','')})
    with (dest/'diagnostics.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=['model','case','trial','status','failed_checks',
                                             *FIELDS,'parity_error','error'])
        writer.writeheader();writer.writerows(rows)
    print('Saved:',dest/'diagnostics.json',flush=True)
    return int(any(x.get('status')=='Blocked' for x in results))

if __name__=='__main__':
    raise SystemExit(main())
