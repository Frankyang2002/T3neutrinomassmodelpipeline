"""Exploratory A(-2) *joint* scalar-UV/Yukawa Sobol search.

Uses unchanged RunT3GeneralScalar48Campaign UV/EFT/C5 calculations and screens.
This is a candidate general-quartic study, NOT full-general-scalar/sequential EFT.
The historical small topology lambdaT3 is fixed at 5e-9 by default.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
import json
import math
from pathlib import Path
import re
import sys

import numpy as np
from scipy.optimize import least_squares
from scipy.stats import qmc

from scripts import RunT3GeneralScalar48Campaign as bulk
from scripts import FitT3BPlus1GeneralThreeStudies as b21

CASE = 'small_uv_running_5e-9'
FORMAT = 't3-joint-scalar-yukawa-exploration-v1'


def quartic_names(cfg, maximum):
    ordinary = cfg['base_state']['ordinary']
    names = []
    for name, raw in ordinary.items():
        if not name.startswith('lambda') or name == 'lambdaT3':
            continue
        if isinstance(raw, dict) and set(raw) == {'real','imag'}:
            z = b21.scalar_complex(raw)
            if abs(z.imag) > 1e-12:
                continue  # real-coefficient scan only, without altering CP conventions
        elif isinstance(raw, (int,float)):
            z = complex(float(raw))
        else:
            continue
        if math.isfinite(z.real):
            names.append(name)
    return tuple(names[:maximum])


def set_quartic(raw, value):
    if isinstance(raw, dict):
        raw['real'] = float(value)
        return raw
    return float(value)


def configure(history, x, factors, qnames, qvals):
    """Return a new history; never modify historical source inputs."""
    if len(factors) != 2 or any(not math.isfinite(v) or v <= 0 for v in factors):
        raise ValueError('Invalid two positive mass factors')
    if len(qnames) != len(qvals):
        raise ValueError('Quartic coordinate mismatch')
    h = dict(history)
    cfg = bulk.update_yukawas(history, np.asarray(x, dtype=float))
    ordinary = cfg['base_state']['ordinary']
    for name, factor in zip(('mS1Sq','mS2Sq'),factors):
        ordinary[name] = float(ordinary[name])*float(factor)**2
    for name,value in zip(qnames,qvals):
        ordinary[name] = set_quartic(ordinary[name],value)
    h['cfg'] = cfg
    return h


def screen_point(history, entry, model, new, payload, q0, x, factors, qnames, qvals, uvbetas):
    h = configure(history,x,factors,qnames,qvals)
    vec,layout,_ = bulk.uv_base(h['cfg'],entry)
    uv, qf, nuv = bulk.run_uv(vec,layout,uvbetas,new,q0)
    eft,masses,qs,neft = bulk.run_eft(uv,qf,payload,model,new,())
    chk = bulk.testscreen(uv,eft,masses,qf,qs)
    return chk, {'uv_nfev':nuv,'eft_nfev':neft,
                 'muF_mS1Sq':float(uv.mS1Sq),'muF_mS2Sq':float(uv.mS2Sq)}


def atomic_json(path, obj):
    path = Path(path)
    tmp = path.with_suffix('.tmp')
    tmp.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    tmp.replace(path)


def bounds_from_seed(seed):
    return b21.search_box(seed,local_fraction=1.5,floor_fraction=.05,
                          abs_yukawa_max=math.sqrt(4*math.pi))


def sample_parameters(row, lo, hi, qnames, qbound, mass_low, mass_high):
    yy = lo + (hi-lo)*row[:18]
    # log-uniform mass scaling, independent for the two scalars
    ff = np.exp(np.log(mass_low)+(np.log(mass_high)-np.log(mass_low))*row[18:20])
    qq = -qbound + 2*qbound*row[20:20+len(qnames)]
    return yy,ff,qq


def parse(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-id',required=True)
    p.add_argument('--beta-campaign',default='overnight_six_v1_bulk48')
    p.add_argument('--case',choices=[CASE],default=CASE)
    p.add_argument('--model',choices=['A_-2'],default='A_-2')
    p.add_argument('--sobol-points',type=int,default=128)
    p.add_argument('--max-rounds',type=int,default=8)
    p.add_argument('--max-nfev',type=int,default=100)
    p.add_argument('--local-seeds',type=int,default=3)
    p.add_argument('--seed',type=int,default=20261030)
    p.add_argument('--mass-factor-min',type=float,default=0.05)
    p.add_argument('--mass-factor-max',type=float,default=5.0)
    p.add_argument('--quartic-bound',type=float,default=1.0)
    p.add_argument('--max-quartics',type=int,default=8)
    p.add_argument('--extra-real',type=float,default=-0.01)
    p.add_argument('--resume',action='store_true')
    return p.parse_args(argv)


def run(args):
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,63}',args.run_id):
        raise ValueError('Invalid run id')
    if (args.sobol_points<2 or args.sobol_points & (args.sobol_points-1)
            or not 1<=args.max_rounds<=1000 or not 0<=args.max_nfev<=1000
            or not 1<=args.local_seeds<=32
            or not 0<args.mass_factor_min<args.mass_factor_max<=100
            or not 0<args.quartic_bound<=4*math.pi
            or not 0<=args.max_quartics<=32
            or not 0<abs(args.extra_real)<=.5):
        raise ValueError('Invalid search budget/bounds')
    root=bulk.OUT/args.run_id
    if root.exists() and not args.resume:
        raise ValueError('Run exists; use --resume or a new run id')
    root.mkdir(parents=True,exist_ok=True)
    statefile=root/'state.json'
    registry=bulk.read(bulk.REGISTRY)
    entry,new,reused=bulk.inventory(registry,args.model)
    if reused: raise ValueError('Special mixed-quadratic models unsupported')
    hist=bulk.historical(args.case,args.model,entry)
    payload=bulk.read(bulk.OUT/args.beta_campaign/'betas'/f'{bulk.safe(args.model)}.json')
    kernel,match_path,match_kind=bulk.matching(args.model)
    from Numerical.orchestration import RunSmallLambdaRefit as refit
    from Numerical.fitting.OscillationFit import OscillationFitTarget
    target=OscillationFitTarget.from_json(refit.TARGET)
    if target.ordering!='NO':raise ValueError('Expected normal ordering')
    vec,layout,_=bulk.uv_base(hist['cfg'],entry)
    uvbetas,uvpath=bulk.uv_betas(args.model,entry,new,layout)
    q0=tuple(complex(args.extra_real) for _ in new)
    qnames=quartic_names(hist['cfg'],args.max_quartics)
    lo,hi=bounds_from_seed(hist['seed'])
    metadata={'format':FORMAT,'model':args.model,'case':args.case,
              'beta_campaign':args.beta_campaign,'qnames':list(qnames),
              'sobol_points':args.sobol_points,'max_rounds':args.max_rounds,
              'max_nfev':args.max_nfev,'local_seeds':args.local_seeds,'seed':args.seed,
              'mass_limits':[args.mass_factor_min,args.mass_factor_max],
              'quartic_bound':args.quartic_bound,'extra_real':args.extra_real,
              'historical_lambdaT3_fixed':hist['lambda_fixed'],
              'matching_kind':match_kind,
              'physics_caveat':'candidate general-quartic; no general quadratic or sequential EFT matching'}
    if statefile.exists():
        data=bulk.read(statefile)
        if data.get('metadata')!=metadata:raise ValueError('Resume settings differ from saved run')
    else:
        data={'metadata':metadata,'rounds':[], 'best_fit':None,'best_screen':None}
        atomic_json(statefile,data)
    if data['best_fit'] is not None and data['best_fit']['chi2']<=5.86:
        print('Already accepted:',data['best_fit']['chi2']);return 0
    for ir in range(len(data['rounds']),args.max_rounds):
        sobol=qmc.Sobol(d=20+len(qnames),scramble=True,seed=args.seed+ir)
        points=sobol.random_base2(int(math.log2(args.sobol_points)))
        # Include historical seed in the first round.
        candidates=[]
        if ir==0:
            candidates.append((hist['seed'],np.ones(2),
                np.array([b21.scalar_complex(hist['cfg']['base_state']['ordinary'][n]).real if isinstance(hist['cfg']['base_state']['ordinary'][n],dict) else float(hist['cfg']['base_state']['ordinary'][n]) for n in qnames]),'historical'))
        for i,z in enumerate(points):
            y,f,q=sample_parameters(z,lo,hi,qnames,args.quartic_bound,
                                     args.mass_factor_min,args.mass_factor_max)
            candidates.append((y,f,q,f'sobol_{i+1}'))
        passing=[]; fail_counts={}; error_count=0
        best_screen=None
        for i,(y,f,q,label) in enumerate(candidates):
            try:
                screen,diag=screen_point(hist,entry,args.model,new,payload,q0,
                                         y,f,qnames,q,uvbetas)
                if screen['screen']:
                    passing.append((y,f,q,label,screen))
                    if best_screen is None:best_screen={'label':label,'screen':screen,
                                 'mass_factors':f.tolist(),'quartics':dict(zip(qnames,q.tolist()))}
                else:
                    key=screen.get('reason','Rejected')
                    fail_counts[key]=fail_counts.get(key,0)+1
            except (ValueError,RuntimeError,ArithmeticError,TypeError,OverflowError,
                    FloatingPointError,np.linalg.LinAlgError) as exc:
                error_count+=1
                k='Exception: '+str(exc)[:100]
                fail_counts[k]=fail_counts.get(k,0)+1
        best_this=None
        # Fit only screen-passing scalar/Yukawa points. Preserve all original screens.
        # Reuse exact candidate Evaluator, which calls C5 and low-energy oscillation fit.
        for j,(y,f,q,label,screen) in enumerate(passing[:args.local_seeds]):
            h=configure(hist,y,f,qnames,q)
            evaluator=bulk.Evaluator(args.case,args.model,entry,new,h,payload,kernel,target,args.extra_real)
            first=evaluator.evaluate(y,label)
            def candidate(row,x):
                nonlocal best_this
                if row.get('chi2') is not None and (best_this is None or row['chi2']<best_this['chi2']):
                    best_this={'chi2':row['chi2'],'prediction':row.get('prediction'),
                               'parameters':dict(zip(bulk.NAMES,map(float,x))),
                               'mass_factors':f.tolist(),
                               'quartics':dict(zip(qnames,map(float,q))),
                               'screen':{k:row.get(k) for k in ('rhoF','scalar_masses_muS_tev','max_yukawa_singular_muF','max_scalar_quartic_abs')},
                               'round':ir+1,'label':label}
            candidate(first,y)
            if args.max_nfev and first.get('chi2') is not None and first['chi2']>5.86:
                z0=np.clip((y-lo)/(hi-lo),1e-8,1-1e-8)
                def residual(z):
                    xx=lo+(hi-lo)*z
                    row=evaluator.evaluate(xx,'local')
                    candidate(row,xx)
                    return np.asarray(row.get('whitened_residual',[1e4]*5),dtype=float)
                least_squares(residual,z0,bounds=(np.zeros(18),np.ones(18)),
                              max_nfev=args.max_nfev,xtol=1e-6,ftol=1e-6,gtol=1e-6)
            if best_this and best_this['chi2']<=5.86:break
        rec={'round':ir+1,'tested':len(candidates),'screen_passed':len(passing),
             'rejected_reasons':fail_counts,'exceptions':error_count,
             'best_chi2':best_this['chi2'] if best_this else None}
        data['rounds'].append(rec)
        if best_screen is not None:data['best_screen']=best_screen
        if best_this and (data['best_fit'] is None or best_this['chi2']<data['best_fit']['chi2']):
            data['best_fit']=best_this
            xx=np.array([best_this['parameters'][k] for k in bulk.NAMES])
            f=best_this['mass_factors'];q=np.array([best_this['quartics'][k] for k in qnames])
            h=configure(hist,xx,f,qnames,q)
            h['cfg']['exploratory_joint_search']={'run_id':args.run_id,'chi2':best_this['chi2'],
                 'candidate_only':True,'full_sequential_matching_certified':False}
            atomic_json(root/'best_fitted_config.json',h['cfg'])
        atomic_json(statefile,data)
        print(f'ROUND {ir+1}/{args.max_rounds} mass-safe {len(passing)}/{len(candidates)} '
              f'best-round={rec["best_chi2"]} best-overall={(data["best_fit"] or {}).get("chi2")}',flush=True)
        if data['best_fit'] and data['best_fit']['chi2']<=5.86:break
    print('Saved:',statefile,flush=True)
    return 0


def main(argv=None):
    return run(parse(argv))

if __name__=='__main__':sys.exit(main())
