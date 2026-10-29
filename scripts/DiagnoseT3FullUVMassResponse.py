"""Full-UV scalar mass boundary response and beta decomposition; no fitting.

Only changes the two UV scalar mass-squared boundary entries. Keeps the
historical Yukawas, RGEs, thresholds, quartics, and screens unchanged.
"""
from __future__ import annotations
import argparse
import json
import math
import re
from pathlib import Path
import numpy as np
from scripts import RunT3GeneralScalar48Campaign as bulk
from scripts.ProbeT3FullUVRunning import BLOCKED_LEGACY, REUSED, CASE

FORMAT='t3-fulluv-scalar-mass-response-v1'

def boundary_vector(vec,layout,factor):
    if not (math.isfinite(factor) and factor>0):
        raise ValueError('mass-squared factor must be positive')
    out=np.array(vec,dtype=float,copy=True)
    for name in ('mS1Sq','mS2Sq'):
        b=layout.block(name)
        if b.size!=1 or out[b.start]<=0:
            raise ValueError(f'Invalid positive scalar mass squared: {name}')
        out[b.start]*=factor
    return out

def scalar_values(vec,layout):
    return {name:float(vec[layout.block(name).start]) for name in ('mS1Sq','mS2Sq')}

def run_one(model,factor,entry,new,history,betas_uv,betas_eft,q0):
    # Use the exact UV+EFT beta and integrators from the 48-case campaign.
    vec,layout,_=bulk.uv_base(history['cfg'],entry)
    vec=boundary_vector(vec,layout,factor)
    names=tuple(new)
    # UV initial beta contributions d(mS^2)/d ln(mu), in GeV^2.
    initial=np.r_[vec,[v for z in q0 for v in (z.real,z.imag)]]
    deriv=bulk.s30.combined_rhs(math.log(1e7),initial,layout=layout,names=names,betas=betas_uv)
    initial_beta=scalar_values(deriv,layout)
    uv,qf,nuv=bulk.run_uv(vec,layout,betas_uv,names,q0)
    eft,masses,qs,neft=bulk.run_eft(uv,qf,betas_eft,model,new,())
    chk=bulk.testscreen(uv,eft,masses,qf,qs)
    return {'factor_mSsq':factor,
      'uv_boundary_mSsq_GeV2':scalar_values(vec,layout),
      'beta_mSsq_at_UV_GeV2_per_lnmu':initial_beta,
      'uv_at_muF_mSsq_GeV2':{'mS1Sq':float(uv.mS1Sq),'mS2Sq':float(uv.mS2Sq)},
      'eft_at_muS_mSsq_GeV2':{'mS1Sq':float(eft.mS1Sq),'mS2Sq':float(eft.mS2Sq)},
      'scalar_masses_muS_tev':sorted((np.asarray(masses)/1000).tolist()),
      'uv_nfev':nuv,'eft_nfev':neft,**chk}

def response(rows):
    valid=[r for r in rows if r.get('status')=='Computed']
    if len(valid)<2:return {'status':'InsufficientPoints'}
    baseline=next((v for v in valid if v['factor_mSsq']==1),None)
    if baseline is None:return {'status':'MissingBaseline'}
    res=[]
    b=np.array(baseline['scalar_masses_muS_tev'])
    for row in valid:
        m=np.array(row['scalar_masses_muS_tev'])
        res.append({'factor_mSsq':row['factor_mSsq'],
                    'delta_mass_TeV':(m-b).tolist(),
                    'mass_window_pass':bool(np.all((m>=.5)&(m<=5)))})
    return {'status':'Computed','comparisons_to_baseline':res}

def diagnose(args):
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,63}',args.run_id):
        raise ValueError('Invalid run ID')
    factors=[float(v) for v in args.factors.split(',')]
    if len(set(factors))!=len(factors) or 1.0 not in factors or len(factors)<2:
        raise ValueError('Factors must be distinct and include 1.0 plus another factor')
    if not all(math.isfinite(v) and 0<v<=100 for v in factors):
        raise ValueError('Factors must be within (0,100]')
    registry=bulk.read(bulk.REGISTRY)
    models=list(registry['models'])
    if args.models:
        requested=[v.strip() for v in args.models.split(',')]
        if set(requested)-set(models):raise ValueError('Unknown model(s)')
        models=[v for v in models if v in requested]
    dest=bulk.OUT/args.run_id
    if dest.exists() and not args.resume:
        raise FileExistsError(f'{dest} exists, specify --resume or new --run-id')
    dest.mkdir(parents=True,exist_ok=True)
    output=dest/'mass_response.jsonl'
    done=set()
    if args.resume and output.exists():
        for line in output.read_text(encoding='utf-8').splitlines():
            if line.strip():done.add(json.loads(line)['model'])
    for model in models:
        if model in done:
            print('SKIP',model,'checkpoint exists',flush=True);continue
        row={'format':FORMAT,'model':model,'case':CASE,'factors':factors,'points':[]}
        try:
            if model in BLOCKED_LEGACY:
                row.update(status='ExcludedLegacyC5Parity')
            elif model in REUSED:
                row.update(status='ExcludedSpecialMixedScalar')
            else:
                entry,new,reused=bulk.inventory(registry,model)
                if reused:raise ValueError('Reused quartic implementation unsupported')
                hist=bulk.historical(CASE,model,entry)
                if not hist['uv_on']:raise ValueError('This is not a full-UV case')
                beta_eft=bulk.read(bulk.OUT/args.source_campaign/'betas'/f'{bulk.safe(model)}.json')
                vec,layout,_=bulk.uv_base(hist['cfg'],entry)
                beta_uv,beta_file=bulk.uv_betas(model,entry,new,layout)
                row.update(status='Computed',uv_beta_path=str(beta_file),
                           historical_uv_boundary_mSsq=scalar_values(vec,layout))
                q0=tuple(complex(args.extra_real) for _ in new)
                for f in factors:
                    try:
                        result=run_one(model,f,entry,new,hist,beta_uv,beta_eft,q0)
                        result['status']='Computed'
                        print(f'{model}: mSsq_factor={f:g} masses_TeV={result["scalar_masses_muS_tev"]} rhoF={result["rhoF"]:.4f} safe={result["screen"]}',flush=True)
                    except Exception as e:
                        result={'factor_mSsq':f,'status':'NumericalFailure',
                                'reason':f'{type(e).__name__}: {e}'[:600]}
                        print(f'{model}: factor={f:g} FAILED {result["reason"]}',flush=True)
                    row['points'].append(result)
                row['response']=response(row['points'])
        except Exception as e:
            row.update(status='InputFailure',reason=f'{type(e).__name__}: {e}'[:600])
        with output.open('a',encoding='utf-8') as stream:
            stream.write(json.dumps(row,allow_nan=False)+'\n')
        print('CHECKPOINT',model,row['status'],flush=True)
    print('COMPLETE',output,flush=True)

def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--run-id',required=True)
    ap.add_argument('--source-campaign',default='overnight_six_v1_bulk48')
    ap.add_argument('--models',default='')
    ap.add_argument('--factors',default='0.25,1,4',help='UV scalar mass-SQUARED multiplicative factors; must include 1')
    ap.add_argument('--extra-real',type=float,default=-.01)
    ap.add_argument('--resume',action='store_true')
    diagnose(ap.parse_args(argv))

if __name__=='__main__':main()
