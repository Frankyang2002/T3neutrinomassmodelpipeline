"""Six-model coupled-RGE sensitivity of the original T3 Yukawa flavor tensor.

P_ij = conj(y1_ir)conj(y2_jr) + conj(y2_ir)conj(y1_jr), summed over r.
This is only the unweighted flavor factor in the ordinary-T3 C5 kernel.
It is NOT the extended matched C5: loop weights, heavy mass diagonalization,
matching phases and new scalar couplings are not included.

Read existing 16 general RGBeta exports and use the already-tested fully
coupled UV integrator. Recompute six seeded and zero-extra trajectories.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys
from typing import Any

import numpy as np

ROOT=Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
from Numerical.experimental.T3GeneralYukawaCoupledUV import (
    EXTRAS,REPORTS,CONFIGS,DynamicLayout,load_case,prepare_state,
    endpoint_for,integrate,matrix,encode,norm
)

OUTPUT=ROOT/'Reports/output/workstream_a/extended_c5_native/kernel_proxy'


def need(ok:bool,msg:str):
    if not ok:raise ValueError(msg)


def original_flavor_tensor(y1:Any,y2:Any)->np.ndarray:
    a=matrix(y1);b=matrix(y2)
    p=np.conjugate(a)@np.conjugate(b).T + np.conjugate(b)@np.conjugate(a).T
    need(bool(np.allclose(p,p.T,atol=1e-23,rtol=1e-12)),
         'Incorrect symmetric Yukawa contraction')
    return p


def run_one(model:str,seed_factor:float=.3)->dict:
    need(model in EXTRAS,'No independent extra Yukawa for model '+model)
    payload,cfg=load_case(model,REPORTS,CONFIGS)
    extra=EXTRAS[model]
    st=prepare_state(model,payload,cfg,seed_factor)
    uv,endpoint=endpoint_for(cfg,st)
    lay,sol,_=integrate(payload,st,uv,endpoint,n_save=3)
    final=lay.unpack(sol.y[:,-1])
    zero={k:(v.copy() if isinstance(v,np.ndarray) else v) for k,v in st.items()}
    zero[extra]=np.zeros((3,3),dtype=complex)
    _,ctrl,_=integrate(payload,zero,uv,endpoint,n_save=3)
    final_zero=lay.unpack(ctrl.y[:,-1])
    pin=original_flavor_tensor(st['y1'],st['y2'])
    pseed=original_flavor_tensor(final['y1'],final['y2'])
    pzero=original_flavor_tensor(final_zero['y1'],final_zero['y2'])
    delta=pseed-pzero
    return {'status':'OriginalT3YukawaFlavorTensorSensitivityOnly',
        'model':model,'extra_yukawa':extra,
        'mu_initial_gev':uv,'mu_final_gev':endpoint,
        'yukawa_basis':'native RGBeta coordinates; not physical Matchete inputs',
        'proxy_definition':'conjugate(y1)@conjugate(y2).T + conjugate(y2)@conjugate(y1).T',
        'proxy_initial':encode(pin),
        'proxy_seeded_endpoint':encode(pseed),
        'proxy_zero_extra_endpoint':encode(pzero),
        'proxy_delta_seeded_minus_zero':encode(delta),
        'proxy_delta_frobenius':norm(delta),
        'proxy_delta_relative_to_control':norm(delta)/max(norm(pzero),1e-300),
        'seeded_extra_endpoint_norm':norm(final[extra]),
        'zero_seed_extra_endpoint_norm':norm(final_zero[extra]),
        'seed_factor':seed_factor,
        'solver_rhs_evaluations_seeded':int(sol.nfev),
        'solver_rhs_evaluations_control':int(ctrl.nfev),
        'no_matchete_rgbeta_action_dictionary_applied':True,
        'not_a_neutrino_mass_prediction':True,
        'not_a_matched_C5':True,
        'loop_weights_not_included':True,
        'hierarchical_thresholds_not_included':True}


def atomic(path:Path,obj:dict):
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_name(path.name+'.tmp')
    tmp.write_text(json.dumps(obj,indent=2)+'\n',encoding='utf-8');tmp.replace(path)


def self_test():
    a=np.eye(3,dtype=complex)*(.01+.02j)
    b=np.diag([.02,-.03j,.04]).astype(complex)
    p=original_flavor_tensor(a,b)
    need(np.allclose(p,2*np.conjugate(a)@np.conjugate(b).T),
         'Diagonal flavor tensor wrong')
    need(np.allclose(original_flavor_tensor(a,0*b),0),'Zero control failed')
    try:original_flavor_tensor(np.eye(2),np.eye(2))
    except ValueError:pass
    else:raise AssertionError('Invalid 2x2 Yukawa accepted')
    need(len(EXTRAS)==6,'Six model catalogue changed')
    print('PASS: symmetric baseline-Yukawa flavor tensor, exact diagonal comparison and negative controls')


def main()->int:
    ap=argparse.ArgumentParser(description=__doc__)
    g=ap.add_mutually_exclusive_group(required=True)
    g.add_argument('--self-test',action='store_true')
    g.add_argument('--all',action='store_true')
    g.add_argument('--model',choices=list(EXTRAS))
    ap.add_argument('--seed-factor',type=float,default=.3)
    args=ap.parse_args()
    if args.self_test:self_test();return 0
    models=[args.model] if args.model else list(EXTRAS)
    rows=[]
    for n,model in enumerate(models,1):
        try:
            row=run_one(model,seed_factor=args.seed_factor)
            atomic(OUTPUT/'models'/f'{model}.json',row)
            rows.append({'model':model,'status':row['status'],
                'relative_delta':row['proxy_delta_relative_to_control'],
                'delta_frobenius':row['proxy_delta_frobenius']})
            print(f'[{n}/{len(models)}] PASS {model}: relative baseline flavor tensor shift '
                  f'{row["proxy_delta_relative_to_control"]:.6e}',flush=True)
        except Exception as exc:
            rows.append({'model':model,'status':'Blocked','reason':str(exc)})
            print(f'[{n}/{len(models)}] BLOCKED {model}: {exc}',flush=True)
        atomic(OUTPUT/'summary.json',{'format':'t3-extended-C5-baseline-flavor-proxy-v1',
            'status':'Partial','rows':rows,'physical_matched_C5':False})
    ok=all(r['status']!='Blocked' for r in rows)
    atomic(OUTPUT/'summary.json',{'format':'t3-extended-C5-baseline-flavor-proxy-v1',
        'status':'AllSixBaselineFlavorProxiesCalculated' if ok else 'PartialOrBlocked',
        'rows':rows,'physical_matched_C5':False,
        'physical_extended_neutrino_mass':False})
    print('Report:',OUTPUT/'summary.json')
    return 0 if ok else 4


if __name__=='__main__':raise SystemExit(main())
