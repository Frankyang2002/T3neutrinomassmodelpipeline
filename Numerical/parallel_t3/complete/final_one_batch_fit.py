"""Numerical 18-real-Yukawa objective for a fully sourced C5 callback.

The objective is model-independent. A caller must supply a genuine complete
source-based C5 pipeline; no surrogate triangle or scalar bubble is inserted.
"""
from __future__ import annotations

from typing import Callable
import numpy as np
from scipy.optimize import least_squares

OBS = ('delta_m21_sq_ev2','delta_m3l_sq_ev2','sin2_theta12',
       'sin2_theta13','sin2_theta23')


def residual5(observables:dict, target:dict)->np.ndarray:
    if observables.get('mixing_defined') is False:
        raise ValueError('PMNS undefined')
    center=target['central_values']
    sigma=target['one_sigma_errors']
    pred=dict(observables)
    pred.setdefault('delta_m21_sq_ev2',pred.get('dm21_sq_eV2'))
    pred.setdefault('delta_m3l_sq_ev2',pred.get('dm31_sq_eV2'))
    z=np.asarray([(float(pred[k])-float(center[k]))/float(sigma[k])
                  for k in OBS],float)
    if not np.isfinite(z).all():
        raise ValueError('Nonfinite fit residual')
    return z


def fit_18_real(initial:np.ndarray, evaluate:Callable[[np.ndarray],dict],
                target:dict,*,max_nfev:int=120,
                max_yukawa:float=0.2, penalty:float=1e4)->dict:
    """Fit all 18 real UV y1/y2 entries against 5 Gaussian observables.

    A bad trial receives a penalty and is never treated as a physical point.
    The caller's `evaluate` must perform all requested threshold stages.
    """
    initial=np.asarray(initial,float)
    if initial.shape != (18,) or not np.isfinite(initial).all():
        raise ValueError('Fit requires exactly 18 finite real UV Yukawas')
    if not 1e-5 < max_yukawa < np.sqrt(4*np.pi):
        raise ValueError('Invalid Yukawa bound')
    if max_nfev<1 or penalty<=100 or not np.isfinite(penalty):
        raise ValueError('Invalid optimiser limits')
    records=[]
    best=None
    cache={}
    def residual(x):
        nonlocal best
        key=np.asarray(x,float).tobytes()
        if key in cache: return cache[key].copy()
        try:
            pred=evaluate(np.asarray(x,float))
            values=residual5(pred,target)
            status='CompleteSourceObjectiveEvaluated'
        except Exception as exc:
            values=np.full(5,penalty,dtype=float)
            status='InvalidTrial:'+type(exc).__name__+':'+str(exc)[:160]
        chi2=float(np.dot(values,values))
        row={'evaluation':len(records)+1,'chi2':chi2,'status':status}
        records.append(row)
        if status=='CompleteSourceObjectiveEvaluated' and (
                best is None or chi2 < best['chi2']):
            best={'chi2':chi2,'parameters':np.asarray(x,float).tolist(),
                  'residuals':values.tolist(),'evaluation':len(records)}
        cache[key]=values
        return values.copy()
    lb=np.full(18,-max_yukawa)
    ub=np.full(18,max_yukawa)
    clipped=np.maximum(lb+1e-11,np.minimum(ub-1e-11,initial))
    sol=least_squares(residual,clipped,bounds=(lb,ub),
                      max_nfev=max_nfev,method='trf',
                      x_scale='jac',ftol=1e-8,xtol=1e-8,gtol=1e-8)
    return {
        'status':'SourceDrivenFitCandidate' if best is not None else 'Blocked_NoValidObjective',
        'best':best,'nfev_solver':int(sol.nfev),
        'unique_objective_calls':len(records),
        'optimizer_success':bool(sol.success),
        'optimizer_message':str(sol.message),
        'trials':records,
        'numerical_target':'uncorrelated Gaussian NuFIT five-observable approximation',
        'finite_native_C5_required_in_callback':True,
    }
