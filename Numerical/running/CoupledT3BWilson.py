"""Standalone one-generation, one-loop RGE-resummed C12/C5 evolution.

Requires an exact SU(2) cross-quartic mapping produced by
validation.DeriveT3BCrossQuarticMap. No changes to production matching.
The background is a JSON array of samples of the already-running scalar EFT;
all Wilson coefficients are dimension -1 and complex. Interpolation in ln(mu).

This is NOT yet the full three-generation physical neutrino-mass pipeline.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

FACTOR=16*np.pi**2

def _complex(z):
    if isinstance(z,(list,tuple)) and len(z)==2:return complex(float(z[0]),float(z[1]))
    if isinstance(z,(int,float,complex)):return complex(z)
    if isinstance(z,str):return complex(z.replace(' ','').replace('I','j'))
    raise ValueError(f'Invalid complex scalar {z!r}; use [Re, Im]')

def _pair(z):
    w=complex(z);return [float(w.real),float(w.imag)]


def build_matrix_evaluator(validation_json,mapping_json):
    check=json.loads(Path(validation_json).read_text(encoding='utf-8'))
    mapping=json.loads(Path(mapping_json).read_text(encoding='utf-8'))
    if check.get('status')!='ChecksPassed':raise ValueError('Symbolic validation not passed')
    if mapping.get('status')!='ExactTensorMap' or mapping.get('exact_tensor_residuals')!=0:
        raise ValueError('Missing exact quartic invariant-basis map')
    if mapping.get('representation')!='T3-B alpha=-1':raise ValueError('Unexpected model')
    names=('gamma_12_12','gamma_12_5','gamma_5_12','gamma_5_5')
    symbolic=[]
    lamH=sp.Symbol('lambdaH')
    lambdas=sp.Matrix([[sp.sympify(x) for x in row] for row in mapping['matrix']])
    linear=sp.expand((lambdas*sp.Matrix(sp.symbols('lambda12 lambda12Adj')))[1])
    for name in names:
        expr=sp.sympify(check['anomalous_dimension_matrix'][name])
        expr=sp.expand(expr.subs(sp.Symbol('λ'),lamH).subs(sp.Symbol('lambda12Inv2'),linear))
        symbolic.append(expr)
    symbols=sorted(set().union(*(e.free_symbols for e in symbolic)),key=lambda s:s.name)
    allowed={'g2','lambdaH','lambda12','lambda12Adj','lambdaT3','yu','yd','ye'}
    unknown={x.name for x in symbols}-allowed
    if unknown:raise ValueError(f'Unmapped symbols: {unknown}')
    # SymPy conjugates of real-coupling symbols are harmless; complex couplings
    # require lambdify to preserve numpy.conjugate.
    f=sp.lambdify(symbols,symbolic,modules='numpy')
    def mat(couplings):
        args=[_complex(couplings[k.name]) for k in symbols]
        g=np.asarray(f(*args),dtype=complex).reshape((2,2))
        if not np.all(np.isfinite(g)):raise ValueError('Nonfinite mixing matrix')
        return g
    return mat,{n:sp.sstr(v) for n,v in zip(names,symbolic)}


def solve(validation_json,mapping_json,background_json,s0,kappa0,rtol=1e-9,atol=1e-12):
    matrix,expressions=build_matrix_evaluator(validation_json,mapping_json)
    data=json.loads(Path(background_json).read_text(encoding='utf-8'))
    samples=data.get('samples',data) if isinstance(data,dict) else data
    if isinstance(data,dict) and data.get('scope') not in (None,'one-generation diagnostic background'):
        raise ValueError('Unsupported background scope')
    if len(samples)<2:raise ValueError('Need at least two background samples')
    # Each sample: {mu_gev,g2,lambdaH,lambda12,lambda12Adj,lambdaT3,yu,yd,ye}.
    samples=sorted(samples,key=lambda p:float(p['mu_gev']))
    mus=np.asarray([float(s['mu_gev']) for s in samples]);ts=np.log(mus)
    if any(not np.isfinite(t) for t in ts) or np.any(np.diff(ts)<=0):
        raise ValueError('Background scales must be positive and distinct')
    required=('g2','lambdaH','lambda12','lambda12Adj','lambdaT3','yu','yd','ye')
    for row in samples:
        for key in required:
            if key not in row:raise ValueError(f'Missing running coupling {key!r}')
    lookup={k:np.asarray([_complex(p[k]) for p in samples]) for k in required}
    def current(t):
        return {k:complex(np.interp(t,ts,a.real),np.interp(t,ts,a.imag))
                for k,a in lookup.items()}
    # Integrate from the highest scale to the lowest scale.
    t0,t1=ts[-1],ts[0]
    init=np.asarray([_complex(s0),_complex(kappa0)],dtype=complex)
    def rhs(t,y):return matrix(current(t))@y/FACTOR
    result=solve_ivp(rhs,(t0,t1),init,method='DOP853',rtol=rtol,atol=atol,
                     dense_output=False)
    if not result.success:raise RuntimeError(result.message)
    leading=init+(t1-t0)*(matrix(current(t0))@init)/FACTOR
    return {'status':'Success','method':'one_generation_coupled_one_loop_RGE',
        'model':'T3-B alpha=-1','mu_high_gev':float(mus[-1]),
        'mu_low_gev':float(mus[0]),'initial':{'s':_pair(init[0]),'kappa':_pair(init[1])},
        'final':{'s':_pair(result.y[0,-1]),'kappa':_pair(result.y[1,-1])},
        'frozen_couplings_fixed_order_LL':{'s':_pair(leading[0]),'kappa':_pair(leading[1])},
        'nfev':int(result.nfev),'matrix_expressions':expressions,
        'trajectory':[{'mu_gev':float(np.exp(t)),'s':_pair(result.y[0,i]),
                       'kappa':_pair(result.y[1,i])} for i,t in enumerate(result.t)],
        'limitations':['one generation','background samples supplied externally',
                       'no scalar-threshold matching','not a production neutrino observable']}

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--validation',type=Path,required=True)
    p.add_argument('--mapping',type=Path,required=True)
    p.add_argument('--background',type=Path,required=True)
    p.add_argument('--s0',required=True,help='[Re,Im] as JSON')
    p.add_argument('--kappa0',required=True,help='[Re,Im] as JSON')
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    result=solve(a.validation,a.mapping,a.background,json.loads(a.s0),json.loads(a.kappa0))
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps({k:result[k] for k in ('status','final','frozen_couplings_fixed_order_LL','nfev')},indent=2))
if __name__=='__main__':main()
