"""Numerical L/eC one-generation gauge-index closure of the project psi2phi2 RGE.

The Wilson coefficients carry six independent *physical* lepton-flavor pairs;
internal SM Yukawa and gauge indices retain the published one-generation RGE
convention.  No full three-generation charged-lepton anomalous dimension is
claimed.  Exactly the seven terms of RGE.general.WilsonTensorRGE are evaluated.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from math import factorial, pi
from pathlib import Path
from typing import Any, Mapping
import json

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

LOOP = 16.0*pi*pi
PAIRS = ((0,0),(0,1),(0,2),(1,1),(1,2),(2,2))
BASE_COUPLINGS = {
    'λ':'lambdaH',
    'lambdaS1Inv1':'lambdaS1', 'lambdaS2Inv1':'lambdaS2',
    'lambdaH1Inv1':'lambdaH1', 'lambdaH1Inv2':'lambdaH1Adj',
    'lambdaH2Inv1':'lambdaH2', 'lambdaH2Inv2':'lambdaH2Adj',
    'lambda12Inv1':'lambda12', 'lambda12Inv2':'lambda12Adj',
    'lambdaT3':'lambdaT3',
}


def require(value: bool, msg: str) -> None:
    if not value:
        raise ValueError('B/C coupled Wilson: '+msg)


def finite_complex(x: Any, name: str) -> complex:
    z=complex(x)
    require(np.isfinite(z.real) and np.isfinite(z.imag), name+' nonfinite')
    return z


def flavor_pairs(masses: np.ndarray, y1: np.ndarray, y2: np.ndarray) -> np.ndarray:
    """F11, F12(ordered Majorana sum), F22, each of shape (6,).

    F12 contains BOTH cross orientations, matching the established project
    hierarchy adapter (WeinbergFlavorMatching.build_majorana_c5_flavor_matrix)
    and DirectWilsonSequential.pair_matrix. The old average introduced an
    exact 1/2 in the direct-source parity check. The same-Yukawa F11/F22
    channels are intentionally unchanged (independent Stage17 normalization).
    This is parity with the established candidate convention, not an
    independent physical matching-normalization certificate.
    """
    m=np.asarray(masses,dtype=float)
    a=np.conjugate(np.asarray(y1,dtype=complex))
    b=np.conjugate(np.asarray(y2,dtype=complex))
    require(m.shape==(3,) and np.all(m>0) and np.isfinite(m).all(), 'three positive heavy masses required')
    require(a.shape==b.shape==(3,3) and np.isfinite(a).all() and np.isfinite(b).all(), 'invalid Yukawa flavor data')
    f11=(a/m[None,:])@a.T
    # Preserve BOTH ordered cross orientations used by the project's existing
    # full-flavor Majorana C5 lift. Do not average here: it halves the direct
    # Wilson beta relative to the separately saved hierarchical calculation.
    f12=(a/m[None,:])@b.T+(b/m[None,:])@a.T
    f22=(b/m[None,:])@b.T
    return np.asarray([[mat[p,q] for p,q in PAIRS] for mat in (f11,f12,f22)],dtype=complex)


def independent_tensor(wilson: Any, factors: np.ndarray, *, nf: int=3, ns: int=12) -> np.ndarray:
    """Lift the published PL tree tensor into full flavor; three polynomial probes.

    A one-generation product conj(y1)*conj(y2) is lifted using the ordered
    Majorana cross sum (F12) in the historical project matching convention.
    """
    require(factors.shape==(3,len(PAIRS)), 'flavor basis mismatch')
    c=np.zeros((nf,nf,ns,ns,len(PAIRS)),dtype=complex)
    for (i,j,aa,bb),expr in wilson.nonzero_items():
        require(1<=i<=nf and 1<=j<=nf and 1<=aa<=ns and 1<=bb<=ns,'seed support outside L/eC closure')
        expression=sp.sympify(expr)
        # ScalarOnlyTensorAdapters._term_prefactor uses Symbol('MF', nonzero=True).
        # SymPy treats it as different from Symbol('MF'), so constructing a
        # fresh unassumed mass symbol leaves 1/MF unevaluated in the probes.
        # Bind by name to the *actual* symbols present in each seed component.
        by_name={}
        for symbol in expression.free_symbols:
            by_name.setdefault(symbol.name,[]).append(symbol)
        require(set(by_name)<= {'y1','y2','MF'} and 'MF' in by_name,
                'unresolved PL Yukawa/fermion mass in seed: '+str(expression))
        require(all(len(symbols)==1 for symbols in by_name.values()),
                'ambiguous symbolic assumptions in PL Wilson seed: '+str(expression))
        a=by_name.get('y1',[sp.Symbol('y1')])[0]
        b=by_name.get('y2',[sp.Symbol('y2')])[0]
        mf=by_name['MF'][0]
        substitutions=[{a:1,b:0,mf:1},{a:0,b:1,mf:1},{a:1,b:1,mf:1}]
        v=[]
        for params in substitutions:
            x=sp.N(expression.subs(params),16)
            require(not x.free_symbols,'unresolved PL Yukawa/fermion mass in seed: '+str(x))
            v.append(complex(x))
        co=np.asarray([v[0],v[2]-v[0]-v[1],v[1]],dtype=complex)
        c[i-1,j-1,aa-1,bb-1]=co@factors
    require(np.allclose(c,c.swapaxes(0,1),rtol=1e-11,atol=1e-25),'seed fermion symmetry')
    require(np.allclose(c,c.swapaxes(2,3),rtol=1e-11,atol=1e-25),'seed scalar symmetry')
    return c


def _symmetric_lookup(lookup: Mapping[tuple[int,int,int,int],Any], ns:int=12) -> np.ndarray:
    tensor=np.zeros((ns,ns,ns,ns),dtype=complex)
    for key,value in lookup.items():
        for perm in set(__import__('itertools').permutations(tuple(k-1 for k in key))):
            tensor[perm]=complex(value)
    return tensor


def _sympy_expression_lambdas(expressions: Mapping[tuple[int,int,int,int],Any]):
    keys=sorted(expressions)
    expr=[sp.sympify(expressions[k]) for k in keys]
    symbols=sorted(set().union(*(e.free_symbols for e in expr)),key=str)
    unknown={str(s) for s in symbols}-set(BASE_COUPLINGS)
    require(not unknown, 'unmapped Matchete baseline quartic parameters '+repr(unknown))
    fn=sp.lambdify(symbols,expr,modules='numpy',cse=True)
    return keys,symbols,fn


def _matchete_complex(text: Any) -> complex:
    from RGE.running.intermediate.IntermediateMatcheteParsing import parse_mathematica_scalar
    return complex(sp.N(parse_mathematica_scalar(str(text)),17))


def extra_matchete_templates(payload: Mapping[str,Any], model: str) -> dict[str, dict[tuple[int,int,int,int],complex]]:
    """Expand saved stage22 CGs in the (R, I) real-component ordering.

    Returns polynomial coefficients for lambda_M * monomial, before H.c.
    """
    require(payload.get('format')=='stage22-exact-tensors-v1','stage22 tensor contract changed')
    tensors=payload.get('models',{}).get(model,{}).get('tensors')
    require(isinstance(tensors,list) and len(tensors)==7,'expected exactly seven additional B/C quartics')
    answer={}
    for row in tensors:
        symbol=row['symbol'];dims=row['dimensions'];species=row['species']
        require(len(species)==len(dims)==4 and all(int(d)==2 for d in dims),'B/C quartic component dimensions changed')
        require(row.get('potential_prefactor')=='1','stage22 quartic prefactor convention changed')
        vals=row['components']
        require(len(vals)==int(np.prod(dims)),'incomplete stage22 tensor data')
        poly={}
        for flat,raw in enumerate(vals):
            z=_matchete_complex(raw)
            if abs(z)<1e-14:continue
            coords=np.unravel_index(flat,dims)
            opts=[]
            for leg,k in zip(species,coords):
                base='H' if leg.startswith('H') else ('S1' if leg.startswith('S1') else 'S2')
                start={'H':0,'S1':4,'S2':8}[base]
                real=start+2*k+1
                imag=real+1
                sign=-1j if leg.endswith('dag') else 1j
                opts.append(((real,1/np.sqrt(2)),(imag,sign/np.sqrt(2))))
            for bits in product((0,1),repeat=4):
                chosen=[opts[i][v] for i,v in enumerate(bits)]
                key=tuple(sorted(t[0] for t in chosen))
                prod=z*np.prod([t[1] for t in chosen])
                poly[key]=poly.get(key,0j)+prod
        require(symbol not in answer,'duplicate stage22 quartic '+symbol)
        answer[symbol]=poly
    return answer


def matchete_from_rgbeta(model: str, handoff:Mapping[str,Any], state:Any) -> dict[str,complex]:
    """Apply the *candidate*, not physically certified, stage29 signed map."""
    require(handoff.get('format')=='stage29-conditional-uv-quartic-handoff-v1','stage29 mapping format changed')
    require(handoff.get('physical_coupling_mapping_certified') is False,'stage29 certification changed: review required')
    meta=handoff['models'][model]
    result={}
    for block in meta['blocks']:
        rhs=[]
        for symbol in block['rgbeta_symbols']:
            if symbol.startswith('lambdaT3General'):
                v=state.extra[symbol]
            else:
                v=getattr(state.baseline,symbol)
            rhs.append(finite_complex(v,symbol))
        raw=block['signed_map_exact']
        mat=np.array([[complex(sp.N(sp.sympify(str(t)),16)) for t in r] for r in raw],dtype=complex)
        require(mat.shape==(len(rhs),len(block['matchete_symbols'])) and mat.shape[0]==mat.shape[1], 'unexpected stage29 map matrix')
        cond=np.linalg.cond(mat)
        require(np.isfinite(cond) and cond<1e8,'ill-conditioned stage29 map')
        values=np.linalg.solve(mat,np.array(rhs))
        for symbol,v in zip(block['matchete_symbols'],values):
            require(symbol not in result,'duplicate stage29 mapped symbol')
            result[symbol]=finite_complex(v,symbol)
    return result


@dataclass
class QuarticEvaluator:
    baseline: Any
    extra_templates: dict
    handoff: Mapping[str,Any]
    model: str
    keys: list
    symbols: list
    function: Any

    @classmethod
    def build(cls,base_tensor:Any, stage22:Mapping[str,Any],handoff:Mapping[str,Any],model:str):
        keys,syms,fn=_sympy_expression_lambdas(dict(base_tensor.nonzero_items()))
        return cls(base_tensor,extra_matchete_templates(stage22,model),handoff,model,keys,syms,fn)

    def evaluate(self,state:Any, *, include_extra:bool=True, omit_active:bool=False)->np.ndarray:
        base=state.baseline
        env={s:finite_complex(getattr(base.sm,'lambdaH') if str(s)=='λ' else getattr(base,BASE_COUPLINGS[str(s)]),str(s)) for s in self.symbols}
        raw=np.atleast_1d(self.function(*[env[s] for s in self.symbols]))
        require(len(raw)==len(self.keys),'baseline quartic evaluation returned wrong size')
        coefficients={k:complex(v) for k,v in zip(self.keys,raw)}
        if include_extra:
            mapped=matchete_from_rgbeta(self.model,self.handoff,state)
            require(set(mapped)==set(self.extra_templates),'stage22/stage29 additional quartic coverage mismatch')
            for coupling,poly in self.extra_templates.items():
                # Certified B/C active two-Higgs quartics are reconstructed
                # from the Stage17 RGBeta potential, not the provisional
                # Matchete-group mapping.
                if coupling in ('lambdaT3ExtraQ020020I1','lambdaT3ExtraQ020200I1'):
                    continue
                # Saved Stage22 tensor is the Matchete Lagrangian tensor.
                # V4=-L4; stage29 map contains the same Lagrangian sign.
                lam=-mapped[coupling]
                for key,value in poly.items():
                    mult=np.prod([factorial(key.count(k)) for k in set(key)])
                    coefficients[key]=coefficients.get(key,0j)+2*np.real(lam*value)*mult
            if not omit_active:
                from Numerical.experimental.GeneralT3ExtraEFT1Mixing import real_quartic_entries
                if not hasattr(self,'_active_functions'):
                    symbols=[sp.Symbol('lambdaHHdagS2S2'),sp.Symbol('lambdaHHdagS1barS1bar')]
                    source=real_quartic_entries()
                    self._active_keys=sorted(source)
                    self._active_functions=sp.lambdify(symbols,
                              [source[k] for k in self._active_keys], modules='numpy',cse=True)
                v=self._active_functions(base.lambdaHHdagS2S2,base.lambdaHHdagS1barS1bar)
                for key,value in zip(self._active_keys,np.atleast_1d(v)):
                    coefficients[key]=coefficients.get(key,0j)+complex(value)
        tensor=_symmetric_lookup(coefficients)
        require(np.max(np.abs(np.imag(tensor)))<2e-9,'quartic potential is not real')
        return np.asarray(np.real(tensor),dtype=float)


def build_static_kinetics(context:Any)->dict:
    """Numerical generators and Yukawa basis from the project's SM master."""
    from RGE.general.GaugeGenerators import quadratic_casimir_matrix
    nf=context.inputs.fermion_dimension
    ns=context.model.total_real_scalar_dimension
    require(ns==12 and nf>=15,'B/C scalar and SM fermion dimensions changed')
    ids=[context.fermion_basis.global_index('L',1)-1,
         context.fermion_basis.global_index('L',2)-1,
         context.fermion_basis.global_index('eC',1)-1]
    require(len(set(ids))==3,'L/eC closure indexing changed')
    Y=[]
    syms=tuple(sp.Symbol(z) for z in ('ye','yu','yd'))
    for symbol in syms:
        array=np.zeros((nf,nf,ns),dtype=complex)
        for i in range(nf):
            for j in range(nf):
                for a in range(ns):
                    val=context.inputs.yukawa(i+1,j+1,a+1)
                    if val!=0:
                        x=sp.N(val.subs({x:1 if x==symbol else 0 for x in syms}),16)
                        require(not x.free_symbols,'unrecognised one-generation SM Yukawa')
                        array[i,j,a]=complex(x)
        Y.append(array)
    gauge=[]
    for sector in context.inputs.gauge_sectors:
        name=str(sector.coupling)
        require(name in ('g1','g2'),'unexpected gauge generator convention '+name)
        ts=[np.asarray(mat,dtype=complex) for mat in sector.scalar_generators]
        tf=[np.asarray(mat,dtype=complex) for mat in sector.fermion_generators]
        cs=np.asarray(quadratic_casimir_matrix(sector.scalar_generators),dtype=complex)
        cf=np.asarray(quadratic_casimir_matrix(sector.fermion_generators),dtype=complex)
        require(all(t.shape==(ns,ns) for t in ts) and all(t.shape==(nf,nf) for t in tf),'gauge tensor dimensions wrong')
        gauge.append((name,ts,tf,cs,cf))
    return {'indices':np.array(ids,dtype=int),'Ybasis':np.asarray(Y),'gauge':gauge,'nf':nf,'ns':ns}


def kernels(state:Any,static:dict,quartic:np.ndarray)->dict[str,np.ndarray]:
    """Background-dependent inputs for the seven-term Wilson master."""
    sm=state.baseline.sm
    # The published gauge-index master has a one-generation SM Yukawa tensor.
    # Match the Higgs anomalous-dimension trace using Frobenius norms, while
    # retaining full lepton flavor on the BSM Wilson coefficient itself.
    sm_norms=[]
    for name in ('ye','yu','yd'):
        matrix=np.asarray(getattr(sm,name),dtype=complex)
        require(matrix.shape==(3,3) and np.isfinite(matrix).all(), 'SM '+name+' 3x3 state invalid')
        sm_norms.append(float(np.linalg.norm(matrix,'fro')))
    ye,yu,yd=sm_norms
    Y=np.einsum('x,xija->ija',np.array([ye,yu,yd]),static['Ybasis'],optimize=True)
    ids=static['indices'];n=3
    y=Y[np.ix_(ids,ids,np.arange(12))]
    y_iall=Y[ids,:,:]
    y_allj=Y[:,ids,:]
    gamma_s=np.zeros((12,12),dtype=complex)
    gamma_f=np.zeros((n,n),dtype=complex)
    scalar_gauge=np.zeros((12,12,12,12),dtype=complex)
    mixed_gauge=np.zeros((n,n,12,12),dtype=complex)
    for name,ts,tf,cs,cf in static['gauge']:
        g=finite_complex(sm.gY if name=='g1' else sm.g2,name)
        g2=g*g
        gamma_s+=-4*g2*cs
        gamma_f+=-3*g2*cf[np.ix_(ids,ids)]
        for t_s,t_f in zip(ts,tf):
            scalar_gauge+=2*g2*np.einsum('ac,bd->abcd',t_s,t_s,optimize=True)
            mixed_gauge+=4*g2*np.einsum('lj,bc->jlbc',t_f[np.ix_(ids,ids)],t_s,optimize=True)
    gamma_s+=0.5*(np.einsum('jia,jib->ab',Y,np.conj(Y),optimize=True)+np.einsum('jib,jia->ab',Y,np.conj(Y),optimize=True))
    gamma_f+=0.5*np.einsum('ika,jka->ij',y_iall,np.conj(y_iall),optimize=True)
    yukmix=2*np.einsum('jkc,klb->jlbc',y_iall,np.conj(y_allj),optimize=True)
    yukmix+=np.einsum('jkb,klc->jlbc',y_iall,np.conj(y_allj),optimize=True)
    return {'y':y,'quartic':quartic,'scalarg':scalar_gauge,
            'mixed':yukmix+mixed_gauge,'gammas':gamma_s,'gammaf':gamma_f}


def rhs_master(C:np.ndarray, k:Mapping[str,np.ndarray])->np.ndarray:
    """All seven terms of RGE.general.WilsonTensorRGE, same pair swaps."""
    c=np.asarray(C,dtype=complex)
    require(c.ndim==5 and c.shape[:4]==(3,3,12,12),'Wilson tensor must be (3,3,12,12,nflavor)')
    y=k['y']
    term=np.einsum('ijd,kld,klabf->ijabf',y,np.conj(y),c,optimize=True)
    term+=np.einsum('abcd,ijcdf->ijabf',k['quartic']-k['scalarg'],c,optimize=True)
    mixed=np.einsum('jlbc,ilacf->ijabf',k['mixed'],c,optimize=True)
    term+=mixed+mixed.swapaxes(0,1)+mixed.swapaxes(2,3)+mixed.swapaxes(0,1).swapaxes(2,3)
    crossed=np.einsum('jla,ikc,lkbcf->ijabf',y,y,c,optimize=True)
    term+=crossed+crossed.swapaxes(0,1).swapaxes(2,3)
    cc=np.conj(c)
    term+=2*np.einsum('ild,jkd,klabf->ijabf',y,y,cc,optimize=True)
    term+=np.einsum('ijd,kld,klabf->ijabf',y,y,cc,optimize=True)
    gs=np.einsum('bd,ijadf->ijabf',k['gammas'],c,optimize=True)
    term+=gs+gs.swapaxes(2,3)
    gf=np.einsum('jk,ikabf->ijabf',k['gammaf'],c,optimize=True)
    term+=gf+gf.swapaxes(0,1)
    require(np.isfinite(term).all(),'Wilson beta nonfinite')
    return term


def weinberg_projection(C:np.ndarray,template:Mapping,ids:np.ndarray,*,threshold:float=2e-5)->tuple[np.ndarray,float]:
    """Check ALL 30 independent LLHH entries against the Weinberg template.

    The template mapping is sparse: entries absent from it must also vanish.
    Testing only nonzero template entries would accept a gauge-breaking tensor.
    """
    from itertools import combinations_with_replacement
    require(all(i-1 in ids for i in (1,2)), 'L indices missing from closure')
    T=[];A=[]
    for i,j in combinations_with_replacement((1,2),2):
        ii=int(np.flatnonzero(ids==i-1)[0]);jj=int(np.flatnonzero(ids==j-1)[0])
        for a,b in combinations_with_replacement((1,2,3,4),2):
            T.append(complex(template.get((i,j,a,b),0)))
            A.append(C[ii,jj,a-1,b-1,:])
    require(len(T)==30,'expected exactly 30 independent LLHH tensors')
    T=np.asarray(T,dtype=complex)
    A=np.asarray(A,dtype=complex)
    denom=np.vdot(T,T)
    require(denom.real>0,'empty Weinberg reference template')
    projected=(np.conj(T)[:,None]*A).sum(axis=0)/denom
    residual=np.linalg.norm(A-T[:,None]*projected[None,:])/max(np.linalg.norm(A),1e-300)
    require(np.isfinite(residual) and residual<=threshold,
            'Higgs output not a pure Weinberg SU(2) tensor: residual='+str(residual))
    return projected,float(residual)


def flavor_to_matrix(values:np.ndarray)->np.ndarray:
    v=np.asarray(values,dtype=complex)
    require(v.shape==(6,),'expected six lepton flavor pairs')
    m=np.zeros((3,3),dtype=complex)
    for k,(p,q) in enumerate(PAIRS):m[p,q]=m[q,p]=v[k]
    return m


def interpolate_kernels(log_mu:np.ndarray, components:list[Mapping[str,np.ndarray]], t:float)->dict:
    """Linear interpolation in log(mu); endpoints fail closed."""
    grid=np.asarray(log_mu,dtype=float)
    require(np.all(np.diff(grid)<0),'source scales must strictly descend')
    require(grid[-1]-1e-10<=t<=grid[0]+1e-10,'requested scale outside EFT1 background')
    t=float(np.clip(t,grid[-1],grid[0]))
    q=np.searchsorted(-grid,-t,side='right')
    if q==0:return components[0]
    if q>=len(grid):return components[-1]
    a=q-1;b=q
    weight=(t-grid[a])/(grid[b]-grid[a])
    return {k:(1-weight)*components[a][k]+weight*components[b][k] for k in components[0]}


def integrate_coupled(initial:np.ndarray, mu:np.ndarray, kernels_on_grid:list,
                      *,rtol:float=5e-7,atol:float=1e-24):
    """Integrate a real-linear complex Wilson system (C and conjugate C)."""
    t=np.log(np.asarray(mu,dtype=float));require(len(t)==len(kernels_on_grid) and len(t)>=9,'invalid EFT1 Wilson background grid')
    shape=initial.shape
    def fun(x,vector):
        flat=vector.reshape((2,)+shape)
        c=flat[0]+1j*flat[1]
        k=interpolate_kernels(t,kernels_on_grid,x)
        b=rhs_master(c,k)/LOOP
        return np.stack([b.real,b.imag]).ravel()
    init=np.stack([initial.real,initial.imag]).ravel()
    output=solve_ivp(fun,(float(t[0]),float(t[-1])),init,method='RK45',rtol=rtol,atol=atol,
                     t_eval=t,max_step=0.15)
    require(output.success and np.isfinite(output.y).all(),'Wilson ODE integration failed: '+str(output.message))
    trajectory=output.y.reshape((2,)+shape+(len(t),))
    return trajectory[0]+1j*trajectory[1],output
