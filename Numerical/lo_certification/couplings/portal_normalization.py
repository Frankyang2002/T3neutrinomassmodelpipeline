"""Independent, exact HHSS portal tensor comparison: native Matchete vs RGBeta.

This is NOT a new RGE.  RGBeta source defines the two renormalisable
portal invariants for every nonsinglet multiplet (one for a singlet):
  L4 = lambda_Hi (H^dag H)(Si^dag Si)
     + lambda_HiAdj (H^dag t^A H)(Si^dag T^A Si).
Their real-field V4 derivatives are constructed algebraically, with the minus
sign from V=-L_int. For each local Matchete quartic column, solve the *exact*
span problem against the RGBeta columns. No candidate aliases are presumed.

The comparison assumes that Matchete's original complex component coordinates
and RGBeta's standard highest-weight SU(2) coordinates are the same. A failure
is recorded, NEVER fitted away with an unknown basis or phase rotation.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from typing import Any

import numpy as np
import sympy as sp


class TensorAuditError(ValueError):
    pass


def _ensure(condition: bool, message: str) -> None:
    if not condition:
        raise TensorAuditError(message)


@lru_cache(maxsize=3)
def su2_generators(dimension: int) -> tuple[sp.Matrix, ...]:
    _ensure(dimension in (1,2,3), 'Only SU(2) multiplets 1,2,3 checked in source')
    j = sp.Rational(dimension-1,2)
    m = [j-i for i in range(dimension)]
    plus = sp.zeros(dimension)
    for col, value in enumerate(m):
        if value+1 in m:
            plus[m.index(value+1), col] = sp.sqrt((j-value)*(j+value+1))
    minus = plus.conjugate().T
    return ((plus+minus)/2, (plus-minus)/(2*sp.I), sp.diag(*m))


def real_bilinear_hessian(generator: sp.Matrix) -> sp.Matrix:
    """For z=(x_R+i x_I)/sqrt2, z^dag G z = x^T K(G) x / 2.

    K is its Hessian in interleaved real coordinates, with NO missing 1/2.
    """
    d=generator.rows
    _ensure(generator.shape==(d,d) and generator==generator.conjugate().T,
            'Only Hermitian bilinear generators supported')
    k=sp.zeros(2*d)
    for p in range(d):
        for q in range(d):
            a=sp.simplify(sp.re(generator[p,q]))
            b=sp.simplify(sp.im(generator[p,q]))
            k[2*p,2*q]=a
            k[2*p,2*q+1]=-b
            k[2*p+1,2*q]=b
            k[2*p+1,2*q+1]=a
    _ensure(k==k.T, 'Real bilinear Hessian must be symmetric')
    return k


def _canonical(a:int,b:int,s:int,t:int) -> tuple[int,int,int,int]:
    return min(a,b),max(a,b),min(s,t),max(s,t)


def _compact(data:dict) -> dict:
    return {k:sp.simplify(v) for k,v in data.items() if sp.simplify(v)!=0}


def rgbeta_portal_columns(dimension: int) -> dict[str,dict[tuple[int,int,int,int],sp.Expr]]:
    """Local heavy multiplet indices only; RGBeta source L4 normalization.

    For V=-lambda (H^dag T H)(S^dag U S),
       d^4 V / dH_A dH_B dS_a dS_b = -lambda K_T[A,B] K_U[a,b].
    """
    kh0=real_bilinear_hessian(sp.eye(2))
    ks0=real_bilinear_hessian(sp.eye(dimension))
    generators=su2_generators(dimension)
    kh=[real_bilinear_hessian(t) for t in su2_generators(2)]
    ks=[real_bilinear_hessian(t) for t in generators]
    result={'singlet':{}, **({'adjoint':{}} if dimension>1 else {})}
    for A in range(4):
        for B in range(A,4):
            for a in range(2*dimension):
                for b in range(a,2*dimension):
                    key=(A,B,a,b)
                    result['singlet'][key]=sp.simplify(-kh0[A,B]*ks0[a,b])
                    if dimension>1:
                        result['adjoint'][key]=sp.simplify(-sum(kh[i][A,B]*ks[i][a,b]
                                                             for i in range(3)))
    return {name:_compact(entries) for name,entries in result.items()}


def native_portal_columns(ch:Any, species:int) -> dict[str,dict[tuple[int,int,int,int],sp.Expr]]:
    _ensure(species in (1,2), 'Species must be 1 or 2')
    d1={'A':1,'B':2,'C':2,'D':3,'E':3}[ch.model[0]]
    offset=0 if species==1 else 2*d1
    dim={'A':(1,3),'B':(2,2),'C':(2,2),'D':(3,1),'E':(3,3)}[ch.model[0]][species-1]
    result={}
    for name,components in ch.by_name.items():
        if not name.startswith(f'lambdaH{species}Inv'):
            continue
        rows={}
        symbol=sp.Symbol(name)
        for key,formula in components.items():
            A,B,s,t=key
            _ensure(offset<=s<offset+2*dim and offset<=t<offset+2*dim,
                    f'{ch.model}/{name}: native tensor crosses scalar multiplets')
            value=sp.simplify(formula.subs(symbol,sp.S.One))
            _ensure(not value.free_symbols and sp.simplify(sp.im(value))==0,
                    f'{ch.model}/{name}: cannot isolate single real native coupling')
            loc=_canonical(A,B,s-offset,t-offset)
            if loc in rows:
                _ensure(sp.simplify(rows[loc]-value)==0,
                        f'{ch.model}/{name}: duplicate inequivalent symmetric components')
            rows[loc]=value
        _ensure(bool(rows),f'{ch.model}/{name}: empty tensor')
        result[name]=_compact(rows)
    return dict(sorted(result.items()))


def _matrix_rows(columns:dict, keys:list[tuple[int,int,int,int]]) -> sp.Matrix:
    names=tuple(columns)
    if not names:
        return sp.zeros(len(keys),0)
    return sp.Matrix([[columns[name].get(key,sp.S.Zero) for name in names] for key in keys])


def _exact_basis_coefficients(rgb:dict,native:dict) -> tuple[sp.Matrix, list[str]]:
    """Returns A with *native-unit tensor* N_j=sum_i RGB_i A_ij.

    This allows independent RGBeta couplings lambda_RGB=A lambda_native,
    but only a square nonsingular A permits reversing the mapping.
    """
    keys=sorted(set().union(*(set(v) for v in (*rgb.values(),*native.values()))))
    R=_matrix_rows(rgb,keys)
    N=_matrix_rows(native,keys)
    if R.rank()!=len(rgb):
        raise TensorAuditError('RGBeta portal tensors unexpectedly linearly dependent')
    pivots=R.T.rref()[1]
    mat=R.extract(pivots,range(R.cols))
    A=mat.inv()*N.extract(pivots,range(N.cols))
    residual=R*A-N
    failures=[]
    for row,key in zip(residual.tolist(),keys):
        if any(sp.simplify(x)!=0 for x in row):
            failures.append(str(key))
    return A,failures


def tensor_to_dense(data:dict,dim:int)->np.ndarray:
    out=np.zeros((4,4,2*dim,2*dim),float)
    for (a,b,s,t),expr in data.items():
        val=float(sp.N(expr,16))
        _ensure(np.isfinite(val), 'Nonfinite exact RGBeta tensor component')
        for A,B in {(a,b),(b,a)}:
            for S,T in {(s,t),(t,s)}:
                out[A,B,S,T]=val
    return out


def real_running_coupling(value:Any,label:str)->float:
    _ensure(type(value) in (int,float) or
            (isinstance(value,dict) and set(value)=={'real','imag'}),
            f'{label}: required real running value absent or invalid')
    if isinstance(value,dict):
        _ensure(abs(float(value['imag']))==0,f'{label}: expected real coupling')
        value=value['real']
    out=float(value)
    _ensure(np.isfinite(out),f'{label}: nonfinite running value')
    return out


def portal_species_audit(ch:Any, species:int)->dict:
    dimension={'A':(1,3),'B':(2,2),'C':(2,2),'D':(3,1),'E':(3,3)}[ch.model[0]][species-1]
    rgb0=rgbeta_portal_columns(dimension)
    prefix=f'lambdaH{species}'
    rgb={prefix if name=='singlet' else prefix+'Adj':v for name,v in rgb0.items()}
    native=native_portal_columns(ch,species)
    A,failures=_exact_basis_coefficients(rgb,native)
    rank=int(A.rank()) if not failures else None
    complete=(not failures and A.rows==A.cols and rank==A.rows)
    status=('exact_full_tensor_span' if complete else
            'exact_partial_tensor_span' if not failures else
            'unresolved_generator_or_complex_scalar_basis')
    transform=[[sp.sstr(sp.simplify(A[i,j])) for j in range(A.cols)] for i in range(A.rows)]
    result={'species':f'S{species}','dimension':dimension,'status':status,
            'basis_assumption':'same highest-weight complex SU2 component frame in RGBeta and Matchete',
            'potential_convention':'V4=lambda_abcd phi_a phi_b phi_c phi_d/4!',
            'rgbeta_source_potential_sign':'V4 = - L4 interaction',
            'native_symbols':list(native),'rgbeta_symbols':list(rgb),
            'native_unit_to_rgbeta_invariant_coefficients':transform,
            'native_tensor_outside_rgbeta_span_components':failures[:12],
            'rgbeta_independent_portal_directions':len(rgb),
            'native_independent_directions_in_rgbeta_span':rank,
            'rgbeta_portal_directions_missing_from_native_span':
                len(rgb)-rank if rank is not None else None,
            'normalization_conversion_physically_usable_for_portal_subsector':complete,
            'full_EFT1_source_completeness_certified':False}
    if complete:
        inv=A.inv()
        result['rgbeta_couplings_to_native_couplings']={
           native_name:{rgb_name:sp.sstr(sp.simplify(inv[j,i]))
                        for i,rgb_name in enumerate(rgb)}
           for j,native_name in enumerate(native)}
    return result


def construct_rgbeta_portal_background(ch:Any, f:dict, s:dict)->tuple[dict,dict]:
    """Generate **portal-only** RGBeta tensors at both scales, in the assumed
    original highest-weight real basis. It omits mixing and extra invariants.
    Not a complete seven-term or joint-scalar solver input.
    """
    species={i:portal_species_audit(ch,i) for i in (1,2)}
    d1,d2={'A':(1,3),'B':(2,2),'C':(2,2),'D':(3,1),'E':(3,3)}[ch.model[0]]
    n=2*(d1+d2)
    arrays={}
    missing={}
    for scale,state in (('muF',f),('muS',s)):
        tensor_sum=np.zeros((4,4,n,n),float)
        individual=[]; labels=[]; missing[scale]=[]
        for i,d in ((1,d1),(2,d2)):
            offset=0 if i==1 else 2*d1
            for basis_name,local_tensor in rgbeta_portal_columns(d).items():
                symbol=f'lambdaH{i}'+('Adj' if basis_name=='adjoint' else '')
                if symbol not in state:
                    missing[scale].append(symbol)
                    continue
                value=real_running_coupling(state[symbol],symbol)
                term=np.zeros((4,4,n,n),float)
                term[:,:,offset:offset+2*d,offset:offset+2*d]=value*tensor_to_dense(local_tensor,d)
                labels.append(symbol);individual.append(term);tensor_sum+=term
        arrays[f'{scale}_portal_channel_tensors']=np.stack(individual,axis=0) if individual else np.zeros((0,4,4,n,n))
        arrays[f'{scale}_portal_sum_tensor']=tensor_sum
        arrays[f'{scale}_portal_labels']=np.asarray(labels,dtype='U64')
    detail={'source_basis':'original S1/S2 interleaved canonical real frame (assumed RGBeta highest weight)',
            'tensor_shape':[4,4,n,n],
            'physical_completeness':False,
            'muF_missing_rgbeta_portals':missing['muF'],
            'muS_missing_rgbeta_portals':missing['muS'],
            'portal_subsectors':list(species.values()),
            'omitted_independent_interactions':['lambdaT3','all extra HHSS symbols in registry'],
            'safe_for_complete_seven_term_or_joint_matching':False}
    return detail,arrays


def numerically_map_full_portal_subsector(ch:Any,species:int,running:dict)->dict:
    """Invert an exactly verified square native-to-RGBeta tensor map.

    No solution is emitted for incomplete basis, even if an underdetermined
    numerical least-squares solution exists. This avoids silently removing
    any independent quartic. No other interaction is assigned a zero value.
    """
    audit=portal_species_audit(ch,species)
    if audit['status']!='exact_full_tensor_span':
        return {'status':'blocked_tensor_basis_not_exact_full_rank',
                'native_couplings':None,
                'missing_numeric_couplings':list(audit['rgbeta_symbols'])}
    mapping=audit['rgbeta_couplings_to_native_couplings']
    needed=audit['rgbeta_symbols']
    absent=sorted(set(needed)-set(running))
    if absent:
        return {'status':'blocked_missing_running_rgbeta_couplings',
                'native_couplings':None,'missing_numeric_couplings':absent}
    values={name:real_running_coupling(running[name],name) for name in needed}
    native={}
    for native_name,terms in mapping.items():
        z=sum(sp.sympify(factor)*values[rgb_name] for rgb_name,factor in terms.items())
        z=complex(sp.N(z,16))
        _ensure(abs(z.imag) <= 1e-12*max(1.,abs(z.real)),
                'Native scalar portal coupling is unexpectedly complex')
        native[native_name]=float(z.real)
    # Independent reconstruction of two V4 tensors and their HHSS fourth
    # derivatives for the *portal subsector only*. Use exactly the same
    # canonical real-field coordinates and source native CG expressions.
    dim=audit['dimension']
    observed=np.zeros((4,4,2*dim,2*dim),float)
    for native_name,column in native_portal_columns(ch,species).items():
        observed+=native[native_name]*tensor_to_dense(column,dim)
    expected=np.zeros_like(observed)
    suffix=('Adj','')
    for name,column in rgbeta_portal_columns(dim).items():
        key=f'lambdaH{species}'+('Adj' if name=='adjoint' else '')
        expected+=values[key]*tensor_to_dense(column,dim)
    scale=max(float(np.max(abs(expected))),1.0)
    residual=float(np.max(abs(expected-observed))/scale)
    _ensure(residual<=1e-10,
            f'Exact portal map fails independent numerical V4 tensor check: {residual}')
    return {'status':'exact_portal_only_numeric_conversion_verified',
            'native_couplings':native,'rgbeta_couplings':values,
            'normalized_real_HHSS_tensor_relative_residual':residual,
            'missing_numeric_couplings':[],
            'full_potential_or_physical_C5_certified':False}
