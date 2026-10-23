"""Exact quartic invariant-basis map for ordinary T3-B (dS1=dS2=2).

Uses the SAME Matchete real-scalar quartic tensor adapter as production.
Compares it to RGBeta's explicit invariant definitions in T3RGBetaModel.wl:
 (S1†S1)(S2†S2) and sum_A(S1†T^A S1)(S2†T^A S2), T=sigma/2.

Does not assume either Matchete invariant is the singlet or adjoint basis.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from itertools import combinations_with_replacement
import sympy as sp


L1,L2=sp.symbols('lambda12Inv1 lambda12Inv2')
N1,N2=sp.symbols('lambda12 lambda12Adj')


def numerical_cross_quartic_tensors():
    """Real 4th derivatives of RGBeta basis polynomials for the two doublets."""
    rr=sp.symbols('r5:13',real=True)
    S1=sp.Matrix([(rr[0]+sp.I*rr[1])/sp.sqrt(2),
                  (rr[2]+sp.I*rr[3])/sp.sqrt(2)])
    S2=sp.Matrix([(rr[4]+sp.I*rr[5])/sp.sqrt(2),
                  (rr[6]+sp.I*rr[7])/sp.sqrt(2)])
    t=(sp.Matrix([[0,1],[1,0]])/2,
       sp.Matrix([[0,-sp.I],[sp.I,0]])/2,
       sp.diag(sp.Rational(1,2),-sp.Rational(1,2)))
    p1=sp.expand((S1.conjugate().T*S1)[0]*(S2.conjugate().T*S2)[0])
    p2=sp.expand(sum((S1.conjugate().T*g*S1)[0]*(S2.conjugate().T*g*S2)[0] for g in t))
    result=[]
    for p in [p1,p2]:
        tensor={}
        polynomial=sp.Poly(p,*rr)
        from math import factorial
        for powers,coefficient in polynomial.terms():
            if sum(powers)!=4:continue
            key=tuple(sorted(k+5 for k,n in enumerate(powers) for _ in range(n)))
            value=sp.simplify(coefficient*sp.prod(factorial(n) for n in powers))
            if value!=0:tensor[key]=value
        result.append(tensor)
    return result


def matchete_cross_quartic_tensors(quartic_seed,rgbeta):
    from RGE.running.intermediate.ScalarOnlyTensorAdapters import load_and_build_scalar_only_quartic_tensor
    full=load_and_build_scalar_only_quartic_tensor(Path(quartic_seed),Path(rgbeta))
    keys=[key for key in combinations_with_replacement(range(5,13),4)
          if sum(5<=x<=8 for x in key)==2 and sum(9<=x<=12 for x in key)==2]
    output=[]
    for coupling in (L1,L2):
        tensor={}
        for key in keys:
            expr=sp.expand(full[key])
            value=sp.simplify(expr.coeff(coupling))
            if value!=0:tensor[key]=value
        if not tensor:raise ValueError(f'Matchete seed has no {coupling} cross-quartic tensor')
        output.append(tensor)
    return output


def solve_map(mat_basis,num_basis):
    """Solve N_j = sum_i M[i,j] T_i over ALL independent real tensor entries."""
    keys=sorted(set().union(*(set(x) for x in (*mat_basis,*num_basis))))
    A=sp.Matrix([[basis.get(key,0) for basis in mat_basis] for key in keys])
    B=sp.Matrix([[basis.get(key,0) for basis in num_basis] for key in keys])
    if A.rank()!=2:raise ValueError('Matchete scalar quartic basis has rank !=2')
    if B.rank()!=2:raise ValueError('RGBeta scalar quartic basis has rank !=2')
    columns=[]
    for j in range(2):
        solution,free=A.gauss_jordan_solve(B[:,j])
        if free.rows:raise ValueError('Quartic mapping is not unique')
        if any(sp.simplify(e)!=0 for e in A*solution-B[:,j]):
            raise ValueError('Nonzero tensor residual in quartic basis mapping')
        columns.append(solution)
    M=sp.Matrix.hstack(*columns)
    return M,len(keys)


def derive(quartic_seed,rgbeta):
    meta=json.loads(Path(rgbeta).read_text(encoding='utf-8'))['metadata']
    if (meta.get('dS1'),meta.get('dS2'),meta.get('alpha'))!=(2,2,-1):
        raise ValueError('This SU(2) invariant mapping is currently restricted to T3-B alpha=-1')
    mat=matchete_cross_quartic_tensors(quartic_seed,rgbeta)
    num=numerical_cross_quartic_tensors()
    M,n=solve_map(mat,num)
    return {'status':'ExactTensorMap','representation':'T3-B alpha=-1',
        'convention':'V4=-(Matchete L4); RGBeta T^A=sigma^A/2',
        'order_matchete':['lambda12Inv1','lambda12Inv2'],
        'order_rgbeta':['lambda12','lambda12Adj'],
        'matrix':[[sp.sstr(M[i,j]) for j in range(2)] for i in range(2)],
        'lambda12Inv2':sp.sstr(sp.expand((M*sp.Matrix([N1,N2]))[1])),
        'exact_tensor_residuals':0,'tested_real_tensor_entries':n}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--data-dir',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    r=derive(args.data_dir/'eft1_after_F_scalar_quartic_seed.json',args.data_dir/'eft1_rgbeta_rge.json')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(r,indent=2),encoding='utf-8')
    print(json.dumps(r,indent=2))
if __name__=='__main__':main()
