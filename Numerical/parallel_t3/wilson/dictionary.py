"""Extend Batch1 numeric polynomial mapping with repeated/reused coupling slots.

A repeated symbol may have several native monomial records; their *sum* is
one coefficient column, not several independent couplings. This prevents
silent rank inflation from the 12 reused RGBeta symbol slots.
"""
from __future__ import annotations
from collections import OrderedDict
from typing import Mapping,Sequence
import math,itertools
import numpy as np


def group_fit(mat_records:Sequence[Mapping],rgb_records:Sequence[Mapping],
              *,samples:int=256,seed:int=1907,tol:float=1e-8) -> dict:
    from Numerical.basis.T3Batch1TensorMap import Polynomial
    pm=[Polynomial.from_record(x) for x in mat_records]
    pr=[Polynomial.from_record(x) for x in rgb_records]
    if not pm or not pr:raise ValueError('Blocked_EmptyNativePolynomialBasis')
    dm={};
    for p in itertools.chain(pm,pr):
        for (spec,_),d in zip(p.species,p.tensor.shape):
            if spec in dm and dm[spec]!=d:raise ValueError('Blocked_ComponentFrameDimensionMismatch')
            dm[spec]=d
    def matrices(polys,fields):
        groups=OrderedDict()
        for p in polys:groups.setdefault(p.symbol,[]).append(p)
        arr=np.array([[sum(p.evaluate(x) for p in cols) for cols in groups.values()]
                      for x in fields],float)
        return arr,list(groups)
    rng=np.random.default_rng(seed)
    def random_fields(n):
        return [{s:(rng.normal(size=d)+1j*rng.normal(size=d))/math.sqrt(2)
                 for s,d in dm.items()} for _ in range(n)]
    train=random_fields(samples)
    vm,mnames=matrices(pm,train);vr,rnames=matrices(pr,train)
    A,_,rank_r,singular=np.linalg.lstsq(vr,vm,rcond=1e-12)
    rank_m=int(np.linalg.matrix_rank(vm,tol=1e-10*max(1.,np.linalg.norm(vm,ord=2))))
    fit=float(np.linalg.norm(vr@A-vm)/max(np.linalg.norm(vm),1e-300))
    test=random_fields(max(32,samples//4))
    mtest,_=matrices(pm,test);rtest,_=matrices(pr,test)
    held=float(np.linalg.norm(rtest@A-mtest)/max(np.linalg.norm(mtest),1e-300))
    complete=(rank_r==len(rnames) and rank_m==len(mnames)
              and rank_r==rank_m and fit<tol and held<tol)
    return {'status':('Candidate_ExactPolynomialSubspace' if complete
                      else 'Blocked_OvercompleteOrDifferentPolynomialSpaces'),
            'matchete_symbols':mnames,'rgbeta_symbols':rnames,
            'matchete_native_record_slots':len(pm),'rgbeta_native_record_slots':len(pr),
            'matchete_reused_slots':len(pm)-len(mnames),
            'rgbeta_reused_slots':len(pr)-len(rnames),
            'map_rgbeta_from_matchete':A.tolist(),
            'rank_matchete':rank_m,'rank_rgbeta':int(rank_r),
            'singular_values_rgbeta':singular.tolist(),
            'fit_relative_residual':fit,'holdout_relative_residual':held,
            'independent_bases':rank_r==len(rnames) and rank_m==len(mnames),
            'native_polynomial_coordinates_assumed_equal':True,
            'all_model_quartics_covered_by_this_comparison':False,
            'physical_matchete_RGBeta_map_certified':False,
            'gamma_cc_spinor_signed_factor_certified':False}
