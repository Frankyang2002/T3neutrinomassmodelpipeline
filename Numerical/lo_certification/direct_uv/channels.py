"""Exact p=0 charge-channel bookkeeping using the unchanged direct UV engine.

Source-field charge projections are performed BEFORE fermion Takagi/scalar
real-Hessian rotations. This is essential: Dirac mass eigenstates need not
carry a unique charge in the doubled-Weyl Takagi basis.
"""
from __future__ import annotations

from collections import defaultdict
from typing import Any
import numpy as np

from Numerical.direct_one_loop.core import (
    _complex, encoded, b0_finite, scalar_real_basis, takagi, LOOP
)
from .certification import charge_bases


def _norm(m):
    return float(np.linalg.norm(m))


def partition_exact(data:dict[str,Any], kind:str, alpha:int, direct:dict[str,Any],
                    *, weak_fraction:float=0.01) -> dict[str,Any]:
    M=_complex(data['fermion_mass_gev'],'M')
    A=_complex(data['scalar_A_gev2'],'A')
    B=_complex(data['scalar_B_gev2'],'B')
    U,mf,_=takagi(M)
    C,ms2,_=scalar_real_basis(A,B)
    mu=float(data['renormalization_scale_gev'])
    qs=charge_bases(kind,alpha)['scalar']
    qset=sorted(set(qs.tolist()))
    nf=len(mf); nr=len(ms2)
    # Charge-tagged mass-eigenstate vertices. A tag identifies the charge of
    # the ORIGINAL scalar component, not that of a mixed real eigenstate.
    g_by_q={}
    for q in qset:
        mask=np.isclose(qs,q)[None,None,:]
        total=np.zeros((3,nf,nr),complex)
        for row in data['vertices']:
            y=_complex(row['field'],'field') if 'field' in row else np.zeros((3,nf,len(qs)),complex)
            z=_complex(row['conjugate'],'conjugate') if 'conjugate' in row else np.zeros((3,nf,len(qs)),complex)
            total+=np.einsum('ais,ik,sr->akr',y*mask,U,C,optimize=True)
            total+=np.einsum('ais,ik,sr->akr',z*mask,U,C.conj(),optimize=True)
        g_by_q[q]=total
    # A decomposition indexed by UNORDERED source-scalar charge pairs.
    # Cross pairs combine both orientations to retain complex symmetry.
    parts={}; by_heavy={}; by_scalar={}
    for qi,q in enumerate(qset):
        for qp in qset[qi:]:
            mat=np.zeros((3,3),complex)
            for k,heavy in enumerate(mf):
                for r,scalar2 in enumerate(ms2):
                    left=g_by_q[q][:,k,r]
                    right=g_by_q[qp][:,k,r]
                    piece=np.outer(left,right)
                    if q!=qp:piece=piece+piece.T
                    piece=LOOP*heavy*b0_finite(float(heavy**2),float(scalar2),mu**2)*piece
                    mat+=piece
                    label=f'{q:g}|{qp:g}'
                    by_heavy[(label,int(k))]=by_heavy.get((label,int(k)),np.zeros((3,3),complex))+piece
                    by_scalar[(label,int(r))]=by_scalar.get((label,int(r)),np.zeros((3,3),complex))+piece
            parts[f'{q:g}|{qp:g}']=mat
    total=_complex(direct['total']['complex_symmetric_matrix_gev'],'direct matrix',(3,3))
    reconstructed=sum(parts.values(),np.zeros((3,3),complex))
    # Independent charge-resolved accumulation and engine accumulation can
    # differ by roundoff after extremely strong cancellations.
    absolute_error=_norm(reconstructed-total)
    largest_part=max((_norm(x) for x in parts.values()),default=0.)
    roundoff_scale=max(_norm(total),1e-12*direct['diagnostics']['sum_individual_frobenius_gev'],1e-50)
    relative_error=absolute_error/roundoff_scale
    if absolute_error>2e-5*roundoff_scale:
        raise ValueError(f'Charge decomposition fails to reconstruct exact UV matrix: {relative_error:g}')
    nonzero=[(n,m) for n,m in parts.items() if _norm(m)>1e-50]
    by_heavy_all=defaultdict(lambda:np.zeros((3,3),complex))
    by_scalar_all=defaultdict(lambda:np.zeros((3,3),complex))
    for (_q,k),m in by_heavy.items():by_heavy_all[k]+=m
    for (_q,r),m in by_scalar.items():by_scalar_all[r]+=m
    source_pair={k:_complex(v,'source pair matrix',(3,3)) for k,v in
                 direct['vertex_pair_matrices_gev'].items()}
    max_pair=max((_norm(x) for x in source_pair.values()),default=0.)
    def record(matrix,scale):
        return {'matrix_gev':encoded(matrix),
                'frobenius_gev':_norm(matrix),
                'relative_to_full_mass':_norm(matrix)/max(_norm(total),1e-80),
                'relative_to_largest_in_partition':_norm(matrix)/max(scale,1e-80),
                'numerically_suppressed':bool(_norm(matrix)<weak_fraction*scale)}
    return {
        'charge_channel_rule':'unordered pair of ORIGINAL scalar component charges; includes interference',
        'charge_channel_matrices':{k:record(m,largest_part) for k,m in parts.items()},
        'active_charge_channels':[k for k,m in nonzero],
        'per_charge_and_fermion':{f'{q}:F{k}':record(m,largest_part) for (q,k),m in by_heavy.items()},
        'per_charge_and_scalar':{f'{q}:S{r}':record(m,largest_part) for (q,r),m in by_scalar.items()},
        'per_fermion_eigenstate':{str(k):record(m,max((_norm(x) for x in by_heavy_all.values()),default=0.))
                                   for k,m in sorted(by_heavy_all.items())},
        'per_scalar_eigenstate':{str(k):record(m,max((_norm(x) for x in by_scalar_all.values()),default=0.))
                                 for k,m in sorted(by_scalar_all.items())},
        'per_vertex_pair':{k:record(m,max_pair) for k,m in source_pair.items()},
        'reconstruction_absolute_error_gev':absolute_error,
        'reconstruction_relative_error':relative_error,
        'reconstruction_relative_error_denominator':'max(coherent norm, 1e-12 * sum individual term norms) to cover exact-zero cancellation',
        'suppression_definition':f'Frobenius norm below {weak_fraction:g} times partition maximum; never discarded',
        'cancellation_metrics': {
            'sum_charge_matrix_norms_gev':sum(map(_norm,parts.values())),
            'coherent_full_norm_gev':_norm(total),
            'charge_coherent_over_incoherent_norm_sum':_norm(total)/max(sum(map(_norm,parts.values())),1e-100),
            'individual_loop_coherent_over_sum':direct['diagnostics']['cancellation_ratio_total_over_sum'],
            'individual_loop_sum_over_coherent':direct['diagnostics']['cancellation_factor_sum_over_total'],
        },
        'warning':'Degenerate Takagi eigenvectors, individual terms and MSbar channel pieces depend on basis/scheme; only the full sum is invariant.',
    }
