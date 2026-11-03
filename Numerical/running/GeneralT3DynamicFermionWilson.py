"""Dynamic full-general F-tree LLSS tensors at the actual 3-family UV endpoint.

14 nonmixed T3 models: source-exact Matchete tree operator factorisation,
heavy-flavour SVD/Takagi, source CG expansion, lepton-index realignment, and
polynomial orbit symmetrisation. The full y1/y2/y1Conj/y2Conj inventory is
retained if it occurs in the current source operators.

B(-1),C(-1): optional explicit provisional reweighting of the pre-existing
mixed-scalar release, fitting its SU2/scalar component geometry against the
corresponding reference heavy-flavour bilinear. The mixed source does not yet
have the same native source/CG factorisation as the 14-model backend.

No GitHub writes, no interpolation of the yellow-region Wilson proxy.
"""
from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Mapping

import numpy as np

from Numerical.running.GeneralT3SequentialLeadingLog import require, decode

MIXED = frozenset(('B_-1','C_-1'))
F_ROOTS = ('output/general_observable/fermion_first_self_fix_v2',
           'output/general_observable/fermion_first_all16_v1')


def _dense_operator(operators: list[dict], rows: list[dict], model: str
                   ) -> tuple[np.ndarray, tuple[int, ...], dict]:
    """Matchete ordered terms -> fully symmetric formal coefficient tensor."""
    from Numerical.running.GeneralT3WilsonPolynomialBatch import (
        source_term_alignment, canonical_polynomial, distinct_orbit)
    from Numerical.running.GeneralT3OrderedTensor import scalar_basis

    basis = scalar_basis(model)
    n = basis.total_dimension-4
    require(n>0,'invalid heavy scalar dimension')
    by_index={int(row['source_index']):row for row in rows}
    tensors=[]; indices=[]; audits=[]
    for operator in operators:
        idx=int(operator['source_index'])
        if idx not in by_index:
            raise ValueError('Matchete source operator absent from aligned term list')
        align=source_term_alignment(operator,by_index[idx])
        perm=align['flavour_axis_permutation_for_lepton_order']
        ordered={}
        for component in operator['components']:
            i,j=component['lepton_su2']
            a,b=component['scalar_real']
            require(5<=a<=n+4 and 5<=b<=n+4,
                    f'{model}: unexpected Higgs leg at F tree LLSS threshold')
            val=decode(component['ordered_flavour_coefficient_gev_inv'])
            for p in range(3):
                for q in range(3):
                    x=complex(val[p,q])
                    if x==0:continue
                    flav=(p+1,q+1)
                    key=(i,flav[perm[0]],j,flav[perm[1]],a,b)
                    ordered[key]=ordered.get(key,0j)+x
        sym,_=canonical_polynomial(ordered)
        C=np.zeros((2,3,2,3,n,n),complex)
        for key,coeff in sym.items():
            for i,p,j,q,a,b in distinct_orbit(key):
                C[i-1,p-1,j-1,q-1,a-5,b-5]=complex(coeff)
        require(np.isfinite(C).all(),'nonfinite normalized F-tree source tensor')
        tensors.append(C)
        indices.append(idx)
        audits.append({'source_index':idx,'raw_ordered_components':len(ordered),
                       'symmetric_orbits':len(sym),
                       'flavour_permutation':perm})
    require(bool(tensors),'Matchete tree source has no nonzero operators')
    return np.stack(tensors),tuple(indices),{
        'sources':audits,'heavy_real_scalar_count':n,
        'operator_symmetrisation':'exact commuting-spinor polynomial orbit average'}


def nonmixed_boundary(model: str, uv_state: Any, *, repo_root: Path,
                      mu_f: float = 1e5, f_roots: tuple[str, ...] = F_ROOTS):
    """Evaluate tree Wilson boundary using ACTUAL fitted muF state, all 14."""
    from scripts.RunT3GeneralWilsonBoundary import (
        _first_threshold, _validate_seed, flavour_matrix, encode_matrix)
    from Numerical.matching.GeneralT3FlavourC5 import prepare_heavy_flavour
    from Numerical.running.GeneralT3OrderedTensor import project_ordered_term
    from Numerical.joint_scalar_matching.fermion_boundary_bridge import SourceBoundary

    require(model not in MIXED,'mixed scalar matching has a separate source')
    roots=tuple(Path(repo_root)/x for x in f_roots)
    _status,seed_path,_certificate = _first_threshold(model,roots)
    doc=json.loads(seed_path.read_text(encoding='utf-8-sig'))
    terms=_validate_seed(doc,model)
    required={f['name'] for row in terms for f in row['factors']}
    masses,yukawas,heavy_diag=prepare_heavy_flavour(uv_state,required)
    operators=[]
    for row in terms:
        axes=[f['external_index'] for f in row['factors']]
        row['flavour_matrix_axes_external_indices']=axes
        row['ordered_flavour_matrix_gev_inv']=encode_matrix(
            flavour_matrix(row,masses,yukawas))
        operators.append(project_ordered_term(row,doc['CGRegistry'],model))
    by_source, indices, audit=_dense_operator(operators,terms,model)
    raw=np.sum(by_source,axis=0)
    require(abs(float(uv_state.mu_gev)-mu_f)<1e-6,
            'actual UV state is not at fermion matching threshold')
    signature=sha256((seed_path.read_bytes()+raw.tobytes())).hexdigest()
    n=raw.shape[-1]
    from Numerical.joint_scalar_matching.special_mass_v6 import (
        SPECIAL, expected_special_labels)
    if model in SPECIAL:
        labels=expected_special_labels(model)
    else:
        d1={'A':1,'B':2,'C':2,'D':3,'E':3}[model[0]]
        labels=tuple(['NewScalar1']*(2*d1)+['NewScalar2']*(n-2*d1))
    source=SourceBoundary(model,seed_path,signature,float(mu_f),raw,by_source,
                          indices,labels,False,'DynamicGeneralFermionTreeCandidate')
    return source,{
        'matched_source':str(seed_path),
        'f_tree_operator_count':len(terms),
        'Yukawa_types':sorted(required),
        'heavy_mass_eigenvalues_gev':np.asarray(masses,float).tolist(),
        'heavy_basis_diagnostics':str(heavy_diag),
        **audit,
        'absolute_matchete_to_weyl_normalization_certified':False}


def _matrix(raw: Any) -> np.ndarray:
    if isinstance(raw,Mapping) and {'real','imag'}<=set(raw):
        return np.asarray(raw['real'],float)+1j*np.asarray(raw['imag'],float)
    if isinstance(raw,list) and len(raw)==3:
        return np.array([[complex(float(v['real']),float(v['imag']))
                          if isinstance(v,Mapping) else complex(v) for v in row]
                         for row in raw],complex)
    return np.asarray(raw,complex)


def _heavy_bilinear(state: Any) -> np.ndarray:
    from Numerical.matching.GeneralT3FlavourC5 import prepare_heavy_flavour
    masses,yukawas,_ = prepare_heavy_flavour(state,{'y1','y2'})
    a,b=yukawas['y1'],yukawas['y2']
    W=(a/masses[None,:])@b.T
    # The Majorana source contains both lepton orientations in its formal
    # commuting-spinor polynomial, so use the symmetrized external bilinear.
    return (W+W.T)/2


def mixed_reference_boundary(model: str, ref_source: Any,
                             reference_uv: Any, current_uv: Any):
    """Explicit B/C candidate flavour lift of CURRENT muF heavy-state inputs.

    Each source's original SU2/scalar component is a least-squares multiple
    of W_ref. Replacing it by W_new retains that component's reference CG.
    This approximation is NOT the yet-missing native mixed-tree factorisation.
    """
    from Numerical.joint_scalar_matching.fermion_boundary_bridge import SourceBoundary
    require(model in MIXED and ref_source.model==model,'wrong mixed source')
    A=_heavy_bilinear(reference_uv)
    B=_heavy_bilinear(current_uv)
    den=float(np.vdot(A,A).real)
    require(den>1e-50 and np.isfinite(den),'reference heavy Yukawa bilinear vanishes')
    tensors=[]; max_res=0.
    for tensor in ref_source.by_source:
        # tensor[i,p,j,q,a,b] -> component[i,j,a,b,p,q]
        axes=np.transpose(tensor,(0,2,4,5,1,3))
        geom=np.einsum('ijabpq,pq->ijab',axes,A.conj(),optimize=True)/den
        estimated=np.einsum('ijab,pq->ijabpq',geom,A,optimize=True)
        den_rel=max(float(np.linalg.norm(axes)),1e-40)
        max_res=max(max_res,float(np.linalg.norm(axes-estimated)/den_rel))
        updated=np.einsum('ijab,pq->ipjqab',geom,B,optimize=True)
        # Preserve exact polynomial symmetry at the boundary.
        updated=(updated+updated.transpose(2,3,0,1,4,5)
                 +updated.transpose(0,1,2,3,5,4)
                 +updated.transpose(2,3,0,1,5,4))/4
        tensors.append(updated)
    by=np.stack(tensors)
    raw=by.sum(axis=0)
    src=SourceBoundary(model,ref_source.path,
                       sha256(ref_source.hash.encode()+raw.tobytes()).hexdigest(),
                       ref_source.mu_f_gev,raw,by,ref_source.source_indices,
                       ref_source.labels,False,
                       'MixedReferenceSourceCG_FittedFlavourCandidate')
    return src,{'status':'provisional_mixed_mass_source_flavour_lift',
                'reference_bilinear_reconstruction_max_relative_residual':max_res,
                'reference_source':str(ref_source.path),
                'current_Yukawa_bilinear_norm':float(np.linalg.norm(B)),
                'spinor_factor_certified':False}


def uv_from_json(model: str, values: Mapping[str, Any], *, mu_f: float=1e5):
    """Minimal UV heavy-basis input from GeneralRunningBenchmark state dump."""
    from types import SimpleNamespace
    extra={x:_matrix(values[x]) for x in ('y1Conj','y2Conj') if x in values}
    required=('MF','y1','y2')
    missing=[x for x in required if x not in values]
    require(not missing,f'{model}: missing UV heavy inputs {missing}')
    return SimpleNamespace(mu_gev=mu_f,MF=_matrix(values['MF']),
             y1=_matrix(values['y1']),y2=_matrix(values['y2']),
             general_couplings=extra,
             representation=SimpleNamespace(majorana_fermion=model in MIXED))
