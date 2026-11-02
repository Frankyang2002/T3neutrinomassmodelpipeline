"""Exact extension of the native HHSS portal tensor basis by RGBeta directions.

Both fields S1 and S2 remain active until the *same* mu_S.  The coefficients
below are algebraically exact in the stated, common highest-weight complex
component frame.  They do not certify that special Matchete boundary fields
or the RGBeta runtime share that frame.

Existing native quartic channels are NEVER dropped or overwritten.  The chosen
complement is a linearly independent RGBeta singlet/adjoint column; it is not
fitted to data, and does not assert a new Matchete coupling was generated.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any
import numpy as np
import sympy as sp

from .portal_normalization import (
    _matrix_rows, native_portal_columns, rgbeta_portal_columns,
    tensor_to_dense, real_running_coupling, TensorAuditError,
)

DIMS = {'A': (1,3), 'B': (2,2), 'C': (2,2), 'D': (3,1), 'E': (3,3)}
FORMAT = 't3-ws2-portal-completed-native-plus-rgbeta-v3'


def exact_portal_completion(ch: Any, species: int) -> dict:
    """Return full-rank transformation: RGBeta tensor = augmented native basis.

    If N is the native tensor-column matrix and R the RGBeta one, construct
    B=[N,R_selected], where each selected R adds an independent direction.
    Then B X=R exactly. Consequently x=X c_RGB is the unique augmented
    coupling vector, and the full real HHSS derivative is unchanged.
    """
    if species not in (1,2):
        raise TensorAuditError('species must be 1 or 2')
    d=DIMS[ch.model[0]][species-1]
    rgb0=rgbeta_portal_columns(d)
    rgb={f'lambdaH{species}'+('Adj' if key=='adjoint' else ''):value
         for key,value in rgb0.items()}
    native=native_portal_columns(ch,species)
    if not native:
        raise TensorAuditError(f'{ch.model}: no native S{species} portal channel')
    keys=sorted(set().union(*(set(v) for v in (*rgb.values(),*native.values()))))
    R=_matrix_rows(rgb,keys)
    N=_matrix_rows(native,keys)
    if R.rank()!=len(rgb):
        raise TensorAuditError('RGBeta portal basis is not independent')
    # This guards against falsely 'completing' an incompatible component frame.
    if (N.row_join(R)).rank()!=R.rank():
        raise TensorAuditError('Native tensor outside RGBeta span; independent phase/basis map required')
    if N.rank()!=N.cols:
        raise TensorAuditError('Native portal channels linearly dependent')
    chosen=[]
    B=N.copy()
    for name in rgb:
        trial=B.row_join(R[:,list(rgb).index(name)])
        if trial.rank()>B.rank():
            chosen.append(name)
            B=trial
    if B.rank()!=R.rank() or B.cols!=R.cols:
        raise TensorAuditError('Failed to complete portal basis without duplicate directions')
    pivots=B.T.rref()[1]
    X=B.extract(pivots,range(B.cols)).inv()*R.extract(pivots,range(R.cols))
    residual=B*X-R
    if any(sp.simplify(x)!=0 for x in residual):
        raise TensorAuditError('Full HHSS tensor reconstruction failed')
    native_names=list(native)
    added_names=[f'WS2_RGBeta_complement_S{species}_{name}' for name in chosen]
    augmented={**native,**{new:rgb[old] for new,old in zip(added_names,chosen)}}
    matrix=[[sp.sstr(sp.factor(X[i,j])) for j in range(X.cols)] for i in range(X.rows)]
    return {
        'format':FORMAT,'model':ch.model,'species':f'S{species}',
        'status':'exact_full_portal_tensor_basis_in_declared_original_component_frame',
        'native_symbols':native_names,'added_complement_symbols':added_names,
        'added_complement_rgbeta_sources':chosen,
        'augmented_symbols':list(augmented),'rgbeta_symbols':list(rgb),
        'augmented_from_rgbeta_exact_matrix':matrix,
        'native_direction_count':len(native),'rgbeta_direction_count':len(rgb),
        'completion_count':len(chosen),
        'basis_identity':'V4=-L4, 4th derivative, standard highest-weight SU2 complex components; no Matchete special-basis certification',
        'source_matchete_native_basis_provisional_for_special_models':True,
        'full_EFT1_potential_certified':False,
    }


def completed_portal_numeric(ch:Any,states:dict[str,dict]) -> tuple[dict,dict]:
    """Evaluate both scales with no silent zeros, preserving every native channel.

    ``states`` is {'muF': EFT1 state(muF), 'muS': EFT1 state(muS)}.
    Every required RGBeta singlet/adjoint coupling must be numeric at its own
    scale.  Return exact-basis and array data in original S1/S2 real indices.
    """
    d1,d2=DIMS[ch.model[0]]
    n=2*(d1+d2)
    audits={species:exact_portal_completion(ch,species) for species in (1,2)}
    arr={};scale_details={};provenance={}
    for scale,state in states.items():
        channels=[];labels=[];check=np.zeros((4,4,n,n),float)
        calculated=np.zeros_like(check)
        values={};missing=[]
        for species,dim,offset in ((1,d1,0),(2,d2,2*d1)):
            cert=audits[species]
            req=cert['rgbeta_symbols']
            absent=sorted(set(req)-set(state))
            if absent:
                missing.extend(absent)
                continue
            c=sp.Matrix([sp.Rational(0) + sp.Float(real_running_coupling(state[s],s),17)
                         for s in req])
            X=sp.Matrix([[sp.sympify(x) for x in row]
                         for row in cert['augmented_from_rgbeta_exact_matrix']])
            numeric=X*c
            native=native_portal_columns(ch,species)
            rgb0=rgbeta_portal_columns(dim)
            added={name:rgb0['adjoint' if symbol.endswith('Adj') else 'singlet']
                   for name,symbol in zip(cert['added_complement_symbols'],
                                          cert['added_complement_rgbeta_sources'])}
            unit={**native,**added}
            for label,coefficient in zip(cert['augmented_symbols'],numeric):
                number=complex(sp.N(coefficient,16))
                if abs(number.imag)>1e-12*max(1.,abs(number.real)):
                    raise TensorAuditError(f'{label}: non-real coefficient')
                contribution=np.zeros_like(check)
                contribution[:,:,offset:offset+2*dim,offset:offset+2*dim]=(
                    number.real*tensor_to_dense(unit[label],dim))
                labels.append(label);channels.append(contribution)
                calculated+=contribution
                values[label]=number.real
            for kind,col in rgb0.items():
                symbol=f'lambdaH{species}'+('Adj' if kind=='adjoint' else '')
                check[:,:,offset:offset+2*dim,offset:offset+2*dim]+=(
                    real_running_coupling(state[symbol],symbol)*tensor_to_dense(col,dim))
        if not missing:
            residual=float(np.max(np.abs(calculated-check)))
            if not np.allclose(calculated,check,rtol=2e-12,atol=1e-12):
                raise TensorAuditError(f'Augmented portal tensors fail {scale} reconstruction: {residual}')
            status='exact_portal_numeric_reconstruction_conditional_component_frame'
        else:
            # No physical partial tensor is exported when one portal is missing.
            status='blocked_missing_same_scale_portal_coupling'
            channels=[];labels=[];calculated=np.zeros_like(check)
            values={}
            residual=None
        if not missing:
            arr[f'{scale}_portal_sum']=calculated
            arr[f'{scale}_channel_tensors']=np.stack(channels,axis=0)
            arr[f'{scale}_channel_labels']=np.asarray(labels,dtype='U128')
        scale_details[scale]={
            'status':status,'missing_rgbeta_couplings':missing,
            'native_and_completion_couplings':values,
            'tensor_relative_reconstruction_residual':residual,
            'channels_included':labels,'input_renormalization_scale':scale}
    return {
       'format':FORMAT,'model':ch.model,
       'status':'baseline_portal_tensor_basis_completed',
       'species':[audits[1],audits[2]],
       'baseline_missing_native_directions_filled':sum(a['completion_count'] for a in audits.values()),
       'scales':scale_details,
       'renormalization_scheme':'MSbar',
       'source_scalar_basis':'original S1/S2 canonical real components, interleaved Re/Im',
       'quartic_normalization':'V4=lambda_abcd phi_a phi_b phi_c phi_d/4!',
       'potential_sign':'V4=-L4 quartic interaction',
       'higher_HHSS_and_lambdaT3_included':False,
       'physical_complete_EFT1_input':False,
       'suitable_for_existing_joint_scalar_matching_as_full_input':False,
       'suitable_for_existing_seven_term_RGE_as_full_input':False,
    },arr
