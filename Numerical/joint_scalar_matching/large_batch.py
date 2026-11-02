"""All-16 source kernel parity + conditional fixed-order EFT1-to-SMEFT response.

There is NO implicit GammaCC/Weyl normalization, no invented running of
heavy scalar couplings, no independent finite matching certification.  The
conditional C5 is per *unit unknown signed normalization* N_GammaCC.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any
import json
import math
import re

import numpy as np

from .adapter import CONVENTION
from .core import LOOP, EigenBasis, check_tensors, project_weinberg, scalar_bubble_matching
from .fermion_boundary_bridge import SourceBoundary, digest, read_json, require

KERNEL_FORMAT = 't3-eft1-conditional-quartic-pivot-kernel-v1'
BACKGROUND_FORMAT = 't3-joint-scalar-conditional-background-v1'
BATCH_FORMAT = 't3-joint-scalar-integrated-conditional-batch-v2'
KERNEL_STATUS = 'ProvisionalLinearQuarticKernel_NotPhysicalC5'
EXPECTED_ORDER = ('H0R','H0R')


def encode_matrix(a: np.ndarray) -> dict:
    a = np.asarray(a, dtype=np.complex128)
    require(a.shape == (3,3) and np.isfinite(a).all(), 'nonfinite or non-3x3 matrix')
    return {'real':a.real.tolist(),'imag':a.imag.tolist()}


def parse_matrix(obj: Any) -> np.ndarray:
    require(isinstance(obj,dict) and set(('real','imag')).issubset(obj),
            'kernel complex matrix must have real/imag fields')
    a=np.asarray(obj['real'],float)+1j*np.asarray(obj['imag'],float)
    require(a.shape==(3,3) and np.isfinite(a).all(), 'nonfinite or wrong kernel matrix shape')
    return a


def _expected_pivot(source: SourceBoundary, idx: int, pair: tuple[int,int]) -> np.ndarray:
    a,b=pair
    n=len(source.labels)
    require(5<=a<=b<5+n,'pivot kernel has out-of-range, noncanonical scalar pair')
    # Existing source kernel's reported normalization: 4*(2-delta_ab)*
    # C_formal[L_su2=1,L_su2=1,a,b].  This is a SOURCE POLYNOMIAL
    # convention, and MUST NOT be silently equated with the full-tensor
    # Weinberg-projected master RGE until the spinor bridge is certified.
    mult=1 if a==b else 2
    return 4*mult*source.by_source[idx,:,:,0,0,a-5,b-5]


def validate_pivot_kernel(path: Path, source: SourceBoundary) -> dict:
    obj=read_json(path)
    require(obj.get('format')==KERNEL_FORMAT and obj.get('status')==KERNEL_STATUS,
            'unexpected source-pivot kernel version/status')
    require(obj.get('model')==source.model and obj.get('source_boundary_sha256')==source.hash,
            'source-pivot kernel is not tied to this fermion release')
    require(obj.get('real_scalar_dimension')==4+len(source.labels) and
            obj.get('H0R_real_scalar_index_one_based')==3,
            'pivot kernel field layout changed')
    norm=obj.get('physical_normalization',{})
    require(norm.get('certified') is False and norm.get('N_GammaCC_to_Weyl') is None,
            'source-pivot kernel must not silently assert a physical spinor normalization')
    require(obj.get('quartic_component_couplings') ==
            'UNSUPPLIED - no numerical C5 is computed',
            'kernel should not claim numerical quartic matching')
    rows=obj.get('per_source')
    require(isinstance(rows,list) and len(rows)==len(source.source_indices),
            'source pivot inventory differs from LLSS boundary')
    summed: dict[tuple[int,int],np.ndarray] = {}
    pair_contribution_scale: dict[tuple[int,int],float] = {}
    max_abs_diff=0.0
    n_terms=0
    for t,row in enumerate(rows):
        require(row.get('matchete_source_index')==source.source_indices[t],
                'Matchete source order changed in pivot kernel')
        entries=row.get('kernel_channels')
        require(isinstance(entries,list),'missing pivot channel list')
        expected_pairs=[]
        n=len(source.labels)
        for a in range(n):
            for b in range(a,n):
                if np.max(np.abs(source.by_source[t,:,:,0,0,a,b]))>0:
                    expected_pairs.append((5+a,5+b))
        seen=set()
        for entry in entries:
            pair=entry.get('heavy_scalar_pair_one_based')
            require(isinstance(pair,list) and len(pair)==2 and
                    all(isinstance(x,int) and not isinstance(x,bool) for x in pair),
                    'invalid pivot pair index')
            key=tuple(pair)
            require(key not in seen, 'duplicate scalar pair in source pivot')
            seen.add(key)
            require(entry.get('multiplicity_in_ordered_ab_sum')==(1 if key[0]==key[1] else 2),
                    'wrong source pivot pair multiplicity')
            require(entry.get('indeterminate_quartic')==f'lambda_[H0R,H0R,{key[0]},{key[1]}]',
                    'source pivot quartic label changed')
            actual=parse_matrix(entry['conditional_beta_C5_per_unit_quartic'])
            expected=_expected_pivot(source,t,key)
            diff=float(np.max(np.abs(actual-expected)))
            max_abs_diff=max(max_abs_diff,diff)
            scale=max(float(np.max(abs(actual))),float(np.max(abs(expected))),1e-30)
            require(diff<=2e-9*scale+1e-28,
                    f'pivot kernel disagrees with independent source coefficient {key}: {diff:g}')
            summed[key]=summed.get(key,np.zeros((3,3),complex))+expected
            pair_contribution_scale[key]=pair_contribution_scale.get(key,0.)+float(np.max(abs(expected)))
            n_terms+=1
        require(set(expected_pairs)==seen,'pivot omitted or added a nonzero scalar pair')
    combined=obj.get('summed_kernel_channels')
    require(isinstance(combined,list),'summed pivot inventory missing')
    # A whole scalar pair can cancel between distinct PL sources.  The upstream
    # pivot worker deliberately omits fully cancelled pairs (e.g. E_-2).
    # Do not demand a summed row for such a pair, but demand that its net
    # coefficient is numerically zero relative to its unsummed sources.
    seen=set()
    for entry in combined:
        key=tuple(entry.get('heavy_scalar_pair_one_based',()))
        require(key in summed and key not in seen,'summed pivot has duplicate/unknown pair')
        seen.add(key)
        actual=parse_matrix(entry['conditional_beta_C5_per_unit_quartic'])
        diff=float(np.max(np.abs(actual-summed[key])))
        max_abs_diff=max(max_abs_diff,diff)
        scale=max(float(np.max(abs(actual))),float(np.max(abs(summed[key]))),1e-30)
        require(diff<=2e-9*scale+1e-28,'summed source pivot coefficients mismatch')
    for key,value in summed.items():
        if key not in seen:
            require(float(np.max(abs(value))) <= 2e-9*pair_contribution_scale[key]+1e-28,
                    f'summed pivot omitted an uncancelled source pair {key}')
    return {'pivot_kernel_sha256':digest(path),'source_count':len(rows),
            'pivot_per_source_channels':n_terms,'distinct_scalar_pairs':len(combined),
            'cancelled_source_pairs_omitted_by_kernel':len(summed)-len(combined),
            'source_vs_kernel_max_abs_error':max_abs_diff,
            'pivot_convention':'4*(2-delta_ab)*C_formal[lepton_SU2=1,lepton_SU2=1,a,b]',
            'physical_normalization_certified':False}


@dataclass(frozen=True)
class ScalarBackground:
    path: Path
    model: str
    mu_f: float
    mu_s: float
    reference_scale: float
    mass_origin: str
    U: np.ndarray
    mass2: np.ndarray
    quartic: np.ndarray
    channels: np.ndarray
    names: tuple[str,...]
    provenance: dict


def load_background(path: Path, source: SourceBoundary) -> ScalarBackground:
    path=Path(path)
    obj=read_json(path)
    require(obj.get('format')==BACKGROUND_FORMAT,'wrong conditional background schema')
    require(obj.get('model')==source.model and obj.get('source_boundary_sha256')==source.hash,
            'conditional background not sourced from correct fermion LLSS release')
    require(obj.get('source_scalar_real_labels')==list(source.labels),
            'background scalar source real basis ordering changed')
    require(obj.get('renormalization_scheme')=='MSbar' and
            obj.get('tensor_convention')==CONVENTION,'background scheme/tensor convention mismatch')
    require(obj.get('no_separate_scalar_threshold') is True and
            obj.get('same_muS_quartic_and_mass') is True,
            'both scalar multiplets must share one scalar threshold and scale')
    require(obj.get('scalar_mass_origin') in ('unbroken_electroweak','EWSB_resummed'),
            'scalar mass origin must be declared; no synthetic physical claims')
    provenance=obj.get('provenance')
    require(isinstance(provenance,dict) and
            isinstance(provenance.get('source'),str) and provenance['source'].strip() and
            provenance.get('quartic_potential_normalization')=='V=lambda_abcd*phi_a*phi_b*phi_c*phi_d/4!',
            'independently supplied provenance and potential normalization required')
    mf=float(obj.get('mu_f_gev',float('nan')))
    ms=float(obj.get('mu_s_gev',float('nan')))
    ref=float(obj.get('reference_scale_gev',ms))
    require(all(map(math.isfinite,(mf,ms,ref))) and 0<ms<=mf and ref>0 and
            abs(mf/source.mu_f_gev-1)<1e-11,'background matching scales mismatch')
    name=obj.get('tensor_npz')
    require(isinstance(name,str) and name and Path(name).name==name and name not in ('.','..'),
            'adjacent NPZ basename required')
    with np.load(path.parent/name,allow_pickle=False) as arr:
        required={'source_to_matching_real_U','mass_squared_gev2_muS',
                  'lambda_HHSS_muS','lambda_HHSS_channels_muS'}
        require(required.issubset(arr.files),f'NPZ missing {sorted(required-set(arr.files))}')
        U=np.array(arr['source_to_matching_real_U'])
        mass2=np.array(arr['mass_squared_gev2_muS'])
        q=np.array(arr['lambda_HHSS_muS'])
        channels=np.array(arr['lambda_HHSS_channels_muS'])
    n=len(source.labels)
    require(U.shape==(n,n) and np.isrealobj(U) and np.isfinite(U).all() and
            np.allclose(U.T@U,np.eye(n),atol=1e-10,rtol=1e-10),
            'invalid source-to-matching orthogonal transform')
    EigenBasis.from_mass2(mass2)
    require(mass2.shape==(n,n),'mass matrix dimension mismatch')
    check_tensors(source.raw,q,n)
    names=tuple(obj.get('quartic_channel_names',()))
    require(channels.ndim==5 and channels.shape[1:]==q.shape and len(channels)==len(names)
            and len(names)>0 and len(set(names))==len(names)
            and all(isinstance(x,str) and x.strip() for x in names),
            'quartic channels missing or mismatched')
    for ch in channels: check_tensors(source.raw,ch,n)
    err=float(np.max(abs(channels.sum(axis=0)-q)))
    require(err<=1e-10*max(float(np.max(abs(q))),1e-30)+1e-28,
            'quartic channels do not sum to the full quartic tensor')
    return ScalarBackground(path,source.model,mf,ms,ref,obj['scalar_mass_origin'],
                            U,mass2,q,channels,names,provenance)


def _rotate(C: np.ndarray,U:np.ndarray)->np.ndarray:
    return np.einsum('fgijab,ak,bl->fgijkl',C,U,U,optimize=True)


def _beta_project(C:np.ndarray,q:np.ndarray,*,strict=False) -> tuple[np.ndarray,float]:
    # For the particular projected convention implemented by the joint bubble
    # kernel, 16pi² dC_LLhh/dln(mu) |quartic = lambda_HHSS*C_LLSS.
    amplitude=np.einsum('fgijab,ABab->fgijAB',C,q,optimize=True)
    return project_weinberg(amplitude,strict=strict)


def _response(C:np.ndarray,q:np.ndarray,mass2:np.ndarray,mf:float,ms:float):
    """Fixed-order LLSS-induced O(hbar) log + single joint-scalar bubble.

    No C_LLSS self-running: its hbar correction multiplied by the scalar
    hbar bubble starts at hbar², outside the claimed order.
    """
    beta,beta_residual=_beta_project(C,q,strict=False)
    generated=math.log(ms/mf)*beta/LOOP
    chunk=scalar_bubble_matching(C,q,mass2,ms,reference_scale=mf)
    bubble,bubble_residual=project_weinberg(chunk['amplitude'],strict=False)
    return {'beta':beta,'running':generated,'bubble':bubble,
            'sum':generated+bubble,'beta_residual':beta_residual,
            'bubble_residual':bubble_residual}


def evaluate_conditional(source:SourceBoundary, background:ScalarBackground)->dict:
    require(source.model==background.model,'model mismatch')
    C=_rotate(source.raw,background.U)
    q=background.quartic
    mass=background.mass2
    mf,ms=background.mu_f,background.mu_s
    whole=_response(C,q,mass,mf,ms)
    # Gauge Ward identity is checked on the FULL amplitude, not individual
    # quartics or PL operators, which need not be Weinberg by themselves.
    require(whole['beta_residual']<=1e-8 and whole['bubble_residual']<=1e-8,
            'full LLhh amplitude violates the Weinberg Ward projection; check scalar basis/quartics')
    # Compare a second reference matching scale while holding source tensors,
    # masses and couplings fixed.  The total O(hbar) result must be muS-invariant.
    other_scale=ms*math.exp(0.03125)
    shift=_response(C,q,mass,mf,other_scale)
    d=shift['sum']-whole['sum']
    size=max(float(np.linalg.norm(whole['running'])),
             float(np.linalg.norm(whole['bubble'])),
             float(np.linalg.norm(whole['sum'])),1e-30)
    r=float(np.linalg.norm(d)/size)
    require(r<=2e-9,'matching-scale log cancellation failed')
    sources=[]
    for t,index in enumerate(source.source_indices):
        resp=_response(_rotate(source.by_source[t],background.U),q,mass,mf,ms)
        sources.append({'matchete_PL_source_index':index,
                        'C5_leading_log_per_N_gev_inverse':encode_matrix(resp['running']),
                        'C5_finite_joint_scalar_per_N_gev_inverse':encode_matrix(resp['bubble']),
                        'C5_total_LLSS_per_N_gev_inverse':encode_matrix(resp['sum']),
                        'individual_weinberg_projection_residual':max(resp['beta_residual'],resp['bubble_residual'])})
    qs=[]
    for name,ch in zip(background.names,background.channels):
        resp=_response(C,ch,mass,mf,ms)
        qs.append({'name':name,
                   'C5_leading_log_per_N_gev_inverse':encode_matrix(resp['running']),
                   'C5_finite_joint_scalar_per_N_gev_inverse':encode_matrix(resp['bubble']),
                   'C5_total_LLSS_per_N_gev_inverse':encode_matrix(resp['sum']),
                   'individual_weinberg_projection_residual':max(resp['beta_residual'],resp['bubble_residual'])})
    def unpack(rows,key):
        return sum((parse_matrix(x[key]) for x in rows),np.zeros((3,3),complex))
    for rows in (sources,qs):
        for key,expected in [('C5_leading_log_per_N_gev_inverse',whole['running']),
                             ('C5_finite_joint_scalar_per_N_gev_inverse',whole['bubble']),
                             ('C5_total_LLSS_per_N_gev_inverse',whole['sum'])]:
            scale=max(float(np.max(abs(expected))),1e-30)
            require(np.max(abs(unpack(rows,key)-expected))<=1e-8*scale+1e-27,
                    'source/quartic decomposition fails to reproduce the total')
    eigen=EigenBasis.from_mass2(mass)
    return {'format':BATCH_FORMAT,'model':source.model,
            'status':'ConditionalOneLoopLLSSResponse_PerUnitUncertifiedSpinorNormalization',
            'model_source_boundary_sha256':source.hash,
            'background_json_sha256':digest(background.path),
            'background_npz_sha256':digest(background.path.parent/
                    read_json(background.path)['tensor_npz']),
            'renormalization_scheme':'MSbar',
            'mu_f_gev':mf,'mu_s_gev':ms,
            'scalar_mass_origin':background.mass_origin,
            'source_to_matching_real_U':background.U.tolist(),
            'scalar_eigenstate_masses_gev':np.sqrt(eigen.masses2).tolist(),
            'C5_LLSS_beta_16pi2_per_N_gev_inverse':encode_matrix(whole['beta']),
            'C5_EFT1_induced_leading_log_per_N_gev_inverse':encode_matrix(whole['running']),
            'C5_joint_scalar_finite_bubble_per_N_gev_inverse':encode_matrix(whole['bubble']),
            'C5_LLSS_total_one_loop_per_N_gev_inverse':encode_matrix(whole['sum']),
            'per_PL_source':sources,'per_independent_quartic':qs,
            'weinberg_beta_projection_residual':whole['beta_residual'],
            'weinberg_bubble_projection_residual':whole['bubble_residual'],
            'relative_matching_scale_cancellation_residual':r,
            'source_to_eft_weyl_normalization_factor':None,
            'physical_C5_certified':False,'full_one_loop_certified':False,
            'missing_coverage':['signed Matchete GammaCC to Weyl operator normalization',
                 'independent same-scale scalar mass/quartic provenance certification',
                 'hard one-loop matching of heavy fermions',
                 'possible additional EFT1 operator insertions and finite counterterms'],
            'interpretation':('The full LLSS-induced O(hbar) response is shown per unit '
                 'unknown signed source-to-Weyl factor N. The LLSS-induced EFT1 log '
                 'and one joint scalar bubble are added ONCE. No hbar LLSS self-running '
                 'is included (it would produce hbar² in this sector). '
                 'Not a physical neutrino mass or complete full-theory one-loop calculation.'),
            'sequence':'UV -> integrate out all F at muF -> EFT1(S1,S2) -> integrate out all scalars together at muS -> SMEFT'}
