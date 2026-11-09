"""Workstream A v8: basis-independent SU(2) triplet portal spectrum audit.

The representation-theoretic eigenvalues are independent of any numerical
RGBeta S2 frame.  Native RGBeta symbolic reduction is tested separately by
scripts/RunT3ParallelAProbeNativeS2ReductionV8.wl.

This module NEVER extracts RGBeta numerical axes, certifies physical coupling
normalization, or authorizes Wilson evolution.
"""
from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
from typing import Any

import numpy as np

from Numerical.parallel_t3.wilson.s2_intertwiner_v7 import (
    S2EvidenceError, read_json, rel, audit_from_paths as v7_revalidate,
    STATUS as V7_STATUS,
)

FORMAT = 't3-parallel-a-native-s2-reduction-and-basis-invariant-spectrum-v8'
STATUS_REDUCED = ('S2BasisInvariantSpectrum14of14_NativeRGBetaFormalReductionObserved_'
                  'NumericalS2FrameStillUnspecifiedConditional')
STATUS_UNEXPOSED = ('S2BasisInvariantSpectrum14of14_NativeRGBetaFormalReductionUnobserved_'
                    'NumericalS2FrameStillUnspecifiedConditional')
NATIVE_FORMAT = 't3-parallel-a-native-rgbeta-s2-reduction-probe-v8'
NATIVE_REDUCED = 'NativeRGBetaS2SymbolicReductionObserved_NumericalFrameUnspecified'
NATIVE_UNEXPOSED = 'NativeRGBetaS2SymbolicReductionNotExposed_NumericalFrameUnspecified'
SOURCE_FILES = ('RGE/running/rgbeta/T3RGBetaModel.wl',
                'scripts/RunT3ParallelAProbeNativeS2ReductionV8.wl')
TOL = 2.e-10
TRIPLET_COUNT = 14
EXPECTED_EIGENVALUES = np.array([-1.,-1.,.5,.5,.5,.5],dtype=float)


def _identity(d: int) -> np.ndarray:
    return np.einsum('ij,ab->ijab',np.eye(2),np.eye(d)).astype(complex)


def require(yes: bool, message: str) -> None:
    if not yes:
        raise S2EvidenceError('Workstream A v8: '+message)


def _sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _spectral(k: np.ndarray) -> dict:
    require(k.shape == (2,2,3,3), 'S2 portal tensor must be (2,2,3,3)')
    # (bar H_i, H_j, bar S_u, S_v) -> H_i S_u | H_j S_v
    km=k.transpose(0,2,1,3).reshape(6,6)
    herm=rel(km,km.conj().T)
    hermitian=(km+km.conj().T)/2
    eig=np.linalg.eigvalsh(hermitian)
    identity=np.eye(6,dtype=complex)
    poly=rel(km@km+.5*km-.5*identity,np.zeros((6,6)))
    out={
        'eigenvalues':eig.tolist(),
        'hermiticity_residual':herm,
        'spectrum_residual':rel(eig,EXPECTED_EIGENVALUES),
        'minimal_polynomial_residual':poly,
        'trace_residual':float(abs(np.trace(km))),
        'trace_square_residual':float(abs(np.trace(km@km)-3.)),
        'operator_norm':float(np.linalg.norm(km)),
    }
    out['max_residual']=max(out[x] for x in (
        'hermiticity_residual','spectrum_residual','minimal_polynomial_residual',
        'trace_residual','trace_square_residual'))
    return out


def select_spectral_sign(raw: list[np.ndarray]) -> tuple[np.ndarray, dict]:
    """Determine adjoint sign using only Matchete native CG and SU2 Casimirs.

    Crucially, v6's selected adjoint coefficient is NOT used to select sign.
    """
    require(len(raw)==2,'Triplet must have two independent Matchete portal CGs')
    ident=_identity(3)
    v=[]
    for t in raw:
        require(t.shape==(2,2,3,3),'Unexpected Matchete portal CG dimensions')
        coeff=np.vdot(ident,t)/np.vdot(ident,ident)
        require(abs(coeff.imag)<TOL,'Matchete singlet projection complex')
        v.append(t-coeff.real*ident)
    pivot=next((a for a in v if np.linalg.norm(a)>1.e-9),None)
    require(pivot is not None,'Missing Matchete triplet adjoint direction')
    # Frobenius norm: sum_a Tr[tH^a tH^a] Tr[tS2^a tS2^a] = 3.
    cand=pivot*np.sqrt(3.)/np.linalg.norm(pivot)
    plus=_spectral(cand)
    minus=_spectral(-cand)
    passing=[sign for sign,x in ((1,plus),(-1,minus)) if x['max_residual']<TOL]
    require(len(passing)==1,
            'No unique sign passes the independent S2 total-Casimir spectrum')
    selected=passing[0]
    checks=plus if selected==1 else minus
    rejected=minus if selected==1 else plus
    require(rejected['spectrum_residual']>.1 and
            rejected['minimal_polynomial_residual']>.1,
            'Opposite sign unexpectedly passes S2 spectral negative control')
    # Both CGs must live in the rank-two singlet/adjoint span.
    worst=0.
    for t in raw:
        a=np.vdot(ident,t)/np.vdot(ident,ident)
        b=np.vdot(selected*cand,t)/np.vdot(cand,cand)
        require(max(abs(a.imag),abs(b.imag))<TOL,
                'Complex coefficient in Hermitian portal span')
        worst=max(worst,rel(t,a.real*ident+b.real*selected*cand))
    require(worst<TOL,'Portal CGs fail independent spectral-singlet span')
    return selected*cand, {
        'sign_relative_to_first_native_traceless_CG':selected,
        'selected_spectral_checks':checks,
        'opposite_sign_spectrum_residual':rejected['spectrum_residual'],
        'opposite_sign_polynomial_residual':rejected['minimal_polynomial_residual'],
        'native_Matchete_portal_span_residual':worst,
    }


def _v6_adjoint(raw: list[np.ndarray], sector: dict) -> np.ndarray:
    """Read v6 after sign selected independently; use only for cross-check."""
    require(sector.get('sector')=='triplet', 'Expected v6 triplet sector')
    key='candidate_rgbeta_from_matchete_real_matrix'  # v6 (NOT v5 key)
    require(key in sector,'v6 triplet coupling map missing')
    coeff=np.asarray(sector[key],float)
    require(coeff.shape==(2,2),'v6 triplet coupling map wrong shape')
    ident=_identity(3)
    candidates=[]
    for j in range(2):
        if abs(coeff[1,j])>1.e-9:
            candidates.append((raw[j]-coeff[0,j]*ident)/coeff[1,j])
    require(candidates,'v6 triplet map has no adjoint coefficient')
    require(max(rel(x,candidates[0]) for x in candidates)<TOL,
            'v6 triplet map inconsistent across CGs')
    require(max(rel(raw[j],coeff[0,j]*ident+coeff[1,j]*candidates[0])
                for j in range(2))<TOL,
            'v6 triplet map does not reconstruct native CGs')
    return candidates[0]


def inspect_witness(witness:dict, root:Path|None=None)->dict:
    require(witness.get('format')==NATIVE_FORMAT,'Unrecognized native v8 format')
    status=witness.get('status')
    require(status in (NATIVE_REDUCED,NATIVE_UNEXPOSED),
            'Unknown native v8 reduction result')
    checks=['fundamental_dimension_and_trace_valid',
            'S2_dimension_trace_and_casimir_valid',
            'native_s2_generator_still_symbolic',
            'native_delS2_still_symbolic',
            'native_S2_symmetric_isometry_valid']
    require(all(witness.get(k) is True for k in checks),
            'Missing required native representation identity')
    reduction=(witness.get('native_s2_generator_reduction_observed') is True and
               witness.get('native_delS2_reduction_observed') is True)
    require((status==NATIVE_REDUCED)==reduction,
            'Native reduction status contradicts measured witness flags')
    forbidden=['native_s2_numeric_components_exposed',
               'native_s2_numeric_frame_certified',
               'native_s2_to_matchete_intertwiner_certified',
               'physical_coupling_dictionary_certified',
               'physical_wilson_evolution_authorized']
    require(all(witness.get(k) is False for k in forbidden),
            'Native v8 claims unsupported numeric/physical certification')
    require(isinstance(witness.get('native_s2_generator_refined_inputform'),str) and
            isinstance(witness.get('native_delS2_refined_inputform'),str),
            'Missing native symbolic refinement records')
    if root is not None:
        hashes=witness.get('source_sha256')
        require(isinstance(hashes,dict) and set(hashes)==set(SOURCE_FILES),
                'Native v8 source hash inventory missing')
        for src in SOURCE_FILES:
            fp=root/src
            require(fp.is_file() and hashes[src]==_sha(fp),
                    'Stale native v8 source: '+src)
    return {'status':status,'symbolic_reduction_observed':reduction,
            'numeric_s2_frame_exposed':False}


def audit(mat:dict, v6:dict, native:dict)->dict:
    from Numerical.parallel_t3.wilson.portal_orientation_v6 import (
        STATUS as V6_STATUS, _projected_portals)
    require(v6.get('status')==V6_STATUS,'v6 status not verified')
    require(v6.get('triplet_sectors_Pauli_algebra_conditional')==TRIPLET_COUNT,
            'v6 triplet inventory incomplete')
    require(v6.get('scope',{}).get('native_RGBeta_triplet_intertwiner_certified') is False,
            'v6 claims native intertwiner')
    w=inspect_witness(native)
    models={};count=0;worst=0.
    for tag in sorted(v6.get('models',{})):
        require(tag in mat.get('models',{}),'Missing Matchete source model '+tag)
        sectors=[]
        for vs in v6['models'][tag]['sectors']:
            if vs.get('sector')!='triplet':continue
            which=vs.get('which')
            source=sorted((p for p in mat['models'][tag]['portals']
                           if p.get('which')==which),key=lambda p:p['invariant'])
            raw=_projected_portals(source,which,3)
            selected,checks=select_spectral_sign(raw)
            from_v6=_v6_adjoint(raw,vs)
            v6agreement=rel(selected,from_v6)
            require(v6agreement<TOL,
                    f'{tag}/S{which}: basis-independent spectrum contradicts v6')
            sector={
                'model':tag,'which':which,
                'test':'total_SU2_Casimir_spectrum_NO_Matchete_fitted_intertwiner',
                'source_CG_sign_selected_without_using_v6_map':True,
                'comparison_to_v6_selected_adjoint_residual':v6agreement,
                'basis_independent_eigenvalues_expected':EXPECTED_EIGENVALUES.tolist(),
                **checks,
                'RGBeta_numeric_triplet_frame_certified':False,
            }
            worst=max(worst,v6agreement,checks['selected_spectral_checks']['max_residual'],
                      checks['native_Matchete_portal_span_residual'])
            count+=1
            sectors.append(sector)
        models[tag]=sectors
    require(len(models)==16 and count==TRIPLET_COUNT,
            'Expected 16 models and 14 triplet sectors')
    require(worst<TOL,'Unacceptable basis-independent portal spectral residual')
    return {
        'format':FORMAT,
        'status':STATUS_REDUCED if w['symbolic_reduction_observed'] else STATUS_UNEXPOSED,
        'native_RGBeta_source_reduction':w,
        'triplet_sectors_basis_invariant_spectrum_verified':count,
        'maximum_spectral_and_v6_comparison_residual':worst,
        'model_sectors':models,
        'certification':{
            'S2_total_Casimir_eigenvalue_pattern_source_grounded':True,
            'native_RGBeta_symbolic_reduction_directly_measured':w['symbolic_reduction_observed'],
            'native_RGBeta_numerical_S2_frame_certified':False,
            'independent_native_RGBeta_to_Matchete_intertwiner_certified':False,
            'physical_complete_scalar_potential_normalization_certified':False,
            'BCminus1_mixed_scalar_mass_diagonalization_certified':False,
            'signed_GammaCC_to_Weyl_certified':False,
            'physical_Wilson_evolution_authorized':False,
            'EFT_order':'UV -> integrate out all fermions -> integrate out all scalars -> SMEFT'
        },
    }


def audit_from_paths(mat_path:Path,rgb_path:Path,v5_path:Path,
                     fierz_path:Path,v6_path:Path,v7_path:Path,
                     native_v8_path:Path,root:Path)->dict:
    """Revalidate live v5-v7 provenance; compare saved v7 before v8."""
    stored_v7=read_json(v7_path)
    require(stored_v7.get('status')==V7_STATUS and
            stored_v7.get('triplet_sectors_conditional_comparison_verified')==14,
            'Saved v7 is not the completed conditional 14-sector audit')
    require(stored_v7.get('independent_native_RGBeta_S2_to_Matchete_intertwiner_certified') is False,
            'Saved v7 promotes unproven native S2 frame')
    v7fresh=v7_revalidate(mat_path,rgb_path,v5_path,fierz_path,v6_path,root)
    require(v7fresh['status']==stored_v7['status'] and
            rel(v7fresh['maximum_comparison_residual'],stored_v7['maximum_comparison_residual'])<TOL and
            set(v7fresh['models'])==set(stored_v7['models']),
            'Saved v7 differs from live source recheck')
    for tag in v7fresh['models']:
        old=stored_v7['models'][tag];new=v7fresh['models'][tag]
        require(len(old)==len(new),'Saved v7 sectors changed '+tag)
        for a,b in zip(old,new):
            require(a.get('which')==b['which'] and
                    rel(a['adjoint_tensor_transport_residual'],b['adjoint_tensor_transport_residual'])<TOL,
                    'Saved v7 sector changed '+tag)
    native=read_json(native_v8_path)
    inspect_witness(native,root)
    # v7 revalidation does not output source Matchete itself: reread native data.
    result=audit(read_json(mat_path),read_json(v6_path),native)
    result['input_sha256']={name:_sha(path) for name,path in (
        ('matchete_native_v5',mat_path),('rgbeta_native_v5',rgb_path),
        ('portal_dictionary_v5',v5_path),('native_fierz_v6',fierz_path),
        ('portal_orientation_v6',v6_path),('s2_intertwiner_v7',v7_path),
        ('native_s2_reduction_v8',native_v8_path))}
    return result
