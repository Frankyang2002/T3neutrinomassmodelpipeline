"""Trace actual Stream-2 LLHH Ward failure to charge sectors / saved endpoints.

Diagnostic only: neither Wilson sources nor SMEFT matching are modified.
Reference scalar convention: H=(h_R + i*h_I)/sqrt(2), interleaved.
The real-to-(H,H*) transform is unitary, so squared Frobenius norms
partition exactly into HH, HHbar + HbarH, and HbarHbar sectors.
"""
from __future__ import annotations

from hashlib import sha256
from pathlib import Path
from typing import Any
import numpy as np

from .matching import extract_lepton_sectors, master_weinberg_template
from .models import n_real, normalize_model
from .pipeline import _load_pair, W_FORMAT
from .rough_projection import diagnose_llhh

FORMAT = 't3-rough-llhh-ward-origin-v5'


def complex_higgs_coordinate_map() -> np.ndarray:
    """Real Higgs phi_A = M[A,k] (H1,H2,H1*,H2*)[k]."""
    r = np.zeros((4, 4), complex)
    for j in range(2):
        r[2*j,j] = r[2*j,j+2] = 1 / np.sqrt(2)
        r[2*j+1,j] = -1j / np.sqrt(2)
        r[2*j+1,j+2] = 1j / np.sqrt(2)
    return r


def real_to_complex_llhh(a: np.ndarray) -> np.ndarray:
    """Convert BOTH real Higgs legs, not just their numerical names."""
    a = np.asarray(a, complex)
    if a.shape != (3,3,2,2,4,4) or not np.isfinite(a).all():
        raise ValueError('Expected finite LLHH [3,3,2,2,4,4]')
    M=complex_higgs_coordinate_map()
    if not np.allclose(M.conj().T @ M, np.eye(4), atol=1e-14):
        raise RuntimeError('Real/complex Higgs map is not unitary')
    return np.einsum('pqijAB,Ak,Bl->pqijkl',a,M,M,optimize=True)


def complex_to_real_llhh(a: np.ndarray) -> np.ndarray:
    """Inverse complex component conversion; tests can verify roundtrip."""
    a=np.asarray(a,complex)
    if a.shape!=(3,3,2,2,4,4):
        raise ValueError('Expected [3,3,2,2,4,4] holomorphic-coordinate LLHH')
    M=complex_higgs_coordinate_map()
    return np.einsum('pqijkl,Ak,Bl->pqijAB',a,M.conj(),M.conj(),optimize=True)


def _weak_ward_hh(hh: np.ndarray) -> dict[str,float]:
    # Same SU(2) convention as existing wilson.basis.spin_generators(2).
    gens={
        'T1':np.array([[0,.5],[.5,0]],complex),
        'T2':np.array([[0,-.5j],[.5j,0]],complex),
        'T3':np.diag([.5,-.5]).astype(complex),
    }
    n=max(float(np.linalg.norm(hh)),1e-100)
    result={}
    for name,g in gens.items():
        v=(np.einsum('ki,pqkjab->pqijab',g,hh,optimize=True)
           +np.einsum('kj,pqikab->pqijab',g,hh,optimize=True)
           +np.einsum('ca,pqijcb->pqijab',g,hh,optimize=True)
           +np.einsum('cb,pqijac->pqijab',g,hh,optimize=True))
        result[name]=float(np.linalg.norm(v)/n)
    return result


def source_charge_analysis(a: np.ndarray) -> dict[str,Any]:
    """Unitary decomposition and independent SU(2)/U(1) identities.

    Only the HH sector carries the Weinberg operator for two ordinary
    hypercharge -1/2 lepton doublets.  U(1) total charges of HH, H H*,
    H*H* are respectively 0, -1, -2.  These are tensor conventions, not
    new finite matching terms and not a projection applied to the source.
    """
    a=np.asarray(a,complex)
    z=real_to_complex_llhh(a)
    norm=float(np.linalg.norm(z))
    h=z[..., :2,:2]
    hb=z[..., :2,2:]
    bh=z[..., 2:,:2]
    bb=z[..., 2:,2:]
    h_norm=float(np.linalg.norm(h))
    hb_norm=float(np.sqrt(np.linalg.norm(hb)**2 + np.linalg.norm(bh)**2))
    bb_norm=float(np.linalg.norm(bb))
    reconstructed=float(np.sqrt(h_norm**2+hb_norm**2+bb_norm**2))
    if abs(reconstructed-norm)>1e-10*max(norm,1e-100):
        raise RuntimeError('Higgs charge sector norms do not reconstruct original LLHH')
    Tz=real_to_complex_llhh(np.einsum('pq,ijAB->pqijAB',np.eye(3),master_weinberg_template()))
    t=Tz[0,0,...,:2,:2]
    c5=np.einsum('pqijab,ijab->pq',h,t.conj(),optimize=True)/np.vdot(t,t)
    c5=(c5+c5.T)/2
    mismatch=float(np.linalg.norm(h-np.einsum('pq,ijab->pqijab',c5,t,optimize=True)))
    from .rough_projection import _weyl_doublet_generators, _real_higgs_generator, _ward
    u1=_ward(a,_weyl_doublet_generators()['U1_Y'],
             _real_higgs_generator(.5*np.eye(2)))
    actual_u1=float(np.linalg.norm(u1)/max(norm,1e-100))
    expected_u1=float(np.sqrt(hb_norm**2+4*bb_norm**2)/max(norm,1e-100))
    if abs(actual_u1-expected_u1)>1e-9*max(1.,actual_u1,expected_u1):
        raise RuntimeError(f'Hypercharge sector decomposition inconsistent with direct Ward ({actual_u1}, {expected_u1})')
    # Since real/complex rotation is unitary, source and transformed norms agree.
    if abs(norm-np.linalg.norm(a))>1e-9*max(norm,1e-100):
        raise RuntimeError('Wrong real-to-complex tensor normalization')
    if norm<1e-100:
        disposition='zero_LLHH_no_relative_source_information'
    elif hb_norm/norm>1e-8 or bb_norm/norm>1e-8:
        disposition='wrong_hypercharge_sector_in_LLHH'
    elif mismatch/max(h_norm,1e-100)>1e-8:
        disposition='correct_HH_hypercharge_but_broken_SU2_tensor'
    else:
        disposition='consistent_with_unique_Weinberg_singlet'
    return {
        'status':disposition,
        'full_real_LLHH_norm_gev_inverse':norm,
        'full_complex_LLHH_norm_gev_inverse':norm,
        'charge_sector_norms_gev_inverse':{'HH':h_norm,'H_Hbar_plus_Hbar_H':hb_norm,'Hbar_Hbar':bb_norm},
        'charge_sector_fractions':{'HH':h_norm/max(norm,1e-100),
                                   'H_Hbar_plus_Hbar_H':hb_norm/max(norm,1e-100),
                                   'Hbar_Hbar':bb_norm/max(norm,1e-100)},
        'independently_reconstructed_U1_Ward_relative':expected_u1,
        'direct_U1_Ward_relative':actual_u1,
        'allowed_HH_SU2_Ward_relative':_weak_ward_hh(h),
        'allowed_HH_singlet_residual_relative':mismatch/max(h_norm,1e-100),
        'projected_C5_from_allowed_HH':{'real':c5.real.tolist(),'imag':c5.imag.tolist()},
        'NO_auto_repair':True,
    }


def _get_llhh(c:np.ndarray, meta:dict, model:str)->np.ndarray:
    h,_=extract_lepton_sectors(c,np.asarray(meta['lepton_indices'],int),
         np.asarray(meta['higgs_source_indices'],int),
         np.asarray(meta['heavy_source_indices'],int),
         number_of_real_scalars=4+n_real(model))
    return h


def audit_pair(manifest:Path, *, model:str) -> dict:
    """Audit a C-owned adapted Stream-2 endpoint with checksum validation.

    No guessed MF boundary: order-zero MF is examined only if actually saved.
    If only an order-zero boundary is present it cannot certify the actual
    one-loop F-threshold or pinpoint which RGE term caused MS contamination.
    """
    manifest=Path(manifest)
    name=normalize_model(model)
    meta,arrays,meta_sha,npz_sha=_load_pair(manifest,W_FORMAT)
    if normalize_model(meta.get('model',''))!=name:
        raise ValueError('Source model mismatch')
    w=_get_llhh(arrays['C_psi2phi2_ms'],meta,name)
    evolved=source_charge_analysis(w)
    purity=diagnose_llhh(w,master_weinberg_template())
    boundary=arrays.get('C_psi2phi2_tree_mf')
    if boundary is None:
        mf={'availability':'missing'}
        inferred='missing_MF_witness_cannot_localise'
    else:
        x=_get_llhh(boundary,meta,name)
        mf={'availability':'C_psi2phi2_tree_mf',**source_charge_analysis(x)}
        if mf['full_real_LLHH_norm_gev_inverse']<1e-100:
            inferred='tree_MF_LLHH_zero;_MS_source_requires_beta_or_hard_boundary_audit'
        elif mf['status']=='consistent_with_unique_Weinberg_singlet':
            inferred='tree_MF_LLHH_singlet;_MS_source_noninvariance_arose_after_saved_tree_endpoint'
        else:
            inferred='tree_MF_LLHH_already_noninvariant;_check_boundary_embedding_or_map'
    if evolved['status']=='consistent_with_unique_Weinberg_singlet':
        inferred='MS_LLHH_singlet;_no_source_Ward_blocker'
    return {'format':FORMAT,'model':name,
            'status':'SourceWardOriginDiagnosticOnly',
            'source_manifest':str(manifest),
            'source_sha256':{'manifest':meta_sha,'npz':npz_sha},
            'source_scalar_basis_id':meta.get('scalar_basis_id'),
            'source_lepton_indices':meta.get('lepton_indices'),
            'source_higgs_indices':meta.get('higgs_source_indices'),
            'synthetic_fixture':bool(meta.get('synthetic_fixture')),
            'available_NPZ_keys':sorted(arrays),
            'MS_full_evolved_LLHH':evolved,
            'MS_Weinberg_purity_relative':purity['relative_Weinberg_purity_residual'],
            'MF_orderzero_LLHH':mf,
            'origin_interpretation_NOT_proof':inferred,
            'does_not_alter_existing_C_preflight':True,
            'no_scalar_matching_or_SMEFT_running_performed':True,
            'omitted_independent_matching_terms_not_inferred':True}
