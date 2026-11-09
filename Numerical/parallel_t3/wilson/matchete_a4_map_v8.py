"""A(-4) exact-source Matchete T3MixCG triplet frame -> canonical SU(2) spin-1 frame.

A fixed, *unitary* basis transformation obtained from the archived Matchete CG,
not a fit to the Wilson boundary and not a Ward-singlet projection. The input
Matchete PL term contains (L L) S1* S2, with NewScalar2 in the native 3-basis.
The initial rough Stream1 export erroneously identifies this native basis with
the canonical spherical spin-1 coordinates. This module only corrects that
explicit scalar index map; it does not certify GammaCC/Weyl normalization,
additional UV hard operators, or the complete scalar coupling dictionary.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import numpy as np

from .basis import build_basis
from .coordinates import embed_formal_llss
from .ward_singlet_v7 import _gauge_norm

MODEL = 'A_-4'
CLASSIFICATION = 'SourceCGMappedA4NumericalCandidateNotCertified'
SOURCE_REL = ('output/general_observable/fermion_first_hermitian_portal_v1/'
              'A_m4/data/eft1_after_F_wilson_seed.json')
# Exact six-entry sparse Mathematica CGRegistry from this source (whitespace removed).
SOURCE_SPARSE = (
    'SparseArray[Automatic,{2,2,3},0,{1,{{0,3,6},'
    '{{1,1},{1,3},{2,2},{1,2},{2,1},{2,3}}},'
    '{3^(-1/4),I/3^(1/4),-3^(-1/4),-3^(-1/4),'
    '-3^(-1/4),I/3^(1/4)}}]'
)
SOURCE_REPS = '{Bar[SU2L[fund]],Bar[SU2L[fund]],T3BSMd3}'


def sha256(path:Path) -> str:
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()


def source_cg(seed_path:Path) -> tuple[np.ndarray,dict]:
    """Read the source file and fail closed unless its exact CG is recognised."""
    path=Path(seed_path)
    data=json.loads(path.read_text(encoding='utf-8-sig'))
    if data.get('IntegratedFields')!=['F'] or data.get('ActiveHeavyFields')!=['S1','S2']:
        raise ValueError('Matchete source does not integrate ALL F with BOTH S active')
    terms=data.get('TreeWilsonTerms',[])
    if len(terms)!=2 or not any(
        t.get('Chirality')=='PL' and 'CG[T3MixCG' in t.get('TermInputForm','')
        and 'Bar[Field[NewScalar1' in t.get('TermInputForm','')
        and 'Field[NewScalar2' in t.get('TermInputForm','')
        for t in terms):
        raise ValueError('Expected Matchete A(-4) PL LL S1* S2 source term is absent')
    matches=[r for r in data.get('CGRegistry',[]) if r.get('Name')=='T3MixCG']
    if len(matches)!=1:
        raise ValueError('Require exactly one native T3MixCG registry entry')
    record=matches[0]
    compact=lambda s: ''.join(str(s).split())
    if compact(record.get('TensorInputForm'))!=SOURCE_SPARSE or \
       compact(record.get('RepsInputForm'))!=SOURCE_REPS:
        raise ValueError('Native T3MixCG changed: cannot reuse source-derived triplet map')
    # SparseArray[Automatic,{2,2,3},0,{1,{{0,3,6},...},values}]:
    # first index=lepton_1, inner pairs=(lepton_2, native_triplet_slot).
    z=3**(-.25)
    g=np.zeros((2,2,3),complex)
    g[0,0,0]=z;g[0,0,2]=1j*z;g[0,1,1]=-z
    g[1,0,1]=-z;g[1,1,0]=-z;g[1,1,2]=1j*z
    return g,{'source_path':str(path),'source_sha256':sha256(path),
              'registry':'T3MixCG','native_reps':record['RepsInputForm'],
              'source_chirality_used':'PL','source_scalar_term':'LL S1_conjugate S2'}


def canonical_singlet_cg() -> np.ndarray:
    """Canonical L L S2(spin-1) invariant, unit Frobenius norm."""
    t=np.zeros((2,2,3),complex)
    t[0,0,2]=1/np.sqrt(3)
    t[0,1,1]=t[1,0,1]=-1/np.sqrt(6)
    t[1,1,0]=1/np.sqrt(3)
    return t


def derive_triplet_unitary(g:np.ndarray) -> tuple[np.ndarray,float]:
    """Solve G_native @ U_native_from_canonical = c * T_spherical.

    This solves for a *unique* 3x3 unitary from the source CG, not from any
    user's numerical LLSS tensor. Requires all of the original CG norm.
    """
    g=np.asarray(g,complex)
    t=canonical_singlet_cg()
    if g.shape!=(2,2,3) or not np.isfinite(g).all():
        raise ValueError('Invalid Matchete T3MixCG tensor')
    coeff,_,rank,_=np.linalg.lstsq(g.reshape(4,3),t.reshape(4,3),rcond=None)
    if rank!=3:
        raise ValueError('Native triplet CG rank is not three')
    left,_,right=np.linalg.svd(coeff)
    u=left@right
    overlap=np.vdot(t,g.reshape(4,3)@u)
    u*=np.exp(-1j*np.angle(overlap))
    scale=float(abs(overlap)/np.vdot(t,t).real)
    rel=float(np.linalg.norm(g.reshape(4,3)@u-scale*t.reshape(4,3)) /
              max(np.linalg.norm(g),1e-100))
    if np.linalg.norm(u.conj().T@u-np.eye(3))>1e-12 or rel>1e-12:
        raise ValueError(f'Source CG not unitarily equivalent to canonical triplet: {rel}')
    return u,scale


def realify_complex_unitary(u:np.ndarray) -> np.ndarray:
    """Real, interleaved R/I action for z_native=U z_canonical."""
    u=np.asarray(u,complex)
    if u.ndim!=2 or u.shape[0]!=u.shape[1] or not np.isfinite(u).all():
        raise ValueError('Expected square finite complex unitary')
    if not np.allclose(u.conj().T@u,np.eye(len(u)),atol=1e-12):
        raise ValueError('Nonunitary source scalar map')
    d=len(u);r=np.zeros((2*d,2*d),float)
    for i in range(d):
        for j in range(d):
            x,y=u[i,j].real,u[i,j].imag
            r[2*i:2*i+2,2*j:2*j+2]=[[x,-y],[y,x]]
    if np.linalg.norm(r.T@r-np.eye(2*d))>1e-11:
        raise ValueError('Native complex U does not realify to an orthogonal map')
    return r


def source_heavy_from_canonical(seed_path:Path) -> tuple[np.ndarray,dict]:
    g,source=source_cg(seed_path)
    u,scale=derive_triplet_unitary(g)
    r=np.eye(8)
    r[2:,2:]=realify_complex_unitary(u)
    source.update({
        'map_definition':'phi_native[a]=sum_A R[a,A] phi_canonical[A]',
        'singlet_S1_real_pair':'identity',
        'triplet_S2_native_from_canonical_U_real':u.real.tolist(),
        'triplet_S2_native_from_canonical_U_imag':u.imag.tolist(),
        'source_CG_to_canonical_scalar_factor':scale,
        'map_orthogonality_error':float(np.linalg.norm(r.T@r-np.eye(8))),
        'reconstruction_relative_error':float(np.linalg.norm(
            g.reshape(4,3)@u-scale*canonical_singlet_cg().reshape(4,3)) /
            np.linalg.norm(g)),
        'no_coefficient_projection':True,
        'no_norm_discarded':True,
    })
    return r,source


def convert_a4_boundary(initial:np.ndarray, raw:np.ndarray, factor:complex,
                        r:np.ndarray, *, ward_tolerance:float=1e-9
                        ) -> tuple[np.ndarray,dict]:
    """Apply the source-derived orthogonal map to *all* native LLSS slots.

    Preserves any original hard C5/non-LLSS blocks; rejects any mismatched
    baseline and does not substitute an invariant template for source data.
    """
    b=build_basis(MODEL)
    c=np.asarray(initial,complex)
    s=np.asarray(raw,complex)
    if c.shape!=(b.nf,b.nf,b.ns,b.ns) or s.shape!=(3,3,2,2,8,8):
        raise ValueError('Incorrect original A(-4) Wilson/formal tensor shapes')
    if np.shape(r)!=(8,8) or np.linalg.norm(np.asarray(r).T@r-np.eye(8))>1e-10:
        raise ValueError('Invalid A(-4) source-to-canonical R')
    identity=embed_formal_llss(s,b,np.eye(8),factor)
    updated=embed_formal_llss(s,b,r,factor)
    ls=[i for pair in b.leptons for i in pair]
    heavy=list(range(4,b.ns))
    ix=np.ix_(ls,ls,heavy,heavy)
    in_norm=float(np.linalg.norm(identity[ix]))
    if in_norm<=1e-100:raise ValueError('Zero A(-4) LLSS boundary')
    discrepancy=float(np.linalg.norm(c[ix]-identity[ix])/in_norm)
    if discrepancy>1e-8:
        raise ValueError(f'Original Stream2 boundary does not match saved Matchete formal LLSS: {discrepancy:g}')
    remainder=c-identity
    if np.linalg.norm(remainder[ix])>1e-8*in_norm:
        raise ValueError('Unexpected extra LLSS terms cannot be silently transformed')
    mapped=updated+remainder
    # This operation is an orthogonal change of scalar frame, not a projection.
    out_norm=float(np.linalg.norm(mapped[ix]))
    norm_delta=abs(out_norm-in_norm)/in_norm
    if norm_delta>1e-10:raise ValueError('Triplet map lost original Wilson tensor norm')
    original_outside=c.copy();original_outside[ix]=0
    mapped_outside=mapped.copy();mapped_outside[ix]=0
    if not np.array_equal(original_outside,mapped_outside):
        raise ValueError('Source map changed a non-LLSS operator')
    original_ward=_gauge_norm(c,b)
    mapped_ward=_gauge_norm(mapped,b)
    if max(mapped_ward.values())>ward_tolerance:
        raise ValueError('Source-CG map does not restore the full gauge Ward identity: '
                         +repr(mapped_ward))
    return mapped,{
        'model':MODEL,'input_identity_ward':original_ward,
        'source_cg_mapped_ward':mapped_ward,
        'source_LLSS_norm_GeV_inverse':in_norm,
        'mapped_LLSS_norm_GeV_inverse':out_norm,
        'LLSS_norm_relative_change':norm_delta,
        'LLSS_norm_fraction_retained':out_norm/in_norm,
        'original_stream2_boundary_reconstruction_error':discrepancy,
        'outside_LLSS_unchanged':True,
        'no_Ward_singlet_projection':True,
        'classification':CLASSIFICATION,
        'signed_GammaCC_to_Weyl_not_independently_certified':True,
    }
