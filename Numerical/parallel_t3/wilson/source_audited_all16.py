"""Source-audited conditional T3 tensor bridge, for the 16-model rough benchmark.

Only CG frame changes derived from the tracked Matchete native tensors are
applied to the *Wilson* boundary. Quartics are built in a canonical SU(2)
representation from the explicitly registered RGBeta invariant structures.
Their native coefficient normalisations remain CONDITIONAL, not certified.
No Ward-singlet projection is performed. No UV/SMEFT dynamics are redefined.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from Numerical.parallel_t3.rough_inputs import builder
from .basis import ALL_MODELS, build_basis, spin_generators
from Numerical.parallel_t3.threshold.models import MIXED, SELF_S1, SELF_S2, _mixed_rotation
from .coordinates import embed_formal_llss
from .matchete_a4_map_v8 import (
    SOURCE_REL, SOURCE_REPS, SOURCE_SPARSE, source_heavy_from_canonical,
    realify_complex_unitary,
)
from .ward_singlet_v7 import _gauge_norm
from Numerical.parallel_t3.threshold.rough_ward_gate_v6 import gauge_ward_report

# The native E-family CG was inspected at GitHub main 9e50ca77e64e9d5fb37079f325bc2b8cab5112e9.
# Source-basis indices are first/second lepton, S1, S2, respectively.
_E_INDICES = (
    (1,1,1,2),(1,1,2,1),(1,1,2,3),(1,1,3,2),
    (1,2,1,3),(1,2,3,1),(2,1,1,3),(2,1,3,1),
    (2,2,1,2),(2,2,2,1),(2,2,2,3),(2,2,3,2),
)
_E_VALUES = (1,-1,-1j,1j,-1j,1j,-1j,1j,1,-1,1j,-1j)
_E_SPARSE = (
    'SparseArray[Automatic,{2,2,3,3},0,{1,{{0,6,12},'
    '{{1,1,2},{1,2,1},{1,2,3},{1,3,2},{2,1,3},{2,3,1},'
    '{1,1,3},{1,3,1},{2,1,2},{2,2,1},{2,2,3},{2,3,2}}},'
    '{1/Sqrt[2],-(1/Sqrt[2]),(-I)/Sqrt[2],I/Sqrt[2],'
    '(-I)/Sqrt[2],I/Sqrt[2],(-I)/Sqrt[2],I/Sqrt[2],'
    '1/Sqrt[2],-(1/Sqrt[2]),I/Sqrt[2],(-I)/Sqrt[2]}}]'
)
_E_REPS = '{Bar[SU2L[fund]],Bar[SU2L[fund]],T3BSMd3,T3BSMd3}'
_BASE_QUARTIC_UNITS = builder.quartic_unit_tensors


def _compact(s: Any) -> str:
    return ''.join(str(s).split())


def hash_file(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def _cg_record(seed: dict, name: str) -> dict:
    records = [r for r in seed.get('CGRegistry', []) if r.get('Name') == name]
    if len(records) != 1:
        raise ValueError(f'Expected one Matchete {name} source CG; got {len(records)}')
    return records[0]


def _e_native_cg() -> np.ndarray:
    g = np.zeros((2, 2, 3, 3), complex)
    for ijab, value in zip(_E_INDICES, _E_VALUES):
        g[tuple(k-1 for k in ijab)] = value / np.sqrt(2)
    return g


def _su2_variation(t: np.ndarray, d1: int, d2: int,
                   *, first_conjugated: bool) -> float:
    """Complex-component SU2 Ward check on LL S1* S2 (or HH S1 S2*)."""
    scale = max(float(np.linalg.norm(t)), 1e-100)
    residual = 0.
    for th, t1, t2 in zip(spin_generators(2), spin_generators(d1), spin_generators(d2)):
        z = (np.einsum('ix,xjab->ijab', th, t)
             + np.einsum('jx,ixab->ijab', th, t))
        if first_conjugated:
            z -= np.einsum('xa,ijxb->ijab', t1, t)
            z += np.einsum('bx,ijax->ijab', t2, t)
        else:
            z += np.einsum('ax,ijxb->ijab', t1, t)
            z -= np.einsum('xb,ijax->ijab', t2, t)
        residual = max(residual, float(np.linalg.norm(z)) / scale)
    return residual


def _triplet_unitary(root: Path) -> tuple[np.ndarray, dict]:
    from .matchete_a4_map_v8 import source_cg, derive_triplet_unitary
    anchor = Path(root) / SOURCE_REL
    g, evidence = source_cg(anchor)
    u, scale = derive_triplet_unitary(g)
    evidence['CG_normalisation_comparison_only'] = scale
    evidence['triplet_map_anchor'] = SOURCE_REL
    return u, evidence


def _self_named_from_interleaved(dim: int) -> np.ndarray:
    """GeneralT3OrderedTensor.scalar_basis packs Minus[all], Plus[all].

    The canonical code uses (Re1, Im1, Re2, Im2, ...). No projection.
    """
    P = np.zeros((2*dim, 2*dim))
    for k in range(dim):
        P[k,2*k] = 1.
        P[dim+k,2*k+1] = 1.
    return P


def _complex_json_matrix(payload: dict) -> np.ndarray:
    return np.asarray(payload['real'],float)+1j*np.asarray(payload['imag'],float)


def _mixed_f_source_rotation(root: Path, model: str, source: dict) -> tuple[np.ndarray,dict]:
    """For B/C(-1) use the exact same-run B source's MF mixing U.

    GeneralT3MassEigenMatchete.wl defines
      (S1, eps conj(S2)) = U (NewScalarX, NewScalarY).
    The B fermion boundary (NOT the old standalone Matchete seed) is already
    expressed in the X/Y basis and stores the U that produced its coefficients.
    """
    if (source.get('format')!='t3-mixed-factorised-LLSS-source-v1'
            or source.get('model')!=model or not source.get('operator_terms')):
        raise ValueError('B/C mixed source is not the archived X/Y ordered source')
    from Numerical.parallel_t3.background.contracts import tag
    base = root/'output/parallel_t3/background/B_all16_qualitative_v1'/tag(model)
    meta = json.loads((base/'fermion_boundary_manifest.json').read_text(encoding='utf-8-sig'))
    if (meta.get('model')!=model or meta.get('source_matchete_seed_sha256')
            !=hash_file(root/meta['source_matchete_seed'].replace('\\','/'))
            or meta.get('mixed_scalar_mass_matrix_scale')!='MF'):
        raise ValueError('Mixed X/Y same-run B source or MF mass-basis provenance mismatch')
    diag=meta.get('mixed_native_diagnostics',{})
    U=_complex_json_matrix(diag['mixing_U_muS'])
    if U.shape!=(2,2) or not np.isfinite(U).all() or not np.allclose(U.conj().T@U,np.eye(2),atol=5e-9):
        raise ValueError('B/C source mixing U is not a unitary 2x2 matrix')
    # models._mixed_rotation: canonical real S1/S2 = R @ native real X/Y.
    # Wilson source index mapping uses native = R.T @ canonical.
    R=_mixed_rotation(U).T
    return R, {'basis':'(S1,epsilon conjugate(S2))=U@(X,Y)',
               'mixing_scale':'MF', 'mixing_U_real':U.real.tolist(),
               'mixing_U_imag':U.imag.tolist(),
               'B_fermion_manifest_sha256':hash_file(base/'fermion_boundary_manifest.json'),
               'epsilon':[[0,1],[-1,0]],
               'original_XY_source_sha256':meta['source_matchete_seed_sha256']}


def source_scalar_map(root: Path, model: str, seed_path: Path) -> tuple[np.ndarray, dict]:
    """Use source-anchored orthogonal Matchete-native -> canonical scalar maps.

    In the six self-bilinear cases, source real slots use a grouped real
    (Minus[1..d],Plus[1..d]) representation, which is NOT interleaved.
    B(-1)/C(-1) use the full same-run mixed scalar U and epsilon intertwiner.
    """
    basis = build_basis(model)
    d1,d2,_ = basis.dims
    p=Path(seed_path); root=Path(root)
    seed=json.loads(p.read_text(encoding='utf-8-sig'))
    R=np.eye(2*d1+2*d2)
    fam=model.split('_')[0]
    evidence:dict[str,Any]={
        'model':model, 'Matchete_source':str(p), 'Matchete_sha256':hash_file(p),
        'source_CG_names':[r.get('Name') for r in seed.get('CGRegistry',[])],
        'definition':'phi_Matchete_native = R @ phi_canonical',
        'source_CG_orientation_verified':False, 'no_Wilson_projection':True,
    }
    if model in MIXED:
        R,mix=_mixed_f_source_rotation(root,model,seed)
        evidence.update(mix)
        evidence['source_CG_orientation_verified']='same_run_MF_XY_unitary_and_epsilon'
    else:
        if seed.get('IntegratedFields')!=['F'] or seed.get('ActiveHeavyFields')!=['S1','S2']:
            raise ValueError('Incorrect source EFT ordering: require ALL F out, BOTH S active')
        if fam in ('B','C'):
            if not seed.get('CGRegistry'):
                raise ValueError('Missing B/C native CG registry')
            reps=[_compact(x.get('RepsInputForm')) for x in seed['CGRegistry']]
            if not any('T3BSMd2' in x for x in reps) or any('T3BSMd3' in x for x in reps):
                raise ValueError('B/C native scalar representation inconsistent')
            evidence['source_CG_orientation_verified']='identity_if_full_Ward_passes'
        else:
            u,anchor=_triplet_unitary(root)
            evidence['triplet_anchor']=anchor
            record=_cg_record(seed,'T3MixCG')
            if fam in ('A','D'):
                if (_compact(record.get('TensorInputForm'))!=SOURCE_SPARSE
                        or _compact(record.get('RepsInputForm'))!=_compact(SOURCE_REPS)):
                    raise ValueError('A/D native triplet CG differs from exact source anchor')
                if fam=='A':R[2*d1:,2*d1:]=realify_complex_unitary(u)
                else:R[:2*d1,:2*d1]=realify_complex_unitary(u)
                evidence['source_CG_orientation_verified']='exact_A4_CG_triplet_representation'
            elif fam=='E':
                if (_compact(record.get('TensorInputForm'))!=_E_SPARSE
                        or _compact(record.get('RepsInputForm'))!=_E_REPS):
                    raise ValueError('E native two-triplet CG differs from archived source')
                g=_e_native_cg()
                canon=np.einsum('ijkl,ka,lb->ijab',g,u.conj(),u)
                residual=_su2_variation(canon,3,3,first_conjugated=True)
                if residual>1e-11:
                    raise ValueError(f'E source CG not SU2 covariant: {residual}')
                R[:6,:6]=realify_complex_unitary(u)
                R[6:,6:]=realify_complex_unitary(u)
                evidence['E_native_CG_canonical_SU2_residual']=residual
                evidence['source_CG_orientation_verified']='two_triplet_CG_exact_check'
            else:
                raise ValueError('Unknown T3 family')
        if model in SELF_S1 or model in SELF_S2:
            neutral_s1=model in SELF_S1
            dim=d1 if neutral_s1 else d2
            start=0 if neutral_s1 else 2*d1
            P=_self_named_from_interleaved(dim)
            R[start:start+2*dim,:] = P@R[start:start+2*dim,:]
            evidence['source_self_bilinear_field_packing']='Minus_all_then_Plus_all'
            evidence['self_bilinear_native_slot_reordering']=P.tolist()
            evidence['self_bilinear_species']='S1' if neutral_s1 else 'S2'
    if np.linalg.norm(R.T@R-np.eye(len(R)))>1e-10:
        raise ValueError('Nonorthogonal source scalar frame map')
    evidence['orthogonality_residual']=float(np.linalg.norm(R.T@R-np.eye(len(R))))
    evidence['map_sha256']=hashlib.sha256(R.tobytes()).hexdigest()
    return R,evidence


def covariance_corrected_mass(root:Path,model:str,original:np.ndarray
                               ) -> tuple[np.ndarray,dict]:
    """Preserve the entire source mass spectrum and the neutral-triplet
    complex bilinear, converting the *Cartesian* self-irrep metric into the
    canonical spherical-frame metric U.T@U. No fit to Wilson data.

    The v1 scalar adapter used three independent real/imag 2x2 blocks for a
    self-conjugate triplet, equivalent to a Cartesian delta. Transforming
    that scalar quadratic to the same canonical spin-1 frame as quartics
    requires the invariant metric U.T@U; the eigenvalues are unchanged.
    """
    b=build_basis(model);d1,d2,_=b.dims
    m=np.asarray(original,float).copy()
    if m.shape!=(2*(d1+d2),)*2 or not np.allclose(m,m.T,atol=1e-7):
        raise ValueError('Invalid original source scalar mass matrix')
    details={'original_mass_eigenvalues_sq_gev2':np.linalg.eigvalsh(m).tolist(),
             'corrected':False}
    if model in SELF_S1 or model in SELF_S2:
        dim=d1 if model in SELF_S1 else d2
        start=0 if model in SELF_S1 else 2*d1
        if dim==3:
            u,anchor=_triplet_unitary(root)
            metric=u.T@u
            if np.linalg.norm(metric.imag)>1e-10:
                raise ValueError('Matchete native spin-1 self-metric is not real')
            metric=metric.real
            old=m[start:start+6,start:start+6].copy()
            if not np.allclose(old,np.kron(np.eye(3),old[:2,:2]),rtol=1e-11,atol=1e-6):
                raise ValueError('Source triplet mass not in expected Cartesian self form')
            avg=float(np.trace(old[:2,:2])/2)
            k=old[:2,:2]-avg*np.eye(2)
            m[start:start+6,start:start+6]=avg*np.eye(6)+np.kron(metric,k)
            details.update({'corrected':True,'species':'S1' if model in SELF_S1 else 'S2',
                           'metric_source':'A(-4) original Matchete Cartesian triplet',
                           'source_CG_sha256':anchor['source_sha256'],
                           'spherical_symmetric_metric':metric.tolist()})
    ev0=np.linalg.eigvalsh(original)
    ev1=np.linalg.eigvalsh(m)
    if not np.allclose(ev0,ev1,atol=1e-6,rtol=1e-10):
        raise ValueError('Mass correction changed scalar eigenvalues')
    w=mass_ward_report(m,b)
    if max(w.values())>1e-9:
        raise ValueError('Canonical source mass is not gauge invariant: '+repr(w))
    details['corrected_mass_eigenvalues_sq_gev2']=ev1.tolist()
    details['corrected_mass_Ward']=w
    return m,details


def canonical_mix_tensor(d1: int, d2: int) -> np.ndarray:
    """Unique HH S1 S2* spin-singlet, using explicit invariant eps/generator CGs.

    This is the SU2-covariant *conditional component* of the RGBeta delS2
    registration, not a certified native RGBeta coefficient calibration.
    """
    g = builder._cg()
    eps = np.array([[0., 1.], [-1., 0.]])
    j = np.array([[0., 0., 1.], [0., -1., 0.], [1., 0., 0.]])
    if (d1,d2) == (1,3):
        t = g.transpose(1,2,0).reshape(2,2,1,3)
    elif (d1,d2) == (3,1):
        t = np.einsum('mij,ma->ija', g, j).reshape(2,2,3,1)
    elif (d1,d2) == (2,2):
        t = (np.einsum('ik,jl->ijkl',eps,np.eye(2))
             + np.einsum('jk,il->ijkl',eps,np.eye(2))) / 2
    elif (d1,d2) == (3,3):
        g2, g3 = spin_generators(2), spin_generators(3)
        t = sum(np.einsum('ij,ba->ijab',eps@g2[k],g3[k]) for k in range(3))
        # The original conditional RGBeta delS2 triple contraction has this
        # invariant's norm: this relative factor is only a convention choice.
        old = np.einsum('mij,mab->ijab',g,builder._triple(g),optimize=True)
        t = t * (np.linalg.norm(old) / np.linalg.norm(t))
    else:
        raise ValueError(f'Unsupported SU2 scalar pair {d1,d2}')
    err = _su2_variation(t,d1,d2,first_conjugated=False)
    if err > 1e-11:
        raise ValueError(f'Canonical lambdaT3 invariant construction failed: {err}')
    return np.asarray(t,complex)


def _holomorphic_bilinear(ns: int, offset_a: int, offset_b: int,
                          matrix: np.ndarray) -> np.ndarray:
    b = np.zeros((ns,ns),complex)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if matrix[i,j] != 0:
                b += matrix[i,j]*np.outer(
                    builder._u(ns,offset_a,i),builder._u(ns,offset_b,j))
    return (b+b.T)/2


def source_invariant_units(model: str) -> dict[str, np.ndarray]:
    """Replace only the known noncovariant conditional RGBeta components."""
    basis = build_basis(model)
    d1,d2,_ = basis.dims
    units = _BASE_QUARTIC_UNITS(model)
    if (d1,d2) != (1,3):
        t = canonical_mix_tensor(d1,d2)
        lam = np.zeros((basis.ns,)*4,complex)
        for i,j,k,l in np.argwhere(np.abs(t)>1e-14):
            lam += t[i,j,k,l] * builder._monomial(
                builder._u(basis.ns,0,int(i)),builder._u(basis.ns,0,int(j)),
                builder._u(basis.ns,4,int(k)),
                builder._u(basis.ns,4+2*d1,int(l),True))
        units['lambdaT3'] = lam
    if (d1,d2) == (3,3):
        # RGBeta source: delta[Bar S1,Bar S2]delta[S1,S2]
        # + delta[Bar S1,S2]delta[S1,Bar S2]. The same-variance
        # spin-1 contraction is J, NOT ordinary identity in spherical axes.
        j = np.array([[0.,0.,1.],[0.,-1.,0.],[1.,0.,0.]])
        a = _holomorphic_bilinear(basis.ns,4,10,j)
        b = builder._bilinear(basis.ns,4,10,np.eye(3))
        units['lambda12Cross'] = (
            builder._product(a,a.conj()) + builder._product(b,b.conj())).real
    return units


def assert_source_rgbeta_definition(root: Path) -> dict:
    """Reject a changed source invariant declaration rather than trusting formulas."""
    rel = 'RGE/running/rgbeta/T3RGBetaModel.wl'
    source = (Path(root)/rel).read_text(encoding='utf-8-sig')
    stripped = _compact(source)
    for phrase in (
        'T3RGBetaMixInvariant[1,3]:=',
        'T3RGBetaMixInvariant[3,1]:=',
        'T3RGBetaMixInvariant[2,2]:=',
        'T3RGBetaMixInvariant[3,3]:=',
        '{H,H,S1,Bar@S2}',
        'del[SU2L@S2,#1,#3]*del[SU2L@S2,#2,#4]+'
        'del[SU2L@S2,#1,#4]*del[SU2L@S2,#2,#3]',
    ):
        if phrase not in stripped:
            raise ValueError('Unrecognized native RGBeta source definition: '+phrase)
    return {'RGBeta_source_path':rel,'RGBeta_source_sha256':hash_file(Path(root)/rel),
            'component_frame':'explicit conditional highest-weight SU2; not native-normalisation-certified'}


def mass_ward_report(mass2: np.ndarray, basis) -> dict[str,float]:
    m = np.asarray(mass2,float)
    nh = basis.ns-4
    if m.shape != (nh,nh) or not np.isfinite(m).all() or not np.allclose(m,m.T,atol=1e-8):
        raise ValueError('Mass matrix must be real, symmetric, and heavy-only')
    values = {}
    for name, pairs in basis.generators().items():
        nr = []
        for _,s in pairs:
            s = s[4:,4:]
            nr.append(float(np.linalg.norm(m@s-s@m)) /
                      max(float(np.linalg.norm(m)*np.linalg.norm(s)),1e-100))
        values[name] = max(nr,default=0.)
    return values


def validated_wilson_map(basis, initial: np.ndarray, raw: np.ndarray,
                         factor: complex, r: np.ndarray, tol: float) -> tuple[np.ndarray,dict]:
    """Orthogonal source frame conversion of the complete LLSS, no projection."""
    ident = embed_formal_llss(raw,basis,np.eye(len(r)),factor)
    mapped = embed_formal_llss(raw,basis,r,factor)
    # Preserve any original hard operator, including complete HH polynomial.
    extra = np.asarray(initial,complex)-ident
    f=[i for pair in basis.leptons for i in pair]
    heav=range(4,basis.ns)
    ix=np.ix_(f,f,heav,heav)
    n=float(np.linalg.norm(ident[ix]))
    if n < 1e-100:
        raise ValueError('Empty original Matchete LLSS boundary')
    if np.linalg.norm(extra[ix]) > 1e-8*n:
        raise ValueError('Native boundary differs from source formal LLSS; cannot transform safely')
    result=mapped+extra
    delta=float(abs(np.linalg.norm(mapped[ix])-n)/n)
    if delta > 1e-10:
        raise ValueError('Source CG conversion lost LLSS tensor norm')
    old=initial.copy(); new=result.copy()
    old[ix]=0;new[ix]=0
    if not np.allclose(old,new,atol=1e-24,rtol=1e-13):
        raise ValueError('Source conversion modified non-LLSS Wilson operators')
    ward=_gauge_norm(result,basis)
    if max(ward.values()) > tol:
        raise ValueError(f'Boundary gauge Ward failed {ward}; coordinate convention unresolved')
    return result, {'norm_ratio':float(np.linalg.norm(mapped[ix])/n),
                    'norm_relative_change':delta,'ward':ward,
                    'non_LLSS_unchanged':True, 'not_a_Ward_projection':True}
