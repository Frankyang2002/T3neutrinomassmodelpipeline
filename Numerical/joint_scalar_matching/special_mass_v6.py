"""Exploratory muS mass-basis extension for eight T3 mixed/self-bilinear models.

Uses pre-existing experimental scalar mass convention code and the actual
1-TeV illustrative running snapshots. The original Matchete real scalar
components are transformed into the *named* special muF boundary components;
identifying the native component CG basis and the muF/muS scalar mass
parameters remains APPROXIMATE; the SOURCE Matchete Minus/Plus field map
is retained unchanged from muF to muS (without an extraneous spherical W).

UV -> all F out at muF -> EFT1(S1,S2) -> all scalars jointly at ONE muS.
"""
from __future__ import annotations

from pathlib import Path
import json
import math
from typing import Any

import numpy as np

from .fermion_boundary_bridge import SourceBoundary, digest, require, slug
from .source_quartic_v4 import QuarticChannels, POTENTIAL_CONVENTION, evaluate_channels
from .numerical_benchmark_v5 import RunningSnapshot, SPECIAL, select_candidate_couplings
from .large_batch import load_background
from .adapter import CONVENTION
from Numerical.experimental.T3SelfBilinearMassBasis import SelfBilinearMassBasis
from Numerical.experimental.GeneralT3MixedScalarMassBasis import MixedDoubletMassBasis

FORMAT = 't3-joint-scalar-exploratory-special-mass-background-v6'
MIXED = frozenset(('B_-1', 'C_-1'))
SELF_S1 = frozenset(('A_+0', 'D_+0', 'E_+0'))
SELF_S2 = frozenset(('A_-2', 'D_-2', 'E_-2'))
MU_PARAMETER = {'self_s1': 'muT3GeneralS1Sq',
                'self_s2': 'muT3GeneralS2Sq',
                'mixed': 'muT3GeneralS1S2'}
DIMS = {'A':(1,3),'B':(2,2),'C':(2,2),'D':(3,1),'E':(3,3)}


def kind_for_model(model: str) -> str:
    require(model in SPECIAL, f'{model} is not one of eight special scalar models')
    if model in MIXED: return 'mixed'
    if model in SELF_S1: return 'self_s1'
    if model in SELF_S2: return 'self_s2'
    raise ValueError(f'missing special scalar category: {model}')


def _finite_complex(value: Any, label: str) -> complex:
    if isinstance(value, dict) and set(value)=={'real','imag'}:
        z=complex(float(value['real']),float(value['imag']))
    elif type(value) in (float,int):
        z=complex(float(value))
    else:
        raise ValueError(f'{label}: expected explicit real/imag numerical value')
    require(math.isfinite(z.real) and math.isfinite(z.imag),f'{label}: nonfinite mass term')
    return z


def expected_special_labels(model: str) -> tuple[str,...]:
    d1,d2=DIMS[model[0]]
    kind=kind_for_model(model)
    if kind=='mixed':
        return ('NewScalarX',)*(2*d1)+('NewScalarY',)*(2*d2)
    if kind=='self_s1':
        return ('T3Scalar1Minus',)*d1+('T3Scalar1Plus',)*d1+('NewScalar2',)*(2*d2)
    return ('NewScalar1',)*(2*d1)+('T3Scalar2Minus',)*d2+('T3Scalar2Plus',)*d2


def _stack_components(z: np.ndarray) -> np.ndarray:
    a=np.asarray(z,dtype=complex)
    return np.column_stack((a.real,a.imag)).reshape(-1)*math.sqrt(2)


def _split_complex_real(flat: np.ndarray) -> np.ndarray:
    a=np.asarray(flat,float)
    require(a.size%2==0,'real scalar coordinates must have even dimension')
    return (a[0::2]+1j*a[1::2])/math.sqrt(2)


def _self_rotation(dim: int) -> np.ndarray:
    """Exact Matchete split-field map, with unchanged native SU(2) indices.

    RunT3GeneralSelfMassBasis.wl defines, component by component,
        S[k] = (T3ScalarMinus[k] + i*T3ScalarPlus[k])/sqrt(2).
    Its 'Minus' field is the Cartesian real part, 'Plus' the imaginary
    part.  In canonical real coordinates (Re, Im) interleaved in each
    triplet slot, this is a pure permutation; NO spherical intertwiner W.

    The unrelated GeneralT3ScalarQuadraticMasses / T3SelfBilinearMassBasis
    prototype uses spherical triplet indices and is not the field map
    of the source Matchete LLSS release.  Multiplying this map by W
    breaks the source Weinberg Ward contraction for d=3.
    """
    require(dim in (1,3),'native self-bilinear dimension must be 1 or 3')
    out=np.zeros((2*dim,2*dim),float)
    for k in range(dim):
        out[2*k,k]=1.0           # Re S[k] = Minus[k]/sqrt(2)
        out[2*k+1,dim+k]=1.0    # Im S[k] = Plus[k]/sqrt(2)
    return out


def _native_self_mass_block(dim:int,mass2:float,mu:complex)->np.ndarray:
    """Unbroken scalar quadratic Hessian in the FIXED native Minus/Plus basis.

    The Matchete source convention is L_int=+mu (S_dag)^2+h.c., hence
    V2 = m^2 S_dag S - [mu (S_dag)^2+h.c.].  With S=(Minus+i Plus)/sqrt2,
    the Hessian is [[(m²-2 Re mu) I,-2 Im mu I],
                    [-2 Im mu I,(m²+2 Re mu) I]].
    For nonreal mu, leave the offdiagonal mass mixing present: the
    existing joint scalar bubble diagonalizes it automatically.  This
    keeps the muF matched scalar labels unchanged at the muS endpoint.
    """
    require(dim in (1,3),'invalid real self-bilinear multiplet dimension')
    a=float(mass2);z=complex(mu)
    require(np.isfinite([a,z.real,z.imag]).all(),'nonfinite self quadratic parameters')
    I=np.eye(dim,dtype=float)
    return np.block([[(a-2*z.real)*I,-2*z.imag*I],
                     [-2*z.imag*I,(a+2*z.real)*I]])


def _mixed_rotation(U: np.ndarray) -> np.ndarray:
    """Old interleaved (S1,S2) = R @ new interleaved (X,Y).

    Reuses MixedDoubletMassBasis.original_fields, including the crucial
    -epsilon conjugation mapping of S2.
    """
    from Numerical.experimental.GeneralT3MixedScalarMassBasis import EPS
    out=np.zeros((8,8),float)
    for j in range(8):
        v=np.zeros(8)
        v[j]=1.
        X=_split_complex_real(v[:4]);Y=_split_complex_real(v[4:])
        s1=U[0,0]*X+U[0,1]*Y
        s2=-EPS@(U[1,0]*X+U[1,1]*Y).conj()
        out[:,j]=np.r_[_stack_components(s1),_stack_components(s2)]
    return out


def build_special_basis(model: str, snapshot: RunningSnapshot,
                        boundary: SourceBoundary) -> dict[str, Any]:
    """Construct R: original S1/S2 real vector = R @ boundary real vector.

    Mass eigenfields are evaluated at muS. The source muF field definitions
    are assumed to agree with these at the intended exploratory accuracy.
    """
    require(model==snapshot.model==boundary.model,'special model mismatch')
    expected=expected_special_labels(model)
    require(tuple(boundary.labels)==expected,
            f'{model}: signed source-boundary scalar ordering differs from mass prototype')
    kind=kind_for_model(model)
    parameter=MU_PARAMETER[kind]
    require(parameter in snapshot.couplings,
            f'{model}: missing {parameter} at muS; zero must not be inferred')
    mu=_finite_complex(snapshot.couplings[parameter],parameter)
    d1,d2=DIMS[model[0]]
    n=2*(d1+d2)
    m1,m2=snapshot.scalar_mass1_sq,snapshot.scalar_mass2_sq
    if kind=='mixed':
        basis=MixedDoubletMassBasis.build(model,m1,m2,mu)
        require(not basis.tachyonic,'mixed scalar squared mass not positive')
        R=_mixed_rotation(basis.U)
        lo,hi=map(float,basis.eigenvalues_sq)
        mass=np.diag([lo]*4+[hi]*4)
        masses={'X_gev2':lo,'Y_gev2':hi}
        convention='(S1,epsilon*conjugate(S2))=U@(X,Y), epsilon12=+1; muS U reused as muF boundary field definition (exploratory)'
    else:
        dim=d1 if kind=='self_s1' else d2
        species='S1' if kind=='self_s1' else 'S2'
        basis=SelfBilinearMassBasis.build(species,dim,m1 if species=='S1' else m2,mu)
        require(not basis.tachyonic,'self-bilinear squared mass not positive')
        # Use the same Matchete source-field map used to generate the muF
        # LLSS release; do not change the definition of Minus/Plus at muS.
        transform=_self_rotation(dim)
        native_mass=_native_self_mass_block(dim,m1 if species=='S1' else m2,mu)
        expected=np.sort(np.repeat(basis.eigen_masses_squared,dim))
        require(np.allclose(np.linalg.eigvalsh(native_mass),expected,
                            rtol=2e-12,atol=1e-8),
                'native Cartesian self mass spectrum differs from analytic split')
        minus=float(basis.eigen_masses_squared[1])
        plus=float(basis.eigen_masses_squared[0])
        if kind=='self_s1':
            R=np.zeros((n,n),float)
            R[:2*d1,:2*d1]=transform
            R[2*d1:,2*d1:]=np.eye(2*d2)
            mass=np.zeros((n,n),float)
            mass[:2*d1,:2*d1]=native_mass
            mass[2*d1:,2*d1:]=m2*np.eye(2*d2)
        else:
            R=np.zeros((n,n),float)
            R[:2*d1,:2*d1]=np.eye(2*d1)
            R[2*d1:,2*d1:]=transform
            mass=np.zeros((n,n),float)
            mass[:2*d1,:2*d1]=m1*np.eye(2*d1)
            mass[2*d1:,2*d1:]=native_mass
        masses={'minus_gev2':minus,'plus_gev2':plus,
                'unsplit_other_gev2':m2 if kind=='self_s1' else m1}
        convention=('Native Matchete fixed component map: S[k]=(Minus[k]+i Plus[k])/sqrt(2), '
                    'same SU2 index k, no spherical W; muS quadratic Hessian includes complex mu mixing; '
                    'muF/muS running coupling map remains exploratory')
    require(R.shape==(n,n) and np.isfinite(R).all() and
            np.allclose(R.T@R,np.eye(n),rtol=1e-11,atol=1e-11),
            'special field change is not orthogonal in canonical real scalar coordinates')
    require(np.isfinite(mass).all() and
            float(np.linalg.eigvalsh(mass)[0])>0,
            'special mass eigenstate squared mass not positive')
    # Strong independent mass check: trace is preserved under real rotation.
    expected_trace=2*d1*m1+2*d2*m2
    require(abs(np.trace(mass)-expected_trace)<=1e-9*max(1.,abs(expected_trace)),
            'special mass basis changes real scalar mass trace')
    return {'R_original_from_boundary':R,'mass2_boundary':mass,'mu2':mu,
            'mass_parameter':parameter,'mass_eigenvalues':masses,
            'quadratic_kind':kind,'mass_basis_assumption':convention,
            'orthogonality_residual':float(np.max(abs(R.T@R-np.eye(n))))}


def build_special_background(ch: QuarticChannels, boundary: SourceBoundary,
                             snapshot: RunningSnapshot, output_dir: Path) -> dict[str,Any]:
    require(ch.model in SPECIAL and ch.model==boundary.model==snapshot.model,
            'special scalar builder requires matching special model/source/snapshot')
    require(ch.source_basis_labels!=boundary.labels,
            'special case unexpectedly has identical native and boundary basis')
    require(abs(snapshot.mu_f/boundary.mu_f_gev-1)<1e-10,'fermion matching scales differ')
    basis=build_special_basis(ch.model,snapshot,boundary)
    values,aliases=select_candidate_couplings(ch,snapshot)
    native,names=evaluate_channels(ch,values)
    R=basis['R_original_from_boundary']
    # quartic ORIGINAL source fields = R @ named matching-boundary fields
    transformed=np.einsum('cABij,ia,jb->cABab',native,R,R,optimize=True)
    summed=transformed.sum(axis=0)
    require(np.isfinite(summed).all() and np.isrealobj(transformed),
            'transformed special quartic must be real and finite')
    target=Path(output_dir)
    target.mkdir(parents=True,exist_ok=True)
    stem=slug(ch.model)
    meta_path=target/f'{stem}.json'; npz_path=target/f'{stem}.npz'
    require(not meta_path.exists() and not npz_path.exists(),
            'refusing to overwrite existing special model background')
    np.savez_compressed(npz_path,source_to_matching_real_U=np.eye(ch.n),
                        mass_squared_gev2_muS=basis['mass2_boundary'],
                        lambda_HHSS_muS=summed,
                        lambda_HHSS_channels_muS=transformed)
    obj={'format':'t3-joint-scalar-conditional-background-v1',
         'model':ch.model,'source_boundary_sha256':boundary.hash,
         'source_scalar_real_labels':list(boundary.labels),
         'renormalization_scheme':'MSbar','tensor_convention':CONVENTION,
         'mu_f_gev':snapshot.mu_f,'mu_s_gev':snapshot.mu_s,
         'reference_scale_gev':snapshot.mu_s,
         'scalar_mass_origin':'unbroken_electroweak',
         'same_muS_quartic_and_mass':True,'no_separate_scalar_threshold':True,
         'tensor_npz':npz_path.name,'quartic_channel_names':list(names),
         'provenance':{
           'source':f'EXPLORATORY mixed/self scalar prototype: {snapshot.source} + {ch.seed_path}',
           'quartic_potential_normalization':POTENTIAL_CONVENTION,
           'running_snapshot_sha256':snapshot.sha256,
           'source_quartic_seed_sha256':ch.seed_sha256,
           'source_boundary_sha256':boundary.hash,
           'candidate_native_to_rgbeta_aliases':aliases,
           'quadratic_parameter_name':basis['mass_parameter'],
           'mass_basis_assumption':basis['mass_basis_assumption'],
           'muF_muS_basis_alignment_certified':False,
           'native_spherical_CG_identification_certified':False,
           'coupling_alias_normalization_certified':False},
         'benchmark_not_physical':True,'physical_C5_certified':False,
         'absolute_GammaCC_to_Weyl_factor':None,
         'operator_inventory_complete':False}
    meta_path.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    bg=load_background(meta_path,boundary)
    return {'background_json':str(meta_path),'background_npz':str(npz_path),
            'mu_s_gev':bg.mu_s,'mass_eigenvalues_gev2':np.linalg.eigvalsh(bg.mass2).tolist(),
            'mass_eigenvalues_gev':np.sqrt(np.linalg.eigvalsh(bg.mass2)).tolist(),
            'quadratic_kind':basis['quadratic_kind'],
            'source_mu2_parameter':basis['mass_parameter'],
            'muS_mass_in_fixed_muF_source_basis':True,
            'source_mu2':{'real':basis['mu2'].real,'imag':basis['mu2'].imag},
            'mass_basis_assumption_NOT_CERTIFIED':basis['mass_basis_assumption'],
            'orthogonality_residual':basis['orthogonality_residual'],
            'coupling_aliases_NOT_CERTIFIED':aliases,
            'physical_C5_certified':False}
