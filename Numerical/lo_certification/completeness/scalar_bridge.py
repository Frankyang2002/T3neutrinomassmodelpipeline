"""Read-only integration of the existing simultaneous-scalar exploratory endpoint.

Only outputs inside this workstream are produced by the CLI.  The old v4
Matchete source-quartic channels and the v5/v6 numerical scalar prototypes are
reused without asserting that their native/RGBeta aliases are physical.

ScalarBackground uses s_SOURCE = U @ s_MATCHING.  The fermion-tree release is
in the SOURCE coordinates.  Pull mass and HHSS tensors back to that basis:
    M_SOURCE = U M_MATCHING U.T,
    lambda_SOURCE[ABij] = U[ia] U[jb] lambda_MATCHING[ABab].
No extra rotation of the LLSS source is made.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping
import math
import numpy as np

from .source_bridge import SourceBridgeError, require, digest

BACKGROUND_FORMAT='t3-joint-scalar-conditional-background-v1'
CONVENTION='T3_REAL_LLSS_HHSS_WEINBERG_v1'


@dataclass(frozen=True)
class ExtractedScalar:
    model: str
    mu_f_gev: float
    mu_s_gev: float
    labels: tuple[str,...]
    lambda_channels_source: np.ndarray
    lambda_source: np.ndarray
    mass2_source: np.ndarray
    channel_names: tuple[str,...]
    metadata_sha256: str
    arrays_sha256: str
    provenance: dict


def parse_scalar_background(meta: Mapping[str,Any], arrays: Mapping[str,np.ndarray], *,
                            model: str, source_sha256: str, labels: tuple[str,...],
                            mu_f_gev: float, meta_hash: str, arrays_hash: str,
                            forbid_physical_claims: bool=True) -> ExtractedScalar:
    """Validate a complete *saved native HHSS source-channel set*, not UV closure.

    All entries are source-derived in the original or explicitly rotated scalar
    basis.  The optional exploratory aliases and mass-basis identifications
    are not thereby certified; no signed GammaCC factor is inserted.
    """
    require(meta.get('format')==BACKGROUND_FORMAT, 'Scalar background format changed')
    require(meta.get('model')==model and meta.get('source_boundary_sha256')==source_sha256,
            'Scalar endpoint source hash/model differs from fermion boundary')
    require(meta.get('source_scalar_real_labels')==list(labels),
            'Scalar basis labels do not match fermion LLSS source')
    require(meta.get('renormalization_scheme')=='MSbar' and meta.get('tensor_convention')==CONVENTION,
            'Scalar scheme/tensor convention mismatch')
    require(meta.get('no_separate_scalar_threshold') is True and
            meta.get('same_muS_quartic_and_mass') is True,
            'Both scalar multiplets must remain until ONE muS with one mass/quartic scale')
    require(meta.get('scalar_mass_origin')=='unbroken_electroweak',
            'Higher-dimensional EWSB-resummed mass input not supported in dim5 hard bridge')
    if forbid_physical_claims:
        require(meta.get('physical_C5_certified') is False and
                meta.get('benchmark_not_physical') is True,
                'Expected explicitly provisional numerical scalar background')
        require(meta.get('absolute_GammaCC_to_Weyl_factor') is None and
                meta.get('operator_inventory_complete') is False,
                'Provisional endpoint must not assert signed or complete physics')
    provenance=meta.get('provenance')
    require(isinstance(provenance,dict) and provenance.get('quartic_potential_normalization')==
            'V=lambda_abcd*phi_a*phi_b*phi_c*phi_d/4!' and provenance.get('source'),
            'Explicit Matchete potential convention/source evidence missing')
    mf,ms=float(meta.get('mu_f_gev',float('nan'))),float(meta.get('mu_s_gev',float('nan')))
    ref=float(meta.get('reference_scale_gev',float('nan')))
    require(all(map(math.isfinite,(mf,ms,ref))) and 0<ms<=mf and
            math.isclose(mf,float(mu_f_gev),rel_tol=1e-11) and
            math.isclose(ref,ms,rel_tol=1e-11),
            'Scalar endpoint scales fail muF/one-muS check')
    names=tuple(meta.get('quartic_channel_names',()))
    require(len(names)>0 and len(names)==len(set(names)) and
            all(isinstance(k,str) and k.strip() for k in names),
            'Complete named quartic-source channels required')
    needed={'source_to_matching_real_U','mass_squared_gev2_muS',
            'lambda_HHSS_muS','lambda_HHSS_channels_muS'}
    require(needed.issubset(arrays),'Missing scalar NPZ entries: '+str(sorted(needed-set(arrays))))
    U=np.asarray(arrays['source_to_matching_real_U'])
    M=np.asarray(arrays['mass_squared_gev2_muS'])
    Q=np.asarray(arrays['lambda_HHSS_muS'])
    channels=np.asarray(arrays['lambda_HHSS_channels_muS'])
    n=len(labels)
    require(n in (8,12), 'Expected both scalar multiplets: n=8 or 12 real components')
    require(U.shape==(n,n) and M.shape==(n,n) and Q.shape==(4,4,n,n) and
            channels.shape==(len(names),4,4,n,n), 'Incorrect complete scalar tensor shapes')
    require(all(np.isrealobj(v) and np.isfinite(v).all() for v in (U,M,Q,channels)),
            'Real scalar tensors must be finite and non-complex')
    U=np.asarray(U,float);M=np.asarray(M,float);Q=np.asarray(Q,float)
    channels=np.asarray(channels,float)
    require(np.allclose(U.T@U,np.eye(n),atol=1e-10,rtol=1e-10),
            'Invalid SOURCE=U@MATCHING orthogonal scalar transform')
    require(np.allclose(M,M.T,rtol=1e-11,atol=1e-8) and
            np.linalg.eigvalsh(M)[0]>0,'Mass2 not positive real symmetric')
    require(np.allclose(Q,Q.swapaxes(0,1),rtol=1e-11,atol=1e-14) and
            np.allclose(Q,Q.swapaxes(2,3),rtol=1e-11,atol=1e-14),
            'Quartic HHSS tensor not symmetrized')
    require(np.allclose(channels,channels.swapaxes(1,2),rtol=1e-11,atol=1e-14) and
            np.allclose(channels,channels.swapaxes(3,4),rtol=1e-11,atol=1e-14),
            'Independent HHSS channels not real symmetric tensors')
    scale=max(float(np.max(abs(Q))),1e-250)
    require(float(np.max(abs(channels.sum(axis=0)-Q)))<=1e-10*scale+1e-28,
            'All source quartic channels must reproduce total HHSS')
    # Rotate both mass AND quartics back into the fermion-release source real basis.
    ms=U@M@U.T
    qs=np.einsum('cABab,ia,jb->cABij',channels,U,U,optimize=True)
    q=np.einsum('ABab,ia,jb->ABij',Q,U,U,optimize=True)
    require(np.allclose(qs.sum(axis=0),q,rtol=1e-10,atol=1e-12*scale+1e-28),
            'Quartic source rotation or channel sum failed')
    require(np.allclose(U.T@ms@U,M,rtol=1e-10,atol=1e-7),
            'Scalar mass basis pullback inconsistent')
    return ExtractedScalar(model,mf,float(meta['mu_s_gev']),labels,qs,q,ms,names,
                           meta_hash,arrays_hash,dict(provenance))


def read_scalar_background(path: Path, source: Any) -> ExtractedScalar:
    """Validate against the pre-existing independent joint matcher first."""
    from Numerical.joint_scalar_matching.large_batch import load_background
    from Numerical.joint_scalar_matching.fermion_boundary_bridge import read_json
    path=Path(path)
    # Existing source-aware parser performs its own channel/orbit checks.
    load_background(path,source)
    raw=read_json(path)
    name=raw.get('tensor_npz')
    require(isinstance(name,str) and name and Path(name).name==name,
            'Scalar NPZ must be adjacent to source metadata')
    npz=path.parent/name
    with np.load(npz,allow_pickle=False) as z:
        arrays={k:np.array(z[k]) for k in z.files}
    return parse_scalar_background(raw,arrays,model=source.model,
         source_sha256=source.hash,labels=tuple(source.labels),
         mu_f_gev=source.mu_f_gev,meta_hash=digest(path),arrays_hash=digest(npz))


def conditional_scalar_payload(data: ExtractedScalar) -> tuple[dict,dict]:
    """Translate to existing Workstream3 strict conditional input contract."""
    payload={
        'model':data.model,'tensor_convention':CONVENTION,'scheme':'MSbar',
        'source_scalar_real_labels':list(data.labels),
        # Means: ALL channels present in native v4 HHSS source file, not an
        # independent proof that every UV quartic is implemented.
        'complete_scalar_HHSS_channels_supplied':True,
        'source_channel_coverage_scope':'complete_cached_native_HHSS_only_NOT_general_UV_closure',
        'mu_s_gev':data.mu_s_gev,
        'quartic_channel_names':list(data.channel_names),
        'physical_C5_certified':False,
        'source_normalization_certified':False,
        'provenance':{
            'source':'v4 Matchete-source HHSS channels plus saved v5/v6 muS scalar benchmark',
            'background_metadata_sha256':data.metadata_sha256,
            'background_npz_sha256':data.arrays_sha256,
            'upstream':data.provenance,
            'candidate_coupling_aliases_certified':False,
            'special_scalar_basis_alignment_certified':False,
            'GammaCC_to_Weyl_signed_factor':None,
            'complete_UV_operator_closure_certified':False,
        },
    }
    arrays={'lambda_HHSS':data.lambda_source,
            'lambda_HHSS_channels':data.lambda_channels_source,
            'scalar_mass_squared_gev2':data.mass2_source}
    return payload,arrays


def build_existing_exploratory_scalar_background(source:Any, running_path:Path,
            v4_cached_json:Path,seed_path:Path,output_dir:Path, *,
            allow_exploratory_map:bool,
            allow_exploratory_special_bases:bool) -> Path:
    """Reuse v4/v5/v6 source routines, writing ONLY at provided owned output.

    The opt-in declares use of unverified native->RGBeta aliases and in
    eight special models unverified muF/muS scalar-basis identifications.
    """
    from Numerical.joint_scalar_matching.numerical_benchmark_v5 import (
        read_v4_cached_channels,load_running_snapshot,build_ordinary_background,
        SPECIAL,REGULAR)
    from Numerical.joint_scalar_matching.special_mass_v6 import build_special_background
    require(allow_exploratory_map,
            'Opt in to UNVERIFIED native Matchete/RGBeta coupling aliases')
    if source.model in SPECIAL:
        require(allow_exploratory_special_bases,
                'Special model requires explicit UNVERIFIED scalar basis opt-in')
    require(source.model in (SPECIAL|REGULAR),'Unknown scalar model class')
    ch=read_v4_cached_channels(Path(v4_cached_json),Path(seed_path),source)
    snap=load_running_snapshot(Path(running_path),source.model)
    output_dir=Path(output_dir)
    if source.model in SPECIAL:
        record=build_special_background(ch,source,snap,output_dir)
    else:
        record=build_ordinary_background(ch,source,snap,output_dir)
    return Path(record['background_json'])
