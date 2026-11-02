"""Source-backed *formal* all-heavy-fermion LLSS decomposition for T3.

Only the verified diagonal-flavour numerical benchmark is reconstructible
from the published source-summed PL coefficients and the stored running state.
The source has no index labelling individual heavy fermion eigenstates; for
these benchmarks the exact p=q=k support together with diagonal Yukawa and
MF matrices permits a faithful decomposition. This does not prove the
GammaCC-to-Weyl signed map or general complex-flavour mass-eigenstate matching.

No scalar mass, quartic or normalization is inferred by this module.
"""
from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path
from typing import Any
import json
import math
import numpy as np


class SourceBridgeError(ValueError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SourceBridgeError(message)


def read_json(path: Path) -> dict:
    val = json.loads(Path(path).read_text(encoding='utf-8-sig'))
    require(isinstance(val,dict),f'Expected JSON object in {path}')
    return val


def digest(path: Path) -> str:
    return sha256(Path(path).read_bytes()).hexdigest()


def complex_running_matrix(obj: Any, label: str) -> np.ndarray:
    """Read the all16_v1 running.json complex-number cell representation."""
    require(isinstance(obj,list) and len(obj)==3 and all(isinstance(r,list) and len(r)==3 for r in obj),
            f'{label}: 3x3 matrix required')
    z=np.zeros((3,3),complex)
    for i in range(3):
        for j in range(3):
            v=obj[i][j]
            require(isinstance(v,dict) and set(v)=={'real','imag'},
                    f'{label}: complex cell must have real/imag')
            z[i,j]=complex(float(v['real']),float(v['imag']))
    require(np.isfinite(z).all(),f'{label}: nonfinite entry')
    return z


def diagonal_fermion_benchmark(running: dict, mu_f_gev: float) -> np.ndarray:
    """Return three *positive real* masses, only for diagonal y1,y2,MF basis.

    This does not derive a general Takagi map. Off-diagonal/complex-phase
    cases must use an independently verified UV vertex/mass basis instead.
    """
    require(running.get('status') is not None,'Running status missing')
    scales=running.get('scales_gev')
    require(isinstance(scales,dict),'Missing saved thresholds')
    # Running states are at the single existing fermion threshold.
    fscale=scales.get('F')
    require(fscale is not None and math.isclose(float(fscale),mu_f_gev,rel_tol=1e-11),
            'Saved running muF differs from LLSS release muF')
    uv=running.get('uv_at_fermion_threshold')
    require(isinstance(uv,dict),'Missing UV-at-muF state')
    y1=complex_running_matrix(uv.get('y1'),'y1(muF)')
    y2=complex_running_matrix(uv.get('y2'),'y2(muF)')
    mf=complex_running_matrix(uv.get('MF'),'MF(muF)')
    for label,a in [('y1',y1),('y2',y2),('MF',mf)]:
        off=a.copy();np.fill_diagonal(off,0)
        scale=max(float(np.max(abs(a))),1e-250)
        require(float(np.max(abs(off))) <= 1e-12*scale,
                f'{label}: non-diagonal flavour basis, no unique p=q=k source decomposition')
    diag=np.diag(mf)
    require(np.max(abs(diag.imag))<1e-11*max(float(np.max(abs(diag))),1.0),
            'Complex heavy-fermion mass phase needs explicit Takagi map')
    require(np.all(diag.real>0), 'Only positive real MF benchmark entries supported')
    require(np.all(np.abs(np.diag(y1))>0) and np.all(np.abs(np.diag(y2))>0),
            'Vanishing Yukawa direction prevents independently checking flavour-to-mass association')
    return diag.real.copy()


@dataclass(frozen=True)
class FormalDecomposition:
    model: str
    c_total: np.ndarray
    c_by_fermion: np.ndarray
    fermion_masses_gev: np.ndarray
    labels: tuple[str,...]
    source_indices: tuple[int,...]
    mu_f_gev: float
    source_hash: str
    running_hash: str


def split_formal_by_flavour(source: Any, running: dict, running_hash: str) -> FormalDecomposition:
    """Respect every PL source, all scalar components, and full p=q=k blocks.

    Valid for the repository's explicitly diagonal benchmarks only. No
    multiplication by the lepton/scalar orbit size or by an unknown Weyl sign.
    """
    raw=np.asarray(source.raw,dtype=complex)
    by=np.asarray(source.by_source,dtype=complex)
    n=len(source.labels)
    require(raw.shape==(3,3,2,2,n,n),'Invalid total LLSS shape')
    require(by.ndim==7 and by.shape[1:]==raw.shape and len(by)==len(source.source_indices),
            'Invalid source-resolved PL LLSS array')
    scale=max(float(np.max(abs(by))),float(np.max(abs(raw))),1e-250)
    require(np.max(abs(by.sum(axis=0)-raw)) <= 1e-9*scale+1e-28,
            'Source sum does not reproduce recorded total')
    for k in range(len(by)):
        off=by[k].copy()
        for p in range(3):
            off[p,p]=0
        require(float(np.max(abs(off))) <= 1e-12*scale,
                'Off-diagonal lepton flavour present: per-fermion attribution ambiguous')
    masses=diagonal_fermion_benchmark(running,float(source.mu_f_gev))
    ck=np.zeros((3,3,3,2,2,n,n),complex)
    for k in range(3):
        # The existing matching source and UV running convention both have
        # y1,y2,MF diagonal in the identical flavour index ordering.
        ck[k,k,k]=raw[k,k]
    require(np.max(abs(ck.sum(axis=0)-raw)) <= 1e-12*scale+1e-28,
            'Flavour-block split failed to reconstruct original LLSS')
    return FormalDecomposition(source.model,raw,ck,masses,tuple(source.labels),
        tuple(source.source_indices),float(source.mu_f_gev),source.hash,running_hash)


def prepare_conditional_input(formal: FormalDecomposition, scalar_json: dict,
                              scalar_data: dict, scalar_sha256: str) -> tuple[dict,dict]:
    """Combine *supplied* matching tensors with formal source, never invent scalar data.

    This produces **per-unknown signed normalization**, not a physical C5.
    Scalar input must explicitly declare a separate provenance and scheme;
    no B/C mixing or self-split mass treatment is guessed here.
    """
    from Numerical.lo_certification.completeness.hard_matching import evaluate_hard_scalar_channel
    m=formal.model
    require(scalar_json.get('model')==m,'Scalar-bundle model mismatch')
    require(scalar_json.get('tensor_convention')=='T3_REAL_LLSS_HHSS_WEINBERG_v1',
            'Explicit real-component master tensor convention required')
    require(scalar_json.get('scheme')=='MSbar','MSbar scheme declaration required')
    require(scalar_json.get('source_scalar_real_labels')==list(formal.labels),
            'Scalar basis labels not identical; explicit source-to-target rotation required')
    require(scalar_json.get('complete_scalar_HHSS_channels_supplied') is True,
            'Incomplete scalar quartic input; no implicit missing zero channels')
    require(isinstance(scalar_json.get('provenance'),dict) and scalar_json['provenance'],
            'External scalar tensor derivation/provenance required')
    mu_s=float(scalar_json.get('mu_s_gev',float('nan')))
    require(math.isfinite(mu_s) and 0<mu_s<=formal.mu_f_gev,'Require one positive muS <= muF')
    names=scalar_json.get('quartic_channel_names')
    require(isinstance(names,list) and names and len(names)==len(set(names)),
            'Every HHSS quartic contribution must have a unique channel name')
    require({'lambda_HHSS','lambda_HHSS_channels','scalar_mass_squared_gev2'} <= set(scalar_data),
            'Explicit mass matrix and all named HHSS channels required')
    q=np.asarray(scalar_data['lambda_HHSS'],float)
    qs=np.asarray(scalar_data['lambda_HHSS_channels'],float)
    ms=np.asarray(scalar_data['scalar_mass_squared_gev2'],float)
    n=len(formal.labels)
    require(q.shape==(4,4,n,n) and qs.ndim==5 and qs.shape[1:]==q.shape and
            qs.shape[0]==len(names) and ms.shape==(n,n),'Scalar component shapes inconsistent')
    # Call the hard kernel for exactly its own strict tensor/Weinberg validations.
    diagnostic=evaluate_hard_scalar_channel(formal.c_by_fermion,formal.fermion_masses_gev,
        q,ms,formal.mu_f_gev,mu_s,C_total_muF=formal.c_total,
        lambda_channels=qs,channel_names=tuple(names))
    require(bool(diagnostic['projection_passed']),
            'Full Weinberg tensor projection failed; numerical input is not physical C5')
    metadata={
        'model':m,'tensor_convention':'T3_REAL_LLSS_HHSS_WEINBERG_v1',
        'scheme':'MSbar','tensor_npz':m.replace('-','m').replace('+','p')+'.npz',
        'mu_f_gev':formal.mu_f_gev,'mu_s_gev':mu_s,
        'quartic_channel_names':names,'source_normalization_certified':False,
        'interpretation':'C5 per unknown signed N_GammaCC_to_Weyl, benchmark diagonal flavour ONLY',
        'physical_C5_certified':False,
        'provenance':{'fermion_source_sha256':formal.source_hash,
                      'running_sha256':formal.running_hash,
                      'scalar_bundle_sha256':scalar_sha256,
                      'scalar_bundle_provenance':scalar_json['provenance']},
    }
    arrays={'C_LLSS_by_fermion':formal.c_by_fermion,
            'C_LLSS_total_muF':formal.c_total,
            'fermion_masses_gev':formal.fermion_masses_gev,
            'lambda_HHSS':q,'lambda_HHSS_channels':qs,
            'scalar_mass_squared_gev2':ms}
    return metadata,arrays
