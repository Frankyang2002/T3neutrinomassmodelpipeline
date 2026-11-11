"""Rebuild same-muS A(-4) scalar HHSS background for conditional Batch17.

Project path: Numerical/experimental/Batch17RebuildBackground.py

Uses the exact Batch14 muS couplings, the Batch16 muS mass Hessian,
source-hash-checked native Matchete HHSS channels, and existing RGBeta
portal/general reference-invariant tensor machinery. The native->RGBeta
coupling map and historical F-tree flavour source are NOT independently
certified against the new illustrative UV benchmark. No physical C5 claim.

Only the ordinary, identical-source-basis A(-4) model is initially enabled.
The six self-bilinear models and B(-1)/C(-1) require separate CG-basis
certification and may not be silently treated as unchanged-field models.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping

import numpy as np

MODEL = 'A_-4'
POTENTIAL = 'V=lambda_abcd*phi_a*phi_b*phi_c*phi_d/4!'


def sha256(path: Path | str) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_json(path: Path | str) -> dict[str, Any]:
    report = json.loads(Path(path).read_text(encoding='utf-8-sig'))
    if not isinstance(report, dict):
        raise ValueError(f'{path}: expected JSON object')
    return report


def scalar_couplings_at_mu_s(endpoint: Mapping[str, Any]) -> dict[str, Any]:
    """Decode every scalar parameter without inventing a zero or alias.

    Matrix-valued Yukawas are not scalar-quartic inputs; they are ignored here
    but remain in the unmodified Batch14 numerical state.
    """
    if endpoint.get('model') != MODEL or endpoint.get('status') != 'CandidateFThresholdAndEFT1Integrated':
        raise ValueError('Expected the successful A_-4 Batch14 EFT1 endpoint')
    raw = endpoint.get('final')
    if not isinstance(raw, Mapping):
        raise ValueError('Batch14 final coupling object missing')
    values: dict[str, Any] = {}
    for name, encoded in raw.items():
        if name in ('yu','yd','ye'):
            continue
        if not isinstance(encoded,Mapping) or set(encoded) != {'real','imag'}:
            raise ValueError(f'{name}: expected explicit real/imag numeric value')
        real = np.asarray(encoded['real'],dtype=float)
        imag = np.asarray(encoded['imag'],dtype=float)
        if real.shape != () or imag.shape != () or not (np.isfinite(real).all() and np.isfinite(imag).all()):
            raise ValueError(f'{name}: scalar coupling must be finite and rank zero')
        x,y = float(real),float(imag)
        # All nonzero imaginary parts must remain explicit. Source-quartic
        # evaluators accept scalar numbers or {real,imag} dictionaries.
        values[name] = {'real':x,'imag':y} if y != 0 else x
    return values


def expected_mass_at_mu_s(mass_audit: Mapping[str, Any]) -> np.ndarray:
    """A_-4 has neither self-bilinear nor mixed-doublet mass operators."""
    if mass_audit.get('model') != MODEL or mass_audit.get('status') != 'CandidateNativeQuadraticHessianComputed':
        raise ValueError('Wrong Batch16 scalar Hessian model/status')
    if mass_audit.get('quadratic_kind') != 'none':
        raise ValueError('A_-4 unexpectedly contains quadratic scalar mixing')
    c = mass_audit.get('conventions',{}).get('rgbeta_potential_plus')
    if not isinstance(c, Mapping) or c.get('tachyon_free_candidate') is not True:
        raise ValueError('No positive candidate muS scalar mass Hessian')
    h = np.asarray(c.get('real_mass_hessian_gev2'),dtype=float)
    if h.shape != (8,8) or not np.isfinite(h).all() or not np.allclose(h,h.T,rtol=0,atol=1e-8):
        raise ValueError('Invalid A_-4 real scalar Hessian')
    return h


def check_native_mass_matrix(mass_audit: Mapping[str,Any], couplings: Mapping[str,Any],
                             *, tolerance: float=1e-10) -> np.ndarray:
    """Compare matrices, not only spectra, in the unchanged A(-4) source basis."""
    from Numerical.experimental.Batch17ConditionalWilson import mass_parity
    h=expected_mass_at_mu_s(mass_audit)
    m1=float(couplings['mS1Sq']);m2=float(couplings['mS2Sq'])
    if not (math.isfinite(m1) and math.isfinite(m2) and min(m1,m2)>0):
        raise ValueError('Invalid muS scalar running masses')
    source_matrix = np.diag([m1]*2+[m2]*6)
    scale=max(abs(m1),abs(m2),1.)
    mismatch=float(np.max(np.abs(h-source_matrix))/scale)
    if mismatch>tolerance:
        raise ValueError(f'Batch16 Hessian and actual Batch14 scalar masses disagree ({mismatch:.4g})')
    mass_parity(mass_audit,source_matrix,tolerance=tolerance)
    return source_matrix


def rebuild_a_m4_background(*, endpoint_path: Path, handoff_path: Path,
                            mass_path: Path, source_path: Path,
                            quartic_cache_path: Path, quartic_seed_path: Path,
                            registry_path: Path, destination: Path) -> dict[str, Any]:
    """Build source-compatible HHSS background from NEW EFT1 endpoint.

    Existing source LLSS and Matchete quartic data are read-only. This
    produces an explicit, provisional cross-package quartic map as in the
    existing GeneralT3SequentialLeadingLog.scalar_tensors implementation.
    """
    from Numerical.experimental.Batch17ConditionalWilson import check_saved_pipeline
    from Numerical.joint_scalar_matching.adapter import CONVENTION
    from Numerical.joint_scalar_matching.fermion_boundary_bridge import load_source
    from Numerical.joint_scalar_matching.large_batch import load_background
    from Numerical.joint_scalar_matching.numerical_benchmark_v5 import (
        RunningSnapshot, read_v4_cached_channels,
    )
    from Numerical.running.GeneralT3SequentialLeadingLog import scalar_tensors
    from Numerical.lo_certification.couplings.extra_hhss_geometry import model_extra_geometry

    endpoint=load_json(endpoint_path)
    handoff=load_json(handoff_path)
    mass_audit=load_json(mass_path)
    mu_f,mu_s,_=check_saved_pipeline(MODEL,handoff,mass_audit,endpoint)
    source=load_source(source_path,expected_model=MODEL)
    if not math.isclose(mu_f,source.mu_f_gev,rel_tol=1e-11):
        raise ValueError('Matchete source and new Batch14 fermion matching scale differ')
    couplings=scalar_couplings_at_mu_s(endpoint)
    h=check_native_mass_matrix(mass_audit,couplings)

    channels=read_v4_cached_channels(quartic_cache_path,quartic_seed_path,source)
    # No unverified original->source transformation is permitted in A(-4).
    if tuple(channels.source_basis_labels) != tuple(source.labels):
        raise ValueError('A_-4 source and Matchete HHSS scalar bases differ')
    registry=load_json(registry_path)
    entry=registry.get('models',{}).get(MODEL)
    if not isinstance(entry,dict):
        raise ValueError('Missing A_-4 general quartic registry entry')
    geometry=model_extra_geometry(MODEL,entry,include_components=True)
    snap=RunningSnapshot(MODEL,Path(endpoint_path),sha256(endpoint_path),mu_f,mu_s,
                         float(couplings['mS1Sq']),float(couplings['mS2Sq']),
                         couplings,{'source':'Batch14 actual evolved muS EFT1 endpoint'})
    tensor=scalar_tensors(MODEL,source,channels,snap,registry,geometry=geometry)
    if tensor.mass2.shape != (8,8) or not np.allclose(tensor.mass2,h,rtol=1e-10,atol=1e-6):
        raise ValueError('Generated tensor matching mass matrix disagrees with Batch16 native Hessian')
    if tensor.channels.shape[1:] != (4,4,8,8):
        raise ValueError('Incorrect HHSS channel tensor shape')
    if not np.allclose(tensor.channels.sum(axis=0),tensor.quartic,rtol=1e-10,atol=1e-12):
        raise ValueError('Individual HHSS channels fail to reproduce total tensor')

    destination=Path(destination)
    output_npz=destination.with_suffix('.npz')
    if destination.exists() or output_npz.exists():
        raise FileExistsError('Refusing to overwrite an existing matching background')
    destination.parent.mkdir(parents=True,exist_ok=True)
    names=list(tensor.names)
    metadata={
       'format':'t3-joint-scalar-conditional-background-v1',
       'model':MODEL,
       'source_boundary_sha256':source.hash,
       'source_scalar_real_labels':list(source.labels),
       'renormalization_scheme':'MSbar',
       'tensor_convention':CONVENTION,
       'mu_f_gev':mu_f,'mu_s_gev':mu_s,
       'reference_scale_gev':mu_s,
       'scalar_mass_origin':'unbroken_electroweak',
       'same_muS_quartic_and_mass':True,
       'no_separate_scalar_threshold':True,
       'tensor_npz':output_npz.name,
       'quartic_channel_names':names,
       'provenance':{
           'source':'Batch14 actual muS RG endpoint + source-hash-checked native Matchete HHSS cache and RGBeta portal/general candidate geometry',
           'quartic_potential_normalization':POTENTIAL,
           'endpoint_sha256':sha256(endpoint_path),
           'batch15_sha256':sha256(handoff_path),
           'batch16_sha256':sha256(mass_path),
           'matchete_source_sha256':source.hash,
           'matchete_scalar_seed_sha256':sha256(quartic_seed_path),
           'matchete_scalar_cache_sha256':sha256(quartic_cache_path),
           'general_registry_sha256':sha256(registry_path),
           'matchete_rgbeta_coupling_normalisation_certified':False,
           'historical_fermion_boundary_matches_new_uv_flavour_benchmark_certified':False,
           'tensor_basis_assumption':'A_-4 native identical scalar real component basis; RGBeta/Matchete quartic crosswalk still conditional',
           'tensor_provenance':tensor.provenance,
       },
       'physical_C5_certified':False,
       'full_one_loop_certified':False,
       'GammaCC_to_Weyl_normalization':None,
       'sm_eft1_wilson_master_complete':False,
    }
    # Write NPZ first; only publish JSON metadata after complete validation.
    try:
        np.savez_compressed(output_npz,
            source_to_matching_real_U=np.eye(8),
            mass_squared_gev2_muS=tensor.mass2,
            lambda_HHSS_muS=tensor.quartic,
            lambda_HHSS_channels_muS=tensor.channels)
        destination.write_text(json.dumps(metadata,indent=2,allow_nan=False)+'\n',encoding='utf-8')
        load_background(destination,source)
    except Exception:
        destination.unlink(missing_ok=True)
        output_npz.unlink(missing_ok=True)
        raise
    return {'model':MODEL,'status':'ConditionalBackgroundRebuiltFromBatch14',
            'background_json':str(destination),'background_npz':str(output_npz),
            'n_quartic_channels':len(names),
            'maximum_quartic_abs':float(np.max(abs(tensor.quartic))),
            'mass_eigenvalues_gev2':np.linalg.eigvalsh(tensor.mass2).tolist(),
            'same_muS_mass_and_quartic':True,
            'physical_C5_certified':False,
            'matchete_rgbeta_coupling_normalisation_certified':False,
            'historical_matchete_F_tree_flavour_boundary_certified_for_Batch13':False}
