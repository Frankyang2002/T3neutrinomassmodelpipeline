"""Batch17: provisional same-endpoint joint-scalar backgrounds for all 16 T3 models.

Intended project path: Numerical/experimental/Batch17All16Background.py

UV -> integrate out ALL F -> EFT1(S1,S2) -> integrate out BOTH scalars
simultaneously -> SMEFT. Uses exact Batch14 endpoint scalar couplings,
Batch16 *candidate* physical eigenmass spectra, and the existing source-backed
Matchete quartic channels / RGBeta portal + extra invariant machinery.

The source-field transforms for six self-bilinears and B_-1/C_-1 come from
special_mass_v6, which documents an *exploratory* CG/component identification.
Its mass eigenvectors, native extra quartic crosswalk, and historical F-tree
flavour boundary are NOT certified. A passed mass-spectrum check is NOT a
certification of the scalar field map. Physical C5 is never claimed.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping

import numpy as np

POTENTIAL = 'V=lambda_abcd*phi_a*phi_b*phi_c*phi_d/4!'


def sha256(path: Path | str) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_json(path: Path | str) -> dict[str, Any]:
    obj=json.loads(Path(path).read_text(encoding='utf-8-sig'))
    if not isinstance(obj,dict):
        raise ValueError(f'{path}: expected JSON object')
    return obj


def scalar_couplings_at_mu_s(endpoint: Mapping[str,Any], model:str) -> dict[str,Any]:
    """Read Batch14 scalar/gauge state. No missing quartics are set to zero."""
    if endpoint.get('model')!=model or endpoint.get('status')!='CandidateFThresholdAndEFT1Integrated':
        raise ValueError(f'{model}: incorrect Batch14 model/status')
    raw=endpoint.get('final')
    if not isinstance(raw,Mapping):
        raise ValueError('Batch14 final couplings missing')
    out={}
    for name,value in raw.items():
        if name in ('yu','yd','ye'):
            continue
        if not isinstance(value,Mapping) or set(value)!={'real','imag'}:
            raise ValueError(f'{model}/{name}: expected explicit real/imag encoding')
        re=np.asarray(value['real'],float);im=np.asarray(value['imag'],float)
        if re.shape!=() or im.shape!=() or not np.isfinite([float(re),float(im)]).all():
            raise ValueError(f'{model}/{name}: scalar field value nonfinite or not rank-zero')
        z=complex(float(re),float(im))
        # Preserve explicitly complex coupling data, do not discard phases.
        out[name]={'real':z.real,'imag':z.imag} if z.imag!=0 else float(z.real)
    for name in ('mS1Sq','mS2Sq'):
        if name not in out:
            raise ValueError(f'{model}: missing {name}')
        if not isinstance(out[name],float) or out[name]<=0:
            raise ValueError(f'{model}: {name} must be real and positive')
    return out


def verify_quadratic_input(mass_audit:Mapping[str,Any], state:Mapping[str,Any],model:str)->None:
    """Require Batch16 quadratic parameter to equal the Batch14 one exactly."""
    kind=mass_audit['quadratic_kind']
    key=mass_audit.get('quadratic_parameter')
    if kind=='none':
        if key is not None:
            raise ValueError(f'{model}: spurious quadratic parameter')
        return
    if key not in state:
        raise ValueError(f'{model}: absent quadratic parameter {key} at muS')
    raw=state[key]
    z=complex(float(raw['real']),float(raw['imag'])) if isinstance(raw,Mapping) else complex(raw)
    reference=mass_audit['quadratic_value']
    target=complex(float(reference['real']),float(reference['imag']))
    scale=max(abs(z),abs(target),1.)
    if abs(z-target)>1e-11*scale:
        raise ValueError(f'{model}: Batch16 quadratic term and Batch14 muS value disagree')


def certify_candidate_mass_parity(mass_audit:Mapping[str,Any], tensor_mass:np.ndarray,
                                  *,tolerance:float=1e-8)->dict:
    """Eigenvalues only for special models; no unverified basis equivalence."""
    mass=np.asarray(tensor_mass,float)
    target=np.asarray(mass_audit['conventions']['rgbeta_potential_plus']
                      ['mass_squared_eigenvalues_gev2'],float)
    if mass.shape!=(target.size,target.size) or not np.isfinite(mass).all():
        raise ValueError('Matching background scalar mass matrix dimension/finiteness mismatch')
    if not np.allclose(mass,mass.T,rtol=1e-11,atol=1e-10):
        raise ValueError('Matching background real scalar mass matrix is not symmetric')
    actual=np.linalg.eigvalsh(mass)
    if actual[0]<=0:
        raise ValueError('Nonpositive scalar mass eigenstate')
    scale=max(float(np.max(np.abs(actual))),float(np.max(np.abs(target))),1.)
    difference=float(np.max(np.abs(actual-target))/scale)
    if difference>tolerance:
        raise ValueError(f'Old scalar matching background is inconsistent with Batch16 muS endpoint '
                         f'(relative mass-spectrum mismatch={difference:.4g})')
    return {'spectrum_relative_max_difference':difference,
            'mass_basis_eigenvectors_identified':False}


def rebuild_background(*,model:str,endpoint_path:Path,handoff_path:Path,
                       mass_path:Path,source_path:Path,quartic_cache_path:Path,
                       quartic_seed_path:Path,registry_path:Path,
                       destination:Path)->dict[str,Any]:
    from Numerical.experimental.Batch17ConditionalWilson import check_saved_pipeline
    from Numerical.joint_scalar_matching.adapter import CONVENTION
    from Numerical.joint_scalar_matching.fermion_boundary_bridge import load_source
    from Numerical.joint_scalar_matching.large_batch import load_background
    from Numerical.joint_scalar_matching.numerical_benchmark_v5 import RunningSnapshot, read_v4_cached_channels
    from Numerical.running.GeneralT3SequentialLeadingLog import scalar_tensors
    from Numerical.lo_certification.couplings.extra_hhss_geometry import model_extra_geometry

    endpoint=load_json(endpoint_path);handoff=load_json(handoff_path)
    mass_audit=load_json(mass_path)
    mu_f,mu_s,_=check_saved_pipeline(model,handoff,mass_audit,endpoint)
    source=load_source(source_path,expected_model=model)
    if not math.isclose(mu_f,source.mu_f_gev,rel_tol=1e-11):
        raise ValueError(f'{model}: fermion matching scale incompatible with source tensor')
    state=scalar_couplings_at_mu_s(endpoint,model)
    verify_quadratic_input(mass_audit,state,model)
    channels=read_v4_cached_channels(quartic_cache_path,quartic_seed_path,source)
    registry=load_json(registry_path)
    entry=registry.get('models',{}).get(model)
    if not isinstance(entry,dict):
        raise ValueError(f'{model}: missing generated general quartic registry entry')
    geometry=model_extra_geometry(model,entry,include_components=True)
    snapshot=RunningSnapshot(model,Path(endpoint_path),sha256(endpoint_path),
                             mu_f,mu_s,float(state['mS1Sq']),float(state['mS2Sq']),
                             state,{'source':'new Batch14 muS EFT1 endpoint'})
    tensor=scalar_tensors(model,source,channels,snapshot,registry,geometry=geometry)
    parity=certify_candidate_mass_parity(mass_audit,tensor.mass2)
    n=len(source.labels)
    if tensor.mass2.shape!=(n,n) or tensor.channels.shape[1:]!=(4,4,n,n):
        raise ValueError(f'{model}: wrong native tensor/mass dimension')
    if not np.allclose(np.sum(tensor.channels,axis=0),tensor.quartic,atol=1e-12,rtol=1e-10):
        raise ValueError(f'{model}: HHSS channel tensor sum failed')
    special=mass_audit['quadratic_kind']!='none'
    if not special:
        # Ordinary models have identical original and source real scalar basis.
        if tuple(channels.source_basis_labels)!=tuple(source.labels):
            raise ValueError(f'{model}: ordinary scalar original/source basis differs')
        d1,d2={'A':(1,3),'B':(2,2),'C':(2,2),'D':(3,1),'E':(3,3)}[model[0]]
        expected=np.diag([state['mS1Sq']]*(2*d1)+[state['mS2Sq']]*(2*d2))
        if not np.allclose(tensor.mass2,expected,rtol=1e-11,atol=1e-6):
            raise ValueError(f'{model}: ordinary mass matrix differs from Batch14 diagonal masses')
    else:
        if tuple(channels.source_basis_labels)==tuple(source.labels):
            raise ValueError(f'{model}: expected special scalar source/boundary basis distinction')
    destination=Path(destination);npz=destination.with_suffix('.npz')
    if destination.exists() or npz.exists():
        raise FileExistsError(f'Existing background would be overwritten: {destination}')
    destination.parent.mkdir(parents=True,exist_ok=True)
    meta={
       'format':'t3-joint-scalar-conditional-background-v1','model':model,
       'source_boundary_sha256':source.hash,
       'source_scalar_real_labels':list(source.labels),
       'renormalization_scheme':'MSbar','tensor_convention':CONVENTION,
       'mu_f_gev':mu_f,'mu_s_gev':mu_s,'reference_scale_gev':mu_s,
       'scalar_mass_origin':'unbroken_electroweak',
       'same_muS_quartic_and_mass':True,'no_separate_scalar_threshold':True,
       'tensor_npz':npz.name,'quartic_channel_names':list(tensor.names),
       'provenance':{
           'source':'Batch14 evolving muS couplings; source-hash-checked Matchete HHSS + RGBeta portal/general geometry',
           'quartic_potential_normalization':POTENTIAL,
           'batch14_sha256':sha256(endpoint_path),
           'batch15_sha256':sha256(handoff_path),
           'batch16_sha256':sha256(mass_path),
           'matchete_source_sha256':source.hash,
           'matchete_scalar_seed_sha256':sha256(quartic_seed_path),
           'matchete_quartic_cache_sha256':sha256(quartic_cache_path),
           'general_registry_sha256':sha256(registry_path),
           'candidate_special_scalar_basis':special,
           'scalar_field_map_certified':False,
           'rgbeta_matchete_coupling_normalisation_certified':False,
           'source_fermion_flavour_boundary_rederived_from_Batch13':False,
           'quartic_tensor_provenance':tensor.provenance,
       },
       'physical_C5_certified':False,'full_one_loop_certified':False,
       'GammaCC_to_Weyl_normalization':None,
       'all_seven_wilson_terms_simultaneously_run':False,
    }
    try:
        np.savez_compressed(npz,
            source_to_matching_real_U=np.eye(n),
            mass_squared_gev2_muS=tensor.mass2,
            lambda_HHSS_muS=tensor.quartic,
            lambda_HHSS_channels_muS=tensor.channels)
        destination.write_text(json.dumps(meta,indent=2,allow_nan=False)+'\n',encoding='utf-8')
        load_background(destination,source)
    except Exception:
        destination.unlink(missing_ok=True);npz.unlink(missing_ok=True)
        raise
    return {
        'model':model,'status':'ConditionalBackgroundRebuiltFromBatch14',
        'background_json':str(destination),'background_npz':str(npz),
        'n_quartic_channels':len(tensor.names),
        'maximum_quartic_abs':float(np.max(np.abs(tensor.quartic))),
        'mass_eigenvalues_gev2':np.linalg.eigvalsh(tensor.mass2).tolist(),
        'batch16_spectrum_relative_difference':parity['spectrum_relative_max_difference'],
        'special_mass_basis_candidate':special,
        'physical_C5_certified':False,
        'matchete_rgbeta_coupling_normalisation_certified':False,
        'historical_fermion_boundary_rederived_from_Batch13':False,
    }
