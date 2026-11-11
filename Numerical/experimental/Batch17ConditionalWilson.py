"""Batch 17: source-explicit candidate EFT1 LLSS -> joint scalar C5 -> SMEFT.

Intended path: Numerical/experimental/Batch17ConditionalWilson.py

This does NOT certify the GammaCC-to-Weyl factor, RGBeta/Matchete scalar
quartic coefficient map, complete EFT1 dim-five operator inventory, or hard
one-loop matching. Those missing inputs are never inferred from mass audits.

Mandatory constraints:
  - the source LLSS release and joint HHSS background identify each other
  - the background mass spectrum agrees with the *new* Batch16 muS Hessian
  - both scalar multiplets are integrated out in ONE joint bubble
  - only the already-implemented LLSS-driven leading log + finite bubble
    enters C5, per unit unknown source-to-Weyl factor N
  - the resulting candidate C5 is evolved with the project's SMEFT runner.

Any inconsistency blocks the candidate calculation; it never silently uses
legacy masses inconsistent with the new renormalisable trajectory.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping

import numpy as np

MODELS = ('A_+0','A_-2','A_-4','B_+1','B_-1','B_-3',
          'C_+1','C_-1','C_-3','D_+0','D_+2','D_-2',
          'E_+0','E_+2','E_-2','E_-4')


def load_json(path: str | Path) -> dict[str,Any]:
    obj=json.loads(Path(path).read_text(encoding='utf-8-sig'))
    if not isinstance(obj,dict):
        raise ValueError(f'JSON must be an object: {path}')
    return obj


def file_hash(path: str | Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _read_value(raw: Any, *, shape: tuple[int,...], name: str) -> np.ndarray:
    if not isinstance(raw,Mapping) or set(raw) != {'real','imag'}:
        raise ValueError(f'{name}: expected explicit real/imag encoding')
    a=np.asarray(raw['real'],dtype=float)
    b=np.asarray(raw['imag'],dtype=float)
    if a.shape!=shape or b.shape!=shape or not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError(f'{name}: expected finite shape {shape}')
    return a+1j*b


def check_saved_pipeline(model: str, handoff: Mapping[str,Any],
                         mass_audit: Mapping[str,Any],
                         eft_endpoint: Mapping[str,Any]) -> tuple[float,float,dict]:
    """Reject mismatched or falsely certified Batch14/15/16 inputs."""
    if model not in MODELS:
        raise ValueError(f'Unsupported model {model!r}')
    expected = (
        (handoff,'CandidateScalarThresholdHandoffAudited'),
        (mass_audit,'CandidateNativeQuadraticHessianComputed'),
        (eft_endpoint,'CandidateFThresholdAndEFT1Integrated'),
    )
    for payload,status in expected:
        if payload.get('model')!=model or payload.get('status')!=status:
            raise ValueError(f'{model}: invalid or crossed Batch14-16 model/status')
        if payload.get('physical_neutrino_mass_prediction') is not False:
            raise ValueError('Input unexpectedly asserts a physical neutrino mass prediction')
    if handoff.get('both_scalars_integrated_out_in_renormalisable_sm_projection') is not True:
        raise ValueError('Expected simultaneous renormalisable scalar projection')
    if eft_endpoint.get('both_scalars_active_in_eft1') is not True:
        raise ValueError('EFT1 did not retain both scalar multiplets')
    mu_f=float(handoff['mu_f_gev']);mu_s=float(handoff['mu_s_gev'])
    for scale,key in ((mu_f,'mu_f_gev'),(mu_s,'mu_s_gev')):
        if not math.isfinite(scale) or scale<=0:
            raise ValueError('Nonfinite or nonpositive threshold scale')
        if not math.isclose(scale,float(eft_endpoint[key]),rel_tol=1e-12):
            raise ValueError(f'Batch14 and Batch15 {key} disagree')
    if mu_s>=mu_f:
        raise ValueError('Invalid fermion-first threshold ordering')
    qkind=mass_audit.get('quadratic_kind')
    if qkind not in ('none','self_s1','self_s2','mixed'):
        raise ValueError('Invalid candidate mass-matrix kind')
    expected_symbols=handoff['scalar_threshold']['additional_quadratic_symbols']
    if expected_symbols != ([] if qkind=='none' else [mass_audit['quadratic_parameter']]):
        raise ValueError('Quadratic inventory mismatch between Batch15 and Batch16')
    conv=mass_audit.get('conventions',{}).get('rgbeta_potential_plus')
    if not isinstance(conv,Mapping):
        raise ValueError('Missing native RGBeta Hessian convention')
    if conv.get('tachyon_free_candidate') is not True:
        raise ValueError('Candidate scalar mass matrix is not positive')
    if conv.get('eigenmass_threshold_hierarchy_candidate') is not True:
        raise ValueError('Candidate physical-eigenmass hierarchy failed')
    if mass_audit.get('scalar_CG_basis_mapping_certified') is not False:
        raise ValueError('Unexpected certified mass-basis mapping')
    state=eft_endpoint.get('final')
    if not isinstance(state,Mapping):
        raise ValueError('Missing evolved scalar EFT1 endpoint')
    sm={}
    for name in ('gY','g2','g3','lambdaH'):
        z=complex(_read_value(state[name],shape=(),name=name))
        if z.imag!=0:
            raise ValueError(f'{name}: nonreal SM parameter')
        sm[name]=float(z.real)
    for name in ('yu','yd','ye'):
        sm[name]=_read_value(state[name],shape=(3,3),name=name)
    return mu_f,mu_s,sm


def mass_parity(mass_audit: Mapping[str,Any], background_mass2: np.ndarray,
                *, tolerance: float=1e-8) -> dict[str,Any]:
    """Compare basis-INVARIANT spectra, never identify unverified eigenvectors."""
    mass=np.asarray(background_mass2,dtype=float)
    target=np.asarray(mass_audit['conventions']['rgbeta_potential_plus']
                      ['mass_squared_eigenvalues_gev2'],dtype=float)
    if mass.ndim!=2 or mass.shape!=(target.size,target.size) or not np.isfinite(mass).all():
        raise ValueError('Background scalar mass matrix shape/finiteness mismatch')
    if not np.allclose(mass,mass.T,rtol=1e-11,atol=1e-10):
        raise ValueError('Background mass matrix must be symmetric')
    actual=np.linalg.eigvalsh(mass)
    if actual[0]<=0:
        raise ValueError('Background has a tachyonic scalar eigenmode')
    scale=max(float(np.max(np.abs(actual))),float(np.max(np.abs(target))),1.)
    worst=float(np.max(np.abs(actual-target))/scale)
    if worst>tolerance:
        raise ValueError(f'Old scalar matching background is inconsistent with Batch16 muS '
                         f'endpoint (relative mass-spectrum mismatch={worst:.4g}). '
                         'Regenerate the joint-scalar HHSS background from the NEW EFT1 endpoint.')
    return {'spectrum_relative_max_difference':worst,
            'background_mass_squared_gev2':actual.tolist(),
            'mass_basis_eigenvectors_identified':False}


def run_conditional(model: str, *, handoff_path: Path, mass_path: Path,
                    endpoint_path: Path, source_path: Path, background_path: Path,
                    mu_low_gev: float=100., vev_gev: float=246.22) -> dict[str,Any]:
    """Compute per-unit-N candidate LLSS contribution, then evolve SMEFT.

    Requires external source and quartic background in the established strict
    joint-scalar convention. Never creates those tensors from coupling names.
    """
    from Numerical.joint_scalar_matching.fermion_boundary_bridge import load_source
    from Numerical.joint_scalar_matching.large_batch import (
        load_background,evaluate_conditional,parse_matrix,
    )
    from Numerical.running.SMWeinbergEvolution import (
        SMWeinbergInitialConditions,evolve_sm_weinberg,
    )
    from physics.NeutrinoMass import numerical_neutrino_mass_matrix

    handoff=load_json(handoff_path)
    mass=load_json(mass_path)
    endpoint=load_json(endpoint_path)
    mu_f,mu_s,sm=check_saved_pipeline(model,handoff,mass,endpoint)
    if not math.isfinite(mu_low_gev) or not 0<mu_low_gev<mu_s:
        raise ValueError('Require 0 < mu_low < muS')
    source=load_source(source_path,expected_model=model)
    background=load_background(background_path,source)
    if not math.isclose(source.mu_f_gev,mu_f,rel_tol=1e-11):
        raise ValueError('F-threshold source and evolved UV endpoint scales differ')
    if not math.isclose(background.mu_s,mu_s,rel_tol=1e-11):
        raise ValueError('Joint scalar matching and evolved EFT1 scales differ')
    mass_check=mass_parity(mass,background.mass2)
    # All scalar states are matched at once in evaluate_conditional; no
    # single-multiplet intermediate EFT is introduced.
    candidate=evaluate_conditional(source,background)
    c5=parse_matrix(candidate['C5_LLSS_total_one_loop_per_N_gev_inverse'])
    if not np.allclose(c5,c5.T,rtol=1e-10,atol=1e-20):
        raise RuntimeError('Candidate Weinberg coefficient lost flavour symmetry')
    initial=SMWeinbergInitialConditions(**sm,K=c5)
    evolved=evolve_sm_weinberg(initial,mu_s,mu_low_gev)
    mnu=numerical_neutrino_mass_matrix(evolved.K,vev_gev=vev_gev)
    # Singular values of complex symmetric Majorana mnu give the Takagi
    # masses, without claiming any measured physical mass prediction.
    masses=np.sort(np.linalg.svd(mnu,compute_uv=False))*1e9
    return {
        'format':'t3-batch17-conditional-llss-smeft-v1',
        'status':'ConditionalLLSSAndSMEFTComputed_NotPhysicalC5',
        'model':model,
        'mu_f_gev':mu_f,'mu_s_gev':mu_s,'mu_low_gev':mu_low_gev,
        'source_boundary_sha256':file_hash(source_path),
        'scalar_background_json_sha256':file_hash(background_path),
        'batch15_sha256':file_hash(handoff_path),
        'batch16_sha256':file_hash(mass_path),
        'batch14_sha256':file_hash(endpoint_path),
        'mass_spectrum_validation':mass_check,
        'joint_scalar_response_per_unit_N':candidate,
        'low_energy_C5_per_unit_N_gev_inverse':{
            'real':evolved.K.real.tolist(),'imag':evolved.K.imag.tolist()},
        'candidate_neutrino_masses_per_unit_N_eV':masses.tolist(),
        'smeft_nfev':int(evolved.nfev),
        'eft1_wilson_running_kind':'LLSS_quartic_source_leading_log_only',
        'all_seven_wilson_terms_run_simultaneously':False,
        'independent_hard_F_loop_matching_included':False,
        'joint_scalar_finite_bubble_included_once':True,
        'both_scalars_integrated_out_together':True,
        'gamma_cc_weyl_normalization_certified':False,
        'matchete_rgbeta_quartic_normalization_certified':False,
        'full_dim5_inventory_certified':False,
        'physical_neutrino_mass_prediction':False,
        'interpretation':(
            'Numerical masses are per unit unknown signed GammaCC-to-Weyl N. '
            'This is only the conditional LLSS-induced leading logarithm plus '
            'one finite joint scalar bubble, followed by SMEFT running. '
            'Background/source CG and quartic normalisations remain uncertified.'),
    }
