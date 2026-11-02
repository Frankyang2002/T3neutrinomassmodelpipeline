"""Workstream 2: all-16 baseline portal completion + exact extra-HHSS reference geometry.

Read only from the existing source release, running outputs, native Matchete
quartic cache and local RGBeta registry. No RGE or scalar matcher is copied.
No output is a certified physical C5 input: 24 extra HHSS sectors remain
unmapped to named couplings and their muF/muS beta trajectories are not supplied.
"""
from __future__ import annotations

import json
from pathlib import Path
import numpy as np

from .batch import (BASELINE_HHSS_REGISTRY, check_source_versions,
                    local_source_recheck, model_paths, read_json, digest,
                    write_json, _registry_hhss_names, observed_running_states)
from .portal_completion import completed_portal_numeric
from .extra_hhss_geometry import model_extra_geometry
from .quartic_audit import rgbeta_potential_declarations, _validate_generated_extra

FORMAT='t3-ws2-complete-baseline-portal-and-extra-reference-batch-v3'


def model_completion(root:Path,model:str,output:Path,release:Path,
                     cache:Path,registry:dict,source_decls:dict)->dict:
    from Numerical.joint_scalar_matching.fermion_boundary_bridge import load_source
    from Numerical.joint_scalar_matching.numerical_benchmark_v5 import read_v4_cached_channels
    paths=model_paths(root,model,release,cache)
    boundary=load_source(paths['boundary'],expected_model=model)
    ch=read_v4_cached_channels(paths['quartic_cache'],paths['quartic_seed'],boundary)
    running=observed_running_states(paths['running'],model)
    mf=float(running['scales_gev']['F']);ms=float(running['scales_gev']['S'])
    if not (mf==boundary.mu_f_gev and 0<ms<mf):
        raise ValueError(f'{model}: F/S grouped thresholds or boundary differ')
    states={'muF':running['eft1_at_fermion_threshold'],
            'muS':running['eft1_at_scalar_threshold']}
    portal,arrays=completed_portal_numeric(ch,states)
    extra=model_extra_geometry(model,registry,include_components=True)
    generated=_validate_generated_extra(root,model,registry)
    expected=len(_registry_hhss_names(registry))
    if extra['registry_HHSS_extra_symbols']!=expected:
        raise ValueError(f'{model}: full extra HHSS source registry inventory lost entries')
    required_portal_count=sum(a['rgbeta_direction_count'] for a in portal['species'])
    both_scale_portals_complete=all(
        len(arrays.get(f'{scale}_channel_labels',[]))==required_portal_count
        for scale in ('muF','muS'))
    # A missing running value blocks its numeric channel, but must NOT hide
    # independent symbolic transformations already completed for this model.
    expected_extra_symbols=sorted(
        {item['rgbeta_symbol'] for sector in extra['monomial_sectors']
         for item in sector['source_registry_symbols']})
    extra_availability=[{
        'rgbeta_symbol':symbol,
        'muF_present':symbol in states['muF'],
        'muS_present':symbol in states['muS'],
        'named_real_tensor_normalization_verified':False,
        'numerical_value_for_physical_eft1':None,
    } for symbol in expected_extra_symbols]
    # Cross-source provenance for the generator parser and actual numerical run.
    source_hashes={label:digest(path) for label,path in paths.items()}
    source_hashes.update({
      'RGBeta_source_wl':source_decls['sha256'],
      'Matchete_scalar_potential_wl':source_decls['matchete_portal_source_sha256'],
      'generated_extra_quartics_wl':generated['sha256']})
    # All outputs are new and exclusively under this workstream's run directory.
    output.mkdir(parents=True,exist_ok=False)
    np.savez_compressed(output/'baseline_native_plus_rgbeta_portals_muF_muS.npz',**arrays)
    write_json(output/'exact_extra_HHSS_reference_tensors.json',extra)
    blockers=[
      'additional_extra_HHSS_named_RGBeta_to_Matchete_CG_normalizations_not_certified',
      'additional_HHSS_muF_muS_numerical_coupling_trajectories_not_supplied',
      'lambdaT3_full_component_phase_and_normalization_not_certified',
      'GammaCC_to_master_Weyl_signed_normalization_not_certified',
      'complete_SM_Weyl_yukawa_gauge_embedding_not_certified',
      'special_scalar_original_to_muF_and_muS_physical_map_not_certified',
      'hard_fermion_one_loop_and_operator_closure_not_certified',
    ]
    if not both_scale_portals_complete:
        blockers.append('missing_complete_muF_or_muS_running_baseline_portal_values')
    if any(not (x['muF_present'] and x['muS_present']) for x in extra_availability):
        blockers.append('missing_extra_HHSS_running_couplings_at_muF_or_muS')
    if registry.get('quadratics_not_implemented',0):
        blockers.append('RG_closed_additional_quadratic_couplings_missing')
    report={
      'format':FORMAT,'model':model,
      'status':('exact_baseline_portal_completed_extra_reference_basis_not_physically_mapped'
          if both_scale_portals_complete else
          'exact_baseline_portal_geometry_completed_numeric_couplings_blocked'),
      'physics_status':'blocked',
      'scales_gev':{'mu_F':mf,'mu_S':ms},'scheme':'MSbar',
      'thresholds':{'mu_F':'all heavy BSM fermions together',
                    'mu_S':'both heavy scalar multiplets together',
                    'SM_fermions':'retained'},
      'baseline_portal_completion':portal,
      'baseline_portal_npz':'baseline_native_plus_rgbeta_portals_muF_muS.npz',
      'full_extra_HHSS_reference_json':'exact_extra_HHSS_reference_tensors.json',
      'generated_rgbeta_extra_registry_source':generated,
      'extra_HHSS_registry_entries':extra['registry_HHSS_extra_symbols'],
      'extra_HHSS_monomial_sectors':len(extra['monomial_sectors']),
      'extra_HHSS_coupling_presence_at_both_scales':extra_availability,
      'source_sha256':source_hashes,
      'source_boundary_sha256':boundary.hash,
      'basis_metadata':{
         'portal_basis':'original S1/S2 complex high-weight converted to canonical real components',
         'source_to_muF_scalar_mass_basis':'not asserted in this batch',
         'source_to_muS_scalar_mass_basis':'not asserted in this batch',
         'special_model_basis_needs_independent_phase_check':model in (
             'A_-2','A_+0','B_-1','C_-1','D_-2','D_+0','E_-2','E_+0'),
      },
      'baseline_portal_exact_tensor_reconstruction_at_both_scales':both_scale_portals_complete,
      'full_HHSS_operator_tensor_and_couplings_certified':False,
      'safe_for_existing_seven_term_EFT1_mixing_as_full_input':False,
      'safe_for_existing_joint_scalar_matching_as_full_input':False,
      'physical_C5_calculated':False,
      'missing_physical_inputs':blockers,
      'interpretation':'The 14 extra baseline portal DIRECTIONS are reconstructed in source convention. '+
        'The 24 symmetry-allowed extra HHSS registry symbols have exact SU2 gauge-invariant '+
        'reference spans, not certified named RGBeta/Matchete coefficient maps. '+
        'No unset physical coupling was set to zero.'}
    write_json(output/'portal_and_extra_status.json',report)
    return {
      'model':model,'status':report['status'], 'physics_status':'blocked',
      'baseline_portal_directions_completed':portal['baseline_missing_native_directions_filled'],
      'extra_HHSS_registry_entries':report['extra_HHSS_registry_entries'],
      'extra_HHSS_reference_sector_count':report['extra_HHSS_monomial_sectors'],
      'numerical_baseline_portals_both_scales_verified':both_scale_portals_complete,
      'physical_C5_calculated':False,'blocking_dependency_count':len(blockers)}


def run_completion_batch(root:Path,output:Path,models:tuple[str,...]|None=None)->dict:
    root=Path(root).resolve();output=Path(output).resolve()
    owner=(root/'output/general_observable/lo_certification_couplings').resolve()
    if output!=owner and owner not in output.parents:
        raise ValueError('Workstream 2 output ownership violation')
    if output.exists():
        raise FileExistsError(f'Refusing to overwrite existing outputs: {output}')
    from Numerical.joint_scalar_matching.adapter import ALL_MODELS
    from Numerical.joint_scalar_matching.fermion_boundary_bridge import load_release_manifest,slug
    chosen=tuple(ALL_MODELS if models is None else models)
    if not chosen or len(set(chosen))!=len(chosen) or not set(chosen)<=set(ALL_MODELS):
        raise ValueError('Invalid model selection')
    release=root/'output/general_observable/lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f'
    cache=root/'output/general_observable/lo_joint_scalar_matching/source_quartic_v4/native_quartic_channels'
    if not set(chosen)<=set(load_release_manifest(release)):
        raise ValueError('Incomplete all-F-first source release')
    registry=read_json(root/'RGE/general/generated/general_rgbeta_quartic_registry.json')
    if registry.get('status')!='GeneratedCandidateBasisNotWolframVerified':
        raise ValueError('RGBeta registry certification state changed; re-audit required')
    for model in chosen:
        missing=set(BASELINE_HHSS_REGISTRY[model])-_registry_hhss_names(registry['models'][model])
        if missing:
            raise ValueError(f'{model}: lost previously inventoried extra HHSS {sorted(missing)}')
    versions=check_source_versions(root)
    source_recheck=local_source_recheck(root,chosen,versions,release,cache,registry)
    declarations=rgbeta_potential_declarations(root)
    output.mkdir(parents=True,exist_ok=False)
    rows=[]
    for model in chosen:
        try:
            rows.append(model_completion(root,model,output/slug(model),release,cache,
                                         registry['models'][model],declarations))
        except Exception as exc:
            path=output/slug(model)
            path.mkdir(exist_ok=True)
            failure={'model':model,'status':'blocked_source_or_tensor_audit_failure',
                     'physics_status':'blocked',
                     'error':f'{type(exc).__name__}: {exc}',
                     'physical_C5_calculated':False}
            write_json(path/'failure.json',failure)
            rows.append(failure)
    summary={
      'format':FORMAT,'models_processed':len(rows),
      'all16_requested':len(chosen)==len(ALL_MODELS),
      'status':'conditional_baseline_portals_completed_extra_HHSS_geometry_unmapped',
      'source_recheck':source_recheck,
      'rgbeta_source_signatures':declarations,
      'all16_model_audits_passed':all(x.get('numerical_baseline_portals_both_scales_verified',False)
                                     for x in rows),
      'baseline_portal_directions_completed':sum(x.get('baseline_portal_directions_completed',0) for x in rows),
      'extra_HHSS_registry_entries':sum(x.get('extra_HHSS_registry_entries',0) for x in rows),
      'physical_C5_calculated':False,
      'models':rows}
    write_json(output/'summary.json',summary)
    return summary
