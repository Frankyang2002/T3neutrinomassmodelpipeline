"""All-16 read-only native-Matchete versus RGBeta HHSS portal normalization audit.

Only Workstream-2 output is written, exclusively. The existing Matchete parser,
CG sources and prior local version check are reused. No fabricated extra-channel
values, group invariants, phase choices or physical matching coefficients.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from .batch import (MODEL_DIMS, BASELINE_HHSS_REGISTRY, check_source_versions,
                    local_source_recheck, model_paths, read_json, digest,
                    write_json, nonzero_sources, _registry_hhss_names)
from .portal_normalization import (portal_species_audit,
                                   construct_rgbeta_portal_background,
                                   numerically_map_full_portal_subsector)

RESULT_FORMAT='t3-ws2-exact-portal-span-audit-v2'


def rgbeta_potential_declarations(root:Path)->dict:
    """Check actual local RGBeta source defines the intended portal polynomials.

    These literal declarations are read and hashed, never constructed from
    approximate beta-functions or candidate numerical aliases.
    """
    p=root/'RGE/running/rgbeta/T3RGBetaModel.wl'
    source=p.read_text(encoding='utf-8')
    import re
    signatures={
      'lambdaH1':r'AddQuartic\[\s*lambdaH1\s*,\s*\{Bar\s*@\s*H,\s*H,\s*Bar\s*@\s*S1,\s*S1\}',
      'lambdaH2':r'AddQuartic\[\s*lambdaH2\s*,\s*\{Bar\s*@\s*H,\s*H,\s*Bar\s*@\s*S2,\s*S2\}',
      'lambdaH1Adj':r'AddQuartic\[\s*lambdaH1Adj\s*,\s*\{Bar\s*@\s*H,\s*H,\s*Bar\s*@\s*S1,\s*S1\}',
      'lambdaH2Adj':r'AddQuartic\[\s*lambdaH2Adj\s*,\s*\{Bar\s*@\s*H,\s*H,\s*Bar\s*@\s*S2,\s*S2\}',
    }
    seen={name:bool(re.search(exp,source)) for name,exp in signatures.items()}
    if not all(seen.values()):
        raise ValueError(f'RGBeta portal source declaration absent/changed: {seen}')
    # Source-level basis factors, not a claim about the RGBeta package's run-time AddQuartic semantics.
    for expression in ('T3RGBetaPortalAdjInvariant[2]',
                       'T3RGBetaPortalAdjInvariant[3]',
                       'tGen[SU2L @ fund', 'tGen[SU2L @ S2',
                       'T3RGBetaPairInvariant[1', 'T3RGBetaPairInvariant[2',
                       'T3RGBetaPairInvariant[3'):
        if expression not in source:
            raise ValueError(f'RGBeta source invariant signature changed: {expression}')
    mat=root/'Lagrangian/interactions/ScalarPotential.wl'
    mat_source=mat.read_text(encoding='utf-8')
    for signature in ('HiggsPortal[which_, d_, coupling_]',
                      'coupling[] Bar[H[i]] H[i] ScalarNorm[which, d, j]',
                      'ScalarNorm[which_Integer, d_Integer?Positive, i_]'):
        if signature not in mat_source:
            raise ValueError('Matchete source portal/potential convention changed: '+signature)
    return {'path':str(p.relative_to(root)).replace('\\','/'),
            'sha256':digest(p),'portal_declarations_present':seen,
            'matchete_portal_source_path':str(mat.relative_to(root)).replace('\\','/'),
            'matchete_portal_source_sha256':digest(mat),
            'verified_scope':'literal Matchete and RGBeta source declarations; no unverified RGBeta runtime tensor factors'}


def _validate_generated_extra(root:Path, model:str, registry:dict)->dict:
    from Numerical.joint_scalar_matching.fermion_boundary_bridge import slug
    path=root/f'RGE/general/generated/{model[0]}_{model[2:]}_extra_quartics.wl'
    if not path.is_file():
        raise FileNotFoundError(f'Missing local generated RGBeta extra invariant source: {path}')
    source=path.read_text(encoding='utf-8')
    names=[]
    for row in registry['new_rgbeta_symbols']:
        name=row['name']
        if f'AddQuartic[{name},' not in source:
            raise ValueError(f'{model}: generated RGBeta invariant {name} not found in {path}')
        if sum(row['counts'][:2])==2 and sum(row['counts'][2:])==2:
            names.append(name)
    return {'relative_path':str(path.relative_to(root)).replace('\\','/'),
            'sha256':digest(path),'HHSS_new_symbol_declarations_verified':sorted(names),
            'native_matchete_tensor_normalization_certified':False}


def _inventory(root:Path, model:str, registry:dict, source, ch, f:dict, s:dict)->dict:
    """Explicitly distinguish portal directions absent from the native seed from
    the generated candidate extras. Neither class may disappear from the inventory.
    """
    native_names=set(ch.by_name)
    extras=nonzero_sources(registry,native_names)
    generated=_validate_generated_extra(root,model,registry)
    missing_values=[]
    for row in extras:
        symbol=row['rgbeta']
        if symbol not in f or symbol not in s:
            missing_values.append({'symbol':symbol,
                                   'at_muF':symbol in f,
                                   'at_muS':symbol in s,
                                   'status':'not_available_as_complete_two-scale_EFT1_input'})
    return {'extra_HHSS_candidates':extras,'generated_extra_source':generated,
            'extra_missing_same_scale_numeric_values':missing_values,
            'physical_extra_HHSS_tensor_normalization_verified':False,
            'quadratic_registry_blocked_count':registry.get('quadratics_not_implemented'),
            'numerical_extra_HHSS_default_zero_used':False}


def audit_model(root:Path,out:Path,model:str,release:Path,cache:Path,reg:dict,
                declarations:dict)->dict:
    from Numerical.joint_scalar_matching.fermion_boundary_bridge import load_source,slug
    from Numerical.joint_scalar_matching.numerical_benchmark_v5 import read_v4_cached_channels
    from Numerical.joint_scalar_matching.source_quartic_v4 import source_seed_path
    from .batch import observed_running_states
    paths=model_paths(root,model,release,cache)
    source=load_source(paths['boundary'], expected_model=model)
    native=read_v4_cached_channels(paths['quartic_cache'],paths['quartic_seed'],source)
    running=observed_running_states(paths['running'],model)
    f=running['eft1_at_fermion_threshold']; s=running['eft1_at_scalar_threshold']
    r1=portal_species_audit(native,1)
    r2=portal_species_audit(native,2)
    bg,arrays=construct_rgbeta_portal_background(native,f,s)
    normalized_numeric={
        'muF':{f'S{i}':numerically_map_full_portal_subsector(native,i,f)
               for i in (1,2)},
        'muS':{f'S{i}':numerically_map_full_portal_subsector(native,i,s)
               for i in (1,2)},
    }
    # The 4-leg tensor is in the RGBeta highest-weight reference basis.
    # It cannot be contracted with the Matchete LLSS source unless a common
    # scalar coordinate frame has been independently validated.
    inventory=_inventory(root,model,reg,source,native,f,s)
    missing_portal_dirs=sum((x['rgbeta_portal_directions_missing_from_native_span'] or 0)
                            for x in (r1,r2))
    unresolved_basis=[x['species'] for x in (r1,r2)
                      if x['status']=='unresolved_generator_or_complex_scalar_basis']
    extra_names=inventory['extra_HHSS_candidates']
    blockers=['GammaCC_to_master_Weyl_signed_normalization_not_certified',
              'full_SM_Weyl_yukawa_gauge_tensor_embedding_not_certified',
              'mixing_lambdaT3_Matchete_RGBeta_phase_normalization_not_certified',
              'full_joint_scalar_background_and_hard_one_loop_closure_not_certified']
    if missing_portal_dirs:
        blockers.append(f'ordinary_RGBeta_portal_span_missing_{missing_portal_dirs}_independent_directions')
    if unresolved_basis:
        blockers.append('portal_Matchete_RGBeta_complex_SU2_component_frame_unverified:'+','.join(unresolved_basis))
    if extra_names:
        blockers.append('extra_HHSS_invariants_not_source_normalized:'+str(len(extra_names)))
    if inventory['quadratic_registry_blocked_count']:
        blockers.append('extra_scalar_quadratic_RG_closure_unimplemented')
    if inventory['extra_missing_same_scale_numeric_values']:
        blockers.append('extra_HHSS_numeric_running_values_missing')
    native_ready=all(x['normalization_conversion_physically_usable_for_portal_subsector']
                     for x in (r1,r2))
    if any(scal['status']!='exact_full_tensor_span' for scal in (r1,r2)):
        native_ready=False
    # No physical values are invented: only the RGBeta portal terms with
    # explicitly available muF/muS numeric couplings are exported.
    out.mkdir(parents=True,exist_ok=False)
    np.savez_compressed(out/'rgbeta_portals_only_muF_muS.npz',**arrays)
    report={'format':RESULT_FORMAT, 'model':model,
        'status':'exact_full_portal_span_only' if native_ready else 'conditional_portal_span_or_index_map',
        'physics_status':'blocked', 'scales_gev':{'mu_F':running['scales_gev']['F'],
                                               'mu_S':running['scales_gev']['S']},
        'scheme':'MSbar', 'two_grouped_matching_thresholds_only':True,
        'original_full_Matchete_seed':str(paths['quartic_seed'].relative_to(root)).replace('\\','/'),
        'source_sha256':{key:digest(value) for key,value in paths.items()},
        'rgbeta_local_potential_source':declarations,
        'native_s1_portal':r1,'native_s2_portal':r2,
        'exact_converted_native_portal_couplings':normalized_numeric,
        'rgbeta_portal_partial_background':bg,
        'partial_background_npz':'rgbeta_portals_only_muF_muS.npz',
        'extra_operator_inventory':inventory,
        'matchete_portal_span_missing_directions':missing_portal_dirs,
        'matchete_baseline_plus_extra_HHSS_complete':False,
        'safe_as_existing_seven_term_solver_input':False,
        'safe_as_existing_joint_scalar_matching_input':False,
        'physical_C5_calculated':False,'blockers':blockers,
        'interpretation':'Exact tensor identities use V=-L4 and 1/4! real-field normalization. '
                         'No source-to-Matchete-GammaCC normalization and no scalar phase rotation inferred.'}
    write_json(out/'portal_normalization.json',report)
    return {'model':model,'status':report['status'],'physical_status':'blocked',
            'RGBeta_portal_directions_missing_from_native_seed':missing_portal_dirs,
            'extra_HHSS_candidate_count':len(extra_names),
            'S1_span_status':r1['status'],'S2_span_status':r2['status'],
            'blockers':blockers,'source_verified':True}


def run_portal_batch(root:Path,output:Path,models:tuple[str,...]|None=None)->dict:
    root=Path(root).resolve();output=Path(output).resolve()
    owner=(root/'output/general_observable/lo_certification_couplings').resolve()
    if owner!=output and owner not in output.parents:
        raise ValueError('Workstream 2 output ownership violation')
    if output.exists():
        raise FileExistsError(f'Run directory exists: {output}')
    from Numerical.joint_scalar_matching.adapter import ALL_MODELS
    from Numerical.joint_scalar_matching.fermion_boundary_bridge import load_release_manifest,slug
    chosen=tuple(ALL_MODELS if models is None else models)
    if not chosen or len(set(chosen))!=len(chosen) or not set(chosen)<=set(ALL_MODELS):
        raise ValueError('Invalid all-16 model selection')
    release=root/'output/general_observable/lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f'
    cache=root/'output/general_observable/lo_joint_scalar_matching/source_quartic_v4/native_quartic_channels'
    manifests=load_release_manifest(release)
    if not set(chosen)<=set(manifests):
        raise ValueError('Fermion source release manifest lacks required models')
    reg=read_json(root/'RGE/general/generated/general_rgbeta_quartic_registry.json')
    if reg.get('status')!='GeneratedCandidateBasisNotWolframVerified':
        raise ValueError('RGBeta registry has changed its certification category')
    for model in chosen:
        existing=_registry_hhss_names(reg['models'][model]);
        missing=set(BASELINE_HHSS_REGISTRY[model])-existing
        if missing:
            raise ValueError(f'{model}: missing previously inventoried HHSS {sorted(missing)}')
    versions=check_source_versions(root)
    recheck=local_source_recheck(root,chosen,versions,release,cache,reg)
    declarations=rgbeta_potential_declarations(root)
    output.mkdir(parents=True,exist_ok=False)
    rows=[]
    for model in chosen:
        path=output/slug(model)
        try:
            rows.append(audit_model(root,path,model,release,cache,reg['models'][model],declarations))
        except Exception as error:
            # Each failed model still counts and remains blocked; no default tensor
            # is substituted and no source input is rewritten.
            path.mkdir(exist_ok=True)
            row={'model':model,'status':'blocked_audit_error','physical_status':'blocked',
                 'source_verified':False,'error':f'{type(error).__name__}: {error}'}
            write_json(path/'failure.json',row)
            rows.append(row)
    summary={'format':RESULT_FORMAT,'models_processed':len(rows),
             'all16_requested':len(chosen)==len(ALL_MODELS),
             'thresholds':['all heavy BSM fermions together at mu_F',
                           'both heavy BSM scalar multiplets together at mu_S'],
             'status':'source_derived_partial_dictionary_not_physical_C5',
             'source_drift_recheck':recheck,'rgbeta_potential_source':declarations,
             'all16_audit_passed':all(x.get('source_verified',False) for x in rows),
             'verified_full_rgbeta_portal_span_models':sum(x.get('status')=='exact_full_portal_span_only' for x in rows),
             'newly_discovered_baseline_portal_deficiency':sum(x.get('RGBeta_portal_directions_missing_from_native_seed',0) for x in rows),
             'extra_HHSS_registry_entries':sum(x.get('extra_HHSS_candidate_count',0) for x in rows),
             'physical_C5_calculated':False,'models':rows}
    write_json(output/'summary.json',summary)
    return summary
