"""V4 source-explicit 24-name RGBeta HHSS index projection audit.

Only owned Workstream-2 outputs. Reads the successful V3 results; validates
all local generator expressions rather than inferring equivalences from
coupling names. No RGE evolution or Wilson matching is duplicated.
"""
from __future__ import annotations

from pathlib import Path
import hashlib
import json
import re

from .batch import read_json, write_json, model_paths, BASELINE_HHSS_REGISTRY, _registry_hhss_names, observed_running_states
from .extra_hhss_geometry import model_extra_geometry
from .named_hhss_projection import FORMAT, named_hhss_sources, projection_report, SourceInvariantError


def _sha(path:Path)->str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _model_slugs(model:str)->str:
    return model.replace('_+', '_p').replace('_-','_m')


def _check_v3_model(v3:Path,model:str,registry:dict)->dict:
    extra_file=v3/'exact_extra_HHSS_reference_tensors.json'
    status_file=v3/'portal_and_extra_status.json'
    prev=read_json(extra_file)
    prev_status=read_json(status_file)
    check=model_extra_geometry(model,registry,include_components=False)
    if prev.get('model')!=model or prev_status.get('model')!=model:
        raise SourceInvariantError('V3 model identity mismatch')
    if prev.get('registry_HHSS_extra_symbols')!=check['registry_HHSS_extra_symbols']:
        raise SourceInvariantError('V3 reference invariant inventory differs from local source')
    prior={(tuple(row['counts_H_BarH_S1_BarS1_S2_BarS2']),
            tuple(sorted(s['rgbeta_symbol'] for s in row['source_registry_symbols'])))
            for row in prev['monomial_sectors']}
    now={(tuple(row['counts_H_BarH_S1_BarS1_S2_BarS2']),
          tuple(sorted(s['rgbeta_symbol'] for s in row['source_registry_symbols'])))
          for row in check['monomial_sectors']}
    if prior!=now:raise SourceInvariantError('Named V3 registry symbols changed')
    return {'v3_reference_sha256':_sha(extra_file),'v3_status_sha256':_sha(status_file),
            'extra_symbols':check['registry_HHSS_extra_symbols']}


def _matchete_name_evidence(seed:dict,entry:dict,category:str)->dict:
    name=entry.get('matchete_basis_name') if category=='new' else None
    terms=seed.get('ScalarQuarticTerms',[])
    if not isinstance(terms,list):raise SourceInvariantError('Unexpected Matchete quartic seed')
    if not name:
        return {'matchete_candidate_pattern':None,'matching_native_terms':[],
                'verified_named_CG_to_RGBeta_map':False}
    pattern=re.escape(name).replace(r'\#',r'\d+')
    hits=[{'term_index':row.get('Index'), 'term_sha256':hashlib.sha256(
              row['TermInputForm'].encode('utf-8')).hexdigest(),
           'CGNames':row.get('CGNames',[])}
          for row in terms if re.search(pattern,row.get('TermInputForm',''))]
    return {'matchete_candidate_pattern':name,'matching_native_terms':hits,
            'verified_named_CG_to_RGBeta_map':False}


def audit_named_model(root:Path,v3root:Path,model:str,registry:dict)->dict:
    v3=v3root/_model_slugs(model)
    prior=_check_v3_model(v3,model,registry)
    sources=named_hhss_sources(root,model,registry)
    # Native seed is the actual symbolic EFT1 source, not an assumed field list.
    release=root/'output/general_observable/lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f'
    cache=root/'output/general_observable/lo_joint_scalar_matching/source_quartic_v4/native_quartic_channels'
    paths=model_paths(root,model,release,cache)
    seed=read_json(paths['quartic_seed'])
    running=observed_running_states(paths['running'],model)
    states={'muF':running['eft1_at_fermion_threshold'],
            'muS':running['eft1_at_scalar_threshold']}
    rows=[]
    bycat={k:{(i['name'] if k=='new' else i['existing_name']):i
             for i in registry.get(k+'_rgbeta_symbols',[])}
           for k in ('new','preexisting')}
    for item in sources['sectors']:
        category=item['category'];decl=item['declaration'];name=decl['symbol']
        report=projection_report(model,item['counts'],decl)
        report['source_category']=category
        report['matchete_seed_link']=_matchete_name_evidence(seed,bycat[category][name],category)
        report['numerical_running_coefficients']={scale:{
          'symbol_present':name in st,
          'value_as_source_data':(st[name] if name in st else None),
          'accepted_as_physical':False} for scale,st in states.items()}
        report['physical_C5_certified']=False
        rows.append(report)
    if len(rows)!=prior['extra_symbols']:
        raise SourceInvariantError('Incomplete 24-sector inventory for model')
    return {'format':FORMAT,'model':model,'status':'source_invariant_projections_conditional',
            'source_projection_count':len(rows),
            'valid_canonical_orientation_symbols':sum(bool(r['valid_canonical_invariant_orientation_count']) for r in rows),
            'orientation_ambiguous_symbols':sum(r['valid_canonical_invariant_orientation_count']!=1 for r in rows),
            'physical_named_Matchete_RGBeta_maps_certified':0,
            'source_evidence':{'source_sha256':sources['source_sha256'],
              'native_matchete_seed_sha256':_sha(paths['quartic_seed']),**prior},
            'muF_gev':running['scales_gev']['F'],
            'muS_gev':running['scales_gev']['S'],
            'entries':rows,
            'physical_C5_calculated':False,
            'full_physical_eft1_quartics_certified':False,
            'blocking_dependencies':[
             'RGBeta_AddQuartic_internal_coefficient_normalization_not_cross_certified',
             'source_index_dual_map_not_physically_fixed_by_Ward_identity_alone',
             'Matchete_CG_to_RGBeta_named_extra_tensor_transform_not_verified',
             'missing_complete_extra_coupling_evolution_muF_to_muS',
             'remaining_spinor_and_scalar_basis_and_operator_closure_dependencies'],
            'two_thresholds_only':True,
            'threshold_sequence':'UV -> all heavy fermions at muF -> EFT1 -> both scalars at muS -> SMEFT'}


def run_named_hhss_batch(root:Path,out:Path,models:tuple[str,...]|None=None,
                         v3_run_name:str='ws2_quartic_completion_v1')->dict:
    root=Path(root).resolve();out=Path(out).resolve()
    owner=(root/'output/general_observable/lo_certification_couplings').resolve()
    if out==owner or owner not in out.parents:
        raise ValueError('Workstream 2 output ownership violation')
    if out.exists():raise FileExistsError(f'Refusing overwrite: {out}')
    from Numerical.joint_scalar_matching.adapter import ALL_MODELS
    chosen=tuple(ALL_MODELS if models is None else models)
    if not chosen or len(chosen)!=len(set(chosen)) or not set(chosen).issubset(ALL_MODELS):
        raise ValueError('Invalid model selection')
    v3root=(owner/v3_run_name).resolve()
    if v3root.parent!=owner or v3_run_name in ('.','..'):
        raise ValueError('V3 dependency must be an existing Workstream 2 run')
    meta=read_json(v3root/'summary.json')
    if not meta.get('all16_model_audits_passed') or meta.get('models_processed')!=16:
        raise ValueError('V3 source audit not successful for all 16 models')
    registry=read_json(root/'RGE/general/generated/general_rgbeta_quartic_registry.json')
    if registry.get('status')!='GeneratedCandidateBasisNotWolframVerified':
        raise ValueError('Registry certification status changed; source re-audit required')
    for model in chosen:
        missing=set(BASELINE_HHSS_REGISTRY[model])-_registry_hhss_names(registry['models'][model])
        if missing:raise ValueError(f'{model}: missing prior HHSS inventory {sorted(missing)}')
    out.mkdir(parents=True,exist_ok=False)
    rows=[]
    for model in chosen:
        slug=_model_slugs(model)
        try:
            report=audit_named_model(root,v3root,model,registry['models'][model])
            path=out/slug;path.mkdir()
            write_json(path/'named_HHSS_source_projection.json',report)
            rows.append({'model':model,'status':report['status'],
                         'source_projection_count':report['source_projection_count'],
                         'valid_canonical_orientation_symbols':report['valid_canonical_orientation_symbols'],
                         'orientation_ambiguous_symbols':report['orientation_ambiguous_symbols'],
                         'physical_named_maps':0,'physical_status':'blocked'})
        except Exception as exc:
            path=out/slug;path.mkdir(exist_ok=True)
            fail={'model':model,'status':'blocked_source_or_projection_failure',
                  'error':f'{type(exc).__name__}: {exc}', 'physical_status':'blocked'}
            write_json(path/'failure.json',fail)
            rows.append(fail)
    summary={'format':FORMAT,'all16_requested':len(chosen)==16,
             'all16_source_projections_passed':all(r.get('status')=='source_invariant_projections_conditional' for r in rows),
             'models_processed':len(rows),
             'extra_HHSS_source_symbols_examined':sum(r.get('source_projection_count',0) for r in rows),
             'valid_canonical_orientation_symbols':sum(r.get('valid_canonical_orientation_symbols',0) for r in rows),
             'physical_named_maps_certified':0,'full_physical_eft1_couplings_certified':False,
             'physical_C5_calculated':False,'input_V3_run_name':v3_run_name,
             'registry_sha256':_sha(root/'RGE/general/generated/general_rgbeta_quartic_registry.json'),
             'models':rows}
    write_json(out/'summary.json',summary)
    return summary
