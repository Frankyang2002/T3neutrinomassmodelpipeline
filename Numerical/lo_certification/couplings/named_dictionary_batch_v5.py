"""V5 all-16 native Matchete CG + RGBeta named-quartic source crosswalk.

Read-only consumers of existing V4 and original native EFT1 seeds. Only writes
under output/general_observable/lo_certification_couplings. Fail closed on
source changes, missing tensor or unparseable actual Mathematica SparseArray.
"""
from __future__ import annotations

from pathlib import Path
import hashlib
import json
import re

from .batch import (read_json,write_json,model_paths,BASELINE_HHSS_REGISTRY,
                    _registry_hhss_names,observed_running_states)
from .named_hhss_projection import (named_hhss_sources,generated_declaration,
                                    _strip_parens,SourceInvariantError)
from .source_cg_v5 import parse_all_cgs,inspect_all_terms,MatcheteCGError
from .matchete_rgbeta_dictionary_v5 import (FORMAT,per_named_HHSS_evidence,
                                             summarize_all_extra_registry)


def _sha(p:Path)->str:return hashlib.sha256(p.read_bytes()).hexdigest()


def _slug(model:str)->str:return model.replace('_+','_p').replace('_-','_m')


def _checked_run_dir(owner:Path,name:str)->Path:
    if not name or name in ('.','..') or '/' in name or '\\' in name:
        raise ValueError('Run name must be a single bare directory name')
    return (owner/name).resolve()


def _source_50_declarations(root:Path,model:str,registry:dict)->dict:
    # Completes V4's 24 HHSS source checks to the 50 general quartic entries.
    generated=root/f'RGE/general/generated/{model}_extra_quartics.wl'
    main=root/'RGE/running/rgbeta/T3RGBetaModel.wl'
    texts={'new':generated.read_text(encoding='utf-8'),
           'preexisting':main.read_text(encoding='utf-8')}
    rows=[]
    for category,key in (('new','name'),('preexisting','existing_name')):
        for item in registry.get(category+'_rgbeta_symbols',[]):
            name=item[key]
            decl=generated_declaration(texts[category],name)
            counts=tuple(decl['fields'].count(s) for s in
                         ('H','Bar @ H','S1','Bar @ S1','S2','Bar @ S2'))
            if counts!=tuple(item['counts']):
                raise MatcheteCGError(f'{model}/{name}: registry and generated fields differ')
            if category=='new':
                norm=lambda s:re.sub(r'\s+','',_strip_parens(s))
                if norm(decl['group_invariant_expression'])!=norm(item['invariant']):
                    raise MatcheteCGError(f'{model}/{name}: registry WL invariant drift')
            rows.append({'name':name,'category':category,'counts':list(counts),
                         'source_self_conjugate_option':decl['self_conjugate_option'],
                         'full_WL_declaration_checked':True})
    if len(rows)!=registry.get('total_additional_quartics'):
        raise MatcheteCGError(f'{model}: generator count not equal registry model total')
    return {'total_source_declarations_verified':len(rows),'entries':rows,
            'file_sha256':{'generated_WL':_sha(generated),'rgbeta_core_WL':_sha(main)}}


def _verify_v4(root:Path,v4root:Path,model:str,registry:dict)->tuple[dict,dict]:
    path=v4root/_slug(model)/'named_HHSS_source_projection.json'
    prev=read_json(path)
    if prev.get('model')!=model or prev.get('status')!='source_invariant_projections_conditional':
        raise MatcheteCGError(f'{model}: V4 source report unsuccessful')
    new=named_hhss_sources(root,model,registry)
    keys={r['rgbeta_symbol'] for r in prev['entries']}
    names={r['declaration']['symbol'] for r in new['sectors']}
    if keys!=names or len(keys)!=len(prev['entries']):
        raise MatcheteCGError(f'{model}: V4 symbol inventory mismatch')
    old=prev['source_evidence']['source_sha256']
    if old!=new['source_sha256']:
        raise MatcheteCGError(f'{model}: RGBeta source changed since V4; rerun V4 first')
    return prev,{'V4_source_report_sha256':_sha(path),
                 'RGBeta_generated_extra_WL_sha256':new['source_sha256']['generated_extra_wl'],
                 'RGBeta_core_WL_sha256':new['source_sha256']['core_rgbeta_wl']}


def audit_v5_model(root:Path,v4root:Path,model:str,registry:dict,release:Path,
                   cache:Path)->dict:
    previous,prior_sources=_verify_v4(root,v4root,model,registry)
    src=named_hhss_sources(root,model,registry)
    paths=model_paths(root,model,release,cache)
    seed_path=paths['quartic_seed']
    seed=read_json(seed_path)
    if previous['source_evidence']['native_matchete_seed_sha256']!=_sha(seed_path):
        raise MatcheteCGError(f'{model}: Matchete seed differs from source used for V4')
    if seed.get('ActiveHeavyFields')!=['S1','S2'] or seed.get('IntegratedFields')!=['F']:
        raise MatcheteCGError(f'{model}: invalid EFT1 threshold source')
    cgs=parse_all_cgs(seed)
    native_terms=inspect_all_terms(seed,cgs)
    all50=_source_50_declarations(root,model,registry)
    src50=summarize_all_extra_registry(model,registry,native_terms)
    if src50['additional_quartic_registry_entries']!=all50['total_source_declarations_verified']:
        raise MatcheteCGError('50-invariant manifest mismatch')
    rows=[per_named_HHSS_evidence(model,item,native_terms,cgs,previous)
          for item in src['sectors']]
    if len(rows)!=sum(item['sector_is_HHSS'] for item in src50['entries']):
        raise MatcheteCGError(f'{model}: HHSS inventory undercount')
    running=observed_running_states(paths['running'],model)
    states={'muF':running['eft1_at_fermion_threshold'],
            'muS':running['eft1_at_scalar_threshold']}
    for row in rows:
        symbol=row['rgbeta_symbol']
        row['running_source_presence_at_scales']={scale:{
            'symbol_present':symbol in state,
            'value_as_source_data':state.get(symbol),
            'accepted_as_physical':False} for scale,state in states.items()}
    # Extra numeric couplings are not filled by assigning missing entries zero.
    summary={'format':FORMAT,'model':model,
      'status':'native_CG_and_RGBeta_source_crosswalk_passed_physical_map_blocked',
      'native_Matchete_CG_tensors_exactly_parsed':len(cgs),
      'native_EFT1_scalar_quartic_terms_checked':len(native_terms),
      'all_additional_RGBeta_quartic_declarations_checked':all50['total_source_declarations_verified'],
      'named_HHSS_channels_crosswalked':len(rows),
      'named_HHSS_direct_native_source_matches':sum(bool(r['matching_named_native_terms']) for r in rows),
      'potential_exact_candidate_CG_maps':sum(bool(r['source_based_candidate_tensor_reconstruction']) for r in rows),
      'physically_certified_named_HHSS_maps':0,
      'physical_C5_calculated':False,'numerical_EFT1_full_coupling_background_exported':False,
      'scales_gev':{'muF':running['scales_gev']['F'],'muS':running['scales_gev']['S']},
      'sources_sha256':{**prior_sources,'Matchete_EFT1_quartic_seed_sha256':_sha(seed_path),
                        'native_running_sha256':_sha(paths['running']),
                        'registry_sha256':_sha(root/'RGE/general/generated/general_rgbeta_quartic_registry.json')},
      'native_CGRegistry':[cg.invariant_summary() for cg in cgs.values()],
      'native_quartic_terms':native_terms,
      'all_extra_quartic_source_declarations':all50,
      'all_extra_quartic_sector_crosswalk':src50,
      'named_extra_HHSS_entries':rows,
      'thresholds':{'fermions':'all heavy fermions at muF','scalars':'both scalars at muS'},
      'blocking_dependencies':[
         'source_Matchete_extra_named_CGs_missing_or_unmapped_for_full_physical_dictionary',
         'Matchete_CG_index_orientation_and_complex_phase_match_not_fully_certified',
         'RGBeta_AddQuartic_internal_coefficient_and_HC_factors_not_cross_certified',
         'full_additional_quartic_trajectory_muF_to_muS_not_certified',
         'independent_GammaCC_to_Weyl_normalization_and_operator_closure_blockers']}
    return summary


def run_v5_batch(root:Path,out:Path,models:tuple[str,...]|None=None,
                 v4_run_name:str='ws2_named_hhss_v1')->dict:
    root=Path(root).resolve();out=Path(out).resolve()
    owner=(root/'output/general_observable/lo_certification_couplings').resolve()
    if out==owner or owner not in out.parents:raise ValueError('Workstream 2 ownership violation')
    if out.exists():raise FileExistsError(f'Refusing to overwrite existing output: {out}')
    v4root=_checked_run_dir(owner,v4_run_name)
    from Numerical.joint_scalar_matching.adapter import ALL_MODELS
    chosen=tuple(ALL_MODELS if models is None else models)
    if not chosen or len(chosen)!=len(set(chosen)) or not set(chosen).issubset(ALL_MODELS):
        raise ValueError('Invalid model selection')
    last=read_json(v4root/'summary.json')
    if last.get('models_processed')!=16 or not last.get('all16_source_projections_passed'):
        raise ValueError('V4 all-16 source projection evidence is required')
    registry_path=root/'RGE/general/generated/general_rgbeta_quartic_registry.json'
    allreg=read_json(registry_path)
    if allreg.get('status')!='GeneratedCandidateBasisNotWolframVerified':
        raise ValueError('RGBeta registry status changed; re-audit first')
    release=root/'output/general_observable/lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f'
    cache=root/'output/general_observable/lo_joint_scalar_matching/source_quartic_v4/native_quartic_channels'
    for model in chosen:
        baseline=set(BASELINE_HHSS_REGISTRY[model]);existing=_registry_hhss_names(allreg['models'][model])
        if baseline-existing:raise ValueError(f'{model}: registry lost previously checked HHSS symbols')
    out.mkdir(parents=True,exist_ok=False)
    rows=[]
    for model in chosen:
        sub=out/_slug(model);sub.mkdir()
        try:
            report=audit_v5_model(root,v4root,model,allreg['models'][model],release,cache)
            write_json(sub/'matchete_rgbeta_exact_source_crosswalk.json',report)
            rows.append({'model':model,'status':report['status'],
                         'native_CG_count':report['native_Matchete_CG_tensors_exactly_parsed'],
                         'native_quartic_term_count':report['native_EFT1_scalar_quartic_terms_checked'],
                         'extra_quartic_declarations':report['all_additional_RGBeta_quartic_declarations_checked'],
                         'HHSS_source_symbols':report['named_HHSS_channels_crosswalked'],
                         'matching_named_native_HHSS':report['named_HHSS_direct_native_source_matches'],
                         'exact_candidate_CG_maps':report['potential_exact_candidate_CG_maps'],
                         'physical_map_count':0,'physical_status':'blocked'})
        except Exception as exc:
            fail={'model':model,'status':'source_or_exact_CG_audit_failed',
                  'error':f'{type(exc).__name__}: {exc}','physical_status':'blocked'}
            write_json(sub/'failure.json',fail)
            rows.append(fail)
    passed=all(r.get('status')=='native_CG_and_RGBeta_source_crosswalk_passed_physical_map_blocked' for r in rows)
    totalHHSS=sum(r.get('HHSS_source_symbols',0) for r in rows)
    total50=sum(r.get('extra_quartic_declarations',0) for r in rows)
    if len(chosen)==16 and (total50!=50 or totalHHSS!=24):
        passed=False  # Changed source inventory requires explicit re-audit.
    summary={'format':FORMAT,'models_processed':len(rows),'all16_requested':len(chosen)==16,
             'all16_source_and_CG_crosswalk_passed':passed and len(chosen)==16,
             'extra_RGBeta_quartic_declarations_checked':total50,
             'extra_HHSS_source_symbols_crosswalked':totalHHSS,
             'native_Matchete_CG_tensors_parsed':sum(r.get('native_CG_count',0) for r in rows),
             'native_EFT1_quartic_terms_parsed':sum(r.get('native_quartic_term_count',0) for r in rows),
             'matching_native_extra_HHSS_channels':sum(r.get('matching_named_native_HHSS',0) for r in rows),
             'physically_certified_HHSS_coupling_maps':0,
             'full_physical_EFT1_couplings_certified':False,'physical_C5_calculated':False,
             'input_V4_run_name':v4_run_name,'RGBeta_registry_sha256':_sha(registry_path),
             'models':rows}
    write_json(out/'summary.json',summary)
    return summary
