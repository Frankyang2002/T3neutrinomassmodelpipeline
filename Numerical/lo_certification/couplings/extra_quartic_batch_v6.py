"""V6 all-16 source-operator completion runner, read-only repository inspection.

Optionally invokes an isolated, live Matchete `GetCGTensor` export (NOT matching).
Never copies UV extra CG directly into EFT1 numerical mixing/matching arrays.
"""
from __future__ import annotations

from pathlib import Path
import json
import os
import shutil
import subprocess

from .extra_quartic_source_v6 import (
    ALL_MODELS, FORMAT, OUTPUT_OWNER, V6SourceError, alpha_tag,
    generator_contract, file_sha256, indexed_seed_evidence, find_seed_files,
    model_registry_quartics, orientation_candidate_ratios,
    parse_live_operator_export, slug, experimental_inputform_evidence,
)


def _read_json(p:Path):
    return json.loads(p.read_text(encoding='utf8'))


def _write(p:Path,doc:dict):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(doc,indent=2,sort_keys=True,default=str)+'\n',encoding='utf8')


def _validate_output(root:Path,out:Path):
    root=root.resolve();out=out.resolve()
    owner=(root/OUTPUT_OWNER).resolve()
    if out==owner or owner not in out.parents:
        raise ValueError('Workstream 2 output ownership violation')
    if out.exists():raise FileExistsError(f'Refusing to overwrite V6 output: {out}')


def _load_v5(root:Path,run_name:str)->dict:
    if not run_name or '/' in run_name or '\\' in run_name or run_name in ('.','..'):
        raise V6SourceError('Invalid source run name')
    p=root/OUTPUT_OWNER/run_name/'summary.json'
    data=_read_json(p)
    if data.get('models_processed')!=16 or not data.get('all16_source_and_CG_crosswalk_passed'):
        raise V6SourceError('Successful V5 all-16 report required')
    return data


def _validate_registry(registry:dict):
    if registry.get('status')!='GeneratedCandidateBasisNotWolframVerified':
        raise V6SourceError('RGBeta registry status changed; cannot trust old inventory')
    if set(registry.get('models',{}))!=set(ALL_MODELS):
        raise V6SourceError('Registry no longer contains exactly 16 expected models')
    entries=sum(len(row['new_rgbeta_symbols'])+len(row['preexisting_rgbeta_symbols'])
                for row in registry['models'].values())
    hhss=sum(sum(item['counts'][:2])==2
             for row in registry['models'].values()
             for item in row['new_rgbeta_symbols']+row['preexisting_rgbeta_symbols'])
    if entries!=50 or hhss!=24:
        raise V6SourceError('The audited registry 50/24 quartic inventory changed')


def _source_scan(root:Path,model:str,expected:list[dict],paths:list[Path]):
    rows=[]
    for path in find_seed_files(root,model,paths):
        try:
            evidence=indexed_seed_evidence(path,model,expected)
            evidence['path']=str(path.relative_to(root))
            evidence['status']='source_read_and_native_CG_validated'
        except Exception as exc:
            evidence={'path':str(path.relative_to(root)),
                      'status':'seed_invalid_or_schema_changed',
                      'error':f'{type(exc).__name__}: {exc}'}
        rows.append(evidence)
    return rows


def _run_wolfram(root:Path,out:Path,model:str,command:str,timeout:int)->dict:
    script=root/'scripts/RunLOExtraQuarticCGExport.wl'
    if not script.exists():raise V6SourceError('Owned V6 Wolfram exporter missing')
    exe=shutil.which(command) if not Path(command).is_file() else command
    if exe is None:return {'executed':False,'status':'wolframscript_not_installed'}
    destination=out/slug(model)/'live_extra_quartic_CG_export.json'
    environment=os.environ.copy()
    environment.update({'T3_EXPERIMENTAL_GENERAL_QUARTICS':'1',
                        'T3_EXPERIMENTAL_GENERAL_QUADRATICS':'0',
                        'T3_EXPERIMENTAL_GENERAL_FERMIONS':'0',
                        'T3_EXPERIMENTAL_FLAVOR_MASSES':'0'})
    args=[str(exe),'-file',str(script),str(destination),model[0],alpha_tag(model)]
    try:
        completed=subprocess.run(args,cwd=str(root),env=environment,capture_output=True,
                                 text=True,timeout=timeout,check=False)
        return {'executed':True,'status':'returned','returncode':completed.returncode,
                'stdout_tail':completed.stdout[-3000:],'stderr_tail':completed.stderr[-3000:],
                'live_file_exists':destination.exists()}
    except subprocess.TimeoutExpired:
        return {'executed':True,'status':'timed_out','timeout_seconds':timeout}


def _v4_source_doc(root:Path,v4_run_name:str,model:str):
    if not v4_run_name or '/' in v4_run_name or '\\' in v4_run_name or v4_run_name in ('.','..'):
        raise V6SourceError('Invalid V4 run name')
    path=root/OUTPUT_OWNER/v4_run_name/slug(model)/'named_HHSS_source_projection.json'
    doc=_read_json(path)
    if doc.get('model')!=model or doc.get('status')!='source_invariant_projections_conditional':
        raise V6SourceError(f'{model}: prior V4 named projection is not validated')
    return doc


def run_v6_batch(root:Path,out:Path,*,models:tuple[str,...]|None=None,
                 previous_v5:str='ws2_matchete_rgbeta_dictionary_v1',
                 v4_run_name:str='ws2_named_hhss_v1',
                 live:bool=False, wolframscript:str='wolframscript',
                 live_export_root:Path|None=None,probe_orientations:bool=False,
                 wolfram_timeout_seconds:int=900)->dict:
    root=Path(root).resolve();out=Path(out).resolve()
    _validate_output(root,out)
    chosen=ALL_MODELS if models is None else tuple(models)
    if not chosen or len(set(chosen))!=len(chosen) or not set(chosen).issubset(ALL_MODELS):
        raise V6SourceError('Invalid or duplicate T3 model selection')
    if live and live_export_root is not None:
        raise V6SourceError('Choose either --run-live or --live-export-root')
    if probe_orientations and not (live or live_export_root):
        raise V6SourceError('Orientation probe requires actual live Matchete CG exports')
    if not 10<=wolfram_timeout_seconds<=7200:
        raise V6SourceError('Invalid per-model Wolfram timeout')
    previous=_load_v5(root,previous_v5)
    registry_path=root/'RGE/general/generated/general_rgbeta_quartic_registry.json'
    registry=_read_json(registry_path)
    _validate_registry(registry)
    generator=generator_contract(root)
    # Scan all native and experimental sources without overwriting anything.
    sources=sorted(p for p in (root/'output').rglob('eft1_after_F_scalar_quartic_seed.json')
                   if OUTPUT_OWNER not in p.relative_to(root).parents)
    experimental_exports=sorted((root/'output').rglob('*.interaction.inputform.txt'))
    out.mkdir(parents=True,exist_ok=False)
    _write(out/'source_generator_contract.json',generator)
    rows=[]
    for model in chosen:
        model_dir=out/slug(model);model_dir.mkdir()
        try:
            expected=model_registry_quartics(model,registry)
            evidence=_source_scan(root,model,expected,sources)
            hhss=[x for x in expected if x['sector_HHSS']]
            experimental=[]
            for path in experimental_exports:
                if not find_seed_files(root,model,[path]):continue
                try:
                    evidence_text=experimental_inputform_evidence(path,model,expected)
                    evidence_text['path']=str(path.relative_to(root))
                    experimental.append(evidence_text)
                except Exception as exc:
                    experimental.append({'path':str(path.relative_to(root)),
                        'status':'invalid_experimental_source',
                        'error':f'{type(exc).__name__}: {exc}'})
            found=sum(any(any(t['matching_named_extra_term_indices']
                              for t in s.get('entries',[]) if t['rgbeta_symbol']==e['rgbeta_symbol'])
                          for s in evidence) for e in hhss)
            live_result={'executed':False,'status':'not_requested'}
            if live:
                live_result=_run_wolfram(root,out,model,wolframscript,wolfram_timeout_seconds)
            if live_export_root is not None:
                live_path=Path(live_export_root).resolve()/slug(model)/'live_extra_quartic_CG_export.json'
                if not live_path.is_file():
                    raise V6SourceError(f'{model}: requested live export does not exist: {live_path}')
                live_result={'executed':False,'status':'external_live_export_provided',
                             'path':str(live_path)}
            else:
                live_path=model_dir/'live_extra_quartic_CG_export.json'
            parsed=None;orientation=[]
            if live_path.is_file():
                parsed=parse_live_operator_export(_read_json(live_path),model,expected)
                live_result.update({'source_validated':True,
                    'file_sha256':file_sha256(live_path),
                    'extra_operators':parsed['live_operator_count'],
                    'native_CGs':parsed['live_CG_count']})
                if probe_orientations:
                    prior=_v4_source_doc(root,v4_run_name,model)
                    orientation=orientation_candidate_ratios(model,expected,parsed,prior)
                    _write(model_dir/'conditional_tensor_ratio_candidates.json',
                           {'model':model,'entries':orientation,
                            'certified_physical_conversion':False})
            # Avoid serializing internal ExactCG Python objects.
            _write(model_dir/'live_operator_validation.json',
                {'model':model,'export':live_result,
                 'validated_native_operator_list':parsed['live_operators'] if parsed else None,
                 'physical_mapping_certified':False})
            status='live_native_extra_CG_export_verified_conditional' if parsed else 'generator_source_verified_live_CG_not_available'
            report={'format':FORMAT,'model':model,'status':status,
                'RGBeta_additional_quartics':expected,
                'extra_HHSS_channels':len(hhss),
                'scanned_EFT1_quartic_seeds':evidence,
                'scanned_experimental_UV_inputform_exports':experimental,
                'native_extra_HHSS_symbols_found_in_any_EFT1_seed':found,
                'live_native_extra_operator_count':parsed['live_operator_count'] if parsed else 0,
                'live_native_extra_CG_tensor_count':parsed['live_CG_count'] if parsed else 0,
                'live_source_export_status':live_result,
                'formal_orientation_ratio_probes':sum(x['candidate_count'] for x in orientation),
                'UV_extra_source_is_not_EFT1_C5_matching_source':True,
                'physical_HHSS_conversion_certified':False,
                'complete_numeric_HHSS_at_muF':False,'complete_numeric_HHSS_at_muS':False,
                'physical_C5_calculated':False,
                'blocking_dependencies':[
                  'EFT1_matching_seed_must_contain_verified_extra_operator_boundary',
                  'Matchete_and_RGBeta_complex_index_phase_and_AddQuartic_HC_normalisation_uncertified',
                  'extra_RGBeta_coupling_beta_functions_and_two_scale_values_not_certified',
                  'independent_GammaCC_spinor_normalization_and_EFT1_operator_closure_blocked'],
                'sources_sha256':{'registry':file_sha256(registry_path),
                    **generator['source_sha256']}}
            _write(model_dir/'V6_source_evidence.json',report)
            rows.append({'model':model,'status':status,
                         'extra_HHSS_channels':len(hhss),'additional_quartics':len(expected),
                         'EFT1_seeds_examined':len(evidence),
                         'experimental_UV_inputform_exports_examined':len(experimental),
                         'named_extra_HHSS_in_EFT1_seeds':found,
                         'live_extra_operators_verified':parsed['live_operator_count'] if parsed else 0,
                         'live_extra_CGs_verified':parsed['live_CG_count'] if parsed else 0,
                         'live_export_status':live_result['status'],
                         'physical_HHSS_maps_certified':0,'physical_status':'blocked'})
        except Exception as exc:
            error={'model':model,'status':'V6_source_validation_failed',
                   'error':f'{type(exc).__name__}: {exc}','physical_status':'blocked'}
            _write(model_dir/'failure.json',error)
            rows.append(error)
    completed=all(r['status']!='V6_source_validation_failed' for r in rows)
    counts=sum(r.get('additional_quartics',0) for r in rows)
    hhss_count=sum(r.get('extra_HHSS_channels',0) for r in rows)
    all16=completed and len(chosen)==16 and counts==50 and hhss_count==24
    all16_live=bool(all16 and sum(r.get('live_extra_operators_verified',0) for r in rows)==50)
    report={'format':FORMAT,'models_processed':len(rows),'all16_source_inventory_passed':all16,
        'all16_live_Matchete_extra_operator_exports_passed':all16_live,
        'additional_quartic_operator_inventory':counts,'additional_HHSS_operator_inventory':hhss_count,
        'live_native_extra_Matchete_operators_verified':sum(r.get('live_extra_operators_verified',0) for r in rows),
        'live_native_extra_Matchete_CGs_verified':sum(r.get('live_extra_CGs_verified',0) for r in rows),
        'native_named_extra_HHSS_in_EFT1_seeds':sum(r.get('named_extra_HHSS_in_EFT1_seeds',0) for r in rows),
        'physically_certified_HHSS_coupling_maps':0,'physical_C5_calculated':False,
        'source_GitHub_snapshot_not_a_runtime_physics_certificate':True,
        'two_matching_thresholds_only':True,
        'previous_v5_run':previous_v5,'previous_v5_sha256':file_sha256(root/OUTPUT_OWNER/previous_v5/'summary.json'),
        'generator_contract_sha256':file_sha256(out/'source_generator_contract.json'),
        'models':rows}
    _write(out/'summary.json',report)
    return report
