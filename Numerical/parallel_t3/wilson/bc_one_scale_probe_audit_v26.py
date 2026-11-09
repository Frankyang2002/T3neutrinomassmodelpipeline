"""v26: inspect previously executed Matchete one-scale all-heavy matching probes.

This DOES NOT run Match[], normalize C5, or promote a simultaneous one-step
probe into the required fermion-first EFT sequence.
"""
from __future__ import annotations
from pathlib import Path
from .bc_native_mass_eigen_reconcile_v25 import (
    BC, GateError, load, digest, require, cp, NATIVE_FORMAT, NATIVE_STATUS)

FORMAT='t3-parallel-a-bc-exploratory-one-scale-source-audit-v26'
STATUS='ExistingMatcheteOneScaleAllHeavyWeinbergProbesAuthenticated_NotFermionFirstEFT'
PROBE_STATUS='OneLoopDirectSymbolicWeinbergProbePassedNumericalC5Pending'
CONVENTION=('L_EFT superset (1/2) C5^(ij) O_ij + h.c., with '
 'O_ij = (L_i^T C eps H)(H^T eps L_j) and symmetric C5; '
 'the ordered Matchete coefficient A5 is mapped by C5^(ij) = A5^(ij) + A5^(ji).')

def audit_model(root:Path,model:str,v25model:dict)->dict:
    stem=model[0]+'_m1_final'
    rootdir=root/'output/general_t3_mass_eigen_matchete'
    cert_path=rootdir/(stem+'_matchete.json')
    probe_dir=rootdir/'one_loop_probe'
    path=probe_dir/(stem+'_one_loop_probe.json')
    cert=load(cert_path);p=load(path)
    require(cert.get('status')==NATIVE_STATUS and cert.get('format')==NATIVE_FORMAT,model+': prior mass certificate changed')
    require(p.get('format')==NATIVE_FORMAT and p.get('status')==PROBE_STATUS and
            p.get('model')==model,model+': native one-loop probe not successful')
    assert_fields=('mass_basis_ordinary_check_passed','heavy_mass_basis_detailed_check_passed',
                   'all_old_scalar_fields_removed','all_17_quartics_present',
                   'both_yukawa_interactions_present')
    for key in assert_fields:
        require(p.get(key) is True and cert.get(key) is True,model+': Matchete basis provenance missing '+key)
    for key in ('matching_to_two_scalar_thresholds_performed','physical_matching_performed',
                'component_basis_map_certified','matchete_rgbeta_quartic_normalization_certified'):
        require(p.get(key) is False and cert.get(key) is False,model+': unsafe matching flag '+key)
    require(p.get('rotation_spurion_values_from_numerical_U')==cert.get('rotation_spurion_values_from_numerical_U'),
            model+': one-loop probe used a different scalar mass basis')
    require(p.get('native_additional_quartic_tensors')==cert.get('native_additional_quartic_tensors') and
            p.get('source_interactions')==cert.get('source_interactions'),
            model+': one-loop probe changed native quartic operator basis')
    require(p.get('one_loop_match_stage_completed') is True and
            p.get('matching_approach')=='SimultaneousAllHeavyOneScaleExploratory' and
            p.get('matching_EFT_order')==5 and p.get('matching_loop_order')==1 and
            p.get('weinberg_extraction_status')=='Success' and
            p.get('weinberg_term_count')==236 and
            p.get('weinberg_holomorphic_term_count')==118 and
            p.get('weinberg_conjugate_term_count')==118 and
            p.get('weinberg_operator_convention')==CONVENTION,
            model+': Weinberg extraction schema changed')
    for key in ('ordered_coefficient_is_canonical_flavour_symmetric_c5',
                'production_sequential_matching_performed','numerical_c5_computed'):
        require(p.get(key) is False,model+': one-loop numerical/physical overclaim '+key)
    files={}
    for key,suffix in (('matched_eft_inputform_file','_one_loop_probe_matched_eft.inputform.txt'),
                       ('weinberg_ordered_coefficient_inputform_file','_one_loop_probe_weinberg_ordered_coefficient.inputform.txt')):
        name=stem+suffix
        require(p.get(key)==name,model+': unexpected inputform artifact name '+key)
        file=probe_dir/name
        require(file.is_file() and 1000<=file.stat().st_size<100_000_000,
                model+': Matchete one-loop source artifact missing/empty '+name)
        files[key]={'file':file.relative_to(root).as_posix(),
                    'byte_count':file.stat().st_size,'sha256':digest(file)}
    require(v25model.get('source_sha256',{}).get('native_matchete_eigenbasis')==digest(cert_path),
            model+': v25 native Matchete certificate changed')
    return {'model':model,'native_exploratory_match_completed':True,
            'weinberg_holomorphic_terms':118,'weinberg_conjugate_terms':118,
            'one_step_all_heavy_not_fermion_first':True,
            'ordered_C5_not_yet_canonical_symmetric':True,
            'actual_EFT1_scalar_threshold_matching_performed':False,
            'source_files':files,
            'source_sha256':{'native_mass_certificate':digest(cert_path),
                             'native_one_scale_probe':digest(path)}}

def audit(root:Path,v25_path:Path)->dict:
    root=Path(root).resolve();v25=load(v25_path)
    from .bc_native_mass_eigen_reconcile_v25 import FORMAT as F25,STATUS as S25
    require(v25.get('format')==F25 and v25.get('status')==S25 and
            v25.get('models_checked')==2 and
            v25.get('certification',{}).get('physical_Wilson_evolution_authorized') is False,
            'v25 native eigenfield bridge not source-conditional')
    require(set(v25.get('models',{}))==set(BC),'v25 missing B/C native model')
    out={m:audit_model(root,m,v25['models'][m]) for m in BC}
    return {'format':FORMAT,'status':STATUS,'models_checked':2,'models':out,
            'one_step_native_Weinberg_holomorphic_total':236,
            'native_one_loop_extraction_run_previously':True,
            'source_sha256':{'v25':digest(v25_path)},
            'certification':{'existing_native_one_step_Weinberg_probe_authenticated':True,
              'physical_fermion_first_matching_performed':False,
              'canonical_flavour_symmetric_C5_recovered':False,
              'one_loop_dim5_ADM':False,
              'physical_Wilson_evolution_authorized':False}}
