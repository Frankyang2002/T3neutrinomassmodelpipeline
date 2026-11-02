"""Workstream3: all-16 exact baseline portal completion + extra HHSS unit responses.

READ ONLY upstream WS2 V3/V4, fermion boundary and saved scalar_hard_01.
WRITES only under output/general_observable/lo_certification_completeness/<new_run>.

The 14 missing native baseline portal directions are evaluated with the
conditional WS2 numerical coefficients (two-scale original scalar frame).
The 24 additional named HHSS couplings have no certified physical mapping:
only independent unit-real/unit-imag reference sensitivities are evaluated.
Do not combine these sensitivities with actual RGBeta coefficients.

UV -> ALL heavy fermions at ONE muF -> EFT1(H,S1,S2,SM fermions)
   -> BOTH scalar multiplets at ONE muS -> SMEFT.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))

from Numerical.lo_certification.completeness.inventory import ALL_MODELS
from Numerical.lo_certification.completeness.source_bridge import (
    digest,read_json,split_formal_by_flavour,prepare_conditional_input)
from Numerical.lo_certification.completeness.scalar_bridge import read_scalar_background
from Numerical.lo_certification.completeness.portal_completion_bridge import (
    read_ws2_completed_portals,original_to_boundary,merge_completed_baseline,
    parse_extra_reference,crosscheck_v4_projections,response,extra_unit_response,
    V3_FORMAT,V4_FORMAT)

OWNED=Path('output/general_observable/lo_certification_completeness')
RELEASE=Path('output/general_observable/lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f')
RUNNING=Path('output/general_observable/all16_v1')
OLD=OWNED/'scalar_hard_01'
WS2_V3=Path('output/general_observable/lo_certification_couplings/ws2_quartic_completion_v1')
WS2_V4=Path('output/general_observable/lo_certification_couplings/ws2_named_hhss_v1')
FORMAT='t3-ws3-conditional-portal-completion-and-reference-response-v1'
SPECIAL=frozenset(('A_-2','A_+0','B_-1','C_-1','D_-2','D_+0','E_-2','E_+0'))


def slug(model:str)->str:
    return model.replace('-', 'm').replace('+','p')


def resolved(p:Path)->Path:
    return (p if p.is_absolute() else ROOT/p).resolve()


def allowed_out(p:Path)->Path:
    path=resolved(p); base=resolved(OWNED)
    if path==base or base not in path.parents:
        raise ValueError(f'Refusing output outside Workstream3: {path}')
    return path


def read_npz(p:Path)->dict[str,np.ndarray]:
    with np.load(p,allow_pickle=False) as obj:
        return {k:np.asarray(obj[k]) for k in obj.files}


def json_new(p:Path,obj:dict)->None:
    p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x',encoding='utf-8') as f:
        json.dump(obj,f,indent=2,ensure_ascii=False,allow_nan=False)
        f.write('\n')


def source_inventory(v3:Path,v4:Path,models:tuple[str,...])->dict:
    a=read_json(v3/'summary.json')
    b=read_json(v4/'summary.json')
    if not (a.get('format')==V3_FORMAT and a.get('all16_model_audits_passed') is True
            and a.get('models_processed')==16
            and a.get('baseline_portal_directions_completed')==14
            and a.get('extra_HHSS_registry_entries')==24
            and a.get('physical_C5_calculated') is False):
        raise ValueError('WS2 V3 all-16 exact portal/reference release is missing or differs')
    if not (b.get('format')==V4_FORMAT and
            b.get('all16_source_projections_passed') is True and
            b.get('models_processed')==16 and
            b.get('extra_HHSS_source_symbols_examined')==24 and
            b.get('physical_named_maps_certified')==0 and
            b.get('physical_C5_calculated') is False):
        raise ValueError('WS2 V4 all-16 named-invariant evidence is not conditional/complete')
    if set(x['model'] for x in a.get('models',[]))!=set(ALL_MODELS) or \
       set(x['model'] for x in b.get('models',[]))!=set(ALL_MODELS):
        raise ValueError('WS2 V3/V4 model inventory has changed')
    return {'v3_sha256':digest(v3/'summary.json'),'v4_sha256':digest(v4/'summary.json'),
            'v3_baseline_directions_added':14,'v3_extra_named_registry_symbols':24,
            'v4_physical_named_maps_certified':0}


def model_input(model:str,release:dict[str,Path],running_root:Path,
                old_root:Path,v3:Path,v4:Path,allow_special:bool):
    from Numerical.joint_scalar_matching.fermion_boundary_bridge import load_source
    from Numerical.joint_scalar_matching.numerical_benchmark_v5 import load_running_snapshot
    from Numerical.joint_scalar_matching.special_mass_v6 import build_special_basis
    name=slug(model)
    source=load_source(release[model],expected_model=model)
    running_path=running_root/name/'running.json'
    run=read_json(running_path)
    formal=split_formal_by_flavour(source,run,digest(running_path))
    old=read_scalar_background(old_root/'scalar_backgrounds'/(name+'.json'),source)
    status_path=v3/name/'portal_and_extra_status.json'
    status=read_json(status_path)
    require_path=status.get('baseline_portal_npz')
    if require_path!='baseline_native_plus_rgbeta_portals_muF_muS.npz':
        raise ValueError(f'{model}: V3 baseline portal NPZ name unexpectedly changed')
    array_path=status_path.parent/require_path
    portal=read_ws2_completed_portals(status,read_npz(array_path),model=model,
          source_hash=source.hash,running_hash=digest(running_path),
          n=len(source.labels),native_names=tuple(x for sp in status['baseline_portal_completion']['species']
                                    for x in sp['native_symbols']))
    # native_names requires ONLY old source channels belonging to the baseline:
    # the parser above checks exact set equality, hence pass native subset.
    if model in SPECIAL:
        snap=load_running_snapshot(running_path,model)
        special=build_special_basis(model,snap,source)
        R=original_to_boundary(model,len(source.labels),
                R=special['R_original_from_boundary'],allow_exploratory_special=allow_special)
        if not np.allclose(special['mass2_boundary'],old.mass2_source,rtol=1e-10,atol=1e-7):
            raise ValueError(f'{model}: existing v6 scalar mass inconsistent with source-basis rotation')
    else:
        R=original_to_boundary(model,len(source.labels),R=None,allow_exploratory_special=False)
    labels,chans,change=merge_completed_baseline(tuple(old.channel_names),
                       old.lambda_channels_source,portal,R)
    extras_path=v3/name/'exact_extra_HHSS_reference_tensors.json'
    extras=read_json(extras_path)
    refs=parse_extra_reference(extras,model=model,n=len(source.labels))
    v4path=v4/name/'named_HHSS_source_projection.json'
    named=read_json(v4path)
    ev=named.get('source_evidence',{})
    if ev.get('v3_reference_sha256')!=digest(extras_path) or \
       ev.get('v3_status_sha256')!=digest(status_path):
        raise ValueError(f'{model}: V4 reference/source provenance SHA drift')
    v4inventory=crosscheck_v4_projections(named,refs,model=model)
    summary={'model':model,'status':'SOURCE_TENSORS_VALIDATED_CONDITIONAL',
       'reference_sha256':digest(extras_path),'v3_status_sha256':digest(status_path),
       'v3_npz_sha256':digest(array_path),'v4_projection_sha256':digest(v4path),
       'old_scalar_background_sha256':old.metadata_sha256,
       'old_scalar_npz_sha256':old.arrays_sha256,
       'fermion_boundary_sha256':source.hash,'running_sha256':digest(running_path),
       'n_scalar_real_components':len(old.labels),'n_heavy_fermions':len(formal.fermion_masses_gev),
       'muF_gev':formal.mu_f_gev,'muS_gev':old.mu_s_gev,
       'source_original_to_boundary_basis':'exploratory_V6' if model in SPECIAL else 'identity',
       'baseline_completion':change,
       'two_scale_portal_tensors_sha256':digest(array_path),
       'extra_reference_basis_count':len(refs),
       'extra_named_inventory':v4inventory,
       'missing_physical_inputs':status.get('missing_physical_inputs',[]),
       'physical_C5_certified':False,'all_F_together_at_muF':True,
       'both_S_together_at_muS':True}
    return formal,old,R,portal,labels,chans,refs,summary


def main(argv=None)->int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,default=OWNED/'portal_completion_01')
    parser.add_argument('--old-scalar-hard-dir',type=Path,default=OLD)
    parser.add_argument('--ws2-v3-dir',type=Path,default=WS2_V3)
    parser.add_argument('--ws2-v4-dir',type=Path,default=WS2_V4)
    parser.add_argument('--release-dir',type=Path,default=RELEASE)
    parser.add_argument('--running-dir',type=Path,default=RUNNING)
    parser.add_argument('--model',choices=ALL_MODELS,default=None)
    parser.add_argument('--allow-exploratory-portal-embedding',action='store_true',
                        help='Approve uncertified Matchete vs RGBeta portal component map')
    parser.add_argument('--allow-exploratory-special-bases',action='store_true',
                        help='Approve unverified original-to-split/mixed scalar map')
    parser.add_argument('--skip-extra-unit-responses',action='store_true',
                        help='Only audit V3/V4 extra inventories, no per-unit reference sensitivities')
    parser.add_argument('--check-only',action='store_true')
    args=parser.parse_args(argv)
    target=allowed_out(args.output_dir)
    if target.exists() and not args.check_only:
        parser.error('Output exists; choose a new directory, no overwrite allowed')
    v3,v4=resolved(args.ws2_v3_dir),resolved(args.ws2_v4_dir)
    old,runroot=resolved(args.old_scalar_hard_dir),resolved(args.running_dir)
    inventory=source_inventory(v3,v4,ALL_MODELS)
    from Numerical.joint_scalar_matching.fermion_boundary_bridge import load_release_manifest
    release=load_release_manifest(resolved(args.release_dir))
    models=(args.model,) if args.model else tuple(ALL_MODELS)
    if not args.check_only and not args.allow_exploratory_portal_embedding:
        parser.error('Numerical baseline portal embedding is uncertified; explicitly opt in')
    if not args.check_only and any(m in SPECIAL for m in models) and not args.allow_exploratory_special_bases:
        parser.error('Special models require explicit exploratory scalar basis opt-in')
    rows=[]; fail=0; total_added=0; total_named=0
    if not args.check_only:target.mkdir(parents=True,exist_ok=False)
    for model in models:
        try:
            formal,old,R,portal,names,channels,refs,record=model_input(model,release,runroot,
                old,v3,v4,args.allow_exploratory_special_bases or args.check_only)
            total_added+=record['baseline_completion']['added_baseline_directions']
            total_named+=record['extra_named_inventory']['extra_named_symbols']
            if not args.check_only:
                computed=response(formal,old.mass2_source,old.mu_s_gev,names,channels)
                # Two recorded WS2 endpoint portal tensors are both inspected.
                # This measures the benchmark's portal coupling evolution;
                # it is NOT a substitute for the actual EFT1 seven-term RGE.
                from Numerical.lo_certification.completeness.portal_completion_bridge import rotate_hhss
                f_portal=rotate_hhss(portal.mu_f_channels_original,R)
                s_portal=rotate_hhss(portal.mu_s_channels_original,R)
                f_kernel=response(formal,old.mass2_source,old.mu_s_gev,portal.names,f_portal)
                s_kernel=response(formal,old.mass2_source,old.mu_s_gev,portal.names,s_portal)
                record['two_scale_baseline_portal_diagnostic']={
                    'muF_portal_only_ward_passed':f_kernel['projection_passed'],
                    'muS_portal_only_ward_passed':s_kernel['projection_passed'],
                    'muF_portal_only_hard_frobenius':f_kernel['hard_norm'],
                    'muS_portal_only_hard_frobenius':s_kernel['hard_norm'],
                    'muF_vs_muS_portal_tensor_relative_shift':float(
                        np.linalg.norm(f_portal.sum(axis=0)-s_portal.sum(axis=0))/
                        max(float(np.linalg.norm(s_portal.sum(axis=0))),1e-300)),
                    'muF_vs_muS_portal_hard_matrix_shift':float(np.linalg.norm(
                        np.asarray(f_kernel['matrices']['hard_muF']['real']) +
                        1j*np.asarray(f_kernel['matrices']['hard_muF']['imag']) -
                        np.asarray(s_kernel['matrices']['hard_muF']['real']) -
                        1j*np.asarray(s_kernel['matrices']['hard_muF']['imag']))),
                    'fixed_input_comparison_not_EFT1_RGE':True,
                    'all_nonportal_couplings_at_muF_available':False,
                }
                record['portal_ward_projection_passed']=computed['projection_passed']
                record['portal_max_projection_residual']=max(computed['residuals'].values())
                record['completed_hard_frobenius']=computed['hard_norm']
                record['completed_dim5_frobenius']=computed['dim5_norm']
                record['conditional_sensitivity_only']=True
                record['physical_C5_certified']=False
                newname=slug(model)
                # Write tensor bundle even if the new nonnative channel fails
                # Ward projection; it must not be promoted as a physical C5.
                npz=target/'completed_tensors'/(newname+'.npz')
                npz.parent.mkdir(parents=True,exist_ok=True)
                np.savez_compressed(npz,lambda_HHSS_channels_muS=channels,
                    lambda_HHSS_completed_muS=channels.sum(axis=0),
                    lambda_HHSS_baseline_muF_original=portal.mu_f_channels_original,
                    lambda_HHSS_baseline_muS_original=portal.mu_s_channels_original,
                    original_from_boundary_R=R,mass_squared_muS_boundary=old.mass2_source,
                    fermion_masses_gev=formal.fermion_masses_gev)
                record['completed_tensor_file']=str(npz.relative_to(target))
                if not computed['projection_passed']:
                    computed['status']='WARD_REJECTED_NOT_A_VALID_WEINBERG_COEFFICIENT'
                    record['status']='WARD_REJECTED_NO_VALID_COMPLETED_C5'
                else:
                    computed['status']='CONDITIONAL_COMPLETED_BASELINE_PORTAL'
                    record['status']='COMPLETED_BASELINE_CONDITIONAL_C5'
                json_new(target/'conditional_portal_response'/(newname+'.json'),
                         {'model':model,'status':computed['status'],
                          'baseline_completion':record['baseline_completion'],
                          'response':computed,'physical_C5_certified':False,
                          'source_normalization_certified':False,
                          'all_independent_HHSS_included':False})
                if not args.skip_extra_unit_responses:
                    units=extra_unit_response(formal,old.mass2_source,old.mu_s_gev,refs,R)
                    json_new(target/'extra_reference_unit_responses'/(newname+'.json'),
                        {'model':model,'status':'REFERENCE_UNIT_SENSITIVITIES_NOT_PHYSICAL',
                         'named_Matchete_RGBeta_maps_certified':0,
                         'physical_coefficients_inserted':False,
                         'linear_response_basis':units,'physical_C5_certified':False})
                    record['unit_reference_responses_computed']=len(units)
                else:record['unit_reference_responses_computed']=0
            else:record['status']='SOURCE_INVENTORY_VERIFIED_CHECK_ONLY'
            print(f'{model}: {record["status"]}; baseline+{record["baseline_completion"]["added_baseline_directions"]}, '
                  f'extra named={record["extra_named_inventory"]["extra_named_symbols"]}',flush=True)
        except (OSError,ValueError,KeyError,TypeError,ArithmeticError,AssertionError) as exc:
            record={'model':model,'status':'BLOCKED_SOURCE_BASIS_OR_PROJECTION',
                    'error':f'{type(exc).__name__}: {exc}', 'physical_C5_certified':False}
            fail+=1
            print(f'{model}: {record["status"]}: {record["error"]}',flush=True)
        rows.append(record)
    result={'format':FORMAT,'model_count':len(models),
       'validated_input_models':len(models)-fail,
       'blocked_source_models':fail,'baseline_directions_verified':total_added,
       'extra_named_symbols_verified':total_named,
       'conditional_baseline_ward_passed':sum(r.get('portal_ward_projection_passed') is True for r in rows),
       'conditional_baseline_ward_rejected':sum(r.get('portal_ward_projection_passed') is False for r in rows),
       'new_reference_unit_response_count':sum(r.get('unit_reference_responses_computed',0) for r in rows),
       'ws2_releases':inventory,'models':rows,
       'physical_C5_certified':False,'complete_one_loop_C5_certified':False,
       'physical_named_HHSS_coupling_maps_certified':0,
       'source_GammaCC_to_Weyl_signed_factor':None,
       'source_original_scalar_component_map_certified':False,
       'extra_reference_coefficients_assigned':False,
       'threshold_sequence':'UV -> ALL BSM fermions @ ONE muF -> EFT1(H,S1,S2,SM fermions) -> BOTH scalars @ ONE muS -> SMEFT',
       'fixed_input_hard_portal_only_not_full_eft1_rge':True}
    if not args.check_only:
        json_new(target/'summary.json',result)
        with (target/'numerical_table.csv').open('x',newline='',encoding='utf-8') as f:
            fields=('model','status','portal_ward_projection_passed','portal_max_projection_residual',
                    'completed_hard_frobenius','completed_dim5_frobenius',
                    'unit_reference_responses_computed')
            w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
            for r in rows:w.writerow({k:r.get(k) for k in fields})
    compact={k:v for k,v in result.items() if k in (
        'model_count','validated_input_models','blocked_source_models','baseline_directions_verified',
        'extra_named_symbols_verified','conditional_baseline_ward_passed',
        'conditional_baseline_ward_rejected','new_reference_unit_response_count',
        'physical_C5_certified')}
    print(json.dumps({'output':str(target) if not args.check_only else None,**compact},indent=2))
    if not args.model and fail==0 and (total_added!=14 or total_named!=24):
        print('INVENTORY MISMATCH: expected 14 baseline directions and 24 extra named symbols',file=sys.stderr)
        return 2
    return 0 if fail==0 else 2

if __name__=='__main__':
    raise SystemExit(main())
