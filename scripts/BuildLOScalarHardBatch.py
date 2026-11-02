"""All-16 source-backed *exploratory* scalar-input + hard-portal matching batch.

Reuse existing v4 Matchete HHSS channel exports, v5 muS renormalizable state,
v6 mixed/split scalar basis, and Workstream3 per-fermion source decomposition.
Nothing is written outside the Workstream3 output subtree, and existing
upstream source/outputs are read-only. UV->ALL F @ one muF->EFT1(S1,S2)
->BOTH scalars @ one muS->SMEFT. No physical C5 certification.

Run from repository root:
python -m scripts.BuildLOScalarHardBatch --allow-exploratory-map \
    --allow-exploratory-special-bases \
    --output-dir output/general_observable/lo_certification_completeness/scalar_hard_01
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path: sys.path.insert(0,str(ROOT))

from Numerical.lo_certification.completeness.source_bridge import (
    digest,read_json,split_formal_by_flavour,prepare_conditional_input)
from Numerical.lo_certification.completeness.scalar_bridge import (
    read_scalar_background,conditional_scalar_payload,
    build_existing_exploratory_scalar_background)
from Numerical.lo_certification.completeness.inventory import ALL_MODELS
from Numerical.joint_scalar_matching.fermion_boundary_bridge import (
    load_release_manifest,load_source,slug)
from scripts.RunLOCompletenessAudit import compute_conditional

OWNED=Path('output/general_observable/lo_certification_completeness')
RELEASE=Path('output/general_observable/lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f')
RUNNING=Path('output/general_observable/all16_v1')
V4=Path('output/general_observable/lo_joint_scalar_matching/source_quartic_v4/native_quartic_channels')
SEEDS=Path('output/full16_symbolic')
BACKGROUND_FORMAT='t3-lo-scalar-hard-source-bridge-v1'


def source_seed_path(root:Path,model:str) -> Path:
    from Numerical.joint_scalar_matching.source_quartic_v4 import source_seed_path
    return source_seed_path(root,model)


def within_owned(target:Path) -> Path:
    dest=target if target.is_absolute() else ROOT/target
    dest=dest.resolve()
    allowed=(ROOT/OWNED).resolve()
    if dest==allowed or allowed not in dest.parents:
        raise ValueError('Output must be a new descendant of '+str(allowed))
    if dest.exists():
        raise FileExistsError('Output already exists: '+str(dest))
    return dest


def main(argv=None)->int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--release-dir',type=Path,default=RELEASE)
    p.add_argument('--running-dir',type=Path,default=RUNNING)
    p.add_argument('--v4-channel-dir',type=Path,default=V4)
    p.add_argument('--source-quartic-dir',type=Path,default=SEEDS)
    p.add_argument('--existing-background-dir',type=Path,default=None,
                   help='Read-only existing v5/v6 numeric_backgrounds (no regeneration)')
    p.add_argument('--allow-exploratory-map',action='store_true',
                   help='Explicit acceptance of UNVERIFIED Matchete/RGBeta quartic aliases')
    p.add_argument('--allow-exploratory-special-bases',action='store_true',
                   help='Explicit acceptance of UNVERIFIED eight split/mixed scalar maps')
    p.add_argument('--check-only',action='store_true',help='Source inventory check; no files written')
    args=p.parse_args(argv)
    dest=within_owned(args.output_dir)
    if args.allow_exploratory_special_bases and not args.allow_exploratory_map:
        p.error('--allow-exploratory-special-bases requires --allow-exploratory-map')
    if args.existing_background_dir is None and not args.allow_exploratory_map and not args.check_only:
        p.error('Auto-build uses uncertified native coupling aliases; set --allow-exploratory-map')
    resolved=lambda x: x if x.is_absolute() else ROOT/x
    release=load_release_manifest(resolved(args.release_dir))
    selected=list(ALL_MODELS)
    if args.check_only:
        inventory=[]
        for model in selected:
            source=load_source(release[model],expected_model=model)
            runpath=resolved(args.running_dir)/slug(model)/'running.json'
            saved=read_json(runpath)
            formal=split_formal_by_flavour(source,saved,digest(runpath))
            cache=resolved(args.v4_channel_dir)/(slug(model)+'.json')
            seed=source_seed_path(resolved(args.source_quartic_dir),model)
            bg=(resolved(args.existing_background_dir)/(slug(model)+'.json') if args.existing_background_dir else None)
            inventory.append({'model':model,'formal_per_k_reconstructed':True,
                'v4_cache_exists':cache.is_file(),'seed_exists':seed.is_file(),
                'background_exists':bg.is_file() if bg else None,
                'n_source_quartic_evidence_files':int(cache.is_file())+int(seed.is_file()),
                'normalization_certified':False})
        print(json.dumps({'models':inventory,'physical_C5_certified':False},indent=2))
        return 0
    dest.mkdir(parents=True,exist_ok=False)
    (dest/'scalar_backgrounds').mkdir()
    (dest/'conditional_inputs').mkdir()
    (dest/'conditional_hard').mkdir()
    (dest/'projection_rejected').mkdir()
    rows=[]
    for model in selected:
        formal=None
        scalar=None
        record={'model':model,'status':'UNRESOLVED','physical_C5_certified':False,
            'signed_GammaCC_normalization':None,'complete_general_scalar_potential_certified':False,
            'hard_one_loop_matching_fully_certified':False,
            'scalar_case':'special' if model in ('A_-2','A_+0','B_-1','C_-1',
                'D_-2','D_+0','E_-2','E_+0') else 'ordinary'}
        try:
            source=load_source(release[model],expected_model=model)
            running_path=resolved(args.running_dir)/slug(model)/'running.json'
            running=read_json(running_path)
            if running.get('model')!=model:
                raise ValueError(f'{model}: running model mismatch')
            formal=split_formal_by_flavour(source,running,digest(running_path))
            record['source_sha256']=source.hash
            record['running_sha256']=formal.running_hash
            record['n_heavy_flavours']=len(formal.fermion_masses_gev)
            record['n_scalar_real_components']=len(formal.labels)
            record['source_formal_per_k_verified']=True
            if args.existing_background_dir is None:
                background_path=build_existing_exploratory_scalar_background(
                    source,running_path,
                    resolved(args.v4_channel_dir)/(slug(model)+'.json'),
                    source_seed_path(resolved(args.source_quartic_dir),model),
                    dest/'scalar_backgrounds',
                    allow_exploratory_map=args.allow_exploratory_map,
                    allow_exploratory_special_bases=args.allow_exploratory_special_bases)
            else:
                background_path=resolved(args.existing_background_dir)/(slug(model)+'.json')
            scalar=read_scalar_background(background_path,source)
            meta,arr=conditional_scalar_payload(scalar)
            metadata,arrays=prepare_conditional_input(formal,meta,arr,scalar.arrays_sha256)
            metadata.update({
                'source_scalar_background_metadata_sha256':scalar.metadata_sha256,
                'source_scalar_background_npz_sha256':scalar.arrays_sha256,
                'source_scalar_case':record['scalar_case'],
                'native_to_rgbeta_quartic_map_certified':False,
                'source_muF_muS_basis_identification_certified':False,
                'source_channel_inventory_scope':'cached native source HHSS, not full UV closure',
                'physical_C5_certified':False})
            name=slug(model)
            npz=dest/'conditional_inputs'/(name+'.npz')
            json_path=dest/'conditional_inputs'/(name+'.json')
            np.savez_compressed(npz,**arrays)
            json_path.write_text(json.dumps(metadata,indent=2,allow_nan=False)+'\n',encoding='utf-8')
            hard=compute_conditional(model,dest/'conditional_inputs')
            (dest/'conditional_hard'/(name+'.json')).write_text(
                json.dumps(hard,indent=2,allow_nan=False)+'\n',encoding='utf-8')
            record.update({'status':'CONDITIONAL_HARD_PORTAL_PER_UNKNOWN_NGAMMA',
                           'hard_input':str(json_path.relative_to(dest)),
                           'hard_diagnostic':'conditional_hard/'+name+'.json',
                           'all_native_channels_preserved':len(scalar.channel_names),
                           'weinberg_projection_passed':hard['projection_passed'],
                           'hard_computation_successful':True})
        except (ValueError,RuntimeError,KeyError,TypeError,FileNotFoundError,AssertionError) as exc:
            record.update({'status':'BLOCKED_SOURCE_OR_WEINBERG_PROJECTION',
                           'hard_computation_successful':False,
                           'blocker':f'{type(exc).__name__}: {exc}'})
            if (formal is not None and scalar is not None and
                    'projection' in str(exc).lower()):
                # Keep the algebraic result for diagnosis without exporting
                # an invalid Weinberg coefficient or a physical hard input.
                try:
                    from Numerical.lo_certification.completeness.hard_matching import (
                        evaluate_hard_scalar_channel)
                    r=evaluate_hard_scalar_channel(formal.c_by_fermion,
                           formal.fermion_masses_gev,scalar.lambda_source,
                           scalar.mass2_source,formal.mu_f_gev,scalar.mu_s_gev,
                           C_total_muF=formal.c_total,
                           lambda_channels=scalar.lambda_channels_source,
                           channel_names=scalar.channel_names)
                    residuals={k:float(v) for k,v in r['residuals'].items()}
                    record['projection_residuals']=residuals
                    diagnostic={'model':model,'status':'REJECTED_WEINBERG_PROJECTION',
                      'projection_residuals':residuals,'physical_C5_certified':False,
                      'projected_matrices_are_NOT_valid_physical_C5':True,
                      'per_pair_terms':r['per_pair_terms'],
                      'matrices_projected_diagnostic_only':{
                         key:{'real':v.real.tolist(),'imag':v.imag.tolist()}
                         for key,v in r['matrices'].items()}}
                    path=dest/'projection_rejected'/(slug(model)+'.json')
                    path.write_text(json.dumps(diagnostic,indent=2,allow_nan=False)+'\n',encoding='utf-8')
                    record['projection_diagnostic']=str(path.relative_to(dest))
                except (ValueError,RuntimeError,TypeError,AssertionError) as second:
                    record['secondary_projection_diagnostic_error']=f'{type(second).__name__}: {second}'
        rows.append(record)
        print(f"{model}: {record['status']}" +
              (f" [{record['blocker']}]" if 'blocker' in record else ''),flush=True)
    computed=sum(x['hard_computation_successful'] for x in rows)
    blocked=len(rows)-computed
    summary={'format':BACKGROUND_FORMAT,'models':rows,'model_count':16,
        'conditional_hard_models':computed,'blocked_models':blocked,
        'physical_C5_certified':False,'complete_leading_one_loop_C5_certified':False,
        'all_fermions_integrated_at_one_muF':True,'both_scalar_multiplets_integrated_at_one_muS':True,
        'SM_fermions_remain':True,'unverified_quartic_mapping_opt_in':args.allow_exploratory_map,
        'unverified_special_scalar_basis_opt_in':args.allow_exploratory_special_bases,
        'source_snapshot':'read-only original Matchete LLSS release and cached HHSS channels',
        'important_caveat':'Hard scalar-portal topology only, per unknown signed GammaCC/Weyl factor; finite renormalized hard matching and other dim5 structures not established'}
    (dest/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({'output':str(dest),'conditional_hard_models':computed,
                       'blocked_models':blocked,'physical_C5_certified':False},indent=2))
    return 0

if __name__=='__main__':
    raise SystemExit(main())
