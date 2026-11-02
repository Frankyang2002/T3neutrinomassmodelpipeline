"""Use existing 1-TeV running snapshots for conditional LLSS joint matching.

Exactly eight ordinary-basis T3 models can be evaluated with the existing
illustrative, diagonal *unbroken* scalar masses. Eight mixed/split-real models
are deliberately blocked without their muS mass/mixing and basis map.

--allow-exploratory-map opts into unverified Matchete->RGBeta alias assumptions.
Every output says NOT a physical C5 or certified neutrino mass prediction.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))

from Numerical.joint_scalar_matching.adapter import ALL_MODELS,normalize_model
from Numerical.joint_scalar_matching.fermion_boundary_bridge import (
    load_release_manifest,load_source,slug,digest)
from Numerical.joint_scalar_matching.source_quartic_v4 import source_seed_path
from Numerical.joint_scalar_matching.numerical_benchmark_v5 import (
    BENCHMARK_FORMAT,REGULAR,SPECIAL,read_v4_cached_channels,load_running_snapshot,
    inspect_model,build_ordinary_background,evaluate_ordinary_response)
from Numerical.joint_scalar_matching.large_batch import load_background,validate_pivot_kernel
from Numerical.joint_scalar_matching.source_su2_batch import (
    validate_su2_response,independent_source_beta_check)

SAFE=ROOT/'output/general_observable/lo_joint_scalar_matching'
RELEASE=ROOT/'output/general_observable/lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f'
V4=SAFE/'source_quartic_v4/native_quartic_channels'
SEEDS=ROOT/'output/full16_symbolic'
RUNNING=ROOT/'output/general_observable/all16_v1'
PIVOTS=ROOT/'output/general_observable/lo_eft1_mixing/source_kernel_v1'
SU2=ROOT/'output/general_observable/lo_eft1_mixing/source_su2_response_v1'


def main(argv=None)->int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--model',help='Only one model; default all 16')
    p.add_argument('--release-dir',type=Path,default=RELEASE)
    p.add_argument('--v4-channel-dir',type=Path,default=V4)
    p.add_argument('--source-quartic-dir',type=Path,default=SEEDS)
    p.add_argument('--running-dir',type=Path,default=RUNNING)
    p.add_argument('--pivot-dir',type=Path,default=PIVOTS)
    p.add_argument('--su2-dir',type=Path,default=SU2)
    p.add_argument('--output-dir',type=Path,default=SAFE/'exploratory_benchmark_v5')
    p.add_argument('--allow-exploratory-map',action='store_true',
                   help='Compute 8 conditional benchmark C5/N matrices using the EXPLICITLY UNVERIFIED quartic aliases')
    args=p.parse_args(argv)
    dest=args.output_dir.resolve()
    safe=SAFE.resolve()
    if dest!=safe and safe not in dest.parents:
        p.error(f'Outputs restricted to {safe}')
    if (dest/'summary.json').exists():
        p.error('Refusing overwrite of existing v5 summary; choose a new output-dir')
    models=[normalize_model(args.model)] if args.model else list(ALL_MODELS)
    try:
        paths=load_release_manifest(args.release_dir)
    except (OSError,ValueError,KeyError,TypeError) as ex:
        p.error(f'Cannot verify 16-source release: {ex}')
    dest.mkdir(parents=True,exist_ok=True)
    results=[]
    for model in models:
        row={'model':model,'status':'blocked','physical_C5_certified':False}
        try:
            src=load_source(paths[model],expected_model=model)
            qpath=source_seed_path(args.source_quartic_dir,model)
            cached=args.v4_channel_dir/f'{slug(model)}.json'
            ch=read_v4_cached_channels(cached,qpath,src)
            running_path=args.running_dir/slug(model)/'running.json'
            snapshot=load_running_snapshot(running_path,model)
            info=inspect_model(ch,src,snapshot)
            row.update(info)
            # Per-source and full SU2 source checks are preserved from v2/v3.
            pivot=validate_pivot_kernel(args.pivot_dir/slug(model)/'quartic_pivot_kernel.json',src)
            su2=validate_su2_response(args.su2_dir/slug(model)/'full_su2_response.json',src)
            row['verified_source_pivot']={'channels':pivot['pivot_per_source_channels'],
                                           'max_abs_error':pivot['source_vs_kernel_max_abs_error']}
            row['verified_source_su2_status']=su2.get('status','verified')
            if model in SPECIAL:
                row['status']='blocked_special_scalar_basis_at_muS'
            elif not args.allow_exploratory_map:
                row['status']='awaiting_explicit_exploratory_alias_opt_in'
            else:
                output=build_ordinary_background(ch,src,snapshot,dest/'numerical_backgrounds')
                bg=load_background(Path(output['background_json']),src)
                parity=independent_source_beta_check(src,bg)
                response=evaluate_ordinary_response(src,Path(output['background_json']))
                response['independent_full_su2_source_parity']=parity
                response['benchmark_sha256']=snapshot.sha256
                response['aliases_not_certified']=output['coupling_aliases_NOT_CERTIFIED']
                target=dest/'C5_per_unknown_NGamma_exploratory'/f'{slug(model)}.json'
                target.parent.mkdir(parents=True,exist_ok=True)
                target.write_text(json.dumps(response,indent=2,allow_nan=False)+'\n',encoding='utf-8')
                row.update({'status':'conditional_numerical_LLSS_C5_per_NGamma_exploratory',
                            'background':output,'conditional_result':str(target),
                            'full_su2_beta_parity':parity})
        except (OSError,ValueError,KeyError,RuntimeError,TypeError,IndexError) as ex:
            row['status']='ERROR_source_or_trajectory_validation'
            row['error']=f'{type(ex).__name__}: {ex}'
        results.append(row)
        print(f'{model}: {row["status"]}',flush=True)
        if 'error' in row:print('  '+row['error'],flush=True)
    error=sum(x['status']=='ERROR_source_or_trajectory_validation' for x in results)
    numeric=sum(x['status']=='conditional_numerical_LLSS_C5_per_NGamma_exploratory' for x in results)
    special=sum(x['status']=='blocked_special_scalar_basis_at_muS' for x in results)
    summary={'format':BENCHMARK_FORMAT,
             'status':'ERRORS_PRESENT' if error else
                'ExploratoryOrdinaryModelsReady_SpecialsBlocked' if numeric else
                'All16TrajectoriesInspected_NoNumericalCalculationRequested',
             'models_checked':len(models),'conditional_ordinary_C5_per_NGamma_count':numeric,
             'special_scalar_models_correctly_blocked':special,'errors':error,
             'physical_C5_certified':False,
             'benchmark_coupling_map_certified':False,
             'absolute_GammaCC_to_Weyl_normalization':None,
             'actual_benchmark':'output/general_observable/all16_v1/*/running.json; 1 TeV illustrative EFT1 endpoint',
             'flags':{'exploratory_alias_map_opt_in':args.allow_exploratory_map,
                      'source_operator_inventory_complete':False,
                      'hard_one_loop_fermion_matching_complete':False,
                      'EFT1_full_Wilson_running_complete':False,
                      'all_scalars_integrated_at_one_muS':True},
             'sequence':'UV -> all F at muF -> EFT1(S1,S2) -> all scalars at ONE muS -> SMEFT',
             'models':results}
    (dest/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(f'Summary: {dest/"summary.json"}',flush=True)
    print(f'Conditional ordinary models: {numeric}; special basis blocked: {special}; source errors: {error}',flush=True)
    return 2 if error else 0

if __name__=='__main__':
    raise SystemExit(main())
