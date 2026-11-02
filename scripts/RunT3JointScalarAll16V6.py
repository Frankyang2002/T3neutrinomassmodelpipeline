"""All-16 conditional T3 LLSS numerical benchmark, including mixed/split scalars.

Rough exploratory 1 TeV scalar endpoint. Native quartic alias, GammaCC
normalization, spherical CG/component and muF/muS special-field reidentification
are NOT certified. The result is NOT full one-loop physical C5 or neutrino masses.
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

from Numerical.joint_scalar_matching.adapter import ALL_MODELS, normalize_model
from Numerical.joint_scalar_matching.fermion_boundary_bridge import (
    load_release_manifest, load_source, slug)
from Numerical.joint_scalar_matching.source_quartic_v4 import source_seed_path
from Numerical.joint_scalar_matching.numerical_benchmark_v5 import (
    REGULAR,SPECIAL,load_running_snapshot,read_v4_cached_channels,
    build_ordinary_background,evaluate_ordinary_response,inspect_model)
from Numerical.joint_scalar_matching.special_mass_v6 import (
    build_special_background,kind_for_model)
from Numerical.joint_scalar_matching.large_batch import (
    load_background,validate_pivot_kernel,parse_matrix)
from Numerical.joint_scalar_matching.source_su2_batch import (
    validate_su2_response,independent_source_beta_check)

SAFE=ROOT/'output/general_observable/lo_joint_scalar_matching'
RELEASE=ROOT/'output/general_observable/lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f'
V4=SAFE/'source_quartic_v4/native_quartic_channels'
SEEDS=ROOT/'output/full16_symbolic'
RUNNING=ROOT/'output/general_observable/all16_v1'
PIVOTS=ROOT/'output/general_observable/lo_eft1_mixing/source_kernel_v1'
SU2=ROOT/'output/general_observable/lo_eft1_mixing/source_su2_response_v1'
FORMAT='t3-exploratory-all16-conditional-c5-benchmark-v6'
RESULT_FORMAT='t3-conditional-c5-per-unknown-gammacc-mixed-split-v6'
VEV_GEV=246.22


def numerical_metrics(response:dict)->dict:
    """Magnitude diagnostics only; NO assumption N_GammaCC = 1 physically.

    SVD singular values of complex symmetric C5 give conditional Takagi
    magnitudes. v^2/2 is the explicitly assumed Weinberg-vev convention.
    """
    o=response['response']
    total=parse_matrix(o['C5_LLSS_total_one_loop_per_N_gev_inverse'])
    leading=parse_matrix(o['C5_EFT1_induced_leading_log_per_N_gev_inverse'])
    finite=parse_matrix(o['C5_joint_scalar_finite_bubble_per_N_gev_inverse'])
    e=np.linalg.svd(total,compute_uv=False)
    ev=(0.5*VEV_GEV**2*1e9*e).tolist()
    return {'conditional_C5_frobenius_per_N_gev_inverse':float(np.linalg.norm(total)),
            'log_C5_frobenius_per_N_gev_inverse':float(np.linalg.norm(leading)),
            'finite_bubble_C5_frobenius_per_N_gev_inverse':float(np.linalg.norm(finite)),
            'conditional_C5_singular_values_per_N_gev_inverse':e.tolist(),
            'conditional_neutrino_mass_singular_values_eV_per_abs_N':ev,
            'neutrino_mass_vev_gev':VEV_GEV,
            'assumed_mass_convention':'mnu = v^2 C5 / 2; v=246.22 GeV; singular values are Takagi magnitudes',
            'absolute_GammaCC_factor_known':False,
            'complete_physical_neutrino_masses':False,
            'projection_residual':max(o['weinberg_beta_projection_residual'],
                                      o['weinberg_bubble_projection_residual']),
            'matching_scale_cancellation_residual':o['relative_matching_scale_cancellation_residual']}


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--model',help='One model; default all 16')
    p.add_argument('--release-dir',type=Path,default=RELEASE)
    p.add_argument('--v4-channel-dir',type=Path,default=V4)
    p.add_argument('--source-quartic-dir',type=Path,default=SEEDS)
    p.add_argument('--running-dir',type=Path,default=RUNNING)
    p.add_argument('--pivot-dir',type=Path,default=PIVOTS)
    p.add_argument('--su2-dir',type=Path,default=SU2)
    p.add_argument('--output-dir',type=Path,default=SAFE/'exploratory_all16_v6')
    p.add_argument('--allow-exploratory-map',action='store_true',
                   help='Opt into UNVERIFIED Matchete/RGBeta quartic aliases')
    p.add_argument('--allow-exploratory-special-bases',action='store_true',
                   help='Opt into UNVERIFIED native spherical CG and muF-to-muS eigenfield alignment')
    args=p.parse_args(argv)
    output=args.output_dir.resolve(); safe=SAFE.resolve()
    if safe!=output and safe not in output.parents:
        p.error(f'write destination must be under {safe}')
    if (output/'summary.json').exists():
        p.error('existing v6 summary will not be overwritten; choose a NEW output directory')
    if args.allow_exploratory_special_bases and not args.allow_exploratory_map:
        p.error('special basis approximation requires --allow-exploratory-map')
    models=[normalize_model(args.model)] if args.model else list(ALL_MODELS)
    try:
        paths=load_release_manifest(args.release_dir)
    except (ValueError,OSError,TypeError,KeyError) as e:
        p.error(f'fermion boundary release validation: {e}')
    output.mkdir(parents=True,exist_ok=True)
    rows=[]; flat=[]
    for model in models:
        row={'model':model,'status':'blocked','physical_C5_certified':False}
        try:
            boundary=load_source(paths[model],expected_model=model)
            source_seed=source_seed_path(args.source_quartic_dir,model)
            ch=read_v4_cached_channels(args.v4_channel_dir/f'{slug(model)}.json',
                                       source_seed,boundary)
            snapshot=load_running_snapshot(args.running_dir/slug(model)/'running.json',model)
            row.update(inspect_model(ch,boundary,snapshot))
            pivot=validate_pivot_kernel(args.pivot_dir/slug(model)/'quartic_pivot_kernel.json',boundary)
            su2=validate_su2_response(args.su2_dir/slug(model)/'full_su2_response.json',boundary)
            row['v2_pivot_verified']=True
            row['v3_full_SU2_verified']=True
            row['v2_pivot_components']=pivot['pivot_per_source_channels']
            if not args.allow_exploratory_map:
                row['status']='awaiting_exploratory_coupling_map_opt_in'
            elif model in SPECIAL and not args.allow_exploratory_special_bases:
                row['status']='awaiting_special_basis_assumption_opt_in'
            else:
                if model in REGULAR:
                    meta=build_ordinary_background(ch,boundary,snapshot,output/'numeric_backgrounds')
                    category='ordinary'
                else:
                    meta=build_special_background(ch,boundary,snapshot,output/'numeric_backgrounds')
                    category=kind_for_model(model)
                bg=load_background(Path(meta['background_json']),boundary)
                parity=independent_source_beta_check(boundary,bg)
                response=evaluate_ordinary_response(boundary,Path(meta['background_json']))
                metrics=numerical_metrics(response)
                response.update({'format':RESULT_FORMAT,'model':model,
                    'exploratory_scalar_case':category,'exploratory_scalar_basis':category!='ordinary',
                    'full_SU2_beta_source_parity':parity,
                    'numerical_metrics':metrics,'physical_C5_certified':False,
                    'native_quartic_alias_mapping_certified':False,
                    'native_special_scalar_basis_map_certified':False,
                    'muF_muS_mass_eigenfield_alignment_certified':False,
                    'complete_one_loop_C5_certified':False})
                result=output/'C5_per_unknown_NGamma_exploratory'/f'{slug(model)}.json'
                result.parent.mkdir(parents=True,exist_ok=True)
                result.write_text(json.dumps(response,indent=2,allow_nan=False)+'\n',encoding='utf-8')
                row.update({'status':'conditional_C5_per_unknown_NGamma_exploratory',
                            'scalar_case':category,'background':meta,
                            'conditional_result':str(result),'numerical_metrics':metrics,
                            'full_su2_beta_source_parity':parity})
                flat.append({'model':model,'scalar_case':category,
                    'C5_frobenius_per_N_GeV_inv':metrics['conditional_C5_frobenius_per_N_gev_inverse'],
                    'C5_log_frobenius_per_N_GeV_inv':metrics['log_C5_frobenius_per_N_gev_inverse'],
                    'C5_finite_frobenius_per_N_GeV_inv':metrics['finite_bubble_C5_frobenius_per_N_gev_inverse'],
                    'mnu_max_eV_per_abs_N':max(metrics['conditional_neutrino_mass_singular_values_eV_per_abs_N']),
                    'mnu_mid_eV_per_abs_N':sorted(metrics['conditional_neutrino_mass_singular_values_eV_per_abs_N'])[1],
                    'mnu_min_eV_per_abs_N':min(metrics['conditional_neutrino_mass_singular_values_eV_per_abs_N']),
                    'weinberg_projection_residual':metrics['projection_residual']})
        except (ValueError,OSError,KeyError,RuntimeError,TypeError,IndexError) as e:
            row.update({'status':'ERROR_unverified_source_or_physics_consistency',
                        'error':f'{type(e).__name__}: {e}'})
        rows.append(row)
        print(f'{model}: {row["status"]}',flush=True)
        if 'error' in row: print('  '+row['error'],flush=True)
    csv_path=output/'all16_conditional_mass_comparison.csv'
    if flat:
        with csv_path.open('w',newline='',encoding='utf-8') as file:
            writer=csv.DictWriter(file,fieldnames=list(flat[0]))
            writer.writeheader();writer.writerows(flat)
    numeric=sum(r['status']=='conditional_C5_per_unknown_NGamma_exploratory' for r in rows)
    errors=sum(r['status'].startswith('ERROR_') for r in rows)
    summary={'format':FORMAT,'status':'All16ConditionalC5Exploratory' if numeric==16 and errors==0 else
                       'PartialConditionalC5Exploratory',
             'models_checked':len(rows),'models_with_conditional_C5':numeric,
             'special_models_conditional_C5':sum(r['status']=='conditional_C5_per_unknown_NGamma_exploratory' and r.get('scalar_case')!='ordinary' for r in rows),
             'errors':errors,'physical_C5_certified':False,
             'full_one_loop_neutrino_mass_calculation':False,
             'absolute_GammaCC_to_Weyl_normalization':None,
             'allow_exploratory_map':args.allow_exploratory_map,
             'allow_exploratory_special_bases':args.allow_exploratory_special_bases,
             'conditional_mass_convention':'mnu=v^2 C5/2 with v=246.22 GeV, PER |N_GammaCC|',
             'source_muS':'illustrative EFT1 renormalisable running at 1000 GeV, not physical certification',
             'known_approximations':['Uncertified Matchete-to-RGBeta quartic aliases',
                 'Uncertified native spherical/component CG identification for mixed/split scalar basis',
                 'Source muF mass eigenfield definitions identified with muS mass eigenfields',
                 'Uncertified GammaCC signed factor',
                 'Hard-fermion loop terms and full EFT1 Wilson anomalous dimensions not included',
                 'Unbroken scalar masses; no EWSB-resummed mass shifts'],
             'all_scalars_integrated_at_one_muS':True,
             'sequence':'UV -> integrate all F at muF -> EFT1(S1,S2) -> integrate all scalars at one muS -> SMEFT',
             'comparison_csv':str(csv_path) if flat else None,
             'models':rows}
    (output/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('Summary:',output/'summary.json')
    if flat:print('Conditional mass comparison:',csv_path)
    print(f'Conditional C5: {numeric}/{len(models)}; errors={errors}')
    return 2 if errors else 0

if __name__=='__main__':
    raise SystemExit(main())
