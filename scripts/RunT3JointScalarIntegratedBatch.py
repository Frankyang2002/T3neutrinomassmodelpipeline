"""All-16 T3: independently check muF source pivots; conditional full LLSS one-loop response.

No external solver; no GitHub writes.  Default audits the 16 real source
boundaries and existing quartic pivot exports, and creates explicit background
input contracts.  --background-dir additionally computes O(hbar) EFT1 running
plus BOTH-scalars-together bubble per unit unknown GammaCC normalization.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))

from Numerical.joint_scalar_matching.adapter import ALL_MODELS,CONVENTION,normalize_model
from Numerical.joint_scalar_matching.fermion_boundary_bridge import (
    digest, load_release_manifest, load_source, slug)
from Numerical.joint_scalar_matching.large_batch import (
    BACKGROUND_FORMAT,BATCH_FORMAT,evaluate_conditional,load_background,validate_pivot_kernel)

SAFE_OUTPUT=(ROOT/'output/general_observable/lo_joint_scalar_matching').resolve()
DEFAULT_RELEASE=ROOT/'output/general_observable/lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f'
DEFAULT_KERNELS=ROOT/'output/general_observable/lo_eft1_mixing/source_kernel_v1'
DEFAULT_OUTPUT=SAFE_OUTPUT/'integrated_batch_v2'


def background_contract(source):
    n=len(source.labels)
    return {'format':BACKGROUND_FORMAT,'model':source.model,
        'source_boundary_sha256':source.hash,
        'source_scalar_real_labels':list(source.labels),
        'renormalization_scheme':'MSbar',
        'tensor_convention':CONVENTION,
        'mu_f_gev':source.mu_f_gev,
        'mu_s_gev':None,
        'reference_scale_gev':None,
        'scalar_mass_origin':None,
        'same_muS_quartic_and_mass':True,
        'no_separate_scalar_threshold':True,
        'tensor_npz':f'{slug(source.model)}.npz',
        'quartic_channel_names':[],
        'provenance':{'source':'REPLACE_WITH_ACTUAL_SCALAR_EFT1_SOURCE',
         'quartic_potential_normalization':'V=lambda_abcd*phi_a*phi_b*phi_c*phi_d/4!'},
        'REQUIRED_NPZ_ARRAYS_WITH_SHAPES':{
          'source_to_matching_real_U':[n,n],
          'lambda_HHSS_muS':[4,4,n,n],
          'mass_squared_gev2_muS':[n,n],
          'lambda_HHSS_channels_muS':['k',4,4,n,n]},
        'STATUS':'INCOMPLETE_CONTRACT_NOT_VALID_INPUT',
        'IMPORTANT':('No mass, quartic, GammaCC sign, scale, or C5 is guessed. '
          'This template cannot be run until replaced by the true scalar EFT1 tensors.')}


def main(argv=None)->int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--release-dir',type=Path,default=DEFAULT_RELEASE)
    p.add_argument('--kernel-dir',type=Path,default=DEFAULT_KERNELS,
                   help='Existing 16 conditional source-quartic pivot kernel files')
    p.add_argument('--background-dir',type=Path,
                   help='Optional model JSON/NPZ containing independently supplied scalar masses/quartics at muS')
    p.add_argument('--output-dir',type=Path,default=DEFAULT_OUTPUT)
    p.add_argument('--model',help='One model (default all sixteen)')
    args=p.parse_args(argv)
    dest=args.output_dir.resolve()
    if dest!=SAFE_OUTPUT and SAFE_OUTPUT not in dest.parents:
        p.error(f'Write only inside {SAFE_OUTPUT}')
    if (dest/'summary.json').exists():
        p.error('Refusing overwrite: choose a fresh --output-dir')
    models=[normalize_model(args.model)] if args.model else list(ALL_MODELS)
    try: paths=load_release_manifest(args.release_dir)
    except (OSError,ValueError,KeyError,TypeError) as err:
        p.error(f'Invalid all-16 source release: {err}')
    dest.mkdir(parents=True,exist_ok=True)
    records=[]
    for model in models:
        record={'model':model,'status':'blocked','physical_C5_certified':False}
        try:
            source=load_source(paths[model],expected_model=model)
            kernel_path=args.kernel_dir/slug(model)/'quartic_pivot_kernel.json'
            k=validate_pivot_kernel(kernel_path,source)
            record.update({'status':'source_and_quartic_kernel_verified',
                           'source_boundary_sha256':source.hash,
                           'mu_f_gev':source.mu_f_gev,
                           'source_matrix_shape':list(source.raw.shape),
                           'source_count':len(source.source_indices),
                           'quartic_pivot_audit':k,
                           'source_GammaCC_normalization_certified':False})
            # Templates are documentary only and are not valid background files.
            target=dest/'scalar_background_contracts'/f'{slug(model)}.json'
            target.parent.mkdir(parents=True,exist_ok=True)
            target.write_text(json.dumps(background_contract(source),indent=2,
                                        allow_nan=False)+'\n',encoding='utf-8')
            record['scalar_background_contract']=str(target)
            if args.background_dir is not None:
                path=args.background_dir/f'{slug(model)}.json'
                bg=load_background(path,source)
                result=evaluate_conditional(source,bg)
                output=dest/'conditional_one_loop_per_gamma'/f'{slug(model)}.json'
                output.parent.mkdir(parents=True,exist_ok=True)
                output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',
                                  encoding='utf-8')
                record.update({'status':'conditional_Ohbar_LLSS_response_per_unknown_GammaCC',
                               'conditional_result':str(output),
                               'mu_s_gev':bg.mu_s,
                               'max_full_Weinberg_projection_residual':max(
                                    result['weinberg_beta_projection_residual'],
                                    result['weinberg_bubble_projection_residual']),
                               'scale_cancellation_residual':
                                    result['relative_matching_scale_cancellation_residual']})
            else:
                record['remaining_blockers']=[
                     'numerical scalar EFT1 quartic tensor at muS not provided',
                     'numerical positive heavy-scalar mass-squared matrix at muS not provided',
                     'absolute signed Matchete GammaCC/Weyl normalization unresolved',
                     'possible hard F and other operator contributions not certified']
        except (OSError,ValueError,KeyError,IndexError,TypeError,RuntimeError) as err:
            record['status']='blocked'
            record['error']=f'{type(err).__name__}: {err}'
        records.append(record)
        print(f'{model}: {record["status"]}',flush=True)
        if 'error' in record:print('  '+record['error'],flush=True)
    n_ok=sum(x['status']=='source_and_quartic_kernel_verified' for x in records)
    n_cond=sum(x['status']=='conditional_Ohbar_LLSS_response_per_unknown_GammaCC' for x in records)
    n_block=sum(x['status']=='blocked' for x in records)
    summary={'format':BATCH_FORMAT,'status':('SourcePivotAuditComplete' if not n_block and args.background_dir is None
                    else 'AllConditionalResponsesComplete' if not n_block else 'ExplicitPerModelFailures'),
        'model_count':len(models),'models_pivot_verified':n_ok+n_cond,
        'conditional_response_models':n_cond,'blocked_models':n_block,
        'source_release_dir':str(args.release_dir.resolve()),
        'source_release_summary_sha256':digest(args.release_dir/'summary.json'),
        'quartic_pivot_dir':str(args.kernel_dir.resolve()),
        'physical_C5_certified':False,'certified_spinor_normalization':False,
        'sequence':'UV -> integrate all F at muF -> EFT1(S1,S2) -> integrate all scalars at one muS -> SMEFT',
        'note':('Source/kernel audits are real repository cross-checks; no quartic/mass inputs are '
            'invented. Optional conditional output is per unknown GammaCC factor and is not '
            'a certified physical C5 or a complete full-theory one-loop result.'),
        'models':records}
    output=dest/'summary.json'
    output.write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('Summary:',output,flush=True)
    return 2 if n_block else 0


if __name__=='__main__':
    raise SystemExit(main())
