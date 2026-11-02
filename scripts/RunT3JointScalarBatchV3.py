"""Single-command T3 all16 source+full-SU2+spinor-evidence readiness audit.

If separately sourced scalar backgrounds at muS are supplied, also computes
conditional leading-log + finite *joint* scalar response, per unknown GammaCC
factor.  NO separate scalar threshold, invented input, or physical C5 claim.
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
    digest,load_release_manifest,load_source,slug)
from Numerical.joint_scalar_matching.large_batch import (
    BATCH_FORMAT,load_background,evaluate_conditional,validate_pivot_kernel)
from Numerical.joint_scalar_matching.source_su2_batch import (
    audit_signed_evidence,independent_source_beta_check,
    validate_su2_manifest,validate_su2_response)

SAFE=ROOT/'output/general_observable/lo_joint_scalar_matching'
RELEASE=ROOT/'output/general_observable/lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f'
PIVOTS=ROOT/'output/general_observable/lo_eft1_mixing/source_kernel_v1'
SU2=ROOT/'output/general_observable/lo_eft1_mixing/source_su2_response_v1'
NATIVE=ROOT/'output/general_observable/lo_fermion_boundary/normalization_probes/native_feynman_a_m4_b669a99ea16d439b8c50e81b16dbd1c8.json'
SIGNED=ROOT/'output/general_observable/lo_fermion_boundary/normalization_probes/signed_zero_momentum_d734e270ed0138fa.json'


def main(argv=None)->int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--release-dir',type=Path,default=RELEASE)
    p.add_argument('--pivot-dir',type=Path,default=PIVOTS)
    p.add_argument('--su2-dir',type=Path,default=SU2)
    p.add_argument('--native-evidence',type=Path,default=NATIVE)
    p.add_argument('--signed-evidence',type=Path,default=SIGNED)
    p.add_argument('--background-dir',type=Path,default=None,
                   help='OPTIONAL source-hashed scalar mass/quartic input JSON+NPZ, at one muS')
    p.add_argument('--model',help='Model name, default full sixteen')
    p.add_argument('--output-dir',type=Path,default=SAFE/'physical_readiness_v3')
    args=p.parse_args(argv)
    dest=args.output_dir.resolve()
    safe=SAFE.resolve()
    if dest!=safe and safe not in dest.parents:
        p.error(f'All writes must be inside {safe}')
    if (dest/'summary.json').exists():
        p.error('Refusing overwrite of summary.json; choose a fresh output-dir')
    models=[normalize_model(args.model)] if args.model else list(ALL_MODELS)
    try:
        paths=load_release_manifest(args.release_dir)
        manifest=validate_su2_manifest(args.su2_dir,args.release_dir,
                                       models=tuple(models))
        normalization=audit_signed_evidence(args.native_evidence,args.signed_evidence,
                                             args.release_dir)
    except (OSError,ValueError,KeyError,TypeError) as ex:
        p.error(f'Invalid source/normalization prerequsite: {ex}')
    dest.mkdir(parents=True,exist_ok=True)
    records=[]
    for model in models:
        row={'model':model,'status':'blocked','physical_C5_certified':False}
        try:
            src=load_source(paths[model],expected_model=model)
            pivot=validate_pivot_kernel(args.pivot_dir/slug(model)/'quartic_pivot_kernel.json',src)
            su2_path=args.su2_dir/slug(model)/'full_su2_response.json'
            su2=validate_su2_response(su2_path,src)
            matching=manifest[model]
            if matching.get('source_count')!=len(src.source_indices) or matching.get('source_boundary_sha256')!=src.hash:
                raise ValueError('full SU2 manifest/model does not match F-boundary source')
            row.update({'status':'source_pivot_and_full_SU2_verified',
                        'fermion_boundary_source_sha256':src.hash,
                        'mu_f_gev':src.mu_f_gev,
                        'source_scalar_real_labels':list(src.labels),
                        'pivot':pivot,'full_su2':su2,
                        'spinor_normalization_certified':False,
                        'background_at_muS_available':False})
            if args.background_dir is not None:
                bg=load_background(args.background_dir/f'{slug(model)}.json',src)
                parity=independent_source_beta_check(src,bg)
                response=evaluate_conditional(src,bg)
                outfile=dest/'conditional_C5_per_GammaCC'/f'{slug(model)}.json'
                outfile.parent.mkdir(parents=True,exist_ok=True)
                result={'model':model,'interpretation':'conditional one-loop LLSS contribution per unit unknown signed GammaCC normalization; not a physical C5',
                        'full_SU2_independent_beta_parity':parity,
                        'conditional_LLN_scalar_response':response,
                        'physical_C5_certified':False}
                outfile.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
                row.update({'status':'conditional_response_at_muS_per_unknown_GammaCC',
                            'conditional_result_path':str(outfile),
                            'background_at_muS_available':True,'mu_s_gev':bg.mu_s,
                            'independent_full_SU2_beta_parity':parity})
            else:
                row['missing_scalar_background_at_muS']={
                    'required_metadata':f'{slug(model)}.json',
                    'required_npz':f'{slug(model)}.npz',
                    'required_content':['source-to-matching real orthogonal basis transform',
                       'positive all-heavy-scalar mass-squared matrix at muS',
                       'full HHSS quartic tensor and independently labelled quartic channels at muS'],
                    'previous_v2_contract':str(SAFE/'integrated_batch_v2/scalar_background_contracts'/f'{slug(model)}.json')}
            row['remaining_physics_blockers']=[
                'signed GammaCC/Weyl source normalization not certified',
                'full renormalized EFT1 same-scale evolution/finite one-loop matching not certified',
                'hard one-loop heavy-F matching and possible additional dim5 operator sectors not certified']
        except (OSError,ValueError,KeyError,TypeError,RuntimeError,IndexError) as ex:
            row['status']='blocked'
            row['error']=f'{type(ex).__name__}: {ex}'
        records.append(row)
        print(f'{model}: {row["status"]}',flush=True)
        if 'error' in row:print('  '+row['error'],flush=True)
    blocked=sum(v['status']=='blocked' for v in records)
    conditional=sum(v['status']=='conditional_response_at_muS_per_unknown_GammaCC' for v in records)
    summary={'format':'t3-joint-scalar-physical-readiness-batch-v3',
         'status':'ExplicitPerModelFailures' if blocked else
                  ('ConditionalSourceResponsesComplete' if args.background_dir is not None else
                   'AllSourcePivotsAndFullSU2Verified_PhysicalInputsPending'),
         'models_checked':len(records),'models_source_pivot_su2_verified':len(records)-blocked,
         'conditional_responses_computed':conditional,'blocked_models':blocked,
         'source_release_summary_sha256':digest(args.release_dir/'summary.json'),
         'signed_normalization_evidence':normalization,
         'normalization_certified':False,'physical_C5_certified':False,
         'notes':('Source resolved F-tree and full-SU2 parity at muF are formal. '
                  'Optional LLSS-only leading log plus finite joint scalar bubble is conditional, '
                  'not the full one-loop SMEFT Wilson coefficient.'),
         'sequence':'UV -> integrate out all F together -> EFT1(S1,S2) -> integrate out all scalars together at single muS -> SMEFT',
         'models':records}
    out=dest/'summary.json'
    out.write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('Summary:',out,flush=True)
    return 2 if blocked else 0


if __name__=='__main__':
    raise SystemExit(main())
