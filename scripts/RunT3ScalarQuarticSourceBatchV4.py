"""All-16 Matchete-native quartic tensor extraction, optional muS background build.

UV -> integrate ALL fermions -> EFT1(S1,S2) -> integrate ALL scalars at ONE muS.
No physical GammaCC/Weyl factor or muS scalar masses are invented.
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
from Numerical.joint_scalar_matching.source_quartic_v4 import (
    FORMAT,build_source_channels,source_seed_path,symbolic_export,build_background)

SAFE=ROOT/'output/general_observable/lo_joint_scalar_matching'
RELEASE=ROOT/'output/general_observable/lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f'
SEEDS=ROOT/'output/full16_symbolic'


def main(argv=None)->int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--model',help='One T3 model (default all 16)')
    p.add_argument('--release-dir',type=Path,default=RELEASE)
    p.add_argument('--source-quartic-dir',type=Path,default=SEEDS)
    p.add_argument('--endpoint-dir',type=Path,default=None,
                   help='Optional source-hashed actual scalar muS endpoints, A_m4.json etc')
    p.add_argument('--output-dir',type=Path,default=SAFE/'source_quartic_v4')
    p.add_argument('--run-v3',action='store_true',
                   help='If ALL requested numerical backgrounds succeed, launch existing v3 conditional matcher')
    args=p.parse_args(argv)
    dest=args.output_dir.resolve()
    safe=SAFE.resolve()
    if dest!=safe and safe not in dest.parents:
        p.error(f'All writes restricted to {safe}')
    if (dest/'summary.json').exists():
        p.error('summary.json already exists; choose a NEW output-dir')
    if args.run_v3 and args.endpoint_dir is None:
        p.error('--run-v3 requires --endpoint-dir with physically sourced muS values')
    models=[normalize_model(args.model)] if args.model else list(ALL_MODELS)
    try:
        path_map=load_release_manifest(args.release_dir)
    except (ValueError,OSError,KeyError,TypeError) as ex:
        p.error(f'F-boundary release failed verification: {ex}')
    dest.mkdir(parents=True,exist_ok=True)
    records=[]
    for model in models:
        record={'model':model,'status':'blocked','physical_C5_certified':False}
        try:
            boundary=load_source(path_map[model],expected_model=model)
            path=source_seed_path(args.source_quartic_dir,model)
            source=build_source_channels(path,model)
            out=dest/'native_quartic_channels'/f'{slug(model)}.json'
            out.parent.mkdir(parents=True,exist_ok=True)
            if out.exists():raise FileExistsError(f'refusing overwrite: {out}')
            info=symbolic_export(source,boundary)
            out.write_text(json.dumps(info,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')
            record.update({'status':'source_HHSS_channels_parsed_muS_numerical_input_pending',
                           'native_quartic_channels_path':str(out),
                           'source_seed_sha256':source.seed_sha256,
                           'source_hhss_term_count':source.seed_hhss_term_count,
                           'independent_quartic_couplings':sorted(source.by_name),
                           'quartic_channels':len(source.by_name),
                           'original_scalar_basis_matches_fermion_boundary':
                               list(source.source_basis_labels)==list(boundary.labels),
                           'source_hhss_operator_completeness_not_certified':True})
            if args.endpoint_dir is not None:
                ep=args.endpoint_dir/f'{slug(model)}.json'
                if not ep.is_file():
                    record['status']='source_HHSS_exported__muS_endpoint_missing'
                    record['endpoint_required']=str(ep)
                else:
                    meta=build_background(source,boundary,ep,dest/'numeric_backgrounds')
                    record.update({'status':'muS_source_quartics_and_mass_background_generated_conditional',
                                   'numeric_background':meta,
                                   'GammaCC_to_Weyl_normalization_certified':False})
        except (OSError,RuntimeError,TypeError,ValueError,KeyError,IndexError) as ex:
            record.update({'status':'blocked',
                           'error':f'{type(ex).__name__}: {ex}'})
        records.append(record)
        print(f'{model}: {record["status"]}',flush=True)
        if 'error' in record: print('  '+record['error'],flush=True)
    nparsed=sum(v['status'].startswith(('source_HHSS_channels_parsed','source_HHSS_exported',
                                        'muS_source_quartics')) for v in records)
    nbackground=sum(v['status'].startswith('muS_source_quartics') for v in records)
    nblocked=sum(v['status']=='blocked' or v['status'].endswith('_missing') for v in records)
    summary={'format':'t3-matchete-source-hhss-quartic-batch-v4',
             'status':'SourceQuarticsParsed_NoNumericalMuSEndpoints' if not nblocked and args.endpoint_dir is None
                      else 'NumericalMuSBackgroundsGenerated_ConditionalOnly' if nbackground==len(models)
                      else 'ExplicitPerModelBlockers',
             'models_checked':len(models),'source_quartic_models_parsed':nparsed,
             'numerical_muS_backgrounds_generated':nbackground,'blocked_or_missing':nblocked,
             'fermion_release_summary_sha256':digest(args.release_dir/'summary.json'),
             'source_quartic_root':str(args.source_quartic_dir.resolve()),
             'potential_normalization':'V=lambda_abcd*phi_a*phi_b*phi_c*phi_d/4!',
             'physical_GammaCC_to_Weyl_normalization_certified':False,
             'physical_all16_C5_certified':False,
             'limitations':[
                'Source tensor is from original (S1,S2) Matchete basis, not automatically the B/C X/Y or split-real fermion basis',
                'Full general UV scalar-potential operator inventory is not certified by a source seed alone',
                'Numerical masses and quartic values at muS require an independently sourced common-scale trajectory',
                'One-loop hard F-matching, additional EFT1 operators and finite counterterms remain unverified'],
             'sequence':'UV -> integrate all F at muF -> EFT1(S1,S2) -> integrate all scalars together at one muS -> SMEFT',
             'models':records}
    out=dest/'summary.json'
    out.write_text(json.dumps(summary,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')
    print('Summary:',out,flush=True)
    if args.run_v3:
        if nbackground!=len(models) or nblocked:
            print('v3 conditional matching NOT launched: incomplete muS backgrounds',flush=True)
            return 2
        from scripts.RunT3JointScalarBatchV3 import main as v3_main
        call=['--background-dir',str(dest/'numeric_backgrounds'),
              '--output-dir',str(dest/'conditional_matching_v3')]
        if args.model:call+=['--model',models[0]]
        return v3_main(call)
    return 2 if nblocked else 0


if __name__=='__main__':
    raise SystemExit(main())
