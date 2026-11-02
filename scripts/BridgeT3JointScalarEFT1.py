"""Bridge all-16 read-only fermion-first source release to joint scalar matching.

Without --endpoint-dir, only exports an *uncertified formal muF tensor audit*.
With --endpoint-dir, requires a separately signed/checked normalized muS EFT1
handoff, then runs the EXISTING joint-scalar engine at one scalar threshold.
No fermion threshold, EFT1 evolution, or scalar EFT2 is executed here.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0,str(ROOT))

from Numerical.joint_scalar_matching.adapter import (
    ALL_MODELS, load_input, normalize_model, save_input)
from Numerical.joint_scalar_matching.fermion_boundary_bridge import (
    BRIDGE_FORMAT, load_eft1_endpoint, load_release_manifest, load_source,
    slug, write_provisional_source)
from scripts.RunT3JointScalarMatching import compute

DEFAULT_RELEASE = (ROOT / 'output/general_observable/lo_fermion_boundary/releases/'
                   'source_021dcc1e8fe8ba0f')
SAFE_OUTPUT = (ROOT / 'output/general_observable/lo_joint_scalar_matching').resolve()
DEFAULT_OUTPUT = SAFE_OUTPUT / 'fermion_boundary_bridge_v1'


def main(argv=None) -> int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--release-dir',type=Path,default=DEFAULT_RELEASE,
                   help='Existing read-only fermion-first 16-model source release')
    p.add_argument('--endpoint-dir',type=Path,
                   help='Optional independent full EFT1 muS endpoint JSON/NPZ handoffs')
    p.add_argument('--output-dir',type=Path,default=DEFAULT_OUTPUT)
    p.add_argument('--model',help='Single model (default all 16)')
    args=p.parse_args(argv)
    selected=[normalize_model(args.model)] if args.model else list(ALL_MODELS)
    dest=args.output_dir.resolve()
    if dest != SAFE_OUTPUT and SAFE_OUTPUT not in dest.parents:
        p.error(f'Output must remain inside {SAFE_OUTPUT}')
    if (dest/'summary.json').exists():
        p.error('Output summary exists. Use a new --output-dir to preserve previous results.')
    # Existing release and endpoints are read-only; never write there.
    try:
        source_files=load_release_manifest(args.release_dir)
    except (OSError,ValueError,KeyError,TypeError) as err:
        p.error(f'Invalid source release: {err}')
    dest.mkdir(parents=True,exist_ok=True)
    reports=[]
    for model in selected:
        row={'model':model,'source_file':str(source_files[model]),
             'status':'blocked','physical_full_one_loop_certified':False}
        try:
            source=load_source(source_files[model],expected_model=model)
            row.update({'source_boundary_sha256':source.hash,
                'mu_f_gev':source.mu_f_gev,
                'n_formal_source_PL_contributions':len(source.source_indices),
                'nonzero_formal_LLSS_components':int((abs(source.raw)>0).sum()),
                'source_GammaCC_Weyl_certified':source.normalization_certified})
            formal=dest/'muF_formal_sources'/f'{slug(model)}.npz'
            write_provisional_source(formal,source)
            row['formal_muF_npz']=str(formal)
            if args.endpoint_dir is None:
                row['status']='formal_muF_exported__muS_endpoint_missing__normalization_open'
                row['blockers']=[
                    'absolute signed GammaCC to Weyl normalization not certified in source release',
                    'no independently renormalized full EFT1 tensors at common muS',
                    'independent complete finite scalar matching not certified']
            else:
                endpoint=Path(args.endpoint_dir).resolve()/f'{slug(model)}.json'
                if not endpoint.is_file():
                    raise ValueError(f'missing EFT1 endpoint metadata {endpoint}')
                inp=load_eft1_endpoint(endpoint,source)
                meta=save_input(dest/'normalized_muS_inputs'/f'{slug(model)}.json',inp)
                # Exercise the actual stable public JSON/NPZ adapter, not a parallel math kernel.
                reloaded=load_input(meta)
                result=compute(reloaded)
                result['joint_scalar_bridge']={
                    'source_boundary_sha256':source.hash,
                    'endpoint_metadata_path':str(endpoint),
                    'source_resolved_PL_terms':len(source.source_indices),
                    'independent_normalization_evidence_declared':True,
                    'physical_full_one_loop_matching_certified':False}
                target=dest/'candidate_C5'/f'{slug(model)}_c5.json'
                target.parent.mkdir(parents=True,exist_ok=True)
                target.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
                row.update({'status':'conditional_joint_scalar_bubble_candidate',
                            'normalized_input':str(meta),'C5_result':str(target),
                            'mu_s_gev':inp.mu_s_gev,
                            'weinberg_projection_residual':result['weinberg_projection_residual'],
                            'blockers':['complete full-theory finite one-loop scalar matching not independently certified']})
        except (OSError,ValueError,KeyError,IndexError,TypeError,RuntimeError) as err:
            row['status']='blocked'
            row['error']=f'{type(err).__name__}: {err}'
        reports.append(row)
        print(f'{model}: {row["status"]}',flush=True)
        if 'error' in row:
            print('  '+row['error'],flush=True)
    n_ready=sum(x['status']=='conditional_joint_scalar_bubble_candidate' for x in reports)
    n_formal=sum('formal_muF_npz' in x for x in reports)
    n_blocked=sum(x['status']=='blocked' for x in reports)
    summary={'format':BRIDGE_FORMAT,'source_release':str(Path(args.release_dir).resolve()),
             'eft1_endpoint_dir':str(Path(args.endpoint_dir).resolve()) if args.endpoint_dir else None,
             'models_requested':len(selected),'source_muF_exported':n_formal,
             'conditional_matching_candidates':n_ready,'blocked_count':n_blocked,
             'physical_C5_certified':False,
             'status':('FormalMuFSourceAuditOnly' if args.endpoint_dir is None and n_blocked==0
                       else 'CandidateMatchingCompleted' if args.endpoint_dir is not None and n_ready==len(selected)
                       else 'ExplicitModelBlockers'),
             'sequence':'UV -> all F at muF -> EFT1(S1,S2) -> both scalars at one muS -> SMEFT',
             'single_scalar_threshold_only':True,
             'models':reports}
    (dest/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('Summary:',dest/'summary.json',flush=True)
    return 2 if n_blocked or (args.endpoint_dir is not None and n_ready!=len(selected)) else 0


if __name__=='__main__':
    raise SystemExit(main())
