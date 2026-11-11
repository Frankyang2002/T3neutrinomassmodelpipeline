"""Rebuild candidate same-muS scalar backgrounds for all 16 T3 models.

Intended project path: scripts/RunBatch17RebuildAll16.py
Preserves prior A_-4 result and all old data; use a NEW --out-dir for new runs.
Physical C5 remains explicitly uncertified; special-model field maps are
exploratory. Block individual models rather than using a wrong tensor basis.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import traceback

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))

from Numerical.experimental.Batch17All16Background import rebuild_background
from Numerical.experimental.Batch17ConditionalWilson import MODELS,run_conditional
from Numerical.joint_scalar_matching.fermion_boundary_bridge import load_release_manifest,slug
from Numerical.joint_scalar_matching.source_quartic_v4 import source_seed_path


def main(argv=None)->int:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--models',help='Comma-separated models, default all16')
    ap.add_argument('--source-release',type=Path,default=ROOT/'output/general_observable/lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f')
    ap.add_argument('--batch14',type=Path,default=ROOT/'output/batch14/eft1')
    ap.add_argument('--batch15',type=Path,default=ROOT/'output/batch15/scalar_threshold')
    ap.add_argument('--batch16',type=Path,default=ROOT/'output/batch16/scalar_mass_audit')
    ap.add_argument('--quartic-cache-dir',type=Path,default=ROOT/'output/general_observable/lo_joint_scalar_matching/source_quartic_v4/native_quartic_channels')
    ap.add_argument('--quartic-seed-dir',type=Path,default=ROOT/'output/full16_symbolic')
    ap.add_argument('--registry',type=Path,default=ROOT/'RGE/general/generated/general_rgbeta_quartic_registry.json')
    ap.add_argument('--out-dir',type=Path,default=ROOT/'output/batch17/rebuilt_all16_v1')
    ap.add_argument('--background-only',action='store_true')
    ap.add_argument('--mu-low',type=float,default=100.)
    args=ap.parse_args(argv)
    models=tuple(x.strip() for x in args.models.split(',')) if args.models else MODELS
    if not models or len(set(models))!=len(models) or not set(models)<=set(MODELS):
        ap.error('Invalid/duplicate --models; names must be in the 16-model registry')
    # Never overwrite output files from previous runs; also refuse to overwrite
    # a prior summary, allowing explicit --out-dir to name a fresh run.
    if (args.out_dir/'summary.json').exists():
        raise FileExistsError(f'Existing run summary: {args.out_dir / "summary.json"}; choose a new --out-dir')
    manifest=load_release_manifest(args.source_release)
    summary={
       'format':'t3-batch17-all16-same-endpoint-conditional-v1',
       'physical_neutrino_mass_prediction':False,
       'interpretation':'All candidate masses are per unit unknown GammaCC-to-Weyl N; candidate scalar coupling map, special mass-field maps and historical F-boundaries are uncertified',
       'models':{},
    }
    args.out_dir.mkdir(parents=True,exist_ok=True)
    for model in models:
        try:
            endpoint=args.batch14/f'{model}_EFT1_running.json'
            handoff=args.batch15/f'{model}_scalar_threshold_handoff.json'
            mass=args.batch16/f'{model}_mass_hessian.json'
            cache=args.quartic_cache_dir/f'{slug(model)}.json'
            seed=source_seed_path(args.quartic_seed_dir,model)
            source=manifest[model]
            files=(endpoint,handoff,mass,cache,seed,source,args.registry)
            missing=[str(p) for p in files if not p.is_file()]
            if missing:raise FileNotFoundError('Missing required files: '+'; '.join(missing))
            background=args.out_dir/f'{slug(model)}.json'
            built=rebuild_background(model=model,endpoint_path=endpoint,
                handoff_path=handoff,mass_path=mass,source_path=source,
                quartic_cache_path=cache,quartic_seed_path=seed,
                registry_path=args.registry,destination=background)
            print(f'[{model}] REBUILT | channels={built["n_quartic_channels"]} | spectrum rel={built["batch16_spectrum_relative_difference"]:.2e} | special={built["special_mass_basis_candidate"]}',flush=True)
            record={'status':'BackgroundRebuiltConditional','background':built}
            if not args.background_only:
                candidate=run_conditional(model,handoff_path=handoff,mass_path=mass,
                    endpoint_path=endpoint,source_path=source,background_path=background,
                    mu_low_gev=args.mu_low)
                out=args.out_dir/f'{model}_conditional_c5.json'
                if out.exists():raise FileExistsError(f'Refusing to overwrite {out}')
                out.write_text(json.dumps(candidate,indent=2,allow_nan=False)+'\n',encoding='utf-8')
                record.update(status='ComputedConditional',
                    candidate_masses_per_unit_N_eV=candidate['candidate_neutrino_masses_per_unit_N_eV'],
                    smeft_nfev=candidate['smeft_nfev'],conditional_c5_result=str(out))
                print(f'[{model}] CONDITIONAL | joint scalar LLSS + SMEFT | physical C5=UNCERTIFIED',flush=True)
            summary['models'][model]=record
        except Exception as exc:
            log=args.out_dir/f'{model}_error.log'
            log.write_text(traceback.format_exc(),encoding='utf-8')
            summary['models'][model]={'status':'Blocked','reason':str(exc),'error_log':str(log)}
            print(f'[{model}] BLOCKED: {exc}',flush=True)
    summary['rebuilt']=sum(x['status']!='Blocked' for x in summary['models'].values())
    summary['conditional']=sum(x['status']=='ComputedConditional' for x in summary['models'].values())
    summary['blocked']=sum(x['status']=='Blocked' for x in summary['models'].values())
    path=args.out_dir/'summary.json'
    path.write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(f'Summary: {path} | rebuilt={summary["rebuilt"]} conditional={summary["conditional"]} blocked={summary["blocked"]} | physical_C5_certified=0',flush=True)
    return 0 if summary['blocked']==0 else 1


if __name__=='__main__':
    raise SystemExit(main())
