"""Batch17: strict conditional LLSS -> joint scalar matching -> SMEFT all16.

Project path: scripts/RunBatch17ConditionalC5All16.py

You must supply a source-resolved Matchete LLSS release and an existing
same-scale joint scalar HHSS background directory, using --source-release
and --background-dir. Without them, print BLOCKED instead of inventing C5.

Each background file is expected at <background-dir>/<MODEL_SLUG>.json
(e.g., A_p0.json -> Ap0.json, A_-2 -> Am2.json). Override with
--background for a single --model.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import traceback

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
from Numerical.experimental.Batch17ConditionalWilson import MODELS,run_conditional


def slug(model:str)->str:
    return model.replace('+','p').replace('-','m').replace('_','')


def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--model',choices=MODELS)
    ap.add_argument('--source-release',type=Path,
                    help='Existing 16-model Matchete source-boundary release directory')
    ap.add_argument('--background-dir',type=Path,
                    help='Existing explicit joint-scalar HHSS backgrounds directory')
    ap.add_argument('--background',type=Path,
                    help='Override single-model background JSON path')
    ap.add_argument('--batch14',type=Path,default=ROOT/'output/batch14/eft1')
    ap.add_argument('--batch15',type=Path,default=ROOT/'output/batch15/scalar_threshold')
    ap.add_argument('--batch16',type=Path,default=ROOT/'output/batch16/scalar_mass_audit')
    ap.add_argument('--out',type=Path,default=ROOT/'output/batch17/conditional_c5')
    ap.add_argument('--mu-low',type=float,default=100.)
    args=ap.parse_args(argv)
    if args.background is not None and args.model is None:
        ap.error('--background requires --model')
    args.out.mkdir(parents=True,exist_ok=True)
    models=(args.model,) if args.model else MODELS
    manifest={}
    if args.source_release is not None:
        from Numerical.joint_scalar_matching.fermion_boundary_bridge import load_release_manifest
        manifest=load_release_manifest(args.source_release)
    summary={'format':'t3-batch17-conditional-summary-v1',
             'physical_neutrino_mass_prediction':False,
             'interpretation':'Conditional LLSS matching per unit unknown source-to-Weyl N; not certified physical C5',
             'models':{}}
    for model in models:
        try:
            handoff=args.batch15/f'{model}_scalar_threshold_handoff.json'
            mass=args.batch16/f'{model}_mass_hessian.json'
            endpoint=args.batch14/f'{model}_EFT1_running.json'
            for path in (handoff,mass,endpoint):
                if not path.is_file():raise ValueError(f'Missing completed Batch14-16 file: {path}')
            if not manifest:
                raise ValueError('Source-resolved LLSS release not supplied (--source-release); cannot construct C5')
            source=manifest[model]
            background=(args.background if args.background is not None else
                        args.background_dir/ f'{slug(model)}.json' if args.background_dir is not None else None)
            if background is None or not background.is_file():
                raise ValueError('Explicit muS HHSS background missing; specify --background-dir or --background')
            report=run_conditional(model,handoff_path=handoff,mass_path=mass,
                endpoint_path=endpoint,source_path=source,background_path=background,
                mu_low_gev=args.mu_low)
            (args.out/f'{model}_conditional_c5.json').write_text(
                json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
            summary['models'][model]={'status':'ComputedConditional',
                'candidate_masses_per_unit_N_eV':report['candidate_neutrino_masses_per_unit_N_eV']}
            print(f'[{model}] CONDITIONAL | joint scalar C5 + SMEFT computed per unit N; '
                  f'physical C5=UNCERTIFIED',flush=True)
        except Exception as exc:
            (args.out/f'{model}_error.log').write_text(traceback.format_exc(),encoding='utf-8')
            summary['models'][model]={'status':'Blocked','reason':str(exc)}
            print(f'[{model}] BLOCKED: {exc}',flush=True)
    summary['computed_conditional']=sum(x['status']=='ComputedConditional' for x in summary['models'].values())
    summary['blocked']=len(models)-summary['computed_conditional']
    (args.out/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(f'Summary: {args.out / "summary.json"} | conditional={summary["computed_conditional"]} '
          f'blocked={summary["blocked"]} | physical_C5_certified=0')
    return 0 if summary['blocked']==0 else 1

if __name__=='__main__':raise SystemExit(main())
