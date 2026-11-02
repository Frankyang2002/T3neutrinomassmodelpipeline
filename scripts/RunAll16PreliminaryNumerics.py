"""Read-only source-backed all-16 one-channel comparison; output only to lo_eft1_mixing."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))

from RGE.running.intermediate.All16PreliminaryNumerics import (
    FORMAT, MissingCertifiedBasis, evaluate_model, output_summary, summary_row, tag,
)
from Numerical.joint_scalar_matching.adapter import ALL_MODELS, normalize_model

SAFE=ROOT/'output/general_observable/lo_eft1_mixing'


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--model',help='One model (default all 16)')
    p.add_argument('--output-root',type=Path,default=SAFE/'all16_lambdaT3_preliminary_v1')
    p.add_argument('--check-only',action='store_true',help='Validate inputs and compute without writing')
    p.add_argument('--require-all',action='store_true',help='Exit code 2 if any model not Ward-valid')
    args=p.parse_args(argv)
    target=args.output_root.resolve(); base=SAFE.resolve()
    if target==base or base not in target.parents:
        p.error(f'Output must be a NEW SUBDIRECTORY of {base}')
    if target.exists() and (not target.is_dir() or any(target.iterdir())):
        p.error('Refusing overwrite/reuse of a nonempty output directory; choose a fresh --output-root')
    selected=[normalize_model(args.model)] if args.model else list(ALL_MODELS)
    docs=[];rows=[]
    for model in selected:
        try:
            result=evaluate_model(ROOT,model)
            docs.append(result)
            row=summary_row(result)
        except MissingCertifiedBasis as exc:
            row={'model':model,'status':'BlockedMissingCertifiedScalarBasis',
                 'reason':str(exc),'physical_mass_certified':False,
                 'masses_ev_per_abs_NGamma':None}
        except (OSError,ValueError,KeyError,IndexError,TypeError,RuntimeError) as exc:
            row={'model':model,'status':'InputError','reason':str(exc),
                 'physical_mass_certified':False,'masses_ev_per_abs_NGamma':None}
        rows.append(row)
        print(f'{model}: {row["status"]}' + (f' ({row["reason"]})' if 'reason' in row else ''),flush=True)
    summary=output_summary(rows)
    if args.check_only:
        print(json.dumps(summary,indent=2,allow_nan=False))
    else:
        target.mkdir(parents=True,exist_ok=True)
        for result in docs:
            path=target/tag(result['model'])/'lambdaT3_one_channel.json'
            path.parent.mkdir(parents=True,exist_ok=True)
            with path.open('x',encoding='utf-8') as fp:
                json.dump(result,fp,indent=2,allow_nan=False)
                fp.write('\n')
        with (target/'summary.json').open('x',encoding='utf-8') as fp:
            json.dump(summary,fp,indent=2,allow_nan=False)
            fp.write('\n')
        print('Summary:',target/'summary.json')
    return 2 if (summary['models_input_errors'] or (args.require_all and
        summary['models_with_Ward_valid_one_channel_results']!=len(rows))) else 0


if __name__=='__main__':
    raise SystemExit(main())
