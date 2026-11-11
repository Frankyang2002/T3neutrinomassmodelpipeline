"""Additive Worker C real-scale all-16 audit and evolved-EFT1 joint threshold."""
from __future__ import annotations

import argparse
from datetime import datetime,timezone
from pathlib import Path
from secrets import token_hex

from Numerical.joint_scalar_matching.worker_c_stream_v2 import (
    MODELS, HandoffError, endpoint_contract_template, read_doc, remaining_work,
    run_endpoint, slug, survey_all, _json_dump,
)

SOURCE_DEFAULT = Path('output/general_observable/lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f')
BACKGROUND_DEFAULT = Path('output/general_observable/lo_joint_scalar_matching/exploratory_all16_v6_native_fix/numeric_backgrounds')
ARCHIVE_DEFAULT = Path('output/general_observable/lo_joint_scalar_matching/exploratory_all16_v6_native_fix/C5_per_unknown_NGamma_exploratory')


def fresh_output(root: Path,kind: str) -> Path:
    stamp=datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S_%f')
    return Path(root)/f'worker_c_{kind}_{stamp}_{token_hex(4)}'


def main(argv=None) -> int:
    ap=argparse.ArgumentParser(description=__doc__)
    modes=ap.add_mutually_exclusive_group(required=True)
    modes.add_argument('--survey-all16', action='store_true',help='16 source+scalar-background conditional re-evaluation at the established scales')
    modes.add_argument('--endpoint',type=Path,help='Actual, evolved EFT1 Wilson+scalar JSON/NPZ handoff')
    modes.add_argument('--contract',choices=MODELS,help='Write a deliberately incomplete handoff specification for one model')
    modes.add_argument('--remaining',action='store_true',help='Print outstanding implementation inventory')
    ap.add_argument('--models',nargs='+',choices=MODELS,help='Subset for --survey-all16, otherwise all 16')
    ap.add_argument('--source-dir',type=Path,default=SOURCE_DEFAULT)
    ap.add_argument('--background-dir',type=Path,default=BACKGROUND_DEFAULT)
    ap.add_argument('--archive-dir',type=Path,default=ARCHIVE_DEFAULT,
                    help='Optional historical same-input regression folder; use --no-archive to disable')
    ap.add_argument('--no-archive',action='store_true')
    ap.add_argument('--output-root',type=Path,default=Path('output/audits'))
    ap.add_argument('--output-dir',type=Path,help='A new (nonexistent) output folder')
    ap.add_argument('--stop-on-error',action='store_true')
    args=ap.parse_args(argv)
    try:
        if args.remaining:
            for x in remaining_work():
                print(f"{x['owner']}: {x['status']}: {x['item']}")
            return 0
        kind=('all16' if args.survey_all16 else 'endpoint' if args.endpoint else 'contract')
        out=args.output_dir if args.output_dir else fresh_output(args.output_root,kind)
        if args.survey_all16:
            result=survey_all(Path.cwd(),args.source_dir,args.background_dir,
                              models=tuple(args.models or MODELS),output_dir=out,
                              keep_going=not args.stop_on_error,
                              archive_dir=None if args.no_archive else args.archive_dir)
            print(f"Worker C: {result['conditionally_evaluated']}/{result['requested_model_count']} conditional models; blocked={result['blocked']}")
            print(f'Summary: {out / "summary.json"}')
            print('PHYSICAL SMEFT C5 CERTIFIED: False')
            return 0 if result['blocked']==0 else 2
        if out.exists():
            raise HandoffError(f'Output directory already exists: {out}')
        if args.contract:
            out.mkdir(parents=True)
            target=out/(slug(args.contract)+'.json')
            _json_dump(target,endpoint_contract_template(args.contract))
            print(f'INCOMPLETE handoff template: {target}')
            return 0
        evaluated=run_endpoint(args.endpoint)
        out.mkdir(parents=True)
        _json_dump(out/'result.json',evaluated)
        print(f'Joint scalar one-bubble endpoint: {out / "result.json"}')
        print('COMPLETE ONE-LOOP SCALAR MATCHING: False')
        return 0
    except (HandoffError,ValueError,KeyError,ImportError,OSError) as exc:
        print(f'BLOCKED: {exc}')
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
