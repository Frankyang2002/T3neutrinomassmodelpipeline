#!/usr/bin/env python3
"""One-command A(-4)/D(-2) source-backed joint-scalar conditional calculation.

    python -X utf8 -m scripts.RunT3FullGeneralJointScalarBatch9 --all
    python -X utf8 -m scripts.RunT3FullGeneralJointScalarBatch9 --model A_-4

Only the scalar-pair leading log and *joint* finite scalar bubble are added;
no absent F-hard or EFT1 channels are silently approximated. Reads the saved
full-general native RG and Matchete source, and writes only a new output tree.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import sys
import traceback

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from Numerical.parallel_t3.complete.two_model_joint_scalar_batch9 import (
    B1_SUBDIR, B2_SUBDIR, MODELS, FORMAT, IntegrationBlocked,
    execute_model, quartic_source_path,
)

DEST = Path('output/parallel_t3/thesis_full_general_joint_scalar_batch9_v1')


def write_json(path: Path, obj: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix('.json.partial')
    tmp.write_text(json.dumps(obj, indent=2, sort_keys=True,
                              ensure_ascii=False, allow_nan=False) + '\n',
                   encoding='utf-8')
    tmp.replace(path)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    select = parser.add_mutually_exclusive_group(required=True)
    select.add_argument('--all', action='store_true')
    select.add_argument('--model', choices=MODELS)
    parser.add_argument('--project-root', type=Path, default=ROOT)
    parser.add_argument('--output', type=Path, default=DEST,
                        help='New output directory relative to project root')
    parser.add_argument('--preflight', action='store_true',
                        help='List required sources, without running physics')
    parser.add_argument('--debug', action='store_true')
    args = parser.parse_args(argv)
    root = args.project_root.resolve()
    output = (root / args.output).resolve()
    mandatory = (root / 'output/parallel_t3').resolve()
    if not output.is_relative_to(mandatory):
        parser.error('Output must remain in the project output/parallel_t3 tree')
    if output.exists() and not args.preflight:
        parser.error(f'Refusing overwrite: {output}; provide a fresh --output')
    models = MODELS if args.all else (args.model,)
    if args.preflight:
        for model in models:
            required = [root/B1_SUBDIR/model/'manifest.json',
                        root/B1_SUBDIR/model/'trajectories.npz',
                        root/B1_SUBDIR/model/'fermion_matching_input.json',
                        root/B1_SUBDIR/model/'scalar_matching_input.json',
                        root/B2_SUBDIR/model/'native_matching/native_status.json',
                        root/B2_SUBDIR/model/'native_matching/data/eft1_after_F_wilson_seed.json',
                        quartic_source_path(root, model)]
            for path in required:
                print(('FOUND ' if path.is_file() else 'MISSING ') + str(path))
        return int(any(not path.is_file() for model in models for path in (
            root/B1_SUBDIR/model/'manifest.json',
            root/B1_SUBDIR/model/'trajectories.npz',
            root/B2_SUBDIR/model/'native_matching/native_status.json',
            quartic_source_path(root, model))))
    rows = []
    for model in models:
        try:
            record = execute_model(root, model, output/model)
            code = 0
        except (IntegrationBlocked, ValueError, FileNotFoundError,
                KeyError, ImportError, OSError, RuntimeError, TypeError) as exc:
            record = {'format': FORMAT, 'model': model,
                      'status': 'Blocked_ExistingSourceMappingOrNumericalHandoff',
                      'reason': f'{type(exc).__name__}: {exc}',
                      'physical_C5_computed': False}
            if args.debug:
                record['traceback_tail'] = traceback.format_exc(limit=10)
            code = 1
        write_json(output/model/'result.json', record)
        rows.append({'model': model, 'status': record['status'],
                     'reason': record.get('reason'),
                     'Frobenius_per_unit_GeV_inverse':
                         record.get('LLSS_response_Frobenius_GeV_inverse')})
        print(model, record['status'], flush=True)
        if code:
            print(record['reason'], file=sys.stderr, flush=True)
    write_json(output/'summary.json', {
        'format': FORMAT+'-summary', 'rows': rows,
        'conditional_success_count': sum(x['status'].startswith(
            'SourceBasedLLSSScalarLoopComputed') for x in rows),
        'physical_C5_or_18_real_fit_computed': False,
        'original_physics_files_overwritten': False,
    })
    print('Summary:', output/'summary.json')
    return 0 if all(x['status'].startswith('SourceBasedLLSSScalarLoopComputed')
                    for x in rows) else 2


if __name__ == '__main__':
    raise SystemExit(main())
