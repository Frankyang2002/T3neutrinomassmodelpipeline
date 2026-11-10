#!/usr/bin/env python3
"""Batch 10 A(-4)/D(-2): source-based formal F hard + full-general HHSS.

python -X utf8 -m scripts.RunT3HardGeneralBatch10 --all
python -X utf8 -m scripts.RunT3HardGeneralBatch10 --all --illustrative-smeft

No renormalized physical C5 is asserted. Reads saved Batch9 success, keeps old
results immutable, and writes source-resolved numerical one-loop diagnostics.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import traceback

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from Numerical.parallel_t3.complete.two_model_hard_general_batch10 import (
    B9, B1_SUBDIR, B2_SUBDIR, EXTRA_D, REGISTRY,
    FORMAT, MODELS, run_one,
)

DEST = Path('output/parallel_t3/thesis_hard_general_hhss_batch10_v1')


def write_json(path: Path, obj: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix('.json.partial')
    tmp.write_text(json.dumps(obj, indent=2, ensure_ascii=False,
                              sort_keys=True, allow_nan=False) + '\n', encoding='utf-8')
    tmp.replace(path)


def required_files(root: Path, model: str) -> list[Path]:
    src = root / B1_SUBDIR / model
    raw = root / B2_SUBDIR / model
    tag = 'A_m4' if model == 'A_-4' else 'D_m2'
    files = [root/B9/model/'result.json', root/B9/model/(tag+'.json'),
             root/B9/model/(tag+'.npz'), src/'manifest.json',
             src/'trajectories.npz', src/'fermion_matching_input.json',
             src/'scalar_matching_input.json',
             raw/'one_loop_operator_ledger.json',
             raw/'native_matching/data/eft1_after_F_wilson_seed.json']
    if model == 'D_-2':
        files += [root/EXTRA_D,root/REGISTRY]
    return files


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    selection = p.add_mutually_exclusive_group(required=True)
    selection.add_argument('--all', action='store_true')
    selection.add_argument('--model', choices=MODELS)
    p.add_argument('--project-root', type=Path, default=ROOT)
    p.add_argument('--output', type=Path, default=DEST,
                   help='New path under output/parallel_t3; never overwrite old runs')
    p.add_argument('--preflight', action='store_true')
    p.add_argument('--skip-a4-seven-term', action='store_true',
                   help='Skip optional A(-4) source-candidate full 7-channel contraction')
    p.add_argument('--illustrative-smeft', action='store_true',
                   help='Run numeric SM/Weinberg/PMNS under explicitly UNVERIFIED '
                        'common signed factor N_GammaCC=+1; NOT physical predictions')
    p.add_argument('--debug', action='store_true')
    args = p.parse_args(argv)
    root = args.project_root.resolve()
    output = (root / args.output).resolve()
    safe = (root / 'output/parallel_t3').resolve()
    if safe not in output.parents or output == safe:
        p.error('Output path must be a new child of output/parallel_t3')
    selected = MODELS if args.all else (args.model,)
    if args.preflight:
        missing = False
        for model in selected:
            for path in required_files(root, model):
                present = path.is_file()
                print(('FOUND ' if present else 'MISSING ') + str(path))
                missing |= not present
        return 2 if missing else 0
    if output.exists():
        p.error('Refusing to overwrite pre-existing outputs: '+str(output)+
                '; select a fresh --output')
    rows = []
    for model in selected:
        try:
            result = run_one(root, model, output/model,
                             include_a4_seven=not args.skip_a4_seven_term,
                             illustrative_smeft=args.illustrative_smeft)
            success = result['status'].startswith('HardAndGeneralScalarSourceComponentsComputed')
        except (ValueError, FileNotFoundError, OSError, TypeError,
                KeyError, ImportError, RuntimeError, FloatingPointError) as exc:
            success = False
            result = {'format':FORMAT, 'model':model,
                      'status':'Blocked_SourceOrConventionGap',
                      'reason':type(exc).__name__+': '+str(exc),
                      'physical_full_MSbar_C5_computed':False}
            if args.debug:
                result['traceback_tail'] = traceback.format_exc(limit=14)
        write_json(output/model/'result.json',result)
        rows.append({'model':model, 'status':result['status'],
                     'formal_source_sum_frobenius':result.get('formal_total_Frobenius'),
                     'reason':result.get('reason')})
        print(model,result['status'], flush=True)
        if not success:
            print(result['reason'],file=sys.stderr,flush=True)
    n = sum(x['status'].startswith('HardAndGeneralScalarSourceComponentsComputed')
            for x in rows)
    write_json(output/'summary.json', {
        'format': FORMAT+'-summary', 'models':rows,
        'success_count':n, 'expected_models':len(rows),
        'complete_physical_matching':False, 'physical_fit':False,
        'source_inputs_overwritten':False,
    })
    print('Summary:',output/'summary.json')
    return 0 if n == len(rows) else 2


if __name__ == '__main__':
    raise SystemExit(main())
