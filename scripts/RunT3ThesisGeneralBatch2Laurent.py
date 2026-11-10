#!/usr/bin/env python3
"""Batch 2 next step: evaluate the saved *complete* F one-loop Matchete EFT.

Does NOT rerun Batch 1 RG or F-threshold Match. Evaluates Matchete loop
functions on all operators, extracts explicit epsilon pole/regular pieces when
possible, and records algebraic flavour-mass orientation diagnostics.

NOT a renormalised MSbar Wilson matching and NOT EFT1 Wilson evolution.
Uncertified results are never marked as complete numerical Wilson coefficients.

PowerShell (project root):
 python -X utf8 -m scripts.RunT3ThesisGeneralBatch2Laurent --all --preflight
 python -X utf8 -m scripts.RunT3ThesisGeneralBatch2Laurent --all \
     --run-id thesis_full_general_batch2_laurent_v1 --accept-formal-audit-only
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path
import os
import re
import shutil
import subprocess
import sys
import traceback

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from Numerical.parallel_t3.wilson.general_batch2_laurent import (
    MODELS, PIPELINE_FORMAT, SourceAuditError, check_audit_report,
    check_saved_raw, matrix_spurion_diagnostics, read_json, require, save_json, sha,
)

B1_ROOT = ROOT / 'output/parallel_t3/full_general_batch1/thesis_full_general_batch1_v1'
B2_RAW = ROOT / 'output/parallel_t3/full_general_batch2/thesis_full_general_batch2_raw_v1'
OUTPUT_ROOT = ROOT / 'output/parallel_t3/full_general_batch2'
AUDITOR = ROOT / 'Lagrangian/RunT3Batch2LaurentAudit.wl'
ENV_GENERAL = (
    'T3_EXPERIMENTAL_GENERAL_FERMIONS',
    'T3_EXPERIMENTAL_GENERAL_QUARTICS',
    'T3_EXPERIMENTAL_GENERAL_QUADRATICS',
    'T3_EXPERIMENTAL_FLAVOR_MASSES',
)


def safe_id(name: str) -> str:
    if re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,63}', name) is None:
        raise argparse.ArgumentTypeError('run-id must be 1..64 alphanumeric/underscore/hyphen')
    return name


def select_models(args: argparse.Namespace) -> tuple[str, ...]:
    if args.all:
        return MODELS
    if args.model:
        return (args.model,)
    raise ValueError('Select --all or --model')


def prepare_model(model: str, args: argparse.Namespace) -> tuple[dict, dict]:
    checked = check_saved_raw(model, args.batch1_root / model,
                             args.raw_root / model)
    numeric = matrix_spurion_diagnostics(checked['numeric'], model)
    record = {
        'format': PIPELINE_FORMAT,
        'status': 'SourceAndBatch1RawMatchingPreflightPassed',
        'model': model,
        'created_utc': datetime.now(timezone.utc).isoformat(),
        'raw_input_run_id': args.raw_root.name,
        'matching_scale_GeV': 1e5,
        'scalar_threshold_scale_GeV': 1e3,
        'source_path_relative_to_project': str(checked['one_loop_source'].resolve().relative_to(ROOT.resolve())).replace('\\', '/'),
        'batch1_manifest_sha256': checked['batch1_manifest_sha256'],
        'raw_manifest_sha256': checked['raw_manifest_sha256'],
        'raw_ledger_sha256': checked['raw_ledger_sha256'],
        'native_matching_status_sha256': checked['native_matching_status_sha256'],
        'raw_one_loop_source_sha256': checked['one_loop_source_sha256'],
        'raw_one_loop_operator_count': checked['raw_one_loop_term_count'],
        'operator_sector_counts_unfiltered': checked['candidate_sector_counts'],
        'audit_source_sha256': sha(AUDITOR),
        'native_fermion_flavour_diagnostics_file': 'native_flavour_diagnostics.json',
        'native_to_Matchete_physical_dictionary_verified': False,
        'numerical_Wilson_initial_conditions_produced': False,
        'EFT1_Wilson_RGE_executed': False,
        'full_one_loop_matching_completed': False,
        'scope': 'full raw F one-loop Matchete, no operator-selection or analytic-triangle projection',
        'physical_blockers': [
            'Matchete-to-RGBeta scalar and fermion coupling/CG normalization dictionary',
            'Explicit UV-hard minus EFT1-soft subtraction with counterterm/pole cancellation',
            'Numerical flavour/operator projection of all required finite EFT1 coefficients',
            'Complete EFT1 Wilson mixing with all active scalars and running parameters',
        ],
    }
    return record, numeric


def launch_wolfram(source: Path, dest: Path, args: argparse.Namespace, log: Path) -> None:
    env = dict(os.environ)
    env.update({key: '1' for key in ENV_GENERAL})
    cmd = [args.wolframscript, '-file', str(AUDITOR), str(source), str(dest),
           str(args.symbolic_seconds)]
    try:
        result = subprocess.run(cmd, cwd=ROOT, env=env, check=False,
                                capture_output=True, text=True,
                                timeout=args.timeout_sec)
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError('Matchete loop evaluation exceeded wall-clock timeout') from exc
    log.write_text('command: ' + ' '.join(cmd) + '\nexit: ' + str(result.returncode)
                   + '\nstdout:\n' + result.stdout + '\nstderr:\n' + result.stderr,
                   encoding='utf-8')
    if result.returncode:
        status_file = dest / 'laurent_audit.json'
        detail = read_json(status_file) if status_file.is_file() else {}
        raise RuntimeError('Wolfram loop audit failed: ' + str(detail.get('reason', 'unknown'))
                           + '; see ' + str(log))


def run(args: argparse.Namespace) -> list[dict]:
    require(AUDITOR.is_file(), 'missing Mathematica audit source ' + str(AUDITOR))
    models = select_models(args)
    require(not (args.output_root / args.run_id).exists() or args.preflight,
            'output run-id already exists; choose a new run-id (old results preserved)')
    planned: list[tuple[str, dict, dict]] = []
    for model in models:
        record, numeric = prepare_model(model, args)
        planned.append((model, record, numeric))
        if args.preflight:
            print(f'{model}: {record["status"]} '
                  f'({record["raw_one_loop_operator_count"]} raw operators)')
    if args.preflight:
        return [r for _, r, _ in planned]

    root = args.output_root / args.run_id
    root.mkdir(parents=True)
    rows = []
    for model, record, numeric in planned:
        dest = root / model
        dest.mkdir()
        save_json(dest / 'native_flavour_diagnostics.json', numeric)
        save_json(dest / 'manifest.json', record)
        source = args.raw_root / model / 'native_matching/data/eft1_after_F_one_loop_inputform.txt'
        try:
            launch_wolfram(source, dest, args, dest / 'wolfram_laurent.log')
            audit = read_json(dest / 'laurent_audit.json')
            result = check_audit_report(audit, record['raw_one_loop_source_sha256'], dest)
            record.update({'status': 'FullOneLoopFormalLaurentAudited_TreeMassRemainderPreserved_NotFiniteWilsons',
                           'symbolic_result': result,
                           'matchete_formal_loop_functions_evaluated':
                               audit.get('matchete_loop_function_evaluator_called') is True,
                           'physical_matching_completed': False})
        except Exception as exc:
            audit_file = dest / 'laurent_audit.json'
            try:
                audit = read_json(audit_file) if audit_file.is_file() else {}
            except (ValueError, OSError):
                audit = {}
            native_status = audit.get('status')
            record.update({
                'status': (native_status if isinstance(native_status, str)
                           and native_status.startswith('Blocked_')
                           else 'Blocked_' + type(exc).__name__),
                'reason': str(exc),
                'native_audit_reason': audit.get('reason'),
                'native_audit_status': native_status,
                'traceback_tail': traceback.format_exc(limit=4),
            })
        save_json(dest / 'manifest.json', record)
        print(f'{model}: {record["status"]}', flush=True)
        rows.append(record)
        save_json(root / 'summary.json', {
            'format': PIPELINE_FORMAT + '-summary',
            'run_id': args.run_id,
            'model_count': len(rows),
            'formal_audit_success_count': sum(x['status'].startswith('FullOneLoopFormal')
                                              for x in rows),
            'finite_MSbar_Wilson_matching_complete': False,
            'full_EFT1_Wilson_RGE_executed': False,
            'models': [{k: m.get(k) for k in ('model', 'status', 'reason')} for m in rows],
        })
    return rows


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    grp = p.add_mutually_exclusive_group(required=True)
    grp.add_argument('--all', action='store_true')
    grp.add_argument('--model', choices=MODELS)
    p.add_argument('--batch1-root', type=Path, default=B1_ROOT)
    p.add_argument('--raw-root', type=Path, default=B2_RAW)
    p.add_argument('--output-root', type=Path, default=OUTPUT_ROOT)
    p.add_argument('--run-id', type=safe_id, default='thesis_full_general_batch2_laurent_v1')
    p.add_argument('--wolframscript', default='wolframscript')
    p.add_argument('--symbolic-seconds', type=int, default=1200)
    p.add_argument('--timeout-sec', type=int, default=1800)
    p.add_argument('--preflight', action='store_true')
    p.add_argument('--accept-formal-audit-only', action='store_true',
                   help='Allow exit 0 for raw Laurent audit (NOT complete finite Wilson matching)')
    args = p.parse_args(argv)
    if not 30 <= args.symbolic_seconds <= 7200 or args.timeout_sec <= args.symbolic_seconds:
        p.error('Require 30 <= symbolic-seconds <= 7200 and timeout-sec > symbolic-seconds')
    if not args.preflight and not (shutil.which(args.wolframscript)
                                   or Path(args.wolframscript).is_file()):
        p.error('wolframscript is not available')
    try:
        rows = run(args)
    except (SourceAuditError, OSError, ValueError) as exc:
        p.error(str(exc))
    if args.preflight:
        print('Existing numerical Batch 1 and full raw one-loop Batch 2: PASS')
        return 0
    if any(not r['status'].startswith('FullOneLoopFormal') for r in rows):
        print('Some loop-function/Laurent audits blocked; see per-model reports')
        return 2
    print('Full raw one-loop Matchete loop-function audit: OK')
    print('Finite MSbar Wilson coefficients and complete EFT1 Wilson RGE: NOT CALCULATED')
    return 0 if args.accept_formal_audit_only else 3


if __name__ == '__main__':
    raise SystemExit(main())
