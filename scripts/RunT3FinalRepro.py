#!/usr/bin/env python3
"""One-command complete conditional T3 reproducibility run.

Default: run the existing full Batch-7 tests and resume Batch 7 from saved
Batches 1/2/4/5, then validate and consolidate every saved branch. Existing
physics outputs are never altered by the final reporting module.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from Numerical.parallel_t3.wilson.final_repro_batch8 import (
    DEST, B7, ReproError, _load_json, verify_report, write_report,
)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, default=ROOT,
                   help='Project root; defaults to this script parent directory')
    p.add_argument('--out', type=Path, default=DEST, help='Separate final report directory')
    p.add_argument('--collect-only', action='store_true',
                   help='Use the finished 16/16 Batch-7 output; do not invoke runner')
    p.add_argument('--verify-only', action='store_true',
                   help='Check saved inputs and generated report without changing files')
    p.add_argument('--skip-tests', action='store_true', help='Skip initial local pytest')
    p.add_argument('--mu-low', type=float, default=None,
                   help='Override low scale when the Batch-7 runner is invoked')
    p.add_argument('--vev', type=float, default=None,
                   help='Override Higgs vev when Batch-7 runner is invoked')
    args = p.parse_args(argv)
    if args.verify_only and args.collect_only:
        p.error('Select --verify-only or --collect-only, not both')
    root = args.root.resolve()
    if not (root / 'scripts/RunT3PoleClosureBatch7.py').is_file():
        p.error('Missing installed Batch-7 runner in selected project root')
    try:
        if args.verify_only:
            manifest = verify_report(root, args.out)
            print('VERIFIED: 16/16 models; 144 scenarios; all input and report hashes match')
            print('Manifest:', (root/args.out).resolve()/'manifest.json')
            return 0
        if not args.skip_tests:
            testpaths = ('tests/test_t3_pole_closure_batch7.py',
                         'tests/test_t3_batch7_resume_cli.py',
                         'tests/test_t3_final_repro_batch8.py')
            print('Running Batch-7 and final report regression tests...', flush=True)
            subprocess.run([sys.executable, '-m', 'pytest', '-q', *testpaths],
                           cwd=root, check=True)
        if not args.collect_only:
            summary = root/B7/'summary.json'
            prior = _load_json(summary) if summary.is_file() else {}
            mu_low = args.mu_low if args.mu_low is not None else prior.get('mu_low_gev', 173.0)
            vev = args.vev if args.vev is not None else prior.get('vev_gev', 246.22)
            print('Resuming existing all-16 source-resolved pipeline...', flush=True)
            subprocess.run([sys.executable, 'scripts/RunT3PoleClosureBatch7.py',
                            '--resume', '--mu-low', str(mu_low), '--vev', str(vev)],
                           cwd=root, check=True)
        print('Checking immutable saved inputs and assembling final report...', flush=True)
        result = write_report(root, args.out)
        verify_report(root, args.out)
        print(f"COMPLETE: {result['completed_models']}/16 models; "
              f"{result['total_scenario_rows']} scenario rows; 0 blocked")
        print('Manifest:', (root/args.out).resolve()/'manifest.json')
        print('Summary:', (root/args.out).resolve()/'all16_overview.csv')
        print('All scenarios:', (root/args.out).resolve()/'all16_scenarios.csv')
        return 0
    except (ReproError, ValueError, FileNotFoundError, KeyError, OSError,
            subprocess.CalledProcessError) as exc:
        print(f'FINAL REPRODUCIBILITY FAILED: {type(exc).__name__}: {exc}',
              file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
