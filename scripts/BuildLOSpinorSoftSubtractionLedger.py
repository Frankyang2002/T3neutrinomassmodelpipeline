#!/usr/bin/env python3
"""Produce all-16 source-resolved EFT1-soft / CT *requirements* and scalar B0 masters.

Does not execute Matchete, perform soft subtraction, or certify C5.
Only writes new immutable files below output/general_observable/lo_certification_spinor/.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import re
import sys


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--native-run', type=Path, default=Path(
        'output/general_observable/lo_certification_spinor/native_hard_e515d82d118284ba_native01'))
    parser.add_argument('--inventory-run', type=Path, default=Path(
        'output/general_observable/lo_certification_spinor/raw_llhh_e515d82d118284ba_inventory01'))
    parser.add_argument('--bridge-run', type=Path, default=Path(
        'output/general_observable/lo_certification_spinor/pole_log_e515d82d118284ba_bridge01'))
    parser.add_argument('--run-id', default='soft01')
    parser.add_argument('--check-only', action='store_true')
    args = parser.parse_args(argv)
    root = args.repo_root.resolve()
    sys.path.insert(0, str(root))
    from Numerical.lo_certification.spinor.boundary import (
        MODELS, BoundaryError, checked_write, require, tag)
    from Numerical.lo_certification.spinor.soft_counterterm import (
        validate_bridge_snapshot, build_readiness)
    try:
        require(re.fullmatch(r'[A-Za-z0-9_-]{1,40}', args.run_id) is not None,
                'Unsafe run-id')
        owned = root / 'output/general_observable/lo_certification_spinor'
        require(owned.resolve().is_dir(), 'Owned output directory missing')

        def owned_input(p: Path) -> Path:
            found = (p if p.is_absolute() else root / p).resolve()
            require(found.is_relative_to(owned.resolve()),
                    f'Source escaped the spinor-owned output directory: {found}')
            return found

        native, inventory, bridge = [owned_input(p) for p in (
            args.native_run, args.inventory_run, args.bridge_run)]
        summary, models, native_summary = validate_bridge_snapshot(
            root, native, inventory, bridge)
        native_sha = sha256((native/'summary.json').read_bytes()).hexdigest()
        bridge_sha = sha256((bridge/'summary.json').read_bytes()).hexdigest()
        report, detailed = build_readiness(root, summary, models,
                                           native_summary, native_sha, bridge_sha)
        output = owned / ('soft_ct_' + report['source_release_sha256'][:16] + '_' + args.run_id)
        require(not output.exists(), 'Immutable output directory already exists: ' + str(output))
        if args.check_only:
            print(json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True))
            print('CHECK ONLY: no output written')
            return 0

        output.mkdir(parents=False, exist_ok=False)
        for model in MODELS:
            d = output / tag(model)
            checked_write(d / 'pole_and_soft_ct_requirements.json', detailed['ledgers'][model])
            checked_write(d / 'scalar_B0_master.json', detailed['scalar_masters'][model])
        checked_write(output / 'summary.json', report)
        print('CREATED:', output.relative_to(root))
        print('VERIFIED: 16/16 source-linked models, 176/176 native pole kernels')
        for model in ('B_-1', 'C_-1'):
            master = detailed['scalar_masters'][model]
            print(model, 'B0 finite at mu_F for XX,XY,YY:',
                  ', '.join(format(x['B0_MSbar_finite_at_mu_F'], '.12g')
                            for x in master['scalar_pairs']))
        print('BLOCKED: actual EFT1 soft loop and CT; Match hard/soft semantics; GammaCC phase; physical C5')
        return 0
    except (BoundaryError, ValueError, KeyError, OSError, TypeError, AssertionError) as exc:
        print('BLOCKED:', type(exc).__name__ + ':', str(exc), file=sys.stderr)
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
