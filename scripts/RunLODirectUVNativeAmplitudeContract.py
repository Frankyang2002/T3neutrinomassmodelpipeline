#!/usr/bin/env python3
"""Read-only Workstream 4 native source-amplitude status witness.

Uses the CURRENT local repository files; does not call Wolfram or perform EFT
matching. Source mismatches cause a nonzero exit without a success report.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

PROJECT = Path(__file__).resolve().parents[1]
if str(PROJECT) not in sys.path:
    sys.path.insert(0, str(PROJECT))
from Numerical.lo_certification.direct_uv.native_amplitude_contract import (
    report, ContractError, audit_frozen_inputs)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--project-root', type=Path, default=PROJECT)
    ap.add_argument('--output-dir', type=Path, default=None)
    ap.add_argument('--check-frozen-inputs', action='store_true',
                    help='Audit the existing all-16 component input shapes and orientations')
    args = ap.parse_args()
    root = args.project_root.resolve()
    allowed = (root / 'output/general_observable/lo_certification_direct_uv').resolve()
    out = (args.output_dir or allowed / 'native_amplitude_contract_01').resolve()
    try:
        if not out.is_relative_to(allowed):
            raise ContractError('Report must remain under owned Workstream 4 output path')
        target = out / 'native_amplitude_source_contract.json'
        if target.exists():
            raise ContractError('Refusing to overwrite existing report; select new output dir')
        result = report(root)
        if args.check_frozen_inputs:
            result['actual_restricted_input_audit'] = audit_frozen_inputs(root)
        out.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n',
                          encoding='utf-8')
    except (ContractError, OSError, ValueError, TypeError, KeyError) as exc:
        print('Native amplitude SOURCE CONTRACT BLOCKED: ' + str(exc), file=sys.stderr)
        return 2
    print('All-16 source-operator dictionaries: 16/16; representative mass witnesses: 2/2')
    print('Physical Matchete/RGBeta-to-Weyl normalisation: 0/16 certified')
    if args.check_frozen_inputs:
        print('Actual restricted frozen input tensor structures: 16/16 verified')
    print('No one-loop matrices or EFT thresholds recalculated')
    print('Report: ' + str(target))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
