#!/usr/bin/env python3
"""Fail-closed all-16 lambdaT3 mixing audit against live registered Matchete CG.

Does not overwrite existing output, input benchmarks or EFT data.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

PROJECT = Path(__file__).resolve().parents[1]
if str(PROJECT) not in sys.path:
    sys.path.insert(0,str(PROJECT))

from Numerical.lo_certification.direct_uv.native_mixing_dictionary import (  # noqa: E402
    audit, MixingAuditError,
)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--project-root',type=Path,default=PROJECT)
    ap.add_argument('--native-json',type=Path,required=True)
    ap.add_argument('--output-dir',type=Path,default=Path(
        'output/general_observable/lo_certification_direct_uv/native_mixing_dictionary_01'))
    ap.add_argument('--no-frozen-inputs',action='store_true',
                    help='Verify native source and five tensors only')
    args = ap.parse_args()
    root = args.project_root.resolve()
    native = (args.native_json if args.native_json.is_absolute() else root/args.native_json).resolve()
    out = (args.output_dir if args.output_dir.is_absolute() else root/args.output_dir).resolve()
    owned = (root/'output/general_observable/lo_certification_direct_uv').resolve()
    try:
        if not out.is_relative_to(owned):
            raise MixingAuditError('Output must be within Workstream 4 output tree')
        if out.exists() and any(out.iterdir()):
            raise MixingAuditError('Output directory not empty; use a new output directory')
        result = audit(root,native,check_frozen=not args.no_frozen_inputs)
        out.mkdir(parents=True,exist_ok=True)
        dest = out/'native_mixing_scalar_A_dictionary.json'
        dest.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
        print('Live registered Matchete mixing CG classes: 5/5')
        print('Native HH exchange and non-HC/HC source structure: verified')
        if result['actual_frozen_inputs']:
            print('Live native CG vs original component source: 16/16')
            print('Frozen scalar-A off-diagonal mixing: 16/16')
            print('Independent real Cartesian scalar Hessians: 16/16')
            print('Maximum conditional mixing residual:',
                  result['actual_frozen_inputs']['max_scalar_A_mixing_relative_residual'])
        print('Full physical Matchete/RGBeta-to-Weyl normalisation: 0/16 certified')
        print('Report:',dest)
        return 0
    except (MixingAuditError,KeyError,ValueError,OSError,TypeError,RuntimeError,
            json.JSONDecodeError) as exc:
        print('Native lambdaT3 mixing dictionary BLOCKED: '+str(exc),file=sys.stderr)
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
