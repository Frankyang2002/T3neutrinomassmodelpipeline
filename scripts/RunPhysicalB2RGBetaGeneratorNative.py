#!/usr/bin/env python3
r"""B2 G2: launch/validate a fresh native RGBeta *symbolic generator* probe.

From project root (new output only):
  .\.venv\Scripts\python.exe scripts\RunPhysicalB2RGBetaGeneratorNative.py

Installed RGBeta is required. This does not resolve numerical SU(2) axes or
certify a physical Matchete/RGBeta coupling dictionary.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'Numerical/lo_certification/physical_normalisation'))
from b2_gauge_coupling import AuditBlocked
from b2_rgbeta_generator_native import OUTPUT_DEFAULT, run_native


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, default=OUTPUT_DEFAULT)
    p.add_argument('--wolfram-exe', default='wolframscript')
    args = p.parse_args(argv)
    output = args.output if args.output.is_absolute() else ROOT/args.output
    try:
        report = run_native(ROOT, output, wolfram=args.wolfram_exe)
    except (AuditBlocked, OSError, ValueError, KeyError, TypeError) as exc:
        print('BLOCKED_B2_NATIVE_SYMBOLIC_GENERATOR:', exc, file=sys.stderr)
        return 3
    print(f"{report['status']}; native families={report['families_registered']}/5; "
          f"numeric SU2 axes=NO; physical=0/16; output={output}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
