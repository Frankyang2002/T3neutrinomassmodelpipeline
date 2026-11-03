"""Render the saved-source-verified Chat A analytical synthesis; no UV recomputation."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from Numerical.lo_patterns_16.final_synthesis import generate, DEST, SynthesisError


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output-dir', type=Path, default=None,
                   help='Must be a fresh child of the Chat A analysis tree')
    args = p.parse_args()
    try:
        d = generate(ROOT, output_dir=args.output_dir)
    except (SynthesisError, OSError, ValueError, KeyError, json.JSONDecodeError) as e:
        print('FINAL SYNTHESIS BLOCKED: ' + str(e), file=sys.stderr)
        return 1
    print('Final synthesis: 16/16 conditional results, 8/8 batch evidence sets; physical certification 0/16')
    print('Baseline heavy-coefficient worst residual: '
          f"{d['scenario_maxima'][0]['max_residual']:.6e}")
    print('Charge-splitting worst residual: '
          f"{d['scenario_maxima'][-1]['max_residual']:.6e}")
    print('Derivative finite-difference maximum error: '
          f"{d['max_exact_derivative_finite_difference_error']:.6e}")
    print('Output:', args.output_dir or ROOT / DEST)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
