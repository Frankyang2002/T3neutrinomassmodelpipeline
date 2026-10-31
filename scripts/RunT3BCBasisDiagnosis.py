"""B/C RGBeta covariance diagnosis; uses saved sources by default."""
from __future__ import annotations
import argparse
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from Numerical.experimental.T3BCBasisDiagnosis import run


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--live', action='store_true', help='Re-evaluate native exported UV beta functions at an independent probe point')
    parser.add_argument('--output', type=Path, default=ROOT / 'output/studies/t3_bc_basis_diagnostic_v1/summary.json')
    args = parser.parse_args()
    result = run(ROOT, args.output, live=args.live)
    print(result['status'])
    for name, row in result['models'].items():
        print(name, 'CConj-failed:', row['CConj_orientation_failed_blocks_at_1e-6'],
              'same-U max:', row['same_U_covariance_max_relative'])
    print('Output:', args.output)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
