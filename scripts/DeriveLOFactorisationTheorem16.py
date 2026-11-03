"""Chat A Batch 7: source-conditional analytical synthesis, read-only on upstream data."""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0,str(ROOT))
from Numerical.lo_patterns_16.factorization_theorem import analyse, TARGET


def main()->None:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,help='Fresh output directory inside Chat A result tree')
    args=p.parse_args()
    data=analyse(ROOT,output_dir=args.output)
    print(f'Batch 7: {data["model_count"]}/16 source-conditional CG+theorem analyses; physical certification 0/16')
    print(f'Max baseline coefficient residual: {data["max_baseline_residual"]:.6e}')
    print(f'Max |candidate CG magnitude ratio - measured |kappa||: '
          f'{data["max_cg_kappa_magnitude_disagreement"]:.6e}')
    for row in data['correction_tests']:
        print(f'{row["model"]:5} split-first-order={row["charge_split_linearized_relative_error"]:.3e} '
              f'large-mix-linear={row["large_mix_linear_relative_error"]:.3e} '
              f'large-mix-cubic={row["large_mix_cubic_relative_error"]:.3e}')
    print('New output: '+str(args.output or TARGET))


if __name__=='__main__':main()
