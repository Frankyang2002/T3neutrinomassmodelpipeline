"""Run Chat A Batch 8 on saved conditional source inputs, no UV or EFT rerun."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Works from the repository checkout independently of the calling cwd.
PROJECT_ROOT=Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0,str(PROJECT_ROOT))
from Numerical.lo_patterns_16.splitting_convergence import analyse


def main() -> None:
    parser=argparse.ArgumentParser(description='All-16 exact scalar splitting O(t²) convergence')
    parser.add_argument('--output-dir',type=Path,help='New, empty directory inside the Chat A analysis root')
    args=parser.parse_args()
    p=analyse(PROJECT_ROOT,output_dir=args.output_dir)
    print('Batch 8: 16/16 conditional splitting derivatives; '+
          f'{len(p["sampled_t"])} t points/model; physical certification 0/16')
    print('Max exact derivative central-difference error:',
          f'{p["max_finite_difference_relative_error"]:.4e}')
    for row in p['model_summary']:
        model=row['model']
        t1=next(r for r in p['points'] if r['model']==model and r['t']==1.)
        order=row['exact_remainder_order']
        print(f'{model:<5} remainder_order={order if order is not None else "noise_floor"} '
              f'transverse_order={row["transverse_signal_order"] if row["transverse_signal_order"] is not None else "noise_floor"} '
              f'endpoint_factorization_residual={t1["nonproportional_residual"]:.4e}')
    from Numerical.lo_patterns_16.splitting_convergence import TARGET
    print('Output:',args.output_dir or TARGET)


if __name__=='__main__':
    main()
