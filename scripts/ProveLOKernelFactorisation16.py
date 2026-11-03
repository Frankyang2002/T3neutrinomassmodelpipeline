"""Run Chat A Batch 6: analytic CG-heavy-kernel factorisation stress tests (conditional)."""
from __future__ import annotations
import argparse
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from Numerical.lo_patterns_16.kernel_factorization import FactorizationError, analyse


def main() -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project-root',type=Path,default=ROOT)
    parser.add_argument('--batch05-dir',type=Path)
    parser.add_argument('--output-dir',type=Path,help='Fresh directory inside Chat A output root')
    parser.add_argument('--no-plots',action='store_true')
    opts=parser.parse_args()
    root=opts.project_root.resolve()
    try:
        result=analyse(root,
            batch05_dir=root/opts.batch05_dir if opts.batch05_dir else None,
            output_dir=root/opts.output_dir if opts.output_dir else None,
            make_plots=not opts.no_plots)
    except (FactorizationError,ValueError,KeyError,TypeError,RuntimeError,OSError) as exc:
        print(f'{type(exc).__name__}: {exc}',file=sys.stderr)
        return 2
    print('Batch 6: 16/16 conditional Yukawa-independent kernels; 5 scenarios; 80 group comparisons; physical certification 0/16')
    summary=result['model_summary']
    for record in summary:
        print(f"{record['model']:<5} baseline B_g residual={record['baseline_coefficient_residual']:.3e} "
              f"heavy-scale spread={record['baseline_heavy_scale_spread']:.3e} "
              f"max stress residual={record['largest_stress_coefficient_residual']:.3e}")
    print('New analysis directory:',opts.output_dir or 'output/general_observable/lo_pattern_study16/analysis_batch06_01')
    return 0

if __name__=='__main__':
    raise SystemExit(main())
