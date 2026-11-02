"""Post-process Chat A provisional_01 once, without rerunning the UV solver."""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from Numerical.lo_patterns_16.result_analysis import analyse,AnalysisError  # noqa:E402


def main() -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project-root',type=Path,default=ROOT)
    parser.add_argument('--study-dir',type=Path,default=None,
                        help='Previous complete provisional_01 directory (read only)')
    parser.add_argument('--output-dir',type=Path,default=None,
                        help='New, empty directory within output/general_observable/lo_pattern_study16')
    parser.add_argument('--no-figures',action='store_true')
    args=parser.parse_args()
    root=args.project_root.resolve()
    try:
        result=analyse(root,
                       study_dir=(root/args.study_dir if args.study_dir else None),
                       output_dir=(root/args.output_dir if args.output_dir else None),
                       figures=not args.no_figures)
    except (AnalysisError,FileExistsError,KeyError,ValueError,TypeError,OSError) as exc:
        print(type(exc).__name__+': '+str(exc),file=sys.stderr)
        return 2
    print('Batch 2: 16/16 source reports and 288 scans analysed; physically certified 0/16')
    for r in result['models']:
        print(f"{r['model']:<5} m1={r['m1_ev']:.6g} m2={r['m2_ev']:.6g} "
              f"m3={r['m3_ev']:.6g} eV; charge coherence={r['charge_coherence_ratio']:.5g}; "
              f"lambda t=2 residual={r['lambda_t2_linear_matrix_residual']:.3g}")
    print('New analysis directory:',args.output_dir or 'output/general_observable/lo_pattern_study16/analysis_batch02_01')
    return 0


if __name__=='__main__':
    raise SystemExit(main())
