"""Chat A batch 3: explain frozen direct-UV patterns using an exact analytic kernel."""
from __future__ import annotations
import argparse
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from Numerical.lo_patterns_16.analytic_explanation import explain,ExplanationError

def main()->int:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--project-root',type=Path,default=ROOT)
    ap.add_argument('--study-dir',type=Path,default=None)
    ap.add_argument('--batch02-dir',type=Path,default=None)
    ap.add_argument('--output-dir',type=Path,default=None)
    a=ap.parse_args()
    root=a.project_root.resolve()
    try:
        result=explain(root,study_dir=root/a.study_dir if a.study_dir else None,
                        batch02_dir=root/a.batch02_dir if a.batch02_dir else None,
                        output_dir=root/a.output_dir if a.output_dir else None)
    except (ExplanationError,ValueError,TypeError,KeyError,RuntimeError,OSError) as ex:
        print(f'{type(ex).__name__}: {ex}',file=sys.stderr)
        return 2
    print(f'Analytic scalar-kernel reconstructions: {result["model_count"]}/16; physically certified: 0/16')
    for x in result['models']:
        print(f'{x["model"]:<5} charges={len(x["channels"])} '
              f'm1/m3={x["m1_over_m3"]:.4g} m2/m3={x["m2_over_m3"]:.4g} '
              f'exact relative error={x["exact_reconstruction_relative_error"]:.3g} '
              f'leading error={x["small_mixing_relative_error"]:.3g} '
              f'charge cosine={x["charge_cosine"]}')
    print('New output: ', a.output_dir or 'output/general_observable/lo_pattern_study16/analysis_batch03_01')
    return 0

if __name__=='__main__':raise SystemExit(main())
