"""Batch 5 analytical source-tensor, Yukawa-texture, and heavy-mode synthesis."""
from __future__ import annotations
import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from Numerical.lo_patterns_16.texture_synthesis import SynthesisError, analyse


def main() -> int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root',type=Path,default=ROOT)
    p.add_argument('--study-dir',type=Path)
    p.add_argument('--batch03-dir',type=Path)
    p.add_argument('--batch04-dir',type=Path)
    p.add_argument('--output-dir',type=Path,help='Fresh path within output/general_observable/lo_pattern_study16')
    p.add_argument('--no-plots',action='store_true')
    args=p.parse_args()
    root=args.project_root.resolve()
    try:
        result=analyse(root,
            **{n:(root/getattr(args,n) if getattr(args,n) is not None else None)
               for n in ('study_dir','batch03_dir','batch04_dir','output_dir')},
            make_plots=not args.no_plots)
    except (SynthesisError,ValueError,TypeError,KeyError,RuntimeError,OSError) as exc:
        print(f'{type(exc).__name__}: {exc}',file=sys.stderr)
        return 2
    print(f'Batch 5: {result["model_count"]}/16 conditional texture+phasor analyses; physical certification 0/16')
    for entry in result['models']:
        obs=entry['observed'];contr=entry['common_mixing_controlled']
        print(f'{entry["model"]:<5} {entry["fermion_type"]:<8} '
              f'orig m1/m3={obs["m1_over_m3"]:.5g} '
              f'C1={obs["m1_coherence_C"]:.5g} '
              f'C_+1 texture m1/m3={contr["C_+1"]["m1_over_m3"]:.5g} '
              f'C_-1 texture m1/m3={contr["C_-1"]["m1_over_m3"]:.5g}')
    print('New output:',args.output_dir or 'output/general_observable/lo_pattern_study16/analysis_batch05_01')
    return 0


if __name__=='__main__':
    raise SystemExit(main())
