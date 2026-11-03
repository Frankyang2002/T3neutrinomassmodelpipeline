"""Chat A Batch 4: explain all 16 frozen neutrino spectra via heavy-state flavour geometry."""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from Numerical.lo_patterns_16.flavour_geometry import analyse,FlavourError


def main()->int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root',type=Path,default=ROOT)
    p.add_argument('--study-dir',type=Path)
    p.add_argument('--batch02-dir',type=Path)
    p.add_argument('--batch03-dir',type=Path)
    p.add_argument('--output-dir',type=Path,
                   help='Fresh directory inside output/general_observable/lo_pattern_study16')
    args=p.parse_args()
    root=args.project_root.resolve()
    try:
        result=analyse(root,**{k:(root/getattr(args,k) if getattr(args,k) is not None else None)
                               for k in ('study_dir','batch02_dir','batch03_dir','output_dir')})
    except (FlavourError,ValueError,TypeError,KeyError,RuntimeError,OSError) as exc:
        print(f'{type(exc).__name__}: {exc}',file=sys.stderr)
        return 2
    print(f'Batch 4: {result["model_count"]}/16 conditional flavour decompositions; physically certified 0/16')
    for x in result['models']:
        l=x['lightest_mode'];c=l['coherence_fraction']
        print(f'{x["model"]:<5} {x["fermion_type"]:<8} groups={x["number_heavy_mass_groups"]} '
              f'm1/m3={x["m1_over_m3"]:.5g} A1={l["geometric_projection_fraction"]:.5g} '
              f'C1={c if c is not None else "unresolved"} '
              f'mechanism={l["mechanism_diagnostic"]}')
    print('New output:',args.output_dir or 'output/general_observable/lo_pattern_study16/analysis_batch04_01')
    return 0


if __name__=='__main__':raise SystemExit(main())
