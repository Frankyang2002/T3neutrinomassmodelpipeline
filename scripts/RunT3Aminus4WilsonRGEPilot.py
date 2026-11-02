"""Audit real-linear master-RGE structure on the successful A_-4 v3 boundary.

No physical Gamma is invented. Optionally solve an independently supplied,
exact-basis A C + B C* source for conditional, not physical evolution.

  python -X utf8 -m scripts.RunT3Aminus4WilsonRGEPilot --mode audit
"""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

from Numerical.running.GeneralT3WilsonEvolution import (
    aggregate_source_coordinates, boundary_coordinates, digest,
    read_json, save_json, tag, verify_batch,
)
from Numerical.running.GeneralT3WilsonRGEPilot import audit_a_minus_4
from Numerical.running.GeneralT3RealLinearWilsonFlow import evolve_real_linear

ROOT=Path(__file__).resolve().parents[1]
DEFAULT_BOUNDARY=Path('output/general_observable/ordered_wilson_tensor_v3')
DEFAULT_OUTPUT=Path('output/general_observable/a_minus4_wilson_rge_pilot_v1')


def absolute(path:Path)->Path:
    return path.resolve() if path.is_absolute() else (ROOT/path).resolve()


def main(argv:list[str]|None=None)->int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=('audit','evolve'), default='audit')
    parser.add_argument('--boundary-root',type=Path,default=DEFAULT_BOUNDARY)
    parser.add_argument('--output-root',type=Path,default=DEFAULT_OUTPUT)
    parser.add_argument('--gamma',type=Path,
                        help='Required for conditional evolution; explicit complete A and B samples')
    parser.add_argument('--points',type=int,default=41)
    args=parser.parse_args(argv)
    if args.mode=='evolve' and args.gamma is None:
        parser.error('--mode evolve requires --gamma: no A/B anomalous dimensions inferred')
    try:
        root=absolute(args.boundary_root)
        inspected=verify_batch(root)
        info=inspected['A_-4']
        path=root/tag('A_-4')/'ordered_tensor.json'
        source=read_json(path)
        audit=audit_a_minus_4(source)
        audit['ordered_tensor_sha256']=digest(path)
        out=absolute(args.output_root)
        save_json(out/'A_m4'/'basis_and_master_rge_audit.json', audit)
        print('A_-4:',audit['status'],
              '; ordered coordinates:',audit['summed_ordered_flavour_coordinate_count'],flush=True)
        if args.mode=='evolve':
            labels,values=aggregate_source_coordinates(*boundary_coordinates(source,'A_-4')[:2])
            gamma=read_json(absolute(args.gamma))
            trajectory=evolve_real_linear(labels,values,gamma,
                                          model='A_-4',tensor_sha256=digest(path),
                                          mu_f=float(info['mu_F_gev']),points=args.points)
            save_json(out/'A_m4'/'conditional_real_linear_trajectory.json',trajectory)
            print('A_-4:',trajectory['status'],
                  '; relative change:',trajectory['relative_change'],flush=True)
        else:
            print('No A/B anomalous-dimension samples supplied; physical EFT1 running NOT CERTIFIED.',flush=True)
        return 0
    except (OSError,ValueError,TypeError,KeyError,RuntimeError) as exc:
        print('BLOCKED: '+f'{type(exc).__name__}: {exc}',file=sys.stderr)
        return 2

if __name__=='__main__':
    raise SystemExit(main())
