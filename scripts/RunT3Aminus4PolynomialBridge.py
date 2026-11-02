"""Read-only A_-4 polynomial-coordinate bridge for the 16-model v3 boundary.

Default output: output/general_observable/a_minus4_polynomial_bridge_v1
No previous matching, ordered tensor or RGE audit file is overwritten.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

from Numerical.running.GeneralT3WilsonEvolution import (
    verify_batch, read_json, save_json, digest, tag,
)
from Numerical.running.GeneralT3WilsonPolynomialBridge import bridge_a4

ROOT=Path(__file__).resolve().parents[1]


def absolute(path:Path)->Path:
    return path.resolve() if path.is_absolute() else (ROOT/path).resolve()


def main(argv:list[str]|None=None)->int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--boundary-root',type=Path,
        default=Path('output/general_observable/ordered_wilson_tensor_v3'))
    p.add_argument('--output-root',type=Path,
        default=Path('output/general_observable/a_minus4_polynomial_bridge_v1'))
    p.add_argument('--skip-master-example-check',action='store_true',
        help='Only for isolated fixture tests without the repository RGE package')
    args=p.parse_args(argv)
    try:
        root=absolute(args.boundary_root)
        verify_batch(root)
        source_file=root/tag('A_-4')/'ordered_tensor.json'
        source=read_json(source_file)
        source_ref=Path(source['source_boundary'])
        source_ref=absolute(source_ref)
        if not source_ref.is_file() or digest(source_ref)!=source.get('source_boundary_sha256'):
            raise ValueError('Stale/missing original A_-4 Matchete flavour source')
        source_boundary=read_json(source_ref)
        report=bridge_a4(source,source_boundary,check_master=not args.skip_master_example_check)
        report['source_ordered_tensor_sha256']=digest(source_file)
        report['original_matchete_boundary_sha256']=digest(source_ref)
        out=absolute(args.output_root)/'A_m4'/'polynomial_bridge.json'
        if out.resolve()==source_file.resolve():
            raise ValueError('Refusing to overwrite matching boundary')
        save_json(out,report)
        print('A_-4:',report['status'])
        print('Ordered:',report['source_ordered_coordinate_count'],
              'symmetric polynomial entries:',report['symmetric_nonzero_coordinate_count'])
        print('Flavour-axis permutation:',
              report['lepton_dummy_alignment']['flavour_axis_for_lepton_field_order_zero_based'])
        print('Polynomial relative residual:',report['polynomial_identity']['relative_monomial_difference'])
        print('Master Weinberg reconstruction:',
              report['master_example_convention']['master_reconstruction_residual']
              if report['master_example_convention'] else 'not checked')
        print('Physical Matchete normalization / EFT1 RGE: NOT CERTIFIED')
        print('Output:',out)
        return 0
    except (OSError,ValueError,TypeError,KeyError,RuntimeError) as exc:
        print('BLOCKED:',f'{type(exc).__name__}: {exc}',file=sys.stderr)
        return 2

if __name__=='__main__':
    raise SystemExit(main())
