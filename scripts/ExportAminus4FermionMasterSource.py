"""Export the *actual A(-4) F-threshold* source in the existing FORMAL format.

The UV-state JSON must contain MF, y1 and y2, each a numerical 3x3 matrix
(as nested lists or {"real":matrix,"imag":matrix}). Use the state evaluated at
mu_F, NOT a Yukawa state from some other scale. Relevant Matchete seed exports
must already exist in the local repository's standard F_ROOTS directories.

This does not identify Matchete GammaCC with the master Weyl convention.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from Numerical.joint_scalar_matching.fermion_boundary_bridge import write_provisional_source
from Numerical.running.GeneralT3DynamicFermionWilson import (
    nonmixed_boundary, uv_from_json,
)

ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--uv-state', type=Path, required=True)
    ap.add_argument('--mu-f', type=float, required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    if args.output.suffix != '.npz':
        ap.error('--output must have .npz extension')
    if args.output.exists() or args.output.with_suffix('.json').exists():
        ap.error('Output already exists: refusing to overwrite')
    uv_values = json.loads(args.uv_state.read_text(encoding='utf-8-sig'))
    state = uv_from_json('A_-4', uv_values, mu_f=args.mu_f)
    source, audit = nonmixed_boundary('A_-4', state, repo_root=ROOT,
                                      mu_f=args.mu_f)
    write_provisional_source(args.output, source)
    report = {
        'model':'A_-4', 'status':'ActualMatchedFSource_AbsoluteWeylNormalizationOpen',
        'mu_F_gev':args.mu_f,'source_tensor_shape':list(source.raw.shape),
        'matching_audit':audit,
        'full_physical_Wilson_boundary_certified':False,
    }
    (args.output.parent / (args.output.stem + '_audit.json')).write_text(
        json.dumps(report, indent=2, default=str) + '\n', encoding='utf8')
    print('Actual A(-4) fermion matching source exported as FORMAL tensor:', args.output)
    print('The signed GammaCC/Weyl factor remains unresolved.')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
