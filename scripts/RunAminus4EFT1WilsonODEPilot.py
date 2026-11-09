"""Run a strictly PROVISIONAL A(-4) master-Weyl Wilson ODE pilot.

Input --source: formal SourceBoundary NPZ + same-stem JSON from
Numerical.joint_scalar_matching.fermion_boundary_bridge.write_provisional_source.

Input --couplings: A(-4) sampled EFT1 tensors NPZ + nearby basis JSON from
Numerical.running.Aminus4EFT1Couplings.write_sampled_npz.

No physical neutrino mass is calculated. The explicit --formal-identity-map
switch admits an unverified basis map and spinor factor ONLY for diagnostics.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import numpy as np

from Numerical.running.T3EFT1WilsonODETransport import transport


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--couplings', type=Path, required=True)
    parser.add_argument('--couplings-json', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--formal-identity-map', action='store_true')
    parser.add_argument('--basis-map-npy', type=Path)
    parser.add_argument('--signed-spinor-factor', type=float)
    parser.add_argument('--rtol', type=float, default=1e-6)
    parser.add_argument('--atol', type=float, default=1e-18)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Output exists: refusing to overwrite historical results')
    if args.formal_identity_map and args.basis_map_npy is not None:
        parser.error('Formal identity flag and explicit map are mutually exclusive')
    R = None if args.basis_map_npy is None else np.load(args.basis_map_npy,
                                                       allow_pickle=False)
    report, C = transport(
        args.source, args.couplings,
        couplings_json=args.couplings_json,
        heavy_source_from_rgbeta=R,
        spinor_signed_factor=args.signed_spinor_factor,
        allow_formal_identity=args.formal_identity_map,
        rtol=args.rtol, atol=args.atol)
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output/'pilot_summary.json').write_text(
        json.dumps(report, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    np.savez_compressed(args.output/'formal_master_C_at_muS.npz', C=C)
    print(report['status'])
    print(args.output/'pilot_summary.json')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
