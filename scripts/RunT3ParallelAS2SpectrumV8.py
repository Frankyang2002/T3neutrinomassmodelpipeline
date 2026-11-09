"""Workstream A v8: source-guarded native S2 reduction and independent Casimir spectrum.

Execute from repository root: python -m scripts.RunT3ParallelAS2SpectrumV8
Read-only; requires fresh native symbolic witness and completed v7 output.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

from Numerical.parallel_t3.wilson.s2_spectral_v8 import audit_from_paths
from Numerical.parallel_t3.wilson.s2_intertwiner_v7 import S2EvidenceError


def main()->int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root',type=Path,default=Path('.'))
    p.add_argument('--matchete',type=Path,default=Path('output/parallel_t3/wilson/matchete_portals_native_v5.json'))
    p.add_argument('--rgbeta',type=Path,default=Path('output/parallel_t3/wilson/rgbeta_portals_native_v5.json'))
    p.add_argument('--v5',type=Path,default=Path('output/parallel_t3/wilson/portal_dictionary_conditional_v5.json'))
    p.add_argument('--fierz',type=Path,default=Path('output/parallel_t3/wilson/rgbeta_fundamental_fierz_native_v6.json'))
    p.add_argument('--v6',type=Path,default=Path('output/parallel_t3/wilson/portal_orientation_conditional_v6.json'))
    p.add_argument('--v7',type=Path,default=Path('output/parallel_t3/wilson/s2_intertwiner_conditional_v7.json'))
    p.add_argument('--native-v8',type=Path,default=Path('output/parallel_t3/wilson/rgbeta_s2_reduction_native_v8.json'))
    p.add_argument('--output',type=Path,default=Path('output/parallel_t3/wilson/s2_spectrum_conditional_v8.json'))
    args=p.parse_args()
    root=args.project_root.resolve()
    resolve=lambda path:(path if path.is_absolute() else root/path).resolve()
    output=resolve(args.output)
    try:
        if output.exists():
            raise S2EvidenceError('Refusing to overwrite '+str(output))
        report=audit_from_paths(resolve(args.matchete),resolve(args.rgbeta),
                 resolve(args.v5),resolve(args.fierz),resolve(args.v6),
                 resolve(args.v7),resolve(args.native_v8),root)
        output.parent.mkdir(parents=True,exist_ok=True)
        with output.open('x',encoding='utf-8') as fp:
            json.dump(report,fp,indent=2,allow_nan=False)
            fp.write('\n')
        print(f"v8: {report['triplet_sectors_basis_invariant_spectrum_verified']}/14 "
              f"triplet spectra, maximum residual "
              f"{report['maximum_spectral_and_v6_comparison_residual']:.3e}")
        print(report['status'])
        print('Native numerical RGBeta S2 frame NOT certified; NO physical Wilson evolution.')
        print(output)
        return 0
    except (S2EvidenceError,ValueError,TypeError,KeyError,IndexError,
            OSError,ArithmeticError) as exc:
        print('Workstream A v8 BLOCKED:',exc,file=sys.stderr)
        return 2


if __name__=='__main__':
    raise SystemExit(main())
