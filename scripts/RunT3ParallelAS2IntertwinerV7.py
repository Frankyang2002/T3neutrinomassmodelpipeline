"""Read-only conditional S2 audit. Run as python -m scripts.RunT3ParallelAS2IntertwinerV7."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys

from Numerical.parallel_t3.wilson.s2_intertwiner_v7 import (
    S2EvidenceError,audit_from_paths)


def main()->int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root',type=Path,default=Path('.'))
    p.add_argument('--matchete',type=Path,default=Path('output/parallel_t3/wilson/matchete_portals_native_v5.json'))
    p.add_argument('--rgbeta',type=Path,default=Path('output/parallel_t3/wilson/rgbeta_portals_native_v5.json'))
    p.add_argument('--v5',type=Path,default=Path('output/parallel_t3/wilson/portal_dictionary_conditional_v5.json'))
    p.add_argument('--fierz',type=Path,default=Path('output/parallel_t3/wilson/rgbeta_fundamental_fierz_native_v6.json'))
    p.add_argument('--v6',type=Path,default=Path('output/parallel_t3/wilson/portal_orientation_conditional_v6.json'))
    p.add_argument('--native-rgbeta-witness',type=Path,default=None,
                   help='Optional historical B2 native SYMBOLIC witness; never a numeric frame')
    p.add_argument('--generator-candidate',type=Path,default=None,
                   help='Optional historical Matchete fixture provenance (not independent live data)')
    p.add_argument('--output',type=Path,
                   default=Path('output/parallel_t3/wilson/s2_intertwiner_conditional_v7.json'))
    args=p.parse_args()
    root=args.project_root.resolve()
    def resolve(x):return None if x is None else (x if x.is_absolute() else root/x).resolve()
    output=resolve(args.output)
    try:
        if output.exists():raise S2EvidenceError('Refusing to overwrite '+str(output))
        report=audit_from_paths(resolve(args.matchete),resolve(args.rgbeta),
                 resolve(args.v5),resolve(args.fierz),resolve(args.v6),root,
                 resolve(args.native_rgbeta_witness),resolve(args.generator_candidate))
        output.parent.mkdir(parents=True,exist_ok=True)
        with output.open('x',encoding='utf-8') as f:
            json.dump(report,f,indent=2,allow_nan=False)
            f.write('\n')
        print(f"v7: {report['triplet_sectors_conditional_comparison_verified']}/14 S2-to-Matchete conditional comparisons; "
              f"max residual {report['maximum_comparison_residual']:.3e}")
        print(report['status'])
        print('Native RGBeta numerical S2 frame NOT measured. No physics evolution authorized.')
        print(output)
        return 0
    except (S2EvidenceError,ValueError,TypeError,KeyError,IndexError,
            OSError,ArithmeticError) as exc:
        print('Workstream A v7 BLOCKED:',exc,file=sys.stderr)
        return 2

if __name__=='__main__':
    raise SystemExit(main())
