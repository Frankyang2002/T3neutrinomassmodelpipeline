"""Run Workstream A v9 portal real-action algebra audit (read-only / conditional).

Native evidence v5-v8 must already be present; never launches physical running.
Only writes a new v9 JSON. Refuses to overwrite existing files.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

from Numerical.parallel_t3.wilson.portal_real_action_v9 import (
    audit_from_paths, V9EvidenceError,
)


def main(argv=None) -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project-root',type=Path,default=Path('.'))
    defaults={
        'matchete':'matchete_portals_native_v5.json',
        'rgbeta':'rgbeta_portals_native_v5.json',
        'v5':'portal_dictionary_conditional_v5.json',
        'fierz':'rgbeta_fundamental_fierz_native_v6.json',
        'v6':'portal_orientation_conditional_v6.json',
        'v7':'s2_intertwiner_conditional_v7.json',
        'native-v8':'rgbeta_s2_reduction_native_v8.json',
        'v8':'s2_spectrum_conditional_v8.json',
        'output':'portal_real_action_conditional_v9.json',
    }
    for key,value in defaults.items():
        parser.add_argument('--'+key,type=Path,default=Path('output/parallel_t3/wilson')/value)
    parser.add_argument('--samples',type=int,default=5)
    args=parser.parse_args(argv)
    root=args.project_root.resolve()
    path=lambda p: (p if p.is_absolute() else root/p).resolve()
    out=path(args.output)
    try:
        if out.exists():raise V9EvidenceError('Refusing to overwrite '+str(out))
        report=audit_from_paths(path(args.matchete),path(args.rgbeta),
            path(args.v5),path(args.fierz),path(args.v6),path(args.v7),
            path(args.native_v8),path(args.v8),root,samples=args.samples)
        out.parent.mkdir(parents=True,exist_ok=True)
        with out.open('x',encoding='utf-8') as f:
            json.dump(report,f,indent=2,allow_nan=False)
            f.write('\n')
        print(f"v9: {report['native_Matchete_portal_coupling_slots']}/58 portal slots, "
              f"{report['triplet_sectors_spectral_normalization_verified']}/14 triplets, "
              f"max conditional real-action residual "
              f"{report['max_constructed_real_action_relative_residual']:.3e}")
        print(report['status'])
        print('Native RGBeta Lam projector NOT executed; physical Wilson evolution BLOCKED.')
        print(out)
        return 0
    except (V9EvidenceError,ValueError,TypeError,KeyError,IndexError,
            OSError,ArithmeticError) as e:
        print('Workstream A v9 BLOCKED: '+str(e),file=sys.stderr)
        return 2


if __name__=='__main__':
    raise SystemExit(main())
