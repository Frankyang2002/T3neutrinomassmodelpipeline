"""Native doublet Fierz and conditional triplet Lie-sign portal audit (v6).

Use --matchete, --rgbeta, --v5, --fierz and a new --output path.
Does not change any source, matching, running or prior numerical result.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

from Numerical.parallel_t3.wilson.portal_orientation_v6 import audit_from_paths
from Numerical.parallel_t3.wilson.matchete_rgbeta_dictionary import DictionaryEvidenceError


def main() -> int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--matchete',type=Path,required=True)
    p.add_argument('--rgbeta',type=Path,required=True)
    p.add_argument('--v5',type=Path,required=True)
    p.add_argument('--fierz',type=Path,required=True)
    p.add_argument('--project-root',type=Path,default=Path('.'))
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    try:
        root=args.project_root.resolve()
        path=lambda x: x.resolve() if x.is_absolute() else (root/x).resolve()
        output=path(args.output)
        if output.exists():
            raise DictionaryEvidenceError('Refusing to overwrite '+str(output))
        report=audit_from_paths(path(args.matchete),path(args.rgbeta),
                                path(args.v5),path(args.fierz),root)
        output.parent.mkdir(parents=True,exist_ok=True)
        with output.open('x',encoding='utf-8') as f:
            json.dump(report,f,indent=2,allow_nan=False)
            f.write('\n')
        print(f"{report['doublet_sectors_Fierz_verified']}/12 doublet Fierz, "
              f"{report['triplet_sectors_Pauli_algebra_conditional']}/14 "
              f"triplet Pauli-Lie, {report['singlet_sectors']}/6 singlets; "
              f"residual={report['max_relative_component_and_potential_residual']:.3e}")
        print(report['status'],'->',output)
        print('NOT native RGBeta S2 component-frame or complete scalar normalization.')
        return 0
    except (DictionaryEvidenceError,ValueError,TypeError,KeyError,
            OSError,ArithmeticError,json.JSONDecodeError) as exc:
        print('Portal v6 BLOCKED:',exc,file=sys.stderr)
        return 2


if __name__=='__main__':
    raise SystemExit(main())
