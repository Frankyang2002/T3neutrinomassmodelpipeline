#!/usr/bin/env python3
"""Audit native portal CG operators and all-16 frozen scalar benchmark inputs.

Requires prior native Wolfram tensor export. Never changes an input file.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from Numerical.lo_certification.direct_uv.native_portal_dictionary import (  # noqa:E402
    audit, PortalDictionaryError,
)


def main()->int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root',type=Path,default=ROOT)
    p.add_argument('--native-json',type=Path,required=True)
    p.add_argument('--output-dir',type=Path,
        default=Path('output/general_observable/lo_certification_direct_uv/native_portal_dictionary_01'))
    p.add_argument('--no-frozen-inputs',action='store_true',
        help='Only verify exact native portal tensor identities; skip all-16 frozen input blocks')
    a=p.parse_args()
    root=a.project_root.resolve()
    native=a.native_json if a.native_json.is_absolute() else root/a.native_json
    dest=a.output_dir if a.output_dir.is_absolute() else root/a.output_dir
    path=dest/'native_portal_operator_dictionary.json'
    try:
        if path.exists():
            raise PortalDictionaryError(f'Refusing to overwrite {path}')
        result=audit(root,native,not a.no_frozen_inputs)
        dest.mkdir(parents=True,exist_ok=True)
        path.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
        print('Live registered portal CG classes verified: 5/5')
        print('Basic RGBeta singlet in native real-Hermitian portal span: 10/10 operator sectors')
        print('Unique native portal coupling dictionaries: '
              f'{10 - len(result["nonunique_operator_subspaces"])}/10')
        for row in result['nonunique_operator_subspaces']:
            print(f'  NON-UNIQUE {row["class"]}/{row["scalar"]}: '
                  f'rank {row["real_Hermitian_rank"]}/'
                  f'{row["native_invariant_count"]}; '
                  f'zero Hermitian projections {row["zero_projected_indices"]}')
        if result['frozen_input_audit']:
            print('Frozen scalar-A basic portal blocks: 16/16 verified')
        print('Full physical Matchete/RGBeta-to-Weyl normalisation: 0/16 certified')
        print('Report:',path)
        return 0
    except (PortalDictionaryError,OSError,KeyError,ValueError,TypeError,
            json.JSONDecodeError) as exc:
        print('Native portal dictionary BLOCKED: '+str(exc),file=sys.stderr)
        return 2


if __name__=='__main__':
    raise SystemExit(main())
