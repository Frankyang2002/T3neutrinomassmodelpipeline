#!/usr/bin/env python3
"""Run read-only all-16 native-source MF Ward and conditional phase-spurion audit.

python scripts/RunLODirectUVSpinorMassDictionary.py --project-root . \\
  --output-dir output/general_observable/lo_certification_direct_uv/spinor_mass_01
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0,str(ROOT))
from Numerical.lo_certification.direct_uv.spinor_mass_dictionary import audit, DictionaryError


def main()->int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root',type=Path,default=ROOT)
    p.add_argument('--output-dir',type=Path,default=None)
    a=p.parse_args()
    root=a.project_root.resolve()
    allowed=(root/'output/general_observable/lo_certification_direct_uv').resolve()
    output=(a.output_dir or allowed/'spinor_mass_01').resolve()
    try:
        if not output.is_relative_to(allowed):
            raise DictionaryError('Output outside the dedicated owned tree')
        path=output/'all_16_spinor_mass_dictionary.json'
        if path.exists():
            raise DictionaryError('Refusing overwrite: '+str(path)+'; choose a new --output-dir')
        result=audit(root)
        output.mkdir(parents=True,exist_ok=True)
        path.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    except (DictionaryError,OSError,ValueError,KeyError,TypeError) as exc:
        print('Source mass/spinor audit BLOCKED: '+str(exc),file=sys.stderr)
        return 2
    print('All-16 source mass Ward and conditional spurion identities: PASS')
    print('Dirac: {Dirac}/16; Majorana: {Majorana}/16; optional extra Yukawas: {extra_yukawa_allowed}/16'.format(**result['model_counts']))
    print('Physical 4-spinor/Weyl/RGBeta dictionary certified: 0/16')
    print('Report: '+str(path))
    return 0

if __name__=='__main__':raise SystemExit(main())
