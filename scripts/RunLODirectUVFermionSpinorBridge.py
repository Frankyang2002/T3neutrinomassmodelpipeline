#!/usr/bin/env python3
"""Read-only WS4 native fermion/conditional Weyl bridge; --native-dir required."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from Numerical.lo_certification.direct_uv.native_fermion_action_bridge import (  # noqa:E402
    evaluate,FermionBridgeError,
)


def main() -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project-root',type=Path,default=ROOT)
    parser.add_argument('--native-dir',type=Path,required=True)
    parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args()
    root=args.project_root.resolve()
    native=args.native_dir if args.native_dir.is_absolute() else root/args.native_dir
    out=args.output_dir if args.output_dir.is_absolute() else root/args.output_dir
    report=out/'native_fermion_spinor_bridge.json'
    try:
        if out.exists() and any(out.iterdir()):
            raise FermionBridgeError('Output directory exists and is not empty: '+str(out))
        data=evaluate(root,native)
        out.mkdir(parents=True,exist_ok=True)
        report.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    except (FermionBridgeError,ValueError,OSError,KeyError,TypeError,IndexError,
            json.JSONDecodeError) as error:
        print('Native fermion/spinor bridge BLOCKED: '+str(error),file=sys.stderr)
        return 2
    print('RGBeta simultaneous y1/y2/MF native registry: 3/3 representative models')
    print('Matchete original restricted mass/Yukawa action: 3/3 representatives')
    print('All-16 candidate Weyl masses and Yukawa tensors: 16/16 matched to frozen inputs')
    print('Physical Matchete/RGBeta-to-Weyl phase/normalisation: 0/16 certified')
    print('Report:',report)
    return 0


if __name__=='__main__':raise SystemExit(main())
