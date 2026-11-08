"""Generate a 1-loop Weyl contraction topology audit for Batch 3--4.

Topology is independent of published beta-coefficient implementations.
It does NOT certify A(-4) physics or replace a full counterterm derivation.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0,str(ROOT))
from RGE.running.intermediate.GeneralEFT1OneLoopTopologyAudit import classify


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path)
    args=p.parse_args(argv)
    report=classify()
    if args.output:
        if args.output.exists():
            raise FileExistsError('Refusing overwrite: '+str(args.output))
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf8')
    for k,v in report['sectors'].items():
        print(f"{k}: chirality-allowed connected contractions={v['number_of_wick_contractions']}, 1PI={v['number_1pi_wick_contractions']}")
    print(report['status'])
    print('Physical A(-4) certification: NO')
    return 0

if __name__=='__main__':
    raise SystemExit(main())
