"""Export exact sparse Cartesian quartic candidates from native Matchete evidence.

This is NOT an RGBeta normalization certificate and cannot unlock production.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from Numerical.running.T3MatcheteCartesianResponses import convert_matchete_export

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REG = ROOT/'RGE/general/generated/general_rgbeta_quartic_registry.json'
DEFAULT_SOURCE = ROOT/'output/general_t3_stage22/matchete_tensors.json'
DEFAULT_OUT = ROOT/'output/general_t3_stage22/matchete_cartesian_candidates.json'


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--matchete',type=Path,default=DEFAULT_SOURCE)
    p.add_argument('--registry',type=Path,default=DEFAULT_REG)
    p.add_argument('--output',type=Path,default=DEFAULT_OUT)
    args = p.parse_args(argv)
    if not args.matchete.is_file():
        p.error(f'Native Matchete export missing: {args.matchete}. Run wolframscript exporter first.')
    source_bytes = args.matchete.read_bytes()
    registry_bytes = args.registry.read_bytes()
    payload = convert_matchete_export(json.loads(source_bytes),json.loads(registry_bytes))
    payload['source_sha256'] = hashlib.sha256(source_bytes).hexdigest()
    payload['registry_sha256'] = hashlib.sha256(registry_bytes).hexdigest()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')
    print('Matchete Cartesian candidate response export: 16 models; 50 quartic slots')
    print('WARNING: NOT a certified RGBeta normalization map; do not use as production tensor')
    print('Output:',args.output)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
