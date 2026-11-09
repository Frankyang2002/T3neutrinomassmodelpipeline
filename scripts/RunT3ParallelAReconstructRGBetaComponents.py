"""Build a CONDITIONAL scalar-invariant component reference from RGBeta native declarations.

Never identifies Matchete and RGBeta axes or coupling normalizations.
Use: python -m scripts.RunT3ParallelAReconstructRGBetaComponents --output FILE.json
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from Numerical.parallel_t3.wilson.rgbeta_components import load_reconstruct


def main(argv=None) -> int:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--declarations',type=Path,
       default=Path('output/general_t3_stage22/rgbeta_native_declarations.json'))
    ap.add_argument('--registry',type=Path,
       default=Path('RGE/general/generated/general_rgbeta_quartic_registry.json'))
    ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args(argv)
    if a.output.exists():
        ap.error('Output exists; use a fresh versioned path (no silent overwrite)')
    result=load_reconstruct(a.declarations,a.registry)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(f'{result["quartic_slot_count"]} / 50 conditional components across '
          f'{result["model_count"]} models -> {a.output}')
    print('NOT normalized; epsilon and S2 component frames are chosen, not native certified.')
    return 0


if __name__=='__main__':
    raise SystemExit(main())
