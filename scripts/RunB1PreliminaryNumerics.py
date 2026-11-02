"""Numerical B_-1 one-channel conditional source-backed EFT1 -> joint S pilot."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path: sys.path.insert(0,str(ROOT))
from RGE.running.intermediate.B1PreliminaryNumerics import produce

SAFE=ROOT/'output/general_observable/lo_eft1_mixing'
DEFAULT=SAFE/'b1_lambdaT3_source_numerics_v1'

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo-root',type=Path,default=ROOT)
    p.add_argument('--output-root',type=Path,default=DEFAULT)
    p.add_argument('--check-only',action='store_true',help='Compute and print summary without writing')
    x=p.parse_args(argv)
    root=x.repo_root.resolve()
    dest=x.output_root.resolve()
    safe=(root/'output/general_observable/lo_eft1_mixing').resolve()
    if dest==safe or safe not in dest.parents: p.error(f'Writes only to NEW subfolder of {safe}')
    final=dest/'b1_lambdaT3_one_loop.json'
    if not x.check_only and (dest.exists() or final.exists()):
        p.error('Refusing to reuse output directory; choose a fresh --output-root')
    data=produce(root)  # all read-only validation occurs before any writes
    if x.check_only:
        print(json.dumps(data,indent=2,allow_nan=False))
    else:
        dest.mkdir(parents=True,exist_ok=False)
        with final.open('x',encoding='utf-8') as h:
            json.dump(data,h,indent=2,allow_nan=False)
            h.write('\n')
        print('Wrote conditional B_-1 source numerical contribution to',final)
    print('SU2 Ward pass:',data['response']['full_SU2_Weinberg_projection_pass'])
    print('Physical C5 computed: false; unknown GammaCC/Matchete-RGBeta mapping')
    return 0 if data['response']['full_SU2_Weinberg_projection_pass'] else 2

if __name__=='__main__':
    raise SystemExit(main())
