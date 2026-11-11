"""Check actual 16-model Matchete candidate/registry mapping without physics runs."""
import argparse
from pathlib import Path
from Numerical.running.T3ExploratoryFullQuartic import build_overlay
ROOT=Path(__file__).resolve().parents[1]

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--candidate',type=Path,default=ROOT/'output/general_t3_stage22/matchete_cartesian_candidates.json')
    parser.add_argument('--registry',type=Path,default=ROOT/'RGE/general/generated/general_rgbeta_quartic_registry.json')
    args=parser.parse_args()
    import json
    from scripts import RunT3SimultaneousEFT1All16 as general
    reg=json.loads(args.registry.read_text(encoding='utf8'))
    total=0
    for model in sorted(reg['models']):
        names={x['name'] for x in reg['models'][model]['new_rgbeta_symbols']}
        names|={x['existing_name'] for x in reg['models'][model]['preexisting_rgbeta_symbols']}
        def q_provider(_state):
            raise RuntimeError('This is a structural mapping check, not an RGE run')
        basis=build_overlay(model,args.candidate,args.registry,general.general_quartic,q_provider,
                            required_quartics=names)
        count=len(basis.directions)//2
        total+=count
        print(f'{model:5} exploratory additional non-HHSS = {count}')
    if total!=26:
        raise RuntimeError(f'Expected 26 non-HHSS Matchete slots; found {total}')
    print('PASS: 16 models, 26 non-HHSS slots, explicit provisional mapping; no production running changed')
if __name__=='__main__':main()
