"""Stage17: B/C(-1) extra-quartic EFT1 -> Weinberg operator mixing audit."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from Numerical.experimental.GeneralT3ExtraEFT1Mixing import run_file


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output-root',type=Path,default=Path('output'))
    p.add_argument('--stage15',type=Path,default=Path('output/general_t3_stage15/source_lagrangian_convention.json'))
    p.add_argument('--out',type=Path,default=Path('output/general_t3_stage17'))
    args=p.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    items=[]; failures=[]
    for model, folder in [('B_m1','T3_B_alpha_m1'),('C_m1','T3_C_alpha_m1')]:
        base=args.output_root/folder/'data'
        try:
            result=run_file(model,base/'eft1_after_F_wilson_seed.json',base/'eft1_rgbeta_rge.json',args.stage15)
            path=args.out/f'{model}_extra_eft1_scalar_pair_beta.json'
            path.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
            print(f'{model}: {result["status"]}, {result["nonzero_weinberg_real_components"]} nonzero Higgs-real components')
            print('Report:',path)
            items.append({'model':model,'status':result['status'],
                          'nonzero_weinberg_real_components':result['nonzero_weinberg_real_components']})
        except Exception as e:
            # Leave the other model running and write the concrete failure.
            print(f'BLOCKED {model}: {type(e).__name__}: {e}',file=sys.stderr)
            failures.append({'model':model,'type':type(e).__name__,'error':str(e)})
    summary={'status':'ExtraScalarPairBetasDerived' if not failures else 'ReviewRequired',
             'models':items,'failures':failures,'complete_extra_eft1_matching':False,
             'production_code_modified':False}
    (args.out/'stage17_summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    return 0 if not failures else 2


if __name__=='__main__':
    raise SystemExit(main())
