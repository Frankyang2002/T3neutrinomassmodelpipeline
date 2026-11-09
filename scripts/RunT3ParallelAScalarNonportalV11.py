"""Audit mathematical SU2 nonportal invariant rank and special-quartic actions."""
from __future__ import annotations
import argparse,json,sys
from pathlib import Path
from Numerical.parallel_t3.wilson.scalar_nonportal_rank_v11 import V11Blocked,audit

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root',type=Path,default=Path('.'))
    p.add_argument('--samples',type=int,default=128)
    p.add_argument('--seed',type=int,default=21551)
    p.add_argument('--output',type=Path,default=Path('output/parallel_t3/wilson/nonportal_scalar_invariants_conditional_v11.json'))
    a=p.parse_args(argv);root=a.project_root.resolve()
    out=(a.output if a.output.is_absolute() else root/a.output).resolve()
    try:
        if out.exists():raise V11Blocked('Refusing to overwrite '+str(out))
        data=audit(seed=a.seed,samples=a.samples,root=root)
        out.parent.mkdir(parents=True,exist_ok=True)
        with out.open('x',encoding='utf-8') as f:
            json.dump(data,f,indent=2,allow_nan=False);f.write('\n')
        print('v11: {} source quartic registrations, ranks {}, special action residual {:.3e}'.format(
            data['source_registration_slots_total'],data['rank_sector_counts'],
            data['max_special_constructed_action_error']))
        print(data['status'])
        print('Independent mathematical SU2 frames, NOT numerical native RGBeta or Matchete normalization.')
        print(out)
        return 0
    except (V11Blocked,ValueError,KeyError,IndexError,TypeError,OSError,ArithmeticError) as e:
        print('Workstream A v11 BLOCKED: '+str(e),file=sys.stderr)
        return 2
if __name__=='__main__':raise SystemExit(main())
