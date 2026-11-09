"""Audit all16 declared-frame lambdaT3 topology complex scalar actions."""
from __future__ import annotations
import argparse,json,sys
from pathlib import Path
from Numerical.parallel_t3.wilson.scalar_topology_action_v12 import V12Blocked,audit

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root',type=Path,default=Path('.'))
    p.add_argument('--output',type=Path,default=Path('output/parallel_t3/wilson/topology_scalar_action_conditional_v12.json'))
    p.add_argument('--samples',type=int,default=5)
    a=p.parse_args(argv);root=a.project_root.resolve()
    out=(a.output if a.output.is_absolute() else root/a.output).resolve()
    try:
        if out.exists():raise V12Blocked('Refusing to overwrite '+str(out))
        result=audit(samples=a.samples,root=root)
        out.parent.mkdir(parents=True,exist_ok=True)
        with out.open('x',encoding='utf-8') as f:
            json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
        print('v12: {}/16 lambdaT3 complex source actions; worst {:.3e}'.format(
            result['lambdaT3_complex_quartic_slots_constructed'],
            result['max_constructed_complex_action_residual']))
        print(result['status']);print(out)
        print('No native RGBeta complex Lam / Matchete component-map / Wilson physics certification.')
        return 0
    except (V12Blocked,ValueError,KeyError,TypeError,OSError,ArithmeticError) as e:
        print('Workstream A v12 BLOCKED: '+str(e),file=sys.stderr);return 2
if __name__=='__main__':raise SystemExit(main())
