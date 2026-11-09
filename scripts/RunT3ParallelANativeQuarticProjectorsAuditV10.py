"""Audit source-guarded native RGBeta Lam/Projector response (v10)."""
from __future__ import annotations
import argparse,json,sys
from pathlib import Path
from Numerical.parallel_t3.wilson.native_quartic_projectors_v10 import (
    V10Blocked,audit_from_paths,
)

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root',type=Path,default=Path('.'))
    p.add_argument('--native',type=Path,default=Path('output/parallel_t3/wilson/rgbeta_quartic_projectors_native_v10.json'))
    p.add_argument('--v9',type=Path,default=Path('output/parallel_t3/wilson/portal_real_action_conditional_v9.json'))
    p.add_argument('--output',type=Path,default=Path('output/parallel_t3/wilson/native_quartic_projectors_conditional_v10.json'))
    a=p.parse_args(argv);root=a.project_root.resolve()
    path=lambda x:(x if x.is_absolute() else root/x).resolve()
    output=path(a.output)
    try:
        if output.exists():raise V10Blocked('Refusing to overwrite '+str(output))
        report=audit_from_paths(root,path(a.native),path(a.v9))
        output.parent.mkdir(parents=True,exist_ok=True)
        with output.open('x',encoding='utf-8') as f:
            json.dump(report,f,indent=2,allow_nan=False);f.write('\n')
        print('v10: {registered_quartic_slots}/178 native registrations; '
              '{native_real_projector_diagonals_verified}/150 real projector diagonals; '
              '{native_portal_response_matrix_elements}/110 portal matrix elements'.format(**report))
        print(report['status'])
        print('Native symbolic Lam/projector VERIFIED; no numerical component Lam; NO Wilson evolution.')
        print(output)
        return 0
    except (V10Blocked,ValueError,KeyError,IndexError,TypeError,OSError,
            ArithmeticError,ImportError) as exc:
        print('Workstream A v10 BLOCKED: '+str(exc),file=sys.stderr)
        return 2
if __name__=='__main__':raise SystemExit(main())
