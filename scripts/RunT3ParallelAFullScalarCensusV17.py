"""Revalidate v9, v14 and v16 evidence; assemble a conditional 150-slot census."""
from __future__ import annotations
import argparse,json,sys
from pathlib import Path
from Numerical.parallel_t3.wilson.full_scalar_census_v17 import audit_from_paths
from Numerical.parallel_t3.wilson.nonportal_native_dictionary_v16 import read

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root',type=Path,default=Path('.'))
    p.add_argument('--output',type=Path,
                   default=Path('output/parallel_t3/wilson/full_scalar_census_conditional_v17.json'))
    p.add_argument('--resume',action='store_true',help='Verify existing census without rewriting')
    a=p.parse_args(argv)
    root=a.project_root.resolve();base=root/'output/parallel_t3/wilson'
    out=a.output if a.output.is_absolute() else root/a.output
    if Path(str(out)+'.partial').exists():p.error('refusing existing partial: '+str(out))
    if out.exists() and not a.resume:p.error('refusing existing output: '+str(out))
    try:
        r=audit_from_paths(root,
             base/'matchete_nonportal_native_v16.json',
             base/'matchete_portals_native_v5.json',
             base/'rgbeta_quartic_projectors_native_v14.json',
             base/'nonportal_dictionary_conditional_v16.json',
             base/'portal_real_action_conditional_v9.json',
             base/'native_quartic_projectors_conditional_v14.json')
        if out.exists():
            old=read(out)
            if (old.get('format')!=r['format'] or old.get('status')!=r['status'] or
                old.get('input_sha256')!=r['input_sha256'] or
                old.get('certification')!=r['certification'] or
                old.get('models')!=r['models']):
                raise ValueError('existing v17 output differs from freshly checked source')
            print('v17: existing full scalar census independently VERIFIED, not overwritten: '+str(out))
            return 0
        out.parent.mkdir(parents=True,exist_ok=True)
        tmp=Path(str(out)+'.partial');tmp.write_text(json.dumps(r,indent=2,allow_nan=False)+'\n',encoding='utf-8')
        if json.loads(tmp.read_text(encoding='utf-8'))['native_RGBeta_real_projector_unit_diagonals']!=150:
            raise ValueError('readback census not 150')
        tmp.replace(out)
    except (OSError,ValueError,KeyError,TypeError) as exc:
        print('Blocked_FullScalarCensusV17: '+str(exc),file=sys.stderr);return 7
    print('v17 150/150 real scalar registry: Higgs 16 + portals 58 + nonportal 76')
    print('RGBeta native real Lam/projector diagonals 150/150, portal responses 110/110')
    print(r['status'])
    print('Complex 28 exploratory; physical Wilson evolution and mixed-mass matching blocked')
    print(out)
    return 0
if __name__=='__main__':raise SystemExit(main())
