"""Source-guarded all-16 Matchete nonportal self/cross real-coupling dictionary."""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path
from Numerical.parallel_t3.wilson.nonportal_native_dictionary_v16 import (
    audit_from_paths, V16Blocked, read, relative,
)


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root',type=Path,default=Path('.'))
    p.add_argument('--matchete',type=Path,default=Path('output/parallel_t3/wilson/matchete_nonportal_native_v16.json'))
    p.add_argument('--portal',type=Path,default=Path('output/parallel_t3/wilson/matchete_portals_native_v5.json'))
    p.add_argument('--native-v14',type=Path,default=Path('output/parallel_t3/wilson/rgbeta_quartic_projectors_native_v14.json'))
    p.add_argument('--output',type=Path,default=Path('output/parallel_t3/wilson/nonportal_dictionary_conditional_v16.json'))
    p.add_argument('--samples',type=int,default=5)
    p.add_argument('--resume',action='store_true',help='Revalidate an existing v16 result without rewriting it')
    a=p.parse_args(argv)
    root=a.project_root.resolve()
    def pp(path):return path if path.is_absolute() else root/path
    out=pp(a.output)
    if Path(str(out)+'.partial').exists():p.error('Refusing existing partial: '+str(out))
    if out.exists() and not a.resume:p.error('Refusing existing output: '+str(out))
    try:
        result=audit_from_paths(pp(a.matchete),pp(a.portal),pp(a.native_v14),root,samples=a.samples)
        if out.exists():
            old=read(out)
            if (old.get('format')!=result['format'] or old.get('status')!=result['status'] or
                old.get('input_sha256')!=result['input_sha256'] or
                old.get('nonportal_coupling_slots')!=76 or
                old.get('nonportal_sectors')!=48 or
                old.get('certification')!=result['certification']):
                raise V16Blocked('saved v16 output does not match freshly validated input')
            for tag in result['models']:
                for fam in ('self1','self2','cross'):
                    x=old['models'][tag][fam];y=result['models'][tag][fam]
                    if (x.get('matchete_couplings')!=y['matchete_couplings'] or
                        x.get('rgbeta_couplings')!=y['rgbeta_couplings'] or
                        relative(x.get('conditional_rgbeta_from_matchete_real_matrix'),
                                 y['conditional_rgbeta_from_matchete_real_matrix'])>1.e-10):
                        raise V16Blocked(f'saved v16 map altered: {tag}/{fam}')
            print('v16: existing output independently rechecked; no overwrite: '+str(out))
            return 0
        out.parent.mkdir(parents=True,exist_ok=True)
        tmp=Path(str(out)+'.partial')
        tmp.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
        if json.loads(tmp.read_text(encoding='utf-8'))['nonportal_coupling_slots'] != 76:
            raise V16Blocked('v16 output readback failed')
        tmp.replace(out)
    except (V16Blocked,ValueError,OSError,KeyError,TypeError) as exc:
        print('Blocked_NonportalV16: '+str(exc),file=sys.stderr)
        return 7
    print(f"v16: {result['nonportal_coupling_slots']}/76 native scalar nonportal CGs; "
          f"{result['nonportal_sectors']}/48 sectors; max residual {result['max_relative_residual']:.3e}")
    print(result['status'])
    print('No native RGBeta triplet component frame or mixed B/C mass-basis certification.')
    print(out)
    return 0

if __name__=='__main__':raise SystemExit(main())
