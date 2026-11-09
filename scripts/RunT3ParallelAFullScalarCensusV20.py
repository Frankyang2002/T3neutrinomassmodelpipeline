"""Full EFT1 quartic source-conditioned census: 150 real + 28 complex = 178 slots."""
from __future__ import annotations
import argparse
from pathlib import Path
import json
import sys
from Numerical.parallel_t3.wilson.full_scalar_census_v20 import audit,require
from Numerical.parallel_t3.wilson.complex_scalar_bridge_v19 import load,V19Blocked
BASE=Path('output/parallel_t3/wilson')

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root',type=Path,default=Path('.'))
    p.add_argument('--output',type=Path,default=BASE/'full_scalar_census_conditional_v20.json')
    p.add_argument('--resume',action='store_true')
    a=p.parse_args(argv);root=a.project_root.resolve()
    def f(name):return root/BASE/name
    try:
        report=audit(root,f('full_scalar_census_conditional_v17.json'),
            f('complex_scalar_bridge_conditional_v19.json'),
            f('rgbeta_quartic_projectors_native_v14.json'),
            f('matchete_topology_native_v18.json'),f('matchete_bc_extra_native_v18.json'),
            f('rgbeta_complex_conjugate_native_v15.json'),
            f('matchete_native_tensors.json'),f('rgbeta_abstract_components_conditional_v2.json'),
            f('additional_quartic_dictionary_conditional_v3.json'))
        out=a.output if a.output.is_absolute() else root/a.output
        if out.exists():
            require(a.resume,'refusing to overwrite '+str(out))
            old=load(out)
            require(old.get('format')==report['format'] and old.get('status')==report['status'] and
                    old.get('input_sha256')==report['input_sha256'] and
                    old.get('certification')==report['certification'] and
                    old.get('models')==report['models'],
                    'existing v20 census differs from current independent recheck')
            print('v20: existing report revalidated and preserved:',out)
        else:
            require(not (out.parent/(out.name+'.partial')).exists(),'stale v20 .partial exists')
            out.parent.mkdir(parents=True,exist_ok=True)
            tmp=out.parent/(out.name+'.partial')
            tmp.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
            require(load(tmp)['format']==report['format'],'v20 output JSON roundtrip failed')
            tmp.replace(out)
        print('v20: all16 178/178 native RGBeta scalar quartic slots; 150 real + 28 complex')
        print(report['status'])
        print('B/C mixed scalar basis, physical coupling map, GammaCC/Weyl and ADM remain blocked')
        print(out)
        return 0
    except (V19Blocked,ValueError,KeyError,TypeError,IndexError,OSError) as exc:
        print('v20 BLOCKED:',exc,file=sys.stderr);return 7
if __name__=='__main__':raise SystemExit(main())
