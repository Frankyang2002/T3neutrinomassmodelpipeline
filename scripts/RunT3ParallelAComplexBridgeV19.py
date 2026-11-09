"""Source-guarded v19: 16 native topology CGs + B/C 14 extra CGs + 28 native RGBeta channels."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys
from Numerical.parallel_t3.wilson.complex_scalar_bridge_v19 import audit,load,V19Blocked,require

BASE=Path('output/parallel_t3/wilson')

def path(root:Path,p:Path)->Path:
    return p if p.is_absolute() else root/p

def resume(out:Path, fresh:dict)->None:
    old=load(out)
    require(old.get('format')==fresh['format'] and
            old.get('status')==fresh['status'] and
            old.get('input_sha256')==fresh['input_sha256'] and
            old.get('certification')==fresh['certification'] and
            old.get('baseline_EFT1_complex_slots')==28 and
            old.get('native_Matchete_BC_extra_CGs_rechecked')==14 and
            float(old.get('max_native_topology_6_over_24_real_action_relative_residual',1.))<1e-9,
            'existing v19 result stale, altered or unsafe')
    for tag in fresh['models']:
        require(old.get('models',{}).get(tag,{}).get('baseline_complex_quartic_slots')==
                fresh['models'][tag]['baseline_complex_quartic_slots'],
                'saved v19 model census changed: '+tag)
    for tag in ('B_-1','C_-1'):
        a=old.get('BC_baseline_six_of_seven_conditional_source_maps',{}).get(tag, {})
        b=fresh['BC_baseline_six_of_seven_conditional_source_maps'][tag]
        require(a==b,'saved v19 B/C 6/7 coupling map changed: '+tag)

def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--project-root',type=Path,default=Path('.'))
    ap.add_argument('--topology',type=Path,default=BASE/'matchete_topology_native_v18.json')
    ap.add_argument('--bc',type=Path,default=BASE/'matchete_bc_extra_native_v18.json')
    ap.add_argument('--native-rgbeta',type=Path,default=BASE/'rgbeta_quartic_projectors_native_v14.json')
    ap.add_argument('--v15',type=Path,default=BASE/'rgbeta_complex_conjugate_native_v15.json')
    ap.add_argument('--mat50',type=Path,default=BASE/'matchete_native_tensors.json')
    ap.add_argument('--rgb50',type=Path,default=BASE/'rgbeta_abstract_components_conditional_v2.json')
    ap.add_argument('--v3',type=Path,default=BASE/'additional_quartic_dictionary_conditional_v3.json')
    ap.add_argument('--output',type=Path,default=BASE/'complex_scalar_bridge_conditional_v19.json')
    ap.add_argument('--resume',action='store_true')
    a=ap.parse_args(argv);root=a.project_root.resolve()
    try:
        files=[path(root,x) for x in (a.topology,a.bc,a.native_rgbeta,a.v15,a.mat50,a.rgb50,a.v3)]
        report=audit(root,*files)
        out=path(root,a.output)
        if out.exists():
            require(a.resume,'refusing to overwrite existing '+str(out))
            resume(out,report)
            print('v19: existing report independently revalidated and preserved:',out)
        else:
            require(not (out.parent/(out.name+'.partial')).exists(),'stale .partial output exists')
            out.parent.mkdir(parents=True,exist_ok=True)
            part=out.parent/(out.name+'.partial')
            part.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
            require(load(part)['format']==report['format'],'v19 JSON roundtrip failed')
            part.replace(out)
        print('v19: 16/16 native Matchete topology CGs; B/C 14/14 extra CGs verified')
        print('v19: 12 baseline B/C special subspace + 16 topology = 28 native RGBeta complex channels')
        print('v19:',report['status'],'residual',report['max_native_topology_6_over_24_real_action_relative_residual'])
        print('NO physical complex coupling map or Wilson evolution authorized')
        return 0
    except (V19Blocked,ValueError,KeyError,TypeError,IndexError,OSError) as exc:
        print('v19 BLOCKED:',exc,file=sys.stderr);return 7
if __name__=='__main__':raise SystemExit(main())
