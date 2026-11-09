"""Produce/recheck saved v18 Matchete CGs in existing B/C complex mass eigenbases."""
from __future__ import annotations
import argparse,json,sys
from pathlib import Path
from Numerical.parallel_t3.wilson.bc_mass_eigen_native_bridge_v21 import (
    BCBridgeBlocked, audit, require, read_json,
)

BASE=Path('output/parallel_t3/wilson')

def resolved(root:Path,name:Path)->Path:
    return name if name.is_absolute() else root/name

def write_or_verify(path:Path,doc:dict,*,resume:bool)->str:
    if path.exists():
        require(resume,'refusing to overwrite '+str(path))
        previous=read_json(path)
        require(previous==doc,'existing source-conditioned v21 differs; refusing overwrite')
        return 'revalidated and retained'
    part=path.with_name(path.name+'.partial')
    require(not part.exists(),'stale .partial output '+str(part))
    path.parent.mkdir(parents=True,exist_ok=True)
    part.write_text(json.dumps(doc,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    require(read_json(part)==doc,'v21 JSON roundtrip mismatch')
    part.replace(path)
    return 'generated'

def main(argv=None)->int:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--project-root',type=Path,default=Path('.'))
    ap.add_argument('--topology',type=Path,default=BASE/'matchete_topology_native_v18.json')
    ap.add_argument('--bc',type=Path,default=BASE/'matchete_bc_extra_native_v18.json')
    ap.add_argument('--v19',type=Path,default=BASE/'complex_scalar_bridge_conditional_v19.json')
    ap.add_argument('--v20',type=Path,default=BASE/'full_scalar_census_conditional_v20.json')
    ap.add_argument('--B-mass-basis',type=Path,default=Path('output/general_t3_mass_basis/B_m1_final_mass_basis.json'))
    ap.add_argument('--C-mass-basis',type=Path,default=Path('output/general_t3_mass_basis/C_m1_final_mass_basis.json'))
    ap.add_argument('--output',type=Path,default=BASE/'bc_native_mass_eigen_complex_conditional_v21.json')
    ap.add_argument('--trials',type=int,default=8)
    ap.add_argument('--resume',action='store_true')
    a=ap.parse_args(argv);root=a.project_root.resolve()
    try:
        report=audit(root,topology_path=resolved(root,a.topology),bc_path=resolved(root,a.bc),
                     v19_path=resolved(root,a.v19),v20_path=resolved(root,a.v20),
                     mass_paths={'B_-1':resolved(root,a.B_mass_basis),
                                 'C_-1':resolved(root,a.C_mass_basis)},trials=a.trials)
        out=resolved(root,a.output)
        result=write_or_verify(out,report,resume=a.resume)
        print(f'v21: {result} {out}')
        print('v21: v19 SHA256 inputs {}/7 verified; legacy Windows separator keys {}'.format(
            report['v19_input_provenance']['source_files_sha256_rechecked'],
            report['v19_input_provenance']['legacy_windows_separator_keys_accepted']))
        print('v21: 16/16 B/C native Matchete complex quartic CGs rotated; '+'residual '+str(report['max_complex_quartic_action_transport_relative_residual']))
        print('v21: 12 baseline B/C complex source slots, optional seventh preserved; no physical matching')
        print(report['status'])
        return 0
    except (BCBridgeBlocked,ValueError,KeyError,IndexError,TypeError,OSError) as exc:
        print('v21 BLOCKED:',exc,file=sys.stderr);return 7
if __name__=='__main__':raise SystemExit(main())
