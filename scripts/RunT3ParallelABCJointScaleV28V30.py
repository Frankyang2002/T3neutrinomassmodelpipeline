#!/usr/bin/env python3
"""Run B/C(-1) all-scalar joint-scale candidate EFT1 replay and guarded handoff.

No scalar threshold Match[], no fermion Wilson boundary, no SMEFT calculation.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0,str(ROOT))

from Numerical.parallel_t3.wilson import bc_joint_scalar_scale_v28 as v28
from Numerical.parallel_t3.wilson import bc_joint_scalar_scale_audit_v29 as v29
from Numerical.parallel_t3.wilson import bc_fermion_first_common_handoff_v30 as v30

OUTPUT=ROOT/'output/parallel_t3/wilson'


def emit(path:Path,payload:dict) -> str:
    """Atomic create or bytewise revalidation; never overwrite saved evidence."""
    encoded=json.dumps(payload,indent=2,sort_keys=True,allow_nan=False)+'\n'
    path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():
        current=v28.read(path)
        if current != payload:
            raise v28.JointScaleBlocked(f'Existing evidence differs, refusing overwrite: {path}')
        return 'REUSED'
    tmp=path.with_name(path.name+'.tmp')
    if tmp.exists():
        raise v28.JointScaleBlocked('Stale temporary file requires manual inspection: '+str(tmp))
    try:
        tmp.write_text(encoded,encoding='utf-8')
        tmp.replace(path)
    finally:
        if tmp.exists():tmp.unlink()
    return 'CREATED'


def run(*,root:Path=ROOT,outdir:Path=OUTPUT,backend:str='auto',scheme:str='geometric',
        v27path:Path|None=None,qeft_dir:Path|None=None,mass_dir:Path|None=None,
        stage33_dir:Path|None=None,mu2_fraction:float=.01,mu2_phase:float=.6,
        higgs_mass_squared:float=-7800.) -> dict:
    root=Path(root).resolve();outdir=Path(outdir)
    v27path=Path(v27path or root/'output/parallel_t3/wilson/bc_fermion_first_scalar_threshold_ledger_v27.json')
    if not v27path.is_file():
        raise v28.JointScaleBlocked('Missing v27 ledger: '+str(v27path)+'. First run v25-v27.')
    v27data=v28.read(v27path)
    from Numerical.parallel_t3.wilson.bc_threshold_order_ledger_v27 import FORMAT as F27,STATUS as S27
    v28.need(v27data.get('format')==F27 and v27data.get('status')==S27,
             'v27 source ledger format/status invalid')
    rows={}
    paths={}
    for model in v28.MODELS:
        name='bc_joint_scalar_candidate_state_v28_'+model.replace('-','m')+'.json'
        path=outdir/name
        if path.exists():
            # Save computational time. The independent v29 audit rehashes
            # *every* live file before this evidence can be reused.
            d=v28.read(path)
            v29.audit_model(root,d)
            v28.need(d['model']==model and d.get('common_scalar_matching_scale',{}).get('scheme')==scheme,
                     'saved v28 model/scale prescription differs')
            if backend!='auto':
                v28.need(d.get('backend')==backend,'saved v28 backend differs from selection')
            action='REVALIDATED'
        else:
            d=v28.run_model(root,model,backend=backend,scheme=scheme,
                qeft_dir=qeft_dir,mass_dir=mass_dir,stage33_dir=stage33_dir,
                higgs_mass_squared=higgs_mass_squared,mu2_fraction=mu2_fraction,
                mu2_phase=mu2_phase)
            action=emit(path,d)
        rows[model]=d;paths[model]=path
        m=d['common_scalar_matching_scale']
        print(f"v28 {model}: {action}; backend={d['backend']}; joint muS={m['matching_scale_gev']:.9g} GeV; "
              f"running X/Y={m['mX_gev']:.9g}/{m['mY_gev']:.9g} GeV; "
              f"max |ln(muS/m)|={m['max_absolute_scalar_matching_log']:.5g}",flush=True)
    checked=v29.audit(root,rows['B_-1'],rows['C_-1'])
    checked['source_sha256']={tag:v28.sha(path) for tag,path in paths.items()}
    checked['source_paths']={tag:path.relative_to(root).as_posix() for tag,path in paths.items()}
    p29=outdir/'bc_joint_scalar_candidate_readiness_v29.json'
    action29=emit(p29,checked)
    print('v29: '+action29+'; both common EFT1 scalar-scale states reverified',flush=True)
    final=v30.audit(v27path,p29,root=root)
    p30=outdir/'bc_fermion_first_joint_scalar_handoff_v30.json'
    action30=emit(p30,final)
    print('v30: '+action30+'; joint scalar candidate inputs prepared; Wilson and physical matching BLOCKED',flush=True)
    return final


def main()->int:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root',type=Path,default=ROOT)
    ap.add_argument('--output-dir',type=Path,default=None)
    ap.add_argument('--v27',type=Path,default=None)
    ap.add_argument('--backend',choices=['auto','five_mass','quadratic_only'],default='auto')
    ap.add_argument('--scheme',choices=['geometric','upper'],default='geometric')
    ap.add_argument('--qeft-dir',type=Path,default=None)
    ap.add_argument('--mass-dir',type=Path,default=None)
    ap.add_argument('--stage33-dir',type=Path,default=None)
    ap.add_argument('--higgs-mass-squared',type=float,default=-7800.)
    ap.add_argument('--mu2-fraction',type=float,default=.01)
    ap.add_argument('--mu2-phase',type=float,default=.6)
    args=ap.parse_args()
    root=args.root.resolve()
    try:
        run(root=root,outdir=args.output_dir or root/'output/parallel_t3/wilson',
            v27path=args.v27 or root/'output/parallel_t3/wilson/bc_fermion_first_scalar_threshold_ledger_v27.json',
            backend=args.backend,scheme=args.scheme,qeft_dir=args.qeft_dir,
            mass_dir=args.mass_dir,stage33_dir=args.stage33_dir,
            higgs_mass_squared=args.higgs_mass_squared,mu2_fraction=args.mu2_fraction,
            mu2_phase=args.mu2_phase)
    except (v28.JointScaleBlocked,ValueError,RuntimeError,OSError,KeyError,TypeError) as ex:
        print('v28-v30 BLOCKED: '+str(ex),file=sys.stderr)
        return 7
    return 0

if __name__=='__main__':raise SystemExit(main())
