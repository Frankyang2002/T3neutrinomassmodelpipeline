"""Read-only v25-v27 integration: native mass basis + old one-step probe + EFT ordering.

Run: python -m scripts.RunT3ParallelABCMatcheteThresholdBridgeV25V27 --project-root . --resume
Exit 7 on failed source-integrity check. No external processes or physics matching.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys
from Numerical.parallel_t3.wilson import (
    bc_native_mass_eigen_reconcile_v25 as v25,
    bc_one_scale_probe_audit_v26 as v26,
    bc_threshold_order_ledger_v27 as v27,
)
BASE=Path('output/parallel_t3/wilson')
TARGETS=(('bc_native_matchete_mass_eigen_concordance_v25.json',v25),
         ('bc_existing_one_scale_matchete_probe_v26.json',v26),
         ('bc_fermion_first_scalar_threshold_ledger_v27.json',v27))

def safe_read(path:Path)->dict:
    return json.loads(path.read_text(encoding='utf-8-sig'),
        parse_constant=lambda s: (_ for _ in ()).throw(ValueError('nonfinite JSON constant')))

def write_or_verify(path:Path,report:dict, *,resume:bool,verify_only:bool)->str:
    if path.exists():
        if not resume and not verify_only:raise v25.GateError('output exists: '+str(path))
        if safe_read(path)!=report:raise v25.GateError('existing output differs; refusing overwrite '+str(path))
        return 'VERIFIED existing'
    if verify_only:raise v25.GateError('verify-only output missing '+str(path))
    path.parent.mkdir(parents=True,exist_ok=True)
    # Refuse overwriting an existing file. New files only, no rewriting witnesses.
    with path.open('x',encoding='utf-8') as f:
        json.dump(report,f,indent=2,allow_nan=False,sort_keys=True)
        f.write('\n')
    return 'CREATED'

def run(root:Path, *,resume:bool=True,verify_only:bool=False)->list[dict]:
    root=root.resolve();base=root/BASE
    path22=base/'bc_quartic_mass_readiness_conditional_v22.json'
    path21=base/'bc_native_mass_eigen_complex_conditional_v21.json'
    reports=[]
    for name,mod in TARGETS:
        if mod is v25: report=mod.audit(root,path21,path22)
        elif mod is v26: report=mod.audit(root,base/TARGETS[0][0])
        else: report=mod.audit(root,path22,base/TARGETS[0][0],base/TARGETS[1][0])
        path=base/name
        act=write_or_verify(path,report,resume=resume,verify_only=verify_only)
        print(f'{name}: {act}: {report["status"]}',flush=True)
        reports.append(report)
    for tag in v25.BC:
        d=reports[-1]['models'][tag]
        m=d['scalar_mass_eigenstates']
        print(f'{tag}: mX={m["X_gev"]:.8f} GeV, mY={m["Y_gev"]:.8f} GeV; '
              f'ln(mY/mX)={d["log_my_over_mx"]:.8f}; '
              f'EFT1 endpoint={d["saved_candidate_EFT1_scale_gev"]:.2f} GeV (BELOW both)',flush=True)
    print('v25-v27: native eigenfield checks linked; one-step probe observational only; '
          'fermion-first EFT1 matching, physical ADM and SMEFT BLOCKED.',flush=True)
    return reports

def main(argv=None)->int:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--project-root',type=Path,default=Path('.'))
    ap.add_argument('--resume',action='store_true')
    ap.add_argument('--verify-only',action='store_true')
    args=ap.parse_args(argv)
    try:
        run(args.project_root,resume=args.resume or args.verify_only,verify_only=args.verify_only)
        return 0
    except (v25.GateError,ValueError,KeyError,TypeError,FileNotFoundError,RuntimeError,OverflowError,AssertionError) as e:
        print('v25-v27 BLOCKED:',e,file=sys.stderr)
        return 7
if __name__=='__main__': raise SystemExit(main())
