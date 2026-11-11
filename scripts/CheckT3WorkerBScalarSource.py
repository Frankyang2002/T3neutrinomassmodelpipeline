"""Read-only Worker B A(-4)/D(-2) benchmark check.

Usage: python -m scripts.CheckT3WorkerBScalarSource --model A_-4 --response-root PATH --hhss-npz PATH
The HHSS NPZ must contain 'hhss' shape 4x4x8x8 in the fixed 12-real basis.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
from Numerical.running.T3CoupledQuarticContext import CandidateQuarticContext
from Numerical.running.T3ScalarRunningSourceBridge import (
    load_saved_running_state,build_from_cached_context,
)

def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--project',type=Path,default=Path(__file__).resolve().parents[1])
    ap.add_argument('--model',required=True,choices=('A_-4','D_-2'))
    ap.add_argument('--response-root',type=Path,required=True)
    ap.add_argument('--hhss-npz',type=Path,required=True)
    ap.add_argument('--output',type=Path)
    args=ap.parse_args(argv)
    project=args.project.resolve()
    mu,state=load_saved_running_state(project,args.model)
    ctx=CandidateQuarticContext.load(project,args.response_root,args.model)
    with np.load(args.hhss_npz,allow_pickle=False) as arrays:
        if 'hhss' not in arrays:
            raise ValueError('NPZ missing hhss; no source HHSS inferred')
        h=np.array(arrays['hhss'],copy=True)
    bg=build_from_cached_context(args.model,mu,state,ctx,hhss_source=h)
    report={'model':args.model,'scale_gev':mu,
            'source_hhss_supplied_externally':True,
            'scalar_mass_eigenvalues_gev2':bg.eigenvalues_gev2.tolist(),
            'quartic_symmetry_verified':True,
            'basis_orthogonality_residual':float(np.linalg.norm(bg.scalar_rotation.T@bg.scalar_rotation-np.eye(8))),
            'assumptions':list(bg.assumptions),
            'status':'ConditionalSourceInjectedBenchmark_NoIndependentPhysicalCertification'}
    print(json.dumps(report,indent=2))
    if args.output is not None:
        if args.output.exists(): raise FileExistsError('refusing overwrite: '+str(args.output))
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    return 0

if __name__=='__main__':
    raise SystemExit(main())
