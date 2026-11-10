#!/usr/bin/env python3
"""Fit 18 real UV Yukawa seeds in the A(-4) fixed-background rough EFT model.

Example:
  python -X utf8 -m scripts.RunT3A4Rough18RealFit
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import traceback
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0,str(ROOT))

from Numerical.parallel_t3.background.general_batch1_api import read_json
from Numerical.parallel_t3.complete.a4_rough_18fit import (
    FitWorkspace,RoughFitError,fit_18_real,
)

B1=ROOT/'output/parallel_t3/full_general_batch1/thesis_full_general_batch1_v1/A_-4'
CONFIG=ROOT/'configs/thesis_numerics/t3_thesis_A4_Dm2_pilot.json'
OUT=ROOT/'output/parallel_t3/a4_rough_18real_fit_v1/A_-4'


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--batch1-dir',type=Path,default=B1)
    p.add_argument('--config',type=Path,default=CONFIG)
    p.add_argument('--out',type=Path,default=OUT)
    p.add_argument('--max-nfev',type=int,default=180)
    p.add_argument('--regularization',type=float,default=.02)
    p.add_argument('--operator-factor',type=float,default=1.)
    args=p.parse_args(argv)
    if args.out.exists():
        p.error('Output folder exists; use --out with a new directory to retain previous fit')
    try:
        f=read_json(args.batch1_dir/'fermion_matching_input.json')
        s=read_json(args.batch1_dir/'scalar_matching_input.json')
        cfg=read_json(args.config)
        target_file=ROOT/cfg['fitting']['target']
        target=read_json(target_file)
        ws=FitWorkspace.from_sources(f,s,cfg,target,
                   gamma_factor=args.operator_factor,mu_low=100.)
        report,arrays=fit_18_real(ws,regularization=args.regularization,
                                   max_nfev=args.max_nfev)
        report['input_files']={
            'fermion_matching_input':str(args.batch1_dir/'fermion_matching_input.json'),
            'scalar_matching_input':str(args.batch1_dir/'scalar_matching_input.json'),
            'configuration':str(args.config),'target':str(target_file),
        }
        args.out.mkdir(parents=True,exist_ok=False)
        (args.out/'fit.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',
                                         encoding='utf-8')
        np.savez_compressed(args.out/'fitted_matrices.npz',**arrays)
        ob=report['final_fit_observables']
        print('A(-4) rough 18-real-Yukawa numerical fit: COMPUTED')
        print('Optimizer:',report['optimizer_message'],'nfev=',report['optimizer_nfev'])
        print('Uniform initial Yukawa scale:',report['initial_uniform_scale_each_Yukawa'])
        print('Start Gaussian chi2:',report['uniform_scale_start_chi2'])
        print('Final Gaussian chi2:',report['final_fit_chi2'])
        print('Masses at 100 GeV [eV]:',ob['masses_eV'])
        print('Delta m21^2 and m31^2 [eV^2]:',ob['dm21_sq_eV2'],ob['dm31_sq_eV2'])
        print('sin2(theta12,theta13,theta23):',
              ob['sin2_theta12'],ob['sin2_theta13'],ob['sin2_theta23'])
        print('Approximate fixed-background UV->F transfer used; full physical matching NOT computed')
        print('Output:',args.out/'fit.json')
        return 0
    except (OSError,ValueError,RuntimeError,KeyError,TypeError,
            np.linalg.LinAlgError) as exc:
        print('BLOCKED:',exc,file=sys.stderr)
        return 2

if __name__=='__main__':
    raise SystemExit(main())
