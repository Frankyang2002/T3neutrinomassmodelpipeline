#!/usr/bin/env python3
"""Read-only convergence study of corrected D(-2) full-tensor EFT1 evolution.

Uses existing saved MF boundary, native EFT1 trajectory and repository ADM.
No Yukawa optimisation, scalar matching, or production-file modifications.
DOP853 on the full 45x45x12x12 complex Wilson tensor is memory/CPU intensive.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))
from RunDMinus2CorrectedEFT1 import prepare
from DiagnoseDMinus2CanonicalConventions import DEFAULT, ward, llhh
from Numerical.parallel_t3.wilson.workflow import run_wilson
from Numerical.parallel_t3.threshold.matching import project_master


def run_case(label, basis, bg, initial, mf, ms, *, solver, steps, rtol, atol, variant):
    print(f'Running {label}: {solver}, {steps} intervals', flush=True)
    result = run_wilson(
        basis, initial, bg, mf, ms, adm_variant=variant,
        solver=solver, steps=steps, rtol=rtol, atol=atol,
    )
    c5 = np.asarray(result['C5_final'], complex)
    c = np.asarray(result['C_final'], complex)
    projected, residual = project_master(llhh(c, basis), enforce_pure=False)
    if np.linalg.norm(c5-projected) > max(1e-24, 1e-8*np.linalg.norm(c5)):
        raise RuntimeError(f'{label}: independent Weinberg projectors disagree')
    ws = ward(c, basis)
    worst = max(ws['SU2']+ws['U1'])
    print(f'{label}: ||C5(MS)||={np.linalg.norm(c5):.12e} GeV^-1; '
          f'Weinberg residual={residual:.3e}; Ward max={worst:.3e}; '
          f'RHS calls={result["metadata"]["rhs_evaluations"]}', flush=True)
    return c5


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input', type=Path, default=DEFAULT)
    ap.add_argument('--rtol', type=float, default=1e-6)
    ap.add_argument('--atol', type=float, default=1e-21)
    ap.add_argument('--dop-steps', type=int, default=4)
    ap.add_argument('--rk-steps', type=int, nargs='+', default=[4, 8])
    ap.add_argument('--variant', choices=['fonseca_eq22_six_term',
                                          'ward_covariant_candidate'],
                    default='fonseca_eq22_six_term')
    ap.add_argument('--skip-dop', action='store_true',
                    help='Run RK2-only sanity check before the expensive DOP853 run')
    args = ap.parse_args()
    if args.rtol<=0 or args.atol<=0 or args.dop_steps<1 or any(s<1 for s in args.rk_steps):
        ap.error('Tolerances and steps must be positive')
    basis, bg, initial, mu, _, _ = prepare(args.input)
    mf, ms = float(mu[-1]), float(mu[0])
    print(f'D(-2) convergence: MF={mf:g} GeV, MS={ms:g} GeV', flush=True)
    print('Corrected mapping: triplet O only; J portal throughout saved background.',flush=True)
    print('Full 45 Weyl x 12 scalar tensor; DOP853 may use substantial RAM.',flush=True)
    current = ward(initial,basis)
    if max(current['SU2']+current['U1'])>1e-9:
        raise RuntimeError('Corrected initial boundary fails Ward test')
    outputs = {}
    for step in dict.fromkeys(args.rk_steps):
        key = f'RK2_{step}'
        outputs[key] = run_case(key,basis,bg,initial,mf,ms,solver='RK2',
                                steps=step,rtol=args.rtol,atol=args.atol,
                                variant=args.variant)
    if not args.skip_dop:
        outputs['DOP853'] = run_case('DOP853',basis,bg,initial,mf,ms,
                                      solver='DOP853',steps=args.dop_steps,
                                      rtol=args.rtol,atol=args.atol,variant=args.variant)
    reference_key='DOP853' if 'DOP853' in outputs else list(outputs)[-1]
    reference=outputs[reference_key]
    print(f'\nMatrix convergence relative to {reference_key}:', flush=True)
    for name, x in outputs.items():
        diff=np.linalg.norm(x-reference)
        rel=diff/max(np.linalg.norm(reference),1e-300)
        norm_ratio=np.linalg.norm(x)/max(np.linalg.norm(reference),1e-300)
        print(f' {name:12s}: relative matrix difference={rel:.12e}; '
              f'norm ratio={norm_ratio:.12g}', flush=True)
    print('\nReference C5(MS) [GeV^-1]:')
    print(np.array2string(reference,precision=9))
    print('No files modified; no threshold or SMEFT calculation performed.')
    print('Candidate ADM and signed GammaCC/Weyl normalization remain uncertified.')

if __name__ == '__main__':
    main()
