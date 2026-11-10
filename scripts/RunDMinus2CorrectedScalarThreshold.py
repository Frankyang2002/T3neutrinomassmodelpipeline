#!/usr/bin/env python3
"""Read-only D(-2) corrected EFT1 -> simultaneous 8-real-scalar threshold.

Reconstructs corrected full Wilson and quartic trajectories from saved benchmark,
runs the existing candidate ADM to MS, then calls the repository's scalar
matching function in rge_improved mode. Matches BOTH S1 (six real) and S2
(two real) at one MS. Prints source/basis/Weinberg and mass consistency checks.
No Yukawa optimization, no production edits, no SMEFT running or mass spectrum.

The existing threshold mass data and identity R_wilson_from_background come
from a previously PROVISIONAL source convention; this script checks SU2/U1
mass covariance but cannot certify phase/normalization of the native singlet,
the signed GammaCC->Weyl map or complete one-loop threshold amplitudes.
"""
from __future__ import annotations
import argparse
import sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
sys.path.insert(0,str(ROOT/'scripts'))

from RunDMinus2CorrectedEFT1 import prepare
from Numerical.parallel_t3.wilson.workflow import run_wilson
from Numerical.parallel_t3.threshold.matching import (
    match_scalar_threshold,scalar_gauge_ward_residuals
)

DEFAULT=Path("output/parallel_t3/rough_inputs/fit16_cc2dcbec4a")


def matrix(label,x):
    print(label+' [GeV^-1]:')
    print(np.array2string(np.asarray(x),precision=9,suppress_small=False))


def main():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument('--input-root',type=Path,default=DEFAULT)
    a.add_argument('--solver',choices=('RK2','DOP853'),default='DOP853')
    a.add_argument('--steps',type=int,default=4)
    a.add_argument('--rtol',type=float,default=1e-6)
    a.add_argument('--atol',type=float,default=1e-21)
    a.add_argument('--variant',
                   choices=('fonseca_eq22_six_term','ward_covariant_candidate'),
                   default='fonseca_eq22_six_term')
    a.add_argument('--ward-tolerance',type=float,default=1e-8)
    args=a.parse_args()
    if args.steps<1:a.error('--steps must be at least 1')
    if args.ward_tolerance<=0:a.error('--ward-tolerance must be positive')

    folder=args.input_root/'a'/'D_m2'
    cb=args.input_root/'c'/'D_m2'/'background.npz'
    basis,bg,initial,mu,qold,qnew=prepare(folder)
    mf=float(mu[-1]);ms=float(mu[0])
    if not cb.is_file():
        raise FileNotFoundError(f'Missing retained scalar-threshold background: {cb}')
    with np.load(cb,allow_pickle=False) as f:
        required=('mass_squared_gev2_ms','R_wilson_from_background','quartic_real_full_ms')
        absent=[key for key in required if key not in f.files]
        if absent:
            raise KeyError(f'{cb}: missing keys {absent}, available {f.files}')
        mass=np.asarray(f['mass_squared_gev2_ms'],float)
        R=np.asarray(f['R_wilson_from_background'],float)
        old_ms_q=np.asarray(f['quartic_real_full_ms'],float)
    n=basis.ns-4
    if mass.shape!=(n,n) or R.shape!=(n,n):
        raise ValueError(f'Mass/R shape mismatch: {mass.shape} {R.shape}')
    if old_ms_q.shape!=(basis.ns,)*4:
        raise ValueError(f'Unexpected stored scalar quartic shape {old_ms_q.shape}')
    if np.linalg.norm(R.T@R-np.eye(n))>1e-10:
        raise ValueError('Nonorthogonal recorded scalar R')
    if np.linalg.norm(R-np.eye(n))>1e-10:
        raise RuntimeError('Existing c background uses a non-identity R; '
                           'its mapping requires independent source alignment')
    # The c saved quartic must match the a saved MS value (prior to correction).
    drift=np.linalg.norm(old_ms_q-qold[0])
    rel=drift/max(np.linalg.norm(qold[0]),1e-300)
    if rel>1e-9:
        raise RuntimeError(f'A/C scalar background mismatch at MS: relative {rel:g}')
    if abs(ms-1000)>1e-5:
        raise ValueError(f'D(-2) expected MS=1000 GeV, got {ms}')
    eig=np.linalg.eigvalsh(mass)
    if np.min(eig)<=0:raise ValueError('Nonpositive scalar mass squared')

    gens=np.asarray([s for group in ('SU2','U1')
                     for _,s in basis.generators()[group]])
    ward=scalar_gauge_ward_residuals(qnew[0],mass,gens)
    print('=== D(-2): corrected EFT1 -> both-scalar threshold ===',flush=True)
    print('Input root:',args.input_root,flush=True)
    print(f'MF={mf:g} GeV; MS={ms:g} GeV; heavy real scalars={n}',flush=True)
    print('A/C scalar quartic relative consistency:',f'{rel:.9e}',flush=True)
    print('Retained mass eigenvalues [GeV]:',np.sqrt(eig),flush=True)
    print('Scalar SU2/U1 Ward diagnostic:',ward,flush=True)
    if not ward['passed'] or max(ward['quartic_invariance_relative'],
                                 ward['mass_commutator_relative'])>=args.ward_tolerance:
        raise RuntimeError('Scalar quartic/mass background fails gauge Ward')
    print(f'Running {args.solver} EFT1 with {args.steps} intervals...',flush=True)
    result=run_wilson(basis,initial,bg,mf,ms,
                      adm_variant=args.variant,solver=args.solver,steps=args.steps,
                      rtol=args.rtol,atol=args.atol)
    print('EFT1 finished: RHS calls',result['metadata']['rhs_evaluations'],flush=True)
    print('||C5_EFT1(MS)||:',f'{np.linalg.norm(result["C5_final"]):.12e}',flush=True)

    ll=np.array(basis.leptons,dtype=int)
    if ll.shape!=(3,2):raise RuntimeError('Unexpected lepton index geometry')
    h=np.array(basis.higgs,dtype=int)
    heavy=np.arange(4,basis.ns,dtype=int)
    print('Matching all eight heavy real scalars simultaneously at MS...',flush=True)
    out=match_scalar_threshold(
        full_wilson_ms=result['C_final'],
        lepton_indices=ll,
        higgs_indices=h,
        heavy_indices=heavy,
        full_quartic_ms=qnew[0],
        scalar_mass2_background=mass,
        R_wilson_from_background=R,
        mu_s=ms,
        mode='rge_improved',
        separately_supplied_c5_ms=result['C5_final'],
        projection_tolerance=args.ward_tolerance,
        enforce_pure_weinberg=True
    )
    print('Number of integrated-out real scalars:',out['scalar_components_removed'])
    print('Number of scalar eigenstate pairs:',out['bubble']['pair_count'])
    print('Weinberg projection residual EFT1:',f'{out["rge_input_C5_residual"]:.12e}')
    print('Weinberg projection residual scalar bubble:',
          f'{out["scalar_bubble_projection_residual"]:.12e}')
    for name,key in (('EFT1 C5','c5_eft1'),
                     ('Scalar bubble delta C5','c5_scalar_bubble'),
                     ('SMEFT C5 at MS','c5_smeft_at_ms')):
        print(f'||{name}|| [GeV^-1]: {np.linalg.norm(out[key]):.12e}')
        matrix(name,out[key])
    print('Bubble pair count expected/actual:',n*(n+1)//2,out['bubble']['pair_count'])
    print('Counting:',out['counting_mode'],out['bubble_Wilson_order_source'])
    print('No scalar kept below MS in this one-scale matcher.')
    print('No saved/production files modified.')
    print('Numerical candidate only: signed GammaCC/Weyl factor, source scalar '
          'phase/mass mapping and full one-loop matching NOT certified.')

if __name__=='__main__':
    main()
