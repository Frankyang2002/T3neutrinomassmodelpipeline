#!/usr/bin/env python3
"""Run corrected D(-2) full 45-Weyl EFT1 Wilson RGE, without refitting Yukawas.

EFT order: UV -> all F integrated out at MF (retained formal boundary)
-> EFT1 with BOTH scalars active to MS -> scalar threshold (NOT done here).

Reuses the installed repository Wilson ADM and saved native EFT1 background.
At every background scale, corrects ONLY the D(-2) lambdaT3 portal component
using the spherical spin-one metric J. The retained formal Wilson boundary
gets the triplet scalar map O only; NO lepton epsilon.

Numerical-candidate results, not physical certification:
GammaCC/Weyl signed factor, complete one-loop ADM, singlet-phase/source mapping
and scalar matching remain uncertain. Saved inputs and production code unchanged.
"""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
sys.path.insert(0,str(ROOT/'scripts'))

from DiagnoseDMinus2CanonicalConventions import (
    DEFAULT, triplet_rotation, corrected_portal, read_inputs, ward, llhh
)
from Numerical.parallel_t3.wilson.basis import build_basis
from Numerical.parallel_t3.wilson.coordinates import embed_formal_llss,embed_c5
from Numerical.parallel_t3.rough_inputs import builder
from Numerical.parallel_t3.wilson.interchange import Background
from Numerical.parallel_t3.wilson.workflow import run_wilson
from Numerical.parallel_t3.threshold.matching import project_master


def prepare(folder):
    d=read_inputs(folder)
    basis=build_basis('D_-2')
    mu=np.asarray(d['mu_gev'],float)
    native_mu=np.asarray(d['mu_eft1_gev'],float)
    if mu.ndim!=1 or native_mu.shape!=mu.shape:
        raise ValueError('Incorrect saved background/native scale shapes')
    if np.max(np.abs(mu-native_mu))/max(np.max(abs(mu)),1.)>1e-10:
        raise ValueError('Saved RGBeta and canonical background scales disagree')
    if len(mu)<2 or np.any(mu<=0) or np.unique(mu).size!=mu.size:
        raise ValueError('Missing/nonunique/invalid EFT1 scales')
    ix=np.argsort(mu)
    mu=mu[ix]
    y=np.asarray(d['Y'],complex)[ix]
    q_old=np.asarray(d['lambda4'],float)[ix]
    g=np.column_stack([np.asarray(d[k],float)[ix] for k in ('g2','gY','g3')])
    z=np.asarray(d['eft1__lambdaT3'])[ix]
    old_portal=builder.quartic_unit_tensors('D_-2')['lambdaT3']
    new_portal=corrected_portal(basis.ns)
    delta=new_portal-old_portal
    q_new=q_old+2*np.real(z[:,None,None,None,None]*delta[None,:,:,:,:])
    if not np.isfinite(q_new).all():raise ValueError('Nonfinite corrected quartic trajectory')
    if not all(np.allclose(q_new,q_new.swapaxes(1,k),atol=1e-11,rtol=1e-11)
               for k in (2,3,4)):
        raise ValueError('Corrected quartic permutation symmetry failure')
    raw=np.asarray(d['C_LLSS_formal_muF'],complex)
    if raw.shape!=(3,3,2,2,8,8):
        raise ValueError(f'Unexpected formal boundary shape {raw.shape}')
    c5hard=np.asarray(d['C5_hard'],complex)
    if c5hard.shape!=(3,3):raise ValueError('Incorrect C5 hard seed')
    c0=embed_formal_llss(raw,basis,triplet_rotation(),1.)
    c0+=embed_c5(c5hard,basis)
    background=Background(
        basis, np.log(mu),y,q_new,g,basis.generators(),
        {'input_path':str(folder),
         'source':'retained background; D(-2) J portal corrected at each saved scale'})
    return basis,background,c0,mu,q_old,q_new


def diagnostics(name,result,basis):
    c0=np.asarray(result['C_initial'])
    c1=np.asarray(result['C_final'])
    c5=np.asarray(result['C5_final'])
    ll=llhh(c1,basis)
    independent,residual=project_master(ll,enforce_pure=False)
    if np.linalg.norm(independent-c5)>max(1e-24,1e-8*np.linalg.norm(c5)):
        raise RuntimeError('Two independent C5 projectors disagree')
    print(f'\n{name}:')
    print('  ||C5(MF)|| [GeV^-1]:',f'{np.linalg.norm(result["C5_initial"]):.12e}')
    print('  ||C5(MS)|| [GeV^-1]:',f'{np.linalg.norm(c5):.12e}')
    print('  ||LLHH(MS)|| [GeV^-1]:',f'{np.linalg.norm(ll):.12e}')
    print('  LLHH Weinberg projection residual:',f'{residual:.12e}')
    print('  ||full C(MF)||:',f'{np.linalg.norm(c0):.12e}')
    print('  ||full C(MS)||:',f'{np.linalg.norm(c1):.12e}')
    print('  Wilson at MS SU2 max-entry:',
          ' '.join(f'{r:.6e}' for r in ward(c1,basis)['SU2']))
    print('  Wilson at MS U1 max-entry:',
          ' '.join(f'{r:.6e}' for r in ward(c1,basis)['U1']))
    print('  RHS evaluations:',result['metadata']['rhs_evaluations'])
    print('  LL first-order C5 norm:',
          f'{np.linalg.norm(result["leading_log_C5_linearized"]):.12e}')
    print('  C5(MS) [GeV^-1]:')
    print(np.array2string(c5,precision=7))
    return c5


def main():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument('--input',type=Path,default=DEFAULT)
    a.add_argument('--solver',choices=['RK2','DOP853'],default='RK2')
    a.add_argument('--steps',type=int,default=4)
    a.add_argument('--rtol',type=float,default=1e-6)
    a.add_argument('--atol',type=float,default=1e-21)
    a.add_argument('--variant',choices=['fonseca_eq22_six_term','ward_covariant_candidate'],
                   default='fonseca_eq22_six_term')
    a.add_argument('--compare-frozen',action='store_true',
                   help='Run second case with background frozen at MF for comparison')
    args=a.parse_args()
    if args.steps<1:a.error('--steps must be >=1')
    basis,bg,initial,mu,qold,qnew=prepare(args.input)
    mf=float(mu[-1]);ms=float(mu[0])
    print('=== D(-2) corrected EFT1: full Wilson evolution ===',flush=True)
    print('Input:',args.input,flush=True)
    print(f'MF={mf:.12g} GeV, MS={ms:.12g} GeV',flush=True)
    print(f'Solver={args.solver}, steps={args.steps}, ADM={args.variant}',flush=True)
    print('Canonical map: triplet O only; J-corrected quartic at ALL saved scales.',flush=True)
    print('Quartic MF change norm:',f'{np.linalg.norm(qnew[-1]-qold[-1]):.12e}',flush=True)
    print('Quartic MS change norm:',f'{np.linalg.norm(qnew[0]-qold[0]):.12e}',flush=True)
    w0=ward(initial,basis)
    print('Initial Ward SU2:', ' '.join(f'{x:.5e}' for x in w0['SU2']),flush=True)
    print('Initial Ward U1:', ' '.join(f'{x:.5e}' for x in w0['U1']),flush=True)
    if max(w0['SU2']+w0['U1'])>1e-9:
        raise RuntimeError('Initial corrected tensor fails canonical Ward')
    print('Running full coupled EFT1 tensor ...',flush=True)
    result=run_wilson(basis,initial,bg,mf,ms,
                      adm_variant=args.variant,solver=args.solver,steps=args.steps,
                      rtol=args.rtol,atol=args.atol)
    fullc5=diagnostics('SAVED RUNNING BACKGROUND',result,basis)
    if args.compare_frozen:
        y0=np.repeat(bg.y[-1:,:,:,:],len(mu),axis=0)
        q0=np.repeat(bg.quartic[-1:,:,:,:,:],len(mu),axis=0)
        g0=np.repeat(bg.gauge_couplings[-1:,:],len(mu),axis=0)
        fb=Background(basis,bg.t.copy(),y0,q0,g0,bg.generators,
                      {'mode':'MF background frozen'})
        print('Running frozen-MF comparison ...',flush=True)
        frozen=run_wilson(basis,initial,fb,mf,ms,
                          adm_variant=args.variant,solver=args.solver,steps=args.steps,
                          rtol=args.rtol,atol=args.atol)
        fc5=diagnostics('FROZEN MF BACKGROUND',frozen,basis)
        print('\nFull/frozen ||C5|| ratio:',
              f'{np.linalg.norm(fullc5)/max(np.linalg.norm(fc5),1e-300):.12g}')
        print('Relative matrix difference:',
              f'{np.linalg.norm(fullc5-fc5)/max(np.linalg.norm(fc5),1e-300):.12g}')
    print('\nNo production or saved files modified.')
    print('Candidate ADM only; GammaCC/Weyl sign, scalar threshold and SMEFT NOT certified.')

if __name__=='__main__':
    main()
