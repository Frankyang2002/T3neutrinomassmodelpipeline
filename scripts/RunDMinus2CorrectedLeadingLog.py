#!/usr/bin/env python3
"""D(-2) corrected MF leading-log EFT1 Weinberg mixing diagnostic.

Read-only original benchmark; no Yukawa fit or RG ODE. Evaluates the installed
repo candidate ADM at MF using the SU2/U1-Ward-covariant TRIPLET_ONLY boundary,
the J-corrected D(-2) quartic, and the unchanged MF Yukawa/gauge background.

C5(MS)_LL = C5(MF) + ln(MS/MF)/(16pi^2) * beta_C5(MF).
Reports each ADM channel, scalar_pair contribution, total candidate beta,
Weinberg template residuals and Ward tests. This first-order estimate is NOT
a substitute for full running; signed matching conventions remain uncertified.
"""
from __future__ import annotations
import argparse
import math
import sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
sys.path.insert(0,str(ROOT/'scripts'))

from DiagnoseDMinus2CanonicalConventions import (
    DEFAULT, triplet_rotation, corrected_portal, read_inputs, llhh, ward
)
from Numerical.parallel_t3.wilson.basis import build_basis
from Numerical.parallel_t3.wilson.coordinates import (
    embed_formal_llss,embed_c5,project_c5
)
from Numerical.parallel_t3.wilson.adm import beta,total_beta
from Numerical.parallel_t3.rough_inputs import builder
from Numerical.parallel_t3.threshold.matching import project_master

CHANNELS=(
 'yukawa_wavefunction','scalar_pair','mixed_yukawa_gauge',
 'crossed_yukawa','conjugate_coefficient_yukawa',
 'scalar_anomalous_dimension','fermion_anomalous_dimension'
)

def check_ward(x,basis,tolerance,label):
    residual=ward(x,basis)
    maxval=max(residual['SU2']+residual['U1'])
    print(f'{label} Ward max-entry SU2:', ' '.join(f'{v:.4e}' for v in residual['SU2']))
    print(f'{label} Ward max-entry U1:', ' '.join(f'{v:.4e}' for v in residual['U1']))
    print(f'{label} Ward max-entry overall: {maxval:.8e}')
    return maxval<tolerance

def summarize(label,x,basis):
    ll=llhh(x,basis)
    c5,res=project_master(ll,enforce_pure=False)
    print(f'{label:37s} ||beta||={np.linalg.norm(x):.8e} '
          f'||LLHH||={np.linalg.norm(ll):.8e} '
          f'||beta_C5||={np.linalg.norm(c5):.8e} '
          f'r_Weinberg={res:.8e}')
    return c5,res

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,default=DEFAULT)
    parser.add_argument('--ms',type=float,default=1e3,
                        help='MS in GeV; affects first-order leading log only')
    parser.add_argument('--tolerance',type=float,default=1e-9)
    parser.add_argument('--variant',default='fonseca_eq22_six_term',
                        choices=['fonseca_eq22_six_term','ward_covariant_candidate'])
    args=parser.parse_args()
    d=read_inputs(args.input)
    basis=build_basis('D_-2')
    t=np.asarray(d['mu_gev'],float)
    tn=np.asarray(d['mu_eft1_gev'],float)
    im=int(np.argmax(t)); inn=int(np.argmax(tn))
    mf=float(t[im])
    if not np.isclose(mf,tn[inn],rtol=1e-10,atol=1e-6):
        raise ValueError('MF mismatch between saved arrays')
    if not (0<args.ms<mf):raise ValueError('Require 0 < MS < MF')
    raw=np.asarray(d['C_LLSS_formal_muF'],complex)
    if raw.shape!=(3,3,2,2,8,8):raise ValueError(f'Unexpected raw shape {raw.shape}')
    C=embed_formal_llss(raw,basis,triplet_rotation(),1.)
    C+=embed_c5(np.asarray(d['C5_hard'],complex),basis)
    y=np.asarray(d['Y'][im],complex)
    z=complex(np.asarray(d['eft1__lambdaT3'][inn]).item())
    units=builder.quartic_unit_tensors('D_-2')
    portal_old=units['lambdaT3']
    portal_new=corrected_portal(basis.ns)
    lam=np.asarray(d['lambda4'][im],float)+2*np.real(z*(portal_new-portal_old))
    gauge=tuple((float(d[key][im]),basis.generators()[group]) for key,group in
                (('g2','SU2'),('gY','U1'),('g3','SU3')))
    print('=== D(-2) corrected scalar-pair / total beta leading log ===')
    print('Input:',args.input)
    print(f'MF={mf:g} GeV, MS={args.ms:g} GeV, ln(MS/MF)={math.log(args.ms/mf):.12g}')
    print('Mapping: triplet O only, NO lepton epsilon; J-corrected quartic.')
    print('Native lambdaT3(MF):',z)
    print('GammaCC/Weyl factor: provisional +1 (not certified)')
    if not check_ward(C,basis,args.tolerance,'Initial full Wilson'):
        raise RuntimeError('Corrected boundary fails repository SU2/U1 Ward')
    if not np.allclose(lam,lam.swapaxes(0,1),atol=1e-12):
        raise RuntimeError('Corrected lambda not symmetric')
    channels,engine=beta(C,y,lam,gauge,variant=args.variant)
    print('ADM engine:',engine,'variant:',args.variant)
    print('All channel units: 16*pi^2*dC/dln(mu), at MF')
    parsed={}
    for key in CHANNELS:
        parsed[key]=summarize(key,channels[key],basis)
    total=total_beta(channels)
    c5_total,res_total=summarize('TOTAL candidate beta',total,basis)
    c5_pair,res_pair=parsed['scalar_pair']
    total_ok=check_ward(total,basis,args.tolerance,'Total candidate beta')
    pair_ok=check_ward(channels['scalar_pair'],basis,args.tolerance,'Scalar-pair beta')
    c5_start,res_start=project_c5(C,basis)
    if res_start>args.tolerance and np.linalg.norm(c5_start)>1e-25:
        print('WARNING: initial C5 boundary LLHH template not pure')
    delta_log=math.log(args.ms/mf)/(16*math.pi**2)
    pair_c5=delta_log*c5_pair
    total_c5=c5_start+delta_log*c5_total
    print('\nLeading-log coefficient:',f'{delta_log:.12e}')
    print('||C5(MF)|| [GeV^-1]:',f'{np.linalg.norm(c5_start):.12e}')
    print('||Delta C5_pair|| [GeV^-1]:',f'{np.linalg.norm(pair_c5):.12e}')
    print('||C5(MS) first-order total|| [GeV^-1]:',
          f'{np.linalg.norm(total_c5):.12e}')
    print('Scalar-pair Delta C5 matrix [GeV^-1]:')
    print(np.array2string(pair_c5,precision=8))
    print('Total first-order C5(MS) matrix [GeV^-1]:')
    print(np.array2string(total_c5,precision=8))
    print('NOTE: first-order local leading log, no running of C, Y, lambda or gauge.')
    print('NOTE: complete candidate ADM and signed GammaCC/Weyl normalization unverified.')
    print('NOTE: no scalar threshold integration or SMEFT evolution; no neutrino masses computed.')
    print('No files modified.')
    if not pair_ok or not total_ok or res_pair>args.tolerance:
        print('WARNING: candidate ADM fails at least one Ward/Weinberg test.')
        raise SystemExit(2)

if __name__=='__main__':
    main()
