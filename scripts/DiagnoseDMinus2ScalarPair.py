#!/usr/bin/env python3
"""Read-only MF scalar-pair ADM test for retained D(-2) benchmark.

Compare the original provisional source embedding + uncorrected quartic against
the epsilon/rotated canonical Wilson and J-corrected quartic. All inputs at MF,
same Yukawa benchmark. No RGE integration, no Yukawa fit, no file writes.

Original scalar_pair channel is evaluated using the installed repository
Numerical.parallel_t3.wilson.adm.beta; physical ADM normalization unverified.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))

from Numerical.parallel_t3.wilson.basis import build_basis
from Numerical.parallel_t3.wilson.coordinates import embed_formal_llss,embed_c5
from Numerical.parallel_t3.rough_inputs import builder
from Numerical.parallel_t3.wilson.adm import beta
from Numerical.parallel_t3.threshold.matching import project_master

DEFAULT=Path('output/parallel_t3/rough_inputs/fit16_cc2dcbec4a/a/D_m2')

def triplet_O():
    # Verified D(-2) complex source-to-canonical intertwiner, with phase fixed.
    s=np.sqrt(0.5)
    u=np.array([[s,0,-1j*s],[0,-1,0],[-s,0,-1j*s]],complex)
    o=np.zeros((6,6),float)
    for a in range(3):
        for b in range(3):
            z=u[a,b]
            o[2*a:2*a+2,2*b:2*b+2]=[[z.real,-z.imag],[z.imag,z.real]]
    if np.linalg.norm(o@o.T-np.eye(6))>1e-12:
        raise RuntimeError('Triplet rotation is not orthogonal')
    return o

def d_portal(ns):
    J=np.array([[0,0,1],[0,-1,0],[1,0,0]],float)
    G=builder._cg()
    T=np.einsum('aij,ak->ijk',G,J,optimize=True)
    lam=np.zeros((ns,)*4,complex)
    for i,j,k in np.argwhere(abs(T)>1e-14):
        lam += T[i,j,k]*builder._monomial(
            builder._u(ns,0,int(i)),builder._u(ns,0,int(j)),
            builder._u(ns,4,int(k)),builder._u(ns,10,0,True))
    return lam

def load(folder):
    keys={
        'boundary.npz':['C_LLSS_formal_muF','C5_hard'],
        'background.npz':['mu_gev','Y','lambda4','g2','gY','g3'],
        'saved_native_trajectory.npz':['mu_eft1_gev','eft1__lambdaT3']}
    data={}
    for file,names in keys.items():
        path=folder/file
        if not path.is_file():raise FileNotFoundError(path)
        with np.load(path,allow_pickle=False) as f:
            for name in names:
                if name not in f:raise KeyError(f'{path}: missing {name}, keys={f.files}')
                data[name]=f[name].copy()
    return data

def to_llhh(c,basis):
    lepton=[a for pair in basis.leptons for a in pair]
    x=c[np.ix_(lepton,lepton,range(4),range(4))]
    return x.reshape(3,2,3,2,4,4).transpose(0,2,1,3,4,5)

def ward_psi(c,basis):
    vals=[]
    for tf,ts in basis.generators()['SU2']:
        d=(np.einsum('ik,kjab->ijab',tf,c,optimize=True)
          +np.einsum('jk,ikab->ijab',tf,c,optimize=True)
          +np.einsum('ac,ijcb->ijab',ts,c,optimize=True)
          +np.einsum('bc,ijac->ijab',ts,c,optimize=True))
        vals.append(np.linalg.norm(d)/max(np.linalg.norm(c),1e-300))
    return vals

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,default=DEFAULT)
    parser.add_argument('--variant',default='fonseca_eq22_six_term',
                        choices=['fonseca_eq22_six_term','ward_covariant_candidate'])
    args=parser.parse_args()
    d=load(args.input)
    basis=build_basis('D_-2')
    ns=basis.ns
    raw=np.asarray(d['C_LLSS_formal_muF'],complex)
    if raw.shape!=(3,3,2,2,8,8):raise ValueError(f'Wrong raw shape {raw.shape}')
    c5=np.asarray(d['C5_hard'],complex)
    if c5.shape!=(3,3):raise ValueError('Wrong C5_hard shape')
    eps=np.array([[0,1],[-1,0]],complex)
    dual_to_fund=np.einsum('iI,jJ,pqIJab->pqijab',eps,eps,raw,optimize=True)
    R=np.eye(8)
    R[:6,:6]=triplet_O().T

    C_old=embed_formal_llss(raw,basis,np.eye(8),1.)+embed_c5(c5,basis)
    C_new=embed_formal_llss(dual_to_fund,basis,R,1.)+embed_c5(c5,basis)

    mu=np.asarray(d['mu_gev'],float)
    idx=int(np.argmax(mu))
    mu_native=np.asarray(d['mu_eft1_gev'],float)
    ni=int(np.argmax(mu_native))
    if not np.isclose(mu[idx],mu_native[ni],rtol=1e-10,atol=1e-6):
        raise ValueError('MF scales do not agree')
    y=np.asarray(d['Y'][idx],complex)
    lam_old=np.asarray(d['lambda4'][idx],float)
    z=complex(np.asarray(d['eft1__lambdaT3'][ni]).item())
    portal_old=builder.quartic_unit_tensors('D_-2')['lambdaT3']
    portal_new=d_portal(ns)
    lam_new=lam_old+2*np.real(z*(portal_new-portal_old))
    gauge=tuple((float(d[k][idx]),basis.generators()[group])
                for k,group in (('g2','SU2'),('gY','U1'),('g3','SU3')))

    print('Input:',args.input,'MF [GeV]:',mu[idx])
    print('Model D_-2, ADM variant:',args.variant)
    print('Wilson norms old/new:',np.linalg.norm(C_old),np.linalg.norm(C_new))
    print('Quartic norms old/new:',np.linalg.norm(lam_old),np.linalg.norm(lam_new))
    print('Native lambdaT3 MF:',z)
    print('Wilson initial SU2 Ward old:',ward_psi(C_old,basis))
    print('Wilson initial SU2 Ward corrected:',ward_psi(C_new,basis))

    for label,C,L in (('PROVISIONAL',C_old,lam_old),
                      ('CORRECTED',C_new,lam_new)):
        channels,backend=beta(C,y,L,gauge,variant=args.variant)
        pair=channels['scalar_pair']
        llhh=to_llhh(pair,basis)
        c5_beta,res=project_master(llhh,enforce_pure=False)
        print('\n'+label+' beta backend: '+backend)
        print('scalar_pair beta norm:',f'{np.linalg.norm(pair):.12e}')
        print('LLHH scalar_pair beta norm:',f'{np.linalg.norm(llhh):.12e}')
        print('LLHH Weinberg C5 beta norm:',f'{np.linalg.norm(c5_beta):.12e}')
        print('LLHH non-Weinberg projection residual:',f'{res:.12e}')
        print('scalar_pair beta SU2 Ward:',*[f'{v:.6e}' for v in ward_psi(pair,basis)])
        print('C5 beta matrix [1/GeV]:')
        print(np.array2string(c5_beta,precision=6,suppress_small=False))
        if label=='CORRECTED':
            good_ward=max(ward_psi(pair,basis)) < 1e-9
            good_projection=res < 1e-9 if np.linalg.norm(llhh)>1e-25 else True
            print('Corrected scalar_pair SU2:', 'PASS' if good_ward else 'FAIL')
            print('Corrected Weinberg tensor:', 'PASS' if good_projection else 'FAIL')
    print('No saved or production files modified.')
    print('Conditional diagnostics only: signed GammaCC/Weyl factor, singlet conventions,')
    print('complete ADM and matching not independently certified.')

if __name__=='__main__':
    main()
