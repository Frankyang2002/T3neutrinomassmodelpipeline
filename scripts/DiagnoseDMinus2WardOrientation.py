#!/usr/bin/env python3
"""D(-2) read-only Ward-index audit for MF Wilson, scalar portal, beta and Weinberg.

Tests the repository's gauge_ward_residual (transposed generator index
placement) against the previous direct-index diagnostic. A conventional
Weinberg template is an independent zero-residual reference.

Uses existing scripts/DiagnoseDMinus2ScalarPair.py as its companion.
No saved inputs are modified, no running or fitting is performed.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
sys.path.insert(0,str(ROOT/'scripts'))
import DiagnoseDMinus2ScalarPair as old
from Numerical.parallel_t3.wilson.basis import build_basis
from Numerical.parallel_t3.wilson.coordinates import embed_formal_llss,embed_c5
from Numerical.parallel_t3.wilson.adm import beta,gauge_ward_residual
from Numerical.parallel_t3.threshold.matching import master_weinberg_template,project_master


def ward(c, f, s, placement):
    if placement=='repository':
        ff=f.T
        ss=s.T
    elif placement=='previous':
        ff=f
        ss=s
    else:
        raise ValueError(placement)
    dx=(np.einsum('ik,kjab->ijab',ff,c,optimize=True)
       +np.einsum('jk,ikab->ijab',ff,c,optimize=True)
       +np.einsum('ac,ijcb->ijab',ss,c,optimize=True)
       +np.einsum('bc,ijac->ijab',ss,c,optimize=True))
    return float(np.linalg.norm(dx)/max(np.linalg.norm(c),1e-300))


def qward(x,s,placement):
    g=s.T if placement=='repository' else s
    d=(np.einsum('ac,cbde->abde',g,x,optimize=True)
      +np.einsum('bc,acde->abde',g,x,optimize=True)
      +np.einsum('dc,abce->abde',g,x,optimize=True)
      +np.einsum('ec,abdc->abde',g,x,optimize=True))
    return float(np.linalg.norm(d)/max(np.linalg.norm(x),1e-300))


def all_ward(name,c,gens):
    print(name,'norm',f'{np.linalg.norm(c):.9e}')
    for group in ('SU2','U1'):
        w_old=[ward(c,f,s,'previous') for f,s in gens[group]]
        w_correct=[ward(c,f,s,'repository') for f,s in gens[group]]
        w_repo=[gauge_ward_residual(c,f,s) for f,s in gens[group]]
        print(f'  {group} old Frobenius       :',*[f'{x:.8e}' for x in w_old])
        print(f'  {group} repository Frobenius:',*[f'{x:.8e}' for x in w_correct])
        print(f'  {group} repo max-entry      :',*[f'{x:.8e}' for x in w_repo])


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input',type=Path,default=old.DEFAULT)
    args=ap.parse_args()
    b=build_basis('D_-2')
    D=old.load(args.input)
    raw=np.asarray(D['C_LLSS_formal_muF'],complex)
    eps=np.array([[0,1],[-1,0]],complex)
    raw_e=np.einsum('iI,jJ,pqIJab->pqijab',eps,eps,raw,optimize=True)
    R=np.eye(8);R[:6,:6]=old.triplet_O().T
    C0=embed_formal_llss(raw,b,np.eye(8),1.)+embed_c5(D['C5_hard'],b)
    C1=embed_formal_llss(raw_e,b,R,1.)+embed_c5(D['C5_hard'],b)
    bgmu=np.asarray(D['mu_gev'],float)
    ni=int(np.argmax(bgmu))
    native_mu=np.asarray(D['mu_eft1_gev'],float)
    nn=int(np.argmax(native_mu))
    if not np.isclose(bgmu[ni],native_mu[nn],rtol=1e-10):
        raise ValueError('Background and native MF mismatch')
    z=complex(np.asarray(D['eft1__lambdaT3'][nn]).item())
    unit0=old.builder.quartic_unit_tensors('D_-2')['lambdaT3']
    unit1=old.d_portal(b.ns)
    lam0=np.asarray(D['lambda4'][ni],float)
    lam1=lam0+2*np.real(z*(unit1-unit0))
    gens=b.generators()
    print('D(-2) Ward convention audit at MF [GeV]:',bgmu[ni])
    print('Conventions: previous = generator [a,c], repo = generator [c,a] on coefficient slots')
    all_ward('PROVISIONAL Wilson',C0,gens)
    all_ward('EPSILON + TRIPLET ROTATION Wilson',C1,gens)

    for nm,l in [('PROVISIONAL quartic',lam0),('J-PORTAL quartic',lam1),
                 ('Original unit portal',unit0),('J unit portal',unit1)]:
        print(nm,'norm',f'{np.linalg.norm(l):.9e}')
        for group in ('SU2','U1'):
            print(' ',group,'previous:',*[f'{qward(l,s,"previous"):.8e}' for _,s in gens[group]])
            print(' ',group,'repo placement:',*[f'{qward(l,s,"repository"):.8e}' for _,s in gens[group]])

    y=np.asarray(D['Y'][ni],complex)
    gauge=tuple((float(D[name][ni]),gens[group]) for name,group in
                (('g2','SU2'),('gY','U1'),('g3','SU3')))
    for nm,C,L in [('PROVISIONAL',C0,lam0),('EPSILON+J',C1,lam1)]:
        channels,engine=beta(C,y,L,gauge,variant='fonseca_eq22_six_term')
        pair=channels['scalar_pair']
        all_ward(nm+' scalar-pair beta ('+engine+')',pair,gens)
        ll=old.to_llhh(pair,b)
        _,res=project_master(ll,enforce_pure=False)
        print(nm,'LLHH Weinberg projection residual:',f'{res:.12e}')

    T=master_weinberg_template()
    lepton=np.array([-0.5,-0.5])
    f=np.diag(lepton).astype(complex)
    print('REFERENCE: repository Weinberg template, U1 Ward')
    for group in ('SU2','U1'):
        for k,(ff,ss) in enumerate(gens[group]):
            l=list(b.leptons[0]); tf=ff[np.ix_(l,l)]
            sc=ss[:4,:4]
            print(group,k+1,'previous',f'{ward(T,tf,sc,"previous"):.12e}',
                  'repository',f'{ward(T,tf,sc,"repository"):.12e}')
    print('No files modified. This audit is diagnostic, not EFT certification.')

if __name__=='__main__':
    main()
