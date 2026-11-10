#!/usr/bin/env python3
"""Diagnose residual LLHH orientation after D(-2) SU2-covariant correction.

Read-only; uses saved MF inputs, the exact existing candidate ADM, and
the prior diagnostic's proven source-to-canonical mapping and J portal.
Reports SU2/U1 Ward residuals, overlap with the repository Weinberg template,
its complex conjugate and other scalar-sector orientations.
Does not refit or run RGEs. Does not certify loop matching.
"""
from __future__ import annotations
import argparse
import sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
sys.path.insert(0,str(ROOT))
import DiagnoseDMinus2ScalarPair as d
from Numerical.parallel_t3.wilson.basis import build_basis
from Numerical.parallel_t3.wilson.coordinates import embed_formal_llss,embed_c5
from Numerical.parallel_t3.wilson.adm import beta
from Numerical.parallel_t3.threshold.matching import (
    master_weinberg_template,project_master
)

def ward(x, tf,ts):
    dx=(np.einsum('ik,kjab->ijab',tf,x,optimize=True)
       +np.einsum('jk,ikab->ijab',tf,x,optimize=True)
       +np.einsum('ac,ijcb->ijab',ts,x,optimize=True)
       +np.einsum('bc,ijac->ijab',ts,x,optimize=True))
    return float(np.linalg.norm(dx)/max(np.linalg.norm(x),1e-300))

def main():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument('--input',type=Path,default=d.DEFAULT)
    args=a.parse_args()
    data=d.load(args.input)
    b=build_basis('D_-2')
    raw=np.asarray(data['C_LLSS_formal_muF'],complex)
    eps=np.array([[0,1],[-1,0]],complex)
    transformed=np.einsum('iI,jJ,pqIJab->pqijab',eps,eps,raw,optimize=True)
    R=np.eye(8);R[:6,:6]=d.triplet_O().T
    C=embed_formal_llss(transformed,b,R,1.)+embed_c5(data['C5_hard'],b)
    mu=np.asarray(data['mu_gev'],float)
    idx=int(np.argmax(mu))
    mun=np.asarray(data['mu_eft1_gev'],float)
    ix=int(np.argmax(mun))
    if not np.isclose(mu[idx],mun[ix],rtol=1e-10,atol=1e-6):
        raise ValueError('MF scale mismatch')
    z=complex(np.asarray(data['eft1__lambdaT3'][ix]).item())
    old=d.builder.quartic_unit_tensors('D_-2')['lambdaT3']
    new=d.d_portal(b.ns)
    lam=np.asarray(data['lambda4'][idx],float)+2*np.real(z*(new-old))
    y=np.asarray(data['Y'][idx],complex)
    gauge=tuple((float(data[name][idx]),b.generators()[group])
               for name,group in [('g2','SU2'),('gY','U1'),('g3','SU3')])
    terms,backend=beta(C,y,lam,gauge,variant='fonseca_eq22_six_term')
    ll=d.to_llhh(terms['scalar_pair'],b)
    T=master_weinberg_template()
    Tstar=T.conj()
    # Alternative orientations have H1/H2 swapped and signs of scalar imaginary
    # components changed; only exploratory diagnostics, not alternate physical definitions.
    flipped=np.array([2,3,0,1])
    templates={
        'repo_Weinberg':T,
        'conjugate_Weinberg':Tstar,
        'weak_flipped_Weinberg':T[:, :, flipped][:, :, :, flipped],
        'weak_flipped_conjugate':Tstar[:,:,flipped][:,:,:,flipped],
    }
    print('D(-2) LLHH representation diagnostic at MF',mu[idx], 'GeV')
    print('ADM backend:',backend)
    print('LLHH scalar-pair beta norm:',f'{np.linalg.norm(ll):.12e}')
    print('Full corrected C SU2/U1 Ward:',
          [f'{ward(C,tf,ts):.6e}' for group in ('SU2','U1')
           for tf,ts in b.generators()[group]])
    print('Full corrected scalar_pair beta SU2/U1 Ward:',
          [f'{ward(terms["scalar_pair"],tf,ts):.6e}' for group in ('SU2','U1')
           for tf,ts in b.generators()[group]])
    # Select one flavour for conventional SU2/U1 Ward on LLHH
    for group in ('SU2','U1'):
        gg=b.generators()[group]
        lepton=list(b.leptons[1])
        print(group,'LLHH residuals:',
              [f'{ward(ll[1,1],tf[np.ix_(lepton,lepton)],ts[:4,:4]):.6e}'
               for tf,ts in gg])
    C5,res=project_master(ll,enforce_pure=False)
    print('Repository Weinberg projector relative residual:',f'{res:.12e}')
    for name,t in templates.items():
        x=np.einsum('pqijab,ijab->pq',ll,t.conj(),optimize=True)/np.vdot(t,t)
        rem=ll-np.einsum('pq,ijab->pqijab',x,t,optimize=True)
        r=np.linalg.norm(rem)/max(np.linalg.norm(ll),1e-300)
        print(f'{name}: projection residual={r:.12e}, coefficient norm={np.linalg.norm(x):.12e}')
    print('Conjugate-W template overlap:',f'{abs(np.vdot(T,Tstar))/(np.linalg.norm(T)*np.linalg.norm(Tstar)):.12e}')
    print('C5 projected [1/GeV]:')
    print(np.array2string(C5,precision=5))
    print('No files modified.')
    print('Interpretation: SU2 covariance alone does not establish Weinberg orientation;')
    print('check U1 and source scalar conjugation/real-component conventions.')

if __name__=='__main__':
    main()
