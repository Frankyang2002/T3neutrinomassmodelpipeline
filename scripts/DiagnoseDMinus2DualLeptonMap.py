#!/usr/bin/env python3
"""Read-only D(-2) formal Wilson source-dual to canonical-fundamental SU(2) check.
No physics files modified. Does not certify Weyl/Matchete normalization or EFT running.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np


def spin(n):
    s=(n-1)/2
    m=[s-i for i in range(n)]
    up=np.zeros((n,n),complex)
    for j,x in enumerate(m):
        if x+1 in m:
            up[m.index(x+1),j]=np.sqrt((s-x)*(s+x+1))
    return ((up+up.conj().T)/2,(up-up.conj().T)/(2j),np.diag(m).astype(complex))


def real_generator(t):
    g=np.zeros((2*len(t),2*len(t)),complex)
    for a in range(len(t)):
        for b in range(len(t)):
            z=t[a,b]
            g[2*a:2*a+2,2*b:2*b+2]=1j*np.array([[z.imag,z.real],[-z.real,z.imag]])
    return g


def triplet_map():
    j=1j
    m=[np.array([[0,0,0],[0,0,-j],[0,j,0]],complex),
       np.array([[0,-j,0],[j,0,0],[0,0,0]],complex),
       np.array([[0,0,j],[0,0,0],[-j,0,0]],complex)]
    source=[-x for x in m]
    canonical=spin(3)
    a=np.vstack([np.kron(np.eye(3),tc)-np.kron(ts.T,np.eye(3)) for tc,ts in zip(canonical,source)])
    _,sv,vh=np.linalg.svd(a)
    u=vh[-1].conj().reshape((3,3),order='F')
    u/=np.sqrt(np.trace(u.conj().T@u).real/3)
    u*=np.exp(-1j*np.angle(u[0,0]))
    o=np.zeros((6,6),float)
    for i in range(3):
        for k in range(3):
            z=u[i,k]
            o[2*i:2*i+2,2*k:2*k+2]=[[z.real,-z.imag],[z.imag,z.real]]
    return source,canonical,o,sv[-1]


def ward(c,l,s):
    d=(np.einsum('ik,pqkjab->pqijab',l,c,optimize=True)
      +np.einsum('jk,pqikab->pqijab',l,c,optimize=True)
      +np.einsum('ac,pqijcb->pqijab',s,c,optimize=True)
      +np.einsum('bc,pqijac->pqijab',s,c,optimize=True))
    return np.linalg.norm(d)/max(np.linalg.norm(c),1e-300)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--npz',type=Path,default=Path('output/parallel_t3/rough_inputs/fit16_cc2dcbec4a/a/D_m2/boundary.npz'))
    parser.add_argument('--key',default='C_LLSS_formal_muF')
    args=parser.parse_args()
    with np.load(args.npz,allow_pickle=False) as f:
        if args.key not in f:
            parser.error(f'Missing {args.key}; keys: {list(f.files)}')
        c=np.asarray(f[args.key],complex)
    if c.shape!=(3,3,2,2,8,8):
        parser.error(f'Expected (3,3,2,2,8,8), got {c.shape}')
    src,can,o,sv=triplet_map()
    o8=np.eye(8)
    o8[:6,:6]=o
    eps=np.array([[0,1],[-1,0]],complex)
    l=spin(2)
    dual=tuple(-t.T for t in l)
    mapped=np.einsum('iI,jJ,Aa,Bb,pqIJab->pqijAB',eps,eps,o8,o8,c,optimize=True)
    print('Input:',args.npz,'key:',args.key)
    print('Source/mapped Wilson norms:',f'{np.linalg.norm(c):.12e}',f'{np.linalg.norm(mapped):.12e}')
    print('Epsilon dual intertwining maximum:',f'{max(np.linalg.norm(t@eps-eps@d) for t,d in zip(l,dual)):.12e}')
    print('Triplet intertwiner smallest singular value:',f'{sv:.12e}')
    fixed=[]
    for k,(ld,lf,ts,tc) in enumerate(zip(dual,l,src,can),1):
        gs=np.zeros((8,8),complex); gs[:6,:6]=real_generator(ts)
        gc=np.zeros((8,8),complex); gc[:6,:6]=real_generator(tc)
        a,b=ward(c,ld,gs),ward(mapped,lf,gc)
        fixed.append(b)
        print(f'T{k}: source dual Ward={a:.12e}; mapped canonical fundamental Ward={b:.12e}')
    print('Exchange residual lepton/scalar:',
          f'{np.linalg.norm(mapped-mapped.transpose(1,0,3,2,4,5))/max(np.linalg.norm(mapped),1e-300):.12e}',
          f'{np.linalg.norm(mapped-mapped.swapaxes(-1,-2))/max(np.linalg.norm(mapped),1e-300):.12e}')
    print('SU(2) covariance:', 'PASS' if max(fixed)<1e-10 else 'FAIL')
    print('Not certified: source index meaning, singlet phase/order, U(1), Weyl normalization, quartic matching.')
    print('No saved files modified.')


if __name__=='__main__':
    main()
