#!/usr/bin/env python3
"""Read-only D(-2) source -> canonical SU(2) triplet Ward diagnostic.

Run from repository root. Optional --npz must contain C_LLSS_formal of shape
(3,3,2,2,8,8), with source real order (triplet R/I interleaved,
T3Scalar2Minus, T3Scalar2Plus). Does not modify physics or files.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np


def spin_generators(n: int):
    s=(n-1)/2
    m=[s-j for j in range(n)]
    plus=np.zeros((n,n),complex)
    for j,x in enumerate(m):
        if x+1 in m:
            plus[m.index(x+1),j]=np.sqrt((s-x)*(s+x+1))
    return ((plus+plus.conj().T)/2,(plus-plus.conj().T)/(2j),np.diag(m).astype(complex))


def real_scalar_generator(t):
    n=t.shape[0]
    out=np.zeros((2*n,2*n),complex)
    for a in range(n):
        for b in range(n):
            z=t[a,b]
            out[2*a:2*a+2,2*b:2*b+2]=1j*np.array([[z.imag,z.real],[-z.real,z.imag]])
    return out


def complex_to_real(u):
    n=u.shape[0]
    o=np.zeros((2*n,2*n),float)
    for a in range(n):
        for b in range(n):
            z=u[a,b]
            o[2*a:2*a+2,2*b:2*b+2]=[[z.real,-z.imag],[z.imag,z.real]]
    return o


def get_u(source, canonical):
    # vec(U) in column-major order; T_can U = U T_src.
    a=np.vstack([np.kron(np.eye(3),tc)-np.kron(ts.T,np.eye(3))
                 for tc,ts in zip(canonical,source)])
    _,singular,vh=np.linalg.svd(a)
    u=vh[-1].conj().reshape((3,3),order='F')
    u=u/np.sqrt(np.trace(u.conj().T@u).real/3)
    # Fix only the arbitrary global phase, make u[0,0] positive real.
    u*=np.exp(-1j*np.angle(u[0,0]))
    return u,singular


def ward_residual(c,lepton,scalar):
    # c[p,q,i,j,a,b]; generators act on all four representation slots.
    d=(np.einsum('ik,pqkjab->pqijab',lepton,c,optimize=True)
       +np.einsum('jk,pqikab->pqijab',lepton,c,optimize=True)
       +np.einsum('ac,pqijcb->pqijab',scalar,c,optimize=True)
       +np.einsum('bc,pqijac->pqijab',scalar,c,optimize=True))
    return np.linalg.norm(d)/max(np.linalg.norm(c),1e-300)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--npz',type=Path,help='Retained formal source tensor .npz')
    parser.add_argument('--key',default='C_LLSS_formal')
    args=parser.parse_args()
    j=1j
    m=[np.array([[0,0,0],[0,0,-j],[0,j,0]],complex),
       np.array([[0,-j,0],[j,0,0],[0,0,0]],complex),
       np.array([[0,0,j],[0,0,0],[-j,0,0]],complex)]
    src=[-x for x in m]
    can=spin_generators(3)
    u,sv=get_u(src,can)
    o=complex_to_real(u)
    r=o.T # phi_source_triplet = r @ phi_canonical_triplet
    print('U: canonical_complex = U @ source_complex')
    print(np.array2string(u,precision=12,suppress_small=True))
    print('O: canonical_interleaved_real = O @ source_interleaved_real')
    print(np.array2string(o,precision=12,suppress_small=True))
    print('R = O.T: source_interleaved_real = R @ canonical_interleaved_real')
    print(np.array2string(r,precision=12,suppress_small=True))
    print('smallest singular value:',sv[-1])
    print('U unitarity Frobenius:',np.linalg.norm(u.conj().T@u-np.eye(3)))
    print('O orthogonality Frobenius:',np.linalg.norm(o.T@o-np.eye(6)))
    for k,(ts,tc) in enumerate(zip(src,can),1):
        print(f'T{k} complex intertwining: {np.linalg.norm(tc@u-u@ts):.6e}')
        print(f'T{k} real intertwining: {np.linalg.norm(real_scalar_generator(tc)@o-o@real_scalar_generator(ts)):.6e}')
    if args.npz is None:
        print('Saved Wilson tensor NOT tested: specify --npz path to formal tensor.')
        return
    with np.load(args.npz,allow_pickle=False) as f:
        if args.key not in f:
            raise KeyError(f'{args.key} not found; available keys: {f.files}')
        c=np.asarray(f[args.key],complex)
    if c.shape != (3,3,2,2,8,8):
        raise ValueError(f'Expected (3,3,2,2,8,8), received {c.shape}')
    # Full heavy scalar map including two SU(2)-singlet real slots.
    o8=np.eye(8)
    o8[:6,:6]=o
    transformed=np.einsum('pqijab,Aa,Bb->pqijAB',c,o8,o8,optimize=True)
    gsrc=[]; gcan=[]
    for ts,tc in zip(src,can):
        gs=np.zeros((8,8),complex); gc=np.zeros((8,8),complex)
        gs[:6,:6]=real_scalar_generator(ts)
        gc[:6,:6]=real_scalar_generator(tc)
        gsrc.append(gs);gcan.append(gc)
    print('Wilson Frobenius norm source/canonical:',np.linalg.norm(c),np.linalg.norm(transformed))
    for k,(l,gs,gc) in enumerate(zip(spin_generators(2),gsrc,gcan),1):
        a=ward_residual(c,l,gs)
        b=ward_residual(transformed,l,gc)
        print(f'T{k}: source Ward={a:.12e}; canonical transformed Ward={b:.12e}; difference={abs(a-b):.3e}')
    print('No files were modified.')


if __name__=='__main__':
    main()
