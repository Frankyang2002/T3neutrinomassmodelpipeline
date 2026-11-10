#!/usr/bin/env python3
"""Read-only analysis of D(-2) formal L L S1 S2 SU(2) singlet structure.

Uses the saved Matchete formal tensor (NOT canonical C_full_provisional).
Source heavy order: triplet (R,I) interleaved [0:6], two real singlets [6:8].
The first two weak axes are tested in the declared L doublet representation.
Alternative sign/dual assignments are diagnostics ONLY, not certified conventions.
No files written; no RG running, refitting, or changes to saved physics.
"""
from __future__ import annotations
import argparse
import numpy as np
from pathlib import Path


def spin(n):
    s=(n-1)/2
    vals=[s-i for i in range(n)]
    up=np.zeros((n,n),complex)
    for j,m in enumerate(vals):
        if m+1 in vals:
            up[vals.index(m+1),j]=np.sqrt((s-m)*(s+m+1))
    return ((up+up.conj().T)/2,(up-up.conj().T)/(2j),np.diag(vals).astype(complex))


def real_generator(t):
    n=len(t)
    g=np.zeros((2*n,2*n),complex)
    for a in range(n):
        for b in range(n):
            z=t[a,b]
            g[2*a:2*a+2,2*b:2*b+2]=1j*np.array([[z.imag,z.real],[-z.real,z.imag]])
    return g


def nullspace_projection(gens):
    # Convention: flatten (weak_i,weak_j,triplet_real_a) in C order.
    eyes=(np.eye(2),np.eye(2),np.eye(6))
    maps=[]
    for l1,l2,t in gens:
        maps.append(np.kron(np.kron(l1,eyes[1]),eyes[2])
                    +np.kron(np.kron(eyes[0],l2),eyes[2])
                    +np.kron(np.kron(eyes[0],eyes[1]),t))
    mat=np.concatenate(maps,axis=0)
    _,sv,vh=np.linalg.svd(mat,full_matrices=True)
    cutoff=max(mat.shape)*np.finfo(float).eps*max(sv[0],1.)*10
    rank=int(np.count_nonzero(sv>cutoff))
    basis=vh[rank:,:].conj().T
    return maps,sv,basis


def evaluate(c, generators, label):
    maps,sv,basis=nullspace_projection(generators)
    # Collect both scalar index orderings, each singlet and all flavours.
    # C[p,q,i,j,triplet_a,singlet_b] -> columns of shape (24,N).
    forward=np.moveaxis(c[:,:,:,:,:6,6:8],-1,0) # (2,3,3,2,2,6)
    reverse=np.moveaxis(c[:,:,:,:,6:8,:6],-2,0) # (2,3,3,2,2,6)
    arr=np.stack([forward,reverse],axis=0)
    vectors=arr.reshape(-1,24).T
    total=np.linalg.norm(vectors)
    proj=basis@(basis.conj().T@vectors)
    bad=vectors-proj
    res=np.linalg.norm(bad)/max(total,1e-300)
    wr=[np.linalg.norm(m@vectors)/max(total,1e-300) for m in maps]
    power=np.sum(abs(vectors)**2,axis=0)
    best=int(np.argmax(power))
    b=proj[:,best]
    relative_best=np.linalg.norm(vectors[:,best]-b)/max(np.linalg.norm(vectors[:,best]),1e-300)
    print(f'\n{label}')
    print(f'  invariant dimension (real-triplet complexification): {basis.shape[1]}')
    print(f'  smallest singular values: {sv[-4:]}')
    print(f'  projection residual (all flavour/singlet/order): {res:.12e}')
    print('  Ward residuals T1/T2/T3:', ' '.join(f'{w:.12e}' for w in wr))
    print(f'  strongest slice projection residual: {relative_best:.12e}')
    return res, basis


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--npz',type=Path,default=Path('output/parallel_t3/rough_inputs/fit16_cc2dcbec4a/a/D_m2/boundary.npz'))
    p.add_argument('--key',default='C_LLSS_formal_muF')
    args=p.parse_args()
    with np.load(args.npz,allow_pickle=False) as f:
        if args.key not in f:
            p.error(f'missing {args.key}; keys: {list(f.files)}')
        c=np.asarray(f[args.key],complex)
    if c.shape!=(3,3,2,2,8,8):
        p.error(f'expected formal shape (3,3,2,2,8,8), got {c.shape}')
    print('Input:',args.npz,'key:',args.key)
    print('Total formal tensor norm:',f'{np.linalg.norm(c):.12e}')
    print('S1-triplet/S2-singlet norm:',f'{np.linalg.norm(c[:,:,:,:,:6,6:8]):.12e}')
    print('S2-singlet/S1-triplet norm:',f'{np.linalg.norm(c[:,:,:,:,6:8,:6]):.12e}')
    print('Scalar exchange residual:',f'{np.linalg.norm(c-c.swapaxes(-1,-2))/max(np.linalg.norm(c),1e-300):.12e}')
    print('Lepton exchange residual:',f'{np.linalg.norm(c-c.transpose(1,0,3,2,4,5))/max(np.linalg.norm(c),1e-300):.12e}')
    m=[np.array([[0,0,0],[0,0,-1j],[0,1j,0]],complex),
       np.array([[0,-1j,0],[1j,0,0],[0,0,0]],complex),
       np.array([[0,0,1j],[0,0,0],[-1j,0,0]],complex)]
    l=spin(2)
    assignments={
        'source triplet = -M, L fundamental':(l,l,[-x for x in m]),
        'source triplet = +M, L fundamental':(l,l,m),
        'source triplet = -M, L dual':(tuple(-x.T for x in l),tuple(-x.T for x in l),[-x for x in m]),
        'source triplet = +M, L dual':(tuple(-x.T for x in l),tuple(-x.T for x in l),m),
    }
    outcomes=[]
    for name,(l1,l2,trip) in assignments.items():
        res,b=evaluate(c,[(a,d,real_generator(t)) for a,d,t in zip(l1,l2,trip)],name)
        outcomes.append((res,name))
    best=sorted(outcomes)[0]
    print('\nLowest projection residual among tested assumptions:',best[1],f'{best[0]:.12e}')
    print('Interpretation: alternative generator assignments are NOT validated by this search.')
    print('No saved files modified.')

if __name__=='__main__':
    main()
