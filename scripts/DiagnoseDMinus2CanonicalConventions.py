#!/usr/bin/env python3
"""Read-only D(-2) canonical representation audit at the saved fermion threshold.

Compare three boundary embeddings under the repository's coefficient Ward identity:
  provisional source-identity;
  scalar triplet O rotation only;
  previous epsilon-on-leptons + scalar O rotation.

For the same saved native lambdaT3, compare the source quartic and its
spin-one-metric-J corrected tensor, then evaluate the existing scalar-pair ADM
and project the resulting LLHH block onto the repository Weinberg template.

This is an MF diagnostic, not an RG run, a Weyl normalization certification,
or proof of the full ADM or scalar-threshold matching.
"""
from __future__ import annotations
import argparse
import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from Numerical.parallel_t3.wilson.basis import build_basis
from Numerical.parallel_t3.wilson.coordinates import embed_formal_llss, embed_c5
from Numerical.parallel_t3.wilson.adm import beta, gauge_ward_residual
from Numerical.parallel_t3.rough_inputs import builder
from Numerical.parallel_t3.threshold.matching import project_master, master_weinberg_template

DEFAULT=Path('output/parallel_t3/rough_inputs/fit16_cc2dcbec4a/a/D_m2')

def triplet_rotation():
    s=1/np.sqrt(2)
    U=np.array([[s,0,-1j*s],[0,-1,0],[-s,0,-1j*s]],complex)
    O=np.zeros((6,6),float)
    for a in range(3):
        for b in range(3):
            z=U[a,b]
            O[2*a:2*a+2,2*b:2*b+2]=((z.real,-z.imag),(z.imag,z.real))
    assert np.linalg.norm(O.T@O-np.eye(6))<1e-12
    R=np.eye(8);R[:6,:6]=O.T
    return R

def corrected_portal(ns):
    G=builder._cg()
    J=np.array([[0,0,1],[0,-1,0],[1,0,0]],float)
    invariant=np.einsum('aij,ak->ijk',G,J,optimize=True)
    portal=np.zeros((ns,)*4,complex)
    for i,j,k in np.argwhere(np.abs(invariant)>1e-14):
        portal+=invariant[i,j,k]*builder._monomial(
            builder._u(ns,0,int(i)),builder._u(ns,0,int(j)),
            builder._u(ns,4,int(k)),builder._u(ns,10,0,True))
    return portal

def read_inputs(folder):
    files={
        'boundary.npz':('C_LLSS_formal_muF','C5_hard'),
        'background.npz':('mu_gev','Y','lambda4','g2','gY','g3'),
        'saved_native_trajectory.npz':('mu_eft1_gev','eft1__lambdaT3')}
    ret={}
    for fname,keys in files.items():
        path=folder/fname
        with np.load(path,allow_pickle=False) as data:
            for k in keys:
                if k not in data:raise KeyError(f'{path}: missing {k}; has {data.files}')
                ret[k]=data[k].copy()
    return ret

def llhh(full,basis):
    leptons=[f for weak in basis.leptons for f in weak]
    block=full[np.ix_(leptons,leptons,range(4),range(4))]
    return block.reshape(3,2,3,2,4,4).transpose(0,2,1,3,4,5)

def ward(full,basis):
    return {
        group:[gauge_ward_residual(full,tf,ts) for tf,ts in basis.generators()[group]]
        for group in ('SU2','U1')
    }

def report_ward(label, w):
    print(label+' SU2 max-entry:', ' '.join(f'{x:.9e}' for x in w['SU2']))
    print(label+' U1 max-entry: ', ' '.join(f'{x:.9e}' for x in w['U1']))

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',type=Path,default=DEFAULT)
    p.add_argument('--variant',choices=['fonseca_eq22_six_term','ward_covariant_candidate'],
                   default='fonseca_eq22_six_term')
    p.add_argument('--tolerance',type=float,default=1e-9)
    args=p.parse_args()
    if args.tolerance<=0:p.error('tolerance must be positive')
    d=read_inputs(args.input)
    basis=build_basis('D_-2')
    raw=np.asarray(d['C_LLSS_formal_muF'],complex)
    if raw.shape!=(3,3,2,2,8,8):raise ValueError(f'Unexpected formal tensor {raw.shape}')
    c5=np.asarray(d['C5_hard'],complex)
    if c5.shape!=(3,3):raise ValueError(f'Unexpected C5_hard {c5.shape}')
    R=triplet_rotation()
    eps=np.array([[0,1],[-1,0]],complex)
    epsraw=np.einsum('iI,jJ,pqIJab->pqijab',eps,eps,raw,optimize=True)
    candidates={
        'PROVISIONAL_identity':embed_formal_llss(raw,basis,np.eye(8),1.),
        'TRIPLET_ONLY':embed_formal_llss(raw,basis,R,1.),
        'EPSILON_AND_TRIPLET':embed_formal_llss(epsraw,basis,R,1.),
    }
    for k in candidates:
        candidates[k]+=embed_c5(c5,basis)
    times=np.asarray(d['mu_gev'],float)
    times_native=np.asarray(d['mu_eft1_gev'],float)
    m=int(np.argmax(times));n=int(np.argmax(times_native))
    if not np.isclose(times[m],times_native[n],atol=1e-6,rtol=1e-10):
        raise ValueError('Saved MF endpoints differ')
    y=np.asarray(d['Y'][m],complex)
    old=np.asarray(d['lambda4'][m],float)
    z=complex(np.asarray(d['eft1__lambdaT3'][n]).item())
    old_unit=builder.quartic_unit_tensors('D_-2')['lambdaT3']
    new_unit=corrected_portal(basis.ns)
    new=old+2*np.real(z*(new_unit-old_unit))
    gauge=tuple((float(d[key][m]),basis.generators()[group]) for key,group in
                (('g2','SU2'),('gY','U1'),('g3','SU3')))

    print('Input:',args.input)
    print('MF [GeV]:',times[m], 'native lambdaT3:',z)
    print('Wilson basis: original (identity), triplet O only, epsilon+triplet O')
    print('Repository Ward uses coefficient-slot transposed generator.')
    print('Old/new quartic norm:',f'{np.linalg.norm(old):.12e}',
          f'{np.linalg.norm(new):.12e}')
    # A physically correct Weinberg seed must pass the same Ward convention.
    T=master_weinberg_template()
    ll_indices=[f for weak in basis.leptons for f in weak]
    for group in ('SU2','U1'):
        ref=[]
        for tf,ts in basis.generators()[group]:
            lf=tf[np.ix_(basis.leptons[0],basis.leptons[0])]
            # Ward coefficient slot indices exactly as in adm.gauge_ward_residual
            var=(np.einsum('ki,kjab->ijab',lf,T,optimize=True)+
                 np.einsum('kj,ikab->ijab',lf,T,optimize=True)+
                 np.einsum('ca,ijcb->ijab',ts[:4,:4],T,optimize=True)+
                 np.einsum('cb,ijac->ijab',ts[:4,:4],T,optimize=True))
            ref.append(np.linalg.norm(var)/np.linalg.norm(T))
        print('Weinberg template '+group+' Ward:',*[f'{v:.4e}' for v in ref])
    if max(ref)>args.tolerance:
        raise RuntimeError('Reference template Ward failed: inconsistent conventions')

    print('\nCandidate results, using existing repository scalar_pair beta')
    results={}
    for name,C in candidates.items():
        # Use the same J-corrected lambda in both transformed alternatives.
        lam=old if name=='PROVISIONAL_identity' else new
        bc,backend=beta(C,y,lam,gauge,variant=args.variant)
        pair=bc['scalar_pair']
        ll=llhh(pair,basis)
        c5_beta, residual=project_master(ll,enforce_pure=False)
        w=ward(C,basis)
        wp=ward(pair,basis)
        print('\n'+name+' | engine: '+backend)
        print('C norm:',f'{np.linalg.norm(C):.12e}')
        report_ward('C',w)
        report_ward('scalar_pair beta',wp)
        print('scalar_pair LLHH norm:',f'{np.linalg.norm(ll):.12e}')
        print('Weinberg projection residual:',f'{residual:.12e}')
        print('projected C5 beta norm:',f'{np.linalg.norm(c5_beta):.12e}')
        results[name]=(w,wp,residual)
    w,wp,r=results['TRIPLET_ONLY']
    gauge_ok=(max(w['SU2']+w['U1']+wp['SU2']+wp['U1'])<args.tolerance)
    c5_ok=r<args.tolerance
    print('\nTRIPLET_ONLY repository-Ward covariance:', 'PASS' if gauge_ok else 'FAIL')
    print('TRIPLET_ONLY scalar_pair LLHH Weinberg projection:', 'PASS' if c5_ok else 'FAIL')
    print('No files modified. The comparison is restricted to MF scalar-pair mixing.')
    print('Not certified: GammaCC/Weyl signed normalization, complete ADM, scalar matching.')
    if not (gauge_ok and c5_ok):raise SystemExit(2)

if __name__=='__main__':main()
