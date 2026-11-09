"""SU(2)-Ward-nullspace synthetic LLSS and portal polynomials for all 16.

This constructs *mathematically gauge-invariant fixtures*, NOT native Matchete
or RGBeta Clebsch--Gordan / coupling normalization. Fermion-indexed LLSS
is built from L L S1^dag S2; quartic portal from H H S1 S2^dag+h.c.
All x=(R+iI)/sqrt(2) and SU(2) irreps use the source-inspected generator
convention. The SU(2) invariant is solved rather than guessed for B,C,D,E.
"""
from __future__ import annotations
from itertools import permutations
import math
import numpy as np
from .basis import ModelBasis,spin_generators
from .coordinates import symmetrize


def invariant_two_doublets(d1:int,d2:int,*,bar_first:bool=False,
                           bar_second:bool=False)->np.ndarray:
    """Return invariant T_{ijab} for 2x2xd1xd2, symmetric in i,j.

    The first two axes are identical SU2 doublets, with both *unbarred*.
    The last two axes have the requested anti-representations. The nullspace
    is determined by the generator Ward identities and exchange symmetry;
    non-unique/null topology fails rather than selecting a guessed CG.
    """
    dims=(2,2,d1,d2);n=math.prod(dims)
    gs=spin_generators(2)
    a=spin_generators(d1);b=spin_generators(d2)
    mat=[]
    for g in range(3):
        ta=(-a[g].conj() if bar_first else a[g])
        tb=(-b[g].conj() if bar_second else b[g])
        cols=[]
        for col in range(n):
            x=np.zeros(n,complex);x[col]=1;x=x.reshape(dims)
            y=(np.einsum('ki,kjab->ijab',gs[g],x,optimize=True)
               +np.einsum('kj,ikab->ijab',gs[g],x,optimize=True)
               +np.einsum('ca,ijcb->ijab',ta,x,optimize=True)
               +np.einsum('cb,ijac->ijab',tb,x,optimize=True))
            cols.append(y.ravel())
        mat.append(np.stack(cols,axis=1))
    cols=[]
    for col in range(n):
        x=np.zeros(n,complex);x[col]=1;x=x.reshape(dims)
        cols.append((x-x.swapaxes(0,1)).ravel())
    mat.append(np.stack(cols,axis=1))
    A=np.concatenate(mat,axis=0)
    _,ss,vh=np.linalg.svd(A,full_matrices=False)
    rank=int(np.count_nonzero(ss>max(ss[0],1)*1e-10))
    if rank!=n-1:raise RuntimeError(f'Blocked_NonuniqueSyntheticInvariant: dims {dims}, nullity {n-rank}')
    z=vh.conj().T[:,-1].reshape(dims)
    z/=np.linalg.norm(z)
    if np.max(abs(A@z.ravel()))>1e-9:raise RuntimeError('Blocked_SyntheticInvariantWard')
    # Fix a deterministic convention, not a physical source phase.
    pivot=np.argmax(np.abs(z))
    z=z*np.exp(-1j*np.angle(z.ravel()[pivot]))
    return z


def _complex_real_map(d:int,*,conjugated:bool)->np.ndarray:
    """Map one complex species to its 2*d real components."""
    out=np.zeros((d,2*d),complex)
    for j in range(d):
        out[j,2*j]=1/math.sqrt(2)
        out[j,2*j+1]=(-1j if conjugated else 1j)/math.sqrt(2)
    return out


def synthetic_llss(basis:ModelBasis,flavour:np.ndarray)->np.ndarray:
    """Gauge-singlet LL S1† S2 source in real-scalars, all 3x3 flavours."""
    if np.shape(flavour)!=(3,3) or not np.allclose(flavour,np.asarray(flavour).T):
        raise ValueError('Require symmetric three-family seed')
    d1,d2,_=basis.dims
    T=invariant_two_doublets(d1,d2,bar_first=True)
    A=_complex_real_map(d1,conjugated=True)
    B=_complex_real_map(d2,conjugated=False)
    pair=np.einsum('ijab,aA,bB->ijAB',T,A,B,optimize=True)
    c=np.zeros((basis.nf,basis.nf,basis.ns,basis.ns),complex)
    for p in range(3):
        for q in range(3):
            for i in range(2):
                for j in range(2):
                    f,h=basis.leptons[p][i],basis.leptons[q][j]
                    c[f,h,4:4+2*d1,4+2*d1:]+=flavour[p,q]*pair[i,j]
    return symmetrize(c)


def synthetic_portal(basis:ModelBasis,strength:float)->np.ndarray:
    """Real quartic for gauge singlet H H S1 S2† + h.c., 4! convention."""
    d1,d2,_=basis.dims
    T=invariant_two_doublets(d1,d2,bar_second=True)
    H=_complex_real_map(2,conjugated=False)
    A=_complex_real_map(d1,conjugated=False)
    B=_complex_real_map(d2,conjugated=True)
    L=np.einsum('ijab,iI,jJ,aA,bB->IJAB',T,H,H,A,B,optimize=True)
    raw=np.zeros((basis.ns,)*4,complex)
    raw[:4,:4,4:4+2*d1,4+2*d1:]=float(strength)*L
    lam=np.zeros_like(raw)
    for p in permutations(range(4)):
        lam+=raw.transpose(p)+raw.conj().transpose(p)
    if np.max(abs(lam.imag))>1e-10:raise RuntimeError('Blocked_PortalNotReal')
    return lam.real
