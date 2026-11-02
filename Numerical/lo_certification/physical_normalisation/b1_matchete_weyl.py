"""B1 restricted T3 Matchete -> conditional Weyl rule algebra.

No claim of installed-Matchete GammaCC, RGBeta epsilon axes, or physical
cross-package normalisation. This is new isolated code, not a production adapter.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from typing import Mapping
import numpy as np
import sympy as sp

CLASSES = {
    'A': (1, 3, 2, (-4, -2, 0)),
    'B': (2, 2, 1, (-3, -1, 1)),
    'C': (2, 2, 3, (-3, -1, 1)),
    'D': (3, 1, 2, (-2, 0, 2)),
    'E': (3, 3, 2, (-4, -2, 0, 2)),
}
FLAGS = {'T3Y1CG': (True, True, False),
         'T3Y2CG': (True, False, True)}

class B1Error(ValueError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise B1Error(message)


def model_specs():
    for label, (d1, d2, df, alphas) in CLASSES.items():
        for alpha in alphas:
            yield f'{label}_{alpha:+d}', (d1, d2, df, alpha)


def is_majorana(df: int, alpha: int) -> bool:
    return (df % 2 == 1) and alpha == -1


def charge(d: int, hypercharge: float, index: int) -> float:
    return (d-1)/2 - index + hypercharge


def complex_array(raw):
    if isinstance(raw, Mapping):
        require(set(raw) == {'real', 'imag'}, 'Unsupported complex encoding')
        z = np.asarray(raw['real'], dtype=float) + 1j*np.asarray(raw['imag'], dtype=float)
    else:
        z = np.asarray(raw, dtype=complex)
    require(np.isfinite(z).all(), 'Nonfinite input tensor')
    return z


def encode(z):
    a = np.asarray(z, complex)
    return {'real': a.real.tolist(), 'imag': a.imag.tolist()}


# Exterior/Grassmann algebra, independent of the production numerical engine.
_N = ('l0','l1','c0','c1','x0','x1','lb0','lb1','cb0','cb1','xb0','xb1')
_P = {name: i for i, name in enumerate(_N)}
_D = dict(zip(_N[:6],_N[6:]))
_D.update(dict(zip(_N[6:],_N[:6])))

def var(label):
    require(label in _P, 'Unknown Grassmann symbol')
    return {(label,): sp.Integer(1)}


def add(*polys):
    out = {}
    for p in polys:
        for m, c in p.items():
            out[m] = sp.expand(out.get(m, 0)+c)
    return {m:c for m,c in out.items() if sp.simplify(c)!=0}


def scale(p, factor):
    return {m:sp.expand(c*factor) for m,c in p.items() if sp.simplify(c*factor)!=0}


def wedge(a,b):
    out = {}
    for ma,ca in a.items():
        for mb,cb in b.items():
            seq = ma+mb
            if len(set(seq))<len(seq):
                continue
            sign = (-1)**sum(_P[seq[i]]>_P[seq[j]] for i in range(len(seq)) for j in range(i+1,len(seq)))
            key = tuple(sorted(seq,key=_P.__getitem__))
            out[key] = sp.expand(out.get(key,0)+sign*ca*cb)
    return {m:c for m,c in out.items() if sp.simplify(c)!=0}


def dagger(poly):
    pieces=[]
    for mon, coefficient in poly.items():
        q={():sp.conjugate(coefficient)}
        for label in reversed(mon):
            q=wedge(q,var(_D[label]))
        pieces.append(q)
    return add(*pieces)


def dot(left,right):
    return add(*[wedge(a,b) for a,b in zip(left,right)])


def mv(mat,vec):
    return [add(*[scale(vec[j],mat[i,j]) for j in range(mat.cols)]) for i in range(mat.rows)]


def bilinear(a,b):
    return add(wedge(a[0],b[1]),scale(wedge(a[1],b[0]),-1))


def _gamma():
    eye=sp.eye(2);zero=sp.zeros(2)
    sigmas=(sp.Matrix([[0,1],[1,0]]),sp.Matrix([[0,-sp.I],[sp.I,0]]),sp.diag(1,-1))
    def block(a,b,c,d):return sp.BlockMatrix([[a,b],[c,d]]).as_explicit()
    g=(block(zero,eye,eye,zero),)+tuple(block(zero,s,-s,zero) for s in sigmas)
    C=sp.simplify(sp.I*g[2]*g[0])
    eps=sp.Matrix([[0,1],[-1,0]])
    require(C==block(eps,zero,zero,-eps), 'Unexpected C')
    for a in range(4):
        for b in range(4):
            require(g[a]*g[b]+g[b]*g[a]==(2 if a==b==0 else -2 if a==b else 0)*sp.eye(4),'Clifford')
    require(C.T==-C and C*C==-sp.eye(4),'C symmetry')
    require(all(C*gmu.T*C.inv()==-gmu for gmu in g),'C intertwiner')
    return g,C


def four_spinor_witness():
    """Exact chosen-basis proof, NOT native Matchete measurement."""
    g,C=_gamma();zero={}
    L=[var('l0'),var('l1'),zero,zero]
    bL=mv(g[0],[dagger(v) for v in L])
    # F=(chi_alpha,-epsilon bar(xi)), barF=(-xi1,xi0,bar(chi)0,bar(chi)1)
    F=[var('c0'),var('c1'),scale(var('xb1'),-1),var('xb0')]
    bF=mv(g[0],[dagger(v) for v in F])
    expected_bF=[scale(var('x1'),-1),var('x0'),var('cb0'),var('cb1')]
    require(bF==expected_bF,'Dirac adjoint')
    CC=mv(C,bF)
    cc_expected=[var('x0'),var('x1'),scale(var('cb1'),-1),var('cb0')]
    require(CC==cc_expected,'CConj field')
    chi=[var('c0'),var('c1')];xi=[var('x0'),var('x1')];l=L[:2]
    w1=bilinear(l,chi);w2=bilinear(l,xi);md=bilinear(chi,xi)
    require(dot(bL,CC)==dagger(w1),'y1 identity')
    require(dot(bL,F)==dagger(w2),'y2 identity')
    require(dot(bF,F)==add(md,dagger(md)),'Dirac mass identity')
    require(dot(bL,mv(-C,bF))!=dagger(w1),'C negative control')
    # Algebraic Majorana substitution xi=chi, with no factorial convention hidden.
    subst={'x0':'c0','x1':'c1','xb0':'cb0','xb1':'cb1'}
    def replace(p):
        out={}
        for mon,coeff in p.items():
            q={():coeff}
            for label in mon:q=wedge(q,var(subst.get(label,label)))
            out=add(out,q)
        return out
    majorana_mass=replace(dot(bF,F))
    majorana_ww=bilinear(chi,chi)
    require(majorana_mass==add(majorana_ww,dagger(majorana_ww)), 'Majorana bilinear')
    require(majorana_ww==scale(wedge(chi[0],chi[1]),2), 'Majorana spinor factor')
    require(replace(dot(bL,CC))==replace(dot(bL,F)), 'Majorana chiral identity')
    return {'chosen_C':'+i gamma^2 gamma^0','epsilon_12':1,
            'original_y1':'Bar[L] CConj[F] = (l epsilon chi)^dagger',
            'original_y2':'Bar[L] F = (l epsilon xi)^dagger',
            'dirac_mass':'Bar[F] F = chi epsilon xi + h.c.',
            'majorana_mass':'Bar[F] F = chi epsilon chi + h.c.',
            'majorana_free_mass_multiplier':'1/2',
            'canonical_YZ':'Z1=-conjugate(y1 * native_CG1); Y2=-conjugate(y2 * native_CG2)',
            'exact_chosen_basis_checked':True,
            'absolute_native_GammaCC_phase_measured':False}


def triplet_metric(U, t3=None):
    U=np.asarray(U,complex)
    require(U.shape==(3,3) and np.linalg.norm(U.conj().T@U-np.eye(3))<1.e-11,'U not unitary triplet')
    G=U.T@U
    require(np.linalg.norm(G-G.T)<1.e-11 and np.linalg.norm(G.conj().T@G-np.eye(3))<1.e-11,'Invalid triplet pairing')
    # Cartesian spin-1 generator T3=-i epsilon_{3bc}, weights ordered +1,0,-1.
    weights=np.array([1.,0.,-1.])
    if t3 is not None:
        require(np.linalg.norm(np.asarray(t3,complex)@U-U@np.diag(weights))<1.e-9,'Incorrect triplet weight ordering')
    for i,j in product(range(3),repeat=2):
        if abs(weights[i]+weights[j])>1.e-10:
            require(abs(G[i,j])<1.e-9,'Triplet charge pairing violated')
    return G


def spurion_charge(majorana=False):
    # Charges under old field rephasings (F,Fc,S1,S2); external neutrino held fixed.
    charges = ({'y1':(1,0,1,0),'y2':(1,0,0,-1),'M':(2,0,0,0),'lam':(0,0,1,-1)}
               if majorana else {'y1':(0,1,1,0),'y2':(1,0,0,-1),'M':(1,1,0,0),'lam':(0,0,1,-1)})
    ans=tuple(charges['y1'][i]+charges['y2'][i]-charges['M'][i]-charges['lam'][i] for i in range(4))
    require(ans==(0,0,0,0),'Incorrect invariant source spurion')
    wrong=tuple(charges['y1'][i]+charges['y2'][i]-charges['M'][i]+charges['lam'][i] for i in range(4))
    require(wrong!=(0,0,0,0),'Incorrect lambda conjugation not detected')
    return {'invariant':'y1*y2*conjugate(MF)*conjugate(lambdaT3)','charge':list(ans),'wrong_lambda_charge':list(wrong)}


def rotate_tensor(vertex,flags,Us):
    dims=tuple(int(d) for d in vertex['dimensions'])
    require(len(dims)==len(flags),'Bad CG tensor rank')
    old=np.zeros(dims,complex)
    for entry in vertex['components']:
        ix=tuple(int(i)-1 for i in entry['indices'])
        require(len(ix)==len(dims) and all(0<=i<d for i,d in zip(ix,dims)),'Bad CG index')
        old[ix]+=complex(sp.N(sp.sympify(entry['cg'],locals={'I':sp.I}),17))
    for axis,(d,conj) in enumerate(zip(dims,flags)):
        U=np.asarray(Us[d],complex)
        require(U.shape==(d,d),'Wrong U size')
        require(np.linalg.norm(U.conj().T@U-np.eye(d))<1.e-10,'Nonunitary generator rotation')
        transformation=U.conj() if conj else U
        old=np.moveaxis(np.tensordot(old,transformation,axes=([axis],[0])),-1,axis)
    return old


def construct_components(case,benchmark,Us,short,*,eta1=1.+0.j,eta2=1.+0.j,mass_factor=1.):
    """Independent candidate construction at fixed B1 chosen Weyl basis.

    eta1/eta2 are the still-unverified Matchete/GammaCC and B2 SU2
    relative map coefficients; they are NOT fitted or silently set to physical.
    mass_factor=1 is the canonical Majorana 1/2 interpretation.
    """
    allowed=dict(model_specs())
    require(short in allowed, 'Unknown restricted T3 model')
    d1,d2,df,a=allowed[short]
    require(tuple(int(case['model'][n]) for n in ('dS1','dS2','dF','alpha')) == (d1,d2,df,a),'CG case/model mismatch')
    rep=benchmark['representation']
    require(tuple(int(rep[n]) for n in ('d_s1','d_s2','d_f','alpha'))==(d1,d2,df,a),'Benchmark/model mismatch')
    o=benchmark['base_state']['ordinary']
    y1,y2,mf=(complex_array(o[name]) for name in ('y1','y2','MF'))
    require(all(v.shape==(3,3) for v in (y1,y2,mf)),'Expected 3x3 flavour tensors')
    majorana=is_majorana(df,a)
    if majorana:require(np.linalg.norm(mf-mf.T)<1.e-10*max(np.linalg.norm(mf),1.),'Majorana MF not symmetric')
    nf=3*df if majorana else 6*df; ns=d1+d2; nc=3*df
    G=np.asarray(Us[df].T@Us[df],complex)
    if df==3:triplet_metric(Us[df])
    if df==1:require(np.linalg.norm(G-1)<1.e-10,'Bad singlet metric')
    M=np.zeros((nf,nf),complex)
    if majorana:
        M=(mass_factor*np.kron(mf,G)).astype(complex)
    else:
        # psi=(xi,chi), M_{xi_r,chi_t}=mf[r,t], M_{chi_t,xi_r}=mf[r,t]
        M[:nc,nc:]=np.kron(mf,np.eye(df));M[nc:,:nc]=M[:nc,nc:].T
    require(np.linalg.norm(M-M.T)<1.e-10*max(1.,np.linalg.norm(M)),'Weyl M not symmetric')
    cg1=rotate_tensor(case['vertices']['T3Y1CG'],FLAGS['T3Y1CG'],Us)
    cg2=rotate_tensor(case['vertices']['T3Y2CG'],FLAGS['T3Y2CG'],Us)
    require(cg1.shape==(2,df,d1) and cg2.shape==(2,df,d2),'Unexpected native CG dimensions')
    Y=np.zeros((3,nf,ns),complex); Z=np.zeros_like(Y)
    for r,k in product(range(3),range(df)):
        q=charge(df,(a+1)/2,k)
        for s in range(d1):
            factor=np.conj(cg1[0,k,s]);
            if abs(factor)<1.e-13:continue
            require(abs(q-charge(d1,a/2,s))<1.e-10, 'y1 charge mismatch')
            column=r*df+k if majorana else nc+r*df+k
            Z[:,column,s]-=eta1*np.conj(y1[:,r])*factor
        for s in range(d2):
            factor=np.conj(cg2[0,k,s]);
            if abs(factor)<1.e-13:continue
            require(abs(q-charge(d2,(a+2)/2,s))<1.e-10,'y2 charge mismatch')
            if majorana:
                Y[:,r*df:(r+1)*df,d1+s]-=eta2*np.conj(y2[:,r])[:,None]*factor*G[k,:][None,:]
            else:
                Y[:,r*df+k,d1+s]-=eta2*np.conj(y2[:,r])*factor
    return {'M':M,'Y':Y,'Z':Z,'fermion_type':'Majorana' if majorana else 'Dirac',
            'metric':G,'eta1':complex(eta1),'eta2':complex(eta2)}


def relative(a,b):
    return float(np.linalg.norm(a-b)/max(np.linalg.norm(a),np.linalg.norm(b),1.e-100))


def congruence(M,Y,Z,D,S):
    """psi_old=D psi_new, phi_old=S phi_new; D,S unitary diagonal."""
    return (D.T@M@D, np.einsum('ais,ij,st->ajt',Y,D,S),
            np.einsum('ais,ij,st->ajt',Z,D,S.conj()))


def effective_tree_kernel(M,V):
    """External-neutrino pair Green coefficient V M^-1 V^T, no 1/2.

    V has (flavour,heavy) axes. Not a T3 one-loop neutrino mass.
    """
    M=np.asarray(M,complex);V=np.asarray(V,complex)
    return V @ np.linalg.solve(M,V.T)


def b0_finite_independent(mf2,ms2,mu2):
    """Dim-reg MSbar B0(0) finite part, independent analytic expression."""
    require(mf2>0 and ms2>0 and mu2>0,'Positive squared masses required')
    if abs(ms2-mf2)<1e-7*mf2:
        # B0(0;m2,m2) = -ln(m2/mu2)
        return -np.log(mf2/mu2)
    return 1. - (mf2*np.log(mf2/mu2)-ms2*np.log(ms2/mu2))/(mf2-ms2)


def diagonal_toy_self_energy(g,fermion_masses,scalar_masses_sq,mu2):
    """Two-component p=0 chirality-flipping Majorana self-energy.

    g[a,k,t] already includes any 1/sqrt(2) from phi=(R+iI)/sqrt(2).
    No Majorana factor 2 or 1/2 is inserted into the loop amplitude.
    """
    g=np.asarray(g,complex)
    masses=np.asarray(fermion_masses,float)
    scalars=np.asarray(scalar_masses_sq,float)
    require(g.ndim==3 and g.shape[1:]==(len(masses),len(scalars)),'Bad toy dimensions')
    out=np.zeros((g.shape[0],g.shape[0]),complex)
    for k,mf in enumerate(masses):
        for t,ms2 in enumerate(scalars):
            out+=np.outer(g[:,k,t],g[:,k,t])*mf*b0_finite_independent(mf*mf,ms2,mu2)/(16*np.pi**2)
    return out
