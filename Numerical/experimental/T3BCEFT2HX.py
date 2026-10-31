"""B/C(-1) conditional EFT2(H,X) inert-doublet scalar and Wilson evolution.

Conventions: real scalar components z=(r+i s)/sqrt(2); X has Y=-1/2;
eta=epsilon conj(X) has Y=+1/2. All quartics in the potential:
  lambda1/2 (H+ H)^2 + lambda2/2 (eta+ eta)^2
  +lambda3 (H+ H)(eta+ eta)+lambda4 |H+ eta|^2
  +Re[lambda5 (H+ eta)^2].
Scalar tensor V4=lambda_abcd phi_a phi_b phi_c phi_d / 4!.
Only the indicated EFT2 operator basis is covered. This is not finite matching.
"""
from __future__ import annotations
from functools import lru_cache
from math import factorial, pi
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

LOOP=16*pi*pi
NAMES=('lambda1','lambda2','lambda3','lambda4','lambda5_real','lambda5_imag')


def require(condition, message):
    if not condition:
        raise ValueError('B/C EFT2: '+message)


def potential_monomials(v):
    """Six real invariants, lambda5 multiplying (H^dag eta)^2 + h.c. / 2."""
    arr=np.asarray(v)
    require(arr.shape==(8,), 'expected eight real scalar coordinates')
    h=np.array([(v[0]+sp.I*v[1])/sp.sqrt(2),
                (v[2]+sp.I*v[3])/sp.sqrt(2)],dtype=object)
    x=np.array([(v[4]+sp.I*v[5])/sp.sqrt(2),
                (v[6]+sp.I*v[7])/sp.sqrt(2)],dtype=object)
    # η = ε X*, epsilon_01=+1, epsilon_10=-1.
    eta=np.array([sp.conjugate(x[1]),-sp.conjugate(x[0])],dtype=object)
    normh=sum(sp.conjugate(z)*z for z in h)
    normx=sum(sp.conjugate(z)*z for z in eta)
    z=sum(sp.conjugate(h[k])*eta[k] for k in range(2))
    return [normh**2/2,normx**2/2,normh*normx,z*sp.conjugate(z),
            (z*z+sp.conjugate(z*z))/2,sp.I*(z*z-sp.conjugate(z*z))/2]


@lru_cache(maxsize=1)
def template_quartics():
    """Unambiguous 4th-derivative real tensor templates for six invariant couplings."""
    v=sp.symbols('r0:8',real=True)
    basis=potential_monomials(v)
    q=np.zeros((6,8,8,8,8),dtype=float)
    from itertools import permutations
    for k,expr in enumerate(basis):
        poly=sp.Poly(sp.expand(expr),*v)
        for powers,coefficient in poly.terms():
            nonzero=[j for j,m in enumerate(powers) for _ in range(m)]
            require(len(nonzero)==4,'not quartic')
            val=float(sp.N(coefficient*sp.prod(sp.factorial(m) for m in powers)))
            for index in set(permutations(nonzero)):
                q[(k,)+index]=val
    return q


def quartic_tensor(params):
    z=np.asarray(params,dtype=float)
    require(z.shape==(6,) and np.isfinite(z).all(),'quartics must be six finite real numbers')
    return np.einsum('k,kabcd->abcd',z,template_quartics(),optimize=True)


def fit_quartics(q, tolerance=1.e-7):
    """Fit ALL 8^4 tensor components, not a small selected invariant subset."""
    Q=np.asarray(q,dtype=float)
    require(Q.shape==(8,)*4 and np.isfinite(Q).all(),'invalid projected quartic tensor')
    B=template_quartics().reshape(6,-1).T
    fit,_,rank,singular=np.linalg.lstsq(B,Q.ravel(),rcond=None)
    require(rank==6 and singular[-1]>1e-8,'inert basis rank deficient')
    residual=float(np.linalg.norm((B@fit)-Q.ravel())/max(np.linalg.norm(Q),1e-28))
    require(residual<tolerance,'quartic tensor is not the inert-doublet basis: '+str(residual))
    return fit,residual


def rotate_potential(q_old,T):
    q=np.asarray(q_old,dtype=float)
    T=np.asarray(T,dtype=float)
    require(q.shape==(12,)*4 and T.shape==(12,12),'wrong EFT1 quartic/T shape')
    full=np.einsum('abcd,ap,bq,cr,ds->pqrs',q,T,T,T,T,optimize=True)
    return np.asarray(full[:8,:8,:8,:8])


def pack_sm(sm):
    """gY,g2,g3,ye,yu,yd in one real vector (3+54)."""
    parts=[float(sm.gY),float(sm.g2),float(sm.g3)]
    for n in ('ye','yu','yd'):
        arr=np.asarray(getattr(sm,n),dtype=complex)
        require(arr.shape==(3,3) and np.isfinite(arr).all(), 'bad '+n)
        parts.extend(arr.real.ravel());parts.extend(arr.imag.ravel())
    return np.array(parts,dtype=float)


def unpack_sm(v):
    v=np.asarray(v,dtype=float)
    require(v.shape==(57,) and np.isfinite(v).all(),'SM numerical state size changed')
    data={};pos=3
    for n in ('ye','yu','yd'):
        data[n]=(v[pos:pos+9]+1j*v[pos+9:pos+18]).reshape(3,3)
        pos+=18
    return dict(gY=v[0],g2=v[1],g3=v[2],**data)


def pack_running(sm,quartics,mx2,mh2):
    q=np.asarray(quartics,dtype=float)
    require(q.shape==(6,),'quartic shape')
    y=np.r_[pack_sm(sm),q,float(mx2),float(mh2)]
    require(y.shape==(65,) and np.isfinite(y).all(),'EFT2 running state invalid')
    return y


def unpack_running(y):
    y=np.asarray(y,dtype=float)
    require(y.shape==(65,) and np.isfinite(y).all(),'EFT2 state must have 65 real entries')
    d=unpack_sm(y[:57]);d.update(quartics=y[57:63],mx2=float(y[63]),mh2=float(y[64]))
    return d


def beta_running(t,y):
    """Full one-loop Z2 inert-doublet renormalisable beta, MSbar convention.

    SM Yukawa/gauge part reuses the repository's SM beta. Adds inert scalar gauge
    shifts +1/6 to bY,b2 and adds inert quartic pieces to beta_lambda1.
    H quadratic initial value is a tree-level proxy because the project lacks
    its UV/EFT1 mass parameter; this is explicitly uncertified.
    """
    from RGE.running.weinberg.WeinbergRGE import sm_weinberg_beta
    d=unpack_running(y)
    l1,l2,l3,l4,r5,i5=d['quartics']
    z5=complex(r5,i5);z52=(z5.conjugate()*z5).real
    gy,g2,g3=d['gY'],d['g2'],d['g3']
    ye,yu,yd=(d[n] for n in ('ye','yu','yd'))
    bgY,bg2,bg3,b1,bye,byu,byd,_=sm_weinberg_beta(
        gY=gy,g2_value=g2,g3=g3,lambdaH_value=l1,ye=ye,yu=yu,yd=yd,
        K=np.zeros((3,3),dtype=complex))
    bgY+=gy**3/6;bg2+=g2**3/6
    T=float(np.vdot(ye,ye).real+3*np.vdot(yu,yu).real+3*np.vdot(yd,yd).real)
    G=3*g2*g2+gy*gy
    gauge=3*(3*g2**4+gy**4+2*g2*g2*gy*gy)/4
    mix=4*l3*l3+4*l3*l4+2*l4*l4+2*z52
    b1+=mix
    b2=12*l2*l2+mix+gauge-3*G*l2
    b3=(2*(l1+l2)*(3*l3+l4)+4*l3*l3+2*l4*l4+2*z52
        +3*(3*g2**4+gy**4-2*g2*g2*gy*gy)/4-3*G*l3+2*T*l3)
    b4=(2*l4*(l1+l2+4*l3+2*l4)+8*z52+3*g2*g2*gy*gy
        -3*G*l4+2*T*l4)
    b5=z5*(2*(l1+l2+4*l3+6*l4)-3*G+2*T)
    # mH2 corresponds to V=mH2 |H|^2+lambda1/2 |H|^4.
    bmx=(6*l2-1.5*G)*d['mx2']+(4*l3+2*l4)*d['mh2']
    bmh=(6*l1-1.5*G+2*T)*d['mh2']+(4*l3+2*l4)*d['mx2']
    sm_part=np.r_[bgY,bg2,bg3,bye.real.ravel(),bye.imag.ravel(),
                  byu.real.ravel(),byu.imag.ravel(),byd.real.ravel(),byd.imag.ravel()]
    out=np.r_[sm_part,b1,b2,b3,b4,b5.real,b5.imag,bmx,bmh]/LOOP
    require(out.shape==(65,) and np.isfinite(out).all(),'EFT2 beta nonfinite')
    return out


def solve_threshold(mu_y,y0,*,lower_fraction=0.55):
    """Integrate EFT2 below Y and find unique self-consistent X crossing."""
    require(mu_y>0 and 0.2<lower_fraction<0.85,'bad threshold search interval')
    initial=unpack_running(y0)
    require(initial['mx2']>0 and initial['mx2']<mu_y**2,
            'X not lighter than Y at threshold')
    t0=np.log(mu_y);tl=np.log(mu_y*lower_fraction)
    sol=solve_ivp(beta_running,(t0,tl),y0,method='DOP853',dense_output=True,
                  max_step=0.045,rtol=2.e-9,atol=1.e-10)
    require(sol.success and np.isfinite(sol.y).all(),'EFT2 scalar RGE failed: '+sol.message)
    f=lambda t: float(sol.sol(t)[63]-np.exp(2*t))
    require(f(t0)<0 and f(tl)>0,'X crossing is outside EFT2 interval')
    from scipy.optimize import brentq
    tx=float(brentq(f,tl,t0,xtol=1e-13))
    yx=np.asarray(sol.sol(tx),dtype=float)
    mx=float(np.sqrt(yx[63]))
    require(abs(mx/np.exp(tx)-1)<1e-10,'running X crossing failed')
    return float(np.exp(tx)),yx,sol


def rotate_static(static,T):
    """Remove heavy Y from gauge generators and SM Yukawa tensors."""
    require(static['ns']==12 and T.shape==(12,12),'incompatible static EFT1 scalar basis')
    result={k:static[k] for k in ('indices','nf')}
    result['ns']=8
    result['Ybasis']=np.einsum('xija,ap->xijp',static['Ybasis'],T[:,:8],optimize=True)
    result['gauge']=[]
    for name,ts,tf,_,cf in static['gauge']:
        t8=[T[:,:8].T @ np.asarray(g) @ T[:,:8] for g in ts]
        # Casimir may be computed from sum generators squared in the same
        # convention as the project; alternatively rotate the certified full
        # source Casimir (preserves all generator-convention details).
        full_casimir=next(q[3] for q in static['gauge'] if q[0]==name)
        c8=T[:,:8].T @ full_casimir @ T[:,:8]
        result['gauge'].append((name,t8,tf,c8,cf))
        for old_g,new_g in zip(ts,t8):
            block=(T.T @ np.asarray(old_g) @ T)[:8,8:]
            require(np.linalg.norm(block)<1e-10,'X/Y gauge generator does not close')
    # Y absent: the only nonzero Yukawa components must have H scalar legs.
    require(np.max(abs(result['Ybasis'][...,4:]))<1e-11,
            'SM Yukawa basis illegally couples to inert X')
    return result


def kernels_eft2(state,static,quartic):
    """Unchanged seven-term project master input for an 8-real-scalar EFT."""
    n=3;ns=8;nf=static['nf'];ids=static['indices']
    require(quartic.shape==(8,)*4,'invalid H/X quartic')
    Yvals=[float(np.linalg.norm(state[x],'fro')) for x in ('ye','yu','yd')]
    Y=np.einsum('x,xija->ija',np.asarray(Yvals),static['Ybasis'],optimize=True)
    y=Y[np.ix_(ids,ids,np.arange(ns))]
    yi=Y[ids,:,:];yj=Y[:,ids,:]
    gs=np.zeros((ns,ns),complex);gf=np.zeros((n,n),complex)
    sg=np.zeros((ns,ns,ns,ns),complex);mixed=np.zeros((n,n,ns,ns),complex)
    for name,ts,tf,cs,cf in static['gauge']:
        g=state['gY'] if name=='g1' else state['g2'];g2=g*g
        gs+=-4*g2*cs;gf+=-3*g2*cf[np.ix_(ids,ids)]
        for st,ft in zip(ts,tf):
            sg+=2*g2*np.einsum('ac,bd->abcd',st,st,optimize=True)
            mixed+=4*g2*np.einsum('lj,bc->jlbc',ft[np.ix_(ids,ids)],st,optimize=True)
    gs+=0.5*(np.einsum('jia,jib->ab',Y,Y.conj(),optimize=True)
              +np.einsum('jib,jia->ab',Y,Y.conj(),optimize=True))
    gf+=0.5*np.einsum('ika,jka->ij',yi,yi.conj(),optimize=True)
    mix=2*np.einsum('jkc,klb->jlbc',yi,yj.conj(),optimize=True)
    mix+=np.einsum('jkb,klc->jlbc',yi,yj.conj(),optimize=True)
    return dict(y=y,quartic=quartic,scalarg=sg,mixed=mix+mixed,
                gammas=gs,gammaf=gf)


def rhs_master_eft2(C,k):
    """Exact same seven tensor contractions as original 12-scalar master."""
    c=np.asarray(C,dtype=complex)
    require(c.shape==(3,3,8,8,6),'EFT2 Wilson tensor shape changed')
    y=k['y']
    term=np.einsum('ijd,kld,klabf->ijabf',y,y.conj(),c,optimize=True)
    term+=np.einsum('abcd,ijcdf->ijabf',k['quartic']-k['scalarg'],c,optimize=True)
    mixed=np.einsum('jlbc,ilacf->ijabf',k['mixed'],c,optimize=True)
    term+=mixed+mixed.swapaxes(0,1)+mixed.swapaxes(2,3)+mixed.swapaxes(0,1).swapaxes(2,3)
    crossed=np.einsum('jla,ikc,lkbcf->ijabf',y,y,c,optimize=True)
    term+=crossed+crossed.swapaxes(0,1).swapaxes(2,3)
    cc=c.conj()
    term+=2*np.einsum('ild,jkd,klabf->ijabf',y,y,cc,optimize=True)
    term+=np.einsum('ijd,kld,klabf->ijabf',y,y,cc,optimize=True)
    gs=np.einsum('bd,ijadf->ijabf',k['gammas'],c,optimize=True)
    term+=gs+gs.swapaxes(2,3)
    gf=np.einsum('jk,ikabf->ijabf',k['gammaf'],c,optimize=True)
    term+=gf+gf.swapaxes(0,1)
    require(np.isfinite(term).all(),'EFT2 beta nonfinite')
    return term


def interpolate_kernels(times,kernels,t):
    ts=np.asarray(times,dtype=float)
    require(np.all(np.diff(ts)<0),'EFT2 interpolation grid not descending')
    require(ts[-1]-1e-10<=t<=ts[0]+1e-10,'EFT2 interpolation out of bounds')
    idx=np.searchsorted(-ts,-t,side='right')
    if idx==0:return kernels[0]
    if idx>=len(ts):return kernels[-1]
    a,b=idx-1,idx;w=(t-ts[a])/(ts[b]-ts[a])
    return {key:(1-w)*kernels[a][key]+w*kernels[b][key] for key in kernels[0]}


def evolve_wilson(C_y,mu_y,mu_x,scalar_solution,static,*,samples=33,rtol=5e-7):
    require(samples>=9 and samples%2==1,'odd >=9 Wilson samples required')
    require(mu_y>mu_x>0,'invalid ordered thresholds')
    ts=np.linspace(np.log(mu_y),np.log(mu_x),samples)
    kernels=[]
    for t in ts:
        st=unpack_running(scalar_solution.sol(float(t)))
        kernels.append(kernels_eft2(st,static,quartic_tensor(st['quartics'])))
    original=np.asarray(C_y,dtype=complex)
    shape=original.shape
    require(shape==(3,3,8,8,6),'bad initial Wilson')
    def rhs(t,raw):
        a=raw.reshape((2,)+shape);c=a[0]+1j*a[1]
        dc=rhs_master_eft2(c,interpolate_kernels(ts,kernels,t))/LOOP
        return np.stack((dc.real,dc.imag)).ravel()
    initial=np.stack([original.real,original.imag]).ravel()
    sol=solve_ivp(rhs,(float(ts[0]),float(ts[-1])),initial,
                  t_eval=ts,method='RK45',max_step=0.045,rtol=rtol,atol=1e-24)
    require(sol.success and np.isfinite(sol.y).all(),'EFT2 Wilson solve failed: '+sol.message)
    arr=sol.y[:,-1].reshape((2,)+shape)
    return arr[0]+1j*arr[1],sol
