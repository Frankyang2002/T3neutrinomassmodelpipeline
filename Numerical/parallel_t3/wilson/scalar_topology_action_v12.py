"""Workstream A v12: all16 lambdaT3 source-CG topology scalar-action check.

Construct numerical *declared-frame* symmetric SU2 doublet-to-triplet CG,
including source's three-triplet epsilon contraction, and apply source's
complex lambdaT3 H H S1 Bar[S2] + H.c. and RGBeta 6*Matrix/4! factor.

This is NOT a native RGBeta numerical component Lam, nor a Matchete portal
matching. B/C doubled-doublet has a pseudoreal Bar[S2] index whose native
anti-fundamental intertwiner is not measured here: no physical Ward claim.
"""
from __future__ import annotations

from hashlib import sha1,sha256
import itertools
from pathlib import Path

import numpy as np

from .native_quartic_projectors_v10 import EXPECTED_MODELS,DIMENSIONS
from .scalar_nonportal_rank_v11 import su2_generators,relative

REVIEWED_T3_RGBETA_MODEL_GIT_BLOB = 'bc9500803f2324ad6504e71d769cb75fbbde9b07'
FORMAT='t3-parallel-a-lambdaT3-topology-constructed-action-v12'
STATUS=('All16LambdaT3DeclaredSymmetricCGComplexActionsVerified_'
        'NativeRGBetaComplexLamAndAntifundamentalFrameUnverifiedConditional')
TOL=2e-10
PERMS=tuple(itertools.permutations(range(4)))

class V12Blocked(ValueError):
    pass


def require(ok,msg):
    if not ok:raise V12Blocked('Workstream A v12: '+msg)


def doublet_triplet_cg()->np.ndarray:
    """Normalized symmetric H⊗H CG in a declared REAL Cartesian triplet frame."""
    s=1./np.sqrt(2.)
    weight=np.zeros((3,2,2),complex)
    weight[0,0,0]=weight[2,1,1]=1.
    weight[1,0,1]=weight[1,1,0]=s
    # Spin-one spherical components -> real Cartesian basis. This gives
    # J^a_{bc}=-i epsilon_{abc}, with real orthogonal SU2 rotations.
    r=np.array([[-s,0,s],[-1j*s,0,-1j*s],[0,1,0]],complex)
    c=np.einsum('am,mij->aij',r,weight,optimize=True)
    require(relative(c.reshape(3,4)@c.reshape(3,4).conj().T,np.eye(3))<TOL,
            'Sym2 CG isometry failed')
    require(relative(c,c.transpose(0,2,1))<TOL,'Sym2 CG not symmetric')
    return c


def triple_invariant(c:np.ndarray)->np.ndarray:
    """Source T3RGBetaS2TripleInvariant[a,b,c] using 3 delS2 and 3 eps.

    Source: epsilon[i,k] epsilon[j,m] epsilon[l,n] times
    delS2[a,i,j] delS2[b,k,l] delS2[c,m,n].
    """
    e=np.array([[0.,1.],[-1.,0.]])
    result=np.einsum('aij,bkl,cmn,ik,jm,ln->abc',c,c,c,e,e,e,optimize=True)
    eps=np.zeros((3,3,3))
    for p in itertools.permutations(range(3)):
        eps[p]=(-1)**sum(p[i]>p[j] for i in range(3) for j in range(i+1,3))
    require(relative(result,-1j/np.sqrt(2)*eps)<TOL,
            'source three-triplet epsilon contraction/phase differs')
    return result


def topology_tensor(d1:int,d2:int)->np.ndarray:
    """Source T3RGBetaMixInvariant[d1,d2] in a declared triplet frame."""
    c=doublet_triplet_cg()
    if (d1,d2)==(1,3):
        return c.transpose(1,2,0)[:,:,None,:].copy()
    if (d1,d2)==(3,1):
        return c.transpose(1,2,0)[:,:,:,None].copy()
    if (d1,d2)==(2,2):
        return np.einsum('aij,akl->ijkl',c,c,optimize=True)
    if (d1,d2)==(3,3):
        triple=triple_invariant(c)
        return np.einsum('aij,abc->ijbc',c,triple,optimize=True)
    raise V12Blocked('unsupported source topology scalar dimensions')


def _projections(d1:int,d2:int):
    dims={'H':2,'S1':d1,'S2':d2}
    offsets={};count=0
    for k in ('H','S1','S2'):
        offsets[k]=count;count+=2*dims[k]
    specs=['H','H','S1','Bar[S2]']
    project=[]
    for s in specs:
        key=s[4:-1] if s.startswith('Bar[') else s
        p=np.zeros((dims[key],count),complex)
        for j in range(dims[key]):
            p[j,offsets[key]+2*j]=1.
            p[j,offsets[key]+2*j+1]=(-1j if s.startswith('Bar[') else 1j)
        project.append(p)
    return project,count,offsets


def _complex_fields(x,d1:int,d2:int):
    proj,n,off=_projections(d1,d2)
    require(x.shape==(n,),'topology scalar real-vector dimension mismatch')
    return {k:(x[o:o+2*d:2]+1j*x[o+1:o+2*d:2])/np.sqrt(2.)
            for k,o,d in (('H',off['H'],2),('S1',off['S1'],d1),('S2',off['S2'],d2))}


def topology_source_complex_action(t:np.ndarray,coup:complex,x:np.ndarray)->float:
    d1,d2=t.shape[2:]
    fields=_complex_fields(x,d1,d2)
    v=np.einsum('ijkl,i,j,k,l->',t,fields['H'],fields['H'],
                fields['S1'],fields['S2'].conj(),optimize=True)
    return float((coup*v+np.conj(coup*v)).real)


def topology_constructed_real_lam(t:np.ndarray,coup:complex)->np.ndarray:
    d1,d2=t.shape[2:]
    p,n,off=_projections(d1,d2)
    raw=np.einsum('ijkl,ia,jb,kc,ld->abcd',t,*p,optimize=True)
    sym=sum(np.transpose(raw,perm) for perm in PERMS)/24
    lam=6*(coup*sym+(coup*sym).conj())
    require(np.linalg.norm(lam.imag)<TOL,'non-real constructed topology Lam')
    return lam.real


def _rotation(d:int,axis:np.ndarray,angle:float)->np.ndarray:
    gs=su2_generators(d)
    mat=np.einsum('a,aij->ij',axis,gs)*angle
    eig,U=np.linalg.eigh(mat)
    return (U*np.exp(1j*eig))@U.conj().T


def algebraic_cg_checks(rng:np.random.Generator)->dict:
    cg=doublet_triplet_cg();triple=triple_invariant(cg)
    errors={'SU2_equivariance':0.,'triple_epsilon_invariance':0.,
            'topology_tensor_virtual_all_fund_Ward':0.}
    for _ in range(16):
        axis=rng.normal(size=3);axis/=np.linalg.norm(axis)
        angle=float(rng.uniform(-1.9,1.9))
        u2=_rotation(2,axis,angle);u3=_rotation(3,axis,angle)
        expected=np.einsum('ab,bij->aij',u3,cg,optimize=True)
        actual=np.einsum('aij,ik,jl->akl',cg,u2,u2,optimize=True)
        errors['SU2_equivariance']=max(errors['SU2_equivariance'],relative(actual,expected))
        trip_rot=np.einsum('abc,ai,bj,ck->ijk',triple,u3,u3,u3,optimize=True)
        errors['triple_epsilon_invariance']=max(errors['triple_epsilon_invariance'],
                                              relative(trip_rot,triple))
        # Four *formal* fundamental slots, not physical Bar[S2] placement.
        t=topology_tensor(2,2)
        rotated=np.einsum('ijkl,ia,jb,kc,ld->abcd',t,u2,u2,u2,u2,optimize=True)
        errors['topology_tensor_virtual_all_fund_Ward']=max(
            errors['topology_tensor_virtual_all_fund_Ward'],relative(rotated,t))
    require(max(errors.values())<TOL,'Sym2 CG covariance or epsilon triple failed')
    return errors


def audit(*,seed:int=3771,samples:int=5,root:Path|None=None)->dict:
    require(type(samples) is int and 1<=samples<=20,'samples outside 1..20')
    src={}
    if root is not None:
        path=root/'RGE/running/rgbeta/T3RGBetaModel.wl'
        require(path.is_file(),'RGBeta topology source missing')
        raw=path.read_bytes().replace(b'\r\n',b'\n')
        git_blob=sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()
        require(git_blob==REVIEWED_T3_RGBETA_MODEL_GIT_BLOB,
                'RGBeta topology source differs from reviewed Git blob; re-review v12 CG')
        content=path.read_text(encoding='utf-8-sig')
        for s in ('T3RGBetaMixInvariant[1, 3]', 'T3RGBetaMixInvariant[3, 1]',
                  'T3RGBetaMixInvariant[2, 2]', 'T3RGBetaMixInvariant[3, 3]',
                  'T3RGBetaS2TripleInvariant', 'SelfConjugate -> False'):
            require(s in content,'RGBeta source topology definition changed '+s)
        src={'T3RGBetaModel.wl_sha256':sha256(path.read_bytes()).hexdigest(),
             'reviewed_T3RGBetaModel_git_blob':git_blob}
    rng=np.random.default_rng(seed);group=algebraic_cg_checks(rng)
    models={};worst=0.;total=0
    for tag in sorted(EXPECTED_MODELS):
        d1,d2,df=DIMENSIONS[tag[0]]
        t=topology_tensor(d1,d2)
        require(t.shape==(2,2,d1,d2) and np.linalg.norm(t)>1e-8,
                tag+': topology source tensor shape/zero')
        errs=[]
        for _ in range(samples):
            x=rng.normal(size=4+2*d1+2*d2)
            coup=complex(rng.normal(),rng.normal())
            a=topology_source_complex_action(t,coup,x)
            l=topology_constructed_real_lam(t,coup)
            b=float(np.einsum('ijkl,i,j,k,l->',l,x,x,x,x,optimize=True)/24)
            res=relative(a,b);errs.append(res);worst=max(worst,res)
            require(res<TOL,tag+': topology complex action not reconstructed')
        models[tag]={'dimensions':[d1,d2,df],
          'quartic_symbol':'lambdaT3','native_self_conjugate':False,
          'fields':['H','H','S1','Bar[S2]'],
          'source_CG_tensor_norm':float(np.linalg.norm(t)),
          'max_constructed_complex_action_residual':max(errs),
          'native_RGBeta_Lam_executed':False,
          'B_C_bar_antifundamental_intertwiner_certified':False}
        total+=1
    require(total==16 and worst<TOL,'incomplete topology quartic action census')
    return {'format':FORMAT,'status':STATUS,'models_checked':16,
        'lambdaT3_complex_quartic_slots_constructed':total,
        'SU2_declared_frame_CG_checks':group,
        'max_constructed_complex_action_residual':worst,
        'models':models,'source_sha256':src,
        'certification':{
         'declared_Sym2_CG_and_triple_epsilon_identity_verified':True,
         'lambdaT3_complex_plus_Hc_constructed_action_verified':True,
         'installed_RGBeta_Lam_and_complex_projector_certified':False,
         'physical_B_C_pseudoreal_conjugate_intertwiner_certified':False,
         'Matchete_native_lambdaT3_CG_component_dictionary_certified':False,
         'complete_scalar_potential_normalization_certified':False,
         'physical_Wilson_evolution_authorized':False,
         'EFT_order':'UV -> integrate out ALL fermions -> integrate out ALL scalars -> SMEFT'
         }}
