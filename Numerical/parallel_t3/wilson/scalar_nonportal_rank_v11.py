"""Workstream A v11: independent SU(2) scalar-invariant rank/action atlas.

Mathematical frames are declared here (Pauli fundamental and *real Cartesian*
spin-one). They are NOT extracted RGBeta numerical axes. Exact RG-closed
quartic shapes are read from T3RGBetaModel.wl, not from physical Matchete CGs.

Coverage:
 * all 16 EFT1 quartic inventories (178 model-coupling registrations);
 * 14 triplet-self adjoint quartic independent directions;
 * 10 S1-S2 adjoint invariant directions;
 * 4 E-class S1-S2 crossed-triplet directions;
 * 12 special B/C(-1) complex quartic real-action checks (6 operators x2).

No numerical native Lam is executed; no physical scalar dictionary is claimed.
"""
from __future__ import annotations

from hashlib import sha1,sha256
import itertools
import json
from pathlib import Path

import numpy as np

from .native_quartic_projectors_v10 import EXPECTED_MODELS,DIMENSIONS,SPECIAL,expected_symbols

REVIEWED_T3_RGBETA_MODEL_GIT_BLOB = 'bc9500803f2324ad6504e71d769cb75fbbde9b07'
FORMAT='t3-parallel-a-nonportal-scalar-invariants-conditional-v11'
STATUS=('NonportalSU2InvariantRankAtlasAnd12SpecialComplexQuarticActionsVerified_'
        'NativeComplexProjectorAndMatcheteNonportalDictionaryUncertifiedConditional')
TOL=2e-10
PERMS=tuple(itertools.permutations(range(4)))

class V11Blocked(ValueError):
    """Source, invariant covariance, action or rank audit blocked."""


def require(ok:bool,msg:str)->None:
    if not ok:raise V11Blocked('Workstream A v11: '+msg)


def relative(a,b)->float:
    a=np.asarray(a);b=np.asarray(b)
    return float(np.linalg.norm(a-b)/max(1.,np.linalg.norm(a),np.linalg.norm(b)))


def su2_generators(d:int)->np.ndarray:
    """Hermitian t^a, trace index T(2)=1/2, T(3)=2, C2(3)=2."""
    if d==2:
        return np.asarray([
          [[0.,.5],[.5,0.]],[[0.,-.5j],[.5j,0.]],
          [[.5,0.],[0.,-.5]]],complex)
    if d==3:
        eps=np.zeros((3,3,3),float)
        for p in itertools.permutations(range(3)):
            parity=sum(p[i]>p[j] for i in range(3) for j in range(i+1,3))
            eps[p]=(-1)**parity
        return -1j*eps
    raise V11Blocked('Unsupported generator dimension '+str(d))


def unitary_rotation(d:int,rng:np.random.Generator)->np.ndarray:
    g=su2_generators(d)
    vec=rng.normal(size=3);vec/=np.linalg.norm(vec)
    angle=float(rng.uniform(-2.6,2.6))
    theta=np.einsum('a,aij->ij',vec,g)*angle
    eigen,vecs=np.linalg.eigh(theta)
    return (vecs*np.exp(1j*eigen))@vecs.conj().T


def source_invariant_values(d1:int,d2:int,s1:np.ndarray,s2:np.ndarray)->dict[str,float]:
    """Real-axis scalar polynomials in a declared SU(2) numerical frame."""
    require(s1.shape==(d1,) and s2.shape==(d2,),'scalar axis dimensions')
    out={'lambda12':float(np.vdot(s1,s1).real*np.vdot(s2,s2).real)}
    if d1>1 and d2>1:
        g1=su2_generators(d1);g2=su2_generators(d2)
        a=np.einsum('i,aij,j->a',s1.conj(),g1,s1)
        b=np.einsum('i,aij,j->a',s2.conj(),g2,s2)
        value=np.dot(a,b)
        require(abs(value.imag)<TOL,'mixed adjoint polynomial not real')
        out['lambda12Adj']=float(value.real)
    if d1==d2==3:
        # delta[i,k]delta[j,l]+delta[i,l]delta[j,k] with fields
        # {Bar[S1], S1, Bar[S2], S2}. A *real Cartesian* S2 basis is
        # necessary for the unstarred delS2 cross-contraction here.
        v=abs(np.dot(s1,s2))**2+abs(np.vdot(s1,s2))**2
        out['lambda12Cross']=float(v)
    return out


def triplet_self_values(s:np.ndarray)->dict[str,float]:
    require(s.shape==(3,),'triplet field wrong dimension')
    j=su2_generators(3)
    norm=np.vdot(s,s)
    adj=np.einsum('i,aij,j->a',s.conj(),j,s)
    out={'self_norm':float((norm**2/2).real),
         'self_adj':float((np.dot(adj,adj)/2).real)}
    require(abs((norm**2/2).imag)<TOL and
            abs((np.dot(adj,adj)/2).imag)<TOL,'complex quartic action')
    # spin-one vector identity: sum_a (s†J_a s)^2
    # = (s†s)^2 - |s^T s|^2 in a REAL Cartesian frame.
    lhs=2*out['self_adj']
    rhs=(norm.real**2 - abs(np.dot(s,s))**2)
    require(abs(lhs-rhs)<TOL*max(1.,abs(lhs),abs(rhs)),
            'triplet-adjoint Fierz relation failed')
    return out


def rank_report(rows:list[list[float]],columns:list[str],label:str)->dict:
    mat=np.asarray(rows,dtype=float)
    require(mat.ndim==2 and mat.shape[1]==len(columns) and
            mat.shape[0]>=4*len(columns) and np.all(np.isfinite(mat)),
            label+': insufficient design matrix')
    colnorm=np.linalg.norm(mat,axis=0)
    require(np.all(colnorm>1.e-12),label+': identically zero source quartic')
    u,sv,vt=np.linalg.svd(mat/colnorm,full_matrices=False)
    r=int(np.sum(sv>1e-8))
    require(r==len(columns) and sv[-1]>1.e-3,
            label+': source scalar invariant rank deficient')
    return {'columns':columns,'verified_real_dimension':r,
            'singular_values_column_normalized':sv.tolist(),
            'minimum_singular_value':float(sv[-1])}


def _epsilon()->np.ndarray:
    return np.asarray([[0.,1.],[-1.,0.]])


def special_invariant_tensor(name:str)->tuple[list[str],np.ndarray]:
    """Six source-exact B/C(-1) GroupInvariant quartic tensor formulas.

    All indices are SU(2) fundamental *formal* slots. The complex conjugate
    action is rebuilt independently below; no Ward/basis matching is claimed.
    """
    e=_epsilon();d=np.eye(2)
    t=np.zeros((2,2,2,2),float)
    if name=='lambdaHHdagS2S2':
        species=['Bar[H]','Bar[H]','S2','S2']
        for i,j,k,l in itertools.product(range(2),repeat=4):
            t[i,j,k,l]=(d[i,k]*d[j,l]+d[i,l]*d[j,k])/2.
    elif name=='lambdaHHdagS1barS1bar':
        species=['Bar[H]','Bar[H]','Bar[S1]','Bar[S1]']
        for i,j,k,l in itertools.product(range(2),repeat=4):
            t[i,j,k,l]=(e[i,l]*e[j,k]+e[i,k]*e[j,l])/2.
    elif name=='lambdaS1bar2S2bar2':
        species=['Bar[S1]','Bar[S1]','Bar[S2]','Bar[S2]']
        for i,j,k,l in itertools.product(range(2),repeat=4):
            t[i,j,k,l]=(e[i,l]*e[j,k]+e[i,k]*e[j,l])/2.
    elif name=='lambdaS1barS2S2bar2':
        species=['Bar[S1]','S2','Bar[S2]','Bar[S2]']
        for i,j,k,l in itertools.product(range(2),repeat=4):
            t[i,j,k,l]=(e[i,l]*d[j,k]+e[i,k]*d[j,l])/2.
    elif name=='lambdaS1S1bar2S2bar':
        species=['S1','Bar[S1]','Bar[S1]','Bar[S2]']
        for i,j,k,l in itertools.product(range(2),repeat=4):
            t[i,j,k,l]=(d[i,k]*e[j,l]+d[i,j]*e[k,l])/2.
    elif name=='lambdaHHdagS1barS2barCross':
        species=['H','Bar[H]','Bar[S1]','Bar[S2]']
        for i,j,k,l in itertools.product(range(2),repeat=4):
            t[i,j,k,l]=d[i,k]*e[j,l]
    else:raise V11Blocked('Unknown native complex quartic '+str(name))
    require(np.linalg.norm(t)>0,'zero complex source tensor '+name)
    return species,t


def _complex_projections(species:list[str])->tuple[list[np.ndarray],int]:
    mapping={};count=0
    for label in species:
        key=label[4:-1] if label.startswith('Bar[') else label
        if key not in mapping:
            mapping[key]=count;count+=4
    ps=[]
    for label in species:
        key=label[4:-1] if label.startswith('Bar[') else label
        p=np.zeros((2,count),complex)
        for i in range(2):
            p[i,mapping[key]+2*i]=1.
            p[i,mapping[key]+2*i+1]= -1j if label.startswith('Bar[') else 1j
        ps.append(p)
    return ps,count


def source_complex_action(t:np.ndarray,species:list[str],c:complex,x:np.ndarray)->float:
    """Native source's complex quartic plus its separate conjugate channel."""
    mapping={};offset=0
    for label in species:
        key=label[4:-1] if label.startswith('Bar[') else label
        if key not in mapping:
            mapping[key]=(x[offset:offset+4:2]+1j*x[offset+1:offset+4:2])/np.sqrt(2)
            offset+=4
    vec=[mapping[label[4:-1]].conj() if label.startswith('Bar[')
         else mapping[label] for label in species]
    z=np.einsum('ijkl,i,j,k,l->',t,*vec,optimize=True)
    return float((c*z+np.conj(c*z)).real)


def constructed_complex_lam(t:np.ndarray,species:list[str],c:complex)->np.ndarray:
    """Constructs a real symmetric tensor from source's 6*Matrix / 4! rule.

    NOT installed native Lam; factor 6 accounts for four complex scalar legs.
    """
    ps,n=_complex_projections(species)
    raw=np.einsum('ijkl,ia,jb,kc,ld->abcd',t,*ps,optimize=True)
    sym=sum(np.transpose(raw,p) for p in PERMS)/24.
    lam=6*(c*sym+(c*sym).conj())
    require(np.linalg.norm(lam.imag)<TOL,'constructed complex quartic Lam not real')
    return lam.real


def complex_action_holdout(name:str,rng:np.random.Generator,samples:int)->dict:
    species,t=special_invariant_tensor(name)
    ps,n=_complex_projections(species)
    worst=0.;nonzero=0
    for i in range(samples):
        x=rng.normal(size=n)
        c=complex(rng.uniform(-1.7,1.7),rng.uniform(-1.7,1.7))
        src=source_complex_action(t,species,c,x)
        lam=constructed_complex_lam(t,species,c)
        cmp=float(np.einsum('ijkl,i,j,k,l->',lam,x,x,x,x,optimize=True)/24.)
        err=relative(src,cmp)
        require(err<TOL,name+': complex action/H.c. or 6*Matrix mismatch')
        worst=max(worst,err)
        nonzero+=int(abs(src)>1.e-8)
    require(nonzero>0,name+': sampled complex source quartic is identically zero')
    return {'fields':species,'raw_tensor_frobenius_norm':float(np.linalg.norm(t)),
            'max_constructed_real_Lam_action_relative_residual':worst,
            'native_RGBeta_Lam_executed':False,'native_Matchete_CG_compared':False}


def _validate_source_contract(root:Path)->dict:
    p=root/'RGE/running/rgbeta/T3RGBetaModel.wl'
    require(p.is_file(),'current RGBeta EFT1 model source missing')
    raw=p.read_bytes()
    # Compare reviewed Git blob after Windows CRLF normalization. Source
    # changes require a fresh review, not a silent reuse of transcribed CGs.
    git_bytes=raw.replace(b'\r\n',b'\n')
    blob_id=sha1(b'blob '+str(len(git_bytes)).encode()+b'\0'+git_bytes).hexdigest()
    require(blob_id==REVIEWED_T3_RGBETA_MODEL_GIT_BLOB,
            'RGBeta model source differs from reviewed Git blob; re-review v11 invariants')
    src=p.read_text(encoding='utf-8-sig')
    for name in (*SPECIAL,'lambdaS1Adj','lambdaS2Adj','lambda12Adj',
                 'lambda12Cross','T3RGBetaAddRGClosedQuartics',
                 'T3RGBetaBuildEFT1','T3RGBetaPortalAdjInvariant'):
        require(name in src,'scalar potential source no longer declares '+name)
    return {'T3RGBetaModel.wl_sha256':sha256(p.read_bytes()).hexdigest(),
            'reviewed_T3RGBetaModel_git_blob':blob_id}


def audit(*,seed:int=21551,samples:int=128,root:Path|None=None)->dict:
    require(type(seed) is int and type(samples) is int and 32<=samples<=300,
            'invalid sampling settings')
    provenance={} if root is None else _validate_source_contract(root)
    rng=np.random.default_rng(seed)
    g2=su2_generators(2);g3=su2_generators(3)
    eps=np.zeros((3,3,3),float)
    for p in itertools.permutations(range(3)):
        inv=sum(p[i]>p[j] for i in range(3) for j in range(i+1,3))
        eps[p]=(-1)**inv
    generator_errors={}
    for d,g,index in ((2,g2,.5),(3,g3,2.)):
        worst=0.
        for a,b in itertools.product(range(3),repeat=2):
            want=1j*sum(eps[a,b,c]*g[c] for c in range(3))
            worst=max(worst,relative(g[a]@g[b]-g[b]@g[a],want))
            want_trace=index if a==b else 0.
            worst=max(worst,abs(np.trace(g[a]@g[b])-want_trace))
        worst=max(worst,relative(sum(m@m for m in g),
            ({2:.75,3:2.}[d])*np.eye(d)))
        require(worst<TOL,'SU2 numerical generators algebra invalid')
        generator_errors[str(d)]=worst
    models={};rank_counts={'triplet_self_sectors':0,'mixed_adjoint_sectors':0,
                           'mixed_crossed_triplet_sectors':0,
                           'special_complex_slots':0}
    worst_action=0.;worst_ward=0.
    for tag in sorted(EXPECTED_MODELS):
        d1,d2,df=DIMENSIONS[tag[0]]
        regs=expected_symbols(tag)
        mode={};
        if d1==3 or d2==3:
            self_table=[]
            for which in (1,2):
                if (d1 if which==1 else d2)!=3:continue
                table=[]
                for _ in range(samples):
                    s=rng.normal(size=3)+1j*rng.normal(size=3)
                    values=triplet_self_values(s)
                    table.append([values['self_norm'],values['self_adj']])
                self_table.append({'which':which,**rank_report(
                    table,['self_norm','self_adj'],tag+f'/S{which} self')})
                rank_counts['triplet_self_sectors']+=1
            mode['triplet_self']=self_table
        cross_table=[]
        for _ in range(samples):
            s1=rng.normal(size=d1)+1j*rng.normal(size=d1)
            s2=rng.normal(size=d2)+1j*rng.normal(size=d2)
            inv=source_invariant_values(d1,d2,s1,s2)
            cross_table.append([inv[k] for k in inv])
            if d1>1 and d2>1 and _<6:
                # Single *same* SU(2) rotation for both scalar reps.
                axis=rng.normal(size=3);axis/=np.linalg.norm(axis)
                angle=float(rng.uniform(-2.,2.))
                u=[]
                for d in (d1,d2):
                    j=su2_generators(d)
                    herm=angle*np.einsum('a,aij->ij',axis,j)
                    eig,v=np.linalg.eigh(herm)
                    u.append((v*np.exp(1j*eig))@v.conj().T)
                rotated=source_invariant_values(d1,d2,u[0]@s1,u[1]@s2)
                worst_ward=max(worst_ward,max(relative(rotated[k],v) for k,v in inv.items()))
        names=list(source_invariant_values(d1,d2,np.ones(d1,complex),np.ones(d2,complex)))
        mode['mixed_S1_S2']=rank_report(cross_table,names,tag+'/mixed')
        if d1>1 and d2>1:rank_counts['mixed_adjoint_sectors']+=1
        if d1==d2==3:rank_counts['mixed_crossed_triplet_sectors']+=1
        special={}
        for name in SPECIAL:
            if name in regs:
                test=complex_action_holdout(name,rng,min(8,samples))
                special[name]=test
                worst_action=max(worst_action,test['max_constructed_real_Lam_action_relative_residual'])
                rank_counts['special_complex_slots']+=1
        mode['special_complex_quartics']=special
        mode['registered_quartic_count']=len(regs)
        mode['source_couplings']=sorted(regs)
        mode['EFT1_all_fermions_integrated_out']=True
        models[tag]=mode
    require(rank_counts=={'triplet_self_sectors':14,'mixed_adjoint_sectors':10,
                          'mixed_crossed_triplet_sectors':4,'special_complex_slots':12},
            'scalar quartic census differs from 16-model source family')
    require(worst_ward<TOL,'nonportal quartic SU2 Ward covariance failed')
    require(worst_action<TOL,'special quartic constructed real action failed')
    return {'format':FORMAT,'status':STATUS,'models_checked':16,
        'source_registration_slots_total':sum(len(expected_symbols(t)) for t in EXPECTED_MODELS),
        'rank_sector_counts':rank_counts,'SU2_generator_errors':generator_errors,
        'max_invariant_SU2_covariance_error':worst_ward,
        'max_special_constructed_action_error':worst_action,
        'models':models,'source_sha256':provenance,
        'certification':{
            'nonportal_real_canonical_SU2_invariant_ranks_verified':True,
            'special_complex_quartic_constructed_Hc_action_verified':True,
            'installed_native_RGBeta_Lam_compared_for_specials':False,
            'all_registered_quartics_matchete_dictionary_certified':False,
            'physical_scalar_potential_normalization_certified':False,
            'native_S2_component_frame_certified':False,
            'BCminus1_mixed_scalar_mass_diagonalization_certified':False,
            'signed_GammaCC_to_Weyl_certified':False,
            'physical_Wilson_evolution_authorized':False,
            'EFT_order':'UV -> integrate out ALL fermions -> integrate out ALL scalars -> SMEFT',
        }}
