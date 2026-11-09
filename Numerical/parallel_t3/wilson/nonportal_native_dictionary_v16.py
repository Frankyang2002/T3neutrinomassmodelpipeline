"""v16 native Matchete S1/S2 self+cross CG -> source-conditioned RGBeta convention.

The target generators are reconstructed from Matchete's *own* phased Higgs
portals. They do NOT constitute an independent native RGBeta component frame.
The native RGBeta v14 Lam/projector check is independently required. Never
promote this conditional component-frame fit to physical EFT authorization.
"""
from __future__ import annotations
from hashlib import sha256
from pathlib import Path
import json
import math
import numpy as np

from .portal_real_action_v9 import (
    EXPECTED_MODELS, DIMENSIONS, TOL, fundamental_fierz,
    spectral_triplet_adjoint, identity, relative,
)
from .portal_tensor_audit import (
    SOURCE_MAT, _sha_files, _component_tensor,
)
from .portal_orientation_v6 import _projected_portals
from .matchete_rgbeta_dictionary import matchete_number
from .s2_intertwiner_v7 import pauli_half

MAT_FORMAT = 't3-parallel-a-native-matchete-nonportal-cg-v16'
FORMAT = 't3-parallel-a-nonportal-76-conditional-dictionary-v16'
STATUS = ('NativeMatcheteNonportal76CGActionAndFullRankSourceFrameVerified_'
          'RGBetaNativeProjectors150Verified_BasisAndMassMixingConditional')
SOURCE_FILES = SOURCE_MAT + ('scripts/RunT3ParallelAExportNonportalMatcheteV16.wl',)
FAMILIES = ('self1','self2','cross')

class V16Blocked(ValueError):
    pass

def require(ok:bool, msg:str)->None:
    if not ok:
        raise V16Blocked('v16 scalar dictionary: '+msg)

def read(path:Path)->dict:
    obj=json.loads(Path(path).read_text(encoding='utf-8-sig'),
                   parse_constant=lambda x: (_ for _ in ()).throw(ValueError('nonfinite JSON '+x)))
    require(isinstance(obj,dict),'not a JSON object '+str(path))
    return obj

def sha(path:Path)->str:
    return sha256(Path(path).read_bytes()).hexdigest()

def counts(d1:int,d2:int)->dict[str,int]:
    require(d1 in (1,2,3) and d2 in (1,2,3),'unsupported dimensions')
    return {'self1':(2 if d1==3 else 1),
            'self2':(2 if d2==3 else 1),
            'cross':(3 if d1==d2==3 else 2 if d1>1 and d2>1 else 1)}

def contract(t:np.ndarray, first:np.ndarray,second:np.ndarray)->complex:
    return complex(np.einsum('ijkl,i,j,k,l->',t,first.conj(),first,
                             second.conj(),second,optimize=True))

def sym_self(t:np.ndarray)->np.ndarray:
    """Commuting identical conjugated pair and unconjugated pair."""
    return (t+t.swapaxes(0,2)+t.swapaxes(1,3)+
            t.transpose(2,3,0,1))/4

def source_record(rec:dict,family:str,k:int,d1:int,d2:int)->np.ndarray:
    require(family in FAMILIES,'bad scalar family')
    field={'self1':1,'self2':2,'cross':0}[family]
    base=('lambdaS'+str(field)) if field else 'lambda12'
    cg=('T3Scalar'+str(field)+'SelfCGInv'+str(k)) if field else 'T3CrossScalarCGInv'+str(k)
    operator=('Scalar'+str(field)+'SelfInv'+str(k)) if field else 'ScalarCrossPortalInv'+str(k)
    dim=d1 if field==1 else d2
    shape=(dim,)*4 if field else (d1,d1,d2,d2)
    bare=tuple(x for x in shape if x>1)
    require(rec.get('family')==family and rec.get('invariant')==k and
            rec.get('coupling')==base+'Inv'+str(k) and
            rec.get('operator')==operator and
            rec.get('dimensions')==list(shape) and
            rec.get('matchete_potential_factor')==('1' if field==0 else '1/2') and
            rec.get('assembled_lagrangian_contains_cg_and_coupling') is True,
            f'{family}/{k}: native Matchete operator provenance invalid')
    if not bare:
        require(rec.get('cg')=='None' and rec.get('native_dimensions')==[] and
                rec.get('native_tensor_accessor')=='Scalar singlet 1',
                'scalar singlet CG evidence invalid')
    else:
        require(rec.get('cg')==cg and rec.get('native_dimensions')==list(bare) and
                rec.get('native_tensor_accessor')=='Matchete`GetCGTensor',
                f'{family}/{k}: native GetCGTensor evidence missing')
    entries=rec.get('components_inputform')
    require(isinstance(entries,list) and len(entries)==math.prod(shape),
            f'{family}/{k}: native component inventory incomplete')
    try:
        raw=np.array([matchete_number(s) for s in entries],complex).reshape(shape)
    except (ValueError,TypeError,KeyError) as exc:
        raise V16Blocked(f'{family}/{k}: malformed native CG entry') from exc
    require(np.all(np.isfinite(raw)) and np.linalg.norm(raw)>1.e-12,
            f'{family}/{k}: native CG nonfinite or zero')
    return raw

def scalar_generators(portal_records:list[dict],which:int,d:int)->tuple[np.ndarray|None,dict]:
    if d==1:
        return None,{'dimension':1,'generator_orientation':'trivial'}
    rows=sorted((x for x in portal_records if x.get('which')==which),
                key=lambda x:x.get('invariant',-1))
    require(len(rows)==2,'missing native portal generators')
    # Full raw Matchete components and 1/2 PlusHc phases are rechecked.
    for i,r in enumerate(rows,1): _component_tensor(r,which,d,i)
    projected=_projected_portals(rows,which,d)
    if d==2:
        k=fundamental_fierz()
        ident=identity(2)
        require(np.linalg.matrix_rank(np.column_stack([x.ravel() for x in (ident,k)]))==2,
                'doublet portal singular')
        q=np.column_stack([x.ravel() for x in (ident,k)])
        for t in projected:
            res=np.linalg.norm(q@np.linalg.lstsq(q,t.ravel(),rcond=None)[0]-t.ravel())
            require(res<1e-10,'Matchete doublet portal misaligned with Fierz')
        sign=None; sres=0.
    else:
        k,sign,sres=spectral_triplet_adjoint(projected)
    th=pauli_half()
    g=2*np.einsum('aij,ijmn->amn',th.conj(),k,optimize=True)
    require(g.shape==(3,d,d),'invalid scalar generator shape')
    for a in range(3):
        require(relative(g[a],g[a].conj().T)<1e-10,'non-Hermitian source-derived generator')
    eps=np.zeros((3,3,3))
    eps[0,1,2]=eps[1,2,0]=eps[2,0,1]=1.
    eps[1,0,2]=eps[2,1,0]=eps[0,2,1]=-1.
    lie=max(relative(g[a]@g[b]-g[b]@g[a],
                     1j*sum(eps[a,b,c]*g[c] for c in range(3)))
            for a in range(3) for b in range(3))
    cas=sum(t@t for t in g)
    c=(d*d-1)/4
    require(lie<1e-10 and relative(cas,c*np.eye(d))<1e-10,
            'source-derived scalar generators fail SU2 algebra')
    trace=np.array([[np.trace(g[a]@g[b]) for b in range(3)] for a in range(3)])
    require(relative(trace,(d*(d*d-1)/12)*np.eye(3))<1e-10,
            'scalar generator trace normalization differs')
    return g,{'dimension':d,'generator_orientation':'Matchete portal-derived, not native RGBeta',
              'max_lie_relative_residual':lie,
              'Casimir':c,'S2_spectral_sign':sign,'S2_spectral_residual':sres}

def invariant_triplet_metric(g:np.ndarray)->np.ndarray:
    """Unit-normalized SU(2)-invariant symmetric spin-1 bilinear metric.

    In a Cartesian (real) triplet basis Q=I. In a spin-1 weight basis,
    Q is the Wigner conjugation anti-diagonal, not I. Its overall phase
    cancels between the barred and unbarred pieces of the crossed quartic.
    """
    require(g.shape==(3,3,3),'triplet invariant metric needs three 3x3 generators')
    basis=[]
    for k in range(3):
        for l in range(3):
            b=np.zeros((3,3),complex);b[k,l]=1
            basis.append(np.concatenate([(a.T@b+b@a).ravel() for a in g]))
    a=np.column_stack(basis)
    _,s,vh=np.linalg.svd(a,full_matrices=False)
    require(s[-1]<1.e-10 and s[-2]>1.e-9,
            'spin-1 invariant symmetric metric is not unique')
    q=vh[-1].conj().reshape(3,3)
    q=(q+q.T)/2
    require(np.linalg.norm(q)>1.e-10,'vanishing triplet invariant metric')
    q=q*np.sqrt(3.)/np.linalg.norm(q)
    first=next((x for x in q.ravel() if abs(x)>1.e-9),None)
    require(first is not None,'no nonzero triplet metric entry')
    q=q/(first/abs(first))
    require(relative(q,q.T)<1.e-10 and
            relative(q.conj().T@q,np.eye(3))<1.e-10 and
            max(relative(t.T@q+q@t,np.zeros_like(q)) for t in g)<1.e-10,
            'triplet invariant metric violates SU2 covariance/normalization')
    return q


def target_tensors(family:str,d1:int,d2:int,g1:np.ndarray|None,
                   g2:np.ndarray|None)->tuple[list[str],list[np.ndarray]]:
    if family=='self1' or family=='self2':
        which=1 if family=='self1' else 2
        d=d1 if which==1 else d2
        g=g1 if which==1 else g2
        ident=np.einsum('ij,kl->ijkl',np.eye(d),np.eye(d)).astype(complex)/2
        names=[f'lambdaS{which}'];tensors=[ident]
        if d==3:
            names.append(f'lambdaS{which}Adj')
            tensors.append(np.einsum('aij,akl->ijkl',g,g,optimize=True)/2)
        return names,[sym_self(t) for t in tensors]
    require(family=='cross','bad family')
    ident=np.einsum('ij,kl->ijkl',np.eye(d1),np.eye(d2)).astype(complex)
    names=['lambda12'];tensors=[ident]
    if d1>1 and d2>1:
        names.append('lambda12Adj')
        tensors.append(np.einsum('aij,akl->ijkl',g1,g2,optimize=True))
    if d1==d2==3:
        # Both triplets share Matchete's single T3BSMd3 representation index.
        require(relative(g1,g2)<1e-10,
                'Matchete triplets use different bases; cannot assert crossed delta')
        q=invariant_triplet_metric(g1)
        # First delta pairs two barred and two unbarred spin-1 legs.
        # When transported from RGBeta's Cartesian frame to Matchete's
        # potentially spherical frame, it is Q* (barred) times Q (unbarred).
        delta=np.eye(3)
        crossed=(np.einsum('ik,jl->ijkl',q.conj(),q)+
                 np.einsum('il,jk->ijkl',delta,delta)).astype(complex)
        names.append('lambda12Cross');tensors.append(crossed)
    return names,tensors

def audit_sector(tag:str,family:str,records:list[dict],d1:int,d2:int,
                 g1:np.ndarray|None,g2:np.ndarray|None,
                 *,rng:np.random.Generator,samples:int=5)->dict:
    target_names,target=target_tensors(family,d1,d2,g1,g2)
    require(len(records)==len(target)==counts(d1,d2)[family],
            f'{tag}/{family}: missing native invariant slots')
    source=[]
    for k,rec in enumerate(sorted(records,key=lambda x:x.get('invariant',-1)),1):
        t=source_record(rec,family,k,d1,d2)
        if family.startswith('self'):
            t=sym_self(t)
        else:
            # Cross is a real coupling with no PlusHc: action must itself be real.
            require(relative(t,t.conj().transpose(1,0,3,2))<1.e-10,
                    f'{tag}/{family}/{k}: real Matchete coupling multiplies non-Hermitian CG')
        factor=.5 if family.startswith('self') else 1.
        source.append(factor*t)
    # Independently check target has correct commuting-field symmetry & rank.
    T=np.stack([x.ravel() for x in target],axis=1)
    M=np.stack([x.ravel() for x in source],axis=1)
    s=np.linalg.svd(T,compute_uv=False)
    require(s[-1]>1.e-10*s[0] and s[0]/s[-1]<1e8,
            f'{tag}/{family}: RGBeta target quartic basis not independent')
    z=np.linalg.lstsq(T,M,rcond=None)[0]
    require(np.linalg.matrix_rank(z,tol=1e-10)==len(records),
            f'{tag}/{family}: Matchete coupling basis loses a direction')
    require(np.linalg.norm(z.imag)<1e-10*max(1.,np.linalg.norm(z)),
            f'{tag}/{family}: no real Matchete-to-RGBeta coefficient map')
    fit=relative(T@z,M)
    require(fit<1.e-10,f'{tag}/{family}: Matchete native CG not in RGBeta invariant span ({fit:.3e})')
    a=z.real
    worst=fit
    n1=d1 if family=='cross' else (d1 if family=='self1' else d2)
    n2=d2 if family=='cross' else n1
    for _ in range(samples):
        v=rng.normal(size=n1)+1j*rng.normal(size=n1)
        w=rng.normal(size=n2)+1j*rng.normal(size=n2)
        if family.startswith('self'): w=v
        c=rng.normal(size=len(records))
        mat=sum(t*x for t,x in zip(source,c))
        rgb=sum(t*x for t,x in zip(target,a@c))
        mact=contract(mat,v,w);ract=contract(rgb,v,w)
        require(abs(mact.imag)<1.e-9*max(1.,abs(mact)) and
                abs(ract.imag)<1.e-9*max(1.,abs(ract)),
                f'{tag}/{family}: real scalar potential not real')
        worst=max(worst,relative(mact,ract))
    require(worst<1.e-10,f'{tag}/{family}: independent complex-field action mismatch')
    return {'family':family,'matchete_couplings':
            [rec['coupling'] for rec in sorted(records,key=lambda x:x['invariant'])],
            'rgbeta_couplings':target_names,
            'conditional_rgbeta_from_matchete_real_matrix':a.tolist(),
            'source_tensor_fit_relative_residual':fit,
            'max_independent_complex_field_action_residual':worst,
            'native_matchete_sector_rank':len(records),
            'constructed_rgbeta_sector_rank':len(target),
            'physical_coupling_dictionary_certified':False}

def audit(mat:dict,portals:dict,*,samples:int=5,seed:int=1618)->dict:
    require(1<=samples<=30 and type(seed) is int,'bad audit settings')
    require(mat.get('format')==MAT_FORMAT and mat.get('nonportal_slots')==76 and
            mat.get('models_checked')==16 and
            set(mat.get('models',{}))==EXPECTED_MODELS and
            set(portals.get('models',{}))==EXPECTED_MODELS,
            'native Matchete nonportal+portal 16-model evidence incomplete')
    require(mat.get('physical_scalar_dictionary_certified') is False and
            mat.get('physical_Wilson_evolution_authorized') is False and
            mat.get('EFT_order')==
             'UV -> integrate out ALL fermions -> integrate out ALL scalars -> SMEFT',
            'unsafe source claim or EFT order')
    rng=np.random.default_rng(seed)
    models={};total=0;worst=0.;mixed=[];genchecks={};sector_count=0
    for tag in sorted(EXPECTED_MODELS):
        d1,d2,df=DIMENSIONS[tag[0]]
        native=mat['models'][tag];old=portals['models'][tag]
        require(native.get('dimensions')==[d1,d2,df] and
                old.get('dimensions')==[d1,d2,df] and
                native.get('alpha')==int(tag.split('_')[1]) and
                native.get('scalar_cg_operator_inventory_verified') is True and
                native.get('native_build_status')==old.get('native_build_status') and
                native.get('matching_mass_basis_ready')==old.get('matching_mass_basis_ready'),
                f'{tag}: inconsistent native Matchete builds/representation')
        if tag in ('B_-1','C_-1'):
            require(native.get('native_build_status')==
                    'RequiresMixedScalarMassDiagonalization' and
                    native.get('mixed_scalar_mass_basis_blocked') is True and
                    native.get('matching_mass_basis_ready') is False,
                    f'{tag}: mixed scalar mass basis inaccurately certified')
            mixed.append(tag)
        else:
            require(native.get('native_build_status')=='Success' and
                    native.get('matching_mass_basis_ready') is True,
                    f'{tag}: unexpected Matchete model build failure')
        g1,c1=scalar_generators(old['portals'],1,d1)
        g2,c2=scalar_generators(old['portals'],2,d2)
        genchecks[tag]={'S1':c1,'S2':c2}
        records=native.get('nonportal')
        require(isinstance(records,list) and len(records)==sum(counts(d1,d2).values()),
                f'{tag}: incorrect native nonportal inventory')
        byfamily={}
        for fam in FAMILIES:
            rr=[x for x in records if x.get('family')==fam]
            row=audit_sector(tag,fam,rr,d1,d2,g1,g2,rng=rng,samples=samples)
            byfamily[fam]=row;total+=len(rr);sector_count+=1
            worst=max(worst,row['max_independent_complex_field_action_residual'])
        models[tag]=byfamily
    require(total==76 and sector_count==48 and len(mixed)==2,
            '76/48/B-C census invalid')
    return {'format':FORMAT,'status':STATUS,'models_checked':16,
            'nonportal_coupling_slots':total,'nonportal_sectors':sector_count,
            'max_relative_residual':worst,'models':models,
            'source_portal_generator_checks':genchecks,
            'mixed_scalar_mass_basis_excluded':mixed,
            'certification':{
                'native_Matchete_76_nonportal_CG_source_registered':True,
                'real_Matchete_nonportal_potential_in_abstract_RGBeta_span':True,
                'native_RGBeta_S2_numerical_component_frame_measured':False,
                'full_physical_Matchete_RGBeta_coupling_map_certified':False,
                'mixed_BC_minus1_mass_basis_certified':False,
                'signed_GammaCC_to_Weyl_certified':False,
                'physical_Wilson_evolution_authorized':False,
                'EFT_order':'UV -> integrate out ALL fermions -> integrate out ALL scalars -> SMEFT'
            }}

def audit_from_paths(mat_path:Path,portal_path:Path,native14_path:Path,root:Path,
                     *,samples:int=5)->dict:
    root=Path(root).resolve()
    mat=read(mat_path);portal=read(portal_path)
    require(mat.get('status')=='NativeMatchete76NonportalCGComponentsCaptured_MassBasisStillConditional',
            'native Matchete v16 exporter status invalid')
    # Revalidate Matchete model sources and scripts against actual local bytes.
    _sha_files(root,mat,SOURCE_FILES)
    _sha_files(root,portal,SOURCE_MAT)
    from .native_context_quartic_v14 import audit_witness
    native14=read(native14_path)
    verified14=audit_witness(native14,root)
    require(verified14.get('status') is not None and
            native14.get('smoke_only') is False and
            native14.get('status')=='NativeRGBetaRealQuarticLamProjectorResponsesObserved_Conditional',
            'full real v14 native Lam/projector witness absent or blocked')
    report=audit(mat,portal,samples=samples)
    report['native_rgbeta_v14_real_projector_witness_validated']=True
    report['input_sha256']={'native_matchete_nonportal_v16':sha(mat_path),
                            'native_matchete_portal_v5':sha(portal_path),
                            'native_rgbeta_quartic_v14':sha(native14_path)}
    return report
