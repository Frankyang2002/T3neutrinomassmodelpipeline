"""Batch 13: source-tree, full-pair-tensor candidate EFT1 Wilson mixing for A(-4), D(-2).

Hierarchy: UV -> all F at 1e5 -> EFT1(S1,S2) -> both S at 1e3 -> SMEFT.

At strict O(1/16pi^2), compute beta[C_tree(muF), renormalizable(muF)]
for the FULL 45-Weyl x 12-real-scalar psi^2 phi^2 tensor. Seven named
channels are retained, though the crossed/gauge-complete anomalous dimension
is NOT an independently certified physical ADM. Source scalar quartics are
expanded with Matchete's exact CGs. D's two *additional* HHSS reference
invariants are included as explicitly conditional candidate mappings.

Replace, NEVER ADD TO, Batch11's pole-derived EFT1 logarithm. Both computations
remain independent diagnostics of the unresolved native-to-master normalization.
This is not a complete physical C5 or EFT1 operator-closure certificate.
"""
from __future__ import annotations

from hashlib import sha256
from itertools import permutations, product
from pathlib import Path
from typing import Any, Mapping
import json
import math

import numpy as np
import sympy as sp

MODELS = ('A_-4', 'D_-2')
FORMAT = 't3-a4-d2-master-wilson-strict-one-loop-batch13-v1'
B9 = Path('output/parallel_t3/thesis_full_general_joint_scalar_batch9_v1')
B11 = Path('output/parallel_t3/thesis_native_laurent_batch11_v1')
B12 = Path('output/parallel_t3/thesis_conditional_uvfit_batch12_v1')
D_EXTRA = Path('output/general_observable/lo_certification_couplings/ws2_quartic_completion_v1/D_m2/exact_extra_HHSS_reference_tensors.json')
REGISTRY = Path('RGE/general/generated/general_rgbeta_quartic_registry.json')
SM_KEYS = ('gY', 'g2', 'g3', 'lambdaH', 'yu', 'yd', 'ye')
LOOP = 1/(16*math.pi**2)
BASIS_SIZES = {'A_-4': (1,3), 'D_-2': (3,1)}
# Explicit, provisional map. A missing entry is a blocker, not an inferred zero.
CANDIDATE_SOURCE_ALIASES = {
    'λ': 'lambdaH',
    'lambdaT3': 'lambdaT3',
    'lambda12Inv1': 'lambda12',
    'lambda12Inv2': 'lambda12Adj',
    'lambdaH1Inv1': 'lambdaH1',
    'lambdaH1Inv2': 'lambdaH1Adj',
    'lambdaH2Inv1': 'lambdaH2',
    'lambdaH2Inv2': 'lambdaH2Adj',
    'lambdaS1Inv1': 'lambdaS1',
    'lambdaS1Inv2': 'lambdaS1Adj',
    'lambdaS2Inv1': 'lambdaS2',
    'lambdaS2Inv2': 'lambdaS2Adj',
}


def require(ok: bool, reason: str) -> None:
    if not ok:
        raise ValueError('Batch13: '+reason)


def digest(path: Path) -> str:
    require(path.is_file(), 'Missing '+str(path))
    return sha256(path.read_bytes()).hexdigest()


def decode_complex(value: Any) -> np.ndarray:
    """Decode Batch1 NPZ trajectory values *or* native serialized JSON values.

    NativeBatch1Trajectory.at_fermion_threshold() returns NumPy arrays,
    including zero-dimensional complex arrays for scalar/gauge couplings.
    general_batch1_api.decode_native() intentionally handles JSON only.
    Retain that strict JSON reader for serialized inputs, and accept the
    already-decoded in-memory numerical state without converting it to JSON.
    """
    if isinstance(value, (np.ndarray, np.generic)):
        arr=np.asarray(value,dtype=np.complex128)
    elif isinstance(value, Mapping) and set(value)=={'real','imag'}:
        re=np.asarray(value['real'],dtype=float)
        im=np.asarray(value['imag'],dtype=float)
        require(re.shape==im.shape,'Mismatched real/imag native coupling shapes')
        arr=re+1j*im
    else:
        from Numerical.parallel_t3.background.general_batch1_api import decode_native
        arr=np.asarray(decode_native(value),dtype=np.complex128)
    require(np.isfinite(arr).all(),'Nonfinite native coupling')
    return arr


def encode_complex(arr: Any) -> dict:
    a=np.asarray(arr,complex)
    require(np.isfinite(a).all(), 'Nonfinite complex output')
    return {'real':a.real.tolist(),'imag':a.imag.tolist()}


def symmetric_matrix(a: Any) -> np.ndarray:
    x=np.asarray(a,complex)
    require(x.shape==(3,3) and np.isfinite(x).all(), 'Expected finite 3x3 matrix')
    return (x+x.T)/2


def replace_eft1_log(native_finite: np.ndarray, native_log: np.ndarray,
                     master_log: np.ndarray, master_to_native_factor: float
                     ) -> np.ndarray:
    """Strict O(hbar): original native hard+soft regular PLUS ONE EFT1 log."""
    mats=[symmetric_matrix(z) for z in (native_finite,native_log,master_log)]
    f=float(master_to_native_factor)
    require(np.isfinite(f) and f!=0, 'Explicit nonzero master/native normalization needed')
    return symmetric_matrix(mats[0]-mats[1]+f*mats[2])


def source_to_master(raw: np.ndarray, fermion_labels: tuple[str,...],
                     heavy_rotation: np.ndarray, model: str) -> np.ndarray:
    """Pure formal tensor embedding, NOT a claimed GammaCC/Weyl conversion.

    raw[p,q,i,j,a,b] has two flavour, two SU2 and two heavy-real coordinates.
    source heavy phi_a = U[a,A] * matching heavy phi_A.
    """
    d1,d2=BASIS_SIZES[model]
    n=2*(d1+d2)
    z=np.asarray(raw,complex)
    U=np.asarray(heavy_rotation,float)
    require(z.shape==(3,3,2,2,n,n) and U.shape==(n,n),
            'Wrong source F-tree or scalar rotation dimensions')
    require(np.isfinite(z).all() and np.isfinite(U).all() and
            np.allclose(U.T@U,np.eye(n),rtol=1e-10,atol=1e-10),
            'Source-to-matching scalar rotation must be real orthogonal')
    # The SourceBoundary is separately symmetric under lepton/scalar exchange.
    require(np.allclose(z,z.transpose(1,0,3,2,4,5),rtol=1e-9,atol=1e-25)
            and np.allclose(z,z.swapaxes(4,5),rtol=1e-9,atol=1e-25),
            'Source F-tree is not a symmetric polynomial tensor')
    moved=np.einsum('pqijab,aA,bB->pqijAB',z,U,U,optimize=True)
    fids={v:i for i,v in enumerate(fermion_labels)}
    indices=[fids[f'L[{p},{i}]'] for p in (1,2,3) for i in (1,2)]
    out=np.zeros((len(fermion_labels),len(fermion_labels),4+n,4+n),complex)
    out[np.ix_(indices,indices,range(4,4+n),range(4,4+n))]=(
        moved.transpose(0,2,1,3,4,5).reshape(6,6,n,n))
    require(np.allclose(out,out.swapaxes(0,1),rtol=1e-9,atol=1e-25) and
            np.allclose(out,out.swapaxes(2,3),rtol=1e-9,atol=1e-25),
            'Formal master Wilson tensor lost operator exchange symmetry')
    return out


def sm_yukawa_tensor(couplings: Mapping[str,Any], fermion_labels: tuple[str,...]) -> np.ndarray:
    """SM left-Weyl Y_ija tensor; exact SU3/weak assignments from Aminus4EFT1Couplings."""
    fi={name:i for i,name in enumerate(fermion_labels)}
    Y=np.zeros((len(fi),len(fi),12),complex)
    values={name:decode_complex(couplings[name]) for name in ('yu','yd','ye')}
    require(all(x.shape==(3,3) and np.isfinite(x).all() for x in values.values()),
            'SM Yukawas must be finite full 3x3 matrices')
    def put(i,j,hcomponent,v,conjugated):
        a=2*hcomponent
        Y[i,j,a]+=v/math.sqrt(2)
        Y[i,j,a+1]+=(-1j if conjugated else 1j)*v/math.sqrt(2)
        Y[j,i,a]=Y[i,j,a]
        Y[j,i,a+1]=Y[i,j,a+1]
    for g,h in product(range(1,4),repeat=2):
        for col in range(1,4):
            for w in range(1,3):
                for hc in range(1,3):
                    eps=1 if (w,hc)==(1,2) else -1 if (w,hc)==(2,1) else 0
                    if eps:
                        put(fi[f'Q[{g},{col},{w}]'],fi[f'uc[{h},{col}]'],hc-1,
                            eps*values['yu'][g-1,h-1],False)
                    if w==hc:
                        put(fi[f'Q[{g},{col},{w}]'],fi[f'dc[{h},{col}]'],hc-1,
                            values['yd'][g-1,h-1],True)
        for w in range(1,3):
            put(fi[f'L[{g},{w}]'],fi[f'ec[{h}]'],w-1,
                values['ye'][g-1,h-1],True)
    return Y


def gauge_sectors(couplings: Mapping[str,Any], model: str,
                  heavy_rotation: np.ndarray) -> tuple:
    """45-Weyl SM generators and 12-real scalar generators, including SU3c."""
    from Numerical.running.Aminus4EFT1Couplings import (
        gauge_generators, real_scalar_generator, su2_generators)
    orig=gauge_generators()
    d1,d2=BASIS_SIZES[model]
    n=2*(d1+d2)
    rot=np.eye(12)
    rot[4:,4:]=np.asarray(heavy_rotation,float)
    heavy_reps=(('H',2,.5,0),('S1',d1,-2. if model=='A_-4' else -1.,4),
                ('S2',d2,-1. if model=='A_-4' else 0.,4+2*d1))
    scalar={'SU3':np.zeros((8,12,12),complex),
            'SU2':np.zeros((3,12,12),complex),
            'U1':np.zeros((1,12,12),complex)}
    for label,d,y,start in heavy_reps:
        ids=np.arange(start,start+2*d)
        for i,t in enumerate(su2_generators(d)):
            scalar['SU2'][i][np.ix_(ids,ids)]=real_scalar_generator(t)
        scalar['U1'][0][np.ix_(ids,ids)]=real_scalar_generator(y*np.eye(d))
    pairs=[]
    for name,sector in (('g2','SU2'),('gY','U1'),('g3','SU3')):
        g=complex(np.asarray(decode_complex(couplings[name])))
        require(abs(g.imag)<1e-10 and g.real>0, 'Invalid gauge coupling '+name)
        term=tuple((np.asarray(orig[sector]['fermion'][i],complex),
                    rot.T@scalar[sector][i]@rot)
                   for i in range(len(scalar[sector])))
        pairs.append((float(g.real),term))
    return tuple(pairs)


def scalar_quartic_from_matchete(seed: Mapping[str,Any], couplings: Mapping[str,Any],
                                 model: str, heavy_rotation: np.ndarray,
                                 native_source_from_reference: np.ndarray | None = None
                                 ) -> tuple[np.ndarray,dict]:
    """CG-exact native topology quartic, supplemented by all baseline RGBeta columns.

    The *archived* native quartic seed does NOT contain all RG-closed portal/self
    invariants. Reusing it alone would silently discard lambdaH2 in A or
    lambdaH1 in D, for example. Here every allowed baseline RGBeta invariant
    enters explicitly. Only the topology lambdaT3 uses native CGs; its
    RGBeta-normalization still requires independent verification.
    """
    from RGE.running.intermediate.ScalarOnlyTensorAdapters import build_eft1_quartic_tensor
    from Numerical.running.Aminus4EFT1Couplings import (
        _quartic_product,su2_generators)
    d1,d2=BASIS_SIZES[model]
    require(seed.get('ScalarQuarticTermCount')==len(seed.get('ScalarQuarticTerms',[])),
            'Native scalar source inventory malformed')
    names={}
    for term in seed['ScalarQuarticTerms']:
        import re
        a=re.findall(r'Coupling\[([^,\]]+),',term['TermInputForm'])
        require(len(a)==1,'Matchete scalar term has an ambiguous coupling')
        name=a[0]
        target=CANDIDATE_SOURCE_ALIASES.get(name)
        if target is None and name.startswith('lambdaT3ExtraQ'):
            target=name.replace('lambdaT3ExtraQ','lambdaT3GeneralQ',1)
        require(target is not None and target in couplings,
                'Missing scalar native-to-RGBeta map for '+name)
        names[name]=target
    mix=[x for x in seed['ScalarQuarticTerms'] if 'Coupling[lambdaT3,' in x['TermInputForm']]
    require(mix,'No native topology quartic lambdaT3')
    one={'CGRegistry':seed['CGRegistry'],'ScalarQuarticTerms':mix}
    sparse=build_eft1_quartic_tensor(one,d_s1=d1,d_s2=d2,shared_scalar=False)
    lam=complex(np.asarray(decode_complex(couplings['lambdaT3'])))
    source_values={sp.Symbol('lambdaT3'):lam}
    source_mix=np.zeros((12,)*4,float)
    for key,expr in sparse.nonzero_items():
        z=complex(sp.N(expr.subs(source_values),16))
        require(np.isfinite([z.real,z.imag]).all() and
                abs(z.imag)<1e-9*max(1.,abs(z.real)),
                'Topology quartic did not produce a real tensor')
        for orbit in set(permutations(tuple(i-1 for i in key))):
            source_mix[orbit]=z.real
    # Build complete baseline real-scalar V4 fourth-derivative columns.
    # Matches the repository Aminus4EFT1Couplings.quartic_columns definitions.
    def hessian(t,offset):
        t=np.asarray(t,complex)
        d=t.shape[0]
        local=np.zeros((2*d,2*d),float)
        for p in range(d):
            for q in range(d):
                x,y=t[p,q].real,t[p,q].imag
                local[2*p:2*p+2,2*q:2*q+2]=[[x,-y],[y,x]]
        x=np.zeros((12,12),float)
        x[offset:offset+2*d,offset:offset+2*d]=local
        return x
    sh=hessian(np.eye(2),0)
    s1=hessian(np.eye(d1),4)
    s2=hessian(np.eye(d2),4+2*d1)
    sm1=[hessian(t,4) for t in su2_generators(d1)]
    sm2=[hessian(t,4+2*d1) for t in su2_generators(d2)]
    hm=[hessian(t,0) for t in su2_generators(2)]
    columns={
       'lambdaH':.5*_quartic_product(sh,sh),
       'lambdaS1':.5*_quartic_product(s1,s1),
       'lambdaS2':.5*_quartic_product(s2,s2),
       'lambdaH1':_quartic_product(sh,s1),
       'lambdaH2':_quartic_product(sh,s2),
       'lambda12':_quartic_product(s1,s2),
    }
    if d1>1:
        columns['lambdaS1Adj']=sum((.5*_quartic_product(x,x) for x in sm1),
                                    np.zeros((12,)*4))
        columns['lambdaH1Adj']=sum((_quartic_product(h,s) for h,s in zip(hm,sm1)),
                                    np.zeros((12,)*4))
    if d2>1:
        columns['lambdaS2Adj']=sum((.5*_quartic_product(x,x) for x in sm2),
                                    np.zeros((12,)*4))
        columns['lambdaH2Adj']=sum((_quartic_product(h,s) for h,s in zip(hm,sm2)),
                                    np.zeros((12,)*4))
    source_baseline=np.zeros((12,)*4,float)
    for name,base in columns.items():
        require(name in couplings,'Missing independent RG-closed scalar term '+name)
        z=complex(np.asarray(decode_complex(couplings[name])))
        require(abs(z.imag)<1e-10*max(1.,abs(z.real)),
                'Real RG-closed scalar coupling contains imaginary phase: '+name)
        source_baseline+=z.real*base
    # Native Matchete T3MixCG is in a Cartesian-like triplet basis, while
    # baseline RGBeta columns above are ALREADY in spherical spin-1 basis.
    # Transform ONLY the native source terms into that reference basis;
    # applying the triplet rotation to the RGBeta reference terms again would
    # silently mix two incompatible coordinate conventions.
    native_R=(np.eye(8) if native_source_from_reference is None
              else np.asarray(native_source_from_reference,float))
    mass_R=np.asarray(heavy_rotation,float)
    require(native_R.shape==mass_R.shape==(8,8) and
            np.allclose(native_R.T@native_R,np.eye(8),atol=1e-10) and
            np.allclose(mass_R.T@mass_R,np.eye(8),atol=1e-10),
            'Native/reference or D mass-basis rotation not orthogonal')
    nrot=np.eye(12);nrot[4:,4:]=mass_R@native_R
    rrot=np.eye(12);rrot[4:,4:]=mass_R
    native_contribution=np.einsum('abcd,aA,bB,cC,dD->ABCD',
                                source_mix,nrot,nrot,nrot,nrot,optimize=True)
    reference_contribution=np.einsum('abcd,aA,bB,cC,dD->ABCD',
                                   source_baseline,rrot,rrot,rrot,rrot,optimize=True)
    result=native_contribution+reference_contribution
    require(np.isfinite(result).all(),'Nonfinite rotated scalar quartic tensor')
    return result,{
       'source_native_to_RGBeta_candidate_alias_inventory':names,
       'baseline_RG_closed_RGBeta_columns':sorted(columns),
       'missing_archived_portal_and_self_columns_completed_from_RGBeta':True,
       'native_topology_cg_tensor_contribution_included':True,
       'source_Cartesian_to_spherical_rotation_applied_only_to_native_quartics':native_source_from_reference is not None,
       'matchete_RGBeta_native_normalization_certified':False,
       'source_scalar_seed_quartic_count':len(seed['ScalarQuarticTerms'])}


def add_d2_extra_scalar_only(full: np.ndarray, couplings: Mapping[str,Any],
                             heavy_rotation: np.ndarray | None = None) -> dict:
    """Three D(-2) scalar-only RGBeta invariants absent from the archived seed.

    This is a *reference polynomial*, not a certified native/Matchete basis.
    The existence of these three couplings is fixed by the model's RGBeta
    registry. Their phases/orientation require a separate source check.
    """
    names=['lambdaT3GeneralQ000004I1','lambdaT3GeneralQ000013I1',
           'lambdaT3GeneralQ001102I1']
    require(full.shape==(12,)*4,'Full D quartic shape invalid')
    reference=np.zeros((12,)*4,float)
    xyz=sp.symbols('u0:8',real=True)
    z=(xyz[6]+sp.I*xyz[7])/sp.sqrt(2)
    zb=sp.conjugate(z)
    s1sq=sum((xyz[k]**2+xyz[k+1]**2)/2 for k in range(0,6,2))
    polys=(zb**4,z*zb**3,s1sq*zb**2)
    for name,P in zip(names,polys):
        require(name in couplings,'Missing D extra scalar-only quartic '+name)
        c=complex(np.asarray(decode_complex(couplings[name])))
        # Reference invariant V4=-cP-c*P*, matching existing D extra-HHSS
        # geometry convention, not assuming any absolute RGBeta->Matchete map.
        pot=sp.expand(-sp.Float(c.real,17)*(P+sp.conjugate(P))-
                      sp.I*sp.Float(c.imag,17)*(P-sp.conjugate(P)))
        poly=sp.Poly(pot,*xyz)
        for a in range(8):
            for b in range(a,8):
                for cc in range(b,8):
                    for d in range(cc,8):
                        coords=(a,b,cc,d)
                        monomial=sp.prod(xyz[i] for i in coords)
                        val=complex(sp.N(poly.coeff_monomial(monomial),15))
                        val*=math.prod(math.factorial(coords.count(k))
                                       for k in set(coords))
                        require(abs(val.imag)<1e-9,'Nonreal D scalar-only reference quartic')
                        if val.real:
                            for orbit in set(permutations(tuple(4+i for i in coords))):
                                reference[orbit]+=val.real
    R=np.eye(8) if heavy_rotation is None else np.asarray(heavy_rotation,float)
    require(R.shape==(8,8) and np.isfinite(R).all() and
            np.allclose(R.T@R,np.eye(8),atol=1e-10),
            'Invalid D split real scalar transform')
    transform=np.eye(12)
    transform[4:,4:]=R
    full+=np.einsum('abcd,aA,bB,cC,dD->ABCD',reference,
                    transform,transform,transform,transform,optimize=True)
    return {'D2_scalar_only_extra_reference_invariants':names,
            'D2_extra_scalar_only_count':3,
            'D2_scalar_only_extra_RGBeta_to_Matchete_map_certified':False}


def add_d2_extra_hhss(full: np.ndarray, root: Path, couplings: Mapping[str,Any],
                       heavy_rotation: np.ndarray) -> dict:
    """Add the two missing D HHSS reference directions, each once per orbit."""
    from Numerical.parallel_t3.complete.two_model_hard_general_batch10 import (
        extra_hhss_from_source)
    doc=json.loads((root/D_EXTRA).read_text(encoding='utf-8-sig'))
    reg=json.loads((root/REGISTRY).read_text(encoding='utf-8-sig'))
    values={k:decode_complex(v).item() for k,v in couplings.items()
            if k.startswith('lambdaT3GeneralQ')}
    channels,audit=extra_hhss_from_source(doc,reg,values)
    for q in channels.values():
        qq=np.einsum('ABij,ia,jb->ABab',q,heavy_rotation,heavy_rotation,
                     optimize=True)
        # Independent sorted H,H,S,S combinations prevent over-counting orbit entries.
        for a in range(4):
            for b in range(a,4):
                for c in range(8):
                    for d in range(c,8):
                        v=float(qq[a,b,c,d])
                        if v==0: continue
                        for orbit in set(permutations((a,b,c+4,d+4))):
                            full[orbit]+=v
    return {'extra_hhss_reference_channels':sorted(channels),
            'extra_hhss_basis_mapping_certified':False,
            'remaining_D2_extra_scalar_only_quartic_count':3,
            'source_extra_reference_sha256':digest(root/D_EXTRA)}


def projection(c: np.ndarray, fermion_labels: tuple[str,...]) -> np.ndarray:
    """Project canonical master coefficient: C5 = 4 C[nu_p,nu_q,h0R,h0R]."""
    nu=[fermion_labels.index(f'L[{p},1]') for p in range(1,4)]
    val=4*np.asarray(c[np.ix_(nu,nu,[2],[2])][:,:,0,0],complex)
    return symmetric_matrix(val)


def _relative(a: np.ndarray,b: np.ndarray)->float:
    return float(np.linalg.norm(a-b)/max(np.linalg.norm(a),np.linalg.norm(b),1e-100))


def run_one(root: Path, model: str, *, master_to_native:float=1.,
            signed_factor:float=1.,beta_mode:str='fonseca_eq22_six_term',
            output:Path|None=None, include_smeft:bool=True)->dict:
    """Compute strict one-loop master beta and substitute its log at muS.

    The master normalization and candidate full ADM are NOT certified. The
    conditional SMEFT run is for a numerical pipeline demonstration only.
    """
    root=Path(root).resolve()
    require(model in MODELS,'Unsupported model')
    require(beta_mode in ('fonseca_eq22_six_term','ward_covariant_candidate'),
            'The source_v1_unstarred beta is known to fail gauge covariance')
    from Numerical.parallel_t3.complete.two_model_joint_scalar_batch9 import load_snapshot
    from Numerical.running.GeneralT3DynamicFermionWilson import uv_from_json,nonmixed_boundary
    from Numerical.joint_scalar_matching.large_batch import load_background
    from Numerical.running.Aminus4EFT1Couplings import F_LABELS
    from Numerical.parallel_t3.wilson.adm import beta, SYMMETRIC_TERMS
    from Numerical.parallel_t3.complete.full_general_completion import existing_evidence

    b9=root/B9/model/'result.json'
    b11=root/B11/model/'result.json'
    a=json.loads(b9.read_text(encoding='utf-8-sig'))
    b=json.loads(b11.read_text(encoding='utf-8-sig'))
    require(a.get('status')=='SourceBasedLLSSScalarLoopComputed_ConditionalNotPhysicalC5'
            and b.get('status')=='NativeHardSoftLaurentAndScaleLogClosed_ConditionalOperatorConvention',
            'Need successful native Batch9 and Batch11 for model '+model)
    snap=load_snapshot(root,model)
    uv=uv_from_json(model,snap.uv_values,mu_f=snap.running.mu_f)
    source,source_info=nonmixed_boundary(model,uv,repo_root=root,mu_f=snap.running.mu_f)
    require(source.raw.shape==(3,3,2,2,8,8),'Wrong source tree dimensions')
    slug=model.replace('_-','_m')
    bg=load_background(root/B9/model/(slug+'.json'),source)
    require(bg.model==model and bg.mu_s==1e3 and bg.mu_f==1e5,
            'Wrong joint scalar background')
    # Ordinary A: same real basis. D: original neutral complex singlet becomes
    # separately split, physical real S2- and S2+ fields; source U is not guessed.
    R=np.eye(8)
    if model=='D_-2':
        from Numerical.joint_scalar_matching.special_mass_v6 import build_special_basis
        from Numerical.joint_scalar_matching.source_quartic_v4 import build_source_channels
        ch=build_source_channels(snap.quartic_seed,model)
        spec=build_special_basis(model,snap.running,source)
        R=np.asarray(spec['R_original_from_boundary'],float)
    from Numerical.parallel_t3.complete.native_triplet_basis_batch15 import (
        heavy_source_from_reference,triplet_cg_from_seed,generators_from_cg,
        relative_operator_ward)
    seed=json.loads(snap.quartic_seed.read_text(encoding='utf-8-sig'))
    source_seed=json.loads(source.seed_path.read_text(encoding='utf-8-sig'))
    native_R,basis_info=heavy_source_from_reference(model,source_seed)
    qcg=generators_from_cg(triplet_cg_from_seed(seed))
    wcg=generators_from_cg(triplet_cg_from_seed(source_seed))
    require(max(np.linalg.norm(a-b) for a,b in zip(qcg,wcg))<1e-10,
            'Native quartic and Wilson tree use different triplet CG bases')
    # The saved Batch 9 original source-to-boundary map is identity for
    # A and the already split D source. The additional native-source triplet
    # map converts these to the SAME spherical basis as the RGBeta columns.
    require(np.allclose(bg.U,np.eye(8),rtol=1e-10,atol=1e-10),
            'Unexpected saved scalar source-to-boundary basis map')
    full_fermions=tuple(F_LABELS)
    c0=source_to_master(source.raw,full_fermions,R@native_R,model)
    fstate=snap.trajectory.at_fermion_threshold()
    # The fixed-order one-loop beta is frozen at muF, not the muS state.
    yukawa=sm_yukawa_tensor(fstate,full_fermions)
    sectors=gauge_sectors(fstate,model,R)
    quartic,quartic_info=scalar_quartic_from_matchete(seed,fstate,model,R,
                                                    native_source_from_reference=native_R)
    d2_info={}
    if model=='D_-2':
        d2_info=add_d2_extra_hhss(quartic,root,fstate,R)
        d2_info.update(add_d2_extra_scalar_only(quartic,fstate,R))
    assert quartic.shape==(12,)*4
    require(all(np.allclose(quartic,quartic.swapaxes(0,i),atol=1e-10)
                for i in (1,2,3)), 'Scalar quartic lacks S4 exchange symmetry')
    tree_ward=max(relative_operator_ward(c0,tf,ts)
                  for _,gens in sectors for tf,ts in gens)
    require(tree_ward<1e-8,
            f'CG-aligned source LLSS still violates SU2/U1/SU3 Ward identity: {tree_ward:.5g}')
    channels,backend=beta(c0,yukawa,quartic,sectors,variant=beta_mode)
    require(set(channels)==set(SYMMETRIC_TERMS),'Candidate seven-channel ADM lost one term')
    dt=math.log(1e3/1e5)
    delta={name:LOOP*dt*term for name,term in channels.items()}
    delta_full=sum(delta.values(),np.zeros_like(c0))
    delta_ward=max(relative_operator_ward(delta_full,tf,ts)
                   for _,gens in sectors for tf,ts in gens)
    quartic_ward=0.0
    for _,gens in sectors:
        for _,ts in gens:
            v=(np.einsum('ea,ebcd->abcd',ts,quartic,optimize=True)
              +np.einsum('eb,aecd->abcd',ts,quartic,optimize=True)
              +np.einsum('ec,abed->abcd',ts,quartic,optimize=True)
              +np.einsum('ed,abce->abcd',ts,quartic,optimize=True))
            quartic_ward=max(quartic_ward,
                float(np.max(abs(v))/max(np.max(abs(quartic)),1e-100)))
    c_s=c0+delta_full
    b5={name:projection(term,full_fermions) for name,term in delta.items()}
    candidate_log=projection(delta_full,full_fermions)
    require(_relative(candidate_log,sum(b5.values(),np.zeros((3,3),complex)))<1e-10,
            'Seven-channel projection failed recomposition')
    native_log=decode_complex(b['source_EFT1_scale_partition_log_GeVinv'])
    native_c5=decode_complex(b['source_finite_sector_candidate_GeVinv'])
    candidate=replace_eft1_log(native_c5,native_log,candidate_log,master_to_native)
    ratio=_relative(master_to_native*candidate_log,native_log)
    result={
      'format':FORMAT,'model':model,
      'status':'FullTensorSevenChannelStrictOneLoopConditional_NotPhysicalC5',
      'hierarchy':[['UV_1e7','all_F_1e5'],['EFT1_S1_S2','both_S_1e3'],['SMEFT_100']],
      'native_B11_source_C5_GeVinv':encode_complex(native_c5),
      'native_B11_pole_derived_EFT1_log_GeVinv':encode_complex(native_log),
      'master_tree_C5_at_muF_GeVinv':encode_complex(projection(c0,full_fermions)),
      'master_candidate_EFT1_log_GeVinv':encode_complex(candidate_log),
      'master_to_native_factor_assumed':float(master_to_native),
      'master_vs_native_log_relative_difference':ratio,
      'replacement_conditional_C5_at_muS_GeVinv':encode_complex(candidate),
      'candidate_ADM':beta_mode,'candidate_ADM_backend':backend,
      'candidate_channel_C5_delta_GeVinv':{k:encode_complex(v) for k,v in b5.items()},
      'source_tree_F_basis_info':source_info,
      'source_CG_derived_triplet_basis':basis_info,
      'SU2_U1_SU3_Ward':{'tree_max':tree_ward,'candidate_ADM_delta_max':delta_ward,
                        'scalar_potential_max':quartic_ward,
                        'all_three_independently_pass':max(delta_ward,quartic_ward,tree_ward)<1e-8},
      'matchete_source_scalar_quartic':quartic_info,
      'D2_extra_reference':d2_info,
      'input_sha256':{'batch9':digest(b9),'batch11':digest(b11),
                      'source_scalar_quartic':digest(snap.quartic_seed)},
      'strict_order':'C(muS)=C_tree(muF)+ln(muS/muF)/(16pi2)*beta_candidate[C_tree,g(muF)]',
      'hard_F_one_loop_re_evolved':False,
      'scalar_bubble_readded':False,
      'native_pole_log_readded':False,
      'full_master_pair_symmetric_coordinate_size':int(c0.size),
      'native_GammaCC_to_Weyl_factor_certified':False,
      'full_EFT1_operator_closure_certified':False,
      'D2_all_extra_scalar_only_quartics_included':True,
      'D2_extra_scalar_only_RGBeta_to_Matchete_map_certified':False,
      'physical_C5_claimed':False,
      'fitting_repeated':False}
    if include_smeft:
        from Numerical.parallel_t3.complete.two_model_conditional_uvfit_batch12 import candidate_c5_as_smeft
        from Numerical.parallel_t3.complete.full_general_completion import finish_smeft
        from Numerical.parallel_t3.complete.two_model_conditional_uvfit_batch12 import BENCHMARK
        b1=root/Path('output/parallel_t3/full_general_batch1/thesis_full_general_batch1_v1')/model
        b2=root/Path('output/parallel_t3/full_general_batch2/thesis_full_general_batch2_raw_v1')/model
        laurent=root/Path('output/parallel_t3/full_general_batch2/thesis_full_general_batch2_laurent_v2')/model
        evidence=existing_evidence(model,b1,b2,laurent)
        cfg=json.loads((root/BENCHMARK).read_text(encoding='utf-8-sig'))
        target=json.loads((root/cfg['fitting']['target']).read_text(encoding='utf-8-sig'))
        out=finish_smeft(evidence,candidate_c5_as_smeft(candidate,signed_factor),
                         vev=float(cfg['vev_GeV']),target=target)
        result['conditional_100GeV_neutrinos']=out
        result['signed_factor_assumed']=float(signed_factor)
    if output is not None:
        output=Path(output)
        output.mkdir(parents=True,exist_ok=True)
        (output/'result.json').write_text(json.dumps(result,indent=2,
                                       sort_keys=True,allow_nan=False,default=str)+'\n',
                                       encoding='utf-8')
        np.savez_compressed(output/'full_tensor_candidate.npz',
                            C_tree_muF=c0,C_candidate_muS=c_s,
                            delta_candidate_total=delta_full,
                            delta_EFT1_log_C5=candidate_log,
                            **{'delta_'+k:v for k,v in delta.items()})
    return result
