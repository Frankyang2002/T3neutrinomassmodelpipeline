"""All-16 source-guarded, *conditional* first-threshold Weyl LLSS amplitudes.

Strict ordering: UV -> integrate ALL heavy fermions at muF -> EFT1 containing
SM fermions and BOTH BSM scalars -> integrate both scalars at muS -> SMEFT.
This is the finite-dimensional tree Gaussian Schur complement ONLY; it omits
hard one-loop matching and is NOT a physical Matchete/RGBeta certificate.

Existing direct_one_loop convention:
    L = -1/2 psi^T M psi - nu_a psi_i (Y[a,i,s] phi_s
        + Z[a,i,s] phi_s^*) + h.c.
    L_eff,tree = +1/2 J^T M^{-1} J,  J_i=sum_a,s nu_a(V^+ phi_s+V^- phi_s^*)
with left-handed spinor contraction implicit. K is the *symmetric* second
source derivative of L_eff; the Lagrangian contains K/2, not K itself.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any
import numpy as np
import sympy as sp

from Numerical.direct_one_loop.core import _complex, encoded, takagi
from .certification import CLASSES, all_models, charge_bases, validate_components, digest, load
from .spinor_mass_dictionary import guard_sources, mass_ward, optional_evidence

FORMAT='t3-direct-uv-first-threshold-gaussian-llss-conditional-v1'
REQUIRED_SOURCE_STATUS='source_mass_ward_plus_declared_weyl_gamma_basis_not_native_certified'
CHANNELS=(('field','field','YY'),('field','conjugate','YZ'),
          ('conjugate','field','ZY'),('conjugate','conjugate','ZZ'))


class AmplitudeError(ValueError):
    pass


def require(condition:bool, why:str)->None:
    if not condition:
        raise AmplitudeError(why)


def _add(*items:dict)->dict:
    out={}
    for b in items:
        for mon,c in b.items():
            out[mon]=sp.expand(out.get(mon,sp.Integer(0))+c)
    return {k:v for k,v in out.items() if sp.simplify(v)!=0}


# Exactly anticommuting spinor variables, with one fixed monomial basis.
_NAMES=('l0','l1','ch0','ch1','xi0','xi1',
        'lb0','lb1','cb0','cb1','xb0','xb1')
_POS={x:i for i,x in enumerate(_NAMES)}
_BAR={_NAMES[i]:_NAMES[i+6] for i in range(6)}
_BAR.update({_NAMES[i+6]:_NAMES[i] for i in range(6)})


def exterior(a:dict,b:dict)->dict:
    out={}
    for ma,ca in a.items():
        for mb,cb in b.items():
            m=ma+mb
            if len(set(m))!=len(m):continue
            swaps=sum(_POS[m[i]]>_POS[m[j]] for i in range(len(m)) for j in range(i+1,len(m)))
            mon=tuple(sorted(m,key=_POS.__getitem__))
            out[mon]=sp.expand(out.get(mon,sp.Integer(0))+(-1)**swaps*ca*cb)
    return {k:v for k,v in out.items() if sp.simplify(v)!=0}


def _v(name:str)->dict:
    require(name in _POS,'Unknown Grassmann component: '+name)
    return {(name,):sp.Integer(1)}


def _scale(a:dict,k:Any)->dict:
    return {m:sp.expand(c*k) for m,c in a.items() if sp.simplify(c*k)!=0}


def dagger(a:dict)->dict:
    out={}
    for mon,coeff in a.items():
        term={():sp.conjugate(coeff)}
        for var in reversed(mon):term=exterior(term,_v(_BAR[var]))
        out=_add(out,term)
    return out


def _dot(a:list[dict],b:list[dict])->dict:
    require(len(a)==len(b),'Incompatible spinor lengths')
    out={}
    for x,y in zip(a,b):out=_add(out,exterior(x,y))
    return out


def _mm(m:sp.Matrix, v:list[dict])->list[dict]:
    return [_add(*[_scale(v[j],m[i,j]) for j in range(len(v))])
            for i in range(len(v))]


def explicit_four_spinor_witness()->dict[str,Any]:
    """One declared 4x4 C and index basis; no native package phase claim.

    F=(chi_alpha, -epsilon_{dot a dot b} bar xi_dot b),
    bar F=( -xi_1, xi_0, bar chi_0, bar chi_1).
    The index-raising convention for bar xi is made explicit and is essential
    to avoid assigning a spurious y2 sign. C=+i gamma^2 gamma^0.
    """
    z=sp.zeros(2); eye=sp.eye(2);eps=sp.Matrix([[0,1],[-1,0]])
    def block(a,b,c,d):return sp.BlockMatrix([[a,b],[c,d]]).as_explicit()
    g0=block(z,eye,eye,z)
    s2=sp.Matrix([[0,-sp.I],[sp.I,0]])
    g2=block(z,s2,-s2,z)
    C=sp.simplify(sp.I*g2*g0)
    require(C==block(eps,z,z,-eps),'Wrong explicitly chosen C orientation')
    require(C.T==-C and C*C==-sp.eye(4),'C must be antisymmetric and square to -1')
    barL=[{}, {},_v('lb0'),_v('lb1')]
    F=[_v('ch0'),_v('ch1'),_scale(_v('xb1'),-1),_v('xb0')]
    barF=[_scale(_v('xi1'),-1),_v('xi0'),_v('cb0'),_v('cb1')]
    from_actual_dirac_dagger=_mm(g0,[dagger(x) for x in F])
    require(barF==from_actual_dirac_dagger, 'bar F must equal F-dagger gamma0')
    L=[_v('l0'),_v('l1'),{},{}]
    require(barL==_mm(g0,[dagger(x) for x in L]), 'bar L must equal L-dagger gamma0')
    ccF=_mm(C,barF)
    expectedCC=[_v('xi0'),_v('xi1'),_scale(_v('cb1'),-1),_v('cb0')]
    require(ccF==expectedCC,'Chosen CConj F does not exchange chi and xi')
    W1=_add(exterior(_v('l0'),_v('ch1')),_scale(exterior(_v('l1'),_v('ch0')),-1))
    W2=_add(exterior(_v('l0'),_v('xi1')),_scale(exterior(_v('l1'),_v('xi0')),-1))
    WM=_add(exterior(_v('ch0'),_v('xi1')),_scale(exterior(_v('ch1'),_v('xi0')),-1))
    y1=_dot(barL,ccF);y2=_dot(barL,F); mass=_dot(barF,F)
    require(y1==dagger(W1),'y1 four-spinor/Weyl identity failed')
    require(y2==dagger(W2),'y2 four-spinor/Weyl identity failed')
    require(mass==_add(WM,dagger(WM)),'Dirac four-spinor/Weyl mass identity failed')
    # Distinguish the physical H.c. of each original Matchete vertex.
    require(y1!=dagger(W2) and y2!=dagger(W1),'Two Dirac Weyl chiralities collapsed')
    # A self-conjugate Majorana spinor is obtained by xi=chi consistently.
    mp={'xi0':'ch0','xi1':'ch1','xb0':'cb0','xb1':'cb1'}
    def subst(poly:dict)->dict:
        r={}
        for mon,c in poly.items():
            term={():c}
            for x in mon:term=exterior(term,_v(mp.get(x,x)))
            r=_add(r,term)
        return r
    majorana_m=subst(mass)
    Wmaj=subst(WM)
    require(majorana_m==_add(Wmaj,dagger(Wmaj)), 'Majorana mass bilinear failed')
    require(subst(y1)==subst(y2),'Majorana y1/y2 use common heavy Weyl spinor in chosen basis')
    require([subst(x) for x in ccF]==[subst(x) for x in F],
            'Majorana CConj F must coincide with F in this chosen convention')
    require(Wmaj==_scale(exterior(_v('ch0'),_v('ch1')),2),
            'Majorana identical-Weyl factor of two wrong')
    # Ensure a false C sign cannot pass the y1 identity.
    require(_dot(barL,_mm(-C,barF))!=dagger(W1),'Wrong-C negative control failed')
    return {
        'status':'explicit_chosen_four_spinor_dirac_majorana_action_equalities_pass',
        'conventions':{
            'gamma0':'[[0,I],[I,0]]','gamma2':'[[0,sigma2],[-sigma2,0]]',
            'C':'i gamma2 gamma0 = diag(epsilon,-epsilon)',
            'epsilon':'[[0,1],[-1,0]]',
            'F_dirac':'(chi_alpha, -epsilon bar(xi)_dot_alpha)',
            'F_charge_conjugate':'(xi_alpha, -epsilon bar(chi)_dot_alpha)',
            'F_majorana':'xi=chi',
        },
        'exact_identities':{
            'BarL_CConjF':'dagger(l epsilon chi)',
            'BarL_F':'dagger(l epsilon xi)',
            'BarF_F':'chi epsilon xi + h.c.',
            'majorana_half_mass':'-MF/2 (chi epsilon chi + h.c.)',
            'yukawa_HC_scalar_orientations':{'y1':'S1^dagger','y2':'S2'},
        },
        'negative_wrong_C_sign_rejected':True,
        'majorana_repeated_grassmann_spinor_factor_two_checked':True,
        'absolute_installed_Matchete_GammaCC_phase_measured':False,
        'native_RGBeta_epsilon_numeric_orientation_measured':False,
        'relative_Matchete_RGBeta_signed_mass_and_yukawa_dictionary_certified':False,
        'candidate_numerical_YZ_from_original_topology_HC':
            'Y1_HC_conjugate and Y2_HC_field coefficients = -conjugate(original Matchete y_i * CG_i) ONLY in the chosen Weyl and scalar basis',
    }


def phase_invariant_source_orientation(majorana:bool)->dict[str,Any]:
    """Cross-check phases in the four-spinor *original-H.c.* orientation.

    h1=conj(y1_M), h2=conj(y2_M), holomorphic mass MF and mixing lambda.
    h1 nu chi S1dag + h2 nu xi S2 + MF chi xi + lambda S1 S2dag.
    """
    charge={
        'h1':(1,0,-1,0), 'h2':(0,1,0,1),
        'MF':(1,1,0,0),'lambda':(0,0,1,-1),
    }
    if majorana:
        charge={'h1':(1,0,-1,0),'h2':(1,0,0,1),
                'MF':(2,0,0,0),'lambda':(0,0,1,-1)}
    def total(powers):return [sum(p*charge[k][i] for k,p in powers.items()) for i in range(4)]
    correct=total({'h1':1,'h2':1,'MF':-1,'lambda':1})
    incorrect=total({'h1':1,'h2':1,'MF':-1,'lambda':-1})
    require(correct==[0]*4 and incorrect!=[0]*4,'Source-orientation spurion covariance failed')
    return {
        'holomorphic_hc_of_native_topology':'h1 nu chi S1dag + h2 nu xi S2 + MF chi xi + lambda HH S1 S2dag',
        'h1_is_candidate_conjugate_native_y1':True,
        'h2_is_candidate_conjugate_native_y2':True,
        'invariant_source_HC_spurion':'h1*h2*conjugate(MF)*lambda',
        'invariant_in_unbarred_Matchete_names':'conjugate(y1_M)*conjugate(y2_M)*conjugate(MF)*lambdaT3_M',
        'sum_charge':correct,'wrong_lambda_conjugation_charge':incorrect,
        'is_a_conditional_common_chirality_and_scalar_basis_statement':True,
        'absolute_physical_native_couplings_certified':False,
    }


def gaussian_tensors(data:dict[str,Any],kind:str,alpha:int,*,include_tensors:bool=True)->dict:
    """All bosonic scalar orientations + every source pair, no truncation.

    External indices are (lepton_a,lepton_b,scalar_s,scalar_t), ordered
    scalar tags are {field,conjugate} squared. The full coefficient K is
    in L_eff=K/2: even for distinct source vertices no extra factor two.
    """
    structural=validate_components(data,kind,alpha)
    M=_complex(data['fermion_mass_gev'],'M'); n=len(M)
    vals=np.linalg.eigvalsh(M.conj().T@M)
    require(vals[0]>max(vals[-1],1.)*1e-26,'Singular heavy fermion mass; no local LLSS tree expansion')
    identity=np.eye(n,dtype=complex)
    Minv=np.linalg.solve(M,identity)
    inv_err=float(np.linalg.norm(M@Minv-identity)/max(np.linalg.norm(identity),1))
    require(inv_err<2e-10,'Heavy Gaussian inverse inaccurate')
    U,masses,takagi_err=takagi(M)
    require(np.min(masses)>0,'Heavy zero Takagi mass')
    second=U@np.diag(1./masses)@U.T
    cross_err=float(np.linalg.norm(Minv-second)/max(np.linalg.norm(Minv),1e-100))
    require(cross_err<3e-9,'Takagi inverse disagrees with direct inverse')
    # The second derivative of the eliminated Gaussian action is M^{-1}.
    rng=np.random.default_rng(20261103)
    J=rng.normal(size=n)+1j*rng.normal(size=n)
    psi=-Minv@J
    action=(-.5*psi.T@M@psi-psi.T@J)
    effective=.5*J.T@Minv@J
    action_err=float(abs(action-effective)/max(abs(effective),1e-100))
    require(action_err<2e-9,'Gaussian EOM/substitution identity failed')
    ns=sum(CLASSES[kind][:2]); vertices=data['vertices']; names=[v['name'] for v in vertices]
    require(len(names)==len(set(names)),'Duplicate source vertex names')
    vfields=[]
    for v in vertices:
        vfields.append({tag:_complex(v[tag],f"{v['name']}.{tag}",(3,n,ns))
                        if tag in v else np.zeros((3,n,ns),complex)
                        for tag in ('field','conjugate')})
    qb=charge_bases(kind,alpha)
    qscalar=qb['scalar']
    records={}; maxviol=0.; maxsym=0.; maxsource=0.
    terms_all={name:np.zeros((3,3,ns,ns),complex) for _,_,name in CHANNELS}
    for tag1,tag2,ch in CHANNELS:
        pieces={}
        for a in range(len(vertices)):
            for b in range(a,len(vertices)):
                va=vfields[a][tag1]; vb=vfields[b][tag2]
                amp=np.einsum('ais,ij,bjt->abst',va,Minv,vb,optimize=True)
                if a!=b:
                    amp+=np.einsum('ais,ij,bjt->abst',vfields[b][tag1],Minv,vfields[a][tag2],optimize=True)
                # Independent Takagi (heavy eigenstate) projection, no fermion dropped.
                ua=np.einsum('ais,ik->aks',va,U,optimize=True)
                ub=np.einsum('ais,ik->aks',vb,U,optimize=True)
                spectral=np.einsum('aks,k,bkt->abst',ua,1./masses,ub,optimize=True)
                if a!=b:
                    spectral+=np.einsum('aks,k,bkt->abst',np.einsum('ais,ik->aks',vfields[b][tag1],U,optimize=True),1./masses,
                                        np.einsum('ais,ik->aks',vfields[a][tag2],U,optimize=True),optimize=True)
                # The sum may vanish by symmetry/cancellation. Do NOT divide a
                # cancellation-level rounding residue by the vanishing sum.
                natural=(np.linalg.norm(va)*np.linalg.norm(vb)+
                         (np.linalg.norm(vfields[b][tag1])*np.linalg.norm(vfields[a][tag2])
                          if a!=b else 0.))*np.linalg.norm(Minv)
                mismatch=float(np.linalg.norm(amp-spectral)/max(np.linalg.norm(amp),
                                 np.linalg.norm(spectral),natural,1e-100))
                require(mismatch<5e-9,'Heavy-eigenstate source term mismatch: '+str(names[a])+' '+str(names[b]))
                maxsource=max(maxsource,mismatch)
                pieces[f'{names[a]} x {names[b]}']=encoded(amp) if include_tensors else {
                    'norm_gev_minus1':float(np.linalg.norm(amp))}
                terms_all[ch]+=amp
        records[ch]={'scalar_tags':[tag1,tag2],
                     'source_pair_tensors_gev_minus1':pieces,
                     'full_tensor_gev_minus1':encoded(terms_all[ch]) if include_tensors else None,
                     'frobenius_gev_minus1':float(np.linalg.norm(terms_all[ch])),
                     'maximum_absolute_gev_minus1':float(np.max(np.abs(terms_all[ch]))) }
    for tag1,tag2,ch in CHANNELS:
        twin=next(name for a,b,name in CHANNELS if (a,b)==(tag2,tag1))
        mat=terms_all[ch]; oth=terms_all[twin].transpose(1,0,3,2)
        err=float(np.linalg.norm(mat-oth)/max(np.linalg.norm(mat),np.linalg.norm(oth),1e-100))
        require(err<2e-9,'LLSS external-leg exchange symmetry violated')
        maxsym=max(maxsym,err)
        sig1=1 if tag1=='field' else -1
        sig2=1 if tag2=='field' else -1
        allowed=np.isclose(sig1*qscalar[:,None]+sig2*qscalar[None,:],0.)
        forbidden=mat*np.logical_not(allowed)[None,None,:,:]
        err=float(np.linalg.norm(forbidden)/max(np.linalg.norm(mat),1e-100))
        require(err<2e-9,'LLSS tensor violates electric charge conservation')
        maxviol=max(maxviol,err)
        records[ch]['external_exchange_conjugate_tag']=twin
        records[ch]['forbidden_charge_tensor_relative_norm']=err
    # Cross-check source-pair sum is exact even under weak cancellations.
    norms={k:float(np.linalg.norm(v)) for k,v in terms_all.items()}
    return {
        'status':'conditional_all_heavy_fermion_tree_LLSS_reconstructed',
        'input_model':data.get('model'),
        'first_threshold':{'formal_algebra':'ALL BSM fermion eigenstates eliminated together at input scale',
                           'retained':'BOTH BSM scalars, all SM leptons and other SM fields',
                           'matched_at_physical_muF':False,
                           'input_scale_gev':float(data['renormalization_scale_gev']),
                           'second_threshold':'BOTH BSM scalars at mu_S; NOT executed by this calculation'},
        'source_vertices':names,'heavy_fermion_state_count':n,'scalar_component_count':ns,
        'heavy_takagi_masses_gev':masses.tolist(),
        'heavy_inverse_direct_vs_takagi_relative_residual':cross_err,
        'heavy_takagi_residual':takagi_err,
        'gaussian_equations_of_motion_relative_residual':action_err,
        'source_pair_heavy_eigenstate_reconstruction_max_error':maxsource,
        'LLSS_exchange_symmetry_max_relative_error':maxsym,
        'LLSS_charge_violation_max_relative_error':maxviol,
        'source_component_charge_and_hessian_status':structural['status'],
        'tensor_channels':records,
        'channel_norms_gev_minus1':norms,
        'total_charge_allowed_coefficient_norm_gev_minus1':float(np.sqrt(sum(v*v for v in norms.values()))),
        'spinor_and_coupling_dictionary_physically_certified':False,
        'fermion_threshold_hard_one_loop_LLHH_included':False,
        'scalar_threshold_matching_or_running_performed':False,
        'physical_C5_in_SMEFT_computed':False,
    }


def _maybe(path:Path)->dict|None:
    return load(path) if path.is_file() else None


def _eft1_boundary_gate(root:Path,short:str,uv_mu:float)->dict[str,Any]:
    """Metadata-only EFT boundary gate; NEVER pretend different scales align."""
    tag=short.replace('_-','_m').replace('_+','_p')
    pattern='output/general_observable/lo_fermion_boundary/releases/source_*/'+tag+'/boundary.json'
    candidates=sorted(root.glob(pattern))
    if not candidates:
        return {'status':'no_source_fermion_threshold_release',
                'compared_LLSS_tensors':False}
    if len(candidates)>1:
        return {'status':'ambiguous_multiple_fermion_boundary_releases',
                'releases':[str(p.relative_to(root)) for p in candidates],
                'compared_LLSS_tensors':False}
    p=candidates[0]
    rec=load(p)
    require(rec.get('model')==tag,'F threshold source release model label disagrees')
    require(rec.get('field_content',{}).get('heavy_fermions_integrated_out')=='all at mu_F',
            'F threshold release incorrectly integrates a partial heavy multiplet')
    require(rec.get('field_content',{}).get('heavy_scalar_multiplets')=='both retained in EFT1',
            'F threshold release scalar ordering changed')
    mu=float(rec['mu_F_gev'])
    spin=rec.get('spinor_normalization',{})
    norm=spin.get('certified') is True
    scale=np.isclose(mu,uv_mu,rtol=1e-12)
    return {'status':'NOT_COMPARED_SCALE_MISMATCH_AND_SPINOR_NORMALIZATION_OPEN'
            if not scale else 'NOT_COMPARED_PHYSICAL_SPINOR_AND_CG_DICTIONARY_OPEN',
            'release':str(p.relative_to(root)),
            'uv_input_mu_gev':uv_mu,'eft1_threshold_mu_F_gev':mu,
            'same_renormalization_scale':bool(scale),
            'native_spinor_normalization_certified_in_release':norm,
            'compared_LLSS_tensors':False,
            'hard_fermion_matching_included':False}


def batch(root:Path,input_dir:Path|None=None,*,include_tensors:bool=True)->dict:
    """Must guard local UV source and frozen benchmark data; never infer fixes."""
    root=root.resolve()
    source=guard_sources(root)
    spinor=explicit_four_spinor_witness()
    native=optional_evidence(root)
    if input_dir is None:
        input_dir=root/'output/general_observable/lo_direct_one_loop/provisional_inputs'
    rows=[]
    for short,key in all_models():
        bench_path=root/'Numerical/benchmarks/fitted_reference'/(key+'.json')
        input_path=input_dir/(short+'.json')
        modelrow={'model':short,'benchmark_path':str(bench_path.relative_to(root)),
                  'input_path':str(input_path),'physically_certified':False}
        try:
            require(bench_path.is_file() and input_path.is_file(),
                    'Missing frozen benchmark or existing source candidate component input')
            b=load(bench_path);inp=load(input_path)
            require(inp.get('model')==short,'Input/model mismatch')
            require(inp.get('provenance',{}).get('benchmark_sha256')==digest(bench_path),
                    'Frozen benchmark SHA256 provenance mismatch')
            rep=b['representation'];kind,alpha=short[0],int(short[2:])
            d1,d2,df,_=CLASSES[kind]
            require(tuple(rep[k] for k in ('d_s1','d_s2','d_f','alpha'))==(d1,d2,df,alpha),
                    'Frozen benchmark representation mismatch')
            require(np.isclose(float(inp['renormalization_scale_gev']),
                               float(b['base_state']['mu_gev']),rtol=1e-12),
                    'Frozen UV input renormalization scale mismatch')
            modelrow['status']='conditional_tree_LLSS_completed'
            modelrow['fermion_type']='Majorana' if charge_bases(kind,alpha)['self_conjugate'] else 'Dirac'
            modelrow['source_F_mass_ward_passed']=mass_ward(df)['all_su2_generators_ward_zero']
            modelrow['chosen_4spinor_weyl_phases']=phase_invariant_source_orientation(
                modelrow['fermion_type']=='Majorana')
            modelrow['frozen_benchmark_sha256']=digest(bench_path)
            modelrow['output']=gaussian_tensors(inp,kind,alpha,include_tensors=include_tensors)
            modelrow['eft1_release_gate']=_eft1_boundary_gate(root,short,
                float(inp['renormalization_scale_gev']))
        except (AmplitudeError,ValueError,KeyError,TypeError,OSError,np.linalg.LinAlgError,RuntimeError) as exc:
            modelrow['status']='blocked'
            modelrow['reason']=str(exc)
        rows.append(modelrow)
    return {
        'format':FORMAT,'status':'conditional_LLSS_source_candidate_only',
        'source_contract_status':REQUIRED_SOURCE_STATUS,
        'local_source_sha256':source['source_sha256'],
        'chosen_four_spinor_basis':spinor,
        'optional_previous_native_evidence':native,
        'tree_amplitude_computed_models':sum(r['status']=='conditional_tree_LLSS_completed' for r in rows),
        'blocked_models':sum(r['status']=='blocked' for r in rows),
        'physical_certified_models':0,
        'EFT_order':['UV','integrate all BSM heavy fermions together at mu_F',
                     'EFT1 retains both BSM scalars and all SM fermions',
                     'integrate both BSM scalars together at mu_S','SMEFT'],
        'not_included':['raw fermion-threshold hard one-loop LLHH matching',
                        'scalar-loop matching at mu_S','EFT1 evolution',
                        'physical Matchete/RGBeta source-to-Weyl phase/CG normalization',
                        'optional extra Yukawas not present in frozen restricted component inputs'],
        'models':rows,
    }
