"""Batch 14: source-based physical closure investigation, A(-4) and D(-2).

This audits THREE independent requirements together, using existing local outputs:
(1) an independent Matchete GammaCC -> physical Weinberg operator anchor;
(2) completeness and finite subtraction of the FULL F/soft matching source;
(3) gauge covariance, seven-term ADM and scalar coupling basis maps.

A passing internal Laurent or Ward test is NOT a certificate of those
independent requirements. No uncertified matrix is renamed physical C5.
Only existing files are read; this module produces a new, separate output.
"""
from __future__ import annotations

from collections import Counter
from hashlib import sha256
from pathlib import Path
import json
import math
import re
from typing import Any

import numpy as np

MODELS=('A_-4','D_-2')
FORMAT='t3-a4-d2-physical-three-gate-closure-audit-batch14-v1'
B1=Path('output/parallel_t3/full_general_batch1/thesis_full_general_batch1_v1')
B2=Path('output/parallel_t3/full_general_batch2/thesis_full_general_batch2_raw_v1')
B9=Path('output/parallel_t3/thesis_full_general_joint_scalar_batch9_v1')
B11=Path('output/parallel_t3/thesis_native_laurent_batch11_v1')
B13=Path('output/parallel_t3/thesis_master_wilson_batch13_cgfix_v3')
SOFT={
 'A_-4':Path('output/parallel_t3/full_general_final_one_batch/thesis_a4_fresh_f_v3/A_-4/A_scalar_soft_only_unrenormalized.txt'),
 'D_-2':Path('output/parallel_t3/full_general_final_one_batch/thesis_physical_scalar_continuation_v1/D_-2/D_scalar_soft_only_unrenormalized.txt'),
}


def insist(condition:bool,message:str)->None:
    if not condition: raise ValueError('Batch14: '+message)


def read_obj(p:Path)->dict:
    insist(p.is_file(),'missing required saved file '+str(p))
    obj=json.loads(p.read_text(encoding='utf-8-sig'))
    insist(isinstance(obj,dict),'expected JSON object '+str(p))
    return obj


def hash_file(p:Path)->str:
    insist(p.is_file(),'missing '+str(p))
    h=sha256()
    with p.open('rb') as stream:
        for b in iter(lambda:stream.read(1<<20),b''):h.update(b)
    return h.hexdigest()


def decoded_matrix(v:Any)->np.ndarray:
    if isinstance(v,dict) and set(v)=={'real','imag'}:
        z=np.asarray(v['real'],dtype=float)+1j*np.asarray(v['imag'],dtype=float)
    else:
        z=np.asarray(v,dtype=complex)
    insist(z.shape==(3,3) and np.isfinite(z).all(),'invalid complex 3x3 source matrix')
    return z


def relative(a:np.ndarray,b:np.ndarray)->float:
    a=np.asarray(a,complex);b=np.asarray(b,complex)
    return float(np.linalg.norm(a-b)/max(np.linalg.norm(a),np.linalg.norm(b),1e-100))


def inferred_real_normalization(master:np.ndarray,native:np.ndarray)->dict:
    """Projection of native onto master; no physical sign inferred.

    A common rescaling of the external LNV operator changes both sides, so
    this measurement CANNOT certify the GammaCC/Weyl sign or absolute factor.
    """
    m=np.asarray(master,complex).reshape(-1)
    n=np.asarray(native,complex).reshape(-1)
    insist(m.shape==n.shape and m.size>0 and np.isfinite(m).all()
           and np.isfinite(n).all(),'invalid log arrays')
    denom=float(np.vdot(m,m).real)
    if denom <= 1e-100:
        return {'status':'Unidentifiable_ZeroMasterLog','factor':None,'residual':None}
    coefficient=float(np.vdot(m,n).real/denom)
    residual=relative(coefficient*m,n)
    return {'status':'EmpiricalRelativeLogFactorOnly_NotGammaCCCertificate',
            'factor':coefficient,'residual':residual,
            'imaginary_projection_ratio':float(np.vdot(m,n).imag/denom),
            'absolute_gamma_CC_to_Weyl_factor_certified':False}


def audit_normalization(ledger:dict, b13:dict)->dict:
    """Verify native H/Hdag spinor structure, never assign its physical sign."""
    ll=[x for x in ledger['terms'] if x['sector']=='CandidateLeptonNumberViolatingLLHH']
    chiral=Counter()
    for row in ll:
        t=row['raw_inputform']
        insist(t.count('GammaCC')==1 and t.count('NCM[')==1,
               'changed native GammaCC LLHH field structure')
        assert 'Proj[1]' in t or 'Proj[-1]' in t
        p='PR' if 'Proj[1]' in t and 'Proj[-1]' not in t else 'PL' if 'Proj[-1]' in t and 'Proj[1]' not in t else None
        insist(p is not None,'ambiguous chiral LLHH native source')
        chiral[p]+=1
    expected=3 if b13['model']=='A_-4' else 6
    insist(chiral==Counter({'PL':expected,'PR':expected}),
           'native LLHH chirality or operator count changed')
    fit=inferred_real_normalization(
        decoded_matrix(b13['master_candidate_EFT1_log_GeVinv']),
        decoded_matrix(b13['native_B11_pole_derived_EFT1_log_GeVinv']))
    return {
        'native_GammaCC_operator_chiralities':dict(chiral),
        'relative_native_master_log_coefficient':fit,
        'absolute_native_GammaCC_to_physical_Weyl_operator_anchored':False,
        'independent_explicit_spinor_amplitude_available':False,
        'independent_source_operator_sign_and_symmetry_factor_fixed':False,
        'reason':'Relative logs and pole residues are invariant under a COMMON operator rescaling. '
                 'An independent normalized on-shell spinor amplitude or a verified '
                 'Matchete basis-operator mapping is required.',
        'gate_passed':False,
    }


def audit_matching(ledger:dict,b11:dict,soft_status:dict,model:str,soft_text:str)->dict:
    """Full archived F operator inventory; distinguish isolated LLHH from all operators."""
    expected={'A_-4':{'OtherEFT1Operator':81,
                       'CandidateLeptonNumberViolatingLLScalarScalar':48,
                       'CandidateLeptonNumberViolatingLLHH':6},
              'D_-2':{'OtherEFT1Operator':836,
                       'CandidateLeptonNumberViolatingLLScalarScalar':352,
                       'CandidateLeptonNumberViolatingLLHH':12}}[model]
    insist(ledger['status']=='UnrenormalisedOneLoopTermsEnumerated'
           and len(ledger['terms'])==sum(expected.values()),'raw F source inventory changed')
    counts=dict(Counter(x['sector'] for x in ledger['terms']))
    insist(counts==expected and ledger['sector_counts']==expected,
           'one-loop sector inventory incomplete or changed')
    pole=decoded_matrix(b11['source_hard_symmetric_pole_GeVinv'])+decoded_matrix(
         b11['source_soft_symmetric_pole_GeVinv'])
    hard=decoded_matrix(b11['source_hard_symmetric_pole_GeVinv'])
    soft=decoded_matrix(b11['source_soft_symmetric_pole_GeVinv'])
    cancellation=float(np.linalg.norm(pole)/max(np.linalg.norm(hard),np.linalg.norm(soft),1e-100))
    insist(cancellation<2e-7,'saved native hard/soft source-sector poles do not cancel')
    match_txt=soft_text.replace(' ','').replace('\n','')
    for token in ('hbar','ϵ','μbar2','T3EFT1W12'):
        insist(token in match_txt,'soft source missing '+token)
    target=10 if model=='A_-4' else 68
    nterms=soft_status.get('A_physical_scalar_soft_Weinberg_terms',
                           soft_status.get('physical_scalar_soft_Weinberg_terms'))
    insist(nterms==target,'native two-scalar soft term inventory changed')
    native_log=decoded_matrix(b11['source_EFT1_scale_partition_log_GeVinv'])
    source_fin=decoded_matrix(b11['source_finite_sector_candidate_GeVinv'])
    hard_reg=decoded_matrix(b11['source_hard_symmetric_regular_GeVinv'])
    soft_reg=decoded_matrix(b11['source_soft_symmetric_regular_GeVinv'])
    insist(relative(hard_reg+soft_reg+native_log,source_fin)<1e-10,
           'native hard + soft + EFT1 log recomposition failed')
    return {
       'full_F_one_loop_operator_count':len(ledger['terms']),
       'full_F_sector_counts':counts,
       'F_LNV_LLSS_one_loop_operator_count_not_dropped_from_ledger':expected['CandidateLeptonNumberViolatingLLScalarScalar'],
       'native_simultaneous_scalar_soft_terms':target,
       'source_channel_hard_soft_pole_residual_relative':cancellation,
       'source_native_one_loop_log_identity_relative':float(b11['native_pole_log_identity_relative']),
       'source_scale_partition_invariance_relative':float(b11['source_scale_partition_invariance_relative']),
       'source_sector_finite_reconstruction_passed':True,
       'other_F_operators_matched_into_complete_EFT1_counterterm_basis':False,
       'finite_MSbar_subtractions_for_all_relevant_operators_certified':False,
       'evanescent_and_FieldRedefinition_scheme_complete':False,
       'gate_passed':False,
       'reason':'Hard/soft LLHH source-sector closure passes, but the full one-loop '
                'renormalized Wilson boundary including counterterms and all '
                'operator mixings is not independently established.',
    }


def operator_ward(c:np.ndarray,t:np.ndarray,theta:np.ndarray)->float:
    """Gauge variation t^T on both Weyl and theta^T on both real scalars."""
    c=np.asarray(c,complex)
    t=np.asarray(t,complex);theta=np.asarray(theta,complex)
    insist(c.ndim==4 and c.shape[0]==c.shape[1] and c.shape[2]==c.shape[3]
           and t.shape==(c.shape[0],)*2 and theta.shape==(c.shape[2],)*2,
           'Ward tensor/generator dimensions incompatible')
    # Algebraically equivalent to GeneralEFT1SevenTermFlow.gauge_ward_residual.
    d=(np.einsum('ki,kjab->ijab',t,c,optimize=True)
       +np.einsum('kj,ikab->ijab',t,c,optimize=True)
       +np.einsum('ca,ijcb->ijab',theta,c,optimize=True)
       +np.einsum('cb,ijac->ijab',theta,c,optimize=True))
    return float(np.max(abs(d))/max(float(np.max(abs(c))),1e-100))


def scalar_quartic_ward(lam:np.ndarray,theta:np.ndarray)->float:
    lam=np.asarray(lam,complex)
    theta=np.asarray(theta,complex)
    insist(lam.ndim==4 and len(set(lam.shape))==1 and
           theta.shape==(lam.shape[0],)*2,'quartic Ward shape invalid')
    d=(np.einsum('ea,ebcd->abcd',theta,lam,optimize=True)
       +np.einsum('eb,aecd->abcd',theta,lam,optimize=True)
       +np.einsum('ec,abed->abcd',theta,lam,optimize=True)
       +np.einsum('ed,abce->abcd',theta,lam,optimize=True))
    return float(np.max(abs(d))/max(float(np.max(abs(lam))),1e-100))


def yukawa_ward(y:np.ndarray,t:np.ndarray,theta:np.ndarray)->float:
    y=np.asarray(y,complex);t=np.asarray(t,complex);theta=np.asarray(theta,complex)
    insist(y.ndim==3 and y.shape[0]==y.shape[1] and
           t.shape==(y.shape[0],)*2 and theta.shape==(y.shape[2],)*2,
           'Yukawa Ward shape invalid')
    d=(np.einsum('ki,kja->ija',t,y,optimize=True)
       +np.einsum('kj,ika->ija',t,y,optimize=True)
       +np.einsum('ba,ijb->ija',theta,y,optimize=True))
    return float(np.max(abs(d))/max(float(np.max(abs(y))),1e-100))


def audit_gauge(root:Path,model:str,npz:Path,b13:dict,*,full_ward:bool=True)->dict:
    """Source-driven full tensor gauge Ward checks; never implies ADM proof."""
    with np.load(npz,allow_pickle=False) as rows:
        required={'C_tree_muF','delta_candidate_total'}
        insist(required<=set(rows.files),'Batch13 full tensor data missing')
        c=np.asarray(rows['C_tree_muF'],complex)
        delta=np.asarray(rows['delta_candidate_total'],complex)
        require_shape=(45,45,12,12)
        insist(c.shape==delta.shape==require_shape,'Batch13 full tensor shape changed')
        channels={k:np.asarray(rows['delta_'+k],complex)
                  for k in ('yukawa_wavefunction','scalar_pair','mixed_yukawa_gauge',
                            'crossed_yukawa','conjugate_coefficient_yukawa',
                            'scalar_anomalous_dimension','fermion_anomalous_dimension')}
    insist(relative(sum(channels.values(),np.zeros_like(delta)),delta)<1e-10,
           'Batch13 seven-term beta failed recomposition')
    from Numerical.parallel_t3.complete.two_model_joint_scalar_batch9 import load_snapshot
    from Numerical.parallel_t3.complete.two_model_master_wilson_batch13 import (
        gauge_sectors, sm_yukawa_tensor, scalar_quartic_from_matchete,
        add_d2_extra_hhss,add_d2_extra_scalar_only,projection)
    from Numerical.running.Aminus4EFT1Couplings import F_LABELS
    snap=load_snapshot(root,model)
    fstate=snap.trajectory.at_fermion_threshold()
    basis=np.eye(8)
    if model=='D_-2':
        from Numerical.running.GeneralT3DynamicFermionWilson import (
            uv_from_json, nonmixed_boundary)
        from Numerical.joint_scalar_matching.special_mass_v6 import build_special_basis
        from Numerical.joint_scalar_matching.source_quartic_v4 import build_source_channels
        uv=uv_from_json(model,snap.uv_values,mu_f=snap.running.mu_f)
        source,_=nonmixed_boundary(model,uv,repo_root=root,mu_f=snap.running.mu_f)
        source_ch=build_source_channels(snap.quartic_seed,model)
        basis=np.asarray(build_special_basis(model,snap.running,source)[
            'R_original_from_boundary'],dtype=float)
    sectors=gauge_sectors(fstate,model,basis)
    Y=sm_yukawa_tensor(fstate,tuple(F_LABELS))
    seed=read_obj(snap.quartic_seed)
    from Numerical.parallel_t3.complete.native_triplet_basis_batch15 import (
        heavy_source_from_reference)
    # Source-only topology terms use native Cartesian triplet; independent
    # RGBeta reference quartic columns use spherical spin-1 coordinates.
    native_R,basis_info=heavy_source_from_reference(model,seed)
    Q,info=scalar_quartic_from_matchete(seed,fstate,model,basis,
                                        native_source_from_reference=native_R)
    if model=='D_-2':
        add_d2_extra_hhss(Q,root,fstate,basis)
        add_d2_extra_scalar_only(Q,fstate,basis)
    ward={}
    for name,(g,gens) in zip(('SU2','U1','SU3'),sectors):
        records=[]
        for a,(tf,ts) in enumerate(gens):
            row={'generator':a+1,
                 'tree_Weinberg_parent':operator_ward(c,tf,ts),
                 'candidate_ADM_delta':operator_ward(delta,tf,ts),
                 'SM_Yukawa':yukawa_ward(Y,tf,ts),
                 'scalar_potential':scalar_quartic_ward(Q,ts)}
            records.append(row)
        ward[name]=records
    maxima={field:max(row[field] for rows in ward.values() for row in rows)
            for field in ('tree_Weinberg_parent','candidate_ADM_delta',
                          'SM_Yukawa','scalar_potential')}
    src_projection=projection(delta,tuple(F_LABELS))
    require_p=decoded_matrix(b13['master_candidate_EFT1_log_GeVinv'])
    insist(relative(src_projection,require_p)<1e-10,
           'full master NPZ and reported LLHH projection disagree')
    # A ward test requires the boundary as well as the beta to be invariant.
    # This cannot by itself establish the correct physical ADM coefficients.
    # Independent controlled limit: remove all BSM scalar interactions and
    # embed a gauge-invariant SM Weinberg operator. Compare the candidate
    # generic psi^2 phi^2 beta to the *separate* WeinbergRGE implementation.
    # This tests the fermion/scalar factors and gauge coefficient in a
    # nontrivial 3-generation limit; it cannot certify arbitrary BSM ADM.
    from Numerical.parallel_t3.wilson.basis import build_basis
    from Numerical.parallel_t3.wilson.coordinates import embed_c5
    from Numerical.parallel_t3.wilson.adm import beta
    from RGE.running.weinberg.WeinbergRGE import sm_weinberg_beta
    from Numerical.parallel_t3.complete.two_model_master_wilson_batch13 import decode_complex
    K=np.array([[.83+.02j,.12+.14j,-.04+.03j],
                [.12+.14j,.37-.05j,.06-.08j],
                [-.04+.03j,.06-.08j,.15+.09j]],complex)
    toy_basis=build_basis(model)
    cw=embed_c5(K,toy_basis)
    lh=float(np.real(decode_complex(fstate['lambdaH']).item()))
    qSM=lh*_baseline_reference_columns(model)['lambdaH']
    variant=b13.get('candidate_ADM','fonseca_eq22_six_term')
    components,infoSM=beta(cw,Y,qSM,sectors,variant=variant)
    numerator=sum(components.values(),np.zeros_like(cw))
    predicted=projection(numerator,tuple(F_LABELS))
    real=lambda k:float(decode_complex(fstate[k]).real)
    expected=sm_weinberg_beta(
        gY=real('gY'),g2_value=real('g2'),g3=real('g3'),
        lambdaH_value=lh,ye=np.asarray(decode_complex(fstate['ye']),complex),
        yu=np.asarray(decode_complex(fstate['yu']),complex),
        yd=np.asarray(decode_complex(fstate['yd']),complex),K=K)[-1]
    sm_limit_residual=relative(predicted,expected)
    sm_weinberg_limit={
        'candidate_vs_independent_SM_Weinberg_beta_relative':sm_limit_residual,
        'tested_operator':'Full 3-generation LLHH embedded with Higgs potential only',
        'reference':'RGE.running.weinberg.WeinbergRGE.sm_weinberg_beta',
        'candidate_backend':infoSM,
        'SM_limit_numerically_consistent':sm_limit_residual<1e-8,
        'full_general_ADM_certified_by_SM_limit':False,
    }
    threshold=1e-8
    passed=all(float(v)<=threshold for v in maxima.values())
    return {
        'source':'saved complete Batch13 C(muF), seven delta channels; '
                 'fresh source-based full SM gauge generators and quartic/Yukawa tensors',
        'max_relative_Ward_residuals':maxima,
        'all_generators':ward if full_ward else {},
        'candidate_ADM_variant':b13.get('candidate_ADM'),
        'source_CG_triplet_basis_info':basis_info,
        'SM_Weinberg_independent_limit':sm_weinberg_limit,
        'full_tensor_candidate_invariance_passed':passed,
        'tree_and_beta_separately_gauge_invariant_at_this_point':passed,
        'independent_gauge_complete_ADM_derivation_or_general_representation_proof':False,
        'full_operator_ADM_certified':False,
        'gate_passed':False,
        'reason':'Gauge covariance at one benchmark is a necessary, not a sufficient, '
                 'condition for the correct seven-term anomalous dimensions.',
    }


def _baseline_reference_columns(model:str)->dict[str,np.ndarray]:
    """Independent RGBeta-potential *candidate* reference, not native CG data."""
    from Numerical.running.Aminus4EFT1Couplings import (
        _quartic_product, su2_generators, quartic_columns)
    if model=='A_-4': return quartic_columns()
    assert model=='D_-2'
    def block(t:np.ndarray,offset:int)->np.ndarray:
        d=t.shape[0]
        out=np.zeros((12,12),float)
        for p in range(d):
            for q in range(d):
                x,y=t[p,q].real,t[p,q].imag
                out[offset+2*p:offset+2*p+2,offset+2*q:offset+2*q+2]=[[x,-y],[y,x]]
        return out
    h=block(np.eye(2),0);s1=block(np.eye(3),4);s2=block(np.eye(1),10)
    adjh=[block(z,0) for z in su2_generators(2)]
    adj1=[block(z,4) for z in su2_generators(3)]
    return {
      'lambdaH':.5*_quartic_product(h,h),
      'lambdaS1':.5*_quartic_product(s1,s1),
      'lambdaS2':.5*_quartic_product(s2,s2),
      'lambdaH1':_quartic_product(h,s1),
      'lambdaH2':_quartic_product(h,s2),
      'lambda12':_quartic_product(s1,s2),
      'lambdaS1Adj':sum((.5*_quartic_product(x,x) for x in adj1),np.zeros((12,)*4)),
      'lambdaH1Adj':sum((_quartic_product(a,b) for a,b in zip(adjh,adj1)),np.zeros((12,)*4)),
    }


def numeric_quartic_source_channel(seed:dict,items:list[dict],model:str,
                                    coupling_name:str)->np.ndarray:
    """Exactly evaluate a native Matchete CG channel with real unit coupling."""
    from RGE.running.intermediate.ScalarOnlyTensorAdapters import (
        build_eft1_quartic_tensor)
    import sympy as sp
    dims={'A_-4':(1,3),'D_-2':(3,1)}[model]
    part={'CGRegistry':seed['CGRegistry'],'ScalarQuarticTerms':items}
    tensor=build_eft1_quartic_tensor(part,d_s1=dims[0],d_s2=dims[1],
                                    shared_scalar=False)
    arr=np.zeros((12,)*4,float)
    x=sp.Symbol(coupling_name)
    from itertools import permutations
    for indices,formula in tensor.nonzero_items():
        z=complex(sp.N(formula.subs({x:1}),18))
        insist(abs(z.imag)<1e-9*max(1.,abs(z.real)),
               'native coupling unit tensor is not real for '+coupling_name)
        for idx in set(permutations(tuple(int(v)-1 for v in indices))):
            arr[idx]=z.real
    return arr


def native_quartic_map_audit(seed:dict,model:str)->dict:
    """Compare archived native CG tensors with independent RGBeta candidate columns.

    This detects normalization/sign differences. Even exact agreement of the
    sampled *archived* channels does not establish the missing D quartics or
    the full RG beta basis transformation.
    """
    from Numerical.parallel_t3.complete.two_model_master_wilson_batch13 import (
        CANDIDATE_SOURCE_ALIASES)
    cols=_baseline_reference_columns(model)
    from Numerical.parallel_t3.complete.native_triplet_basis_batch15 import heavy_source_from_reference
    native_R,basis_meta=heavy_source_from_reference(model,seed)
    full_R=np.eye(12);full_R[4:,4:]=native_R
    grouped={}
    for row in seed['ScalarQuarticTerms']:
        found=re.findall(r'Coupling\[([^,\]]+),',row['TermInputForm'])
        insist(len(found)==1,'source quartic monomial contains multiple couplings')
        grouped.setdefault(found[0],[]).append(row)
    candidates={}
    for name,terms in sorted(grouped.items()):
        alias=CANDIDATE_SOURCE_ALIASES.get(name)
        if alias is None:
            candidates[name]={'status':'NoIndependentlyIdentifiedRGBetaCoordinate'}
            continue
        if alias not in cols:
            candidates[name]={'status':'RGBetaReferenceTopologyNotIndependent',
                              'candidate_alias':alias}
            continue
        src=numeric_quartic_source_channel(seed,terms,model,name)
        # This input is Matchete CARTESIAN T3BSMd3; the independent RGBeta
        # comparison column is already highest-weight SPHERICAL spin-1.
        src=np.einsum('abcd,aA,bB,cC,dD->ABCD',src,
                      full_R,full_R,full_R,full_R,optimize=True)
        ref=np.asarray(cols[alias],complex)
        nr=float(np.vdot(ref,ref).real)
        insist(nr>0,'reference quartic channel unexpectedly empty')
        factor=float(np.vdot(ref,src).real/nr)
        fit=float(np.linalg.norm(src-factor*ref)/max(np.linalg.norm(src),1e-100))
        candidates[name]={'status':'NativeVersusRGBetaPotentialDirectionMeasured',
                          'candidate_alias':alias,'native_over_reference_scale':factor,
                          'relative_direction_residual':fit,
                          'one_real_unit_coupling_only':True,
                          'coordinate_normalization_certified':False}
    return {'archived_native_channel_comparison':candidates,
            'source_cartesian_to_spherical_CG_intertwiner':basis_meta,
            'D2_missing_extra_HHSS_reference_basis_certified':False if model=='D_-2' else None,
            'all_complex_coupling_phases_and_orientations_verified':False,
            'matchete_to_RGBeta_complete_coordinate_map_certified':False,
            'gate_passed':False}


def run_one(root:Path,model:str,*,output:Path|None=None,
            compare_native_quartics:bool=True,include_generator_details:bool=False)->dict:
    insist(model in MODELS,'unsupported model')
    root=Path(root).resolve()
    ledger_file=root/B2/model/'one_loop_operator_ledger.json'
    b11_file=root/B11/model/'result.json'
    b13_file=root/B13/model/'result.json'
    npz=root/B13/model/'full_tensor_candidate.npz'
    source=root/SOFT[model]
    native_status=source.with_name('native_two_threshold_status.json')
    input_files={'full_F_operator_ledger':ledger_file,'B11_native_matching':b11_file,
                 'B13_master_wilson':b13_file,'B13_full_tensor':npz,
                 'native_scalar_soft_expression':source,'native_scalar_status':native_status}
    src=read_obj(ledger_file);b11=read_obj(b11_file);b13=read_obj(b13_file)
    meta=read_obj(native_status)
    insist(b13.get('model')==b11.get('model')==model and
           b13.get('status')=='FullTensorSevenChannelStrictOneLoopConditional_NotPhysicalC5'
           and b11.get('status')=='NativeHardSoftLaurentAndScaleLogClosed_ConditionalOperatorConvention',
           'previous matching/Wilson sources not successful for '+model)
    insist(meta.get('model')==model,'native soft status model mismatch')
    matching=audit_matching(src,b11,meta,model,source.read_text(encoding='utf-8-sig'))
    normalization=audit_normalization(src,b13)
    ward=audit_gauge(root,model,npz,b13,full_ward=include_generator_details)
    snap=None
    quartic={}
    if compare_native_quartics:
        from Numerical.parallel_t3.complete.two_model_joint_scalar_batch9 import load_snapshot
        snap=load_snapshot(root,model)
        quartic=native_quartic_map_audit(read_obj(snap.quartic_seed),model)
        input_files['native_quartic_seed']=snap.quartic_seed
    unresolved={
        'GammaCC_to_Weyl_absolute_normalization':not normalization['gate_passed'],
        'complete_finite_MSbar_matching':not matching['gate_passed'],
        'complete_gauge_ADM_and_scalar_map':not ward['gate_passed'] or
                  bool(quartic and not quartic['gate_passed']),
    }
    result={
       'format':FORMAT,'model':model,
       'status':'ThreePhysicsClosureGatesInvestigated_StillConditionalNotPhysicalC5',
       'hierarchy':[['UV','F_all_3_flavours'],['EFT1','S1_and_S2_together'],['SMEFT_100GeV']],
       'normalization_gate':normalization,'finite_matching_gate':matching,
       'gauge_Wilson_gate':ward,'native_RGBeta_scalar_map_gate':quartic,
       'remaining_unresolved_physics_gates':unresolved,
       'all_three_physics_gates_certified':False,
       'physical_C5_certified':False,
       'B11_B12_B13_outputs_modified':False,
       'additional_yukawa_fit_performed':False,
       'input_sha256':{k:hash_file(p) for k,p in input_files.items()},
       'source_conclusion':'Existing input files do not independently fix the absolute '
           'GammaCC/Weyl sign and normalisation or guarantee complete one-loop '
           'MSbar matching. A successful native/master log and gauge-Ward check '
           'is necessary but does not fix either normalization or completeness.',
    }
    if output is not None:
        dest=Path(output);dest.mkdir(parents=True,exist_ok=True)
        target=dest/'result.json'
        target.write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)
                          +'\n',encoding='utf-8')
    return result


def preflight_paths(root:Path,model:str)->tuple[Path,...]:
    root=Path(root)
    src=root/SOFT[model]
    return (root/B2/model/'one_loop_operator_ledger.json',
            root/B11/model/'result.json',root/B13/model/'result.json',
            root/B13/model/'full_tensor_candidate.npz',src,
            src.with_name('native_two_threshold_status.json'),
            root/B1/model/'manifest.json',root/B9/model/'result.json')
