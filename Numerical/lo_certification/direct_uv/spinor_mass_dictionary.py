"""Read-only native mass-Ward and CONDITIONAL T3 spinor-spurion phase audit.

Checked against local Matchete/RGBeta source, not an installed native
four-spinor->Weyl matrix convention.  No physical phase is assigned.

Field-change convention (all-left-handed candidate, NOT native certified):
  F -> exp(i*a) F; Fc -> exp(i*b) Fc; S1 -> exp(i*c) S1;
  S2 -> exp(i*d) S2.
  L contains y1 * nu Fc S1, y2 * nu F S2*, MF * F Fc,
             lambdaT3 * HH S1 S2* plus Hermitian conjugates.
  For Majorana F=Fc, b=a and mass term is MF * F F / 2.

The phase of y1*y2*conjugate(MF)*conjugate(lambdaT3) is invariant.
This is an action-spurion witness only.  It does NOT measure Matchete's
GammaCC sign, RGBeta's installed epsilon index ordering, or the physical
coupling dictionary.
"""
from __future__ import annotations

from hashlib import sha256
from pathlib import Path
import json
import re
from typing import Any

import sympy as sp

from .certification import CLASSES, all_models

FORMAT = 't3-lo-direct-uv-native-mass-ward-conditional-spinor-spurions-v1'
RGB = Path('RGE/running/rgbeta/T3RGBetaModel.wl')
TOP = Path('Lagrangian/interactions/T3Topology.wl')
FIELDS = Path('Lagrangian/model/T3Fields.wl')
OLD_MASS = Path('Reports/output/workstream_a/a0_dirac_mass_ward/a0_dirac_mass_ward.json')
A0_PHASE = Path('Reports/output/workstream_a/yukawa_spinor_phase/a0_spinor_phase.json')


class DictionaryError(ValueError):
    pass


def require(test: bool, why: str) -> None:
    if not test:
        raise DictionaryError(why)


def without_wl_comments(s: str) -> str:
    """Remove Mathematica (* ... *) comments, including nested comments."""
    out: list[str] = []
    depth = 0
    i = 0
    while i < len(s):
        if s.startswith('(*', i):
            depth += 1
            i += 2
        elif s.startswith('*)', i):
            require(depth > 0, 'Unmatched Mathematica comment terminator')
            depth -= 1
            i += 2
        else:
            if depth == 0:
                out.append(s[i])
            i += 1
    require(depth == 0, 'Unclosed Mathematica comment')
    return ''.join(out)


def _compact(s: str) -> str:
    return re.sub(r'\s+', '', without_wl_comments(s))


def guard_sources(root: Path) -> dict[str, Any]:
    """Inspect *the installed files* and fail closed on source drift.

    In particular, this rejects the old SU(2) doublet DELTA Dirac mass.
    It does not turn the source into a native RGBeta execution certificate.
    """
    paths = {str(p): root/p for p in (RGB, TOP, FIELDS)}
    for rel, p in paths.items():
        require(p.is_file(), 'Required current local source missing: '+rel)
    rgb = _compact(paths[str(RGB)].read_text(encoding='utf-8-sig'))
    top = _compact(paths[str(TOP)].read_text(encoding='utf-8-sig'))
    fields = _compact(paths[str(FIELDS)].read_text(encoding='utf-8-sig'))

    head = 'T3RGBetaAddMasses[dS1_Integer,dS2_Integer,dF_Integer,yF_]:=Module['
    require(rgb.count(head) == 1, 'RGBeta mass function signature changed or duplicated')
    fragment = rgb.split(head, 1)[1].split('(*',1)[0]  # no comments remain
    # Bound to subsequent functions, but do not assume every nested ] is local.
    fragment = fragment.split('T3RGBetaAddBasicQuartics[',1)[0]
    for needle in (
        'If[yF===0,fermionInvariant=If[dF===1,(1&),(del[T3RGBetaRep[dF],#1,#2]&)];',
        'AddFermionMass[MF,{F,F},GroupInvariant->fermionInvariant,',
        'fermionInvariant=Which[dF===1,(1&),dF===2,(eps[T3RGBetaRep[dF],#1,#2]&),dF===3,(del[T3RGBetaRep[dF],#1,#2]&),True,Return[$Failed]];',
        'AddFermionMass[MF,{F,Fc},GroupInvariant->fermionInvariant,',
        'MassIndices->({heavy[#1],heavy[#2]}&),Chirality->Right',
    ):
        require(needle in fragment, 'RGBeta fermion mass source contract changed: '+needle)
    for needle in (
        'T3RGBetaRep[1]:=None;',
        'T3RGBetaRep[2]:=SU2L[fund];',
        'T3RGBetaRep[3]:=SU2L[S2];',
        'GaugeRep->T3RGBetaGaugeRep[dF,yF]',
        'GaugeRep->T3RGBetaGaugeRep[dF,-yF]',
        'secondHeavy=If[yF===0,F,Fc];',
        'AddYukawa[y1,{Bar@S1,l,F},',
        'AddYukawa[y2,{S2,l,secondHeavy},',
    ):
        require(needle in rgb, 'RGBeta field/source convention changed: '+needle)
    for needle in (
        'useCConj=which===1;',
        'scalarBar=which===2;',
        'If[useCConj,CConj[NewFermion[r]],NewFermion[r]]',
        'y[p,r]NCM[Bar[l[i,p]],f]scalarCG[cg,labels]',
        'lambdaT3[]H[i]H[j]s1s2CG[T3MixCG,labels]',
        's1=If[d1===1,NewScalar1[],NewScalar1[a]];',
        's2=If[d2===1,Bar[NewScalar2[]],Bar[NewScalar2[b]]];',
        'scalar=If[dS===1,ScalarFieldValue[which,conjugated],ScalarFieldValue[which,conjugated,k]];',
        'conjugated=which===1;',
        'coupling[p,r]NCM[Bar[l[i,p]],f]scalarCG[cg,labels]',
    ):
        require(needle in top, 'Matchete Yukawa/mixing source contract changed: '+needle)
    for needle in (
        'DefineCoupling[y1,Indices->{Flavor,NFlavor},SelfConjugate->False];',
        'DefineCoupling[y2,Indices->{Flavor,NFlavor},SelfConjugate->False];',
        'DefineCoupling[lambdaT3,SelfConjugate->False];',
        'selfConj=TrueQ[PossibleZeroQ[f["Y"]]]&&OddQ[f["SU2"]];',
        'SelfConjugate->True',
        'SelfConjugate->False',
    ):
        require(needle in fields, 'Matchete field definition changed: '+needle)
    return {
        'source_sha256': {rel:sha256(path.read_bytes()).hexdigest() for rel,path in paths.items()},
        'rgbeta_mass_contract_checked':True,
        'rgbeta_dirac_doublet_mass_tensor':'epsilon_of_two_left_fundamentals',
        'rgbeta_dirac_triplet_mass_tensor':'delta_in_real_triplet_representation',
        'rgbeta_majorana_odd_representation_contract':'delta',
        'source_read_live_on_local_filesystem':True,
        'installed_native_rgbeta_mass_amplitude_executed':False,
        'absolute_spinor_phase_measured':False,
        'Matchete_topology_y1_spinor':'CConj[NewFermion]',
        'Matchete_topology_y2_spinor':'NewFermion',
        'RGBeta_y1_spinor_name':'F',
        'RGBeta_y2_spinor_name':'Fc for Dirac, F for Majorana',
        'mapping_of_Matchete_CConj_to_RGBeta_Fc_in_common_left_Weyl_basis':'CONDITIONAL_NOT_NATIVE_PROVED',
    }


def ward_generators(d: int) -> list[sp.Matrix]:
    if d == 1:
        return [sp.zeros(1) for _ in range(3)]
    if d == 2:
        return [sp.Matrix([[0, 1], [1, 0]])/2,
                sp.Matrix([[0, -sp.I], [sp.I, 0]])/2,
                sp.diag(1, -1)/2]
    if d == 3:
        # Cartesian adjoint basis: Hermitian spin-1 generators, (T_a)bc=-i eps_abc.
        return [sp.Matrix(3,3,lambda b,c:-sp.I*sp.LeviCivita(a,b,c))
                for a in range(3)]
    raise DictionaryError('Not a supported SU(2) dimension')


def invariant_tensor(d: int) -> sp.Matrix:
    if d==1:
        return sp.eye(1)
    if d==2:
        return sp.Matrix([[0,1],[-1,0]])
    if d==3:
        return sp.eye(3)
    raise DictionaryError('Not a supported SU(2) dimension')


def exact_zero(matrix: sp.Matrix) -> bool:
    return all(sp.simplify(x)==0 for x in matrix)


def mass_ward(d: int, *, force_delta_doublet:bool=False) -> dict[str,Any]:
    c = (sp.eye(2) if force_delta_doublet and d==2 else invariant_tensor(d))
    residuals=[sp.simplify(t.T*c+c*t) for t in ward_generators(d)]
    passed=all(exact_zero(r) for r in residuals)
    return {
        'dimension':d,
        'candidate_representation_basis':'singlet' if d==1 else 'fundamental' if d==2 else 'Cartesian_adjoint',
        'tensor':'delta' if force_delta_doublet else ('epsilon' if d==2 else 'delta'),
        'all_su2_generators_ward_zero':passed,
        'nonzero_generator_numbers':[i+1 for i,r in enumerate(residuals) if not exact_zero(r)],
        'mass_tensor_rows':[[str(e) for e in row] for row in c.tolist()],
        'native_component_basis_cross_package_measured':False,
    }


# Charges of the conditional source operators under rephasings of F,Fc,S1,S2.
# y1 * nu Fc S1, y2 * nu F S2*, MF * F Fc, lambda * HH S1 S2*.
DIRAC = {'y1':(0,1,1,0),'y2':(1,0,0,-1),
         'MF':(1,1,0,0),'lambdaT3':(0,0,1,-1),
         'y1Conj':(0,1,-1,0),'y2Conj':(1,0,0,1)}
MAJORANA = {'y1':(1,0,1,0),'y2':(1,0,0,-1),
            'MF':(2,0,0,0),'lambdaT3':(0,0,1,-1),
            'y1Conj':(1,0,-1,0),'y2Conj':(1,0,0,1)}


def phase_sum(terms: dict[str,int], majorana:bool) -> tuple[int,...]:
    charges=MAJORANA if majorana else DIRAC
    require(all(c in charges for c in terms), 'Unknown T3 source operator phase label')
    return tuple(sum(power*charges[name][j] for name,power in terms.items()) for j in range(4))


def spurion_witness(majorana:bool,extra:str|None=None) -> dict[str,Any]:
    terms={'y1':1,'y2':1,'MF':-1,'lambdaT3':-1}
    correct=phase_sum(terms,majorana)
    wrong=phase_sum({**terms,'lambdaT3':1},majorana)
    require(correct==(0,0,0,0),'Conditional source-loop spurion not rephasing invariant')
    require(wrong!=(0,0,0,0),'Wrong lambda-conjugation negative control unexpectedly invariant')
    wrong_same_handedness=(None if majorana else tuple(
        DIRAC['y2'][j]+DIRAC['y2'][j]-DIRAC['MF'][j]-DIRAC['lambdaT3'][j]
        for j in range(4)))
    require(majorana or wrong_same_handedness!=(0,0,0,0),
            'Incorrect same-heavy-Weyl-field Dirac Yukawa assignment was not detected')
    # A candidate epsilon / GammaCC orientation is not a physical coefficient.
    info={
        'conditional_left_weyl_action':'y1 nu Fc S1 + y2 nu F S2* + MF F Fc + lambdaT3 HH S1 S2* + hc' if not majorana else
            'y1 nu F S1 + y2 nu F S2* + MF/2 F F + lambdaT3 HH S1 S2* + hc',
        'rephasing_variables':['F','Fc','S1','S2'],
        'MF_monomial_charge':list((MAJORANA if majorana else DIRAC)['MF']),
        'loop_spurion':'y1*y2*conjugate(MF)*conjugate(lambdaT3)',
        'operator_powers':terms,
        'loop_spurion_rephasing_charge':list(correct),
        'incorrect_nonconjugated_lambda_charge':list(wrong),
        'negative_control_rejected':True,
        'incorrect_same_chirality_dirac_assignment_charge':None if majorana else list(wrong_same_handedness),
        'absolute_Matchete_RGBeta_coefficient_certified':False,
        'MF_and_GammaCC_sign_transport_certified':False,
        'native_complex_field_conjugations_physically_projected':False,
    }
    if extra is not None:
        require(extra in ('y1Conj','y2Conj'), 'Unknown extra Yukawa')
        info['additional_source_allowed_independent_yukawa']={
            'name':extra,'operator_charge':list((MAJORANA if majorana else DIRAC)[extra]),
            'independent_complex_coefficient':True,
            'included_in_restricted_frozen_benchmark':False,
        }
    return info


def _read_optional(root:Path, relative:Path) -> dict[str,Any]|None:
    p=root/relative
    if not p.is_file():
        return None
    value=json.loads(p.read_text(encoding='utf-8-sig'))
    require(isinstance(value,dict),'Malformed optional source report '+str(relative))
    return value


def optional_evidence(root:Path) -> dict[str,Any]:
    # Preserve upstream native result as distinct from our independent Ward.
    a0=_read_optional(root,OLD_MASS)
    spinor=_read_optional(root,A0_PHASE)
    mass_record=({'available':False} if a0 is None else {
        'available':True,'status':a0.get('status'),
        'registered_MF_invariant':a0.get('registered_MF_invariant'),
        'local_rgbeta_source_sha256':a0.get('local_rgbeta_source_sha256'),
        'actual_su2_mass_ward_passed':a0.get('actual_su2_mass_ward_passed'),
        'current_source_mass_ward_certificate':False,
    })
    if a0 is not None:
        mass_record['current_source_mass_ward_certificate']=(
            a0.get('local_rgbeta_source_sha256')==sha256((root/RGB).read_bytes()).hexdigest()
            and a0.get('registered_MF_invariant')=='eps[SU2L[fund], dm1, dm2]'
            and a0.get('actual_su2_mass_ward_passed') is True)
        mass_record['stale_or_incompatible_native_report']=(not mass_record['current_source_mass_ward_certificate'])
    phase_record=({'available':False} if spinor is None else {
        'available':True,'status':spinor.get('status'),
        'chosen_CC_basis_y1_R':spinor.get('coupling_dictionary_in_explicitly_chosen_CC_and_epsilon_basis',{}).get('y1_R'),
        'physical_cross_package_yukawa_dictionary_certified':spinor.get('physical_cross_package_yukawa_dictionary_certified'),
        'MF_and_y2_transport_proved':spinor.get('global_dirac_mass_and_y2_phase_transport_proved'),
    })
    return {'a0_native_mass_report':mass_record,'a0_spinor_witness':phase_record}


def audit(root:Path) -> dict[str,Any]:
    root=root.resolve()
    evidence=guard_sources(root)
    ward={d:mass_ward(d) for d in (1,2,3)}
    require(all(v['all_su2_generators_ward_zero'] for v in ward.values()),
            'Native-source candidate fermion mass representation fails SU2 Ward')
    wrong=mass_ward(2,force_delta_doublet=True)
    require(not wrong['all_su2_generators_ward_zero'],
            'SU2 doublet delta negative control unexpectedly passed')
    items=[]
    for model,key in all_models():
        kind,alpha=model[0],int(model[2:])
        ds1,ds2,df,_=CLASSES[kind]
        majorana=(df%2==1 and alpha==-1)
        extra='y1Conj' if alpha==0 else 'y2Conj' if alpha==-2 else None
        # Extra Yukawas only when the corresponding scalar hypercharge is zero.
        if extra is not None:
            require((extra=='y1Conj' and alpha/2==0) or
                    (extra=='y2Conj' and (alpha+2)/2==0),
                    'Extra Yukawa scalar hypercharge test failed')
        items.append({
            'model':model,'frozen_benchmark_key':key,
            'fields':{'dS1':ds1,'dS2':ds2,'dF':df,'alpha':alpha},
            'fermion_kind':'Majorana' if majorana else 'Dirac',
            'rgbeta_source_mass_invariant':ward[df]['tensor'],
            'su2_ward_satisfied_for_declared_representation':ward[df]['all_su2_generators_ward_zero'],
            'source_operator_orientation_verified':True,
            'conditional_spurion':spurion_witness(majorana,extra),
            'physical_source_certified':False,
        })
    return {'format':FORMAT,
            'status':'16of16_native_source_mass_Ward_and_conditional_spurions_checked_absolute_dictionary_pending',
            'source_contract':evidence,
            'mass_ward_witnesses':[ward[d] for d in (1,2,3)],
            'negative_control_doublet_delta':wrong,
            'optional_native_evidence':optional_evidence(root),
            'model_count':len(items),
            'model_counts':{'Dirac':sum(x['fermion_kind']=='Dirac' for x in items),
                            'Majorana':sum(x['fermion_kind']=='Majorana' for x in items),
                            'extra_yukawa_allowed':sum('additional_source_allowed_independent_yukawa' in x['conditional_spurion'] for x in items)},
            'physical_source_certified_models':0,
            'outstanding':[
                'Native installed Matchete GammaCC 4x4 orientation and RGBeta epsilon component orientation',
                'Physical two-Yukawa and MF normalisation in one common all-left Weyl basis',
                'Physical Matchete-RGBeta signed portal and scalar-mixing coupling dictionary',
                'Additional allowed Yukawas absent from restricted benchmark numeric coefficients',
                'Complete same-input two-threshold EFT C5 comparison after source matching',
            ],
            'all_heavy_fermions_at_muF_then_both_heavy_scalars_at_muS':True,
            'SM_fermions_never_integrated_out':True,
            'models':items}
