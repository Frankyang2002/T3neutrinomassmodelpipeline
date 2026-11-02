"""B2: an exact, assumption-labelled SU(2) Ward witness for RGBeta contracts.

The *literal* field interpretation is fund for H,S1,l and anti-fund for
Bar[S1],Bar[S2]. RGBeta's abstract index/Bar semantics have not been certified
numerically. Consequently a failed literal Ward test is a **diagnostic**, not
proof of a package bug. In particular, it cannot grant or reject a physical
Matchete/RGBeta coupling certificate on its own.

No production imports, writes, or Lagrangian modifications.
"""
from __future__ import annotations

import hashlib
import json
import re
from itertools import product
from pathlib import Path
from typing import Any, Callable

import sympy as sp

from b2_gauge_coupling import (
    AUDITS, CLASSES, FLAGS, SOURCE_BLOBS, AuditBlocked, exact, git_result,
    models, native_tensor, require, source_rotations,
)

FUND = 'fund'
ANTI = 'anti'
EPS = sp.Matrix([[0, 1], [-1, 0]])
PAULI = (
    sp.Matrix([[0, 1], [1, 0]]) / 2,
    sp.Matrix([[0, -sp.I], [sp.I, 0]]) / 2,
    sp.diag(sp.Rational(1, 2), -sp.Rational(1, 2)),
)


def delta(a: int, b: int) -> sp.Integer:
    return sp.Integer(a == b)


def literal_rgbeta_22_mix(indices: tuple[int, int, int, int]) -> sp.Expr:
    """FORMAL delS2*delS2 expansion observed in the native RGBeta export.

    Index order (H,H,S1,Bar[S2]). Not an asserted physical component tensor.
    """
    h, h2, s, sb = indices
    return (delta(h, sb)*delta(h2, s) + delta(h2, sb)*delta(h, s))/2


def repaired_22_mix(indices: tuple[int, int, int, int]) -> sp.Expr:
    """One explicitly convention-chosen singlet of 2⊗2⊗2⊗2bar.

    The physical source action normalization and RGBeta's index interpretation
    have NOT been measured. This is not a production replacement proposal.
    """
    h, h2, s, sb = indices
    return (delta(h, sb)*EPS[h2, s] + delta(h2, sb)*EPS[h, s])/2


def literal_B_y1(indices: tuple[int, int]) -> sp.Expr:
    """Formal source eps on (Bar[S1],l), B singlet heavy fermion."""
    a, b = indices
    return EPS[a, b]


def repaired_B_y1(indices: tuple[int, int]) -> sp.Expr:
    a, b = indices
    return delta(a, b)


def literal_B_y2(indices: tuple[int, int]) -> sp.Expr:
    """Formal source eps on (S2,l), B singlet heavy fermion."""
    a, b = indices
    return EPS[a, b]


def ward_residual(tensor: Callable[[tuple[int, ...]], sp.Expr],
                  reps: tuple[str, ...]) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    """Exact Σ_axis C(...j...) R_axis[j,i]; bars use -T* as generators."""
    require(bool(reps) and set(reps) <= {FUND, ANTI}, 'bad SU2 representations')
    results = []
    for t in PAULI:
        nonzero = []
        for indices in product(range(2), repeat=len(reps)):
            value = sp.S.Zero
            for axis, rep in enumerate(reps):
                representation = t if rep == FUND else -t.conjugate()
                for j in range(2):
                    switched = list(indices)
                    switched[axis] = j
                    value += tensor(tuple(switched)) * representation[j, indices[axis]]
            nonzero.append(sp.simplify(value))
        results.append(sp.simplify(sum(sp.conjugate(z)*z for z in nonzero)))
    return tuple(results)  # type: ignore[return-value]


def norm2(tensor: Callable[[tuple[int, ...]], sp.Expr], n: int) -> sp.Expr:
    return sp.simplify(sum(sp.conjugate(tensor(idx))*tensor(idx)
                           for idx in product(range(2), repeat=n)))


def invariant_identity() -> dict[str, Any]:
    tests = (
        ('rgbeta_22_mix_literal_delta_delta', literal_rgbeta_22_mix,
         (FUND, FUND, FUND, ANTI), False),
        ('convention_chosen_22_mix_delta_epsilon', repaired_22_mix,
         (FUND, FUND, FUND, ANTI), True),
        ('rgbeta_B_y1_literal_epsilon', literal_B_y1,
         (ANTI, FUND), False),
        ('convention_chosen_B_y1_delta', repaired_B_y1,
         (ANTI, FUND), True),
        ('rgbeta_B_y2_literal_epsilon', literal_B_y2,
         (FUND, FUND), True),
    )
    rows = []
    for name, fn, reps, expected in tests:
        residual = ward_residual(fn, reps)
        ok = all(v == 0 for v in residual)
        require(ok is expected, f'{name}: unexpected exact SU2 Ward outcome')
        n2 = norm2(fn, len(reps))
        require(n2 > 0, f'{name}: identically zero tensor')
        rows.append({'name': name, 'ordered_representations': list(reps),
                     'ward_generator_squared_norms': [str(x) for x in residual],
                     'tensor_squared_norm': str(n2), 'ward_exact_zero': ok,
                     'component_axes_are_RGBeta_measured': False})
    # independent negative check: the 2⊗2 Kronecker delta is not an SU2 singlet.
    require(any(x != 0 for x in ward_residual(lambda ab: delta(*ab), (FUND,FUND))),
            'fundamental delta wrongly classified as an invariant')
    require(all(x == 0 for x in ward_residual(lambda ab: EPS[ab[0],ab[1]],
                                            (FUND,FUND))), 'epsilon not invariant')
    for h, h2, s, sb in product(range(2), repeat=4):
        require(repaired_22_mix((h,h2,s,sb)) ==
                repaired_22_mix((h2,h,s,sb)), 'Higgs Bose symmetry failed')
    return {'basis': 'chosen SU2L fundamental (+1/2,-1/2), epsilon_01=+1',
            'bar_generator': '-conjugate(T)',
            'fields_for_22_mix': ['H','H','S1','Bar[S2]'],
            'b1_yukawa_fields': ['Bar[S1]','l','F'],
            'b2_yukawa_fields': ['S2','l','F'],
            'tensor_checks': rows,
            'native_RGBeta_index_semantics_verified': False,
            'native_RGBeta_action_prefactors_verified': False}


def mix_charge_support(alpha: int, tensor: Callable[[tuple[int, ...]], sp.Expr]) -> list[dict]:
    """H0=weight -1/2; nonzero scalar monomials and their two physical Q."""
    rows=[]
    for s,t in product(range(2), repeat=2):
        cg = sp.simplify(tensor((1,1,s,t)))
        if cg == 0:
            continue
        q1=sp.Rational(1,2)-s+sp.Rational(alpha,2)
        q2=sp.Rational(1,2)-t+sp.Rational(alpha+2,2)
        rows.append({'S1_weight_index_1based':s+1,'S2_weight_index_1based':t+1,
                     'coefficient':str(cg),'charge_S1':str(q1), 'charge_S2':str(q2),
                     'electric_charge_preserved': q1==q2})
    return rows


def all16_charge_warning() -> list[dict]:
    rows=[]
    for m in models():
        if m['family'] not in ('B','C'):
            continue
        a=m['alpha']
        bad=mix_charge_support(a,literal_rgbeta_22_mix)
        good=mix_charge_support(a,repaired_22_mix)
        require(bad and good and not all(r['electric_charge_preserved'] for r in bad)
                and all(r['electric_charge_preserved'] for r in good),
                'expected 2x2 electric-charge anomaly not found')
        rows.append({'model':m['label'],'alpha':a,'literal_source_tensor':bad,
                     'explicit_Ward_invariant_candidate':good,
                     'not_an_observed_RGBeta_numerical_tensor':True})
    require(len(rows)==6,'expected B/C six charge samples')
    return rows


def source_provenance(root: Path) -> dict:
    """Allow new Git HEADs only if ALL six tracked source blobs remain exact.

    Reuses Git's path-specific clean filters for CRLF. Checks actual repo root,
    HEAD tree, index, canonical worktree blob and clean path for every file.
    No mutations, no reliance on a specific commit hash.
    """
    root=root.resolve()
    rc, top=git_result(root,'rev-parse','--show-toplevel')
    require(rc==0 and Path(top).resolve()==root,'not exact requested repository root')
    rc, head=git_result(root,'rev-parse','HEAD')
    require(rc==0 and re.fullmatch(r'[0-9a-f]{40}', head) is not None,
            'cannot measure local Git commit')
    rows=[]
    for rel, expected in SOURCE_BLOBS.items():
        path=root/rel
        require(path.is_file(), 'missing source '+rel)
        tree_rc,tree=git_result(root,'ls-tree','HEAD','--',rel)
        idx_rc,index=git_result(root,'ls-files','--stage','--',rel)
        hash_rc,hash_value=git_result(root,'hash-object',f'--path={rel}','--',rel)
        diff_rc,_=git_result(root,'diff','--quiet','HEAD','--',rel)
        tree_fields=tree.split()
        index_fields=index.split()
        tree_blob=tree_fields[2] if tree_rc==0 and len(tree_fields)>=4 else None
        index_blob=index_fields[1] if idx_rc==0 and len(index_fields)>=4 and index_fields[2]=='0' else None
        require(tree_blob==expected and index_blob==expected and hash_rc==0
                and hash_value==expected and diff_rc==0,
                f'changed/unverified tracked source: {rel}; '
                'do not reuse the old tensor interpretation')
        rows.append({'path':rel,'git_blob_sha':expected,
                     'worktree_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
                     'validated_with_git_clean_filter':True})
    # Do not accept matching blobs but a different source/operator contract.
    rgb=(root/'RGE/running/rgbeta/T3RGBetaModel.wl').read_text(encoding='utf-8-sig')
    from b2_gauge_coupling import compact_wl
    stripped=compact_wl(rgb)
    for token in (
        'T3RGBetaMixInvariant[2,2]:=',
        'delS2[SU2L,a,#1,#2]*delS2[SU2L,a,#3,#4]',
        'AddQuartic[lambdaT3,{H,H,S1,Bar@S2}',
        'T3RGBetaYukawaInvariant[2,1]:=(eps[SU2L@fund,#1,#2]&)',
    ):
        require(token.strip() in stripped,
                'literal formal invariant source changed: '+token)
    return {'repository_head_measured':head,'validated_git_blobs':rows,
            'all_six_source_blobs_unchanged':True,
            'repo_head_must_equal_old_reference':False}


def formal_evidence(data: dict) -> dict:
    require(data.get('format')=='b2-native-rgbeta-formal-invariant-probe-v1' and
            data.get('status')=='FormalRegistriesCapturedPhysicalNormalisationPending' and
            data.get('physical_source_certified') is False and
            data.get('numerical_component_axes_verified') is False and
            data.get('full_RGBeta_action_normalisation_verified') is False,
            'native formal evidence incorrect status or overclaims certification')
    records=data.get('five_classes')
    require(isinstance(records,list) and len(records)==5,'not five native classes')
    found={}
    for row in records:
        family=row.get('family')
        require(family in CLASSES and family not in found and
                row.get('status')=='RegistryAndFormalInvariantCaptured' and
                row.get('physical_action_certified') is False and
                row.get('native_numeric_RGBeta_component_axes_certified') is False,
                'class registry missing/overclaimed')
        d1,d2,df,alphas=CLASSES[family]
        require(row.get('dimensions_alpha')[:3]==[d1,d2,df],
                'native family dims changed')
        found[family]=row
    require(set(found)==set(CLASSES),'native family list not A–E')
    for family in ('B','C'):
        row=found[family]
        require(row['registered_yukawas'][0]['coupling']=='y1' and
                row['registered_yukawas'][0]['actual_registered_fields']=='{Bar[S1], l, F}' and
                row['registered_yukawas'][1]['coupling']=='y2' and
                row['registered_yukawas'][1]['actual_registered_fields']=='{S2, l, F}' and
                row['registered_yukawas'][0]['formal_invariant_inputform']==
                ('eps[SU2L[fund], b2a, b2b]' if family=='B' else
                 'delS2[SU2L, b2c, b2a, b2b]'),
                'B/C native y1 registration changed')
        mix=row['mixing_invariant_formal'].replace(' ','')
        require('del[SU2L[fund]' in mix and mix.count('del[SU2L[fund]')==4
                and mix.count('/2')==2,
                'B/C native mixing formal expansion changed')
    return {'five_classes_formal_records_validated':True,
            'B_C_formal_delta_delta_confirmed':True,
            'B_y1_epsilon_on_BarS1_l_confirmed':True,
            'RGBeta_formal_file_is_numeric_tensor_witness':False,
            'RGBeta_formal_file_is_action_witness':False}


def fit_tensor(reference: Callable[[tuple[int, ...]], sp.Expr],
               observed: Callable[[tuple[int, ...]], sp.Expr],
               n: int) -> dict:
    den=norm2(reference,n)
    require(den!=0,'zero candidate tensor')
    num=sum(sp.conjugate(reference(ix))*observed(ix)
            for ix in product(range(2),repeat=n))
    factor=sp.simplify(num/den)
    error=sp.simplify(sum(sp.conjugate(observed(ix)-factor*reference(ix)) *
                          (observed(ix)-factor*reference(ix))
                          for ix in product(range(2),repeat=n)))
    return {'complex_ratio_to_chosen_candidate':str(factor),
            'exact_residual_squared_norm':str(error),
            'is_proportional_exact':bool(error==0),
            'ratio_is_a_physical_RGBeta_coupling_conversion':False}


def matchete_b_c_ratios(root: Path) -> dict:
    """Compare B/C native Matchete mixing CG with chosen candidate, not RGBeta."""
    gene=root/AUDITS/'matchete_generator_analysis.json'
    comp=root/AUDITS/'component_support_audit.json'
    cov=root/AUDITS/'cg_covariance_audit.json'
    if not all(p.is_file() for p in (gene,comp,cov)):
        return {'status':'source_CG_audits_missing','models_checked':0}
    try:
        g=json.loads(gene.read_text(encoding='utf-8-sig'))
        c=json.loads(comp.read_text(encoding='utf-8-sig'))
        v=json.loads(cov.read_text(encoding='utf-8-sig'))
        require(v.get('status')=='cg_generator_covariance_verified' and
                v.get('count')==16 and len(c['cases'])==16,
                'source CG audit incomplete')
        u=source_rotations(g)
        rows=[]
        for case in c['cases']:
            m=case['model']
            if int(m['dS1'])!=2 or int(m['dS2'])!=2:
                continue
            t=native_tensor(case['vertices']['T3MixCG'],FLAGS['T3MixCG'],u)
            observed=lambda ix: sp.simplify(t[ix])
            good=fit_tensor(repaired_22_mix,observed,4)
            bad=fit_tensor(literal_rgbeta_22_mix,observed,4)
            rows.append({'family':'B' if int(m['dF'])==1 else 'C',
                         'alpha':int(m['alpha']),
                         'matchete_vs_chosen_ward_candidate':good,
                         'matchete_vs_literal_source_delta_delta':bad})
        require(len(rows)==6,'not all six B/C Matchete tensor cases')
        return {'status':'conditional_Matchete_CG_comparisons_done',
                'models_checked':6,'models':rows,
                'native_RGBeta_numeric_axes_measured':False,
                'physical_RGBeta_coupling_ratios':None}
    except (ValueError,KeyError,TypeError,OSError) as exc:
        return {'status':'matchete_CG_comparison_blocked','reason':str(exc),
                'models_checked':0}


def compute(root: Path, formal_file: Path | None = None,
            registry_file: Path | None = None, *, validate_sources: bool = True) -> dict:
    result={'format':'b2-rgbeta-literal-ward-witness-v1',
            'physical_source_certified':False,
            'physical_crosspackage_complete':False,
            'native_RGBeta_component_axes_certified':False,
            'full_RGBeta_action_normalisation_verified':False,
            'numerical_RGE_or_matching_performed':False}
    result['ward_basis_witness']=invariant_identity()
    result['B_C_all_six_charge_witnesses']=all16_charge_warning()
    result['B_C_models_checked']=6
    if validate_sources:
        result['source_provenance']=source_provenance(root)
    else:
        result['source_provenance']={'status':'synthetic_only_no_repo_source_check'}
    if formal_file is not None:
        require(formal_file.is_file(),'native formal file missing')
        obj=json.loads(formal_file.read_text(encoding='utf-8-sig'))
        result['native_formal_evidence']=formal_evidence(obj)
        result['native_formal_evidence']['sha256']=hashlib.sha256(formal_file.read_bytes()).hexdigest()
    else:
        result['native_formal_evidence']={'status':'no_native_formal_file_provided'}
    if registry_file is not None:
        require(registry_file.is_file(),'native registry file missing')
        registry=json.loads(registry_file.read_text(encoding='utf-8-sig'))
        require(registry.get('format')=='b2-rgbeta-extended-native-registry-v1' and
                registry.get('physical_source_certified') is False and
                registry.get('status')=='ExtendedFormalRegistryProbeCompletedActionPending' and
                registry.get('numeric_RGBeta_components_certified') is False and
                registry.get('full_action_feynman_rules_verified') is False and
                len(registry.get('five_classes',[])) == 5,
                'native registry file is malformed, partial or overclaims')
        result['native_extended_registry']={
            'status':registry.get('status'),
            'sha256':hashlib.sha256(registry_file.read_bytes()).hexdigest(),
            'family_count':len(registry.get('five_classes',[])),
            'physical_source_certified':False}
    result['matchete_B_C_CG_comparison']=matchete_b_c_ratios(root) if validate_sources else {
        'status':'not_checked_no_local_matchete_audits'}
    result['status']='LITERAL_RGBETA_DOUBLE_DELTA_WARD_MISMATCH_CONDITIONAL'
    result['interpretation']='Mismatch assumes raw del/eps indices are literal fundamental/antifundamental components; this has not been verified for RGBeta.'
    result['source_bug_confirmed']=False
    result['source_change_authorised']=False
    result['crosspackage_yukawa_ratios']=None
    result['crosspackage_mass_ratio']=None
    result['crosspackage_mixing_ratio']=None
    return result
