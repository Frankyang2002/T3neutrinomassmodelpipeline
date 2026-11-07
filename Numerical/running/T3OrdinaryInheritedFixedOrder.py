"""Source-audited inherited F-threshold Weinberg term for D(+2), A(-4).

Basis of this restricted numerical calculation:
  (1) Native scalar-only one-loop Matchete hard C5 template, verified in Batch 3.
  (2) Archived general T3 F-first one-loop EFT1 InputForm for the same model.
  (3) Three nondegenerate diagonal heavy-fermion masses at muF.

The legacy F-first source has six pure-SM 2L2H terms: three holomorphic
finite/epsilon/log terms, and their three hermitian-conjugate counterparts.
The holomorphic coefficient is +2*hbar*lambdaT3/(sqrt(3)*MF), while the
independent scalar hard template has -2*hbar*lambdaT3/(sqrt(3)*MF).
The equality of these pole residues is audited at common fixed-order
couplings; an evolved-coupling combination is reported separately and is
NOT certified pole-free at finite resummed running.

No new Matchete matching is performed. Source archives were produced by an
older full16 symbolic matching pipeline; this does not replace a native
all-three-mass one-loop first-threshold calculation or certify SU(2)/spinor
normalization and full SMEFT matching.
"""
from __future__ import annotations

from hashlib import sha256
import math
from pathlib import Path
from typing import Any

import numpy as np

from Numerical.running.T3OrdinaryHardFlavor import (
    MODELS, TAGS, HardFlavorError, decode_complex, encode_complex,
    file_sha, mass_kernel, metrics, read_json, require,
    running_parameters, source_pole_matrix, verify_hard_stage,
)

ARCHIVE_NAMES = {
    'D_+2': 'T3_D_alpha_p2',
    'A_-4': 'T3_A_alpha_m4',
}
ARCHIVE_FILE = 'data/eft1_after_F_one_loop_inputform.txt'
EXPECTED_ORIENTATION = ('holomorphic', 'anti_holomorphic')


def top_level_terms(source: str) -> list[str]:
    """Split Matchete InputForm Plus without splitting nested sums/products."""
    require(isinstance(source, str) and len(source) > 100, 'One-loop source empty')
    brackets: list[str] = []
    pairs = {')': '(', ']': '[', '}': '{'}
    terms = []
    start = 0
    for idx, c in enumerate(source):
        if c in '([{':
            brackets.append(c)
        elif c in ')]}':
            require(bool(brackets) and brackets.pop() == pairs[c],
                    'Malformed Matchete InputForm brackets')
        elif c in '+-' and not brackets and idx > start:
            part = source[start:idx].strip()
            if part:
                terms.append(part)
                start = idx
    require(not brackets, 'Unbalanced Matchete InputForm brackets')
    final = source[start:].strip()
    if final:
        terms.append(final)
    return terms


def _compact(text: str) -> str:
    return ''.join(text.split())


def audit_inherited_source(path: Path, model: str) -> dict[str, Any]:
    """Require exact +2/sqrt(3) F-threshold finite/pole/log operator triplets.

    This is a term-structure audit; no symbolic evaluation/Matchete Match.
    The archived F-first term has a *single unindexed MF*. Only the
    explicitly diagonal mass-eigenstate lift is used downstream.
    """
    require(model in MODELS, f'Unsupported model {model}')
    path = Path(path)
    require(path.is_file(), f'Missing inherited one-loop source: {path}')
    content = path.read_text(encoding='utf-8-sig')
    parts = top_level_terms(content)
    weinberg = [t for t in parts
                if t.count('Field[l,') == 2 and t.count('Field[H,') == 2
                and 'hbar' in t
                and not any(f'Field[{n},' in t for n in ('NewScalar1', 'NewScalar2', 'NewFermion'))]
    require(len(weinberg) == 6,
            f'{model}: expected exactly six pure-SM two-lepton/two-Higgs terms, found {len(weinberg)}')
    groups: dict[str, dict[str, str]] = {k: {} for k in EXPECTED_ORIENTATION}
    log_token = 'Log[μbar2/Coupling[MF,{},0]^2]*'
    for raw in weinberg:
        t = _compact(raw)
        require(t.startswith('+') or not t.startswith('-'),
                f'{model}: negative inherited Weinberg term')
        t = t.lstrip('+')
        barred = all(f'Bar[Coupling[{n},' in t for n in ('y1', 'y2'))
        direct = all(f'Coupling[{n},' in t for n in ('y1', 'y2'))
        if barred and 'Coupling[lambdaT3,{},0]' in t and 'Bar[Coupling[lambdaT3' not in t:
            orientation = 'holomorphic'
        elif direct and 'Bar[Coupling[lambdaT3,{},0]]' in t and not barred:
            orientation = 'anti_holomorphic'
        else:
            raise HardFlavorError(f'{model}: ambiguous Yukawa/complex-conjugation sector')
        require(t.count('hbar') == 1 and 'GammaCC' in t and
                'Coupling[MF,{},0]' in t and
                'Coupling[MF,{Index[' not in t and
                '2*hbar*' in t and 'Sqrt[3]' in t,
                f'{model}: unsupported inherited source normalization or heavy-mass structure')
        has_eps = 'ϵ' in t
        has_log = log_token in t
        require(not (has_eps and has_log), f'{model}: pole-log mixed term')
        category = 'pole' if has_eps else 'log' if has_log else 'finite'
        denom = ('/(Sqrt[3]*ϵ*Coupling[MF,{},0])' if has_eps
                 else '/(Sqrt[3]*Coupling[MF,{},0])')
        require(t.endswith(denom), f'{model}: unexpected inherited Weinberg denominator {category}')
        require(category not in groups[orientation], f'{model}: duplicate {orientation} {category}')
        groups[orientation][category] = t
    for key, group in groups.items():
        require(set(group) == {'finite', 'pole', 'log'},
                f'{model}: incomplete {key} finite/pole/log triplet')
        normalized = {k: v.replace('ϵ*', '').replace(log_token, '')
                      for k, v in group.items()}
        require(len(set(normalized.values())) == 1,
                f'{model}: {key} finite/pole/log coefficients differ algebraically')
        require(normalized['finite'].count('Coupling[y1,') == 1 and
                normalized['finite'].count('Coupling[y2,') == 1 and
                'Coupling[lambdaT3,{},0]' in normalized['finite'],
                f'{model}: source has extra Yukawa/λ factors')
    return {
        'status': 'ArchivedFOneLoopWeinbergTripletsVerified',
        'model': model,
        'archive_path': str(path),
        'archive_sha256': file_sha(path),
        'inputform_additive_term_count': len(parts),
        'pure_sm_weinberg_term_count': len(weinberg),
        'holomorphic_triplet': ['finite', 'pole', 'log'],
        'anti_holomorphic_triplet': ['finite', 'pole', 'log'],
        'per_eigenvalue_inherited_pole_coefficient': '+2*hbar*lambdaT3/(sqrt(3)*MF) times barred y1*y2',
        'per_eigenvalue_hard_pole_coefficient': '-2*hbar*lambdaT3/(sqrt(3)*MF) times barred y1*y2',
        'pole_cancellation_at_common_couplings': True,
        'legacy_unindexed_MF_requires_diagonal_heavy_mass_lift': True,
        'native_first_threshold_three_mass_one_loop_matching_recomputed': False,
        'complete_counterterm_and_scheme_certified': False,
    }


def _sum_parts(parts: list[np.ndarray]) -> np.ndarray:
    return sum(parts, np.zeros((3, 3), complex))


def matching_kernel_at_heavy_mass(m1sq: float, m2sq: float, mf: float) -> float:
    """Exact algebraic combination of hard finite + inherited finite + log.

      [1 + weighted_log(muS^2/ms^2)] - 2log(muS/muF)
       - [1 + log(muF^2/MF^2)] = weighted_log(MF^2/ms^2).
    """
    return mass_kernel(m1sq, m2sq, mf) - 1.0


def assemble_fixed_order(*, model: str, inherited_path: Path,
                         hard_stage_folder: Path, running_doc: dict,
                         batch3_result: dict, vev_gev: float = 246.22) -> dict:
    """Compute the source-backed inherited Weinberg sector without rematching.

    The fixed-coupling result is the controlled O(hbar) number; the
    muS-evolved combination is a separately flagged higher-order diagnostic.
    """
    require(model in MODELS, f'Unsupported model {model}')
    hard_source = verify_hard_stage(hard_stage_folder, model)
    inherited_source = audit_inherited_source(inherited_path, model)
    p = running_parameters(running_doc, model)
    require(batch3_result.get('model') == model and
            batch3_result.get('status') == 'HardThreeFlavorAndDirectLogComputed_NotFullC5',
            f'{model}: missing verified Batch 3 numerical result')
    require(batch3_result.get('source_template', {}).get('source_expression_sha256_normalized')
            == hard_source['source_expression_sha256_normalized'],
            f'{model}: Batch 3 hard expression differs')
    require(batch3_result.get('source_template', {}).get('source_hard_report_sha256')
            == hard_source['source_hard_report_sha256'],
            f'{model}: Batch 3 hard-stage report differs')
    require(batch3_result.get('source_template', {}).get('source_provenance_sha256')
            == hard_source['source_provenance_sha256'],
            f'{model}: Batch 3 provenance differs')
    require(math.isfinite(vev_gev) and vev_gev > 0, 'Positive finite vev required')
    require(abs(float(batch3_result.get('vev_gev', -1)) - vev_gev) < 1e-10,
            f'{model}: vev differs from saved Batch 3')

    ps, evolved_poles = source_pole_matrix(p['y1'],p['y2'],p['mf'],p['lam'])
    pf, frozen_poles = source_pole_matrix(p['y1'],p['y2'],p['mf'],p['lam_f'])
    archived_pole = decode_complex(batch3_result['raw_unsubtracted_hard_1_over_epsilon_residue_gev_inverse'],
                                  'Batch3 scalar pole')
    archived_subtotal = decode_complex(batch3_result['c5_hard_finite_plus_direct_log_subset_gev_inverse'],
                                      'Batch3 hard+log subtotal')
    require(np.linalg.norm(archived_pole-ps) <= max(float(np.linalg.norm(ps))*1e-12,1e-28),
            f'{model}: saved Batch3 pole no longer matches running inputs')
    expected_evolved = ps*(mass_kernel(p['m1sq'],p['m2sq'],p['mu_s'])
                           -2*math.log(p['mu_s']/p['mu_f']))
    require(np.linalg.norm(archived_subtotal-expected_evolved)
            <= max(float(np.linalg.norm(expected_evolved))*1e-12,1e-28),
            f'{model}: saved Batch3 subtotal changed; refusing double count')

    # The source's holomorphic inherited pole is minus the hard one in the
    # SAME renormalisable parameter convention. It is not an ad-hoc counterterm.
    inherited_pole = [-z for z in frozen_poles]
    inherited_finite = [(-z)*(1.0+2*math.log(p['mu_f']/m))
                        for z, m in zip(frozen_poles,p['mf'])]
    inh = _sum_parts(inherited_finite)
    fixed_hard = pf * mass_kernel(p['m1sq_f'],p['m2sq_f'],p['mu_s'])
    fixed_log = (-2*math.log(p['mu_s']/p['mu_f']))*pf
    fixed_subtotal = fixed_hard + fixed_log
    c5_fixed = fixed_subtotal + inh
    shortcut = _sum_parts([z*matching_kernel_at_heavy_mass(p['m1sq_f'],p['m2sq_f'],m)
                           for z,m in zip(frozen_poles,p['mf'])])
    denominator = max(float(np.linalg.norm(c5_fixed)),1e-40)
    identity_residual = float(np.linalg.norm(shortcut-c5_fixed)/denominator)
    require(identity_residual < 3e-12,
            f'{model}: fixed-order matching-log identity did not close')
    pole_residue = pf + _sum_parts(inherited_pole)
    pole_cancel_residual = float(np.linalg.norm(pole_residue) /
                                 max(float(np.linalg.norm(pf)),1e-40))
    require(pole_cancel_residual < 1e-13,
            f'{model}: common-coupling hard/inherited poles did not cancel')
    # Source coefficients are exact +/-2/sqrt(3), but numerical y/MF are
    # floating point. The matrix-level equality is verified separately.
    test_scale = p['mu_s']*math.exp(.046875)
    shifted = (pf*mass_kernel(p['m1sq_f'],p['m2sq_f'],test_scale)
               -2*pf*math.log(test_scale/p['mu_f']) + inh)
    scale_residual = float(np.linalg.norm(shifted-c5_fixed)/denominator)
    require(scale_residual < 3e-12,
            f'{model}: common-coupling muS variation failed')

    # This combines endpoints with different couplings; the unequal pole
    # residue is NOT 'renormalized away'. It estimates the effect of RG
    # improvement and missing higher-order consistency in this truncation.
    c5_evolved = archived_subtotal + inh
    evolved_pole_residual = ps + _sum_parts(inherited_pole)
    evolved_relative = float(np.linalg.norm(evolved_pole_residual)/
                             max(float(np.linalg.norm(ps)),float(np.linalg.norm(pf)),1e-40))

    return {
        'format':'t3-ordinary-inherited-fixed-order-source-audit-v1',
        'model':model,
        'status':'InheritedWeinbergAndCommonCouplingPoleCancellationComputed_Conditional',
        'basis':'holomorphic Weinberg C5 with symmetric lepton-index projection',
        'source_hard':hard_source,
        'source_legacy_fermion_one_loop':inherited_source,
        'matching_scales_gev':{'muF':p['mu_f'],'muS':p['mu_s']},
        'fermion_masses_gev':p['mf'].tolist(),
        'fixed_order_common_renormalisable_couplings':{
            'location':'muF',
            'lambdaT3':{'real':p['lam_f'].real,'imag':p['lam_f'].imag},
            'mS1sq_gev2':p['m1sq_f'],'mS2sq_gev2':p['m2sq_f'],
            'hbar':'1/(16*pi^2)',
            'matching_scheme':'source InputForm finite terms with muBarF=muF, muBarS=muS',
            'common_coupling_vs_RG_improved':'One-loop fixed order; scalar beta effects formally higher order',
        },
        'common_hard_pole_residue_gev_inverse':encode_complex(pf),
        'common_inherited_pole_residue_gev_inverse':encode_complex(_sum_parts(inherited_pole)),
        'common_total_pole_residue_gev_inverse':encode_complex(pole_residue),
        'common_hard_finite_gev_inverse':encode_complex(fixed_hard),
        'common_direct_eft1_log_gev_inverse':encode_complex(fixed_log),
        'inherited_F_one_loop_finite_gev_inverse':encode_complex(inh),
        'common_fixed_order_inherited_inclusive_c5_gev_inverse':encode_complex(c5_fixed),
        'common_fixed_order_simplified_kernel_gev_inverse':encode_complex(shortcut),
        'common_fixed_order_masses':metrics(c5_fixed,vev_gev),
        'RG_improved_mixed_endpoint_c5_gev_inverse_not_pole_free':encode_complex(c5_evolved),
        'RG_improved_mixed_endpoint_masses_diagnostic':metrics(c5_evolved,vev_gev),
        'RG_improved_mixed_endpoint_uncancelled_pole_residue_gev_inverse':encode_complex(evolved_pole_residual),
        'per_heavy_eigenstate':[{
            'r':r+1, 'MF_gev':float(p['mf'][r]),
            'common_hard_pole':encode_complex(frozen_poles[r]),
            'common_inherited_pole':encode_complex(inherited_pole[r]),
            'inherited_F_finite':encode_complex(inherited_finite[r]),
            'fixed_order_total':encode_complex(frozen_poles[r]*matching_kernel_at_heavy_mass(
                p['m1sq_f'],p['m2sq_f'],p['mf'][r])),
        } for r in range(3)],
        'checks':{
            'native_hard_scalar_matchete_template_verified':True,
            'archived_fermion_one_loop_weinberg_triplet_verified':True,
            'hard_and_inherited_pole_opposite_at_common_couplings':True,
            'common_pole_cancellation_relative_residual':pole_cancel_residual,
            'fixed_order_log_kernel_identity_relative_residual':identity_residual,
            'common_scalar_matching_scale_variation_relative_residual':scale_residual,
            'mixed_endpoint_uncancelled_pole_relative_residual':evolved_relative,
            'three_nondegenerate_eigenvalue_sum':True,
            'no_Matchete_matching_recomputed':True,
            'archived_LLSS_not_added_again':True,
            'independent_native_F_one_loop_three_mass_recheck':False,
            'all_first_threshold_operator_classes_EOM_reduced':False,
            'full_MSbar_renormalisation_certified':False,
            'physical_SU2_spinor_normalisation_certified':False,
            'full_general_sequential_C5_certified':False,
        },
        'interpretation':(
            'The exact ±2/sqrt(3) Weinberg poles cancel at common fixed-order '
            'couplings and the inherited F-threshold finite+logarithmic term '
            'is included for each heavy-mass eigenstate. This is a source-backed '
            'restricted holomorphic Weinberg-sector calculation, not a native '
            'independent NF=3 F-threshold one-loop re-match, an all-operator '
            'SMEFT/EOM proof, or a physically normalized full-model C5. '
            'The evolved-scalar endpoint remains a distinct non-pole-free diagnostic.'),
    }
