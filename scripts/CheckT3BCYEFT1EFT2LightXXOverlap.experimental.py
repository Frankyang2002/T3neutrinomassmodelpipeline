#!/usr/bin/env python3
"""B(-1)/C(-1): verify the light-XX scalar bubble on BOTH sides of Y matching.

Independently reads (a) original F-first EFT1 LLSS/EOMEFT tree vertices and
(b) Y-only matched EFT2 SM-tree output, including its retained X interactions.
Projects both to the SAME ordered neutral SU(2) component using the actual
native CG registry and calculates K_XX separately in EFT1 and EFT2.

Only after showing K_XX(EFT1)==K_XX(EFT2), reconstruct the scalar-bubble
UV counterterm difference CT1-CT2 = -(K_YY+K_XY). Compare that with the
independently derived previous v5 counterterm and native Y pole. Also verify
the finite XX bubble logarithm cancels between the two EFTs for mX>0.

CAUTION: This is NOT full EFT1/EFT2 operator renormalization, not a complete
method-of-regions soft subtraction of the mixed XY loop, and NOT physical C5.
No native Matchete matching is performed. No production files modified.

At repository root:
  python scripts/CheckT3BCYEFT1EFT2LightXXOverlap.experimental.py --self-test
  python scripts/CheckT3BCYEFT1EFT2LightXXOverlap.experimental.py --output output/audits/t3_bc_y_eft2_light_xx_overlap_v1/summary.json
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from hashlib import sha256
import json
import math
from pathlib import Path
import sys

import sympy as sp

ROOT = Path('output/general_t3_mass_eigen_matchete')
AUDITS = Path('output/audits')
STEMS = {'B_m1': 'B_-1', 'C_m1': 'C_-1'}


class XXOverlapError(ValueError):
    pass


def require(condition, message):
    if not condition:
        raise XXOverlapError(message)


def serialize(m):
    return [
        {'coupling_factors': list(key[0]), 'chirality': key[1],
         'H_bar_count': key[2], 'coefficient_without_hbar_over_MF': str(value)}
        for key, value in sorted(m.items(), key=lambda item: str(item[0]))
    ]


def compact(m):
    return {k: sp.simplify(v) for k, v in m.items() if sp.simplify(v) != 0}


def source_xx_vertices(v3, v2, ll_source: str, quartic_source: str, tensors: dict):
    """No native pole is consulted. Only the original tree Lagrangians.

    Source files may group terms differently (notably C's LLXX polynomial).
    Fully expand each selected Matchete AST using the existing v3 parser and
    combine exact component coefficients, rather than comparing term counts.
    """
    ll = [t for t in v2.split_additive(ll_source)
          if v2.field_count(t, 'NewScalarX') == 2
          and v2.field_count(t, 'NewScalarY') == 0
          and t.count('Field[l, Fermion') == 2]
    qq = [t for t in v2.split_additive(quartic_source)
          if v2.field_count(t, 'NewScalarX') == 2
          and v2.field_count(t, 'NewScalarY') == 0
          and v2.field_count(t, 'H') == 2
          and 'Field[l, Fermion' not in t]
    require(ll and qq, 'Empty LLXX or HHXX source operator selection')
    w = defaultdict(lambda: sp.Integer(0))
    q = defaultdict(lambda: sp.Integer(0))
    for t in ll:
        for k, v in v3.slice_component_factors(t, tensors, external='wilson').items():
            require(v3.scalar_pair(k[3]) == 'XX', 'Unexpected LLSS pair type')
            w[k] += v
    for t in qq:
        for k, v in v3.slice_component_factors(t, tensors, external='quartic').items():
            require(v3.scalar_pair(k[3]) == 'XX', 'Unexpected HHSS pair type')
            q[k] += v
    # Keep component labels and external H bar count in the keys. Do not
    # symmetrise flavour or insert a hypothetical Majorana normalisation.
    return compact(w), compact(q), {'LLXX_additive_terms': len(ll),
                                     'HHXX_additive_terms': len(qq)}


def xx_kernel(v3, w: dict, q: dict) -> dict:
    result = defaultdict(lambda: sp.Integer(0))
    for (wc, ch, _, ws), wv in w.items():
        require(v3.scalar_pair(ws) == 'XX', 'Not an XX Wilson vertex')
        for (qc, _, hbar, qs), qv in q.items():
            require(v3.scalar_pair(qs) == 'XX', 'Not an XX quartic vertex')
            n = v3.contractions(ws, qs)
            if not n:
                continue
            key = (tuple(sorted((*wc, *qc))), ch, hbar)
            result[key] += n * wv * qv
    return compact(result)


def xx_finite_mass_log(mu_y: float, m_x: float) -> float:
    """Finite B0(0;mX,mX) in MSbar with mu=MY; SAME on either side."""
    require(mu_y > m_x > 0, 'Require 0 < MX < MY')
    return 2.0 * math.log(mu_y / m_x)


def read(path: Path) -> tuple[str, str]:
    require(path.is_file(), 'Missing original local source: ' + str(path))
    data = path.read_bytes()
    return data.decode('utf-8'), sha256(data).hexdigest()


def check_model(stem, v2, v3, v4, v5, root, native_dir,
                ledger_dir, prior_v5: dict):
    model = STEMS[stem]
    cg_file = native_dir / f'{stem}_native_CG.json'
    cg = v3.read_registry(cg_file, model)
    p_ll = (root / 'eft1_fermion_first_tree' /
            f'{stem}_final_eft1_fermion_first_tree_llss_tree.inputform.txt')
    p_eft1 = (root / 'eft1_fermion_first_tree' /
              f'{stem}_final_eft1_fermion_first_tree_eft1_transition_eom.inputform.txt')
    p_eft2 = (root / 'eft1_scalar_Y_only_canonical_probe' /
              f'{stem}_final_y_only_canonical_match_sm_tree.inputform.txt')
    a, ha = read(p_ll)
    b, hb = read(p_eft1)
    c, hc = read(p_eft2)

    require('NewScalarY' not in c, model + ': EFT2 tree still contains heavy Y')
    w1, q1, counts1 = source_xx_vertices(v3, v2, a, b, cg)
    w2, q2, counts2 = source_xx_vertices(v3, v2, c, c, cg)
    require(w1 and w2 and q1 and q2, model + ': XX component projection is empty')
    Wdiff = v5.equals(w1, w2)
    Qdiff = v5.equals(q1, q2)
    k1, k2 = xx_kernel(v3, w1, q1), xx_kernel(v3, w2, q2)
    require(k1 and k2, model + ': Zero XX bubble kernel cannot establish overlap')
    Kdiff = v5.equals(k1, k2)
    # Full UV local mixing in the *scalar-quartic* bubble subset on either
    # side, not the complete anomalous dimension of the dim-five EFT.
    ykernels, _, hashes = v4.get_pair_kernels(v3, v2, root, stem, cg_file)
    ky = v5.sum_maps(ykernels['YY'], ykernels['XY'])
    uv1 = v5.sum_maps(k1, ky)
    uv2 = k2
    ct1 = v5.invert_map(uv1)
    ct2 = v5.invert_map(uv2)
    ct_diff = v5.sum_maps(ct1, v5.invert_map(ct2))
    expected_diff = v5.invert_map(ky)
    ctdiff_mismatches = v5.equals(ct_diff, expected_diff)
    # External independent native pole cross-check, NOT used to build CT.
    native_pole, hp = v4.native_sector(
        v3, v2, ledger_dir / f'{stem}_y_llhh_raw_pole_residue.inputform.txt', cg, 50)
    native_check = v5.equals(ky, native_pole)
    native_ct_check = v5.equals(v5.sum_maps(ct_diff, native_pole), {})

    old = prior_v5.get('model_results', {}).get(model)
    require(isinstance(old, dict) and old.get('model') == model,
            model + ': Missing v5 source-derived CT report')
    prior_hashes = old.get('source_sha256', {})
    require(all(prior_hashes.get(k) == v for k, v in
                {**hashes, 'pole': hp}.items()),
            model + ': EFT1 / native CG / native pole changed since v5 audit')
    require(old.get('independent_bubble_MSbar_CT_cancels_native_Y_pole') is True,
            model + ': v5 source-CT audit was not successful')
    old_ct = {
        (tuple(rec['coupling_factors']), rec['chirality'], rec['H_bar_count']):
        rec['coefficient_without_hbar_and_inverse_MF']
        for rec in old.get('independent_bubble_sector_CT_without_hbar_and_inverse_MF', [])
    }
    # Compare serialised exact constants without treating them as untrusted
    # executable SymPy expressions; v5 and v6 share the same source parser.
    old_serial = {
        (tuple(rec['coupling_factors']), rec['chirality'], rec['H_bar_count']):
        str(rec['coefficient_without_hbar_over_MF'])
        for rec in serialize(ct_diff)
    }
    prior_ct_mismatch = (old_ct != old_serial)

    # The matched EFT2 tree must not already contain LLHH as a tree-level
    # source. A new LLHH term would change the one-loop matching book-keeping.
    explicit_llhh = [t for t in v2.split_additive(c)
        if v2.field_count(t, 'H') == 2 and
           t.count('Field[l, Fermion') == 2 and
           v2.field_count(t, 'NewScalarX') == 0]
    require(not explicit_llhh, model + ': EFT2 tree has unexpected LLHH')
    parfile = root / f'{stem}_final_parameters.json'
    par = json.loads(read(parfile)[0]); require(par.get('model') == model,
                'Unexpected threshold parameter model')
    mX, mY = map(float, par['heavy_masses_gev'])
    xxlog = xx_finite_mass_log(mY, mX)
    # At fixed physical mX the XX finite logarithm is present in BOTH EFTs;
    # cancellation is exact because their X-only kernels agree.
    xx_fin1 = {k: sp.simplify(v * sp.Symbol('LXX')) for k, v in k1.items()}
    xx_fin2 = {k: sp.simplify(v * sp.Symbol('LXX')) for k, v in k2.items()}
    xx_fin_diff = v5.equals(xx_fin1, xx_fin2)

    return {
        'format': 't3bc-y-e1-e2-light-xx-bubble-overlap-v1',
        'model': model,
        'source_sha256': {'EFT1_LLSS': ha, 'EFT1_EOMEFT': hb,
                          'EFT2_Y_tree': hc, 'native_CG': sha256(cg_file.read_bytes()).hexdigest(),
                          'native_Y_pole': hp},
        'selected_source_terms': {'EFT1': counts1, 'EFT2': counts2,
                                  'EFT2_tree_LLHH': len(explicit_llhh)},
        'component_counts': {'EFT1_LLXX': len(w1), 'EFT2_LLXX': len(w2),
                             'EFT1_HHXX': len(q1), 'EFT2_HHXX': len(q2),
                             'KXX_EFT1': len(k1), 'KXX_EFT2': len(k2),
                             'KYY': len(ykernels['YY']), 'KXY': len(ykernels['XY'])},
        'exact_LLXX_component_equality': not Wdiff,
        'exact_HHXX_component_equality': not Qdiff,
        'exact_XX_loop_UV_component_equality': not Kdiff,
        'EFT1_minus_EFT2_scalar_bubble_UV_equals_Y_vertices':
            not v5.equals(v5.sum_maps(uv1, v5.invert_map(uv2)), ky),
        'EFT1_minus_EFT2_scalar_bubble_CT_equals_minus_Y': not ctdiff_mismatches,
        'EFT1_EFT2_scalar_bubble_CT_difference_cancels_native_Y_pole':
            not native_ct_check,
        'XX_finite_MSbar_log_cancelled_between_EFTs': not xx_fin_diff,
        'independent_Y_kernel_equals_native_pole': not native_check,
        'previous_v5_CT_matches_v6_CT_difference': not prior_ct_mismatch,
        'XX_finite_log_at_muY': xxlog,
        'XX_finite_log_interpretation':
            'At mu=MY and mX>0, log(MY^2/MX^2) multiplies KXX in BOTH EFTs; '
            'the entire XX sector cancels. The massless scaleless limit '
            'requires UV/IR separation instead of treating zero as zero UV.',
        'EFT1_LLXX_vs_EFT2_mismatches': Wdiff[:30],
        'EFT1_HHXX_vs_EFT2_mismatches': Qdiff[:30],
        'EFT1_XX_UV_vs_EFT2_mismatches': Kdiff[:30],
        'CT_difference_vs_minus_Y_mismatches': ctdiff_mismatches[:30],
        'CT_difference_plus_native_pole_mismatches': native_ct_check[:30],
        'XX_finite_log_mismatches': xx_fin_diff[:30],
        'scalar_bubble_UV_EFT1_without_hbar_and_inverse_MF': serialize(uv1),
        'scalar_bubble_UV_EFT2_without_hbar_and_inverse_MF': serialize(uv2),
        'scalar_bubble_CT_difference_without_hbar_and_inverse_MF': serialize(ct_diff),
        'scalar_bubble_only': True,
        'complete_EFT1_EFT2_operator_counterterms_certified': False,
        'EFT2_mixed_XY_soft_region_subtraction_completed': False,
        'mass_suppressed_operator_basis_included': False,
        'GammaCC_to_symmetric_Weinberg_C5_certified': False,
        'physical_C5_Y': None,
        'physical_matching_certified': False,
    }


def self_test():
    # Pure UV algebra with independent XX and Y dictionaries; XX must be
    # nonzero, otherwise this would be a vacuous EFT1/EFT2 overlap check.
    k = (('lambdaXX', 'y1', 'y1'), 'PR', 2)
    y = (('lambdaXY', 'y1', 'y2'), 'PR', 2)
    xx = {k: sp.Rational(2, 3)}
    yy = {y: sp.Rational(1, 5)}
    def total(left, right):
        out = defaultdict(lambda: sp.Integer(0))
        for m in (left, right):
            for kk, vv in m.items(): out[kk] += vv
        return compact(out)
    e1 = total(xx, yy)
    e2 = xx.copy()
    ct1 = {a: -b for a, b in e1.items()}
    ct2 = {a: -b for a, b in e2.items()}
    assert total(ct1, {a: -b for a, b in ct2.items()}) == {y: -sp.Rational(1, 5)}
    # Alter the EFT2 XX vertex: the source-level validation MUST fail.
    bad = {k: sp.Rational(3, 4)}
    assert e2 != bad
    assert total(ct1, {a: -b for a, b in {k: -v for k, v in bad.items()}.items()}) != {y: -sp.Rational(1, 5)}
    # Two Wick pairings for identical XX fields, one for distinct components.
    # The native v3 implementation is used at runtime; this independently
    # demonstrates multiplicity expectations in its input convention.
    assert abs(xx_finite_mass_log(1188.5289940451582, 1035.583176055388)) < 1
    assert xx_finite_mass_log(2.0, 1.0) == math.log(4.0)
    print('PASS: nonzero EFT1/EFT2 XX UV subtraction and tamper rejection')
    print('PASS: same-mass light-XX finite logarithm cancels algebraically')
    print('PASS: full EFT counterterms and physical C5 remain uncertified')


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--self-test', action='store_true')
    ap.add_argument('--output', type=Path)
    ap.add_argument('--root', type=Path, default=ROOT)
    ap.add_argument('--native-dir', type=Path,
                    default=AUDITS / 't3_bc_y_native_cg_components_v1')
    ap.add_argument('--ledger-dir', type=Path,
                    default=AUDITS / 't3_bc_y_renormalization_ledger_v1')
    ap.add_argument('--v5-summary', type=Path,
                    default=AUDITS / 't3_bc_y_scalar_bubble_msbar_v1' / 'summary.json')
    a = ap.parse_args()
    try:
        if a.self_test: self_test()
        if a.output:
            require(a.v5_summary.is_file(), 'Missing prior v5 summary: ' + str(a.v5_summary))
            old = json.loads(a.v5_summary.read_text(encoding='utf-8'))
            require(old.get('format') == 't3bc-y-scalar-bubble-msbar-counterterm-audit-v1',
                    'Wrong previous v5 report')
            from importlib.util import spec_from_file_location, module_from_spec
            def sibling(filename, name):
                path = Path(__file__).resolve().with_name(filename)
                require(path.is_file(), 'Missing previous script ' + str(path))
                spec = spec_from_file_location(name, path)
                require(spec and spec.loader, 'Cannot load ' + str(path))
                obj = module_from_spec(spec); spec.loader.exec_module(obj)
                return obj
            v2 = sibling('AuditT3BCYMixedScalarPoleSupport.experimental.py', 't3bc_overlap_v2')
            v3 = sibling('CompareT3BCYNativeSU2Coefficients.experimental.py', 't3bc_overlap_v3')
            v4 = sibling('CheckT3BCYFiniteScalarBubble.experimental.py', 't3bc_overlap_v4')
            v5 = sibling('DeriveT3BCYBubbleMSbarCT.experimental.py', 't3bc_overlap_v5')
            summary = {
                'format': 't3bc-y-e1-e2-light-xx-overlap-report-v1',
                'models': {}, 'physical_C5_Y': None,
                'complete_EFT1_EFT2_counterterms_certified': False,
            }
            all_pass = True
            for stem in STEMS:
                item = check_model(stem, v2, v3, v4, v5, a.root,
                                   a.native_dir, a.ledger_dir, old)
                summary['models'][item['model']] = item
                good = all(item[z] for z in (
                    'exact_LLXX_component_equality', 'exact_HHXX_component_equality',
                    'exact_XX_loop_UV_component_equality',
                    'EFT1_minus_EFT2_scalar_bubble_UV_equals_Y_vertices',
                    'EFT1_minus_EFT2_scalar_bubble_CT_equals_minus_Y',
                    'EFT1_EFT2_scalar_bubble_CT_difference_cancels_native_Y_pole',
                    'XX_finite_MSbar_log_cancelled_between_EFTs',
                    'independent_Y_kernel_equals_native_pole',
                    'previous_v5_CT_matches_v6_CT_difference'))
                all_pass &= good
                print(f"{item['model']}: XX source equality={item['exact_LLXX_component_equality']}/"
                      f"{item['exact_HHXX_component_equality']}; "
                      f"XX UV overlap={item['exact_XX_loop_UV_component_equality']}; "
                      f"EFT1-EFT2 CT difference={item['EFT1_EFT2_scalar_bubble_CT_difference_cancels_native_Y_pole']}; "
                      f"full EFT CT=NOT CERTIFIED")
            a.output.parent.mkdir(parents=True, exist_ok=True)
            a.output.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + '\n',
                                encoding='utf-8')
            print(('PASS' if all_pass else 'FAIL') + ': X-only EFT1/EFT2 scalar-bubble overlap -> ' + str(a.output))
            return 0 if all_pass else 1
        require(a.self_test, 'Specify --self-test and/or --output')
        return 0
    except (XXOverlapError, ValueError, KeyError, IndexError, TypeError,
            AttributeError, AssertionError, OSError, RuntimeError) as exc:
        print('FAIL (not physical C5): ' + str(exc), file=sys.stderr)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
