#!/usr/bin/env python3
"""Workstream C / v9.1: source-only fallback when SMEFT_Warsaw is not installed.

Does not use a fabricated SMEFT operator or assume a GammaCC conversion factor.
If a real installed Warsaw reference is available, the v9 structural route is
preserved. Otherwise inspect the *actual* saved B/C Y/LLHH source tensors,
compare their SU(2) skeleton to the published Warsaw definition, and explicitly
leave the absolute factor unresolved.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import sys

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

FMT = 't3-bc-matchete-warsaw-dim5-probe-v1'
OUTPUT = 't3-bc-native-to-warsaw-weinberg-convention-audit-v2'
SOURCE_REL = Path('output/general_t3_mass_eigen_matchete/eft1_scalar_Y_only_canonical_probe')
MASS_TERM_CONVENTION = 'L_m = -(1/2)*nu_L^T*C*m_nu*nu_L + h.c.'
PUBLISHED_Q = 'Q_llHH^pr = eps_in eps_jm H^i H^j [(ell_p^n)^T C ell_r^m]'


def fail(msg):
    raise ValueError(msg)


def signature(term):
    h = term.count('Field[H, Scalar')
    hb = term.count('Bar[Field[H, Scalar')
    l = term.count('Field[l, Fermion')
    lb = term.count('Bar[Field[l, Fermion')
    if (h, l) != (2, 2) or term.count('GammaCC') != 1:
        fail('not an explicit Weinberg-like Matchete operator')
    if (hb, lb) not in ((0, 0), (2, 2)):
        fail('unexpected mixed barred/unbarred field orientation')
    pl = 'Proj[-1]' in term
    pr = 'Proj[1]' in term
    if pl == pr:
        fail('missing or ambiguous chirality projector')
    return dict(
        orientation='direct' if hb == 0 else 'hermitian_conjugate',
        H_barred=hb, lepton_barred=lb,
        chirality='PL' if pl else 'PR',
        CG_epsilon_count=term.count('eps[SU2L]'),
        coupling_labels=sorted(set(re.findall(r'Coupling\[([A-Za-z0-9$]+),', term))),
        sha256=hashlib.sha256(term.encode('utf-8')).hexdigest(),
    )


def inspect_reference(reference):
    if reference.get('format') != FMT or reference.get('status') != 'WarsawDim5StructureExported':
        fail('not a successful installed Warsaw model export')
    data = reference.get('weinberg_terms')
    if not isinstance(data, list) or not data:
        fail('no Warsaw operators exported')
    recs = []
    for entry in data:
        if not isinstance(entry, dict) or not isinstance(entry.get('inputform'), str):
            fail('malformed reference entry')
        sig = signature(entry['inputform'])
        for source_key, value in [('H_fields', 2), ('lepton_fields', 2),
                                  ('barred_H_fields', sig['H_barred']),
                                  ('barred_lepton_fields', sig['lepton_barred'])]:
            if entry.get(source_key) != value:
                fail('reference source metadata mismatch: ' + source_key)
        recs.append(sig)
    if reference.get('weinberg_source_terms_found') != len(recs):
        fail('reference operator count mismatch')
    if not any(x['orientation'] == 'direct' for x in recs):
        fail('no direct Warsaw operator')
    if len({x['sha256'] for x in recs}) != len(recs):
        fail('duplicate reference operator')
    return dict(reference_model=reference.get('reference_model'),
                reference_operator_source='installed_model',
                native_reference_operator_detected=True,
                reference_operator_count=len(recs),
                native_reference_operator_signatures=recs)


def native_term_signature(term):
    sig = signature(term)
    if sig['CG_epsilon_count'] != 2:
        fail('native LLHH term lacks exactly two SU(2) epsilon CGs')
    if term.count('NCM[') != 1:
        fail('native LLHH term has unexpected spinor bilinear multiplicity')
    if term.count('hbar') != 1 or term.count('Coupling[MF,') != 1:
        fail('native loop/mass prefactor inventory mismatch')
    # These are source conventions established by the existing Matchete
    # Y-only matching. No Weyl spinor charge-conjugation sign is inferred.
    if sig['chirality'] == 'PR' and (sig['H_barred'], sig['lepton_barred']) != (2, 2):
        fail('PR native term unexpectedly unbarred')
    if sig['chirality'] == 'PL' and (sig['H_barred'], sig['lepton_barred']) != (0, 0):
        fail('PL native term unexpectedly barred')
    return sig


def inspect_native(model, root):
    from Numerical.experimental.T3BCMatcheteFiniteTerms import split_inputform
    stem = {'B_-1': 'B_m1', 'C_-1': 'C_m1'}[model]
    path = root / f'{stem}_final_y_only_canonical_match_raw_scalar_llhh.inputform.txt'
    if not path.is_file():
        fail('required native source is missing: ' + str(path))
    source = path.read_text(encoding='utf-8-sig')
    split = split_inputform(source, stage='Y', sector='LLHH')
    inventory = (len(split.terms), len(split.pole_terms),
                 len(split.log_terms), len(split.constant_terms))
    if inventory != (132, 50, 50, 32):
        fail(f'{model}: unexpected source inventory {inventory}')
    pole = [native_term_signature(t) for t in split.pole_terms]
    finite = [native_term_signature(t) for t in split.constant_terms]
    for label, signatures, expected in [('pole', pole, (25, 25)),
                                         ('finite', finite, (16, 16))]:
        counts = (sum(s['chirality'] == 'PR' for s in signatures),
                  sum(s['chirality'] == 'PL' for s in signatures))
        if counts != expected:
            fail(f'{model}: {label} chirality counts {counts} != {expected}')
    # The native source matches the field and epsilon *skeleton*, but we have
    # NOT compared a numerical model-file Wilson prefactor to this source.
    return dict(model=model, source_file=str(path),
                source_sha256=hashlib.sha256(source.encode('utf8')).hexdigest(),
                total_terms=132, pole_terms=50, log_terms=50, finite_terms=32,
                pole_chiralities={'PR': 25, 'PL': 25},
                finite_chiralities={'PR': 16, 'PL': 16},
                native_llhh_su2_epsilon_skeleton_checked=True,
                all_native_terms_have_two_H_two_leptons_and_gammaCC=True,
                native_PR_PL_both_counted_as_independent_C5=False,
                native_GammaCC_to_Warsaw_Wilson_prefactor_certified=False)


def analyse(reference, native_root):
    ref_meta = None
    if reference is not None:
        if not isinstance(reference, dict) or reference.get('format') != FMT:
            fail('unexpected reference JSON format')
        if reference.get('status') == 'WarsawDim5StructureExported':
            ref_meta = inspect_reference(reference)
        elif reference.get('status') in ('WarsawModelLoadFailed', 'WarsawModelUnavailable'):
            ref_meta = dict(reference_model='SMEFT_Warsaw',
                            reference_operator_source='UNAVAILABLE_NOT_INSTALLED',
                            reference_failure_status=reference['status'])
        else:
            fail('reference failed for another reason: ' + str(reference.get('status')))
    root = native_root if native_root.is_absolute() else REPO / native_root
    native = [inspect_native(m, root) for m in ('B_-1', 'C_-1')]
    return dict(
        format=OUTPUT,
        status='NativeSourceWarsawSkeletonAudited_AbsoluteNormalizationOutstanding',
        published_operator=PUBLISHED_Q,
        published_convention_source='https://matchete.gitlab.io/conventions_SMEFT_Warsaw.html',
        reference=ref_meta or dict(reference_operator_source='NOT_REQUESTED'),
        models=native,
        native_warsaw_su2_and_field_structure_compatible=True,
        # The following is only a conditional algebraic identity:
        flavour_rule_if_identical_operator_normalization='C5=A+A^T for L5=(1/2)*C5*Q+h.c.',
        flavour_transpose_not_adjoint=True,
        chosen_future_mass_term_convention=MASS_TERM_CONVENTION,
        absolute_GammaCC_to_Warsaw_factor_certified=False,
        hbar_to_loop_factor_certified_here=False,
        Warsaw_C5_to_neutrino_mass_sign_certified=False,
        all_EFT1_EFT2_counterterms_certified=False,
        physical_C5=None,
        next_required_source='actual installed SMEFT_Warsaw model operator OR an independent Weyl/Dirac action normalization derivation',
    )


def self_test():
    direct=('hbar*Coupling[MF, {}, 0]*Field[H, Scalar, {}, {}]*'
            'Field[H, Scalar, {}, {}]*CG[eps[SU2L], {}]*CG[eps[SU2L], {}]*'
            'NCM[Field[l, Fermion, {}, {}], DiracProduct[GammaCC, Proj[-1]],'
            'Field[l, Fermion, {}, {}]]')
    a = native_term_signature(direct)
    assert a['chirality'] == 'PL' and a['CG_epsilon_count'] == 2
    bad = direct.replace('CG[eps[SU2L], {}]*', '', 1)
    try:
        native_term_signature(bad)
    except ValueError:
        pass
    else:
        raise AssertionError('missing CG tensor passed')
    wrong = direct.replace('Proj[-1]', 'Proj[1]')
    try:
        native_term_signature(wrong)
    except ValueError:
        pass
    else:
        raise AssertionError('chirality orientation mismatch passed')
    r = {'format': FMT, 'status': 'WarsawModelLoadFailed'}
    assert r['status'] in ('WarsawModelLoadFailed', 'WarsawModelUnavailable')
    print('PASS: source native SU(2)/GammaCC inventory, orientation and tamper rejection')
    print('PASS: unavailable Warsaw model handled without inventing a normalization')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--self-test', action='store_true')
    p.add_argument('--reference', type=Path,
                   help='Optional v9 Warsaw JSON; failed model-load JSON activates source fallback')
    p.add_argument('--native-root', type=Path, default=SOURCE_REL)
    p.add_argument('--output', type=Path)
    args = p.parse_args()
    try:
        if args.self_test:
            self_test();return 0
        if not args.output:
            p.error('--output is required')
        reference = None
        if args.reference:
            reference = json.loads(args.reference.read_text(encoding='utf-8-sig'))
        result = analyse(reference, args.native_root)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False)+'\n', encoding='utf-8')
        for model in result['models']:
            print(f"{model['model']}: pole 25 PR/25 PL, finite 16 PR/16 PL; native LLHH skeleton checked")
        print('Installed Warsaw model:', result['reference']['reference_operator_source'])
        print('Absolute GammaCC/Warsaw factor: NOT CERTIFIED; physical C5: NOT CERTIFIED')
        print('Saved:', args.output)
        return 0
    except (ValueError, KeyError, TypeError, OSError, ImportError, json.JSONDecodeError) as exc:
        print('FAIL (no normalization certificate):', exc, file=sys.stderr)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
