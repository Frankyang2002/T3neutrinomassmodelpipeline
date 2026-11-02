"""Regression tests for exact, fail-closed Matchete raw-LLHH classification."""
from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

from Numerical.lo_certification.spinor.boundary import BoundaryError, MODELS, tag
from Numerical.lo_certification.spinor.raw_llhh import (
    split_top_level, classify, analyze_text, sha_bytes,
    validate_batch_directory, create_inventory)

# Close to the real InputForm, with nested SU2 brackets and charged scalar HC.
H = 'Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]'
L = 'Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}]'
PL = f'NCM[Transp[{L}], DiracProduct[GammaCC, Proj[-1]], {L}]'
PR = f'NCM[Bar[{L}], DiracProduct[GammaCC, Proj[1]], Transp[Bar[{L}]]]'


def terms():
    pl = f'(hbar*Bar[{H}]*{H}*Coupling[y1, {{}}, 0]*{PL})/(ϵ*Coupling[MF, {{}}, 0])'
    pr = f'(hbar*Bar[{H}]*{H}*Coupling[y1, {{}}, 0]*{PR})/Coupling[MF, {{}}, 0]'
    return pl, pr


class RawLLHHTests(unittest.TestCase):
    def test_signed_nested_expression_exact_split(self):
        a, b = terms()
        source = f'{a} - ({b}) + (2*(hbar+1))/(ϵ)'
        got = split_top_level(source)
        self.assertEqual(got, [a, f'- ({b})', '+ (2*(hbar+1))/(ϵ)'])
        self.assertEqual(' '.join(got), source)

    def test_reject_unclosed_bracket(self):
        with self.assertRaises(BoundaryError):
            split_top_level('(hbar*Field[H, Scalar, {}, {}]')

    def test_quote_with_plus_is_not_operator(self):
        self.assertEqual(split_top_level('f["a + b"] + c'), ['f["a + b"]', '+ c'])

    def test_one_loop_chiral_higgs_only(self):
        a, b = terms()
        x, y = classify(a), classify(b)
        self.assertTrue(x['is_LLHH'])
        self.assertTrue(x['explicit_epsilon_pole'])
        self.assertEqual(x['chirality'], 'PL')
        self.assertEqual(x['fields'], {'H': 2, 'l': 2})
        self.assertEqual(x['couplings'], ['MF', 'y1'])
        self.assertEqual(y['chirality'], 'PR')
        self.assertFalse(y['explicit_epsilon_pole'])

    def test_extra_scalar_rejected(self):
        a, _ = terms()
        s = a + '*Field[NewScalarY, Scalar, {}, {}]'
        self.assertFalse(classify(s)['is_LLHH'])

    def test_wrong_chirality_not_selected(self):
        a, _ = terms()
        self.assertFalse(classify(a.replace('Proj[-1]', 'Proj[1]')
                                  .replace('GammaCC', 'GammaCC, GammaM[mu]'))['is_LLHH'])

    def test_explicit_power_counted(self):
        a, _ = terms()
        self.assertTrue(classify(a.replace(f'Bar[{H}]*{H}', f'{H}^2'))['is_LLHH'])

    def test_nonhbar_llhh_kept_separate(self):
        a, b = terms()
        x = classify(a.replace('hbar*', ''))
        self.assertTrue(x['is_LLHH'])
        self.assertFalse(x['is_one_loop_LLHH'])
        self.assertEqual(x['chirality'], 'PL')

    def test_unpaired_chiral_terms_block(self):
        with self.assertRaisesRegex(BoundaryError, 'PL/PR'):
            analyze_text('A_-4', terms()[0])

    def test_poles_and_tree_contamination_not_renormalized(self):
        a, b = terms()
        x = analyze_text('A_-4', 'Coupling[MS1, {}, 1]^2 + ' + a + ' + ' + b)
        self.assertEqual(x['PL_holomorphic_terms'], 1)
        self.assertEqual(x['PR_conjugate_terms_excluded_from_PL'], 1)
        self.assertEqual(x['explicit_epsilon_pole_terms'], 1)
        self.assertEqual(x['non_hbar_terms_in_raw_export'], 1)
        self.assertFalse(x['MSbar_subtraction_certified'])
        self.assertIsNone(x['physical_C5'])

    def test_validate_all16_raw_batch_checks_hashes(self):
        a, b = terms()
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            models = []
            for m in MODELS:
                d = root / tag(m)
                d.mkdir()
                is_mixed = m in ('B_-1', 'C_-1')
                f = (d / f'{tag(m)}_native_one_loop_raw_hard_llhh_sector.inputform.txt') if is_mixed \
                    else (d/'data/eft1_after_F_one_loop_inputform.txt')
                f.parent.mkdir(parents=True, exist_ok=True)
                f.write_text(a+' + '+b, encoding='utf8')
                report = d / ('native_one_loop.json' if is_mixed else 'native_status.json')
                report.write_text('{"dummy":1}', encoding='utf8')
                (d/'wolfram.log').write_text('mock native log', encoding='utf8')
                rec = dict(model=m, native_exit_code=0,
                           hard_sector_renormalized=False, physical_wilson_coefficient=None,
                           native_report_file=report.name,
                           native_report_sha256=sha_bytes(report.read_bytes()),
                           native_log_sha256=sha_bytes((d/'wolfram.log').read_bytes()))
                if is_mixed:
                    rec.update(status='MATCHeTE_RAW_MIXED_LOOP_LLHH_EXTRACTED_NOT_MSbar',
                               raw_hard_llhh_file=f.name,
                               raw_hard_llhh_sha256=sha_bytes(f.read_bytes()),
                               explicit_hard_LLHH_candidate_terms=2)
                else:
                    rec.update(status='MATCHeTE_F_ONLY_LOOP_ORDER_ONE_RETURNED_RAW_SECTOR_UNRESOLVED',
                               raw_loop_inputform_file='data/'+f.name,
                               raw_loop_inputform_sha256=sha_bytes(f.read_bytes()))
                models.append(rec)
                (d/'inspection.json').write_text(json.dumps(rec), encoding='utf8')
            summary = dict(format='t3-f-only-raw-one-loop-native-evidence-v1',
                           status='RAW_F_FIRST_ONE_LOOP_EVIDENCE_COLLECTED_NOT_PHYSICAL',
                           n_models_with_raw_native_one_loop_evidence=16,
                           models_requested=list(MODELS), model_records=models,
                           single_simultaneous_fermion_threshold=True,
                           both_scalar_multiplets_retained_at_F_threshold=True,
                           SM_fermions_remain_active=True,
                           n_models_with_MSbar_hard_sector_certified=0,
                           native_GammaCC_to_master_Weyl_phase=None,
                           physical_master_wilson_boundary=None,
                           mu_F_gev=1e5)
            (root/'summary.json').write_text(json.dumps(summary), encoding='utf8')
            _, inputs = validate_batch_directory(root)
            self.assertEqual(len(inputs), 16)
            self.assertTrue(all(create_inventory(m, inputs[m])['selected_llhh_term_count'] == 2
                                for m in MODELS))
            # A source-input corruption must abort before generating any output.
            (root/'A_m4/data/eft1_after_F_one_loop_inputform.txt').write_text('0')
            with self.assertRaisesRegex(BoundaryError, 'hash mismatch'):
                validate_batch_directory(root)


if __name__ == '__main__':
    unittest.main()
