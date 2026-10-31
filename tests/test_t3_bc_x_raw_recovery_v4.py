"""Read-only v3 X-only matching report-recovery regression tests."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from scripts.RecoverT3BCXRawMatching import (
    ARTIFACTS, FORMAT, MODELS, PROGRESS_STAGES, RELATIVE_V3,
    STATUS, V3_WL_SHA256, count_top_level_terms, inspect_model, recover,
)


class SyntaxAndSource(unittest.TestCase):
    def test_top_level_terms(self):
        self.assertEqual(count_top_level_terms('-(a+b) + (a*b) - (c[d]+e)'), 3)
        self.assertEqual(count_top_level_terms('x*y'), 1)
        self.assertEqual(count_top_level_terms('(a + b) - c'), 2)

    def test_invalid_unbalanced(self):
        for s in ('', 'a+(b', 'x+[a}', '"unterminated'):
            with self.subTest(expr=s), self.assertRaises(RuntimeError):
                count_top_level_terms(s)

    def test_full_wolfram_report_fix(self):
        source = (Path(__file__).resolve().parents[1] /
                  'Lagrangian/RunT3BCXFromSavedY.experimental.wl').read_text(encoding='utf8')
        self.assertIn('report=Join[report,<|"status"->"XFollowupYTreeOnlyXRawMatchingNotCertified"', source)
        self.assertIn('report=Join[report,<|"status"->"XFollowupTreeMatchingFailed"', source)
        self.assertIn('report=Join[report,<|"status"->"XFollowupOneLoopMatchingFailed"', source)
        self.assertNotIn('AssociateTo[report,"status"->"XFollowupYTreeOnlyXRawMatchingNotCertified"', source)
        self.assertIn('saved_fermion_EOMEFT_equivalence_certified"->False', source)
        self.assertIn('"y_loop_inherited_in_x_matching"->False', source)
        self.assertIn('XFollowupFinalReportExportFailed', source)


def fixture(root: Path, model: str):
    tag = model.replace('+', 'p').replace('-', 'm')
    base = root / RELATIVE_V3
    base.mkdir(parents=True, exist_ok=True)
    stem = tag + '_final_x_from_saved_y_match'
    prefix = base / stem
    checks = {s: True for s in ('ordinary_check', 'hermiticity', 'canonical',
                                'heavy_mass_basis', 'gauge_invariant', 'all_objects_defined')}
    x = {'model': model, 'status': 'XFollowupOneLoopMatchingInProgress',
         'x_tree_matching_completed': True, 'physical_matching_certified': False,
         'matching_input_selection': 'heavy_free_rebased_y_tree_eomeft',
         'rebased_check_passed': True, 'selected_x_input_full_matchete_check_passed': True,
         'source_is_y_tree_EOMEFT': True, 'y_loop_inherited_in_x_matching': False,
         'running_threshold_couplings_applied': False,
         'saved_fermion_EOMEFT_equivalence_certified': False,
         'rebased_x_heavy_eft2_detailed_check': checks,
         'light_x_eft2_detailed_check': checks}
    y = {'model': model, 'status': 'SavedYTreeEOMConsumedForFreshXNoYRematch',
         'x_followup_status': x['status'], 'physical_sequential_matching_certified': False}
    provenance = {'model': model, 'format': 't3-bc-fresh-x-only-native-attempt-v3',
                  'mode': 'match', 'exit_code': 0, 'fresh_wl_sha256': V3_WL_SHA256,
                  'fresh_x_matchete_matching_completed': False,
                  'physical_matching_certified': False}
    progress = {'model': model, 'stage_history': list(PROGRESS_STAGES),
                'last_stage': PROGRESS_STAGES[-1], 'physical_matching_performed': False}
    cross = {'status': 'XHeavyRegistryVersusRebaseComparisonCompleted',
             'rebased_heavy_check_passed': True, 'matched_eft_fallback_used': False}
    for suffix, obj in (('_x_followup.json', x), ('.json', y),
                        ('.provenance.json', provenance),
                        ('.progress.json', progress),
                        ('_x_followup_rebase_crosschecks.json', cross)):
        (base / (stem + suffix)).write_text(json.dumps(obj), encoding='utf8')
    (base / (stem + '.log')).write_text('exit=0\n' +
        'Matching pipeline completed successfully.\n' * 2 +
        'AssociateTo::argrx: 6 arguments\n', encoding='utf8')
    expression = (' + '.join('hbar*GammaCC*Field[H, Scalar, {a}]*'
                             'Field[l, Fermion, {i}]' for _ in range(36)))
    for suffix in ARTIFACTS:
        if suffix == '_x_raw_llhh.inputform.txt':
            s = expression
        elif suffix == '_y_tree_eomeft.inputform.txt':
            s = 'Field[NewScalarX, Scalar, {i}]*' + ('a' * 650)
        elif suffix == '_y_full_eomeft.inputform.txt':
            s = 'Field[NewScalarX, Scalar, {i}]*' + ('b' * 650)
        else:
            s = 'hbar*Field[H, Scalar, {i}]*' + ('a' * 650)
        (base / (stem + '_x_followup' + suffix)).write_text(s, encoding='utf8')
    return prefix


class Recovery(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.prefix = fixture(self.root, 'B_-1')

    def test_success_and_no_old_file_mutation(self):
        original = (self.prefix.parent / (self.prefix.name + '_x_followup.json'))
        data = original.read_bytes()
        report, files = inspect_model('B_-1', self.root)
        self.assertEqual(report['status'], STATUS)
        self.assertEqual(report['x_raw_llhh_term_count'], 36)
        self.assertFalse(report['physical_sequential_matching_certified'])
        summary = recover(('B_-1',), root=self.root)
        self.assertEqual(summary['format'], FORMAT)
        self.assertEqual(original.read_bytes(), data)
        out = self.root / 'output/audits/t3_bc_x_raw_matching_recovery_v4'
        self.assertTrue((out / 'recovered_x_native_raw_evidence_v4.zip').is_file())

    def test_reject_wrong_wolfram_source_hash(self):
        p = self.prefix.with_suffix('.provenance.json')
        data = json.loads(p.read_text())
        data['fresh_wl_sha256'] = '0' * 64
        p.write_text(json.dumps(data))
        with self.assertRaises(RuntimeError):
            inspect_model('B_-1', self.root)

    def test_reject_uncertified_stage(self):
        p = self.prefix.with_suffix('.progress.json')
        data = json.loads(p.read_text())
        data['stage_history'].remove('x-followup-one-loop-matching-returned')
        p.write_text(json.dumps(data))
        with self.assertRaises(RuntimeError):
            inspect_model('B_-1', self.root)

    def test_reject_partial_llhh(self):
        p = self.prefix.parent / (self.prefix.name + '_x_followup_x_raw_llhh.inputform.txt')
        p.write_text('hbar*GammaCC')
        with self.assertRaises(RuntimeError):
            inspect_model('B_-1', self.root)

    def test_reject_unproven_fallback(self):
        p = self.prefix.parent / (self.prefix.name + '_x_followup.json')
        data = json.loads(p.read_text())
        data['matching_input_selection'] = 'MatchedEFT'
        p.write_text(json.dumps(data))
        with self.assertRaises(RuntimeError):
            inspect_model('B_-1', self.root)

    def test_no_partial_writes_when_other_model_bad(self):
        fixture(self.root, 'C_-1')
        p = self.root / RELATIVE_V3 / 'C_m1_final_x_from_saved_y_match.provenance.json'
        d = json.loads(p.read_text()); d['exit_code'] = 3; p.write_text(json.dumps(d))
        with self.assertRaises(RuntimeError):
            recover(MODELS, root=self.root)
        self.assertFalse((self.root/'output/audits/t3_bc_x_raw_matching_recovery_v4').exists())


if __name__ == '__main__':
    unittest.main()
