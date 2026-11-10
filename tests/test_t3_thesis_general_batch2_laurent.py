"""Offline unit tests: no Wolfram required; does NOT certify loop integrals."""
from __future__ import annotations

import argparse
import tempfile
import unittest
from pathlib import Path

import numpy as np

from Numerical.parallel_t3.background.general_batch1_api import encode_native
from Numerical.parallel_t3.wilson.general_batch2_laurent import (
    SourceAuditError, REPORT_FORMAT, check_audit_report,
    matrix_spurion_diagnostics, sha,
)
from scripts.RunT3ThesisGeneralBatch2Laurent import safe_id


def input_record(model: str = 'A_-4', *, offdiag: bool = False) -> dict:
    MF = np.diag([1e5, 1.1e5, 1.2e5]).astype(complex)
    if offdiag:
        MF[0, 1] = 1e3
    a = np.array([[1., 2., 3.], [4., 5., 6.], [7., 8., 9.]])*1e-4
    b = np.array([[9., 3., 1.], [2., 8., 4.], [5., 6., 7.]])*1e-4
    values = {'MF': encode_native(MF), 'y1': encode_native(a),
              'y2': encode_native(b)}
    if model == 'D_-2':
        values['y2Conj'] = encode_native(np.zeros((3, 3)))
        values['muT3GeneralS2Sq'] = encode_native(complex(-4338.314613978704,0))
        for k in ('000004', '000013', '001102', '020101', '110002'):
            values['lambdaT3GeneralQ' + k + 'I1'] = encode_native(.01)
    return {'model': model, 'native_parameter_values': values}


class SpuriousWilsonProtectionTests(unittest.TestCase):
    def test_diagonal_mf_matches_source_diagonal_formula(self):
        x = matrix_spurion_diagnostics(input_record(), 'A_-4')
        self.assertFalse(x['physical_Wilson_coefficient'])
        self.assertLess(x['baseline_Y1_Y2']['orientation_max_abs_difference_GeVinv'], 1e-27)
        y = np.asarray(input_record()['native_parameter_values']['y1'])
        self.assertEqual(len(y), 3)
        self.assertEqual(x['MF_off_diagonal_relative_Frobenius'], 0.)

    def test_nondiagonal_mf_orientation_not_silently_assumed(self):
        x = matrix_spurion_diagnostics(input_record(offdiag=True), 'A_-4')
        self.assertGreater(x['baseline_Y1_Y2']['orientation_max_abs_difference_GeVinv'], 0.)
        self.assertGreater(x['MF_off_diagonal_relative_Frobenius'], 0.)
        self.assertFalse(x['native_to_Matchete_heavy_chirality_map_verified'])

    def test_d2_general_inventory_preserved(self):
        x = matrix_spurion_diagnostics(input_record('D_-2'), 'D_-2')
        self.assertEqual(len(x['D_five_extra_quartics_at_F']), 5)
        self.assertAlmostEqual(x['D_self_bilinear_at_F_GeV2']['real'], -4338.314613978704)
        self.assertIn('extra_Y1_Y2Conj', x)

    def test_d2_extra_coupling_missing_rejected(self):
        d = input_record('D_-2')
        d['native_parameter_values'].pop('muT3GeneralS2Sq')
        with self.assertRaises(SourceAuditError):
            matrix_spurion_diagnostics(d, 'D_-2')

    def test_noninvertible_mf_rejected(self):
        d = input_record()
        d['native_parameter_values']['MF'] = encode_native(np.zeros((3,3)))
        with self.assertRaises(SourceAuditError):
            matrix_spurion_diagnostics(d, 'A_-4')

    def _report(self, path: Path, *, tree_terms: int = 4):
        directory = path.parent
        files = {}
        for name in ('raw_hbar1_coefficient.wxf',
                     'unreconciled_hbar0_mass_remainder.wxf',
                     'evaluated_full_hbar1_coefficient.wxf',
                     'explicit_pole_residue_unsubtracted.wxf',
                     'epsilon_zero_regular_unsubtracted.wxf'):
            dst = directory / name
            dst.write_bytes(('mock WXF: ' + name).encode('utf-8'))
            files[name] = sha(dst)
        hbar1 = 131 if tree_terms == 4 else 1194
        return {
            'format': REPORT_FORMAT,
            'status': 'LoopFunctionsEvaluated_ExplicitLaurentSeparated_TreeRemainderPreserved_AuditOnly',
            'source_sha256': 'native-source-hash',
            'one_loop_hbar_linearity_exact': tree_terms == 0,
            'hbar_higher_order_absent': True,
            'hbar1_nonzero': True,
            'full_expression_hbar_decomposition_exact': True,
            'unreconciled_hbar0_mass_remainder_present': tree_terms > 0,
            'unreconciled_hbar0_mass_term_count': tree_terms,
            'hbar0_only_recognised_quadratic_mass_terms': True,
            'raw_hbar1_term_count': hbar1,
            'original_source_term_count': hbar1 + tree_terms,
            'unreconciled_tree_mass_coordinate_identity_verified': False,
            'tree_mass_remainder_dropped_from_source': False,
            'tree_mass_remainder_included_in_one_loop_coefficient': False,
            'explicit_laurent_reconstruction_exact': True,
            'all_operators_preserved': True,
            'native_fermions_absent': True,
            'renormalised_MSbar_Wilson_boundary_computed': False,
            'source_scale_substitution_performed': False,
            'output_sha256': files,
        }

    def test_successful_audit_is_still_not_physical_matching(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)/'raw_hbar1_coefficient.wxf'
            report = self._report(p)
            rec = check_audit_report(report, 'native-source-hash', Path(tmp))
            self.assertFalse(rec['numerical_finite_Wilson_matching_complete'])
            self.assertTrue(rec['pole_expression_formal_only'])
            self.assertFalse(rec['original_expression_hbar_linear'])
            self.assertEqual(rec['unreconciled_hbar0_mass_term_count'], 4)
            self.assertEqual(rec['raw_hbar1_term_count'], 131)

    def test_audit_source_checksum_mismatch_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)/'raw_hbar1_coefficient.wxf'
            report = self._report(p)
            with self.assertRaises(SourceAuditError):
                check_audit_report(report, 'changed-source', Path(tmp))

    def test_audit_output_checksum_mismatch_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)/'raw_hbar1_coefficient.wxf'
            report = self._report(p)
            p.write_bytes(b'corrupted wxf')
            with self.assertRaises(SourceAuditError):
                check_audit_report(report, 'native-source-hash', Path(tmp))

    def test_blocked_laurent_cannot_be_called_success(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)/'raw_hbar1_coefficient.wxf'
            report = self._report(p)
            report['status'] = 'Blocked_UnresolvedLoopFunctions'
            with self.assertRaises(SourceAuditError):
                check_audit_report(report, 'native-source-hash', Path(tmp))

    def test_d2_1200_terms_preserved_separately(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)/'raw_hbar1_coefficient.wxf'
            report = self._report(p, tree_terms=6)
            rec = check_audit_report(report, 'native-source-hash', Path(tmp))
            self.assertEqual(rec['unreconciled_hbar0_mass_term_count'], 6)
            self.assertEqual(rec['raw_hbar1_term_count'], 1194)
            self.assertFalse(rec['numerical_finite_Wilson_matching_complete'])

    def test_unreconciled_tree_wxf_cannot_be_omitted(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)/'raw_hbar1_coefficient.wxf'
            report = self._report(p)
            report['output_sha256'].pop('unreconciled_hbar0_mass_remainder.wxf')
            with self.assertRaises(SourceAuditError):
                check_audit_report(report, 'native-source-hash', Path(tmp))

    def test_tree_mass_identity_not_assumed(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)/'raw_hbar1_coefficient.wxf'
            report = self._report(p)
            report['unreconciled_tree_mass_coordinate_identity_verified'] = True
            with self.assertRaises(SourceAuditError):
                check_audit_report(report, 'native-source-hash', Path(tmp))

    def test_tree_mass_inventory_required(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)/'raw_hbar1_coefficient.wxf'
            report = self._report(p)
            report['hbar0_only_recognised_quadratic_mass_terms'] = False
            with self.assertRaises(SourceAuditError):
                check_audit_report(report, 'native-source-hash', Path(tmp))

    def test_hbar2_contamination_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)/'raw_hbar1_coefficient.wxf'
            report = self._report(p)
            report['hbar_higher_order_absent'] = False
            with self.assertRaises(SourceAuditError):
                check_audit_report(report, 'native-source-hash', Path(tmp))

    def test_source_terms_must_be_accounted(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)/'raw_hbar1_coefficient.wxf'
            report = self._report(p)
            report['original_source_term_count'] -= 1
            with self.assertRaises(SourceAuditError):
                check_audit_report(report, 'native-source-hash', Path(tmp))

    def test_laurent_pole_wxf_must_be_present(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)/'raw_hbar1_coefficient.wxf'
            report = self._report(p)
            report['output_sha256'].pop('explicit_pole_residue_unsubtracted.wxf')
            with self.assertRaises(SourceAuditError):
                check_audit_report(report, 'native-source-hash', Path(tmp))

    def test_run_id_is_validation_only(self):
        self.assertEqual(safe_id('thesis_batch2_laurent_v1'), 'thesis_batch2_laurent_v1')
        for bad in ('../bad','bad name','', 'Z'*65):
            with self.assertRaises(argparse.ArgumentTypeError):
                safe_id(bad)


if __name__ == '__main__':
    unittest.main()
