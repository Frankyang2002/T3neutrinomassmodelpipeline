"""Offline synthetic fixtures only; does not claim native RGBeta execution."""
from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import tempfile
import unittest

from scripts import InspectAminus4RGBetaSymbolicContractions as audit


class TestAminus4SymbolicContractions(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.source = self.root / audit.SOURCE_REL
        self.source.parent.mkdir(parents=True)
        self.source.write_text('(* synthetic test source, NOT production *)\n',
                               encoding='utf8')
        self.digest = sha256(self.source.read_bytes()).hexdigest()
        pkg_root = self.root / 'InstalledRGBeta'
        (pkg_root / 'Kernel').mkdir(parents=True)
        self.pkg = pkg_root / 'Kernel' / 'init.m'
        self.pkg.write_text('(* synthetic RGBeta init *)\n', encoding='utf8')
        (pkg_root / 'GroupInvariants.m').write_text(
            'tGen[r_,A_,i_,j_] := nativeRule[r,A,i,j]\n'
            'CheckPotentialCouplingMixing[x_] := x\n'
            'BetaTerm[Quartic,1] := projected\n', encoding='utf8')
        self.pkg_input = json.dumps(str(self.pkg))
        self.previous = {
            'format': audit.OLD_FORMAT,
            'status': 'NATIVE_SYMBOLIC_OPERATOR_TENSOR_CAPTURED_COMPONENT_MAP_PENDING',
            'model': 'A_-4', 'stage': 'EFT1',
            'dimensions_alpha': [1, 3, 2, -4],
            'source_model_path': audit.SOURCE_REL.as_posix(),
            'source_model_sha256': self.digest,
            'rgbeta_package_file': self.pkg_input,
            'native_operator_registry_complete': True,
            'native_Lam_formal_tensor_available': True,
            'matchete_physical_coupling_map_certified': False,
            'native_to_reference_component_axes_certified': False,
        }
        q = {
            'lambdaS2': 'del[SU2L[S2],aci,acj]*del[SU2L[S2],ack,acl]/2',
            'lambdaS2Adj': 'tGen[SU2L[S2],acA,aci,acj]*tGen[SU2L[S2],acA,ack,acl]/2',
            'lambdaH2': 'del[SU2L[fund],aci,acj]*del[SU2L[S2],ack,acl]',
            'lambdaH2Adj': 'tGen[SU2L[fund],acA,aci,acj]*tGen[SU2L[S2],acA,ack,acl]',
            'lambdaT3': 'delS2[SU2L,acl,aci,acj]',
        }
        group = [
            {'id': name, 'status': 'SYMBOLICALLY_UNRESOLVED_OR_DIFFERENT',
             'source_expression_inputform': 'HoldComplete[tGen[SU2L[S2],a,i,j]*tGen[SU2L[S2],a,j,k]]',
             'native_evaluated_inputform': 'tGen[SU2L[S2],a,i,j]*tGen[SU2L[S2],a,j,k]',
             'target_inputform': '2*del[SU2L[S2],i,k]',
             'native_difference_inputform': 'unreduced',
             'native_exact_target_match': False}
            for name in audit.GROUP_IDS
        ]
        self.native = {
            'format': audit.FORMAT, 'status': audit.STATUS,
            'model': 'A_-4', 'stage': 'EFT1',
            'dimensions_alpha': [1, 3, 2, -4],
            'source_model_path': audit.SOURCE_REL.as_posix(),
            'source_model_sha256': self.digest,
            'rgbeta_package_file': self.pkg_input,
            'generator_contraction_probes': group,
            'quartic_invariant_samples': [
                {'name': name, 'captured': True,
                 'native_invariant_inputform': expr,
                 'field_order_inputform': 'recorded in previous witness'}
                for name, expr in q.items()],
            'native_dummy_index_double_call_probes': [
                {'id': name,
                 'first_call_inputform': 'tGen[SU2L[S2],A,i,j]',
                 'second_call_inputform': 'tGen[SU2L[S2],A,k,l]',
                 'same_literal_A_present_in_both_raw_calls': True}
                for name in ('S2Adj_times_S2Adj', 'H2Adj_times_H2Adj',
                             'S2Adj_times_H2Adj')
            ] + [{
                'id': 'S2Singlet_times_S2Singlet_control',
                'first_call_inputform': 'del[SU2L[S2],i,j]',
                'second_call_inputform': 'del[SU2L[S2],k,l]',
                'same_literal_A_present_in_both_raw_calls': False,
            }],
            'native_projector_related_symbol_names': [
                'RGBeta`PackageScope`ProjectQuartics'],
            'native_projector_symbol_definitions': [
                {'symbol': 'RGBeta`PackageScope`ProjectQuartics',
                 'definition_prefix_inputform': 'DownValues[...]',
                 'full_definition_length': 15, 'truncated': False}],
            'native_group_reductions_certified': False,
            'quartic_projection_certified': False,
            'physical_matchete_coupling_map_certified': False,
            'native_to_reference_axes_certified': False,
            'source_modified': False,
        }

    def test_valid_native_report_uncertified(self):
        result = audit.analyze(self.native, self.previous, self.root,
                               source_scan=False)
        self.assertEqual(result['generator_probes_evaluated_to_reference_targets'], 0)
        self.assertFalse(result['scalar_beta_polynomial_certified'])
        self.assertFalse(result['native_to_Matchete_component_map_certified'])
        self.assertEqual(result['native_operator_registry_verified'], 9)
        self.assertEqual(result['native_raw_dummy_index_reuse_count'], 3)
        self.assertFalse(result['native_projector_alpha_renaming_verified'])

    def test_one_exact_native_group_match(self):
        x = self.native['generator_contraction_probes'][1]
        x['status'] = 'EVALUATED_MATCH'
        x['native_exact_target_match'] = True
        x['native_difference_inputform'] = '0'
        result = audit.analyze(self.native, self.previous, self.root,
                               source_scan=False)
        self.assertEqual(result['generator_probes_evaluated_to_reference_targets'], 1)
        self.assertFalse(result['physical_quartic_normalization_certified'])

    def test_stale_source_blocks(self):
        self.source.write_text('(* changed *)', encoding='utf8')
        with self.assertRaisesRegex(audit.ContractionAuditError, 'stale'):
            audit.analyze(self.native, self.previous, self.root,
                          source_scan=False)

    def test_stale_previous_witness_blocks(self):
        self.previous['source_model_sha256'] = '0'*64
        with self.assertRaises(audit.ContractionAuditError):
            audit.analyze(self.native, self.previous, self.root,
                          source_scan=False)

    def test_different_installed_rgbeta_blocks(self):
        self.previous['rgbeta_package_file'] = '"OtherRGBeta"'
        with self.assertRaises(audit.ContractionAuditError):
            audit.analyze(self.native, self.previous, self.root,
                          source_scan=False)

    def test_no_fabricated_generator_match(self):
        x = self.native['generator_contraction_probes'][0]
        x['native_exact_target_match'] = True
        with self.assertRaisesRegex(audit.ContractionAuditError, 'inconsistent'):
            audit.analyze(self.native, self.previous, self.root,
                          source_scan=False)

    def test_no_fabricated_status(self):
        self.native['generator_contraction_probes'][0]['status'] = 'CERTIFIED_PHYSICAL'
        with self.assertRaises(audit.ContractionAuditError):
            audit.analyze(self.native, self.previous, self.root,
                          source_scan=False)

    def test_incomplete_native_registration_blocks(self):
        self.previous['native_operator_registry_complete'] = False
        with self.assertRaises(audit.ContractionAuditError):
            audit.analyze(self.native, self.previous, self.root,
                          source_scan=False)

    def test_unsafe_certification_flag_blocks(self):
        self.native['quartic_projection_certified'] = True
        with self.assertRaises(audit.ContractionAuditError):
            audit.analyze(self.native, self.previous, self.root,
                          source_scan=False)

    def test_invariant_wrong_representation_blocks(self):
        self.native['quartic_invariant_samples'][1]['native_invariant_inputform'] = \
            'tGen[SU2L[fund],a,i,j]'
        with self.assertRaises(audit.ContractionAuditError):
            audit.analyze(self.native, self.previous, self.root,
                          source_scan=False)

    def test_duplicate_symbol_blocks(self):
        self.native['native_projector_related_symbol_names'] *= 2
        with self.assertRaises(audit.ContractionAuditError):
            audit.analyze(self.native, self.previous, self.root,
                          source_scan=False)

    def test_native_source_scan(self):
        result = audit.analyze(self.native, self.previous, self.root,
                               source_scan=True)
        scan = result['installed_package_source']
        self.assertEqual(scan['status'], 'INSTALLED_PACKAGE_SOURCE_SCANNED')
        self.assertEqual(len(scan['source_files']), 2)
        self.assertTrue(any('CheckPotentialCouplingMixing' in x['code_excerpt']
                            for x in scan['relevant_source_snippets']))

    def test_missing_package_source_not_physical_certificate(self):
        self.pkg.unlink()
        result = audit.analyze(self.native, self.previous, self.root,
                               source_scan=True)
        self.assertEqual(result['installed_package_source']['status'],
                         'INSTALLED_PACKAGE_NOT_ACCESSIBLE')
        self.assertFalse(result['scalar_beta_polynomial_certified'])

    def test_safe_quoted_path(self):
        got = audit.parse_installed_path(json.dumps(str(self.pkg)))
        self.assertEqual(got, self.pkg)

    def test_injected_path_not_evaluated(self):
        with self.assertRaises(audit.ContractionAuditError):
            audit.parse_installed_path('"RGBeta"+__import__("os").system("echo PWNED")')

    def test_dummy_index_probe_missing_blocks(self):
        self.native['native_dummy_index_double_call_probes'].pop()
        with self.assertRaises(audit.ContractionAuditError):
            audit.analyze(self.native, self.previous, self.root,
                          source_scan=False)

    def test_double_call_without_generator_blocks(self):
        row = self.native['native_dummy_index_double_call_probes'][0]
        row['first_call_inputform'] = 'del[SU2L[S2],i,j]'
        with self.assertRaises(audit.ContractionAuditError):
            audit.analyze(self.native, self.previous, self.root,
                          source_scan=False)

    def test_missing_quartic_sample_blocks(self):
        self.native['quartic_invariant_samples'].pop()
        with self.assertRaises(audit.ContractionAuditError):
            audit.analyze(self.native, self.previous, self.root,
                          source_scan=False)

    def test_no_replacement_beta_file_emitted(self):
        self.assertNotIn('compatible_beta',
                         audit.analyze(self.native, self.previous, self.root,
                                       source_scan=False))


if __name__ == '__main__':
    unittest.main()
