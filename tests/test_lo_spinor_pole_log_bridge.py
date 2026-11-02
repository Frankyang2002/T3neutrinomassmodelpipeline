"""Exact raw-native Laurent identities, no physical C5 promotion."""
from __future__ import annotations

import copy
import tempfile
from pathlib import Path
import subprocess
import sys
import unittest

from Numerical.lo_certification.spinor.boundary import BoundaryError, MODELS
from Numerical.lo_certification.spinor.raw_llhh import analyze_text, classify, sha_bytes, split_top_level
from Numerical.lo_certification.spinor.pole_log_bridge import (
    normalize_native_kernel, validate_and_group, expressions_for,
    EXPECTED_SNAPSHOT, FORMAT,
)

H = 'Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]'
L = 'Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}]'
MF = 'Coupling[MF, {Index[d$$1, NFlavor]}, 0]'
PL = f'NCM[Transp[{L}], DiracProduct[GammaCC, Proj[-1]], {L}]'
PR = f'NCM[Bar[{L}], DiracProduct[GammaCC, Proj[1]], Transp[Bar[{L}]]]'


def triple(chirality='PL', coefficient=2, first_implicit=False):
    spinor = PL if chirality == 'PL' else PR
    numerator = f'{coefficient}*hbar*Coupling[y1, {{}}, 0]*{H}^2*{spinor}'
    base = f'({numerator})/{MF}'
    pole = f'({numerator})/(ϵ*{MF})'
    logged = f'({numerator.replace("*" + spinor, "*Log[μbar2/" + MF + "^2]*" + spinor)})/{MF}'
    if first_implicit:
        return base, '+ ' + pole, '+ ' + logged
    return '+ ' + base, '+ ' + pole, '+ ' + logged


def make_inventory(model='A_-4'):
    # Real Matchete mixes constant/pole first, followed by all logarithms.
    a = triple('PR', first_implicit=True)
    b = triple('PL')
    text = ' '.join([a[0], a[1], b[0], b[1], a[2], b[2]])
    result = analyze_text(model, text)
    result['native_raw_inputform_sha256'] = sha_bytes(text.encode('utf8'))
    result['raw_source_file'] = 'test.inputform.txt'
    result['native_raw_inputform_bytes'] = len(text.encode('utf8'))
    return result


class PoleLogBridgeTests(unittest.TestCase):
    def test_normalize_pole_and_finite_and_log_exact(self):
        base, pole, logged = triple('PL', first_implicit=True)
        n1, k1, arg1 = normalize_native_kernel(base)
        n2, k2, arg2 = normalize_native_kernel(pole)
        n3, k3, arg3 = normalize_native_kernel(logged)
        self.assertEqual(n1, n2)
        self.assertEqual(n1, n3)
        self.assertEqual((k1,k2,k3), ('constant','pole','scale_log'))
        self.assertEqual(arg1, '')
        self.assertEqual(arg2, '')
        self.assertEqual(arg3, 'Log[μbar2/' + MF + '^2]')

    def test_nontrivial_real_cg_prefactor_and_mass_denominator(self):
        spinor = PL
        core = f'2*Sqrt[2/3]*hbar*Coupling[y1, {{}}, 0]*{H}^2*{spinor}'
        constant = f'- ({core})/(Sqrt[3]*{MF})'
        pole = f'- ({core})/(Sqrt[3]*ϵ*{MF})'
        log = f'- ({core.replace("*" + spinor, "*Log[μbar2/" + MF + "^2]*" + spinor)})/(Sqrt[3]*{MF})'
        self.assertEqual(len({normalize_native_kernel(q)[0] for q in (constant,pole,log)}), 1)

    def test_source_provenance_and_exact_pl_pr_groups(self):
        got = validate_and_group('A_-4', make_inventory())
        self.assertEqual(got['status'], 'EXACT_NATIVE_POLE_CONSTANT_LOG_TRIPLETS_VERIFIED_NOT_RENORMALIZED')
        self.assertEqual(got['exact_kernel_groups'], 2)
        self.assertEqual((got['PL_groups'], got['PR_groups']), (1,1))
        self.assertEqual((got['explicit_1_over_epsilon_terms'],
                          got['finite_constant_terms'],got['finite_scale_log_terms']), (2,2,2))
        self.assertTrue(got['raw_pole_residue_equals_finite_constant_equals_log_coefficient'])
        self.assertIsNone(got['physical_MSbar_hard_coefficient'])
        self.assertIsNone(got['master_Weyl_phase'])
        self.assertFalse(got['EFT1_soft_subtraction_performed'])
        self.assertFalse(got['counterterm_subtraction_performed'])

    def test_every_model_synthetic_checked_with_no_phase(self):
        for m in MODELS:
            with self.subTest(model=m):
                o = validate_and_group(m, make_inventory(m))
                self.assertEqual(o['exact_kernel_groups'], 2)
                self.assertEqual(o['format'], FORMAT)
                self.assertFalse(o['counterterm_subtraction_performed'])

    def test_finite_constant_equals_residue(self):
        o = validate_and_group('A_-4', make_inventory())
        ex = expressions_for(o, 'PL')
        self.assertEqual(ex['pole_residue_raw'], ex['finite_constant_raw'])
        self.assertIn('Log[μbar2/' + MF + '^2]', ex['scale_log_raw'])
        self.assertIn('Log[', ex['formal_pole_drop_NOT_MSbar'])
        self.assertNotIn('ϵ', ex['formal_pole_drop_NOT_MSbar'])

    def test_exported_formal_sum_is_parseable_and_exact_chirality(self):
        result = validate_and_group('A_-4', make_inventory())
        for ch in ('PL', 'PR'):
            out = expressions_for(result, ch)
            self.assertEqual(len(split_top_level(out['pole_residue_raw'])), 1)
            self.assertEqual(len(split_top_level(out['scale_log_raw'])), 1)
            self.assertEqual(len(split_top_level(out['formal_pole_drop_NOT_MSbar'])), 2)
            for item in split_top_level(out['formal_pole_drop_NOT_MSbar']):
                self.assertEqual(classify(item)['chirality'], ch)

    def test_mixed_mass_basis_scalar_single_MF_factor(self):
        a = triple('PL')
        b = triple('PR')
        self.assertEqual(normalize_native_kernel(a[0])[0], normalize_native_kernel(a[1])[0])
        self.assertEqual(normalize_native_kernel(b[0])[0], normalize_native_kernel(b[2])[0])

    def test_differing_sign_is_not_concealed(self):
        inv = make_inventory()
        for term in inv['selected_terms']:
            if term['chirality'] == 'PL' and term['explicit_epsilon_pole']:
                term['inputform'] = '- ' + term['inputform'][2:]
                term['term_sha256'] = sha_bytes(term['inputform'].encode('utf8'))
        with self.assertRaisesRegex(BoundaryError, 'unmatched'):
            validate_and_group('A_-4', inv)

    def test_differing_flavour_tensor_is_not_concealed(self):
        inv = make_inventory()
        for term in inv['selected_terms']:
            if term['chirality'] == 'PL' and term['explicit_epsilon_pole']:
                term['inputform'] = term['inputform'].replace('y1', 'y2')
                term['term_sha256'] = sha_bytes(term['inputform'].encode('utf8'))
                term['couplings'] = ['MF','y2']
        with self.assertRaisesRegex(BoundaryError, 'unmatched'):
            validate_and_group('A_-4', inv)

    def test_no_double_count_of_hc(self):
        inv = make_inventory()
        o = validate_and_group('A_-4', inv)
        for chirality in ('PL', 'PR'):
            ex = expressions_for(o, chirality)
            self.assertIn('Proj[-1]' if chirality=='PL' else 'Proj[1]', ex['finite_constant_raw'])
            self.assertNotIn('Proj[1]' if chirality=='PL' else 'Proj[-1]', ex['finite_constant_raw'])

    def test_unmatched_pole_rejected(self):
        inv = make_inventory()
        inv['selected_terms'] = [x for x in inv['selected_terms'] if not
             (x['chirality'] == 'PL' and x['contains_log'])]
        inv['selected_llhh_term_count'] -= 1
        inv['nonpole_terms_not_yet_MSbar_subtracted'] -= 1
        inv['terms_with_Log'] -= 1
        inv['PL_holomorphic_terms'] -= 1
        with self.assertRaisesRegex(BoundaryError, 'unmatched'):
            validate_and_group('A_-4', inv)

    def test_multiple_regulator_atoms_rejected(self):
        x = triple()[1].replace('ϵ*', 'ϵ*ϵ*')
        with self.assertRaisesRegex(BoundaryError, 'Multiple epsilon'):
            normalize_native_kernel(x)

    def test_mixed_epsilon_and_log_rejected(self):
        x = triple()[1].replace('hbar*', 'hbar*Log[μbar2/'+MF+'^2]*')
        inv = make_inventory()
        p = next(t for t in inv['selected_terms'] if t['chirality']=='PL' and t['explicit_epsilon_pole'])
        p['inputform'] = x
        p['term_sha256'] = sha_bytes(x.encode('utf8'))
        p['contains_log'] = True
        with self.assertRaisesRegex(BoundaryError, 'classification inconsistent'):
            validate_and_group('A_-4', inv)

    def test_bogus_log_mass_denominator_rejected(self):
        x = triple()[2].replace('Coupling[MF, {Index[d$$1, NFlavor]}, 0]^2]', 'Coupling[MS1, {}, 0]^2]')
        with self.assertRaisesRegex(BoundaryError, 'Unsupported mass log'):
            normalize_native_kernel(x)

    def test_nonmultiplicative_log_rejected(self):
        x = triple()[2].replace('*Log[', '+Log[')
        with self.assertRaisesRegex(BoundaryError, 'multiplicative'):
            normalize_native_kernel(x)

    def test_unexpected_physical_certificate_rejected(self):
        inv = make_inventory()
        inv['MSbar_subtraction_certified'] = True
        with self.assertRaisesRegex(BoundaryError, 'unrenormalized'):
            validate_and_group('A_-4', inv)

    def test_mutated_inputform_hash_rejected(self):
        inv = make_inventory()
        inv['selected_terms'][0]['inputform'] += ' '
        with self.assertRaisesRegex(BoundaryError, 'SHA256'):
            validate_and_group('A_-4', inv)

    def test_source_label_duplicate_rejected(self):
        inv = make_inventory()
        inv['selected_terms'][1]['source_additive_term_index'] = inv['selected_terms'][0]['source_additive_term_index']
        with self.assertRaisesRegex(BoundaryError, 'Duplicate'):
            validate_and_group('A_-4', inv)

    def test_invalid_unsupported_regulator_rejected(self):
        with self.assertRaisesRegex(BoundaryError, 'Unsupported regulator'):
            normalize_native_kernel(triple()[0] + '*ε')

    def test_batch_snapshots_declared(self):
        self.assertEqual(EXPECTED_SNAPSHOT[:16], 'e515d82d118284ba')
        self.assertEqual(len(MODELS), 16)


if __name__ == '__main__':
    unittest.main()
