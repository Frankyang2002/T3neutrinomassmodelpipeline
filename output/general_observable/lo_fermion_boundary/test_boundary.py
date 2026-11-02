"""Offline unit tests and optional read-only live 16-model input check.

Project destination: output/general_observable/lo_fermion_boundary/test_boundary.py
"""
import unittest
from pathlib import Path
import subprocess
import sys

import build_boundary as b


class BoundaryTests(unittest.TestCase):
    def test_all_model_count_and_labels(self):
        self.assertEqual(len(b.MODELS), 16)
        self.assertEqual(len(set(b.MODELS)), 16)
        self.assertIn('B_-1', b.MODELS)
        self.assertIn('C_-1', b.MODELS)

    def test_orbit_cardinalities(self):
        self.assertEqual(len(b.orbit((1, 1, 1, 1, 5, 5))), 1)
        self.assertEqual(len(b.orbit((1, 1, 1, 1, 5, 6))), 2)
        self.assertEqual(len(b.orbit((1, 1, 2, 3, 5, 5))), 2)
        self.assertEqual(len(b.orbit((1, 1, 2, 3, 5, 6))), 4)
        self.assertEqual(b.canon((2, 3, 1, 1, 6, 5)), (1, 1, 2, 3, 5, 6))

    def test_source_polynomial_does_not_double_identical_legs(self):
        raw = {(1, 1, 1, 1, 5, 6): 2j}
        rows = b.rows_from_aligned(raw)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]['orbit_size'], 2)
        self.assertEqual(rows[0]['coefficient_imag'], 1)
        self.assertEqual(rows[0]['coefficient_real'], 0)
        b.compare_rows(rows, rows, 'fixture')

    def test_distinct_sources_preserved_before_summation(self):
        a = b.rows_from_aligned({(1, 1, 1, 2, 5, 6): 1 + 2j})
        c = b.rows_from_aligned({(1, 1, 1, 2, 5, 6): -0.25 - 1j})
        self.assertEqual(len(a), len(c))
        summed = b.rows_from_aligned({(1, 1, 1, 2, 5, 6): 0.75 + 1j})
        self.assertAlmostEqual(b.as_map(a)[(1, 1, 1, 2, 5, 6)]
                               + b.as_map(c)[(1, 1, 1, 2, 5, 6)],
                               b.as_map(summed)[(1, 1, 1, 2, 5, 6)])

    def test_gamma_cc_projectors_detected_without_guessing_sign(self):
        text = 'Transp[Field[l, Fermion, {Index[i, SU2L[fund]], Index[p, Flavor]}, 0]] GammaCC Proj[-1]'
        self.assertEqual(b.chirality(text), 'PL')
        conj = 'Transp[Bar[Field[l, Fermion, {Index[i, SU2L[fund]], Index[p, Flavor]}, 0]]] GammaCC Proj[1]'
        self.assertEqual(b.chirality(conj), 'PR')
        with self.assertRaises(ValueError):
            b.chirality(text.replace('Proj[-1]', 'Proj[1]'))
        with self.assertRaises(ValueError):
            b.chirality(text + ' Proj[1]')

    def test_flavour_axis_transposition(self):
        pre = 'Field[l, Fermion, {Index[xx, SU2L[fund]], Index[p, Flavor]}, 0]'
        post = 'Field[l, Fermion, {Index[yy, SU2L[fund]], Index[q, Flavor]}, 0]'
        source = {'flavour_matrix_axes_external_indices': ['q', 'p'],
                  'factors': [{'external_index': 'q'}, {'external_index': 'p'}],
                  'operator_template_unit_flavor_inputform':
                  'Transp[' + pre + '] GammaCC Proj[-1] ' + post}
        operator = {'flavour_matrix_axes_external_indices': ['q', 'p']}
        self.assertEqual(b.source_alignment(operator, source), [1, 0])

    def test_no_unverified_spinor_factor_in_code(self):
        # No argument exists that can cause an unknown physical factor to be assigned.
        source = Path(b.__file__).read_text(encoding='utf-8')
        self.assertIn("'gamma_cc_to_master_weyl_factor': None", source)
        self.assertIn("'gamma_cc_to_master_weyl_sign': None", source)
        self.assertIn("'one_loop_F_hard_matching_included': False", source)
        self.assertIn('not destination.exists()', source)


if __name__ == '__main__':
    live = '--live' in sys.argv
    sys.argv = [v for v in sys.argv if v != '--live']
    result = unittest.main(exit=False, verbosity=2)
    if not result.result.wasSuccessful():
        raise SystemExit(1)
    if live:
        rc = subprocess.call([sys.executable, str(Path(__file__).with_name('build_boundary.py')),
                              '--check-only'])
        raise SystemExit(rc)
