"""Focused tests for the additive Batch 3 comparison; no project source mutation."""
from __future__ import annotations

import copy
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from Numerical.lo_certification.benchmark import uv_eft_compare as c

ROOT = Path(__file__).resolve().parents[1]


class TestBatch3Comparison(unittest.TestCase):
    def fixture(self):
        return c.read_json(ROOT / c.SNAPSHOT_REL)

    def test_snapshot_reference_validates(self):
        r = c.compare_all(ROOT)
        self.assertEqual(r['status'], 'COMPARISON_ONLY_NOT_PHYSICAL_MATCHING')
        self.assertEqual(r['model_count'], 2)
        self.assertIsNone(r['physical_mnu_ev'])
        self.assertIsNone(r['physical_GammaCC_normalization'])

    def test_known_A_m4_reference_residual(self):
        r = c.compare_all(ROOT)['models'][0]
        self.assertEqual(r['model'], 'A_-4')
        self.assertAlmostEqual(r['residual_after_best_complex_scale'], 0.9089276327726868, places=11)
        self.assertAlmostEqual(r['uv_offdiagonal_frobenius_fraction'], 0.8031116076976612, places=11)
        self.assertAlmostEqual(r['uv_neutrino_mass_singular_values_ev_desc'][0], 0.043600644036855436, places=11)

    def test_known_B_m1_reference_residual(self):
        r = c.compare_all(ROOT)['models'][1]
        self.assertEqual(r['model'], 'B_-1')
        self.assertAlmostEqual(r['residual_after_best_complex_scale'], 0.886969835439317, places=11)
        self.assertAlmostEqual(r['uv_offdiagonal_frobenius_fraction'], 0.7947401905417403, places=11)

    def test_global_complex_factor_witness(self):
        x = np.array([[1+2j, 0.25j, 0], [0.25j, 2-1j, 0], [0, 0, 3]], complex)
        z = 2.5*np.exp(1j*0.3)
        rec = {'origin': 'synthetic', 'uv_matrix_gev': {'real': (z*x/1e9).real.tolist(), 'imag': (z*x/1e9).imag.tolist()},
               'eft_matrix_ev_per_N': {'real': x.real.tolist(), 'imag': x.imag.tolist()},
               'uv_masses_ev': np.linalg.svd(z*x, compute_uv=False).tolist(),
               'eft_masses_ev_per_abs_N': np.linalg.svd(x, compute_uv=False).tolist(),
               'uv_input_scale_gev': 1e7, 'mu_f_gev': 1e5, 'mu_s_gev': 1e3,
               'fermion_boundary_flavour_support': {'total_independent_LLSS_coordinates': 2,
                        'lepton_flavour_offdiagonal_coordinates': 1,
                        'lepton_flavour_diagonal_coordinates': 1}}
        r = c.compare_one('A_-4', rec)
        self.assertLess(r['residual_after_best_complex_scale'], 1e-14)
        self.assertAlmostEqual(complex(**r['best_complex_scale_reported_axes']).real, z.real, places=12)

    def test_pinned_boundary_flavour_support(self):
        results = c.compare_all(ROOT)['models']
        self.assertEqual([r['fermion_boundary_flavour_support']['total_independent_LLSS_coordinates']
                          for r in results], [60, 108])
        self.assertTrue(all(r['fermion_boundary_contains_only_same_flavour_lepton_pairs']
                            for r in results))

    def test_local_flavour_support_reconstructs(self):
        source = {'model': 'A_-4',
                  'status': 'NumericalLLSSBySourceVerified_SpinorNormalizationOpen',
                  'spinor_normalization': {'certified': False},
                  'summed_independent_symmetric_coordinates': [
                      {'coordinate': [1, 1, 1, 1, 5, 7]},
                      {'coordinate': [1, 1, 2, 2, 5, 9]}]}
        counts = c.flavour_support_from_boundary(source, 'A_-4')
        self.assertEqual(counts['lepton_flavour_offdiagonal_coordinates'], 1)

    def test_wrong_spectrum_rejected(self):
        d = self.fixture()['models']['A_-4']
        r = c._record(d['uv'], d['eft'], origin='test')
        r['uv_masses_ev'] = [1.0, 1.0, 1.0]
        with self.assertRaisesRegex(c.ComparisonError, 'singular masses disagree'):
            c.compare_one('A_-4', r)

    def test_wrong_source_certification_rejected(self):
        d = copy.deepcopy(self.fixture()['models']['A_-4'])
        d['eft']['physical_C5_certified'] = True
        with self.assertRaisesRegex(c.ComparisonError, 'status invalid'):
            c._record(d['uv'], d['eft'], origin='test')

    def test_non_symmetric_rejected(self):
        d = self.fixture()['models']['A_-4']
        r = c._record(d['uv'], d['eft'], origin='test')
        r['uv_matrix_gev'] = copy.deepcopy(r['uv_matrix_gev'])
        r['uv_matrix_gev']['real'][0][1] *= -1
        with self.assertRaisesRegex(c.ComparisonError, 'not complex symmetric'):
            c.compare_one('A_-4', r)

    def test_local_mode_rechecks_full_format_reports(self):
        snap = self.fixture()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            summary_path = root / c.LOCAL_EFT_REL
            summary_path.parent.mkdir(parents=True)
            summary = {'format': 't3-conditional-all16-consolidated-batch2-v1',
                       'status': 'All16ConditionalC5Exploratory',
                       'physical_C5_certified': False,
                       'absolute_GammaCC_to_Weyl_factor': None,
                       'complete_one_loop_neutrino_mass_prediction': False,
                       'models': []}
            for model in c.MODELS:
                rec = snap['models'][model]
                short = rec['eft']
                full = {'model': model, 'status': short['status'],
                        'physical_C5_certified': False,
                        'absolute_GammaCC_to_Weyl_factor': None,
                        'complete_one_loop_prediction': False,
                        'mu_f_gev': short['mu_f_gev'], 'mu_s_gev': short['mu_s_gev'],
                        'conditional_mnu_per_N_eV': short['matrix_ev_per_N'],
                        'mass_singular_values_eV_per_abs_N': short['singular_masses_ev_per_abs_N']}
                summary['models'].append(full)
                uv = rec['uv']
                uv_path = root / c.LOCAL_UV_DIR / f'{model}.json'
                uv_path.parent.mkdir(parents=True, exist_ok=True)
                uv_path.write_text(json.dumps({
                    'model': model, 'status': uv['status'],
                    'physical_source_normalisation_certified': False,
                    'complex_symmetric_matrix_gev': uv['matrix_gev'],
                    'takagi_masses_ev_from_independent_singular_values': uv['singular_masses_ev'],
                    'input_scale_gev': uv['input_scale_gev']}))
                support = snap['fermion_boundary_flavour_support'][model]
                length = support['total_independent_LLSS_coordinates']
                boundary_path = root / c.BOUNDARY_REL / c.BOUNDARY_SLUG[model] / 'boundary.json'
                boundary_path.parent.mkdir(parents=True, exist_ok=True)
                boundary_path.write_text(json.dumps({
                    'model': model, 'status': 'NumericalLLSSBySourceVerified_SpinorNormalizationOpen',
                    'spinor_normalization': {'certified': False},
                    'summed_independent_symmetric_coordinates': [
                        {'coordinate': [1, 1, 1, 1, 5, 7]} for _ in range(length)]}))
            summary_path.write_text(json.dumps(summary))
            records = c.compare_all(root, 'local')
            self.assertEqual(records['sources']['source_mode'], 'local')
            self.assertAlmostEqual(records['models'][0]['residual_after_best_complex_scale'],
                                   0.9089276327726868, places=11)

    def test_missing_local_sources_fails_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(c.ComparisonError, 'Missing source file'):
                c.compare_all(Path(tmp), 'local')

    def test_immutable_output(self):
        data = c.compare_all(ROOT)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            dest = c.emit(root, data, 'test01')
            self.assertEqual(json.loads((dest / 'summary.json').read_text())['physical_certification'], False)
            self.assertIn('DIAGNOSTIC ONLY', (dest / 'comparison.md').read_text())
            with self.assertRaisesRegex(c.ComparisonError, 'already exists'):
                c.emit(root, data, 'test01')

    def test_directory_traversal_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(c.ComparisonError, 'Unsafe run ID'):
                c.emit(Path(tmp), c.compare_all(ROOT), '../escape')

    def test_wrong_scale_rejected(self):
        d = copy.deepcopy(self.fixture()['models']['B_-1'])
        d['eft']['mu_s_gev'] = 2e5
        with self.assertRaisesRegex(c.ComparisonError, 'Invalid scales'):
            c._record(d['uv'], d['eft'], origin='test')

    def test_wrong_model_inventory_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path = root / c.SNAPSHOT_REL
            path.parent.mkdir(parents=True)
            obj = self.fixture()
            del obj['models']['B_-1']
            path.write_text(json.dumps(obj))
            with self.assertRaisesRegex(c.ComparisonError, 'inventory'):
                c.compare_all(root)


if __name__ == '__main__':
    unittest.main()
