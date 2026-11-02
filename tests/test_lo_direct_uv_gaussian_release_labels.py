"""Regression: release folder uses B_m1 but boundary JSON uses B_-1.

These fixtures reproduce the documented release metadata convention for all
16 model points.  No UV/EFT matching or source normalization is inferred.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from Numerical.lo_certification.direct_uv.certification import CLASSES, all_models
from Numerical.lo_certification.direct_uv.gaussian_amplitude import (
    AmplitudeError, _eft1_boundary_gate, batch,
)


RELEASE = 'source_021dcc1e8fe8ba0f'
F_THRESHOLD_GEV = 1.0e5
UV_INPUT_GEV = 1.0e7


def make_release(root: Path, model: str, *, recorded_model: str | None = None,
                 n_f: str = 'all at mu_F', n_s: str = 'both retained in EFT1') -> Path:
    """Source-shaped record; directory slug intentionally differs from model."""
    slug = model.replace('_-', '_m').replace('_+', '_p')
    p = (root / 'output' / 'general_observable' / 'lo_fermion_boundary'
         / 'releases' / RELEASE / slug / 'boundary.json')
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps({
        'format': 't3-lo-fermion-boundary-source-resolved-provisional-v1',
        'model': model if recorded_model is None else recorded_model,
        'status': 'NumericalLLSSBySourceVerified_SpinorNormalizationOpen',
        'mu_F_gev': F_THRESHOLD_GEV,
        'field_content': {
            'heavy_fermions_integrated_out': n_f,
            'heavy_scalar_multiplets': n_s,
        },
        'spinor_normalization': {
            'certified': False,
            'gamma_cc_to_master_weyl_factor': None,
            'gamma_cc_to_master_weyl_sign': None,
        },
    }) + '\n', encoding='utf-8')
    return p


class BoundaryReleaseLabels(unittest.TestCase):
    def test_real_release_label_convention_all16(self):
        with tempfile.TemporaryDirectory() as t:
            root = Path(t)
            points = [m for m, _ in all_models()]
            self.assertEqual(len(points), 16)
            for model in points:
                path = make_release(root, model)
                self.assertNotEqual(path.parent.name, model)
            for model in points:
                with self.subTest(model=model):
                    gate = _eft1_boundary_gate(root, model, UV_INPUT_GEV)
                    self.assertEqual(gate['status'],
                        'NOT_COMPARED_SCALE_MISMATCH_AND_SPINOR_NORMALIZATION_OPEN')
                    self.assertFalse(gate['same_renormalization_scale'])
                    self.assertFalse(gate['compared_LLSS_tensors'])
                    self.assertFalse(gate['native_spinor_normalization_certified_in_release'])
                    self.assertEqual(gate['uv_input_mu_gev'], UV_INPUT_GEV)
                    self.assertEqual(gate['eft1_threshold_mu_F_gev'], F_THRESHOLD_GEV)

    def test_wrong_old_slug_inside_release_rejected(self):
        with tempfile.TemporaryDirectory() as t:
            root = Path(t)
            make_release(root, 'B_-1', recorded_model='B_m1')
            with self.assertRaisesRegex(AmplitudeError, 'model label disagrees'):
                _eft1_boundary_gate(root, 'B_-1', UV_INPUT_GEV)

    def test_swapped_signed_model_record_rejected(self):
        with tempfile.TemporaryDirectory() as t:
            root = Path(t)
            make_release(root, 'A_-4', recorded_model='A_-2')
            with self.assertRaisesRegex(AmplitudeError, 'model label disagrees'):
                _eft1_boundary_gate(root, 'A_-4', UV_INPUT_GEV)

    def test_threshold_ordering_remains_fail_closed(self):
        with tempfile.TemporaryDirectory() as t:
            root = Path(t)
            make_release(root, 'C_-1', n_f='one heavy fermion only')
            with self.assertRaisesRegex(AmplitudeError, 'partial heavy multiplet'):
                _eft1_boundary_gate(root, 'C_-1', UV_INPUT_GEV)
            make_release(root, 'C_-1', n_s='one scalar integrated already')
            with self.assertRaisesRegex(AmplitudeError, 'scalar ordering'):
                _eft1_boundary_gate(root, 'C_-1', UV_INPUT_GEV)

    def test_equal_scales_do_not_claim_physical_match(self):
        with tempfile.TemporaryDirectory() as t:
            root = Path(t)
            make_release(root, 'E_+2')
            gate = _eft1_boundary_gate(root, 'E_+2', F_THRESHOLD_GEV)
            self.assertTrue(gate['same_renormalization_scale'])
            self.assertFalse(gate['compared_LLSS_tensors'])
            self.assertEqual(gate['status'],
                'NOT_COMPARED_PHYSICAL_SPINOR_AND_CG_DICTIONARY_OPEN')

    def test_full_all16_batch_reaches_correct_release_gate(self):
        """Test the formerly failing branch, not just the gate in isolation.

        The Gaussian tensor evaluator is stubbed because its mathematical
        correctness is tested separately by the existing 60-test suite.
        Input hashes, model labels, representation and scales remain real
        checks. This catches the original all-16 false-negative regression.
        """
        with tempfile.TemporaryDirectory() as t:
            root = Path(t)
            inputs = (root / 'output/general_observable/lo_direct_one_loop'
                      / 'provisional_inputs')
            inputs.mkdir(parents=True, exist_ok=True)
            for model, benchmark_key in all_models():
                make_release(root, model)
                kind, alpha = model[0], int(model[2:])
                d1, d2, df, _ = CLASSES[kind]
                bench = (root / 'Numerical/benchmarks/fitted_reference'
                         / (benchmark_key + '.json'))
                bench.parent.mkdir(parents=True, exist_ok=True)
                bench.write_text(json.dumps({
                    'representation': {'d_s1': d1, 'd_s2': d2,
                                       'd_f': df, 'alpha': alpha},
                    'base_state': {'mu_gev': UV_INPUT_GEV},
                }), encoding='utf-8')
                (inputs / (model + '.json')).write_text(json.dumps({
                    'model': model,
                    'renormalization_scale_gev': UV_INPUT_GEV,
                    'provenance': {
                        'benchmark_sha256': hashlib.sha256(bench.read_bytes()).hexdigest(),
                    },
                }), encoding='utf-8')
            with patch(
                'Numerical.lo_certification.direct_uv.gaussian_amplitude.guard_sources',
                return_value={'source_sha256': {}},
            ), patch(
                'Numerical.lo_certification.direct_uv.gaussian_amplitude.optional_evidence',
                return_value={},
            ), patch(
                'Numerical.lo_certification.direct_uv.gaussian_amplitude.gaussian_tensors',
                return_value={'status': 'formal_fixture_only'},
            ) as gaussian:
                result = batch(root, include_tensors=False)
            self.assertEqual(gaussian.call_count, 16)
            self.assertEqual(result['tree_amplitude_computed_models'], 16)
            self.assertEqual(result['blocked_models'], 0)
            self.assertEqual(result['physical_certified_models'], 0)
            self.assertEqual(len(result['models']), 16)
            for row in result['models']:
                with self.subTest(model=row['model']):
                    self.assertEqual(row['status'], 'conditional_tree_LLSS_completed')
                    self.assertEqual(row['eft1_release_gate']['status'],
                        'NOT_COMPARED_SCALE_MISMATCH_AND_SPINOR_NORMALIZATION_OPEN')
                    self.assertFalse(row['eft1_release_gate']['compared_LLSS_tensors'])

    def test_missing_boundary_is_not_a_hard_gaussian_failure(self):
        with tempfile.TemporaryDirectory() as t:
            gate = _eft1_boundary_gate(Path(t), 'D_+2', UV_INPUT_GEV)
            self.assertEqual(gate['status'], 'no_source_fermion_threshold_release')
            self.assertFalse(gate['compared_LLSS_tensors'])


if __name__ == '__main__':
    unittest.main()
