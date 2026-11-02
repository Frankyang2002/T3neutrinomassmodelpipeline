"""Local-checkout all16 exact native-to-completed portal and extra SU2 tests."""
from __future__ import annotations
from pathlib import Path
import unittest

from Numerical.lo_certification.couplings.portal_completion import (
    exact_portal_completion,completed_portal_numeric)
from Numerical.lo_certification.couplings.extra_hhss_geometry import model_extra_geometry
from Numerical.lo_certification.couplings.batch import (
    BASELINE_HHSS_REGISTRY,read_json,observed_running_states,_registry_hhss_names)

ROOT=Path(__file__).resolve().parents[1]
READY=(ROOT/'Numerical/joint_scalar_matching/source_quartic_v4.py').is_file()

@unittest.skipUnless(READY,'Requires full local repository with source-backed caches')
class TestAll16RealCompletion(unittest.TestCase):
    @staticmethod
    def models():
        from Numerical.joint_scalar_matching.adapter import ALL_MODELS
        return ALL_MODELS

    @staticmethod
    def load(model):
        from Numerical.joint_scalar_matching.fermion_boundary_bridge import load_source,slug
        from Numerical.joint_scalar_matching.numerical_benchmark_v5 import read_v4_cached_channels
        from Numerical.joint_scalar_matching.source_quartic_v4 import source_seed_path
        s=slug(model)
        release=ROOT/'output/general_observable/lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f'
        boundary=load_source(release/s/'boundary.json',expected_model=model)
        seed=source_seed_path(ROOT/'output/full16_symbolic',model)
        cache=ROOT/'output/general_observable/lo_joint_scalar_matching/source_quartic_v4/native_quartic_channels'
        ch=read_v4_cached_channels(cache/f'{s}.json',seed,boundary)
        running=observed_running_states(ROOT/'output/general_observable/all16_v1'/s/'running.json',model)
        return ch,running

    def test_all16_native_plus_added_portals_complete_full_exact_tensors(self):
        total=0
        for model in self.models():
            with self.subTest(model=model):
                ch,run=self.load(model)
                a=[exact_portal_completion(ch,i) for i in (1,2)]
                self.assertEqual(sum(x['rgbeta_direction_count'] for x in a),
                                 sum(len(x['augmented_symbols']) for x in a))
                total+=sum(x['completion_count'] for x in a)
        self.assertEqual(total,14)

    def test_all16_same_scale_numerical_portals_without_missing_values(self):
        for model in self.models():
            with self.subTest(model=model):
                ch,run=self.load(model)
                doc,arr=completed_portal_numeric(ch,{
                    'muF':run['eft1_at_fermion_threshold'],
                    'muS':run['eft1_at_scalar_threshold']})
                for scale in ('muF','muS'):
                    if doc['scales'][scale]['missing_rgbeta_couplings']:
                        self.assertNotIn(f'{scale}_portal_sum',arr)
                        self.assertEqual(doc['scales'][scale]['status'],
                                         'blocked_missing_same_scale_portal_coupling')
                    else:
                        self.assertIn(f'{scale}_portal_sum',arr)
                        self.assertIsNotNone(doc['scales'][scale]['tensor_relative_reconstruction_residual'])
                self.assertFalse(doc['physical_complete_EFT1_input'])

    def test_all16_extra_registry_rank_and_full_derivative_export(self):
        reg=read_json(ROOT/'RGE/general/generated/general_rgbeta_quartic_registry.json')
        total=0
        for model in self.models():
            with self.subTest(model=model):
                row=reg['models'][model]
                audit=model_extra_geometry(model,row,include_components=True)
                n=audit['registry_HHSS_extra_symbols']
                self.assertEqual(n,len(_registry_hhss_names(row)))
                self.assertGreaterEqual(n,len(BASELINE_HHSS_REGISTRY[model]))
                total+=n
                for sector in audit['monomial_sectors']:
                    self.assertIsNone(sector['physical_mapping_to_named_RGBeta_coefficients'])
                    self.assertEqual(len(sector['exact_reference_real_component_tensor_basis']),
                                     sector['independent_complex_invariant_count'])
        self.assertGreaterEqual(total,24)

if __name__=='__main__':
    unittest.main()
