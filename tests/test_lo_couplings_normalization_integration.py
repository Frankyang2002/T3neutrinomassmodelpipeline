"""All-16 real-checkout tensor/registry integration checks, read-only.

Skipped outside a complete local repository. The run CLI performs the additional
source-drift structural recheck before writing, but this test does not overwrite
existing outputs or invoke either RGE solver.
"""
from pathlib import Path
import unittest

from Numerical.lo_certification.couplings.portal_normalization import (
    portal_species_audit,construct_rgbeta_portal_background)
from Numerical.lo_certification.couplings.quartic_audit import rgbeta_potential_declarations

ROOT=Path(__file__).resolve().parents[1]
REQUIRED=ROOT/'Numerical/joint_scalar_matching/source_quartic_v4.py'

@unittest.skipUnless(REQUIRED.is_file(),'requires current full local T3 repository')
class TestRealLocalAll16PortalSpans(unittest.TestCase):
    def test_all16_native_cached_quartics_have_explicit_portal_span_status(self):
        from Numerical.joint_scalar_matching.adapter import ALL_MODELS
        from Numerical.joint_scalar_matching.fermion_boundary_bridge import load_source,slug
        from Numerical.joint_scalar_matching.numerical_benchmark_v5 import read_v4_cached_channels
        from Numerical.joint_scalar_matching.source_quartic_v4 import source_seed_path
        from Numerical.lo_certification.couplings.batch import (read_json, observed_running_states)
        release=ROOT/'output/general_observable/lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f'
        cache=ROOT/'output/general_observable/lo_joint_scalar_matching/source_quartic_v4/native_quartic_channels'
        for model in ALL_MODELS:
            with self.subTest(model=model):
                source=load_source(release/slug(model)/'boundary.json',expected_model=model)
                seed=source_seed_path(ROOT/'output/full16_symbolic',model)
                ch=read_v4_cached_channels(cache/f'{slug(model)}.json',seed,source)
                result=[portal_species_audit(ch,s) for s in (1,2)]
                self.assertEqual(len(result),2)
                for row in result:
                    self.assertIn(row['status'],(
                        'exact_full_tensor_span','exact_partial_tensor_span',
                        'unresolved_generator_or_complex_scalar_basis'))
                    self.assertIsInstance(row['native_unit_to_rgbeta_invariant_coefficients'],list)
                data=observed_running_states(
                    ROOT/f'output/general_observable/all16_v1/{slug(model)}/running.json',model)
                bg,arrays=construct_rgbeta_portal_background(
                    ch,data['eft1_at_fermion_threshold'],data['eft1_at_scalar_threshold'])
                self.assertFalse(bg['safe_for_complete_seven_term_or_joint_matching'])
                self.assertEqual(arrays['muS_portal_sum_tensor'].shape[0:2],(4,4))

    def test_local_matchete_rgbeta_source_declarations(self):
        doc=rgbeta_potential_declarations(ROOT)
        self.assertTrue(all(doc['portal_declarations_present'].values()))
        self.assertEqual(len(doc['matchete_portal_source_sha256']),64)

if __name__=='__main__':unittest.main()
