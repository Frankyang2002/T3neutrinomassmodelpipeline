"""Strict native hard-threshold evidence tests: no Wolfram needed."""
import json
from pathlib import Path
import tempfile
import unittest

from Numerical.lo_certification.spinor.hard_threshold import (
    model_command, inspect_raw, _owned_filename, blocked_record, MIXED,
    inspect_native_phase, PHASE_KEYS, native_environment, phase_environment, FLAGS)
from Numerical.lo_certification.spinor.boundary import BoundaryError


class HardThresholdTests(unittest.TestCase):
    def test_14_general_2_mixed(self):
        from Numerical.lo_certification.spinor.boundary import MODELS
        self.assertEqual(len(MODELS), 16)
        self.assertEqual(len([m for m in MODELS if m not in MIXED]), 14)
        self.assertEqual(len(MIXED), 2)

    def test_general_one_loop_command(self):
        with tempfile.TemporaryDirectory() as t:
            root = Path(t)
            cmd = model_command(root, 'A_-4', root / 'out/A_m4')
            self.assertEqual(cmd[-2:], ['A_-4', '1'])
            self.assertEqual(cmd[0:2], ['wolframscript', '-file'])
            self.assertIn('RunT3GeneralFThreshold.wl', cmd[2])

    def test_mixed_command_requires_actual_native_mass_input(self):
        with tempfile.TemporaryDirectory() as t:
            with self.assertRaises(BoundaryError):
                model_command(Path(t), 'B_-1', Path(t)/'out/B_m1')

    def test_general_native_validated_not_certified(self):
        with tempfile.TemporaryDirectory() as t:
            d = Path(t)
            obj = dict(model='E_-2', loop_order=1,
                       status='FThresholdTreeWilsonExported',
                       raw_one_loop_F_matching_included=True,
                       hard_one_loop_MSbar_certified=False,
                       integrated_fields=['F'], light_bsm_fields=['S1', 'S2'],
                       heavy_flavour_count=3)
            (d/'native_status.json').write_text(json.dumps(obj))
            (d/'data').mkdir()
            (d/'data/eft1_after_F_one_loop_inputform.txt').write_text('hbar*GammaCC*Field[H, Scalar, {}, {}]')
            r = inspect_raw('E_-2', d)
            self.assertIsNone(r['physical_wilson_coefficient'])
            self.assertIsNone(r['explicit_hard_LLHH_candidate_terms'])
            self.assertTrue(r['raw_loop_nonzero'])
            self.assertTrue(r['raw_loop_contains_GammaCC'])
            self.assertFalse(r['hard_sector_renormalized'])

    def test_general_tree_level_rejected(self):
        with tempfile.TemporaryDirectory() as t:
            d = Path(t)
            obj = dict(model='E_-2', loop_order=0,
                       status='FThresholdTreeWilsonExported',
                       raw_one_loop_F_matching_included=False,
                       hard_one_loop_MSbar_certified=False,
                       integrated_fields=['F'], light_bsm_fields=['S1', 'S2'],
                       heavy_flavour_count=3)
            (d/'native_status.json').write_text(json.dumps(obj))
            with self.assertRaises(BoundaryError):
                inspect_raw('E_-2', d)

    def test_mixed_raw_zero_does_not_promote_physical(self):
        with tempfile.TemporaryDirectory() as t:
            d = Path(t)
            model = 'B_-1'
            data = dict(model=model,
                        format='t3-mass-eigen-fermion-first-eft1-oneloop-v1',
                        status='FermionFirstOneLoopRawSectorsExtractedNotRenormalized',
                        matching_EFT_order=5, matching_loop_order=1,
                        matching_approach='FermionFirstTreeAndOneLoopScalarsRemainLight',
                        integrated_heavy_fields=['NewFermion'],
                        scalar_light_fields=['NewScalarX', 'NewScalarY'],
                        one_loop_matching_completed=True,
                        loop_difference_hbar_present=True,
                        loop_difference_fermion_absent=True,
                        raw_hard_llhh_scheme_renormalized=False,
                        raw_hard_llhh_pole_rge_certified=False,
                        physical_C5_computed=False,
                        raw_hard_llhh_candidate_term_count=0,
                        one_loop_difference_file='loop.inputform.txt',
                        raw_hard_llhh_file='hard.inputform.txt')
            (d/'loop.inputform.txt').write_text('hbar*C')
            (d/'hard.inputform.txt').write_text('0')
            (d/'B_m1_native_one_loop.json').write_text(json.dumps(data))
            r = inspect_raw(model, d)
            self.assertEqual(r['explicit_hard_LLHH_candidate_terms'], 0)
            self.assertFalse(r['hard_sector_renormalized'])
            self.assertIsNone(r['physical_wilson_coefficient'])

    def test_mixed_hard_sector_source_nonzero(self):
        with tempfile.TemporaryDirectory() as t:
            d = Path(t)
            data = dict(model='C_-1',
                        format='t3-mass-eigen-fermion-first-eft1-oneloop-v1',
                        status='FermionFirstOneLoopRawSectorsExtractedNotRenormalized',
                        matching_EFT_order=5, matching_loop_order=1,
                        matching_approach='FermionFirstTreeAndOneLoopScalarsRemainLight',
                        integrated_heavy_fields=['NewFermion'],
                        scalar_light_fields=['NewScalarX','NewScalarY'],
                        one_loop_matching_completed=True,
                        loop_difference_hbar_present=True,
                        loop_difference_fermion_absent=True,
                        raw_hard_llhh_scheme_renormalized=False,
                        raw_hard_llhh_pole_rge_certified=False,
                        physical_C5_computed=False,
                        raw_hard_llhh_candidate_term_count=2,
                        one_loop_difference_file='loop.inputform.txt',
                        raw_hard_llhh_file='hard.inputform.txt')
            (d/'loop.inputform.txt').write_text('hbar*C')
            (d/'hard.inputform.txt').write_text('GammaCC*H*H')
            (d/'C_m1_native_one_loop.json').write_text(json.dumps(data))
            r = inspect_raw('C_-1', d)
            self.assertEqual(r['explicit_hard_LLHH_candidate_terms'], 2)
            self.assertFalse(r['hard_sector_renormalized'])

    def test_unsafe_native_file_path_rejected(self):
        with tempfile.TemporaryDirectory() as t:
            for x in ('../escape.txt', '/outside/file', '..\\escape.txt', 'C:\\foo\\bad.txt'):
                with self.assertRaises(BoundaryError):
                    _owned_filename(Path(t), x)

    def test_native_driver_environment_isolated_by_model(self):
        parent = {key: 'bad' for key in FLAGS}
        ordinary = native_environment('A_-4', parent)
        mixed = native_environment('B_-1', parent)
        phase = phase_environment(parent)
        self.assertEqual([ordinary[k] for k in FLAGS], ['1']*4)
        self.assertEqual([mixed[k] for k in FLAGS], ['0']*4)
        self.assertEqual([phase[k] for k in FLAGS], ['1','1','0','1'])
        self.assertEqual([parent[k] for k in FLAGS], ['bad']*4)

    def test_native_abstract_CConj_never_sets_absolute_sign(self):
        with tempfile.TemporaryDirectory() as t:
            d = Path(t)/'native_A_p0.json'
            obj = {
                'format': 't3-workstream-a-a0-native-spinor-algebra-v1',
                'model': 'A_p0',
                'status': 'NativeMatcheteAbstractChargeConjugationAndHcPassed',
                'checks': {k: True for k in PHASE_KEYS},
            }
            d.write_text(json.dumps(obj))
            r = inspect_native_phase(d)
            self.assertEqual(r['native_checks_passed'], 6)
            self.assertIsNone(r['GammaCC_to_master_Weyl_factor'])
            obj['checks']['cconj_flips_left_projector'] = False
            d.write_text(json.dumps(obj))
            with self.assertRaises(BoundaryError):
                inspect_native_phase(d)

    def test_blocked_result_is_never_physical(self):
        r = blocked_record('B_-1','mocked native failure',5)
        self.assertEqual(r['status'],'BLOCKED')
        self.assertIsNone(r['physical_wilson_coefficient'])
        self.assertEqual(r['native_exit_code'],5)


if __name__ == '__main__':
    unittest.main()
