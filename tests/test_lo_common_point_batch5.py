"""Batch 5: source-controlled UV versus conditional EFT diagnostic tests."""
from __future__ import annotations
from copy import deepcopy
from pathlib import Path
import tempfile
import unittest
import sys

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0,str(ROOT))
from Numerical.lo_certification.benchmark import common_point_batch5 as m


def complex_matrix(array):
    a=np.asarray(array,float)
    return {'real':a.tolist(),'imag':np.zeros_like(a).tolist()}

def encoded(array):
    return [[{'real':float(z.real),'imag':float(z.imag)} for z in row]
            for row in np.asarray(array,complex)]

def fixture():
    y1=np.array([[1,2,3],[-2,1,2],[1,3,2]],float)*1e-5
    y2=np.array([[2,1,-3],[1,2,1],[-1,2,3]],float)*1e-5
    mf=np.diag([100000,110000,120000]).astype(float)
    ord={'y1':complex_matrix(y1),'y2':complex_matrix(y2),
         'MF':complex_matrix(mf),'mS1Sq':1306122.,'mS2Sq':1e6,
         'lambdaH1':0.02,'lambdaH2':0.03,'lambdaT3':{'real':0.1,'imag':0.0},
         'lambdaH2Adj':0.,'lambdaS2Adj':0.}
    frozen={'base_state':{'mu_gev':1e7,'ordinary':ord},
            'numerical':{'vev_gev':246.22}}
    scalar_f={'mS1Sq':1400000.,'mS2Sq':1200000.,'lambdaH1':.04,
              'lambdaH2':.05,'lambdaT3':{'real':.12,'imag':0.},
              'lambdaH2Adj':.07,'lambdaS2Adj':.05}
    scalar_s={'mS1Sq':1550000.,'mS2Sq':1450000.,'lambdaH1':.041,
              'lambdaH2':.052,'lambdaT3':{'real':.13,'imag':0.},
              'lambdaH2Adj':.03,'lambdaS2Adj':.02}
    run={'uv_at_fermion_threshold':{
         **scalar_f,'y1':encoded(y1*1.15),'y2':encoded(y2*.92),'MF':encoded(mf*1.08)},
         'eft1_at_scalar_threshold':scalar_s}
    return frozen,run

class TestCommonPointBatch5(unittest.TestCase):
    def test_01_convert_evolved_matrix_without_guessing(self):
        f,r=fixture()
        value=m.to_frozen_style(r['uv_at_fermion_threshold']['y1'],'y1',
                                   f['base_state']['ordinary']['y1'])
        self.assertEqual(value['real'][0][0],1.15e-5)
        self.assertEqual(len(value['imag']),3)

    def test_02_wrong_matrix_encoding_rejected(self):
        with self.assertRaises((ValueError,TypeError,KeyError)):
            m.from_encoded([[{'real':4.}]],'y')

    def test_03_complex_scalar_round_trip(self):
        val=m.to_frozen_style({'real':.3,'imag':.05},'lambdaT3',
                               {'real':.1,'imag':0.})
        self.assertEqual(val,{'real':.3,'imag':.05})

    def test_04_real_scalar_complex_rejected(self):
        with self.assertRaisesRegex(m.Batch5Error,'real scalar expected'):
            m.to_frozen_style({'real':.2,'imag':.8},'mS1Sq',1e6)

    def test_05_frozen_scalar_stage_is_unchanged(self):
        f,r=fixture()
        x=m.state_benchmark(f,r,scalar_stage='UV')
        self.assertEqual(x['benchmark']['base_state']['ordinary']['mS1Sq'],1306122.)
        self.assertEqual(x['benchmark']['base_state']['ordinary']['lambdaT3']['real'],.1)
        self.assertEqual(x['omitted_portals'],{})
        self.assertEqual(x['fermion_parameter_scale_gev'],1e5)
        self.assertEqual(f['base_state']['ordinary']['mS1Sq'],1306122.)

    def test_06_F_scalar_stage_uses_F_values(self):
        f,r=fixture()
        x=m.state_benchmark(f,r,scalar_stage='F')
        self.assertEqual(x['benchmark']['base_state']['ordinary']['mS1Sq'],1.4e6)
        self.assertEqual(x['benchmark']['base_state']['ordinary']['lambdaT3']['real'],.12)
        self.assertEqual(x['scalar_parameter_scale_gev'],1e5)

    def test_07_S_scalar_stage_uses_eft1_muS_not_muF(self):
        f,r=fixture()
        x=m.state_benchmark(f,r,scalar_stage='S')
        self.assertEqual(x['benchmark']['base_state']['ordinary']['mS1Sq'],1.55e6)
        self.assertEqual(x['benchmark']['base_state']['ordinary']['lambdaT3']['real'],.13)
        self.assertEqual(x['scalar_parameter_scale_gev'],1e3)
        self.assertEqual(x['benchmark']['base_state']['mu_gev'],1e3)
        self.assertNotEqual(x['benchmark']['base_state']['ordinary']['mS1Sq'],1.4e6)

    def test_08_unsupported_nonzero_portal_is_recorded_and_removed_from_uv_constructor(self):
        f,r=fixture()
        x=m.state_benchmark(f,r,scalar_stage='S')
        self.assertEqual(x['omitted_portals']['lambdaH2Adj']['real'],.03)
        self.assertEqual(x['omitted_portals']['lambdaS2Adj']['real'],.02)
        self.assertEqual(x['benchmark']['base_state']['ordinary']['lambdaH2Adj'],0.)
        self.assertEqual(x['benchmark']['base_state']['ordinary']['lambdaS2Adj'],0.)
        self.assertEqual(f['base_state']['ordinary']['lambdaS2Adj'],0.)

    def test_09_fermion_matrix_same_in_all_three_scenarios(self):
        f,r=fixture()
        mats=[m.z3(m.state_benchmark(f,r,scalar_stage=s)['benchmark']
                      ['base_state']['ordinary']['MF'],'MF') for s in ('UV','F','S')]
        for z in mats[1:]:self.assertTrue(np.array_equal(z,mats[0]))
        self.assertEqual(mats[0][0,0],108000.)

    def test_10_best_complex_factor_exact_for_proportional_matrix(self):
        e=np.array([[.05,.01,0],[.01,.02,.007],[0,.007,.01]],complex)
        n=1.2-2.5j
        out=m.comparison(e*n/1e9,e)
        self.assertLess(out['residual_after_one_complex_multiplier'],1e-13)
        self.assertAlmostEqual(out['best_complex_multiplier_of_EFT_per_N_to_UV']['real'],n.real)
        self.assertAlmostEqual(out['best_complex_multiplier_of_EFT_per_N_to_UV']['imag'],n.imag)
        self.assertFalse(out['matrix_agreement_certified'])

    def test_11_mismatched_flavour_shape_residual_nonzero(self):
        e=np.diag([.05,.02,.01]).astype(complex)
        u=np.array([[.05,.01,0],[.01,.02,0],[0,0,.01]],complex)/1e9
        out=m.comparison(u,e)
        self.assertGreater(out['residual_after_one_complex_multiplier'],.1)
        self.assertAlmostEqual(out['offdiagonal_eft_norm_fraction'],0.)
        self.assertGreater(out['offdiagonal_uv_norm_fraction'],.1)

    def test_12_nonzero_but_wrong_shape_rejected(self):
        with self.assertRaisesRegex(m.Batch5Error,'invalid'):
            m.comparison(np.eye(2),np.eye(3))

    def test_13_zero_mass_matrix_not_interpreted_as_finite_normalization(self):
        with self.assertRaisesRegex(m.Batch5Error,'Zero comparison'):
            m.comparison(np.zeros((3,3)),np.eye(3))

    def test_14_nonfinite_matrix_rejected(self):
        a=np.eye(3);a[0,0]=float('nan')
        with self.assertRaisesRegex(m.Batch5Error,'invalid'):
            m.comparison(a,np.eye(3))

    def test_15_keeps_signed_factor_uncertified_by_design(self):
        e=np.diag([.05,.02,.01]).astype(complex)
        out=m.comparison(e*2/1e9,e)
        self.assertNotIn('gamma_cc_physical_factor',out)
        self.assertFalse(out['matrix_agreement_certified'])

    def test_16_all_F_first_both_S_simultaneous_only(self):
        code=Path(m.__file__).read_text(encoding='utf-8')
        self.assertIn('all_BSM_fermions_integrated_first',code)
        self.assertIn('both_scalars_integrated_together',code)
        self.assertNotIn('scalar1_separate_threshold',code)

    def test_17_input_scenario_invalid_stage_rejected(self):
        f,r=fixture()
        with self.assertRaisesRegex(m.Batch5Error,'Unknown scalar stage'):
            m.state_benchmark(f,r,scalar_stage='NEVER')

    def test_18_unicode_line_endings_sha_pure(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t)/'x.json';p.write_bytes(b'{\r\n"a":1\r\n}\r\n')
            q=Path(t)/'y.json';q.write_bytes(b'{\n"a":1\n}\n')
            self.assertNotEqual(m.digest(p),m.digest(q))
            self.assertEqual(m.sha256_lf(p),m.sha256_lf(q))

    def test_19_output_exclusive_owned_workstream(self):
        self.assertIn('lo_certification_benchmark',str(m.OWNED))
        self.assertNotIn('fermion_boundary/releases',str(m.OWNED))
        self.assertNotIn('joint_scalar_matching/',str(m.OWNED))

    def test_20_local_integration_preflight_or_skip(self):
        if not (ROOT/m.b4.BENCH).is_file():
            self.skipTest('Full Windows checkout not mounted in build environment')
        if not (ROOT/m.OWNED/'batch4_common_input_source01'/'summary.json').is_file():
            self.skipTest('Successful Batch4 source01 output missing; run --batch4-run-id explicitly')
        self.assertEqual(m.preflight(ROOT,'source01')['status'],
                         'PREFLIGHT_OK_NO_OUTPUT_WRITTEN')

if __name__=='__main__':
    unittest.main()
