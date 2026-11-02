"""Batch6 isolated algebra and contractual regression tests; no physics claims."""
from __future__ import annotations
import sys
from pathlib import Path
import unittest
from tempfile import TemporaryDirectory
from copy import deepcopy
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
from Numerical.lo_certification.benchmark import native_portal_batch6 as m


def synthetic():
    q=np.zeros((4,4,8,8),float)
    # Diagonal positive scalar correction to one complex S2 component.
    q[2,2,2,2]=2.;q[2,2,3,3]=2.
    return q


def encoded(x):
    return {'real':float(np.real(x)), 'imag':float(np.imag(x))}


def state():
    return {'mS1Sq':1.4e6,'mS2Sq':1.2e6,'lambdaH1':.04,
            'lambdaH2':.05,'lambdaH2Adj':.07,'lambdaS2Adj':.02,
            'lambdaT3':encoded(.12)}

class TestNativePortalBatch6(unittest.TestCase):
    def test_01_higgs_vev_correct_neutral_real_axis_and_half(self):
        q=synthetic();h,delta=m.mass_hessian_from_quartic(q,1.4e6,1.2e6,246.)
        self.assertAlmostEqual(delta[2,2],246.**2)
        self.assertAlmostEqual(delta[3,3],246.**2)
        self.assertEqual(h[0,0],1.4e6)

    def test_02_hermitian_complex_A_recovered(self):
        a=np.diag([1.3e6,1.4e6,1.5e6,1.6e6]).astype(complex)
        a[0,2]=2.7e3+142j;a[2,0]=a[0,2].conjugate()
        b=np.zeros((4,4),complex)
        h=np.block([[a.real+b.real,-a.imag-b.imag],[a.imag-b.imag,a.real-b.real]])
        perm=np.ravel(np.column_stack((np.arange(4),np.arange(4,8))))
        aa,bb=m.hessian_to_complex(h[np.ix_(perm,perm)])
        self.assertLess(np.linalg.norm(aa-a),1e-9)
        self.assertLess(np.linalg.norm(bb),1e-9)

    def test_03_holomorphic_B_recovered(self):
        a=np.diag([1e6]*4).astype(complex)
        b=np.zeros((4,4),complex);b[1,2]=230.+520j;b[2,1]=230.+520j
        h=np.block([[a.real+b.real,-a.imag-b.imag],[a.imag-b.imag,a.real-b.real]])
        perm=np.ravel(np.column_stack((np.arange(4),np.arange(4,8))))
        aa,bb=m.hessian_to_complex(h[np.ix_(perm,perm)])
        self.assertLess(np.linalg.norm(aa-a),1e-9)
        self.assertLess(np.linalg.norm(bb-b),1e-9)

    def test_04_channel_real_scalar_source(self):
        c=m.candidate_couplings(state(),'adjoint_alias')
        self.assertEqual(c['lambdaH2Inv2'],.07)
        self.assertEqual(c['lambdaH1Inv1'],.04)
        self.assertEqual(c['lambdaT3'],{'real':.12,'imag':0.})

    def test_05_alternative_alias_not_silent(self):
        c=m.candidate_couplings(state(),'ordinary_alias')
        self.assertEqual(c['lambdaH2Inv2'],.05)
        self.assertNotEqual(c['lambdaH2Inv2'],m.candidate_couplings(state(),'adjoint_alias')['lambdaH2Inv2'])

    def test_06_unknown_alias_fails(self):
        with self.assertRaises(m.Batch6Error):
            m.candidate_couplings(state(),'my_guess')

    def test_07_scalar_pure_quartic_not_portal(self):
        data=m.classify_omissions({'lambdaH2Adj':encoded(.03),'lambdaS2Adj':encoded(.01)})
        self.assertTrue(data['pure_scalar_zero_background_tree_hessian'])
        self.assertIn('heavy scalar legs',data['note'])
        self.assertEqual(data['candidate_higgs_scalar_portal']['real'],.03)

    def test_08_missing_higgs_portal_allowed(self):
        self.assertIsNone(m.classify_omissions({'lambdaS2Adj':encoded(.01)})['candidate_higgs_scalar_portal'])

    def test_09_unknown_omission_fails_closed(self):
        with self.assertRaisesRegex(m.Batch6Error,'Unresolved additional'):
            m.classify_omissions({'lambdaHHdagS2S2':encoded(.01)})

    def test_10_negative_scalar_spectrum_blocked(self):
        q=synthetic();q[2,2,2,2]=-2.;q[2,2,3,3]=-2.
        with self.assertRaisesRegex(m.Batch6Error,'Nonpositive'):
            m.mass_hessian_from_quartic(q,1e6,1.,246.)

    def test_11_bad_symmetric_quartic_fails(self):
        q=synthetic();q[1,2,0,0]=.1
        with self.assertRaisesRegex(m.Batch6Error,'Higgs leg symmetry'):
            m.mass_hessian_from_quartic(q,1e6,1e6,246.)

    def test_12_bad_scalar_exchange_fails(self):
        q=synthetic();q[2,2,0,1]=.2
        with self.assertRaisesRegex(m.Batch6Error,'Heavy scalar leg symmetry'):
            m.mass_hessian_from_quartic(q,1e6,1e6,246.)

    def test_13_scalar_stage_UV_F_S(self):
        st=state();frozen={'base_state':{'ordinary':deepcopy(st)}}
        running={'uv_at_fermion_threshold':dict(st,mS1Sq=1.6e6),
                 'eft1_at_scalar_threshold':dict(st,mS1Sq=1.8e6)}
        self.assertEqual(m.scalar_stage(frozen,running,'UV')['mS1Sq'],1.4e6)
        self.assertEqual(m.scalar_stage(frozen,running,'F')['mS1Sq'],1.6e6)
        self.assertEqual(m.scalar_stage(frozen,running,'S')['mS1Sq'],1.8e6)

    def test_14_missing_evolved_mass_fails(self):
        s=state();s.pop('mS2Sq')
        with self.assertRaisesRegex(m.Batch6Error,'incomplete scalar running state'):
            m.scalar_stage({'base_state':{'ordinary':s}}, {'uv_at_fermion_threshold':s,'eft1_at_scalar_threshold':s},'F')

    def test_15_unsupported_complex_scalar_encoding_fails(self):
        with self.assertRaisesRegex(m.Batch6Error,'Unrecognised'):
            m.value('0.5','lambda')

    def test_16_negative_mass_inputs_fails(self):
        with self.assertRaisesRegex(m.Batch6Error,'Invalid scalar mass'):
            m.mass_hessian_from_quartic(synthetic(),-100.,1e6,246.)

    def test_17_aliased_projection_never_certifies(self):
        code=Path(m.__file__).read_text(encoding='utf-8')
        self.assertIn("'native_quartic_alias_certified':False",code)
        self.assertIn("'gammaCC_to_Weyl_signed_factor':None",code)
        self.assertIn("'full_hard_fermion_one_loop_matching_included':False",code)

    def test_18_only_one_scalar_threshold_used(self):
        code=Path(m.__file__).read_text(encoding='utf-8')
        self.assertIn("'both_scalars_integrated_out_simultaneously_at_muS_in_Batch4':True",code)
        self.assertNotIn('separate_S1_threshold',code)

    def test_19_output_isolated(self):
        self.assertIn('lo_certification_benchmark',str(m.OWNED))
        self.assertNotIn('fermion_boundary/releases',str(m.OWNED))

    def test_20_optional_windows_preflight(self):
        if not (ROOT/m.b4.BENCH).is_file():
            self.skipTest('Full project checkout not mounted here')
        if not (ROOT/m.OWNED/'batch4_common_input_source01'/'summary.json').is_file():
            self.skipTest('Missing Batch4 source01 in local checkout')
        if not (ROOT/m.OWNED/'batch5_uv_eft_source01'/'summary.json').is_file():
            self.skipTest('Missing Batch5 source01 in local checkout')
        self.assertEqual(m.preflight(ROOT,'source01','source01')['status'],
                         'PREFLIGHT_OK_NO_OUTPUT_WRITTEN')


if __name__=='__main__':unittest.main()
