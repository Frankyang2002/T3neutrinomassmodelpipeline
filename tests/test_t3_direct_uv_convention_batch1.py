"""Isolated numerical tests for the additive integrated UV convention audit."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import numpy as np

SOURCE = (Path(__file__).resolve().parents[1] / 'scripts' /
          'AuditT3DirectUVConventionBatch1.py')
spec = importlib.util.spec_from_file_location('direct_uv_convention_batch1', SOURCE)
assert spec is not None and spec.loader is not None
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


class MassConventions(unittest.TestCase):
    def setUp(self):
        self.mf = np.array([[10+2j,1-3j,2],[3,12-1j,-4j],[0.5,1j,17]],complex)

    def test_dirac_doubled_weyl_mass_is_symmetric_for_general_flavour_mf(self):
        m = mod.mass_tensor(self.mf,2,False,np.eye(2))
        self.assertEqual(m.shape,(12,12))
        self.assertTrue(np.allclose(m,m.T,atol=1e-15))
        self.assertTrue(np.allclose(m[:6,6:],np.kron(self.mf,np.eye(2))))
        self.assertTrue(np.allclose(m[6:,:6],np.kron(self.mf.T,np.eye(2))))
        self.assertTrue(np.allclose(np.diag(m),0))

    def test_majorana_charge_pairing_metric_nontrivial(self):
        mf=(self.mf+self.mf.T)/2
        pair=np.array([[0,0,1],[0,-1,0],[1,0,0]],complex)
        m=mod.mass_tensor(mf,3,True,pair)
        self.assertEqual(m.shape,(9,9))
        self.assertTrue(np.allclose(m,m.T))
        self.assertAlmostEqual(m[0,2].real,mf[0,0].real)
        self.assertAlmostEqual(m[1,1].real,-mf[0,0].real)
        with self.assertRaisesRegex(ValueError,'not symmetric'):
            mod.mass_tensor(self.mf,3,True,pair)

    def test_negative_majorana_antisymmetric_pairing_rejected(self):
        pair=np.array([[0,1],[-1,0]],complex)
        with self.assertRaisesRegex(ValueError,'Nonsymmetric'):
            mod.mass_tensor(np.eye(3),2,True,pair)


class ScalarConventions(unittest.TestCase):
    def setUp(self):
        self.mix=np.array([[0.6+0.8j,-0.3],[0.2j,1.2-0.5j]])
        self.v=246.0
        self.o={'mS1Sq':800**2,'mS2Sq':920**2,'lambdaH1':0.13,
                'lambdaH2':-0.07,'lambdaT3':{'real':0.002,'imag':0.001}}
        self.a=mod.scalar_tensor(self.o['mS1Sq'],self.o['mS2Sq'],
                    self.o['lambdaH1'],self.o['lambdaH2'],
                    mod._number(self.o['lambdaT3']),self.v,self.mix)

    def test_complex_scalar_action_and_mixing_orientation(self):
        b=np.zeros_like(self.a)
        r=mod.scalar_action_witness(self.a,b,2,self.o,self.v,self.mix)
        self.assertEqual(r['independent_complex_field_witnesses'],2)
        self.assertLess(r['maximum_absolute_source_potential_discrepancy_gev4'],1e-8)
        self.assertTrue(np.allclose(self.a,self.a.conj().T))
        value=-(mod._number(self.o['lambdaT3'])*self.v**2/2)*self.mix[0,1]
        self.assertAlmostEqual(abs(self.a[3,0]-value),0,delta=1e-10)

    def test_sign_half_and_portal_negative_controls(self):
        checks=mod.candidate_negative_controls(
            self.a,self.mix,self.v,mod._number(self.o['lambdaT3']),
            self.o['lambdaH1'],self.o['lambdaH2'],2)
        self.assertEqual(checks['nonzero_controls_checked'],3)
        self.assertTrue(all(checks['controls'].values()))
        wrong=self.a.copy()
        wrong[2:,:2]*=-1
        wrong[:2,2:]=wrong[2:,:2].conj().T
        with self.assertRaisesRegex(ValueError,'not represented'):
            mod.scalar_action_witness(wrong,np.zeros_like(wrong),
                                      2,self.o,self.v,self.mix)

    def test_negative_controls_are_standard_json_booleans(self):
        checks=mod.candidate_negative_controls(
            self.a,self.mix,self.v,mod._number(self.o['lambdaT3']),
            self.o['lambdaH1'],self.o['lambdaH2'],2)
        self.assertEqual(checks['nonzero_controls_checked'],3)
        self.assertTrue(all(type(flag) is bool for flag in checks['controls'].values()))
        # Regression for the original crash: numpy.bool_ cannot be encoded by json.
        self.assertEqual(json.loads(json.dumps(checks)),checks)

    def test_hermiticity_guard_rejects_one_sided_mixing(self):
        bad=self.a.copy()
        bad[0,3]+=1j
        with self.assertRaises(ValueError):
            mod.scalar_action_witness(bad,np.zeros_like(bad),
                                      2,self.o,self.v,self.mix)

    def test_scalar_vev_normalization_doublet_and_singlet(self):
        one=np.array([[1+0j],[0.4j]])
        a=mod.scalar_tensor(100,200,0.2,-0.4,1.5,4,one)
        self.assertAlmostEqual(a[0,0].real,100-0.2*16/2)
        self.assertAlmostEqual(a[2,2].real,200+0.4*16/2)
        self.assertAlmostEqual(a[2,0].real,-1.5*16/2)
        self.assertAlmostEqual(a[2,1].imag,-1.5*16/2*0.4)


class GuardConventions(unittest.TestCase):
    def test_missing_prior_evidence_fails_closed(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            with self.assertRaisesRegex(ValueError,'Missing required evidence'):
                mod._previous_evidence(root,root/'source_audits')

    def test_charge_support_and_wrong_channel_rejected(self):
        c1=np.zeros((2,1,2),complex); c1[0,0,0]=1
        c2=np.zeros((2,1,2),complex); c2[0,0,1]=1
        mix=np.zeros((2,2),complex);mix[0,1]=1
        got=mod.charge_channel_witness(2,2,1,-1,np.eye(1),c1,c2,mix,True)
        self.assertEqual(got['native_component_charge_channels_checked'],4)
        wrong=c1.copy();wrong[0,0,1]=1
        with self.assertRaisesRegex(ValueError,'electric charge'):
            mod.charge_channel_witness(2,2,1,-1,np.eye(1),wrong,c2,mix,True)

    def test_16_unique_model_labels(self):
        pairs=list(mod.models())
        self.assertEqual(len(pairs),16)
        self.assertEqual(len({a for a,_ in pairs}),16)
        self.assertEqual(len({b for _,b in pairs}),16)

    def test_complete_audit_report_is_json_serializable(self):
        # Synthetic source witnesses exercise the full report assembly without
        # requiring a Wolfram installation or all-16 local native CG exports.
        controls=mod.candidate_negative_controls(
            np.eye(2,dtype=complex),np.ones((1,1),complex),246,
            0.003+0.001j,0.2,-0.1,1)
        def fake_case(root,short,key,cases,rot,physical,audits):
            return {'model':short, 'status':mod.CONDITIONAL,
                    'negative_controls':controls,
                    'physical_uv_mass_certified':False}
        with (mock.patch.object(mod,'check_source_contracts',return_value={'fixture':'sha'}),
              mock.patch.object(mod,'_previous_evidence',return_value={
                  'cases':{},'rotations':{},'source_evidence_hashes':{}}),
              mock.patch.object(mod,'prior_spinor_evidence',return_value={
                  'present':False,'cross_package_physical_dictionary_certified':False}),
              mock.patch.object(mod,'audit_case',side_effect=fake_case)):
            report=mod.audit_all(Path('.'),Path('.'),Path('.'))
        self.assertEqual(report['algebraic_candidate_models_passed'],16)
        self.assertEqual(report['physical_models_certified'],0)
        self.assertEqual(len(json.loads(json.dumps(report,allow_nan=False))['models']),16)

    def test_zerocoupling_negative_control_is_not_counted_as_pass(self):
        a=np.eye(2,dtype=complex)
        report=mod.candidate_negative_controls(a,np.ones((1,1)),246,0j,0,0,1)
        self.assertEqual(report['nonzero_controls_checked'],0)


if __name__=='__main__':
    unittest.main()
