"""Worker C V2 tests; all pure tests execute without Wolfram or local pipeline."""
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import types
import unittest
from unittest.mock import patch

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
PATH=ROOT/'Numerical'/'joint_scalar_matching'/'worker_c_stream_v2.py'
spec=importlib.util.spec_from_file_location('worker_c_stream_v2_tested',PATH)
stream=importlib.util.module_from_spec(spec)
spec.loader.exec_module(stream)


class TestT3WorkerCV2(unittest.TestCase):
    def test_16_unique_models(self):
        self.assertEqual(len(stream.MODELS),16)
        self.assertEqual(len(set(stream.MODELS)),16)
        self.assertEqual(stream.slug('D_-2'),'D_m2')
        self.assertEqual(stream.slug('A_+0'),'A_p0')

    def test_scale_rejects_old_rough_benchmark(self):
        with self.assertRaisesRegex(stream.HandoffError,'mu_s_gev'):
            stream._scale(3000,1000,'mu_s_gev')
        with self.assertRaisesRegex(stream.HandoffError,'mu_f_gev'):
            stream._scale(4000,100000,'mu_f_gev')

    def test_complex_matrix_round_trip(self):
        z=np.array([[1,2j,3],[2j,4,5],[3,5,6]],complex)*1e-16
        np.testing.assert_allclose(stream._decode(stream._encode(z)),z)
        self.assertAlmostEqual(stream._symmetry(stream._encode(z)),0)

    def test_nonfinite_c5_fails(self):
        with self.assertRaises(stream.HandoffError):
            stream._encode(np.full((3,3),np.nan))

    def test_no_overwrite(self):
        with tempfile.TemporaryDirectory() as td:
            out=Path(td)/'result.json'
            stream._json_dump(out,{'first':True})
            with self.assertRaises(FileExistsError):
                stream._json_dump(out,{'second':True})
            self.assertEqual(stream.read_doc(out),{'first':True})

    def test_template_fail_closed(self):
        template=stream.endpoint_contract_template('D_-2')
        self.assertFalse(template['C_LLSS_evolved_to_mu_s'])
        self.assertFalse(template['no_double_counted_eft1_log'])
        self.assertEqual(template['mu_s_gev'],1000)
        self.assertEqual(template['yukawa_sector'],stream.RESTRICTED)

    def test_endpoint_contract_no_tree_promotion(self):
        with tempfile.TemporaryDirectory() as td:
            p=self._endpoint_metadata(Path(td))
            with self._mock_adapter():
                with self.assertRaisesRegex(stream.HandoffError,'C_LLSS_evolved_to_mu_s'):
                    stream.validate_endpoint_contract(p)

    @staticmethod
    def _mock_adapter():
        adapter=types.ModuleType('Numerical.joint_scalar_matching.adapter')
        adapter.CONVENTION='T3_REAL_LLSS_HHSS_WEINBERG_v1'
        adapter.normalize_model=lambda name:name
        return patch.dict(sys.modules,{'Numerical':types.ModuleType('Numerical'),
                                       'Numerical.joint_scalar_matching':types.ModuleType('Numerical.joint_scalar_matching'),
                                       'Numerical.joint_scalar_matching.adapter':adapter})

    @staticmethod
    def _endpoint_metadata(td):
        doc=stream.endpoint_contract_template('D_-2')
        doc.update(scalar_mass_origin='unbroken_electroweak',
                   all_fermions_integrated_at_mu_f=True,
                   all_scalars_joint_at_mu_s=True,
                   C5_EFT1_evolved_to_mu_s=True,
                   no_double_counted_eft1_log=True,
                   canonical_real_matching_basis_verified=True)
        p=td/'D_m2.json'
        p.write_text(json.dumps(doc))
        np.savez(td/'D_m2.npz', C5_EFT1=np.zeros((3,3)))
        return p

    def test_endpoint_requires_restricted_yukawa_provenance(self):
        with tempfile.TemporaryDirectory() as td:
            p=self._endpoint_metadata(Path(td))
            doc=stream.read_doc(p)
            doc['C_LLSS_evolved_to_mu_s']=True
            doc['yukawa_sector']='Y1_Y2_y2Conj'
            p.write_text(json.dumps(doc))
            with self._mock_adapter():
                with self.assertRaisesRegex(stream.HandoffError,'Restricted Yukawa'):
                    stream.validate_endpoint_contract(p)

    def test_endpoint_audit_does_not_claim_full_loop(self):
        with tempfile.TemporaryDirectory() as td:
            p=self._endpoint_metadata(Path(td))
            doc=stream.read_doc(p)
            doc['C_LLSS_evolved_to_mu_s']=True
            doc['spinor_normalization_verified']=True
            doc['scalar_quartic_and_mass_provenance_verified']=True
            doc['full_scalar_one_loop_operator_inventory_verified']=True
            p.write_text(json.dumps(doc))
            adapter=types.ModuleType('Numerical.joint_scalar_matching.adapter')
            adapter.CONVENTION='T3_REAL_LLSS_HHSS_WEINBERG_v1'
            adapter.normalize_model=lambda name:name
            with patch.dict(sys.modules,{'Numerical':types.ModuleType('Numerical'),
                                         'Numerical.joint_scalar_matching':types.ModuleType('Numerical.joint_scalar_matching'),
                                         'Numerical.joint_scalar_matching.adapter':adapter}):
                result=stream.validate_endpoint_contract(p)
            self.assertTrue(result['all_physical_checks_passed'])
            self.assertFalse(result['safe_to_claim_full_one_loop'])

    def test_survey_respects_missing_model_and_writes_blocked_reports(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            out=root/'fresh'
            report=stream.survey_all(root,root/'no_source',root/'no_background',
                                     models=('D_-2','A_-4'),output_dir=out)
            self.assertEqual(report['blocked'],2)
            self.assertFalse(report['absolute_physical_C5_available'])
            self.assertTrue((out/'D_m2.json').is_file())
            self.assertTrue((out/'A_m4.json').is_file())
            with self.assertRaisesRegex(stream.HandoffError,'already exists'):
                stream.survey_all(root,root,root,models=('D_-2',),output_dir=out)

    def test_survey_mocked_end_to_end_decomposition(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td); src=root/'source'; bg=root/'background'
            (src/'D_m2').mkdir(parents=True);bg.mkdir()
            (src/'D_m2'/'boundary.json').write_text('{}')
            (bg/'D_m2.json').write_text('{}')
            (bg/'D_m2.npz').write_bytes(b'synthetic-test-bytes')
            (bg/'D_m2.json').write_text(json.dumps({'tensor_npz':'D_m2.npz'}))
            z=np.eye(3)*1e-16
            dataset={'model':'D_-2',
                     'C5_joint_scalar_finite_bubble_per_N_gev_inverse':stream._encode(z),
                     'C5_EFT1_induced_leading_log_per_N_gev_inverse':stream._encode(-z*2),
                     'C5_LLSS_total_one_loop_per_N_gev_inverse':stream._encode(-z),
                     'relative_matching_scale_cancellation_residual':0.0}
            source=types.SimpleNamespace(mu_f_gev=1e5)
            background=types.SimpleNamespace(mu_f=1e5,mu_s=1e3,reference_scale=1e3)
            bridge=types.ModuleType('Numerical.joint_scalar_matching.fermion_boundary_bridge')
            bridge.load_source=lambda path,expected_model:source
            large=types.ModuleType('Numerical.joint_scalar_matching.large_batch')
            large.load_background=lambda path,source:background
            large.evaluate_conditional=lambda source,background:dataset
            with patch.dict(sys.modules,{
                'Numerical':types.ModuleType('Numerical'),
                'Numerical.joint_scalar_matching':types.ModuleType('Numerical.joint_scalar_matching'),
                'Numerical.joint_scalar_matching.fermion_boundary_bridge':bridge,
                'Numerical.joint_scalar_matching.large_batch':large}):
                result=stream.survey_one(src/'D_m2'/'boundary.json',bg/'D_m2.json','D_-2')
            self.assertAlmostEqual(result['diagnostics']['bubble_norm_gev_inv_per_unit_N'],np.linalg.norm(z))
            self.assertFalse(result['restricted_Y1_Y2_source_certified'])

    def test_source_fail_bad_scale(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            p=root/'source.json';p.write_text('{}')
            q=root/'background.json';q.write_text(json.dumps({'tensor_npz':'b.npz'}))
            (root/'b.npz').write_bytes(b'x')
            bridge=types.ModuleType('Numerical.joint_scalar_matching.fermion_boundary_bridge')
            bridge.load_source=lambda path,expected_model:types.SimpleNamespace(mu_f_gev=1e5)
            large=types.ModuleType('Numerical.joint_scalar_matching.large_batch')
            large.load_background=lambda path,source:types.SimpleNamespace(mu_f=1e5,mu_s=3000,reference_scale=3000)
            large.evaluate_conditional=lambda *a: (_ for _ in ()).throw(AssertionError('must not call'))
            with patch.dict(sys.modules,{
                'Numerical':types.ModuleType('Numerical'),
                'Numerical.joint_scalar_matching':types.ModuleType('Numerical.joint_scalar_matching'),
                'Numerical.joint_scalar_matching.fermion_boundary_bridge':bridge,
                'Numerical.joint_scalar_matching.large_batch':large}):
                with self.assertRaisesRegex(stream.HandoffError,'background.mu_s_gev'):
                    stream.survey_one(p,q,'D_-2')

    def test_remaining_work_is_explicit(self):
        r=stream.remaining_work()
        self.assertGreaterEqual(len(r),6)
        self.assertTrue(any('additional independent finite scalar' in row['item'].lower() for row in r))
        self.assertTrue(any(row['owner']=='Worker D' for row in r))


if __name__=='__main__':
    unittest.main()

class TestEndToEndMocked(unittest.TestCase):
    def test_endpoint_inherited_once_channels_and_pairs(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            doc=stream.endpoint_contract_template('D_-2')
            doc.update(scalar_mass_origin='unbroken_electroweak',
                       all_fermions_integrated_at_mu_f=True,
                       all_scalars_joint_at_mu_s=True,
                       C_LLSS_evolved_to_mu_s=True,
                       C5_EFT1_evolved_to_mu_s=True,
                       no_double_counted_eft1_log=True,
                       canonical_real_matching_basis_verified=True)
            p=root/'D_m2.json';p.write_text(json.dumps(doc))
            (root/'D_m2.npz').write_bytes(b'fixture-mocked-only')
            z=np.eye(3,dtype=complex)*1e-18
            bubble=np.eye(3,dtype=complex)*2e-18
            adapter=types.ModuleType('Numerical.joint_scalar_matching.adapter')
            adapter.CONVENTION='T3_REAL_LLSS_HHSS_WEINBERG_v1'
            adapter.normalize_model=lambda name:name
            adapter.scalar_groups=lambda model:['S1']*6+['S2']*2
            adapter.load_input=lambda path:types.SimpleNamespace(
                model='D_-2',c5_eft1=z,c_llss=np.zeros((3,3,2,2,8,8)),
                mass_squared=np.eye(8),lambda_hhss=np.zeros((4,4,8,8)),
                mu_s_gev=1000,scalar_mass_origin='unbroken_electroweak',
                provenance={'source':'mock evolving EFT1'},
                quartic_channels=None,quartic_channel_names=())
            worker=types.ModuleType('Numerical.joint_scalar_matching.worker_c_joint')
            worker.match_joint_scalars=lambda ew,sb,mu:{
                'c5_inherited':z,'c5_finite_scalar':bubble,'c5_total':z+bubble,
                'weinberg_projection_residual':0.0,
                'scalar_eigenvalues_sq_gev2':np.ones(8),
                'contribution_breakdown':{'eigenstate_pairs':{'e00_e01':bubble},
                                          'quartic_channels':{}}}
            with patch.dict(sys.modules,{
                'Numerical':types.ModuleType('Numerical'),
                'Numerical.joint_scalar_matching':types.ModuleType('Numerical.joint_scalar_matching'),
                'Numerical.joint_scalar_matching.adapter':adapter,
                'Numerical.joint_scalar_matching.worker_c_joint':worker}):
                output=stream.run_endpoint(p)
            self.assertEqual(output['status'],'EvolvedEFT1_OneJointScalarBubble_Evaluated')
            self.assertAlmostEqual(output['c5_norms_gev_inv']['total'],np.linalg.norm(z+bubble))
            self.assertFalse(output['complete_scalar_matching_certified'])
            self.assertFalse(output['safe_to_add_as_second_scalar_bubble'])

    def test_b_c_mixed_basis_explicit(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            doc=stream.endpoint_contract_template('B_-1')
            doc.update(scalar_mass_origin='unbroken_electroweak',
                       all_fermions_integrated_at_mu_f=True,
                       all_scalars_joint_at_mu_s=True,
                       C_LLSS_evolved_to_mu_s=True,
                       C5_EFT1_evolved_to_mu_s=True,
                       no_double_counted_eft1_log=True,
                       canonical_real_matching_basis_verified=True)
            p=root/'B_m1.json';p.write_text(json.dumps(doc))
            (root/'B_m1.npz').write_bytes(b'mock')
            adapter=types.ModuleType('Numerical.joint_scalar_matching.adapter')
            adapter.CONVENTION='T3_REAL_LLSS_HHSS_WEINBERG_v1'
            adapter.normalize_model=lambda name:name
            with patch.dict(sys.modules,{'Numerical':types.ModuleType('Numerical'),
                         'Numerical.joint_scalar_matching':types.ModuleType('Numerical.joint_scalar_matching'),
                         'Numerical.joint_scalar_matching.adapter':adapter}):
                with self.assertRaisesRegex(stream.HandoffError,'mixed-scalar'):
                    stream.validate_endpoint_contract(p)

class TestArchiveRegression(unittest.TestCase):
    def test_same_input_archive_matches(self):
        with tempfile.TemporaryDirectory() as td:
            z=stream._encode(np.eye(3)*1e-16)
            fields={key:z for key in ('C5_joint_scalar_finite_bubble_per_N_gev_inverse',
                       'C5_EFT1_induced_leading_log_per_N_gev_inverse',
                       'C5_LLSS_total_one_loop_per_N_gev_inverse')}
            data={'model':'D_-2','model_source_boundary_sha256':'a',
                  'background_json_sha256':'b','background_npz_sha256':'c',
                  'mu_f_gev':100000.0,'mu_s_gev':1000.0,**fields}
            p=Path(td)/'archive.json'
            p.write_text(json.dumps({'model':'D_-2','response':data}))
            current={'model':'D_-2','conditional_response_per_unit_unknown_normalization':data}
            result=stream.compare_archived_conditional(current,p)
            self.assertEqual(result['status'],'IdenticalInputsNumericalRegressionPassed')
            self.assertEqual(result['max_relative_matrix_difference'],0.)

    def test_changed_input_cannot_be_numerically_compared(self):
        with tempfile.TemporaryDirectory() as td:
            old={'model_source_boundary_sha256':'old','mu_s_gev':1000.0}
            now={'model_source_boundary_sha256':'new','mu_s_gev':3000.0}
            p=Path(td)/'archive.json'
            p.write_text(json.dumps({'model':'D_-2','response':old}))
            result=stream.compare_archived_conditional(
                {'model':'D_-2','conditional_response_per_unit_unknown_normalization':now},p)
            self.assertEqual(result['status'],'IncomparableDifferentInputs')
            self.assertIn('mu_s_gev',result['mismatched_fields'])
