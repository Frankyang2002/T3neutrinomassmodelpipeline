"""Isolated algebra, provenance, mass-basis and fail-closed tests for next batch."""
from __future__ import annotations

import json
from math import log, pi
from pathlib import Path
from tempfile import TemporaryDirectory
from types import ModuleType, SimpleNamespace
from unittest import TestCase
from unittest.mock import patch
import sys

import numpy as np

from RGE.running.intermediate.All16PreliminaryNumerics import (
    FORMAT, MIXED_MODELS, SPLIT_MODELS, MissingCertifiedBasis,
    _ordinary_scalar_point, _mixed_scalar_point, _mass_c5, output_summary,
    resolve_scalar_point, tag,
)
from RGE.running.intermediate.B1PreliminaryNumerics import real_lift_b1


class AlgebraTests(TestCase):
    def test_model_tags(self):
        self.assertEqual(tag('A_-4'),'A_m4')
        self.assertEqual(tag('E_+2'),'E_p2')
        self.assertEqual(tag('C_-1'),'C_m1')

    def test_exact_model_partition(self):
        models={'A_-4','A_-2','A_+0','B_-3','B_-1','B_+1','C_-3',
                'C_-1','C_+1','D_-2','D_+0','D_+2','E_-4','E_-2','E_+0','E_+2'}
        self.assertEqual(len(models),16)
        self.assertEqual(len(SPLIT_MODELS),6)
        self.assertEqual(len(MIXED_MODELS),2)
        self.assertTrue(SPLIT_MODELS.isdisjoint(MIXED_MODELS))
        self.assertEqual(len(models-SPLIT_MODELS-MIXED_MODELS),8)

    def test_split_real_field_cannot_be_silently_identified(self):
        for model in SPLIT_MODELS:
            with self.subTest(model=model), self.assertRaisesRegex(MissingCertifiedBasis,'not|No independent'):
                resolve_scalar_point(Path('.'),model,None,None)

    def test_real_lift_complex_eps(self):
        alpha=.113
        U=np.array([[np.cos(alpha),np.sin(alpha)*np.exp(1j*.2)],
                    [-np.sin(alpha)*np.exp(-1j*.2),np.cos(alpha)]])
        R=real_lift_b1(U)
        self.assertEqual(R.shape,(8,8))
        self.assertLess(np.linalg.norm(R.T@R-np.eye(8)),1e-12)

    def test_muS_ordinary_illustrative_provenance_and_all_12_components(self):
        with TemporaryDirectory() as temp:
            root=Path(temp)
            p=root/'output/general_observable/all16_v1/E_m4/running.json'
            p.parent.mkdir(parents=True)
            j={'model':'E_-4','status':'NumericalRenormalisableRunningOnly',
               'illustrative_benchmark':True,'general_matching_performed':False,
               'physical_mass_parameters_available':False,
               'eft1_wilson_running_performed':False,
               'scales_gev':{'F':100000,'S':1000},
               'eft1_at_scalar_threshold':{'mS1Sq':1.7e6,'mS2Sq':2.3e6,
                    'lambdaT3':{'real':.012,'imag':-.003}}}
            p.write_text(json.dumps(j))
            source=SimpleNamespace(mu_f_gev=100000.,labels=('NewScalar1',)*6+('NewScalar2',)*6)
            ch=SimpleNamespace(source_basis_labels=source.labels)
            r=_ordinary_scalar_point(root,'E_-4',source,ch)
            self.assertEqual(r.mass2_boundary_muS.shape,(12,12))
            self.assertTrue(np.array_equal(np.diag(r.mass2_boundary_muS),[1.7e6]*6+[2.3e6]*6))
            self.assertEqual(r.lambda_candidate,.012-.003j)
            self.assertFalse(r.inputs['muF_matched_to_this_running_benchmark'])
            j['physical_mass_parameters_available']=True
            p.write_text(json.dumps(j))
            with self.assertRaises(ValueError):_ordinary_scalar_point(root,'E_-4',source,ch)

    def test_ordinary_rejects_special_basis_even_with_mass(self):
        with TemporaryDirectory() as temp:
            root=Path(temp);path=root/'output/general_observable/all16_v1/A_m4/running.json'
            path.parent.mkdir(parents=True)
            path.write_text(json.dumps({'model':'A_-4','status':'NumericalRenormalisableRunningOnly',
                'illustrative_benchmark':True,'general_matching_performed':False,
                'physical_mass_parameters_available':False,'eft1_wilson_running_performed':False,
                'scales_gev':{'F':1e5,'S':1e3},'eft1_at_scalar_threshold':{'mS1Sq':1e6,'mS2Sq':2e6,'lambdaT3':.01}}))
            source=SimpleNamespace(mu_f_gev=100000.,labels=('NewScalarX',)*8)
            ch=SimpleNamespace(source_basis_labels=('NewScalar1',)*8)
            with self.assertRaisesRegex(ValueError,'DIFFER'):
                _ordinary_scalar_point(root,'A_-4',source,ch)

    def test_mixed_basis_source_handoff_accepts_both_slash_styles(self):
        with TemporaryDirectory() as temp:
            root=Path(temp)
            name='C_m1';model='C_-1'
            f=root/f'output/general_observable/mixed_f_synchronised_v1/{name}/mass_basis/{name}_final_mass_basis.json'
            s=root/f'output/general_t3_mass_basis/{name}_final_mass_basis.json'
            o=root/f'output/general_observable/ordered_wilson_tensor_v4/{name}/ordered_tensor.json'
            fp=root/f'output/general_observable/mixed_f_synchronised_v1/{name}/matchete/{name}_final_parameters.json'
            for p in (f,s,o,fp):p.parent.mkdir(parents=True,exist_ok=True)
            eigen=[1.6e6,2.5e6]
            U={'real':[[1.,0.],[0.,1.]],'imag':[[0.,0.],[0.,0.]]}
            A={'real':[[eigen[0],0.],[0.,eigen[1]]],'imag':[[0.,0.],[0.,0.]]}
            def basis(mu):
                return {'format':'t3-mixed-doublet-mass-basis-transport-v1','model':model,
                        'scale_gev':mu,'status':'ExperimentalMassBasisTensorTransportNoPhysicalMatching',
                        'quadratic_hessian_agreement':True,'max_su2_covariance_absolute_error':0.,
                        'basis':{'basis':'(S1, epsilon*conjugate(S2)) = U @ (X,Y)',
                        'epsilon':[[0.,1.],[-1.,0.]],'unitarity_verified':True,
                        'diagonal_mass_verified':True,'unitary_flavour_rotation_U':U,
                        'mass_squared_flavour_matrix_gev2':A,'masses_squared_gev2':eigen},
                        'quartics':{'lambdaT3':{'source_rgbeta_coupling':{'real':.01,'imag':.002}}}}
            f.write_text(json.dumps(basis(100000)))
            s.write_text(json.dumps(basis(1000)))
            o.write_text(json.dumps({'format':'t3-mixed-native-matchete-ordered-tensor-v1',
                'model':model,'mu_F_gev':100000,'source_mass_basis_sha256':'mock-sha'}))
            fp.write_text(json.dumps({'format':'t3-experimental-mass-eigen-matchete-parameters-v1',
                'model':model,'source_eft1_scale_gev':100000,'source_sha256':'mock-sha',
                'source_path':str(f.relative_to(root)).replace('/','\\'),'unitary_U':U}))
            source=SimpleNamespace(mu_f_gev=100000.,labels=('NewScalarX',)*4+('NewScalarY',)*4)
            ch=SimpleNamespace(n=8)
            point=_mixed_scalar_point(root,model,source,ch)
            self.assertEqual(point.status,'SourceSyncedMixedDoubletCandidate_NotPhysical')
            self.assertEqual(point.mass2_boundary_muS.shape,(8,8))
            self.assertLess(np.max(abs(np.sort(np.linalg.eigvalsh(point.mass2_boundary_muS))
                      -np.sort(np.repeat(eigen,4)))),1e-7)

    def test_complete_summary_explicitly_reports_nonphysical(self):
        rows=[{'model':'A_-4','status':'ConditionalWeinbergLambdaT3OneChannel'},
              {'model':'A_-2','status':'BlockedMissingCertifiedScalarBasis'}]
        s=output_summary(rows)
        self.assertEqual(s['models_with_Ward_valid_one_channel_results'],1)
        self.assertEqual(s['models_blocked_special_scalar_basis'],1)
        self.assertFalse(s['physical_neutrino_masses_computed'])
        self.assertEqual(s['source_to_Weyl_factor'],None)

    def test_fixed_order_one_loop_scale_cancellation_and_neutrino_mass(self):
        # Deterministic mock tests implementation-level order counting and
        # master/joint factor two; the actual repo core is used by the CLI.
        mod=ModuleType('Numerical.joint_scalar_matching.core')
        class Box:pass
        def check(C,q,n):
            self.assertEqual(C.shape,(3,3,2,2,n,n))
            self.assertEqual(q.shape,(4,4,n,n))
            return C,q
        def project(amp,strict=False):
            a=np.asarray(amp)
            return a[:,:,0,0,2,2],0.0
        def scalar(C,q,mass,mu,reference_scale=None):
            b=np.einsum('fgijab,ABab->fgijAB',C,q)
            return {'amplitude':-log(mu/reference_scale)*b/(16*pi*pi)}
        mod.check_tensors=check;mod.project_weinberg=project;mod.scalar_bubble_matching=scalar
        C=np.zeros((3,3,2,2,12,12),complex)
        q=np.zeros((4,4,12,12),float)
        C[0,0,0,0,0,0]=1e-11+2e-12j
        q[2,2,0,0]=.3
        with patch.dict(sys.modules,{'Numerical.joint_scalar_matching.core':mod}):
            result=_mass_c5(C,q,np.diag([1e6]*6+[2e6]*6),100000.,1000.)
        self.assertTrue(result['full_SU2_Weinberg_projection_pass'])
        self.assertLess(result['matching_scale_cancellation_relative'],1e-10)
        self.assertGreater(result['neutrino_mass_singular_values_ev_per_abs_NGamma'][2],0.)
        expected=2*.3*(1e-11+2e-12j)
        got=result['beta_C5_master_16pi2_per_N_gev_inverse']
        self.assertAlmostEqual(got['real'][0][0]/expected.real,1.0,places=12)

    def test_non_Ward_blocks_neutrino_mass(self):
        mod=ModuleType('Numerical.joint_scalar_matching.core')
        mod.check_tensors=lambda C,q,n:(C,q)
        mod.project_weinberg=lambda a,strict=False:(np.zeros((3,3),complex),.25)
        mod.scalar_bubble_matching=lambda C,q,m,mu,reference_scale=None:{'amplitude':
            np.zeros((3,3,2,2,4,4),complex)}
        with patch.dict(sys.modules,{'Numerical.joint_scalar_matching.core':mod}):
            out=_mass_c5(np.zeros((3,3,2,2,8,8)),np.zeros((4,4,8,8)),
                        np.eye(8)*1e6,1e5,1e3)
        self.assertEqual(out['status'],'NonWeinbergTensorBlockedMassInterpretation')
        self.assertIsNone(out['neutrino_mass_singular_values_ev_per_abs_NGamma'])

if __name__=='__main__':
    import unittest
    unittest.main()


class CLITests(TestCase):
    @staticmethod
    def _load_cli():
        import importlib.util
        from pathlib import Path
        path=Path(__file__).resolve().parents[1]/'scripts/RunAll16PreliminaryNumerics.py'
        fake=ModuleType('Numerical.joint_scalar_matching.adapter')
        fake.ALL_MODELS=('A_-4','A_-2')
        fake.normalize_model=lambda x:x if x in fake.ALL_MODELS else (_ for _ in ()).throw(ValueError('unknown model'))
        with patch.dict(sys.modules,{'Numerical.joint_scalar_matching.adapter':fake}):
            sp=importlib.util.spec_from_file_location('all16_cli_test_loader',path)
            module=importlib.util.module_from_spec(sp)
            sp.loader.exec_module(module)
        return module

    def test_preflight_check_only_and_exclusive_output(self):
        cli=self._load_cli()
        with TemporaryDirectory() as temp:
            root=Path(temp)
            cli.ROOT=root
            cli.SAFE=root/'output/general_observable/lo_eft1_mixing'
            destination=cli.SAFE/'test_output'
            def one_model(_,model):
                if model=='A_-2': raise MissingCertifiedBasis('split basis not sourced')
                return {'model':model,'status':'ConditionalWeinbergLambdaT3OneChannel',
                    'response':{'full_SU2_Weinberg_projection_pass':True,
                    'beta_SU2_residual_relative':1e-15,
                    'matching_scale_cancellation_relative':1e-14,
                    'neutrino_mass_singular_values_ev_per_abs_NGamma':[.01,.02,.03]}}
            with patch.object(cli,'evaluate_model',side_effect=one_model):
                self.assertEqual(cli.main(['--output-root',str(destination),'--check-only']),0)
                self.assertFalse(destination.exists())
                self.assertEqual(cli.main(['--output-root',str(destination)]),0)
                self.assertTrue((destination/'summary.json').is_file())
                self.assertTrue((destination/'A_m4/lambdaT3_one_channel.json').is_file())
                with patch.object(cli.argparse.ArgumentParser,'error',side_effect=ValueError('refuse')):
                    with self.assertRaisesRegex(ValueError,'refuse'):
                        cli.main(['--output-root',str(destination)])

    def test_invalid_inputs_are_nonzero_but_known_split_basis_not_fatal(self):
        cli=self._load_cli()
        with TemporaryDirectory() as temp:
            root=Path(temp);cli.ROOT=root
            cli.SAFE=root/'output/general_observable/lo_eft1_mixing'
            dest=cli.SAFE/'test_error'
            with patch.object(cli,'evaluate_model',side_effect=ValueError('bad tensor')):
                self.assertEqual(cli.main(['--output-root',str(dest)]),2)
            s=json.loads((dest/'summary.json').read_text())
            self.assertEqual(s['models_input_errors'],2)
