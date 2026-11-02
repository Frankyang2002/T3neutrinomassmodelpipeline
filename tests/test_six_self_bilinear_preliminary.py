"""Real-component and failure-path checks for all six source-split scalar cases."""
from __future__ import annotations

import json
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace, ModuleType
from unittest import TestCase
from unittest.mock import patch
import importlib.util
import sys

import numpy as np

from RGE.running.intermediate.SixSelfBilinearPreliminary import (
    SPLIT_SHAPES, native_cartesian_lift, native_original_mass2, split_scalar_point,
    evaluate_split_model,
)
from RGE.running.intermediate.All16PreliminaryNumerics import SPLIT_MODELS


def write_fixture(root:Path,model:str,*,muF:complex=15000.,muS:complex=1200.,
                  m1:float=1.8e6,m2:float=2.9e6):
    d1,d2,selected=SPLIT_SHAPES[model]
    tag=model.replace('-','m').replace('+','p')
    target=root/'output/general_observable/all16_v1'/tag/'running.json'
    target.parent.mkdir(parents=True,exist_ok=True)
    key='muT3GeneralS1Sq' if selected==1 else 'muT3GeneralS2Sq'
    cj=lambda z:{'real':complex(z).real,'imag':complex(z).imag}
    data={'model':model,'status':'NumericalRenormalisableRunningOnly',
        'illustrative_benchmark':True,'general_matching_performed':False,
        'physical_mass_parameters_available':False,'eft1_wilson_running_performed':False,
        'scales_gev':{'UV':1e7,'F':1e5,'S':1e3},
        'eft1_at_fermion_threshold':{key:cj(muF)},
        'eft1_at_scalar_threshold':{key:cj(muS),'mS1Sq':m1,'mS2Sq':m2,
                                    'lambdaT3':cj(.014)}}
    target.write_text(json.dumps(data))
    met=root/'output/general_observable/fermion_first_self_fix_v2'/tag/'self_mass_basis.json'
    met.parent.mkdir(parents=True,exist_ok=True)
    meta={'mass_basis_kind':'RealSelfBilinearSplitLightFields',
          'split_species':selected,'complex_self_bilinear_phase_supported':False,
          'field_basis_normalization_physically_certified':False,
          'mass_relation':'Mminus^2=MS^2-2 Re(kappa); Mplus^2=MS^2+2 Re(kappa)',
          'minus_field':f'T3Scalar{selected}Minus','plus_field':f'T3Scalar{selected}Plus'}
    met.write_text(json.dumps(meta))
    fst=(('NewScalar1',)*(2*d1) if selected==2 else
         ('T3Scalar1Minus',)*d1+('T3Scalar1Plus',)*d1)
    snd=(('NewScalar2',)*(2*d2) if selected==1 else
         ('T3Scalar2Minus',)*d2+('T3Scalar2Plus',)*d2)
    source=SimpleNamespace(mu_f_gev=1e5,labels=fst+snd)
    ch=SimpleNamespace(source_basis_labels=('NewScalar1',)*(2*d1)+('NewScalar2',)*(2*d2))
    return source,ch,target


class AlgebraTests(TestCase):
    def test_six_exact_shapes_and_real_orthogonality(self):
        self.assertEqual(set(SPLIT_SHAPES),set(SPLIT_MODELS))
        for model,(d1,d2,split) in SPLIT_SHAPES.items():
            with self.subTest(model=model):
                R=native_cartesian_lift(d1,d2,split)
                self.assertEqual(R.shape,(2*(d1+d2),2*(d1+d2)))
                self.assertLess(np.linalg.norm(R.T@R-np.eye(len(R))),1e-13)
                # A source real/imaginary scalar field has unit kinetic norm;
                # no erroneous 1/sqrt(2) in the *real* change of coordinates.
                self.assertTrue(np.allclose(np.sum(R*R,axis=0),1.))

    def test_native_cartesian_rotation_maintains_complex_scalar(self):
        rng=np.random.default_rng(5)
        for model,(d1,d2,split) in SPLIT_SHAPES.items():
            for phase in (1.,np.exp(.35j)):
                R=native_cartesian_lift(d1,d2,split,phase)
                real_source=rng.normal(size=2*(d1+d2))
                x=R@real_source
                offset=0 if split==1 else 2*d1
                d=(d1,d2)[split-1]
                for k in range(d):
                    z=(x[offset+2*k]+1j*x[offset+2*k+1])/np.sqrt(2)
                    alt=phase*(real_source[offset+d+k]+1j*real_source[offset+k])/np.sqrt(2)
                    self.assertLess(abs(z-alt),1e-13)

    def test_cartesian_native_vs_spherical_basis_not_conflated(self):
        # Native Matchete replacement preserves its SU(2) index. A spherical
        # real intertwiner would be an EXTRA, unsupported transformation.
        R=native_cartesian_lift(3,3,1)
        for k in range(3):
            self.assertEqual(R[2*k,3+k],1.)
            self.assertEqual(R[2*k+1,k],1.)

    def test_mass_splitting_and_complex_hessian_all_six(self):
        for model,(d1,d2,split) in SPLIT_SHAPES.items():
            for mu in (14000.,-11083.,2000.+3900j):
                with self.subTest(model=model,mu=mu):
                    m1,m2=2.1e6,2.8e6
                    mass=native_original_mass2(d1,d2,split,m1,m2,mu)
                    spectrum=np.linalg.eigvalsh(mass)
                    m=(m1,m2)[split-1];d=(d1,d2)[split-1]
                    spectator=(m2 if split==1 else m1)
                    self.assertTrue(np.allclose(np.sort(spectrum),sorted([m-2*abs(mu)]*d+
                        [m+2*abs(mu)]*d+[spectator]*2*((d2 if split==1 else d1)))))
                    R=native_cartesian_lift(d1,d2,split)
                    self.assertTrue(np.allclose(np.linalg.eigvalsh(R.T@mass@R),spectrum))

    def test_negative_muS_swaps_mass_order_without_swapping_F_fields(self):
        with TemporaryDirectory() as temp:
            root=Path(temp);model='E_-2'
            source,ch,_=write_fixture(root,model,muF=10000.,muS=-11000.)
            point=split_scalar_point(root,model,source,ch)
            # E_-2 self species S2, source coordinate order: minus x3, plus x3
            diag=np.diag(point.mass2_boundary_muS)
            self.assertTrue(np.allclose(diag[6:9],[2.9e6+22000]*3))
            self.assertTrue(np.allclose(diag[9:12],[2.9e6-22000]*3))
            self.assertTrue(point.inputs['F_to_S_self_bilinear_sign_changed'])
            self.assertFalse(point.inputs['muF_matched_to_this_running_benchmark'])

    def test_floating_phase_unsupported_at_F(self):
        with TemporaryDirectory() as temp:
            root=Path(temp);model='A_-2'
            source,ch,_=write_fixture(root,model,muF=1000.+2j)
            with self.assertRaisesRegex(ValueError,'positive-real'):
                split_scalar_point(root,model,source,ch)
            source,ch,_=write_fixture(root,model,muF=-1000.)
            with self.assertRaisesRegex(ValueError,'positive-real'):
                split_scalar_point(root,model,source,ch)

    def test_source_field_order_mismatch_fails(self):
        with TemporaryDirectory() as temp:
            root=Path(temp);model='D_+0'
            source,ch,_=write_fixture(root,model)
            source.labels=tuple(reversed(source.labels))
            with self.assertRaisesRegex(ValueError,'source real field order'):
                split_scalar_point(root,model,source,ch)

    def test_eigenvalue_instability_tachyon_fails(self):
        with self.assertRaisesRegex(ValueError,'tachyonic'):
            native_original_mass2(1,3,2,1000.,500.,1000.)

    def test_missing_split_metadata_fails(self):
        with TemporaryDirectory() as temp:
            root=Path(temp);model='A_+0'
            source,ch,_=write_fixture(root,model)
            metadata=root/'output/general_observable/fermion_first_self_fix_v2/A_p0/self_mass_basis.json'
            metadata.unlink()
            with self.assertRaises(FileNotFoundError):
                split_scalar_point(root,model,source,ch)


class RunnerTests(TestCase):
    @staticmethod
    def loader():
        path=Path(__file__).resolve().parents[1]/'scripts/RunAll16PreliminaryNumericsV2.py'
        fake=ModuleType('Numerical.joint_scalar_matching.adapter')
        fake.ALL_MODELS=('A_-4','A_-2','E_+0')
        fake.normalize_model=lambda x:x if x in fake.ALL_MODELS else (_ for _ in ()).throw(ValueError('model'))
        with patch.dict(sys.modules,{'Numerical.joint_scalar_matching.adapter':fake}):
            spec=importlib.util.spec_from_file_location('new_all16_cli',path)
            module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        return module

    def test_v2_all_results_and_no_overwrites(self):
        cli=self.loader()
        with TemporaryDirectory() as temp:
            root=Path(temp);cli.ROOT=root
            cli.SAFE=root/'output/general_observable/lo_eft1_mixing'
            dest=cli.SAFE/'v2'
            def evaluate(_,model):
                return {'model':model,'status':'ConditionalWeinbergLambdaT3OneChannel',
                    'response':{'full_SU2_Weinberg_projection_pass':True,
                    'beta_SU2_residual_relative':2e-16,
                    'matching_scale_cancellation_relative':3e-16,
                    'neutrino_mass_singular_values_ev_per_abs_NGamma':[.01,.02,.03]}}
            with patch.object(cli,'evaluate_model',side_effect=evaluate),\
                 patch.object(cli,'evaluate_split_model',side_effect=evaluate):
                self.assertEqual(cli.main(['--output-root',str(dest),'--check-only']),0)
                self.assertFalse(dest.exists())
                self.assertEqual(cli.main(['--output-root',str(dest),'--require-all']),0)
                j=json.loads((dest/'summary.json').read_text())
                self.assertEqual(j['models_with_Ward_valid_one_channel_results'],3)
                self.assertEqual(j['six_self_bilinear_cartesian_hypotheses'],2)
                self.assertFalse(j['physical_neutrino_masses_computed'])
                self.assertTrue((dest/'summary.csv').is_file())
                for m in ('A_m4','A_m2','E_p0'):
                    self.assertTrue((dest/m/'lambdaT3_one_channel.json').is_file())
                with patch.object(cli.argparse.ArgumentParser,'error',side_effect=ValueError('overwrite')):
                    with self.assertRaisesRegex(ValueError,'overwrite'):
                        cli.main(['--output-root',str(dest)])

    def test_errors_and_nonweinberg_results_are_explicit(self):
        cli=self.loader()
        with TemporaryDirectory() as temp:
            root=Path(temp);cli.ROOT=root
            cli.SAFE=root/'output/general_observable/lo_eft1_mixing'
            def fail(_,model):
                if model=='A_-2': raise ValueError('source phase unverified')
                return {'model':model,'status':'NonWeinbergTensorBlockedMassInterpretation',
                    'response':{'full_SU2_Weinberg_projection_pass':False,
                    'beta_SU2_residual_relative':.2,
                    'matching_scale_cancellation_relative':0.,
                    'neutrino_mass_singular_values_ev_per_abs_NGamma':None}}
            dest=cli.SAFE/'v2'
            with patch.object(cli,'evaluate_model',side_effect=fail),\
                 patch.object(cli,'evaluate_split_model',side_effect=fail):
                self.assertEqual(cli.main(['--output-root',str(dest),'--require-all']),2)
                summary=json.loads((dest/'summary.json').read_text())
                self.assertEqual(summary['models_input_errors'],1)
                self.assertEqual(summary['models_non_Weinberg'],2)
                self.assertEqual(summary['models_with_Ward_valid_one_channel_results'],0)


if __name__=='__main__':
    import unittest
    unittest.main()
