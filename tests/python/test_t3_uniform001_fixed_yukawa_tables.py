"""Isolated tests: no UV running, scalar matching or refit performed here."""
from __future__ import annotations

import csv
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import sys
import types

import numpy as np

from scripts import BuildT3Uniform001FixedYukawaTables as report


def enc(a):
    a=np.asarray(a,dtype=complex)
    return {'real':a.real.tolist(),'imag':a.imag.tolist()}


def fake_result(model,ref, c5, *, chi2=1.0):
    return {'schema':report.SCHEMA,'model':model,'yukawa_reference':ref,
            'status':'Success','extra_scalar_quartics_uv':0.01,
            'yukawa_conjugate_uv_zero':True, 'C5_low_gev_inv':enc(c5),
            'neutrino_masses_ev':[0.001,0.0087,0.05],
            'chi2':chi2,'prediction':{
                'delta_m21_sq_ev2':7.537e-5,'delta_m3l_sq_ev2':2.511e-3,
                'sin2_theta12':.3088,'sin2_theta13':.02248,'sin2_theta23':.47}}


def write_pinned_fixture(folder: Path, *, model: str='A_-2') -> dict:
    """Synthetic canonical export, with original Windows machine paths."""
    folder.mkdir(parents=True, exist_ok=True)
    lag={'model':model}
    for key,filename in report.PINNED_INPUTFORMS.items():
        lag[key]=rf'C:\thesistemp\T3neutrinomassmodelpipeline\output\studies\general_scalar_extended\complete_scalar_symbolic\full16_general_v1\A_m2\{filename}'
        (folder/filename).write_text('Nonempty inputform: '+filename,encoding='utf-8')
    (folder/'lagrangian.json').write_text(json.dumps(lag),encoding='utf-8')
    (folder/'uv_rgbeta.json').write_text('UV beta fake',encoding='utf-8')
    (folder/'eft1_rgbeta.json').write_text('EFT1 beta fake',encoding='utf-8')
    (folder/'provenance.json').write_text(json.dumps({'model':model,'source_sha256':{'example.wl':'pinned-hash'}}),encoding='utf-8')
    sources={sha_name:report.sha(folder/filename) for filename,sha_name in report.PINNED_STATUS_SHA256.items()}
    (folder/'status.json').write_text(json.dumps({'sources':sources}),encoding='utf-8')
    return lag


class Tests(unittest.TestCase):
    def test_complex_proportionality_exact(self):
        c=np.array([[1,2j,3],[2j,4,5j],[3,5j,6]],complex)*1e-16
        k,epsilon=report.proportion((-2+3j)*c,c)
        self.assertAlmostEqual(k.real,-2,12)
        self.assertAlmostEqual(k.imag,3,12)
        self.assertLess(epsilon,1e-14)

    def test_orthogonal_c5(self):
        a=np.diag([1,0,0]).astype(complex)*1e-16
        b=np.diag([0,1,0]).astype(complex)*1e-16
        k,eps=report.proportion(b,a)
        self.assertEqual(k,0j)
        self.assertAlmostEqual(eps,1)

    def test_complex_C5_matrix_shape_and_symmetry(self):
        c=np.array([[2,1j,0],[1j,3,0],[0,0,1]],complex)*1e-16
        np.testing.assert_allclose(report.c5_from_json(enc(c)),c)
        with self.assertRaises(ValueError):
            report.c5_from_json(enc(np.ones((2,2))))
        with self.assertRaises(ValueError):
            report.c5_from_json(enc(np.array([[1,2,0],[0,1,0],[0,0,1]])*1e-16))

    def test_reference_diagonal_reuse(self):
        base={'chi2':.005,'prediction':dict(fake_result('A_-4','A_-4',np.eye(3))['prediction']),
              'neutrino_masses_ev':[1e-3,8e-3,5e-2],
              'numerical_details':{'c5_at_low_scale':enc(np.eye(3)*1e-16)}}
        sig={'lambdaT3':.1}
        result=report.reuse_diagonal(base,'A_-4','A_-4',sig)
        self.assertEqual(result['status'],'ReusedReference')
        self.assertEqual(result['chi2'],.005)
        self.assertEqual(result['C5_low_gev_inv']['real'][0][0],1e-16)

    def test_report_all_four_tables(self):
        with tempfile.TemporaryDirectory() as root:
            root=Path(root)
            data=root/'data';out=root/'reports'
            for model in report.MODELS:
                for ref in report.SOURCE_REFS:
                    scale=1.0 if ref=='D_-2' else 2.0
                    c=scale*np.diag([1.,2.,3.])*1e-16
                    # Intentionally perturb model-specific C5 with a second tensor.
                    if model!=ref:
                        c=c+np.diag([0.,0.,0.02])*1e-16
                    p=data/report.tag(model)/report.tag(ref)/'result.json'
                    p.parent.mkdir(parents=True,exist_ok=True)
                    p.write_text(json.dumps(fake_result(model,ref,c)),encoding='utf-8')
            meta=report.build_tables(data,out)
            self.assertEqual(meta['points'],32)
            for name in ('table_C1_Dm2','table_C2_Am4','table_C3_PMNS',
                         'table_C4_Weinberg_proportionality'):
                for ext in ('.csv','.md','.tex'):
                    self.assertTrue((out/(name+ext)).exists())
            for name,num in [('table_C1_Dm2',16),('table_C2_Am4',16),
                             ('table_C3_PMNS',32),('table_C4_Weinberg_proportionality',16)]:
                with (out/(name+'.csv')).open() as f:
                    self.assertEqual(len(list(csv.DictReader(f))),num)
            with (out/'table_C5_Weinberg_matrix_entries.csv').open() as f:
                self.assertEqual(len(list(csv.DictReader(f))),32*9)
            with (out/'table_C1_Dm2.csv').open() as f:
                rows={x['model']:x for x in csv.DictReader(f)}
            self.assertAlmostEqual(float(rows['D(-2)']['epsilon_C5']),0,places=12)
            with (out/'table_C2_Am4.csv').open() as f:
                rows={x['model']:x for x in csv.DictReader(f)}
            self.assertAlmostEqual(float(rows['A(-4)']['epsilon_C5']),0,places=12)
            with (out/'table_C4_Weinberg_proportionality.csv').open() as f:
                w={r['model']:r for r in csv.DictReader(f)}
            self.assertAlmostEqual(float(w['D(-2)']['k_Dm2_real']),1,places=12)
            self.assertAlmostEqual(float(w['A(-4)']['k_Am4_real']),1,places=12)
            self.assertAlmostEqual(float(w['D(-2)']['epsilon_Dm2']),0,places=12)
            self.assertTrue((out/'provenance.json').exists())

    def test_reject_missing_transfer_no_fake_rows(self):
        with tempfile.TemporaryDirectory() as root:
            with self.assertRaises(FileNotFoundError):
                report.build_tables(Path(root),Path(root)/'reports')

    def test_new_uv_policy_rejected_if_wrong(self):
        with tempfile.TemporaryDirectory() as root:
            d=Path(root)
            model='D_-2';ref='D_-2'
            p=d/report.tag(model)/report.tag(ref)/'result.json'
            p.parent.mkdir(parents=True)
            row=fake_result(model,ref,np.eye(3)*1e-16)
            row['extra_scalar_quartics_uv']=0.05
            p.write_text(json.dumps(row))
            with self.assertRaisesRegex(ValueError,'Wrong benchmark'):
                report.finished_rows(d,(model,),(ref,))

    def test_transfer_yukawas_unchanged_without_minimisation(self):
        frozen=np.arange(1,19,dtype=float)*1e-5
        seen={}
        module=types.ModuleType('scripts.RunT3QualitativeScalarRefits')
        module.fitted_vector=lambda cfg: frozen.copy()
        module.yukawa_matrices=lambda x:{'y1':enc(x[:9].reshape(3,3)),
                                          'y2':enc(x[9:].reshape(3,3))}
        module.parameter_dict=lambda x:{f'y{i//9+1}_{(i%9)//3}{i%3}':float(x[i])
                                        for i in range(18)}
        module.fit_prediction=lambda fit: fake_result('A_-4','D_-2',np.eye(3))['prediction']
        module.low_masses=lambda details,cfg,order:[.001,.008,.05]
        class FakeEvaluator:
            def evaluate(self,params,*,detail):
                self.detail=detail
                seen.update(params)
                fit=types.SimpleNamespace(chi2=7.0,ordering='NO')
                return fit,{'c5_at_low_scale':enc(np.eye(3)*1e-16)}
        class Factory:
            def build(self,lam):
                seen['lambda_t3']=lam
                cfg={'base_state':{'ordinary':{'y1':'old','y2':'old'}}}
                return FakeEvaluator(),cfg,None,{'extra_quartic_value':.01}
        with patch.dict(sys.modules,{'scripts.RunT3QualitativeScalarRefits':module}):
            row=report.evaluate_fixed(Factory(),'A_-4','D_-2',{'fitted':'x'},
                                      {'lambdaT3':.1})
        self.assertEqual(row['status'],'Success')
        self.assertEqual(row['chi2'],7.0)  # A fixed transfer can be a bad oscillation fit.
        self.assertEqual(seen['lambda_t3'],.1)
        self.assertEqual(seen['y1_00'],float(frozen[0]))
        self.assertEqual(seen['y2_22'],float(frozen[17]))
        self.assertEqual(row['uv_transferred_yukawas']['y1']['real'][0][0],float(frozen[0]))

    def test_frozen_symbolic_baseline_reproduces_saved_matrix(self):
        """The pinned cache must reproduce both reference observables AND C5."""
        expected=fake_result('A_-2','A_-2',np.eye(3)*1e-16)
        source={'prediction':expected['prediction'],
                'numerical_details':{'c5_at_low_scale':expected['C5_low_gev_inv']}}
        base=types.ModuleType('scripts.RunT3QualitativeScalarRefits')
        base.real_lambda=lambda cfg:0.1
        base.fitted_vector=lambda cfg:np.arange(18,dtype=float)*1e-5
        base.parameter_dict=lambda x:{'vector':x}
        base.fit_prediction=lambda fit:expected['prediction']
        class E:
            def evaluate(self,params,*,detail):
                self.assert_detail=detail
                return types.SimpleNamespace(chi2=0.0),{'c5_at_low_scale':enc(np.eye(3)*1e-16)}
        class F:
            model='A_-2'
            def build(self,lam):
                self.lam=lam
                return E(),None,None,None
        with patch.dict(sys.modules,{'scripts.RunT3QualitativeScalarRefits':base}):
            report.verify_pinned_baseline(F(),source,{'fitted':'y'})

    def test_frozen_symbolic_baseline_rejects_nonreproducible_c5(self):
        expected=fake_result('A_-2','A_-2',np.eye(3)*1e-16)
        source={'prediction':expected['prediction'],
                'numerical_details':{'c5_at_low_scale':expected['C5_low_gev_inv']}}
        base=types.ModuleType('scripts.RunT3QualitativeScalarRefits')
        base.real_lambda=lambda cfg:0.1
        base.fitted_vector=lambda cfg:np.arange(18,dtype=float)*1e-5
        base.parameter_dict=lambda x:{'vector':x}
        base.fit_prediction=lambda fit:expected['prediction']
        class E:
            def evaluate(self,params,*,detail):
                return types.SimpleNamespace(chi2=0.0),{'c5_at_low_scale':enc(np.eye(3)*2e-16)}
        class F:
            model='A_-2'
            def build(self,lam):return E(),None,None,None
        with patch.dict(sys.modules,{'scripts.RunT3QualitativeScalarRefits':base}):
            with self.assertRaisesRegex(ValueError,'C5 baseline differs'):
                report.verify_pinned_baseline(F(),source,{'fitted':'y'})

    def test_pinned_factory_uses_cached_source_stamp_and_retains_validator(self):
        """A changed live source must not be relabelled as matched."""
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            folder=root/'full16_general_v1'/'A_m2'
            folder.mkdir(parents=True)
            stamp={'model':'A_-2','source_sha256':{'example.wl':'pinned-hash'}}
            write_pinned_fixture(folder)
            # Keep the fixed source-stamp used by this factory unit test.
            (folder/'provenance.json').write_text(json.dumps(stamp))
            seen=[]
            canonical=types.ModuleType('scripts.RunT3GeneralSymbolicStudy')
            canonical.DEFAULT_DATA_ROOT=root
            canonical.model_tag=lambda m:m.replace('-','m').replace('+','p')
            canonical.load_json=lambda p: (stamp if p.name=='provenance.json' else
                    {'models':{'A_-2':{'model':'A_-2'}}} if p.name=='registry.json' else
                    {'report_betas':{'y1':'x'}} if p.name=='uv.json' else
                    {'report_betas':{'lambdaS1':'x'}} if p.name=='eft.json' else
                    {'lambdaS1':0.01})
            canonical.fingerprint=lambda:{'example.wl':'edited-live-hash'}
            canonical.expected_yukawa=lambda m:None
            def validate_rebased(model, folder, entry, kind, hashes):
                lag=report.load(folder/'lagrangian.json')
                for key in report.PINNED_INPUTFORMS:
                    self.assertTrue(Path(lag[key]).is_file())
                seen.append((model,dict(hashes)))
            canonical.validate_cached=validate_rebased
            scalar=types.ModuleType('scripts.RunT3CompleteGeneralScalarSymbolic')
            scalar.REGISTRY=root/'registry.json'
            scalar.MANIFEST=root/'manifest.json'
            general=types.ModuleType('scripts.RunT3GeneralAll16FromReference')
            general.ensure_general_matching=lambda *args,**kwargs:root/'matching.txt'
            class E:
                def __init__(self,*args,**kwargs):self.args=args
            general.NonMixedEvaluator=E
            uniform=types.ModuleType('scripts.RunT3UniformScalarBenchmark')
            uniform.QUARTIC=0.01
            uniform.uniform_config=lambda cfg,lam: (
                {'base_state':{'ordinary':{'lambdaT3':{'real':lam,'imag':0.0}}}},
                 ['lambdaS1'])
            uniform.uniform_overrides=lambda old:({'lambdaS1':0.01},['lambdaS1'],[])
            audits=types.ModuleType('scripts.AuditGeneralT3QuadraticBasis')
            audits.audit=lambda:{'models':{'A_-2':True}}
            quad=types.ModuleType('Numerical.experimental.GeneralT3ScalarQuadraticMasses')
            quad.load_manifest=lambda path:{}
            quad.quadratic_kind=lambda data,model:'none'
            physics=types.ModuleType('Numerical.matching.GeneralT3FlavourC5')
            class K:
                def __init__(self,*args,**kwargs):self.required_yukawas={'y1','y2'}
            physics.GeneralT3FlavourC5=K
            mods={
                'scripts.RunT3GeneralSymbolicStudy':canonical,
                'scripts.RunT3CompleteGeneralScalarSymbolic':scalar,
                'scripts.RunT3GeneralAll16FromReference':general,
                'scripts.RunT3UniformScalarBenchmark':uniform,
                'scripts.AuditGeneralT3QuadraticBasis':audits,
                'Numerical.experimental.GeneralT3ScalarQuadraticMasses':quad,
                'Numerical.matching.GeneralT3FlavourC5':physics}
            paths={'uv_beta':root/'uv.json','eft1_beta':root/'eft.json',
                   'matching':root/'matching.txt','overrides':root/'overrides.json'}
            with patch.dict(sys.modules,mods):
                factory=report.PinnedCanonicalNonMixedFactory(
                    'A_-2',{},paths,object(),'full16_general_v1')
                evaluator,changed,fixed,boundary=factory.build(0.1)
            self.assertEqual(seen,[('A_-2',{'example.wl':'pinned-hash'})])
            self.assertEqual(factory.changed_symbolic_sources,['example.wl'])
            self.assertFalse(boundary['live_symbolic_source_fingerprint_matches'])
            self.assertEqual(boundary['symbolic_provenance_mode'],'pinned-canonical-export')
            self.assertEqual(evaluator.args[4]['lambdaT3']['real'],0.1)
            self.assertEqual(fixed['lambdaS1'],0.01)

    def test_pinned_absolute_inputform_paths_relocated_without_editing_cache(self):
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d)/'full16_general_v1'/'A_m2'
            original=write_pinned_fixture(folder)
            validator=types.SimpleNamespace()
            seen=[]
            def check(model, shadow, entry, kind, pinned):
                metadata=report.load(shadow/'lagrangian.json')
                self.assertEqual(model,'A_-2')
                self.assertEqual(kind,'none')
                self.assertEqual(pinned,{'example.wl':'pinned-hash'})
                for key,name in report.PINNED_INPUTFORMS.items():
                    self.assertEqual(Path(metadata[key]).read_text(encoding='utf-8'),
                                     'Nonempty inputform: '+name)
                    self.assertEqual(Path(metadata[key]).name,name)
                self.assertNotEqual(shadow,folder)
                seen.append(True)
            validator.validate_cached=check
            report.validate_relocated_pinned_cache(validator,'A_-2',folder,{},
                                                   'none',{'example.wl':'pinned-hash'})
            self.assertEqual(seen,[True])
            self.assertEqual(report.load(folder/'lagrangian.json'),original)

    def test_pinned_inputform_hash_mismatch_still_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d)
            write_pinned_fixture(folder)
            (folder/'lagrangian.json.interaction.inputform.txt').write_text('changed')
            validator=types.SimpleNamespace(validate_cached=lambda *args: self.fail('must not validate'))
            with self.assertRaisesRegex(ValueError,'SHA256 mismatch'):
                report.validate_relocated_pinned_cache(validator,'A_-2',folder,{},
                                                       'none',{'example.wl':'pinned-hash'})

    def test_pinned_file_pointer_cannot_be_redirected(self):
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d)
            lag=write_pinned_fixture(folder)
            lag['interaction_inputform_file']=r'C:\other-machine\unrelated.inputform.txt'
            (folder/'lagrangian.json').write_text(json.dumps(lag),encoding='utf-8')
            # Hash the deliberately tampered metadata to exercise pointer check.
            status=report.load(folder/'status.json')
            status['sources']['lagrangian_sha256']=report.sha(folder/'lagrangian.json')
            (folder/'status.json').write_text(json.dumps(status),encoding='utf-8')
            validator=types.SimpleNamespace(validate_cached=lambda *args: self.fail('must not validate'))
            with self.assertRaisesRegex(ValueError,'unexpected pinned inputform filename'):
                report.validate_relocated_pinned_cache(validator,'A_-2',folder,{},
                                                       'none',{'example.wl':'pinned-hash'})

    def test_pinned_missing_inputform_remains_blocked(self):
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d)
            write_pinned_fixture(folder)
            (folder/'lagrangian.json.bsm.inputform.txt').unlink()
            validator=types.SimpleNamespace(validate_cached=lambda *args: self.fail('must not validate'))
            with self.assertRaisesRegex(FileNotFoundError,'missing pinned symbolic inputform'):
                report.validate_relocated_pinned_cache(validator,'A_-2',folder,{},
                                                       'none',{'example.wl':'pinned-hash'})

    def test_pinned_mode_is_explicit_not_default(self):
        import inspect
        body=inspect.getsource(report.main)
        self.assertIn("'--frozen-symbolic'",body)
        self.assertIn('campaign.validated_reference',body)
        self.assertIn('verify_pinned_baseline',body)
        self.assertIn('PinnedCanonicalNonMixedFactory',body)

    def test_source_markdown_has_fixed_yukawa_warnings(self):
        self.assertIn('fixed-Yukawa',report.__doc__)
        self.assertIn('conditional',report.__doc__)
        self.assertIn('ReusedReference',report.reuse_diagonal.__code__.co_consts)


if __name__=='__main__':
    unittest.main()
