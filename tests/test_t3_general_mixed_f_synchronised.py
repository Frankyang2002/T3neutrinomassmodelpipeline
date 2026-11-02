"""Regression tests for the NEW current-running B/C mass-basis handoff."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import types
import unittest
from unittest.mock import patch

import numpy as np

SOURCE = Path(__file__).resolve().parents[1] / 'scripts' / 'RunT3GeneralMixedFSynchronised.py'
spec = importlib.util.spec_from_file_location('t3_mixed_sync_under_test', SOURCE)
mod = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = mod
spec.loader.exec_module(mod)


def fixture(model='B_-1', mix=5+3j, m1=2.0e6, m2=3.0e6):
    state = {'mS1Sq':m1, 'mS2Sq':m2,
             mod.MIXING_NAME: {'real': mix.real, 'imag': mix.imag},
             'y1': [[0.01]*3]*3, 'y2': [[0.02]*3]*3,
             'MF': [[1e5,0,0],[0,1.2e5,0],[0,0,1.5e5]],
             'lambdaH':0.13}
    for k in ('lambdaS1','lambdaS2','lambdaH1','lambdaH2','lambda12',
              'lambdaT3','lambdaH1Adj','lambdaH2Adj','lambda12Adj',
              'lambdaHHdagS2S2','lambdaHHdagS1barS1bar',
              'lambdaS1bar2S2bar2','lambdaS1barS2S2bar2',
              'lambdaS1S1bar2S2bar','lambdaHHdagS1barS2barCross',
              'lambdaT3GeneralQ110101I1'):
        state[k] = .01
    return {'model':model,'status':'NumericalRenormalisableRunningOnly',
            'scales_gev':{'UV':1e7,'F':1e5,'S':1e3},
            'uv_at_fermion_threshold':state}


class TestBoundary(unittest.TestCase):
    def test_model_tag(self):
        self.assertEqual(mod.tag('B_-1'),'B_m1')
        self.assertEqual(mod.tag('C_-1'),'C_m1')
        with self.assertRaises(ValueError): mod.tag('A_-2')

    def test_number(self):
        self.assertEqual(mod.number({'real':5,'imag':-2},'mu'),5-2j)
        with self.assertRaises(ValueError): mod.number({'real':1},'mu')
        with self.assertRaises(ValueError): mod.number(float('nan'),'mu')

    def test_positive_real(self):
        self.assertEqual(mod.positive_real(150.0,'mass'),150.0)
        with self.assertRaises(ValueError): mod.positive_real(2j,'mass')
        with self.assertRaises(ValueError): mod.positive_real(-2.,'mass')

    def test_running_reads_only_fermion_endpoint(self):
        p=fixture()
        p['eft1_at_scalar_threshold']={'mS1Sq':-42}
        s,scale=mod.running_state(p,'B_-1')
        self.assertEqual(s['mS1Sq'],2e6)
        self.assertEqual(scale,1e5)

    def test_model_mismatch_fails(self):
        with self.assertRaisesRegex(ValueError,'unexpected'):
            mod.running_state(fixture(), 'C_-1')

    def test_missing_mu_f_fails(self):
        p=fixture();p['scales_gev'].pop('F')
        with self.assertRaisesRegex(ValueError,'matching scale'):
            mod.running_state(p,'B_-1')

    def test_missing_mixing_fails(self):
        p=fixture();p['uv_at_fermion_threshold'].pop(mod.MIXING_NAME)
        with self.assertRaisesRegex(ValueError,'missing'):
            mod.running_state(p,'B_-1')

    def test_file_selection_keeps_per_model(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            (root/'one'/'B_m1').mkdir(parents=True)
            (root/'one'/'B_m1'/'running.json').write_text('{}')
            with patch.object(mod,'ROOT',root):
                a=mod.running_file('B_-1',(Path('one'),))
                self.assertEqual(a,root/'one'/'B_m1'/'running.json')
                with self.assertRaises(FileNotFoundError):
                    mod.running_file('C_-1',(Path('one'),))

    def test_real_tensor_transport_exact_state_handoff_mock(self):
        """Mock the existing physics APIs, not a mass mixing or beta formula."""
        a=types.ModuleType('Numerical.experimental.GeneralT3MixedScalarMassBasis')
        class Basis:
            @classmethod
            def build(cls, model, m1, m2, mu):
                obj=cls()
                A=np.array([[m1,mu],[np.conj(mu),m2]])
                obj.eigenvalues_sq, obj.U=np.linalg.eigh(A)
                obj.tachyonic=bool(min(obj.eigenvalues_sq)<0)
                obj.degenerate=False
                obj.info=lambda: {'masses_gev':np.sqrt(obj.eigenvalues_sq).tolist(),
                                  'unitary_flavour_rotation_U':{'real':obj.U.real.tolist(),
                                                               'imag':obj.U.imag.tolist()}}
                return obj
        a.MixedDoubletMassBasis=Basis
        a.cjson=lambda z:{'real':complex(z).real,'imag':complex(z).imag}
        a.transport_operator=lambda fields,U:[{'fields':list(fields)}]
        a.verify_su2_covariance=lambda basis:0.
        a.verify_tensor_transport=lambda fields,U,trials:0.
        a.quadratic_mass_tensor_transport=lambda basis:{'verified':True}
        b=types.ModuleType('Numerical.experimental.GeneralT3ScalarQuadraticMasses')
        def hessian(manifest,model,m1,m2,mu):
            eig=np.linalg.eigvalsh([[m1,mu],[np.conj(mu),m2]])
            return {'mass_squared_eigenvalues_gev2':np.repeat(eig,4)}
        b.quadratic_mass_hessian=hessian
        old=types.ModuleType('scripts.RunGeneralT3MixedScalarMassBasis')
        old.STATUS='ExperimentalMassBasisTensorTransportNoPhysicalMatching'
        old.YUKAWA_SCALAR_LEGS={'y1':'S1dag','y2':'S2'}
        keys=list(fixture()['uv_at_fermion_threshold'])
        quartics=[k for k in keys if k.startswith('lambda')]
        self.assertEqual(len(quartics),17)
        old.quartic_operators=lambda manifest,registry,model:{k:('H','Hdag','S1dag','S1') for k in quartics}
        old._state_value=lambda state,name,sidecar:state[name]
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); (root/'source.json').write_text('test')
            (root/'registry.json').write_text('{}')
            (root/'manifest.json').write_text('{}')
            modules={'scripts': types.ModuleType('scripts'),
                     'scripts.RunGeneralT3MixedScalarMassBasis':old,
                     'Numerical': types.ModuleType('Numerical'),
                     'Numerical.experimental': types.ModuleType('Numerical.experimental'),
                     'Numerical.experimental.GeneralT3MixedScalarMassBasis':a,
                     'Numerical.experimental.GeneralT3ScalarQuadraticMasses':b,
                     'scripts.RunGeneralT3MixedScalarMassBasis':old}
            with patch.dict(sys.modules, modules):
                out=mod.make_mass_basis('B_-1',fixture(),root/'source.json',
                                        root/'manifest.json',root/'registry.json')
            self.assertEqual(out['scale_gev'],1e5)
            self.assertEqual(len(out['quartics']),17)
            self.assertEqual(len(out['yukawa_scalar_leg_transforms']),2)
            self.assertFalse(out['MatcheteRGBetaPhysicalCouplingMapVerified'])
            self.assertFalse(out['running_boundary_provenance']['historical_mixed_EFT1_endpoint_reused'])

    def test_prepare_only_does_not_call_mismatch_tree(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            f=root/'run'/'B_m1'/'running.json'; f.parent.mkdir(parents=True)
            f.write_text(json.dumps(fixture()),encoding='utf-8')
            mb={'scale_gev':1e5, 'basis':{'masses_gev':[1414.,1732.]}}
            fake_mat=types.ModuleType('scripts.RunGeneralT3MassEigenMatchete')
            calls=[]
            def prep(model, which, folder, outdir):
                calls.append((model,which,folder,outdir));outdir.mkdir(parents=True,exist_ok=True)
                target=outdir/'B_m1_final_parameters.json';target.write_text('{}')
                return {'parameters':target}
            fake_mat.prepare_model=prep
            fake_tree=types.ModuleType('scripts.RunGeneralT3FermionFirstEFT1Tree')
            fake_tree.run_model=lambda **kw: self.fail('Tree called in prepare-only mode')
            modules={'scripts':types.ModuleType('scripts'),
                     'scripts.RunGeneralT3MassEigenMatchete':fake_mat,
                     'scripts.RunGeneralT3FermionFirstEFT1Tree':fake_tree}
            with patch.dict(sys.modules,modules),patch.object(mod,'ROOT',root),\
                    patch.object(mod,'make_mass_basis',return_value=mb):
                result=mod.process('B_-1',(Path('run'),),root/'out',prepare_only=True)
            self.assertEqual(result['status'],'CurrentRunningMassBasisPrepared')
            self.assertEqual(len(calls),1)
            self.assertEqual(calls[0][1],'final')
            self.assertTrue((root/'out'/'B_m1'/'mass_basis'/'B_m1_final_mass_basis.json').is_file())

    def test_live_tree_requires_current_matchete_certificate(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            f=root/'run'/'C_m1'/'running.json';f.parent.mkdir(parents=True)
            f.write_text(json.dumps(fixture('C_-1')),encoding='utf-8')
            mb={'scale_gev':1e5, 'basis':{'masses_gev':[1400.,1600.]}}
            mat=types.ModuleType('scripts.RunGeneralT3MassEigenMatchete')
            tree=types.ModuleType('scripts.RunGeneralT3FermionFirstEFT1Tree')
            calls=[]
            def live(model,which,folder,outdir,wolframscript,timeout_sec):
                calls.append(('mass',model,folder,outdir))
                outdir.mkdir(parents=True,exist_ok=True)
                (outdir/'C_m1_final_matchete.json').write_text('{}')
                return {'parameters':outdir/'C_m1_final_parameters.json'}
            def tree_run(model,which,folder,certificate_dir,outdir,wolframscript,timeout_sec):
                calls.append(('tree',model,folder,certificate_dir))
                self.assertTrue((certificate_dir/'C_m1_final_matchete.json').is_file())
                return {'llss_file':'current_llss.txt','matched_eft_file':'current_tree.txt',
                        'llss_term_count':24}
            mat.run_model=live;tree.run_model=tree_run
            modules={'scripts':types.ModuleType('scripts'),
                     'scripts.RunGeneralT3MassEigenMatchete':mat,
                     'scripts.RunGeneralT3FermionFirstEFT1Tree':tree}
            with patch.dict(sys.modules,modules),patch.object(mod,'ROOT',root),\
                    patch.object(mod,'make_mass_basis',return_value=mb):
                record=mod.process('C_-1',(Path('run'),),root/'out')
            self.assertEqual([a[0] for a in calls],['mass','tree'])
            self.assertEqual(calls[1][3],calls[0][3])
            self.assertEqual(record['status'],'CurrentRunningMassBasisFThresholdTreeLLSSExported')
            self.assertFalse(record['numerical_flavor_wilson_matrix_computed'])
            self.assertEqual(record['llss_term_count'],24)

    def test_relative_output_root_normalised_before_native_matching(self):
        """Reproduces Windows relative_to(ROOT) failure after Wolfram success."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            f = root/'run'/'C_m1'/'running.json'
            f.parent.mkdir(parents=True)
            f.write_text(json.dumps(fixture('C_-1')), encoding='utf-8')
            mb = {'scale_gev': 1e5, 'basis': {'masses_gev': [1400., 1600.]}}
            mat = types.ModuleType('scripts.RunGeneralT3MassEigenMatchete')
            tree = types.ModuleType('scripts.RunGeneralT3FermionFirstEFT1Tree')
            def live(model, which, folder, outdir, wolframscript, timeout_sec):
                # Actual upstream validation: Path.relative_to does not resolve.
                (outdir/'C_m1_final_matchete.json').relative_to(root)
                return {'parameters': outdir/'C_m1_final_parameters.json'}
            def tree_run(model, which, folder, certificate_dir, outdir, wolframscript, timeout_sec):
                (certificate_dir/'C_m1_final_matchete.json').relative_to(root)
                (outdir/'C_m1_final_eft1_fermion_first_tree.json').relative_to(root)
                return {'llss_file': 'test_llss.txt', 'matched_eft_file': 'test_eft.txt',
                        'llss_term_count': 24}
            mat.run_model = live
            tree.run_model = tree_run
            modules = {'scripts': types.ModuleType('scripts'),
                       'scripts.RunGeneralT3MassEigenMatchete': mat,
                       'scripts.RunGeneralT3FermionFirstEFT1Tree': tree}
            with patch.dict(sys.modules, modules), patch.object(mod, 'ROOT', root), \
                    patch.object(mod, 'make_mass_basis', return_value=mb):
                record = mod.process('C_-1', (Path('run'),), Path('out'))
            self.assertEqual(record['status'], 'CurrentRunningMassBasisFThresholdTreeLLSSExported')
            self.assertTrue((root/'out'/'C_m1'/'status.json').is_file())

    def test_recover_existing_requires_certificates_and_same_running(self):
        """Recovery never reruns Wolfram or accepts a stale UV boundary."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            f = root/'run'/'B_m1'/'running.json'
            f.parent.mkdir(parents=True)
            f.write_text(json.dumps(fixture()), encoding='utf-8')
            basis = {'unitary_flavour_rotation_U': {'real': [[1.,0.],[0.,1.]],
                                                   'imag': [[0.,0.],[0.,0.]]},
                     'mass_squared_flavour_matrix_gev2': {'real': [[2e6,0.],[0.,3e6]],
                                                         'imag': [[0.,0.],[0.,0.]]},
                     'masses_gev': [1400., 1600.]}
            scalar = {'scale_gev': 1e5, 'basis': basis,
                      'source_sha256': {}, 'quartics': {'lambdaH': {'source_rgbeta_coupling':
                          {'real':.13,'imag':0.}}},
                      'running_boundary_provenance': {'running_sha256': mod.sha(f),
                                                       'boundary': 'uv_at_fermion_threshold',
                                                       'boundary_scale_gev': 1e5}}
            out = root/'out'/'B_m1'
            mod.save(out/'mass_basis'/'B_m1_final_mass_basis.json', scalar)
            mod.save(out/'matchete'/'B_m1_final_matchete.json', {'model':'B_-1'})
            mod.save(out/'fermion_tree'/'B_m1_final_eft1_fermion_first_tree.json',
                     {'llss_term_count': 20})
            mat = types.ModuleType('scripts.RunGeneralT3MassEigenMatchete')
            tree = types.ModuleType('scripts.RunGeneralT3FermionFirstEFT1Tree')
            mat.validate_wolfram_output = lambda obj, model, U: self.assertEqual(model, 'B_-1')
            mat.run_model = lambda **kwargs: self.fail('Live Matchete called during recovery')
            tree._stem = lambda model, which: 'B_m1_final_eft1_fermion_first_tree'
            tree.validate_tree_report = lambda obj, model, folder, stem: (folder/'tree.txt',folder/'llss.txt')
            tree.run_model = lambda **kwargs: self.fail('Live tree matching called during recovery')
            modules = {'scripts': types.ModuleType('scripts'),
                       'scripts.RunGeneralT3MassEigenMatchete': mat,
                       'scripts.RunGeneralT3FermionFirstEFT1Tree': tree}
            with patch.dict(sys.modules, modules), patch.object(mod, 'ROOT', root), \
                    patch.object(mod, 'make_mass_basis', return_value=scalar):
                result = mod.process('B_-1', (Path('run'),), Path('out'), recover_existing=True)
                self.assertTrue(result['reused_existing_validated_wolfram_outputs'])
                self.assertEqual(result['llss_term_count'], 20)
                # A changed running source cannot silently reuse old certificates.
                f.write_text(json.dumps(fixture(mix=12+2j)), encoding='utf-8')
                with self.assertRaisesRegex(ValueError, 'previous mass-basis input differs'):
                    mod.process('B_-1', (Path('run'),), Path('out'), recover_existing=True)

    def test_summary_failure_is_not_success(self):
        with tempfile.TemporaryDirectory() as tmp:
            with patch.object(mod,'process',side_effect=ValueError('unsafe basis')):
                ret=mod.main(['--all','--output-root',tmp])
            self.assertEqual(ret,1)
            summary=mod.read(Path(tmp)/'summary.json')
            self.assertEqual(summary['live_tree_success_count'],0)
            self.assertTrue(all(v['status']=='Blocked' for v in summary['results']))


if __name__=='__main__':
    unittest.main()
