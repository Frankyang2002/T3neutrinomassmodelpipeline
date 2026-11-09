"""Tests for the v10 clean-process resume pathway (not a Wolfram simulation)."""
from __future__ import annotations

from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path
import tempfile
import unittest

from scripts.RunT3ParallelAScalarBatchCheckV10V12 import (
    ResumeBlocked, compare_saved_to_fresh, verify_smoke,
)


class TestSavedWitnessRevalidation(unittest.TestCase):
    def test_nested_source_data_passes(self):
        doc={'format':'v11','source_sha256':{'model':'deadbeef'},
             'models':{'A_-4':{'rank':2,'values':[1.,2.,{'x':3.}]}},
             'physical_Wilson_evolution_authorized':False}
        compare_saved_to_fresh(deepcopy(doc),doc)

    def test_modified_saved_physics_rejected(self):
        good={'models':{'A_-4':{'rank':2,'error':0.0}}}
        bad=deepcopy(good);bad['models']['A_-4']['rank']=1
        with self.assertRaisesRegex(ResumeBlocked,'saved numeric evidence differs'):
            compare_saved_to_fresh(bad,good)

    def test_source_sha_tampering_rejected(self):
        good={'source_sha256':{'model':'one'}}
        with self.assertRaisesRegex(ResumeBlocked,'saved evidence differs'):
            compare_saved_to_fresh({'source_sha256':{'model':'two'}},good)

    def test_nan_and_physics_authorization_rejected(self):
        good={'status':True,'number':1.}
        with self.assertRaisesRegex(ResumeBlocked,'NaN'):
            compare_saved_to_fresh({'status':True,'number':float('nan')},good)
        with self.assertRaisesRegex(ResumeBlocked,'value type differs'):
            compare_saved_to_fresh({'status':1,'number':1.},good)

    def test_v10_smoke_rejects_wrong_status(self):
        with self.assertRaisesRegex(ResumeBlocked,'did not pass'):
            verify_smoke({'format':'t3-parallel-a-native-rgbeta-quartic-projectors-v10',
                          'smoke_only':True,'status':'Blocked_NativeRGBetaQuarticProjectorResponses'},
                         Path('.'))

    def test_native_smoke_source_hash_and_response_guard(self):
        from Numerical.parallel_t3.wilson.native_quartic_projectors_v10 import (
            FORMAT, SOURCE_FILES, DIMENSIONS, expected_symbols,
        )
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            hashes={}
            for name in SOURCE_FILES:
                path=root/name
                path.parent.mkdir(parents=True,exist_ok=True)
                path.write_bytes(name.encode())
                hashes[name]=sha256(path.read_bytes()).hexdigest()
            package=root/'installed/RGBeta.m'
            package.parent.mkdir(parents=True,exist_ok=True)
            package.write_bytes(b'fake immutable package fixture')
            models={}
            for tag in ('A_-4','B_-1','E_-2'):
                d1,d2,df=DIMENSIONS[tag[0]]
                symbols=expected_symbols(tag)
                inventory=[{
                  'symbol':name,'native_self_conjugate':is_real,
                  'has_installed_RGBeta_projector':True,
                  'native_unique_arrangement_count':1,
                  'native_coupling_inputform':f'6*Matrix[{name}][##1] &'
                } for name,is_real in symbols.items()]
                diag=[{'symbol':name,'native_response':{'ok':True,'value':1.}}
                      for name,is_real in symbols.items() if is_real]
                blocks=[]
                for which,d in ((1,d1),(2,d2)):
                    names=[f'lambdaH{which}']+([f'lambdaH{which}Adj'] if d>1 else [])
                    blocks.append({'which':which,'dimension':d,'symbols':names,
                       'native_response_matrix':[
                         [{'ok':True,'value':float(i==j)} for j in range(len(names))]
                         for i in range(len(names))]})
                models[tag]={'ok':True,'dimensions':[d1,d2,df],
                             'alpha':int(tag.split('_')[1]),
                             'registered_quartics_count':len(symbols),
                             'registry_exact_source_inventory':True,
                             'inventory':inventory,'native_real_self_projectors':diag,
                             'native_portal_response_blocks':blocks,
                             'native_real_self_projectors_all_unit':True,
                             'native_portal_response_blocks_evaluated':True,
                             'native_numeric_S2_components_measured':False}
            witness={'format':FORMAT,
                     'status':'NativeRGBetaQuarticProjectorSmoke3ModelsPassed',
                     'smoke_only':True,'registered_model_count':3,
                     'native_Lam_evaluated':True,'native_quartic_projectors_evaluated':True,
                     'physical_Wilson_evolution_authorized':False,
                     'source_sha256':hashes,
                     'installed_RGBeta_main_file':str(package),
                     'installed_RGBeta_main_file_sha256':sha256(package.read_bytes()).hexdigest(),
                     'models':models}
            self.assertEqual(verify_smoke(witness,root)['status'],witness['status'])
            broken=deepcopy(witness)
            broken['models']['B_-1']['inventory'][0]['native_coupling_inputform']='12*Matrix[lambdaH][##1] &'
            with self.assertRaisesRegex(ResumeBlocked,'registration changed'):
                verify_smoke(broken,root)
            broken=deepcopy(witness)
            broken['models']['B_-1']['native_portal_response_blocks'][0][
              'native_response_matrix'][1][1]['value']=0.
            with self.assertRaisesRegex(ResumeBlocked,'response singular'):
                verify_smoke(broken,root)
            package.write_bytes(b'tampered fake package')
            with self.assertRaisesRegex(ResumeBlocked,'entrypoint changed'):
                verify_smoke(witness,root)

    def test_runner_enforces_scoped_environment_and_resume(self):
        runner=(Path(__file__).resolve().parents[1] /
                'scripts/RunT3ParallelAScalarBatchV10V12.ps1').read_text(encoding='utf-8')
        wolfram=(Path(__file__).resolve().parents[1] /
                 'scripts/RunT3ParallelANativeQuarticProjectorsV10.wl').read_text(encoding='utf-8')
        self.assertIn('finally {',runner)
        self.assertIn("$old[$flag]",runner)
        self.assertIn("[Environment]::SetEnvironmentVariable($flag, $null, 'Process')",runner)
        self.assertIn("[Environment]::SetEnvironmentVariable($flag, $old[$flag], 'Process')",runner)
        self.assertIn('Verify-Existing',runner)
        self.assertIn('v10-smoke',runner)
        self.assertIn('v10-native',runner)
        self.assertIn('v10-final',runner)
        self.assertIn('T3A10EnabledFlags',wolfram)
        self.assertIn('Exit[3]',wolfram)


if __name__ == '__main__':
    unittest.main()
