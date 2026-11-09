"""Regression tests for fail-closed native v10 smoke forensics (v13)."""
from __future__ import annotations

from copy import deepcopy
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

from Numerical.parallel_t3.wilson.native_quartic_projectors_v10 import (
    FORMAT, SOURCE_FILES, expected_symbols, DIMENSIONS,
)

ROOT=Path(__file__).resolve().parents[1]
SPEC=importlib.util.spec_from_file_location(
    'v13forensic',ROOT/'scripts/RunT3ParallelAQuarticForensicsV13.py')
assert SPEC is not None and SPEC.loader is not None
M=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)


def fixture(root:Path):
    paths={}
    for i,name in enumerate(SOURCE_FILES):
        p=root/name;p.parent.mkdir(parents=True,exist_ok=True)
        p.write_bytes(('source text '+str(i)).encode())
        paths[name]=sha256(p.read_bytes()).hexdigest()
    entry=root/'RGBeta_main.m'
    entry.write_bytes(b'RGBeta v1.2.0 fixture')
    rows={}
    for tag in M.EXPECTED:
        expected=expected_symbols(tag)
        dims=DIMENSIONS[tag[0]]
        diag=[{'symbol':name,'native_response':{'ok':True,'value':1.}}
              for name,real in expected.items() if real]
        blocks=[]
        for which,d in enumerate(dims[:2],1):
            ss=[f'lambdaH{which}']+([f'lambdaH{which}Adj'] if d>1 else [])
            matrix=[[{'ok':True,'value':1. if i==j else 0.}
                     for j in range(len(ss))] for i in range(len(ss))]
            blocks.append({'which':which,'dimension':d,'symbols':ss,
                           'native_response_matrix':matrix})
        rows[tag]={
            'ok':False,
            'registry_exact_source_inventory':True,
            'registered_quartics_count':len(expected),
            'inventory':[{'symbol':k} for k in expected],
            'native_real_self_projectors':diag,
            'native_portal_response_blocks':blocks,
            'native_real_self_projectors_all_unit':False,
            'native_portal_response_blocks_evaluated':True,
        }
    rows['A_-4']['native_real_self_projectors'][0]['native_response']={
        'ok':False,'reason':'unresolved_native_symbolic_coefficient',
        'coefficient_inputform':'del[SU2L[S2], i, j]'}
    rows['B_-1']['native_real_self_projectors'][0]['native_response']={
        'ok':True,'value':2.}
    rows['E_-2']['native_portal_response_blocks'][0]['native_response_matrix'][0][0]={
        'ok':False,'reason':'native_refinement_failed_or_timed_out'}
    return {
        'format':FORMAT,'status':M.BLOCKED,'smoke_only':True,
        'registered_model_count':3,
        'physical_Wilson_evolution_authorized':False,
        'native_numeric_S2_frame_measured':False,
        'full_scalar_potential_dictionary_certified':False,
        'source_sha256':paths,
        'installed_RGBeta_main_file':str(entry),
        'installed_RGBeta_main_file_sha256':sha256(entry.read_bytes()).hexdigest(),
        'models':rows,
    }


class TestNativeV13FailureForensics(unittest.TestCase):
    def test_classifies_without_authorizing_physics(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td); doc=fixture(root)
            res=M.analyze(doc,root,'a'*64)
            self.assertEqual(res['status'],M.STATUS)
            self.assertEqual(res['input_smoke_sha256'],'a'*64)
            self.assertFalse(res['native_Lam_projector_certified'])
            self.assertFalse(res['physical_Wilson_evolution_authorized'])
            self.assertEqual(len(res['minimal_reprobe_targets']),3)
            self.assertEqual(res['response_class_counts']['native_evaluation_failed'],2)
            self.assertGreaterEqual(res['response_class_counts']['non_unit_self_response'],1)
            self.assertEqual(res['attempted_by_type']['diagonal'],29)
            self.assertEqual(res['attempted_by_type']['portal'],
                             sum((1 if d==1 else 4) for t in M.EXPECTED
                                 for d in DIMENSIONS[t[0]][:2]))

    def test_wrong_registered_symbols_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);doc=fixture(root)
            doc['models']['E_-2']['inventory'].pop()
            with self.assertRaisesRegex(M.ForensicBlocked,'inventory'):
                M.analyze(doc,root,'a'*64)

    def test_modified_source_hash_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);doc=fixture(root)
            (root/SOURCE_FILES[0]).write_text('changed')
            with self.assertRaisesRegex(M.ForensicBlocked,'hash'):
                M.analyze(doc,root,'a'*64)

    def test_modified_installed_rgbeta_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);doc=fixture(root)
            (root/'RGBeta_main.m').write_text('changed')
            with self.assertRaisesRegex(M.ForensicBlocked,'installed RGBeta'):
                M.analyze(doc,root,'a'*64)

    def test_false_physics_certification_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);doc=fixture(root)
            doc['physical_Wilson_evolution_authorized']=True
            with self.assertRaisesRegex(M.ForensicBlocked,'physics'):
                M.analyze(doc,root,'a'*64)

    def test_status_must_be_blocked(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);doc=fixture(root)
            doc['status']='NativeRGBetaQuarticProjectorSmoke3ModelsPassed'
            with self.assertRaisesRegex(M.ForensicBlocked,'blocked'):
                M.analyze(doc,root,'a'*64)

    def test_cli_writes_separate_idempotent_and_immutable_output(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);doc=fixture(root)
            saved=root/'output/parallel_t3/wilson/rgbeta_quartic_projectors_native_v10_smoke.json'
            saved.parent.mkdir(parents=True,exist_ok=True)
            original=json.dumps(doc,indent=2).encode()
            saved.write_bytes(original)
            self.assertEqual(M.main(['--project-root',str(root)]),0)
            h=sha256(original).hexdigest()[:12]
            forensic=saved.parent/f'rgbeta_quartic_forensics_v13_{h}.json'
            self.assertTrue(forensic.is_file())
            self.assertEqual(saved.read_bytes(),original)
            self.assertEqual(M.main(['--project-root',str(root)]),0)
            forensic.write_text('{"bad":true}')
            self.assertEqual(M.main(['--project-root',str(root)]),2)

    def test_runner_preserves_original_failed_smoke(self):
        text=(ROOT/'scripts/RunT3ParallelAScalarBatchV10V12.ps1').read_text()
        for snippet in ('Verify-Existing','v10-smoke','v10-native','v10-final',
              'Blocked_NativeRGBetaQuarticProjectorResponses',
              'scripts.RunT3ParallelAQuarticForensicsV13',
              'RunT3ParallelAProbeQuarticResponseV13.wl',
              'finally {',
              "[Environment]::SetEnvironmentVariable($flag, $null, 'Process')",
              "[Environment]::SetEnvironmentVariable($flag, $old[$flag], 'Process')"):
            self.assertIn(snippet,text)
        self.assertIn('exit 7',text)

    def test_reprobe_source_contains_two_independent_routes_and_hard_blocks(self):
        text=(ROOT/'scripts/RunT3ParallelAProbeQuarticResponseV13.wl').read_text()
        for snippet in ('Block[{$quartics=KeyTake[$quartics,{source}]}',
                       'full_registry_then_other_couplings_zero',
                       'T3Stage22NativeGet','T3RGBetaBuildEFT1',
                       'physical_Wilson_evolution_authorized',
                       'Exit[3]', 'Exit[6]', 'source_sha256'):
            self.assertIn(snippet,text)


if __name__=='__main__':unittest.main()

class TestNativeReprobeComparison(unittest.TestCase):
    @staticmethod
    def fixture_pair():
        forensic={
            'status':'NativeQuarticSmokeFailureClassified_NotCertified',
            'physical_Wilson_evolution_authorized':False,
            'installed_RGBeta_main_file_sha256':'9'*64,
            'minimal_reprobe_targets':[{'model':'B_-1','target':'lambdaH1',
                'source':'lambdaH1', 'reason':'unresolved_native_symbolic_coefficient'}],
        }
        same={'native_coefficient_numeric':True,
              'native_coefficient_value':1.,'linear_in_source':True}
        reprobe={'format':'t3-parallel-a-native-rgbeta-quartic-response-reprobe-v13',
          'status':'NativeResponseForensicsCaptured_NoPhysicsAuthorization',
          'forensic_input_sha256':'a'*64,
          'physical_Wilson_evolution_authorized':False,
          'installed_RGBeta_main_file_sha256':'9'*64,
          'models':{'B_-1':{'probes':[{'target':'lambdaH1','source':'lambdaH1',
                        'isolated_registry':dict(same),
                        'full_registry_then_other_couplings_zero':dict(same)}]}}}
        return forensic,reprobe

    def test_both_routes_numeric_classified_not_certified(self):
        from scripts.RunT3ParallelACompareQuarticReprobeV13 import compare
        forensic,reprobe=self.fixture_pair()
        result=compare(forensic,reprobe,'a'*64)
        self.assertEqual(result['results'][0]['comparison'],
           'both_native_routes_numeric_check_old_v10_failure')
        self.assertFalse(result['native_150_diagonals_certified'])
        self.assertFalse(result['physical_Wilson_evolution_authorized'])

    def test_different_routes_and_abstract_evidence_classified(self):
        from scripts.RunT3ParallelACompareQuarticReprobeV13 import compare
        forensic,reprobe=self.fixture_pair()
        first=reprobe['models']['B_-1']['probes'][0]
        first['isolated_registry']['native_coefficient_numeric']=False
        r=compare(forensic,reprobe,'a'*64)
        self.assertEqual(r['results'][0]['comparison'],
                         'isolation_and_full_registry_differ')

    def test_wrong_probe_source_hash_or_unsafe_physics_rejected(self):
        from scripts.RunT3ParallelACompareQuarticReprobeV13 import compare
        forensic,reprobe=self.fixture_pair()
        from scripts.RunT3ParallelAQuarticForensicsV13 import ForensicBlocked as SourceBlock
        with self.assertRaisesRegex(SourceBlock,'match'):
            compare(forensic,reprobe,'b'*64)
        reprobe['physical_Wilson_evolution_authorized']=True
        with self.assertRaisesRegex(SourceBlock,'physics'):
            compare(forensic,reprobe,'a'*64)

    def test_unknown_probe_targets_rejected(self):
        from scripts.RunT3ParallelACompareQuarticReprobeV13 import compare
        forensic,reprobe=self.fixture_pair()
        reprobe['models']['B_-1']['probes'][0]['source']='lambdaH2'
        from scripts.RunT3ParallelAQuarticForensicsV13 import ForensicBlocked as SourceBlock
        with self.assertRaisesRegex(SourceBlock,'mismatch'):
            compare(forensic,reprobe,'a'*64)
