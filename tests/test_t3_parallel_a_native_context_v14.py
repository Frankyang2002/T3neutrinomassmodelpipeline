"""v14 PackageScope native RGBeta quartic source guards and response tests.

Synthetic inputs here are ONLY schema fixtures. They cannot substitute for the
installed Wolfram witness; live native integration is executed on Windows.
"""
from __future__ import annotations
from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest

from Numerical.parallel_t3.wilson.native_context_quartic_v14 import (
    FORMAT, SMOKE, FULL, BLOCKED, SOURCE_FILES, API,
    V14Blocked, read, file_hash, audit_witness,
)
from Numerical.parallel_t3.wilson.native_quartic_projectors_v10 import expected_symbols
from test_t3_parallel_a_native_quartic_projectors_v10 import fixture_native

ROOT=Path(__file__).resolve().parents[1]


def fixture(root:Path,*,smoke:bool=True,blocked:bool=False):
    base=fixture_native()
    if smoke:base['models']={k:v for k,v in base['models'].items() if k in {'A_-4','B_-1','E_-2'}}
    base['format']=FORMAT
    base['status']=BLOCKED if blocked else (SMOKE if smoke else FULL)
    base['smoke_only']=smoke
    base['registered_model_count']=len(base['models'])
    base['native_package_scoped_api_verified']=True
    base['resolved_native_package_symbols']=deepcopy(API)
    base['native_complex_projector_28_channels_certified']=False
    base['EFT_order']='UV -> integrate out ALL fermions -> integrate out ALL scalars -> SMEFT'
    p=root/'rgbeta_entrypoint.m';p.write_text('RGBeta v1.2.0 fake fixture')
    base['installed_RGBeta_main_file']=str(p)
    base['installed_RGBeta_main_file_sha256']=file_hash(p)
    base['source_sha256']={}
    for i,name in enumerate(SOURCE_FILES):
        q=root/name;q.parent.mkdir(parents=True,exist_ok=True)
        q.write_text('Source fixture '+str(i))
        base['source_sha256'][name]=file_hash(q)
    for tag,model in base['models'].items():
        model['native_complex_conjugate_projector_probes']=[
            {'symbol':s,'native_complex_conjugate_response':{'ok':True,
                'direct_channel':1.,'conjugate_channel':0.,'unit_direct_zero_conjugate':True}}
            for s,real in expected_symbols(tag).items() if not real]
    if blocked:
        base['native_Lam_evaluated']=False
        base['native_quartic_projectors_evaluated']=False
        base['models']['B_-1']['native_real_self_projectors'][0]['native_response']={
            'ok':False, 'reason':'native_tensor_contraction_remains_abstract'}
        base['models']['B_-1']['native_real_self_projectors_all_unit']=False
        base['models']['B_-1']['ok']=False
    return base


class TestNativeContextV14(unittest.TestCase):
    def test_source_has_correct_package_scoped_implementations(self):
        wl=(ROOT/'scripts/RunT3ParallelANativeQuarticProjectorsV14.wl').read_text()
        self.assertIn('RGBeta`PackageScope`Lam[t14a,t14b,t14c,t14d]',wl)
        self.assertIn('RGBeta`PackageScope`RefineGroupStructures[',wl)
        self.assertIn('RGBeta`PackageScope`TStructure',wl)
        self.assertIn('DownValues[RGBeta`PackageScope`Lam]',wl)
        self.assertIn('native_tensor_contraction_remains_abstract',wl)
        self.assertIn('complex_conjugate_channels_not_linear',wl)
        self.assertNotIn(' *\n          Lam[t14a',wl)
        self.assertNotIn('FullSimplify[RefineGroupStructures[',wl)
        self.assertIn('"physical_Wilson_evolution_authorized"->False',wl)

    def test_3_model_synthetic_smoke_counts_and_complex_routes(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);native=fixture(root)
            out=audit_witness(native,root)
            self.assertEqual(out['source_census']['quartic_registered'],38)
            self.assertEqual(out['source_census']['real_diagonal_unit'],29)
            self.assertEqual(out['source_census']['portal_response_attempted'],21)
            self.assertEqual(out['complex_probe_census']['complex_attempted'],9)
            self.assertTrue(out['real_projector_conditionally_verified'])
            self.assertTrue(out['native_complex_exploratory_complete'])
            self.assertFalse(out['physical_Wilson_evolution_authorized'])

    def test_full_synthetic_all16_counts(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);native=fixture(root,smoke=False)
            out=audit_witness(native,root)
            self.assertEqual(out['source_census']['quartic_registered'],178)
            self.assertEqual(out['source_census']['real_diagonal_unit'],150)
            self.assertEqual(out['source_census']['portal_response_numeric'],110)
            self.assertEqual(out['complex_probe_census']['complex_attempted'],28)
            self.assertTrue(out['real_projector_conditionally_verified'])

    def test_old_inert_package_api_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);native=fixture(root)
            native['resolved_native_package_symbols']['Lam']=['Global`Lam']
            with self.assertRaisesRegex(V14Blocked,'PackageScope'):
                audit_witness(native,root)

    def test_old_source_hash_or_installed_package_change_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);native=fixture(root)
            (root/SOURCE_FILES[2]).write_text('tampered script')
            with self.assertRaisesRegex(V14Blocked,'source changed'):
                audit_witness(native,root)
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);native=fixture(root)
            (root/'rgbeta_entrypoint.m').write_text('different installed package')
            with self.assertRaisesRegex(V14Blocked,'entrypoint'):
                audit_witness(native,root)

    def test_native_failure_inspectable_but_not_certified(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);native=fixture(root,blocked=True)
            out=audit_witness(native,root,allow_blocked=True)
            self.assertEqual(out['status'],'Blocked_NativeRGBetaRealProjectors_PackageScopeProbe')
            self.assertEqual(out['source_census']['real_diagonal_unit'],28)
            self.assertFalse(out['real_projector_conditionally_verified'])
            with self.assertRaisesRegex(V14Blocked,'status is blocked'):
                audit_witness(native,root)

    def test_complex_failure_is_separate_from_real_projector_acceptance(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);native=fixture(root)
            row=native['models']['A_-4']['native_complex_conjugate_projector_probes'][0]
            row['native_complex_conjugate_response']={'ok':False,'reason':'symbolic_conjugate_channel'}
            out=audit_witness(native,root)
            self.assertTrue(out['real_projector_conditionally_verified'])
            self.assertFalse(out['native_complex_exploratory_complete'])
            self.assertEqual(out['complex_probe_census']['complex_direct_unit_conj_zero'],8)
            self.assertFalse(out['complex_28_certified'])

    def test_invalid_real_diagonal_or_portal_singularity_detected(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);native=fixture(root,blocked=True)
            native['models']['B_-1']['native_real_self_projectors'][0]['native_response']={'ok':True,'value':1.}
            native['models']['B_-1']['native_portal_response_blocks'][0]['native_response_matrix'][0][1]['value']=1.
            native['models']['B_-1']['native_portal_response_blocks'][0]['native_response_matrix'][1][0]['value']=1.
            out=audit_witness(native,root,allow_blocked=True)
            self.assertTrue(any(x['type']=='portal_response_invalid' for x in out['failures']))

    def test_unsafe_physics_claim_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);native=fixture(root)
            native['physical_Wilson_evolution_authorized']=True
            with self.assertRaisesRegex(V14Blocked,'unsupported physics'):
                audit_witness(native,root)

    def test_v15_complex_channel_9_and_28_source_guarded(self):
        from Numerical.parallel_t3.wilson.native_complex_channel_v15 import audit_from_native
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);native=fixture(root)
            path=root/'native_v14_smoke.json';path.write_text(json.dumps(native))
            report=audit_from_native(root,path)
            self.assertEqual(report['counts']['attempted'],9)
            self.assertEqual(report['counts']['unit_direct_zero_conjugate'],9)
            self.assertFalse(report['physical_Wilson_evolution_authorized'])
            self.assertTrue(report['native_complex_direct_conjugate_symbolic_verified'])
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);native=fixture(root,smoke=False)
            path=root/'native_v14_full.json';path.write_text(json.dumps(native))
            report=audit_from_native(root,path)
            self.assertEqual(report['counts']['attempted'],28)
            self.assertTrue(report['native_complex_direct_conjugate_symbolic_verified'])

    def test_v15_complex_failure_report_is_not_a_false_pass(self):
        from Numerical.parallel_t3.wilson.native_complex_channel_v15 import audit_from_native
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);native=fixture(root,smoke=False)
            native['models']['C_-1']['native_complex_conjugate_projector_probes'][0][
                'native_complex_conjugate_response']={'ok':False, 'reason':'unresolved_conjugate_channel'}
            path=root/'native_v14_full.json';path.write_text(json.dumps(native))
            report=audit_from_native(root,path)
            self.assertEqual(report['counts']['attempted'],28)
            self.assertEqual(report['counts']['unit_direct_zero_conjugate'],27)
            self.assertFalse(report['native_complex_direct_conjugate_symbolic_verified'])
            self.assertFalse(report['native_complex_action_physical_normalization_certified'])

    def test_v15_tampered_native_source_or_physics_claim_rejected(self):
        from Numerical.parallel_t3.wilson.native_complex_channel_v15 import audit_from_native
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);native=fixture(root)
            native['physical_Wilson_evolution_authorized']=True
            path=root/'native_v14_smoke.json';path.write_text(json.dumps(native))
            with self.assertRaisesRegex(V14Blocked,'unsupported physics'):
                audit_from_native(root,path)

    def test_live_native_full_if_successful_witness_available(self):
        from Numerical.parallel_t3.wilson.native_context_quartic_v14 import audit_from_paths
        root=ROOT
        base=root/'output/parallel_t3/wilson'
        n=base/'rgbeta_quartic_projectors_native_v14.json'
        if not n.is_file():
            self.skipTest('Native all16 v14 export not yet executed')
        native=read(n)
        if native.get('status') != FULL:
            self.skipTest('Native all16 v14 witness exists but remains blocked')
        v9=base/'portal_real_action_conditional_v9.json'
        if not v9.is_file():
            self.skipTest('Native v14 witness present but v9 input unavailable')
        out=audit_from_paths(root,n,v9)
        self.assertEqual(out['source_census']['real_diagonal_unit'],150)
        self.assertEqual(out['source_census']['portal_response_numeric'],110)
        self.assertFalse(out['physical_Wilson_evolution_authorized'])

    def test_immutable_cli_and_runner_safety(self):
        ps=(ROOT/'scripts/RunT3ParallelAQuarticNativeContextV14.ps1').read_text()
        self.assertIn('finally {',ps)
        self.assertIn('RunT3ParallelAQuarticNativeContextAuditV14',ps)
        self.assertIn('NativeRGBetaQuarticProjectorSmoke3ModelsPassed',ps)
        self.assertIn('RunT3ParallelANativeQuarticProjectorsV14.wl',ps)
        self.assertNotIn('Remove-Item',ps)
        from scripts.RunT3ParallelAQuarticNativeContextAuditV14 import write_immutable
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'witness.json';write_immutable(p,{'x':1})
            write_immutable(p,{'x':1})
            with self.assertRaisesRegex(V14Blocked,'overwrite'):
                write_immutable(p,{'x':2})

if __name__=='__main__':unittest.main()
