"""Offline tests for X-heavy free-Lagrangian rebase diagnostic and fail-closed matching gate.

These are source/driver tests; they do not execute Matchete and cannot certify
sequential Wilson matching.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace
import tempfile
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
WL=ROOT/'Lagrangian/RunT3BCXYNativeFollowup.experimental.wl'
PYS=ROOT/'scripts/RunT3BCXYNativeFollowup.py'
spec=importlib.util.spec_from_file_location('xy_rebase_stage_v2',PYS)
stage=importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(stage)

FLAGS={x:True for x in stage.X_FLAGS}


class TestMatcheteCheckSelection(unittest.TestCase):
    def test_correct_full_six_flags(self):
        self.assertTrue(stage.full_native_check(FLAGS))
        for label in FLAGS:
            with self.subTest(label=label):
                bad=dict(FLAGS,**{label:False})
                self.assertFalse(stage.full_native_check(bad))
                missing=dict(FLAGS)
                missing.pop(label)
                self.assertFalse(stage.full_native_check(missing))
        self.assertFalse(stage.full_native_check({**FLAGS,'other':True}))
        self.assertFalse(stage.full_native_check({**FLAGS,'canonical':1}))

    def test_rebased_valid(self):
        report={'matching_input_selection':'heavy_free_rebased_y_tree_eomeft',
                'selected_x_input_full_matchete_check_passed':True,
                'rebased_x_heavy_eft2_detailed_check':dict(FLAGS),
                'free_lagrangians_algebraically_identical':False}
        self.assertTrue(stage.certified_x_input_selection(report))
        report['rebased_x_heavy_eft2_detailed_check']['hermiticity']=False
        self.assertFalse(stage.certified_x_input_selection(report))

    def test_unchanged_direct_requires_free_equality(self):
        report={'matching_input_selection':'unchanged_y_tree_eomeft_proven_equal_free_lagrangians',
                'selected_x_input_full_matchete_check_passed':True,
                'direct_x_heavy_eft2_detailed_check':dict(FLAGS),
                'rebased_x_heavy_eft2_detailed_check':{'ordinary_check':False},
                'free_lagrangians_algebraically_identical':False}
        self.assertFalse(stage.certified_x_input_selection(report))
        report['free_lagrangians_algebraically_identical']=True
        self.assertTrue(stage.certified_x_input_selection(report))
        report['direct_x_heavy_eft2_detailed_check']['heavy_mass_basis']=False
        self.assertFalse(stage.certified_x_input_selection(report))

    def test_no_unqualified_fallback(self):
        report={'matching_input_selection':'matched_eft',
                'selected_x_input_full_matchete_check_passed':True,
                'direct_x_heavy_eft2_detailed_check':dict(FLAGS),
                'rebased_x_heavy_eft2_detailed_check':dict(FLAGS),
                'free_lagrangians_algebraically_identical':True}
        self.assertFalse(stage.certified_x_input_selection(report))
        report['matching_input_selection']='heavy_free_rebased_y_tree_eomeft'
        del report['selected_x_input_full_matchete_check_passed']
        self.assertFalse(stage.certified_x_input_selection(report))


class TestGeneratedMatchete(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.wl=WL.read_text(encoding='utf8')

    def test_full_source_retained(self):
        self.assertGreater(len(self.wl),110000)
        self.assertIn('T3MEVerifiedRenormalizableReconstruction',self.wl)
        self.assertIn('nativeExtra=T3MEBuildNativeExtraQuartics',self.wl)
        self.assertIn('rawTree=yTree["EOMEFT"]',self.wl)
        self.assertIn('"y_loop_inherited_in_x_matching"->False',self.wl)

    def test_two_candidate_checks_and_delta(self):
        for token in ('directHeavyCheck=T3BCXFollowupCheck[rawTree]',
                      'check=T3BCXFollowupCheck[rebased]',
                      'freeDelta=CheckAbort[Expand[freeHeavy-freeLight]',
                      'freeDeltaZero=TrueQ[freeDelta===0]',
                      '"direct_heavy_checks"->directHeavyCheck',
                      '"rebased_heavy_checks"->check',
                      '"free_lagrangians_algebraically_identical"->freeDeltaZero',
                      '"matched_eft_fallback_used"->False'):
            with self.subTest(token=token):self.assertIn(token,self.wl)

    def test_heavy_mass_reclassification(self):
        self.assertIn('RemoveField[NewScalarX]',self.wl)
        self.assertIn('Mass->{Heavy,MSX}',self.wl)
        self.assertIn('RunT3Matching[matchingInput,5,0]',self.wl)
        self.assertIn('RunT3Matching[matchingInput,5,1]',self.wl)
        self.assertNotIn('RunT3Matching[yTree["MatchedEFT"]',self.wl)
        self.assertIn('If[directHeavyOK && freeDeltaZero,',self.wl)

    def test_brackets_mathematica_basic(self):
        src=self.wl
        stack=[];i=0;state='normal';depth=0
        while i<len(src):
            t=src[i:i+2];c=src[i]
            if state=='string':
                if c=='\\':i+=2;continue
                if c=='"':state='normal'
            elif depth:
                if t=='(*':depth+=1;i+=2;continue
                if t=='*)':depth-=1;i+=2;continue
            else:
                if t=='(*':depth=1;i+=2;continue
                if c=='"':state='string'
                elif c in '[{(':stack.append(c)
                elif c in ']})':
                    self.assertTrue(stack,f'extra {c} at {i}')
                    self.assertEqual({'[':']','{':'}','(':')'}[stack.pop()],c)
            i+=1
        self.assertEqual(depth,0)
        self.assertEqual(state,'normal')
        self.assertEqual(stack,[])

    def test_physical_flags_false(self):
        for token in ('"saved_fermion_EOMEFT_equivalence_certified"->False',
                      '"running_threshold_couplings_applied"->False',
                      '"physical_matching_certified"->False',
                      '"two_threshold_physical_prediction"->False'):
            self.assertIn(token,self.wl)


class TestDiagnosticsDriver(unittest.TestCase):
    def test_diagnose_records_failure_but_does_not_certify(self):
        with tempfile.TemporaryDirectory() as td:
            out=Path(td)
            def dummy_process(command,**kwargs):
                target=Path(command[-1]);target.parent.mkdir(parents=True,exist_ok=True)
                stem=target.with_suffix('')
                target.write_text(json.dumps({'status':stage.Y_STATUS,'physical_matching_performed':False}))
                Path(str(stem)+'_x_followup.json').write_text(json.dumps({
                    'status':'XFollowupRebasedEFT2CheckFailed',
                    'direct_x_heavy_eft2_detailed_check':dict(FLAGS),
                    'rebased_x_heavy_eft2_detailed_check':{**FLAGS,'hermiticity':False},
                    'free_lagrangians_algebraically_identical':False,
                    'physical_matching_certified':False}))
                Path(str(stem)+'_x_followup_rebase_crosschecks.json').write_text(json.dumps({
                    'direct_heavy_check_passed':True,
                    'rebased_heavy_check_passed':False,
                    'free_lagrangians_algebraically_identical':False}))
                return SimpleNamespace(returncode=0,stdout='Matchete diagnostic',stderr='')
            fakepaths={'parameters':Path(td)/'params.json','first_fermion_eomeft':Path(td)/'source.txt'}
            with patch.object(stage,'preflight',return_value={'script':'hash'}), \
                 patch.object(stage,'input_paths',return_value=fakepaths), \
                 patch.object(stage.subprocess,'run',side_effect=dummy_process):
                result=stage.run_model('B_-1',mode='diagnose',output_dir=out)
                self.assertFalse(result['all_x_checks_passed'])
                self.assertEqual(result['status'],'XFollowupRebasedEFT2CheckFailed')
                self.assertFalse(result['physical_c5_certified'])
                self.assertTrue((out/'B_m1_final_xy_native_diagnose_x.provenance.json').exists())
                with self.assertRaisesRegex(RuntimeError,'rebased checks'):
                    stage.run_model('B_-1',mode='check',output_dir=out)

if __name__=='__main__':unittest.main()
