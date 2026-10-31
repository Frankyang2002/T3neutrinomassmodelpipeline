"""Offline source and failure-path tests for the experimental fresh-kernel X probe."""
from __future__ import annotations
from importlib.util import spec_from_file_location,module_from_spec
import json
from pathlib import Path
from types import SimpleNamespace
import tempfile
import unittest
from unittest.mock import patch
from zipfile import ZipFile

ROOT=Path(__file__).resolve().parents[1]

def load(name:str,path:Path):
    spec=spec_from_file_location(name,path)
    mod=module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

stage=load('fresh_stage',ROOT/'scripts/RunT3BCXFreshKernel.py')
collector=load('fresh_collector',ROOT/'scripts/CollectT3BCXFreshEvidence.py')
FLAGS={key:True for key in stage.FLAGS}

class TestSource(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.wl=(ROOT/'Lagrangian/RunT3BCXFromSavedY.experimental.wl').read_text(encoding='utf8')
    def test_native_full_source_retained(self):
        self.assertGreater(len(self.wl),115000)
        for token in ('T3MEVerifiedRenormalizableReconstruction',
                      'T3MEBuildNativeExtraQuartics',
                      'T3BCXFollowup[<|"EOMEFT"->parsedTree|>',
                      'If[Environment["T3_BC_X_FRESH_EOM_TREE"] =!= "",',
                      'T3MECheckpoint["x-fresh-saved-y-eomeft-parsed"]',
                      'Return[Module[{pathTree',
                      '"fresh_kernel_without_y_matching"->True',
                      '"y_one_loop_correction_inherited_in_x_source"->False',
                      '"saved_fermion_eomeft_equivalence_certified"->False'):
            self.assertIn(token,self.wl)
    def test_snapshot_immediately_before_tree_and_loop(self):
        a=self.wl.index('"XFollowupTreeMatchingInProgress"')
        b=self.wl.index('xTree=CheckAbort[RunT3Matching[matchingInput,5,0]')
        self.assertLess(a,b)
        self.assertIn('Export[stem<>".json",report,"RawJSON"]',self.wl[a:b])
        a=self.wl.index('"XFollowupOneLoopMatchingInProgress"')
        b=self.wl.index('xLoop=CheckAbort[RunT3Matching[matchingInput,5,1]')
        self.assertLess(a,b)
        self.assertIn('Export[stem<>".json",report,"RawJSON"]',self.wl[a:b])
    def test_no_silent_matched_eft_substitution(self):
        self.assertIn('matchingInput=rebased;',self.wl)
        self.assertNotIn('matchingInput=yTree["MatchedEFT"]',self.wl)
        self.assertIn('freeDeltaZero=TrueQ[freeDelta===0]',self.wl)
        self.assertIn('"SavedYTreeEOMConsumedForFreshXNoYRematch"',self.wl)
    def test_delimiters(self):
        stack=[];mode='normal';depth=0;i=0;s=self.wl
        while i<len(s):
            pair=s[i:i+2];c=s[i]
            if mode=='string':
                if c=='\\':i+=2;continue
                if c=='"':mode='normal'
            elif depth:
                if pair=='(*':depth+=1;i+=2;continue
                if pair=='*)':depth-=1;i+=2;continue
            else:
                if pair=='(*':depth=1;i+=2;continue
                if c=='"':mode='string'
                elif c in '[{(':stack.append(c)
                elif c in ']})':
                    self.assertTrue(stack, f'extra {c} at {i}')
                    self.assertEqual({'[':']','{':'}','(':')'}[stack.pop()],c)
            i+=1
        self.assertEqual(depth,0);self.assertEqual(stack,[]);self.assertEqual(mode,'normal')

class TestChecks(unittest.TestCase):
    def test_six_conditions(self):
        self.assertTrue(stage.passed(FLAGS))
        for key in FLAGS:
            with self.subTest(key=key):
                self.assertFalse(stage.passed({**FLAGS,key:False}))
        self.assertFalse(stage.passed({**FLAGS,'new_flag':True}))
        self.assertFalse(stage.passed({**FLAGS,'canonical':1}))
        self.assertFalse(stage.passed({}))
    def test_classification_exact(self):
        y={'status':stage.FRESH_Y_STATUS}
        x={'status':stage.MATCH_STATUS,
           'rebased_x_heavy_eft2_detailed_check':dict(FLAGS),
           'matching_input_selection':'heavy_free_rebased_y_tree_eomeft',
           'physical_matching_certified':False}
        self.assertEqual(stage.classify(y,x,0),(stage.MATCH_STATUS,True))
        self.assertFalse(stage.classify(y,x,1)[1])
        self.assertFalse(stage.classify(y,{**x,'physical_matching_certified':True},0)[1])
        self.assertFalse(stage.classify(y,{**x,'matching_input_selection':'matched_eft'},0)[1])
        self.assertFalse(stage.classify({},x,0)[1])
        self.assertFalse(stage.classify(y,{**x,'rebased_x_heavy_eft2_detailed_check':{}},0)[1])
    def test_driver_handles_interrupted_stage(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            def abort_stub(args,**kwargs):
                target=Path(args[-1]);stem=target.with_suffix('')
                target.write_text(json.dumps({'status':stage.FRESH_Y_STATUS}))
                (Path(str(stem)+'_x_followup.json')).write_text(json.dumps({
                    'status':'XFollowupTreeMatchingInProgress',
                    'rebased_x_heavy_eft2_detailed_check':dict(FLAGS),
                    'matching_input_selection':'heavy_free_rebased_y_tree_eomeft'}))
                (Path(str(stem)+'.progress.json')).write_text(json.dumps({
                    'last_stage':'x-followup-tree-matching-started',
                    'stage_history':['x-followup-tree-matching-started']}))
                return SimpleNamespace(returncode=5,stdout='tree attempted',stderr='kernel aborted')
            fake={k:root/k for k in ('parameters','fermion_eom','prior_y_tree','prior_y_full','fresh_wl')}
            with patch.object(stage,'preflight',return_value={'model':'B_-1'}),\
                 patch.object(stage,'model_paths',return_value=fake),\
                 patch.object(stage.subprocess,'run',side_effect=abort_stub):
                res=stage.run_one('B_-1',mode='match',out=root/'out')
                self.assertEqual(res['status'],'XFollowupTreeMatchingInProgress')
                self.assertFalse(res['fresh_x_matchete_matching_completed'])
                self.assertEqual(res['checkpoint'],'x-followup-tree-matching-started')
                self.assertFalse(res['physical_matching_certified'])
    def test_driver_saves_success_status(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            def success_stub(args,**kwargs):
                target=Path(args[-1]);stem=target.with_suffix('')
                target.write_text(json.dumps({'status':stage.FRESH_Y_STATUS}))
                (Path(str(stem)+'_x_followup.json')).write_text(json.dumps({
                    'status':stage.MATCH_STATUS,'physical_matching_certified':False,
                    'rebased_x_heavy_eft2_detailed_check':dict(FLAGS),
                    'matching_input_selection':'heavy_free_rebased_y_tree_eomeft'}))
                return SimpleNamespace(returncode=0,stdout='complete',stderr='')
            fake={k:root/k for k in ('parameters','fermion_eom','prior_y_tree','prior_y_full','fresh_wl')}
            with patch.object(stage,'preflight',return_value={'model':'B_-1'}),\
                 patch.object(stage,'model_paths',return_value=fake),\
                 patch.object(stage.subprocess,'run',side_effect=success_stub):
                res=stage.run_one('B_-1',mode='match',out=root/'out')
                self.assertTrue(res['fresh_x_matchete_matching_completed'])
                self.assertFalse(res['physical_matching_certified'])
                self.assertTrue((root/'out/B_m1_final_x_from_saved_y_match.provenance.json').is_file())

class TestEvidence(unittest.TestCase):
    def test_partials_allowed(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            fake_mat=root/'mat';fake_mat.mkdir()
            (fake_mat/'B_m1_final_parameters.json').write_text('{}')
            with patch.object(collector,'ROOT',root),\
                 patch.object(collector,'MAT',fake_mat),\
                 patch.object(collector,'V2',root/'v2'),\
                 patch.object(collector,'V3',root/'v3'),\
                 patch.object(collector,'Y',root/'y'):
                result=collector.collect(('B_-1','C_-1'),root/'out.zip')
                self.assertFalse(result['physical_sequential_matching_certified'])
                self.assertGreater(result['missing_reports'],0)
                with ZipFile(root/'out.zip') as z:
                    obj=json.loads(z.read('MANIFEST.json'))
                    self.assertFalse(obj['physical_sequential_matching_certified'])
                    self.assertEqual(obj['models'],['B_-1','C_-1'])
    def test_destination_not_overwritten(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);(root/'params').write_text('x')
            with patch.object(collector,'ROOT',root),\
                 patch.object(collector,'MAT',root),\
                 patch.object(collector,'V2',root),\
                 patch.object(collector,'V3',root),\
                 patch.object(collector,'Y',root):
                with self.assertRaisesRegex(ValueError,'invalid model'):
                    collector.collect(('UNKNOWN',),root/'o.zip')

if __name__=='__main__':unittest.main()
