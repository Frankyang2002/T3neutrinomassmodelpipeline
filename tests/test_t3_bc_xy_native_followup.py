"""Offline structural and contract tests for experimental XY native followup."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from scripts import RunT3BCXYNativeFollowup as stage

WL=Path(__file__).resolve().parents[1]/'Lagrangian/RunT3BCXYNativeFollowup.experimental.wl'

class TestXYFollowup(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source=WL.read_text(encoding='utf-8')

    def test_generated_wl_exists_and_size(self):
        self.assertGreater(len(self.source),100000)

    def test_expected_source_invariants(self):
        for s in ('T3BCXFollowup[treeResult,loopResult,outPath]',
                  'XFollowupRebasedEFT2CheckFailed',
                  'XFollowupYTreeOnlyXRawMatchingNotCertified',
                  'T3_BC_X_FOLLOWUP_MATCH'):
            with self.subTest(s=s):self.assertIn(s,self.source)

    def test_only_y_tree_eomeft_used_for_x_match(self):
        self.assertIn('rawTree=yTree["EOMEFT"]',self.source)
        self.assertIn('rebased=CheckAbort[Expand[rawTree-freeLight+freeHeavy]',self.source)
        self.assertIn('"y_loop_inherited_in_x_matching"->False',self.source)

    def test_no_matched_eft_fallback(self):
        self.assertNotIn('rebased=yTree["MatchedEFT"]',self.source)
        self.assertIn('!KeyExistsQ[yTree,"EOMEFT"]',self.source)

    def test_full_hermiticity_and_heavy_basis_gate(self):
        for token in ('"Hermiticity"','"CanonicallyNormalized"',
                      '"HeavyMassBasis"','"GaugeInvariant"','"AllObjectsDefined"'):
            self.assertIn(token,self.source)
        self.assertIn('If[!rebasedOK,',self.source)
        self.assertIn('If[directHeavyOK && freeDeltaZero,',self.source)

    def test_delimiter_balance_ignoring_comments_strings(self):
        # Basic lexer: no claims of genuine Mathematica parsing or evaluation.
        s=self.source; stack=[]; i=0; mode='normal'; comment=0
        while i<len(s):
            two=s[i:i+2];ch=s[i]
            if mode=='string':
                if ch=='\\':i+=2;continue
                if ch=='"':mode='normal'
            elif comment:
                if two=='(*':comment+=1;i+=2;continue
                if two=='*)':comment-=1;i+=2;continue
            else:
                if two=='(*':comment=1;i+=2;continue
                if ch=='"':mode='string'
                elif ch in '[{(' :stack.append((ch,i))
                elif ch in ']})':
                    self.assertTrue(stack,f'unmatched {ch} at {i}')
                    a,_=stack.pop();self.assertEqual({'[':']','{':'}','(':')'}[a],ch)
            i+=1
        self.assertEqual(comment,0)
        self.assertEqual(mode,'normal')
        self.assertEqual(stack,[])

    def test_invalid_model(self):
        with self.assertRaises(RuntimeError):stage.input_paths('A_+0')

    def test_status_constants(self):
        self.assertEqual(stage.CHECK_STATUS,'XFollowupEFT2XHeavyCheckPassedOnly')
        self.assertEqual(stage.MATCH_STATUS,'XFollowupYTreeOnlyXRawMatchingNotCertified')

    def test_paths(self):
        for model,stem in [('B_-1','B_m1'),('C_-1','C_m1')]:
            self.assertIn(stem,stage.input_paths(model)['parameters'].name)
            self.assertIn(stem,stage.input_paths(model)['first_fermion_eomeft'].name)

    def test_missing_inputs_fail(self):
        with patch.object(stage,'WL',Path('nonexistent.wl')):
            with self.assertRaises(RuntimeError):stage.preflight('B_-1')

    def test_flags_not_promoted(self):
        for x in ('"physical_matching_certified"->False',
                  '"two_threshold_physical_prediction"->False',
                  '"running_threshold_couplings_applied"->False'):
            self.assertIn(x,self.source)

    def test_x_mass_reclassified(self):
        self.assertIn('RemoveField[NewScalarX]',self.source)
        self.assertIn('Mass->{Heavy,MSX}',self.source)

if __name__=='__main__':unittest.main()

class TestEvidenceCollector(unittest.TestCase):
    def test_no_prior_evidence_fails_closed(self):
        from scripts import CollectT3BCXYNativeEvidence as collector
        with tempfile.TemporaryDirectory() as d:
            with patch.object(collector,'ROOT',Path(d)),\
                 patch.object(collector,'DIR',Path(d)/'matchete'),\
                 patch.object(collector,'YDIR',Path(d)/'matchete'/'y'),\
                 patch.object(collector,'XDIR',Path(d)/'matchete'/'x'):
                with self.assertRaises(RuntimeError):
                    collector.collect(('B_-1',),Path(d)/'evidence.zip')

    def test_collector_flags_uncertified(self):
        from scripts import CollectT3BCXYNativeEvidence as collector
        import inspect
        src=inspect.getsource(collector.collect)
        self.assertIn("'physical_sequential_matching_certified':False",src)
        self.assertIn("'source_equivalence_certified':False",src)
