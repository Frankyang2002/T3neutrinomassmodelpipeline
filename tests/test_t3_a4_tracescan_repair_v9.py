"""Offline source-contract tests for v9 TraceScan repair; not physics validation."""
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'Lagrangian/RunT3AMinus4FullGeneralPhysicalScalar.wl'


class TraceScanRepairTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = SOURCE.read_text(encoding='utf-8')

    def test_unsupported_tracescan_option_removed(self):
        # TraceOriginal is an option of Trace, not TraceScan.
        self.assertNotIn('TraceOriginal -> True', self.source)
        self.assertIn('TraceInternal -> True]', self.source)
        self.assertIn('TraceScan[t3ACGTraceVisitor,', self.source)

    def test_full_match_still_runs_unmodified_loop_order_one(self):
        self.assertIn('Match[expr, EFTOrder -> 5, LoopOrder -> loop]', self.source)
        self.assertIn('"selected-full-theory-one-loop"', self.source)
        self.assertIn('"selected-full-theory-tree-control"', self.source)
        self.assertIn('t3AOneLoopAttempt["report"]["structurally_valid_for_requested_order"]', self.source)

    def test_abort_and_invalid_instrumentation_are_separate(self):
        for s in ('$Aborted', '$TimedOut', '$Failed',
                  'UnevaluatedTraceScan_InvalidInstrumentation',
                  'DefineCG_trace_count'):
            self.assertIn(s, self.source)

    def test_both_trace_streams_present(self):
        self.assertIn('A_one_loop_DefineCG_calls_DIAGNOSTIC.txt', self.source)
        self.assertIn('A_two_heavy_DefineCG_calls_DIAGNOSTIC.txt', self.source)
        self.assertIn('Close[t3ACGTraceStream]', self.source)
        self.assertIn('Flush[t3ACGTraceStream]', self.source)
        self.assertIn('TraceScan[t3ACGTraceVisitor,\n            CovariantLoop', self.source)

    def test_both_scalars_in_targeted_native_loop(self):
        self.assertIn('CovariantLoop[t3AMatchInput, {NewScalar1, NewScalar2}', self.source)
        self.assertIn('EFTOrder -> 5, Mode -> Matching', self.source)
        self.assertIn('"TwoHeavyScalarCovariantLoopOnly"', self.source)

    def test_three_flavour_physics_and_authority_gates(self):
        for s in ('"fermion_generations_in_W" -> 3',
                  '"F_hard_inherited" -> False',
                  '"physical_C5_certified" -> False',
                  '"full_one_loop_Match_replaced" -> False',
                  'CheckLagrangian',
                  'literal_F_boundary_equivalence_certified',
                  't3ASelectedExploratory'):
            self.assertIn(s, self.source)

    def test_failed_native_match_still_blocks(self):
        self.assertIn('t3Stop["A(-4) raw Matchete Match failed at one loop:', self.source)
        self.assertIn('A_joint_scalar_match_diagnostics.json', self.source)

if __name__ == '__main__':
    unittest.main()
