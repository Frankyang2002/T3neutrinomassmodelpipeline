"""Offline structural contracts for the A(-4) v8 targeted Matchete continuation.

Static tests do not certify CG invariance, Mathematica execution, or physical C5.
"""
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "Lagrangian" / "RunT3AMinus4FullGeneralPhysicalScalar.wl"


class TestTargetedMatcheteContinuation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.src = SRC.read_text(encoding="utf-8")

    def test_no_disabling_invariant_validation(self):
        self.assertIn("TraceScan[t3ACGTraceVisitor", self.src)
        self.assertIn("HoldPattern[DefineCG[___]]", self.src)
        self.assertIn("TraceInternal -> True", self.src)
        self.assertNotIn("Unprotect[DefineCG]", self.src)
        self.assertNotIn("DefineCG::notInvariant = .", self.src)

    def test_trace_follows_real_full_match(self):
        self.assertIn('Environment["T3_A4_TRACE_DEFINECG"]', self.src)
        self.assertIn("A_one_loop_DefineCG_calls_DIAGNOSTIC.txt", self.src)
        self.assertIn("Match[expr, EFTOrder -> 5, LoopOrder -> loop]", self.src)
        self.assertIn('"selected-full-theory-one-loop"', self.src)

    def test_two_heavy_trace_native_api(self):
        self.assertIn('Environment["T3_A4_PROBE_TWO_HEAVY_TRACE"]', self.src)
        self.assertIn("CovariantLoop[t3AMatchInput, {NewScalar1, NewScalar2}", self.src)
        self.assertIn("EFTOrder -> 5, Mode -> Matching", self.src)

    def test_no_substitution_of_diagnostic_for_full_theory(self):
        self.assertIn('"full_one_loop_Match_replaced" -> False', self.src)
        self.assertIn('"physical_C5_certified" -> False', self.src)
        self.assertIn('"F_hard_inherited" -> False', self.src)
        self.assertIn('"EFT1_Wilson_running_included" -> False', self.src)
        self.assertIn('t3Stop["A(-4) raw Matchete Match failed at one loop: "', self.src)

    def test_nonzero_weinberg_extraction_only(self):
        self.assertIn('ExtractWeinbergCoefficient[t3ATraceLinear]', self.src)
        self.assertIn('"Weinberg_term_count"', self.src)
        self.assertIn('A_two_heavy_scalar_Weinberg_DIAGNOSTIC.wxf', self.src)
        self.assertIn('A_two_heavy_scalar_Weinberg_DIAGNOSTIC.txt', self.src)

    def test_no_single_flavour_fallback(self):
        self.assertIn('"fermion_generations_in_W" -> 3', self.src)
        self.assertNotIn('NF=1', self.src)
        self.assertNotIn('NF = 1', self.src)

if __name__ == "__main__":
    unittest.main()
