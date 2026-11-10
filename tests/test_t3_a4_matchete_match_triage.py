"""Offline regression tests for v7 Match-failure classification.

These tests do not execute Matchete, establish EFT equivalence, or certify C5.
"""
from pathlib import Path
import unittest

SRC = Path(__file__).resolve().parents[1] / "Lagrangian" / "RunT3AMinus4FullGeneralPhysicalScalar.wl"


class TestA4MatchTriage(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = SRC.read_text(encoding="utf-8")

    def test_exact_failure_statuses_recorded(self):
        for item in ('$Failed', '$Aborted', '$TimedOut', 'OtherAtomicResult',
                     'return_InputForm_if_atomic', 'return_head'):
            self.assertIn(item, self.source)

    def test_tree_control_before_one_loop(self):
        text = self.source
        tree = text.index('"selected-full-theory-tree-control"')
        loop = text.index('"selected-full-theory-one-loop"')
        self.assertLess(tree, loop)
        self.assertIn('Match[expr, EFTOrder -> 5, LoopOrder -> loop]', text)

    def test_json_checkpoint_before_failure(self):
        text = self.source
        self.assertIn('A_joint_scalar_match_diagnostics.json', text)
        self.assertIn('t3ASaveMatchReport[];', text)
        self.assertIn('"one_loop_return_inputform"', text)

    def test_same_validated_input_used_in_both_orders(self):
        text = self.source
        self.assertIn('t3AMatchAttempt[t3AMatchInput, 0,', text)
        self.assertIn('t3AMatchAttempt[t3AMatchInput, 1,', text)
        self.assertIn('"scalars_integrated_together" -> {"S1", "S2"}', text)
        self.assertIn('"fermion_generations_in_W" -> 3', text)

    def test_wilson_free_probe_not_accepted_as_physics(self):
        text = self.source
        self.assertIn('Environment["T3_A4_PROBE_WILSON_FREE_MATCH"] === "1"', text)
        self.assertIn('"wilson_free_probe_used_for_physics" -> False', text)
        self.assertIn('t3AWithoutW = t3AMatchAttempt[t3AScalarBase, 1,', text)
        self.assertNotIn('t3AMatched = t3AWithoutW', text)

    def test_finite_c5_not_certified(self):
        text = self.source
        self.assertIn('"physical_C5_certified" -> False', text)
        self.assertIn('"F_hard_inherited" -> False', text)
        self.assertIn('"EFT1_Wilson_running_included" -> False', text)
        self.assertIn('t3ASelectedExploratory', text)


def bracket_balance_mathematica(src: str):
    """Checks brackets ignoring strings and nested Mathematica comments."""
    pairs = {')': '(', ']': '[', '}': '{'}
    stack = []
    i, comment_depth, inside_string = 0, 0, False
    while i < len(src):
        if comment_depth:
            if src.startswith('(*', i):
                comment_depth += 1; i += 2; continue
            if src.startswith('*)', i):
                comment_depth -= 1; i += 2; continue
            i += 1; continue
        if inside_string:
            if src[i] == '\\':
                i += 2; continue
            if src[i] == '"':
                inside_string = False
            i += 1; continue
        if src.startswith('(*', i):
            comment_depth = 1; i += 2; continue
        if src[i] == '"':
            inside_string = True; i += 1; continue
        if src[i] in '([{':
            stack.append(src[i])
        elif src[i] in ')]}':
            if not stack or stack.pop() != pairs[src[i]]:
                return False
        i += 1
    return not stack and not comment_depth and not inside_string


class TestMathematicaStructure(unittest.TestCase):
    def test_bracket_balance(self):
        self.assertTrue(bracket_balance_mathematica(SRC.read_text(encoding="utf-8")))


if __name__ == "__main__":
    unittest.main()
