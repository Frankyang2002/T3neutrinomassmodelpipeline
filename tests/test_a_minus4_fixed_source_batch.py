"""Read-only structural regression for the complete corrected RGBeta model file."""
from __future__ import annotations

from hashlib import sha256
import ast
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0,str(ROOT))

from scripts.RunAminus4FreshIndexProductionBatch import (
    FIXED_SOURCE_SHA, ORIGINAL_SOURCE_SHA, MODEL, NONQUARTICS, QUARTICS,
    check_model_source, compatible_native, read, OLD_RAW, CAUSAL_BOTH,
    SCALAR_COORDS,
)
from Numerical.running import Aminus4EFT1Couplings as eft


class PermanentAminus4Source(unittest.TestCase):
    def test_checked_complete_model_source_hash(self):
        report=check_model_source(ROOT)
        self.assertEqual(report['sha256'],FIXED_SOURCE_SHA)
        self.assertTrue(report['other_models_not_numerically_certified'])

    def test_only_two_Aminus4_adjoint_dummy_instances_are_fresh(self):
        text=(ROOT/MODEL).read_text(encoding='utf8')
        self.assertEqual(text.count('Module[{adj},'),2)
        # Ensure Module is INSIDE the pure function, not evaluated once
        # when the invariant registration function is defined.
        portal=text.split('T3RGBetaPortalAdjInvariant[3] :=',1)[1].split(';',1)[0]
        triplet=text.split('lambdaS2Adj,\n            {Bar @ S2, S2, Bar @ S2, S2}',1)[1].split('];',1)[0]
        for expression in (portal,triplet):
            self.assertIn('Module[{adj}',expression)
            self.assertIn('] &',expression)
            self.assertNotRegex(expression,r'tGen\[[^\n]*, A,')

    def test_complete_model_changed_only_two_verified_invariants(self):
        """Replace only the two edited source spans; require exact old-file recovery.

        SequenceMatcher's changed lines need not contain the unmodified
        function declaration or coupling name. Compare the *whole file* after
        reversing the two anchored edits instead.
        """
        old = ROOT / 'RGE/running/rgbeta/T3RGBetaModel.before_fresh_index_20261108.wl'
        self.assertTrue(old.is_file(), 'Make source backup BEFORE extracting ZIP')
        initial = old.read_text(encoding='utf8')
        fixed = (ROOT / MODEL).read_text(encoding='utf8')
        self.assertEqual(sha256(initial.encode()).hexdigest(), ORIGINAL_SOURCE_SHA)
        self.assertEqual(sha256(fixed.encode()).hexdigest(), FIXED_SOURCE_SHA)

        def span(body: str, start: str, end: str) -> str:
            self.assertEqual(body.count(start), 1, f'Expected one start anchor: {start}')
            i = body.index(start)
            j = body.index(end, i + len(start))
            self.assertGreater(j, i)
            return body[i:j]

        # Includes the changed portal body and its added comment, but not
        # the unchanged function that follows it.
        portal_begin = '(* Extracted generators for doublet triplet*)'
        portal_end = 'T3RGBetaAddRGClosedQuartics['
        old_portal = span(initial, portal_begin, portal_end)
        new_portal = span(fixed, portal_begin, portal_end)
        self.assertIn('T3RGBetaPortalAdjInvariant[3] :=', new_portal)
        self.assertIn('Module[{adj},', new_portal)
        self.assertIn('tGen[SU2L @ fund, adj, #1, #2]', new_portal)
        self.assertIn('tGen[SU2L @ S2, adj, #3, #4]', new_portal)
        self.assertNotEqual(new_portal, old_portal)

        # Use the exact four-field source declaration, and check the stop
        # anchor belongs to its AddQuartic registration.
        self_begin = ('            lambdaS2Adj,\n'
                      '            {Bar @ S2, S2, Bar @ S2, S2},')
        self_end = '\n        ];\n    ];'
        old_tail = initial[initial.index(self_begin):]
        new_tail = fixed[fixed.index(self_begin):]
        old_self = old_tail[:old_tail.index(self_end)]
        new_self = new_tail[:new_tail.index(self_end)]
        self.assertIn('Module[{adj},', new_self)
        self.assertIn('tGen[SU2L @ S2, adj, #1, #2]', new_self)
        self.assertIn('tGen[SU2L @ S2, adj, #3, #4] / 2', new_self)
        self.assertNotEqual(new_self, old_self)

        restored = fixed.replace(new_portal, old_portal, 1)
        restored = restored.replace(new_self, old_self, 1)
        self.assertEqual(restored, initial,
                         'Unexpected change outside the two certified invariant regions')

    def test_full_EFT1_parameter_inventory_unchanged(self):
        archived=read(ROOT/OLD_RAW)
        self.assertEqual(set(archived['report_betas']),eft.BETA_NAMES)
        self.assertEqual(set(archived['report_betas']),QUARTICS|NONQUARTICS)
        self.assertEqual(len(SCALAR_COORDS),10)
        compatible=compatible_native(archived)
        self.assertEqual(compatible['metadata']['IntegratedField'],'F')
        self.assertEqual(compatible['metadata']['ActiveBSMFields'],['S1','S2'])

    def test_prior_independent_native_both_witness_preserved(self):
        w=read(ROOT/CAUSAL_BOTH)
        self.assertEqual(w['mode'],'both')
        self.assertEqual(w['source_model_SHA256'],ORIGINAL_SOURCE_SHA)
        self.assertTrue(w['fresh_triplet_adjoint_dummy'])
        self.assertTrue(w['fresh_portal_adjoint_dummy'])
        self.assertEqual(set(w['report_betas']),QUARTICS)

    def test_Wilson_physical_gate_still_present(self):
        """Check the compiled error message and guard, not source-line layout."""
        source = (ROOT / 'RGE/running/intermediate/GeneralEFT1SevenTermFlow.py').read_text(
            encoding='utf8')
        module = ast.parse(source)
        messages = [
            node.args[0].value
            for node in ast.walk(module)
            if isinstance(node, ast.Call)
            and isinstance(node.func, ast.Name)
            and node.func.id == 'ValueError'
            and node.args
            and isinstance(node.args[0], ast.Constant)
            and isinstance(node.args[0].value, str)
        ]
        self.assertTrue(any(
            'Physical A(-4) input is blocked' in message and
            'crossed term loop convention is not independently certified' in message
            for message in messages), 'Physical Wilson fail-closed error was removed')
        # Guard is imposed on either boundary or couplings marked certified.
        self.assertIn('context.authority == "certified"', source)
        self.assertIn('coupling_payload.get("authority") == "certified"', source)



if __name__ == '__main__':
    unittest.main()
