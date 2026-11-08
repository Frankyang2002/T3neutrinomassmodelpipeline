"""Non-circular 2026 one-loop cross-check of the Batch34 gaugeless beta.

Run: python -m unittest discover -s tests -p test_batch34_independent_master.py -v
"""
from __future__ import annotations
import unittest
import numpy as np

from RGE.running.intermediate.GeneralEFT1IndependentOneLoopReference import (
    CHANNELS,equation22_gaugeless,sum_reference)
from RGE.running.intermediate.GeneralEFT1SevenTermFlow import (
    Couplings,master_beta,total_beta)
from RGE.running.intermediate.GeneralEFT1CrossedWardCandidate import candidate_beta
from scripts.CheckBatch34IndependentMaster import diagnostics,fixture


class IndependentMasterTest(unittest.TestCase):
    def test_six_channels_independently_reproduce_eq22(self):
        for seed in (311,42,2026):
            c,ctx=fixture(seed=seed)
            ref=equation22_gaugeless(c,ctx.yukawa,ctx.quartic)
            repo=master_beta(c,ctx,corrected_crossed=True)
            for key in CHANNELS:
                np.testing.assert_allclose(ref[key],repo[key],rtol=2e-13,atol=2e-13,err_msg=key)

    def test_total_equals_six_term_beta_and_neither_seven_term_variant(self):
        c,ctx=fixture()
        ref=sum_reference(equation22_gaugeless(c,ctx.yukawa,ctx.quartic))
        source=master_beta(c,ctx,corrected_crossed=True)
        candidate=candidate_beta(c,ctx)
        np.testing.assert_allclose(ref,total_beta(source)-source['crossed_yukawa'],rtol=2e-13,atol=2e-13)
        self.assertGreater(np.max(np.abs(ref-total_beta(source))),.5)
        self.assertGreater(np.max(np.abs(ref-total_beta(candidate))),.5)

    def test_both_crossed_candidates_are_nonzero_and_unmatched(self):
        report=diagnostics()
        self.assertEqual(report['status'],'BLOCKED_CONFLICTING_ONE_LOOP_FORMULAE')
        self.assertGreater(report['source_four_exchange_crossed_norm'],1.)
        self.assertGreater(report['Cstar_candidate_crossed_norm'],1.)
        self.assertLess(report['non_crossed_max_absolute_discrepancy'],1e-13)
        self.assertFalse(report['operator_basis_or_loop_formula_physically_certified'])

    def test_reference_is_real_linear_not_complex_linear(self):
        c,ctx=fixture(71)
        def beta(z):
            return sum_reference(equation22_gaugeless(z,ctx.yukawa,ctx.quartic))
        a,b=.7,-1.3
        np.testing.assert_allclose(beta(a*c+b*1j*c),a*beta(c)+b*beta(1j*c),rtol=1e-12,atol=1e-12)
        self.assertGreater(np.linalg.norm(beta(1j*c)-1j*beta(c)),1e-3)

    def test_reference_rejects_invalid_tensor_symmetry_and_shapes(self):
        c,ctx=fixture()
        wrong=c.copy(); wrong[0,1,0,1]+=1.
        with self.assertRaisesRegex(ValueError,'separate'):
            equation22_gaugeless(wrong,ctx.yukawa,ctx.quartic)
        with self.assertRaisesRegex(ValueError,'dimension'):
            equation22_gaugeless(c,ctx.yukawa[:2],ctx.quartic)


if __name__=='__main__':
    unittest.main()
