"""Five-monomial diagnostic regression; no physical RGE/coupling certification.

Run: python -m unittest discover -s tests -p test_a_minus4_quartic_five_monomial_batch.py -v
"""
from __future__ import annotations

import unittest
from fractions import Fraction

import numpy as np

from Numerical.running.Aminus4QuarticDiscrepancy import (
    COORDINATES, DELTA_INPUTFORM, EXPECTED_MONOMIALS,
    assert_expected_five,
    eval_polynomial, extract_monomials,
)


class FiveMonomialAlgebra(unittest.TestCase):
    def test_source_discrepancy_evaluates_expected_five_monomials(self):
        ix={x:i for i,x in enumerate(COORDINATES)}
        x=np.zeros(10)
        x[ix['lambdaS2Adj']]=.3
        x[ix['lambdaH2Adj']]=.4
        x[ix['lambdaT3.real']]=.2
        x[ix['lambdaT3.imag']]=-.1
        d=eval_polynomial(x,EXPECTED_MONOMIALS)
        self.assertAlmostEqual(d[ix['lambdaS2']],8*.3**2)
        self.assertAlmostEqual(d[ix['lambdaS2Adj']],14*.3**2+2*.4**2)
        self.assertAlmostEqual(d[ix['lambdaH']],2*.4**2)
        self.assertAlmostEqual(d[ix['lambdaH2']],4*.4**2)
        self.assertAlmostEqual(d[ix['lambdaH2Adj']],10*.3*.4-16/3*(.2**2+.1**2))
        self.assertAlmostEqual(d[ix['lambdaT3.real']],0)
        self.assertAlmostEqual(d[ix['lambdaT3.imag']],0)

    def test_polarization_finds_only_five_monomials_and_34_failing_probes(self):
        fake_reference=np.array([.17*i for i in range(10)],float)
        def probe(names):
            x=np.array([float(name in names) for name in COORDINATES])
            reference=fake_reference*(x@x)+.14*(sum(x)**2)
            native=reference+eval_polynomial(x,EXPECTED_MONOMIALS)
            records={n:{'reference_scalar_loop':float(r),
                        'native_rgbeta_scalar_loop':float(v)}
                     for n,r,v in zip(COORDINATES,reference,native)}
            return {'coupling_beta_comparison':records,
                    'scalar_loop_projection_closure_relative':4e-16,
                    'reference_and_native_scalar_loop_match':bool(
                        np.max(abs(native-reference))<1e-10)}
        result=extract_monomials(probe)
        self.assertEqual(result['probe_count'],55)
        self.assertEqual(result['failing_probes'],34)
        self.assertEqual(len(result['monomials']),5)
        self.assertLess(assert_expected_five(result['monomials']),3e-13)
        self.assertLess(result['max_probed_reconstruction_abs'],3e-13)

    def test_rejects_unexpected_new_discrepancy(self):
        altered={k:dict(v) for k,v in EXPECTED_MONOMIALS.items()}
        altered[('lambdaS1','lambdaS1')]={'lambdaH':Fraction(1)}
        with self.assertRaisesRegex(ValueError,'support changed'):
            assert_expected_five(altered)
        altered={k:dict(v) for k,v in EXPECTED_MONOMIALS.items()}
        altered[('lambdaS2Adj','lambdaS2Adj')]['lambdaS2']=Fraction(9)
        with self.assertRaisesRegex(ValueError,'coefficient changed'):
            assert_expected_five(altered)

    def test_full_complex_quartic_is_norm_squared_not_real_squared(self):
        try:
            from Numerical.core.RGBetaEvaluator import evaluate_inputform_expression
        except ModuleNotFoundError:
            self.skipTest('Full T3 project RGBeta InputForm evaluator not available in minimal build environment')
        for re,im in ((.3,.4),(.4,0),(-.1,.2),(.0,-.3)):
            env={'lambdaS2Adj':.2,'lambdaH2Adj':.3,'lambdaT3':re+1j*im}
            observed=evaluate_inputform_expression(DELTA_INPUTFORM['lambdaH2Adj'],env)
            expected=10*.2*.3-16/3*(re*re+im*im)
            self.assertAlmostEqual(complex(observed).real,expected,places=12)
            self.assertAlmostEqual(complex(observed).imag,0,places=12)

    def test_known_a_minus4_signed_weyl_evidence_remains_unmeasured(self):
        # This is a project-file regression, not a proof that a sign is +1/-1.
        from pathlib import Path
        import json
        root=Path(__file__).resolve().parents[1]
        path=root/'output/general_observable/a_minus4_boundary_batch1/declared_convention_boundary.json'
        if not path.is_file():self.skipTest('native A(-4) source witness not in this checkout')
        witness=json.loads(path.read_text(encoding='utf-8-sig'))
        self.assertIsNone(witness['native_matchete_to_declared_Weyl_factor'])
        self.assertFalse(witness['physical_certification'])
        self.assertEqual(len(witness['declared_rho_plus_one_independent_coordinates']),36)
        self.assertEqual(len(witness['rho_minus_one_fixed_native_y1_alternative_coordinates']),36)

if __name__=='__main__':unittest.main()
