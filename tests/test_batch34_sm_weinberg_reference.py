"""Independent SM Weinberg-RGE limiting-case regression for Batch 3--4.

Uses the gauge-invariant L,eC,H toy fixture from the separately installed
Batch34 WardGate.  Antusch et al., Phys. Lett. B519 (2001) 238, Eq. (21):
  16*pi² beta(kappa) = -3/2 [kappa Y_e†Y_e + (Y_e†Y_e)^T kappa]
                         + lambda_Antusch*kappa - 3*g2²*kappa + 2*T*kappa.
Their scalar-potential convention is V=(lambda_Antusch/4)(H†H)^2;
we use V=lambda_SM (H†H)^2 and lambda_Antusch=4*lambda_SM.
For one charged-lepton generation and no quark Yukawas, this reduces to
  16*pi² beta(kappa) = (-|ye|² + 4*lambda_SM - 3*g2²)*kappa.

This test is a necessary independent SM reduction, NOT a verification of
crossed-Yukawa anomalous dimensions in the general T3 theory: the crossed
term is zero in this SM Weinberg projection.
"""
from __future__ import annotations

import unittest
import numpy as np

from scripts.CheckBatch34GaugeCovariance import example, ward_c
from RGE.running.intermediate.GeneralEFT1CrossedWardCandidate import candidate_beta
from RGE.running.intermediate.GeneralEFT1SevenTermFlow import (
    Couplings, total_beta,
)

YE_ABS_SQUARED = 0.30**2 + 0.10**2
LAMBDA_SM = 0.20
G2 = 0.65
TARGET = -YE_ABS_SQUARED + 4*LAMBDA_SM - 3*G2**2


def projected_ratio(c, couplings):
    beta=total_beta(candidate_beta(c,couplings))
    return beta[0,0,2,2]/c[0,0,2,2], beta


class StandardModelReferenceTests(unittest.TestCase):
    def setUp(self):
        self.c, self.inputs = example()
        self.y=self.inputs.yukawa
        self.lam=self.inputs.quartic
        self.gauge=self.inputs.gauge
        self.zeros_y=np.zeros_like(self.y)
        self.zeros_lam=np.zeros_like(self.lam)

    def assert_weinberg(self, beta, coefficient):
        np.testing.assert_allclose(beta, coefficient*self.c,
                                   atol=5e-13, rtol=5e-12)

    def test_charged_lepton_yukawa_reproduces_minus_one_ye_squared(self):
        c=Couplings(self.y,self.zeros_lam,())
        ratio,beta=projected_ratio(self.c,c)
        self.assertAlmostEqual(ratio.real,-YE_ABS_SQUARED,places=12)
        self.assertAlmostEqual(ratio.imag,0,places=12)
        self.assert_weinberg(beta,-YE_ABS_SQUARED)
        # Yukawa phase does not enter the SM Weinberg RGE.
        y2=self.y*np.exp(0.57j)
        ratio2,beta2=projected_ratio(self.c,Couplings(y2,self.zeros_lam,()))
        np.testing.assert_allclose(beta2,beta,atol=5e-13,rtol=5e-12)
        self.assertAlmostEqual(ratio2.real,ratio.real,places=12)

    def test_scalar_quartic_reproduces_four_lambda_sm(self):
        ratio,beta=projected_ratio(self.c,Couplings(self.zeros_y,self.lam,()))
        self.assertAlmostEqual(ratio.real,4*LAMBDA_SM,places=12)
        self.assert_weinberg(beta,4*LAMBDA_SM)

    def test_gauge_su2_coefficient_and_hypercharge_cancellation(self):
        beta_g=total_beta(candidate_beta(self.c,
            Couplings(self.zeros_y,self.zeros_lam,self.gauge)))
        beta_su2=total_beta(candidate_beta(self.c,
            Couplings(self.zeros_y,self.zeros_lam,(self.gauge[0],))))
        np.testing.assert_allclose(beta_g,beta_su2,atol=5e-13,rtol=5e-12)
        self.assert_weinberg(beta_g,-3*G2**2)

    def test_full_sm_reference_and_gauge_ward(self):
        parts=candidate_beta(self.c,self.inputs)
        beta=total_beta(parts)
        ratio=beta[0,0,2,2]/self.c[0,0,2,2]
        self.assertAlmostEqual(ratio.real,TARGET,places=12)
        self.assert_weinberg(beta,TARGET)
        for _, pairs in self.gauge:
            for tf,ts in pairs:
                self.assertLess(np.max(np.abs(ward_c(beta,tf,ts))),5e-12)
        # The crossed term does not contribute in the SM Weinberg limit.
        self.assertLess(np.linalg.norm(parts['crossed_yukawa']),1e-12)
        self.assertLess(np.linalg.norm(parts['conjugate_coefficient_yukawa']),1e-12)


if __name__ == '__main__':
    unittest.main()
