"""Joint extended EFT1 analytic integration tests (synthetic tensors only)."""
import math
import unittest
import numpy as np
from scripts.TestT3QuarticWilsonBridge import fixture
from Numerical.running.T3ExtendedSimultaneousEFT1 import (
    JointEFT1Adapter, integrate_extended_eft1)
from Numerical.running.T3RunningQuarticTensor import QuarticEvidenceError


def setup(evolving=True):
    ctx,c,basis=fixture()
    mu_f,mu_s=1e5,5e4
    def unpack(ren,mu):
        state={'lambdaS1':float(ren[0]),'lambdaH1':0.3}
        h=basis.evaluate(state)[:4,:4,4:,4:]
        return state,h
    def derivative(t,ren):
        return np.array([0.07 if evolving else 0.0])
    def context(running,hhss,mu):
        return ctx,None
    adapter=JointEFT1Adapter('A_-4',mu_f,mu_s,unpack,derivative,context)
    return adapter,np.array([0.2]),c,basis


class TestFullRunningIVP(unittest.TestCase):
    def test_extended_eft1_runs_one_ivp(self):
        args=setup()
        got=integrate_extended_eft1(*args,rtol=1e-9,atol=1e-11)
        self.assertTrue(got['complete_quartic_response_used'])
        self.assertGreater(got['nfev'],0)
        self.assertAlmostEqual(got['end_ren_vector'][0],
                               0.2+0.07*math.log(0.5),places=9)
        self.assertTrue(np.isfinite(got['end_c']).all())

    def test_extended_and_frozen_differ(self):
        a=integrate_extended_eft1(*setup(True))
        b=integrate_extended_eft1(*setup(False))
        self.assertGreater(np.linalg.norm(a['end_c']-b['end_c']),1e-12)

    def test_fail_closed_missing_coupling(self):
        adapter,ren,c,basis=setup()
        broken=JointEFT1Adapter(adapter.model,adapter.mu_f,adapter.mu_s,
            lambda y,mu:({'lambdaS1':float(y[0])},
                         np.zeros((4,4,8,8))),adapter.derivative,adapter.context)
        with self.assertRaises(QuarticEvidenceError):
            integrate_extended_eft1(broken,ren,c,basis)

    def test_fail_closed_model(self):
        adapter,ren,c,basis=setup()
        wrong=JointEFT1Adapter('B_-1',adapter.mu_f,adapter.mu_s,
                               adapter.unpack,adapter.derivative,adapter.context)
        with self.assertRaises(QuarticEvidenceError):
            integrate_extended_eft1(wrong,ren,c,basis)

    def test_incomplete_quartic_inventory(self):
        adapter,ren,c,basis=setup()
        guarded=JointEFT1Adapter(adapter.model,adapter.mu_f,adapter.mu_s,
            adapter.unpack,adapter.derivative,adapter.context,
            ('lambdaS1','lambdaH1','lambdaExtraMissing'))
        with self.assertRaises(QuarticEvidenceError):
            integrate_extended_eft1(guarded,ren,c,basis)

    def test_invalid_beta_dimension(self):
        adapter,ren,c,basis=setup()
        bad=JointEFT1Adapter(adapter.model,adapter.mu_f,adapter.mu_s,
             adapter.unpack,lambda t,y:np.ones(2),adapter.context)
        with self.assertRaises(QuarticEvidenceError):
            integrate_extended_eft1(bad,ren,c,basis)


if __name__=='__main__':
    unittest.main()
