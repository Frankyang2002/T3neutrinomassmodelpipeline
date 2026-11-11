"""Standalone regression tests for explicit running quartic -> seven-term Wilson coupling."""
from __future__ import annotations
import unittest
from dataclasses import replace
import numpy as np

from Numerical.running.T3RunningQuarticTensor import (
    RunningQuarticBasis, QuarticEvidenceError,
)
from Numerical.running.T3QuarticWilsonBridge import (
    running_context, running_seven_term_beta,
)
from RGE.running.intermediate.GeneralEFT1C5Mixing import (
    Context, calculate_master_beta, TERMS,
)


def fixture():
    n = 12  # A_-4: H(4) + S1(2) + S2(6) Cartesian components
    q = np.zeros((n, n, n, n))
    q[4, 4, 4, 4] = 24.0
    hhss = np.zeros_like(q)
    for k in set(__import__('itertools').permutations((0, 0, 4, 4))):
        hhss[k] = 2.0
    payload = {
        'format': 't3-running-real-quartic-response-v1',
        'model': 'A_-4',
        'convention': 'V4=lambda_abcd phi_a phi_b phi_c phi_d/4!',
        'basis': 'fixed-original-Cartesian',
        'source': 'SYNTHETIC-UNIT-TEST-NOT-PHYSICAL-EVIDENCE',
        'constant': np.zeros_like(q).tolist(),
        'responses': [
            {'symbol': 'lambdaS1', 'component': 'real', 'tensor': q.tolist()},
            {'symbol': 'lambdaH1', 'component': 'real', 'tensor': hhss.tolist()},
        ],
    }
    basis = RunningQuarticBasis.from_payload(payload, {'lambdaS1', 'lambdaH1'})
    c = np.zeros((9,9,n,n), complex)
    c[0,0,4,4] = 0.7 + 0.2j
    ctx = Context(
        model='A_-4', nf=9, ns=n, leptons=((0,1),(2,3),(4,5)),
        higgs=(0,1,2,3), mu_f=1e5, mu_s=1e3,
        yukawa=np.zeros((9,9,n), complex),
        quartic=np.zeros_like(q), gauge=(), operators=(('synthetic',c),),
        authority='synthetic', certificates={},
    )
    return ctx, c, basis


class TestBridge(unittest.TestCase):
    def test_preserves_every_term_and_matches_direct(self):
        ctx,c,basis=fixture()
        state={'lambdaS1':0.2,'lambdaH1':0.3}
        got,terms=running_seven_term_beta(ctx,c,state,basis)
        direct=calculate_master_beta(replace(ctx,quartic=basis.evaluate(state)),c)
        self.assertEqual(set(terms),set(TERMS))
        np.testing.assert_allclose(got,sum(direct.values()))

    def test_scaling_is_instantaneous(self):
        ctx,c,basis=fixture()
        a,_=running_seven_term_beta(ctx,c,{'lambdaS1':0.2,'lambdaH1':0},basis)
        b,_=running_seven_term_beta(ctx,c,{'lambdaS1':0.4,'lambdaH1':0},basis)
        self.assertAlmostEqual((b-a)[0,0,4,4].real,a[0,0,4,4].real)
        self.assertNotEqual(a[0,0,4,4],0)

    def test_hhss_parity_passes_and_fails(self):
        ctx,c,basis=fixture()
        state={'lambdaS1':0.2,'lambdaH1':0.3}
        hhss=basis.evaluate(state)[:4,:4,4:,4:]
        running_context(ctx,state,basis,expected_hhss=hhss)
        with self.assertRaises(QuarticEvidenceError):
            running_context(ctx,state,basis,expected_hhss=hhss+1)

    def test_model_mismatch_fails(self):
        ctx,c,basis=fixture()
        with self.assertRaises(QuarticEvidenceError):
            running_context(replace(ctx,model='C_-1'),{'lambdaS1':0.2,'lambdaH1':0.3},basis)

    def test_missing_running_state_fails(self):
        ctx,c,basis=fixture()
        with self.assertRaises(QuarticEvidenceError):
            running_seven_term_beta(ctx,c,{'lambdaS1':0.2},basis)

    def test_additional_gamma_contribution(self):
        ctx,c,basis=fixture()
        state={'lambdaS1':0.2,'lambdaH1':0.3}
        bare,_=running_seven_term_beta(ctx,c,state,basis)
        gamma=np.zeros((12,12))
        gamma[4,4]=0.17
        modified,_=running_seven_term_beta(
            ctx,c,state,basis,extra_scalar_anomalous_dimension=gamma)
        self.assertAlmostEqual((modified-bare)[0,0,4,4].real,2*0.17*c[0,0,4,4].real)

if __name__=='__main__':
    unittest.main()
