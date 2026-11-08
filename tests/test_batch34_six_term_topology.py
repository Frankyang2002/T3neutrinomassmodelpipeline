"""Targeted 1PI topology and provisional Fonseca Eq.(22) six-term RGE tests."""
from __future__ import annotations

from itertools import permutations
import unittest
import numpy as np

from RGE.running.intermediate.GeneralEFT1OneLoopTopologyAudit import (
    classify, enumerate_three_vertex_diagrams,
)
from RGE.running.intermediate.GeneralEFT1SixTermReferenceCandidate import (
    six_term_beta, six_term_beta_total, MODE,
)
from RGE.running.intermediate.GeneralEFT1IndependentOneLoopReference import (
    equation22_gaugeless, sum_reference,
)
from RGE.running.intermediate.GeneralEFT1SevenTermFlow import (
    Couplings, TensorBasis, evolve, TERMS, total_beta,
)
from RGE.running.intermediate.GeneralEFT1CrossedWardCandidate import candidate_beta


def fixture(seed=121):
    rng=np.random.default_rng(seed)
    nf,ns=3,4
    yy=rng.normal(size=(nf,nf,ns)) + 1j*rng.normal(size=(nf,nf,ns))
    y=(yy+yy.swapaxes(0,1))/2
    cc=rng.normal(size=(nf,nf,ns,ns))+1j*rng.normal(size=(nf,nf,ns,ns))
    c=(cc+cc.swapaxes(0,1)+cc.swapaxes(2,3)+cc.transpose(1,0,3,2))/4
    lam=rng.normal(size=(ns,)*4)
    lam=sum(lam.transpose(p) for p in permutations(range(4)))/24
    return c,Couplings(y,lam,())


class OneLoopTopologyTests(unittest.TestCase):
    def test_unstarred_C_yy_has_no_chiral_massless_contractions(self):
        self.assertEqual(enumerate_three_vertex_diagrams(('C','y','y')),())

    def test_Cstar_yy_crossed_external_attachment_is_all_1pr(self):
        ds=enumerate_three_vertex_diagrams(('Cstar','y','y'))
        crossed=[d for d in ds if d.external[0]==('phi',) and
                 sum('phi' in ext for ext in d.external[1:])==1]
        self.assertGreater(len(crossed),0)
        self.assertTrue(all(not d.one_particle_irreducible for d in crossed))
        self.assertTrue(all(d.external[0]==('phi',) for d in crossed))

    def test_Cstar_yy_noncrossed_Cstar_external_two_phi_has_1pi(self):
        ds=enumerate_three_vertex_diagrams(('Cstar','y','y'))
        triangle=[d for d in ds if d.external == (('phi','phi'),('psi',),('psi',))]
        self.assertGreater(len(triangle),0)
        self.assertTrue(all(d.one_particle_irreducible for d in triangle))

    def test_mixed_yukawa_has_expected_1pi_control(self):
        ds=enumerate_three_vertex_diagrams(('C','y','ystar'))
        self.assertGreater(sum(d.one_particle_irreducible for d in ds),0)
        self.assertTrue(all(len(d.edges)==3 for d in ds))

    def test_structural_report_does_not_certify_physics(self):
        report=classify()
        self.assertTrue(report['crossed_external_scalar_attachment_only_1PR'])
        self.assertTrue(report['unstarred_C_y_y_chirality_forbids_diagram'])
        self.assertFalse(report['is_complete_counterterm_derivation'])
        self.assertFalse(report['is_physical_Aminus4_certification'])


class SixTermCandidateTests(unittest.TestCase):
    def test_independent_reference_gaugeless_eq22(self):
        c,couplings=fixture()
        independent=sum_reference(equation22_gaugeless(c,couplings.yukawa,couplings.quartic))
        parts=six_term_beta(c,couplings)
        self.assertEqual(tuple(parts), TERMS)
        np.testing.assert_array_equal(parts['crossed_yukawa'],np.zeros_like(c))
        np.testing.assert_allclose(total_beta(parts),independent,rtol=8e-13,atol=2e-12)

    def test_six_term_excludes_nonzero_crossed_Cstar_contraction(self):
        c,couplings=fixture()
        seven=candidate_beta(c,couplings)
        six=six_term_beta(c,couplings)
        self.assertGreater(np.linalg.norm(seven['crossed_yukawa']),1e-6)
        np.testing.assert_allclose(total_beta(seven)-total_beta(six),
                                   seven['crossed_yukawa'],rtol=3e-14,atol=3e-14)

    def test_full_symmetric_ode_zeroes_crossed_channel(self):
        # Three synthetic lepton generations, four real scalars.
        basis=TensorBasis(6,4,((0,1),(2,3),(4,5)),(0,1,2,3))
        c=np.zeros((6,6,4,4),complex)
        c[0,0,2,2]=1.1e-10 + .3e-10j
        y=np.zeros((6,6,4),complex)
        y[0,1,2]=y[1,0,2]=.2+.1j
        y[3,3,1]=.07-.04j
        lam=np.zeros((4,)*4,dtype=float)
        lam[2,2,2,2]=.2
        lo=Couplings(y*.9,lam*.9,())
        hi=Couplings(y,lam,())
        res=evolve(basis,c,1e5,2e3,[2e3,1e5],[lo,hi],points=3,
                   rtol=1e-8,atol=1e-21,beta_mode=MODE)
        self.assertEqual(res['beta_mode'],MODE)
        self.assertFalse(res['beta_physics_certified'])
        self.assertTrue(res['solver']['success'])
        self.assertEqual(res['anti_linear_term_names'],['conjugate_coefficient_yukawa'])
        np.testing.assert_array_equal(res['beta_initial']['crossed_yukawa'],np.zeros_like(c))
        np.testing.assert_array_equal(res['beta_final']['crossed_yukawa'],np.zeros_like(c))

    def test_unknown_mode_still_fails_closed(self):
        self.assertEqual(MODE,'fonseca_eq22_six_term')
        c,C=fixture()
        self.assertEqual(set(six_term_beta(c,C)),set(TERMS))


if __name__=='__main__':
    unittest.main()
