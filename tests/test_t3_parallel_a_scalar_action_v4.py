import unittest
import numpy as np
from Numerical.parallel_t3.wilson.scalar_action_audit import (
    _permutations,real_lam,complex_field_from_real,monomial,
    baseline_invariants,audit_baseline_all16,
)
from Numerical.parallel_t3.wilson.matchete_rgbeta_dictionary import DictionaryEvidenceError

class TestT3ParallelAScalarActionV4(unittest.TestCase):
    def test_baseline_six_all16(self):
        result=audit_baseline_all16(samples=2)
        self.assertEqual(result["baseline_norm_coupling_slots"],96)
        self.assertLess(result["max_baseline_norm_action_residual"],1e-12)
    def test_source_permutation_validation(self):
        self.assertEqual(len(_permutations('{{1,2,3,4},{2,1,3,4}}')),2)
        with self.assertRaises(DictionaryEvidenceError):
            _permutations('{{1,2,3,4},{1,2,3,4}}')
    def test_single_complex_field_quartic_action(self):
        t=np.ones((1,1,1,1)); species=['H','Bar[H]','H','Bar[H]']
        x=np.array([1.2,-0.7]);c=.35+.6j
        lam=real_lam(t,species,c,_permutations('{{1,2,3,4},{2,1,3,4}}'))
        lhs=np.einsum('ijkl,i,j,k,l',lam,x,x,x,x)/24
        rhs=2*(c*monomial(t,species,complex_field_from_real(species,[1]*4,x))).real
        self.assertAlmostEqual(lhs,rhs,places=12)
    def test_triple_species_and_complex_epsilon(self):
        rng=np.random.default_rng(12)
        species=['H','Bar[S1]','S2','Bar[S2]'];dims=[2,2,3,3]
        t=rng.normal(size=dims)+1j*rng.normal(size=dims)
        x=rng.normal(size=14); c=.31-.43j
        permutations=_permutations('{{1,2,3,4},{4,2,1,3},{2,3,4,1}}')
        lam=real_lam(t,species,c,permutations)
        lhs=np.einsum('ijkl,i,j,k,l',lam,x,x,x,x)/24
        fields=complex_field_from_real(species,dims,x)
        rhs=2*(c*monomial(t,species,fields)).real
        self.assertAlmostEqual(lhs,rhs,places=10)
    def test_antihermitian_portal_needs_i(self):
        # Analogue of T3HiggsPortalHermitianPhase: purely imaginary
        # antisymmetric Hermitian bilinear is cancelled by an H.c. unless
        # multiplied by I before taking the H.c.
        T=np.array([[0,1],[-1,0]],dtype=complex)
        self.assertTrue(np.allclose(T.conj().T,-T))
        z=np.array([1+2j,3-1j])
        v=z.conj()@T@z
        self.assertAlmostEqual((v+v.conjugate()).real,0)
        self.assertGreater(abs((1j*v+(1j*v).conjugate()).real),1e-6)

if __name__=='__main__':unittest.main()
