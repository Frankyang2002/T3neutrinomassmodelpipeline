"""A(-4): actual Grassmann spinor-operator identity and negative controls."""
from __future__ import annotations

import unittest

import numpy as np

from Numerical.lo_certification.spinor import a_minus4_boundary_verify as b
from Numerical.lo_certification.spinor import a_minus4_master_basis_gate as m
from Numerical.lo_certification.spinor import a_minus4_weyl_operator_dictionary as d


def native_cg() -> np.ndarray:
    g=3**(-.25)
    result=np.zeros((2,2,3),complex)
    result[0,0,0]=g
    result[1,1,0]=-g
    result[0,1,1]=result[1,0,1]=-g
    result[0,0,2]=result[1,1,2]=1j*g
    return result


class TestExactWeylDictionary(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        rng=np.random.default_rng(456)
        cls.g=native_cg()
        cls.y1=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
        cls.y2=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
        cls.masses=np.array([111.,159.,202.])
        cls.w=b.dense_flavour(cls.y1,cls.y2,cls.masses)
        cls.n=m.triplet_intertwiner(cls.g)['triplet_N']
        qmaster=(cls.g.reshape(4,3)@cls.n.conj().T).reshape(2,2,3)
        cls.independent,cls.mono=b.symmetrise(b.raw_components(cls.w,qmaster))
        cls.phi=np.array([0.,0.,0.,0.,1.2,-.3,.4,.7,-.5,.8,1.1,-.2])
        cls.oldphi=cls.phi.copy()
        cls.oldphi[6:12]=m.realification(cls.n).T@cls.phi[6:12]
        cls.native=d.native_pl_operator(cls.w,cls.g,cls.oldphi)

    def test_distinct_spinors_have_symmetric_grassmann_contraction(self):
        for i in range(1,7):
            for j in range(1,7):
                self.assertEqual(d.weyl_epsilon(i,j),d.weyl_epsilon(j,i))

    def test_identical_Weyl_field_carries_two_Grassmann_monomials(self):
        self.assertEqual(d.weyl_epsilon(2,2),{(2,3):2})

    def test_Grassmann_generators_anticommute_and_nilpotent(self):
        self.assertEqual(d.grassmann_product(4,7),{(4,7):1})
        self.assertEqual(d.grassmann_product(7,4),{(4,7):-1})
        self.assertEqual(d.grassmann_product(3,3),{})

    def test_source_derived_general_complex_3x3_full_operator(self):
        result=d.full_grassmann_probe(self.w,self.g,self.n,self.mono)
        self.assertLess(result['max_relative_Grassmann_operator_residual'],3e-12)
        self.assertGreater(result['full_degree_two_Grassmann_component_count'],20)
        self.assertFalse(result['rho_sign_discriminated'])

    def test_full_Grassmann_master_operator_identical_native_PL(self):
        new=d.master_fullsum_operator(self.mono,self.phi)
        self.assertLess(d.diff(self.native,new)[1],2e-13)

    def test_wilson_orbit_factor_two_missing_is_rejected(self):
        wrong=d.master_fullsum_operator(self.mono,self.phi,coefficient_factor=.5)
        self.assertGreater(d.diff(self.native,wrong)[1],0.25)

    def test_opposite_cc_sign_without_action_transport_is_rejected(self):
        wrong=d.master_fullsum_operator(self.mono,self.phi,coefficient_factor=-1)
        self.assertGreater(d.diff(self.native,wrong)[1],1.5)

    def test_both_consistent_sign_choices_are_indistinguishable(self):
        base=d.master_fullsum_operator(self.mono,self.phi)
        for rho in (-1,1):
            l={k:rho*v for k,v in self.native.items()}
            r={k:rho*v for k,v in base.items()}
            self.assertLess(d.diff(l,r)[1],2e-13)

    def test_flavour_antisymmetric_piece_cancels_after_Grassmann_symmetrization(self):
        # Q[a,b,c] is symmetric under a<->b. The antisymmetric part of W
        # contributes only to the ordered source, not to the physical LL term.
        swapped=d.native_pl_operator(self.w.T,self.g,self.oldphi)
        projected=d.native_pl_operator((self.w+self.w.T)/2,self.g,self.oldphi)
        self.assertLess(d.diff(self.native,swapped)[1],2e-13)
        self.assertLess(d.diff(self.native,projected)[1],2e-13)
        self.assertGreater(np.max(np.abs(self.w-self.w.T)),.01)

    def test_complex_Yukawas_must_be_conjugated_in_PL(self):
        wrong=d.native_pl_operator(np.conjugate(self.w),self.g,self.oldphi)
        self.assertGreater(d.diff(self.native,wrong)[1],1e-2)

    def test_wrong_native_triplet_axis_order_detected(self):
        wrong_phi=self.oldphi.copy()
        wrong_phi[6:12]=self.phi[6:12]  # forgetting the nontrivial triplet map
        wrong=d.native_pl_operator(self.w,self.g,wrong_phi)
        self.assertGreater(d.diff(self.native,wrong)[1],1e-2)

    def test_wrong_complex_scalar_conjugation_detected(self):
        wrong_phi=self.oldphi.copy()
        wrong_phi[5]*=-1  # replaces S1* with S1
        wrong=d.native_pl_operator(self.w,self.g,wrong_phi)
        self.assertGreater(d.diff(self.native,wrong)[1],1e-2)

    def test_invalid_lepton_and_scalar_coordinates_fail_closed(self):
        with self.assertRaises(d.DictionaryBlocked):d.weyl_epsilon(7,1)
        with self.assertRaises(d.DictionaryBlocked):d.master_fullsum_operator(
            {(1,2,7,5):1.},self.phi)
        with self.assertRaises(d.DictionaryBlocked):d.native_pl_operator(
            np.eye(2),self.g,self.oldphi)

    def test_no_implicit_plus_one_sign_certificate(self):
        result=d.full_grassmann_probe(self.w,self.g,self.n,self.mono)
        self.assertEqual([x['rho_C0'] for x in result['rho_branches']],[1,-1])
        self.assertFalse(result['rho_sign_discriminated'])


if __name__=='__main__':
    unittest.main()
