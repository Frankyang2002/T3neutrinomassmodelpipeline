"""A(-4) declared Weyl/UV-coupling dictionary, including negative controls."""
from __future__ import annotations

import unittest
import numpy as np

from Numerical.lo_certification.spinor import a_minus4_boundary_verify as b
from Numerical.lo_certification.spinor import a_minus4_master_basis_gate as m
from Numerical.lo_certification.spinor import a_minus4_weyl_operator_dictionary as w
from Numerical.lo_certification.spinor import a_minus4_declared_convention_boundary as d


def cg_native() -> np.ndarray:
    x=3**(-.25)
    g=np.zeros((2,2,3),complex)
    g[0,0,0]=x
    g[1,1,0]=-x
    g[0,1,1]=g[1,0,1]=-x
    g[0,0,2]=g[1,1,2]=1j*x
    return g


class TestDeclaredWeylUVCouplings(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        rng=np.random.default_rng(5003)
        cls.y1=(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))/7
        cls.y2=(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))/8
        cls.mass=np.array([151.,191.,227.])
        cls.cg=cg_native()
        cls.n=m.triplet_intertwiner(cls.cg)['triplet_N']
        cls.alpha=np.array([0.7,-1.3,2.1])
        cls.beta=np.array([-2.3,1.1,0.4])
        cls.sigma1=.38
        cls.sigma2=-.81
        cls.phi=np.array([0.2,-0.3,0.4,-0.5,0.6,0.7,0.8,-0.9,1.1,0.2,-1.2,0.5])

    def bridge(self,rho=1):
        return d.uv_couplings_in_declared_basis(
            self.y1,self.y2,self.mass,self.alpha,self.beta,
            self.sigma1,self.sigma2,rho=rho)

    def test_original_w_is_full_complex_3_by_3(self):
        z=d.complex_symmetric_w(self.y1,self.y2,self.mass)
        self.assertEqual(z.shape,(3,3))
        self.assertGreater(np.max(np.abs(z-z.T)),1e-5)
        self.assertGreater(np.max(np.abs(z.imag)),1e-5)

    def test_couplings_transform_into_expected_PL_operator(self):
        y1,y2,mass=self.bridge()
        actual=d.complex_symmetric_w(y1,y2,mass)
        expected=d.expected_pl_phase(
            d.complex_symmetric_w(self.y1,self.y2,self.mass),
            self.alpha,self.sigma1,self.sigma2,rho=1)
        self.assertLess(np.max(np.abs(actual-expected)),1e-15)

    def test_negative_native_CC_branch_changes_source_coupling(self):
        plus_y1,plus_y2,_=self.bridge(1)
        minus_y1,minus_y2,_=self.bridge(-1)
        self.assertLess(np.max(np.abs(plus_y1+minus_y1)),1e-14)
        self.assertLess(np.max(np.abs(plus_y2-minus_y2)),1e-14)
        self.assertGreater(np.max(np.abs(plus_y1-minus_y1)),0.01)

    def test_equivalent_gamma_and_y1_redefinition_same_master_theory(self):
        a1,a2,_=self.bridge(1)
        b1,b2,_=d.uv_couplings_in_declared_basis(
            -self.y1,self.y2,self.mass,self.alpha,self.beta,
            self.sigma1,self.sigma2,rho=-1)
        self.assertLess(np.max(np.abs(a1-b1)),1e-14)
        self.assertLess(np.max(np.abs(a2-b2)),1e-14)

    def test_independent_vectorlike_heavy_phases_cancel_with_positive_mass(self):
        y1,y2,_=self.bridge(1)
        z1,z2,mass=d.uv_couplings_in_declared_basis(
            self.y1,self.y2,self.mass,self.alpha,self.beta+.7,
            self.sigma1,self.sigma2,rho=1)
        self.assertTrue(np.array_equal(mass,self.mass))
        self.assertLess(np.max(np.abs(d.complex_symmetric_w(y1,y2,self.mass)-
                                       d.complex_symmetric_w(z1,z2,mass))),1e-15)
        self.assertGreater(np.max(np.abs(y1-z1)),0.01)

    def test_scalar_native_rotation_has_exact_inverse(self):
        p=d.scalar_native_from_master(self.phi,self.n,self.sigma1,self.sigma2)
        back=p.copy()
        back[4:6]=m.realification(np.array([[np.exp(-1j*self.sigma1)]]))@p[4:6]
        back[6:12]=m.realification(np.exp(-1j*self.sigma2)*self.n)@p[6:12]
        self.assertLess(np.max(np.abs(back-self.phi)),1e-14)

    def test_full_complex_yukawa_all_Grassmann_components(self):
        r=d.grassmann_rephasing_witness(self.y1,self.y2,self.mass,self.cg,self.n,
                                       seed=1123,samples=5)
        self.assertLess(r['max_full_grassmann_relative_residual'],3e-12)
        self.assertLess(r['max_w_matrix_relative_residual'],3e-12)
        self.assertFalse(r['gamma_cc_branch_physically_discriminated'])

    def test_uv_barred_S2_CG_transports_with_N_not_conjugate_N(self):
        g=3**(-.25)
        g1=np.array([[0,1],[-1,0]],complex)
        g2=np.zeros((2,2,3),complex)
        g2[0,0,1]=g
        g2[0,1,0]=g
        g2[0,1,2]=-1j*g
        g2[1,0,0]=g
        g2[1,0,2]=1j*g
        g2[1,1,1]=-g
        match=np.einsum('ah,bhc->abc',g1.conj(),g2.conj())
        self.assertLess(np.max(np.abs(match-self.cg)),1e-14)
        g2_master=np.einsum('bhc,kc->bhk',g2,self.n)
        matched_master=np.einsum('abc,kc->abk',match,self.n.conj())
        result=np.einsum('ah,bhk->abk',g1.conj(),g2_master.conj())
        self.assertLess(np.max(np.abs(result-matched_master)),1e-14)
        wrong=np.einsum('bhc,kc->bhk',g2,self.n.conj())
        wrong_matched=np.einsum('ah,bhk->abk',g1.conj(),wrong.conj())
        self.assertGreater(np.max(np.abs(wrong_matched-matched_master)),0.1)

    def test_scalar_sigma2_opposite_sign_fails_wilson_transport(self):
        y1,y2,mass=self.bridge(1)
        actual=d.complex_symmetric_w(y1,y2,mass)
        wrong=d.expected_pl_phase(d.complex_symmetric_w(self.y1,self.y2,self.mass),
                                   self.alpha,self.sigma1,-self.sigma2,rho=1)
        self.assertGreater(np.max(np.abs(actual-wrong)),1e-5)

    def test_missing_left_lepton_phase_fails_transport(self):
        y1,y2,mass=self.bridge(1)
        actual=d.complex_symmetric_w(y1,y2,mass)
        wrong=d.expected_pl_phase(d.complex_symmetric_w(self.y1,self.y2,self.mass),
                                   np.zeros(3),self.sigma1,self.sigma2,rho=1)
        self.assertGreater(np.max(np.abs(actual-wrong)),1e-5)

    def test_wrong_rho_on_one_side_cannot_match(self):
        w_native=d.complex_symmetric_w(self.y1,self.y2,self.mass)
        expected=d.expected_pl_phase(w_native,self.alpha,self.sigma1,self.sigma2,rho=-1)
        actual=d.complex_symmetric_w(*self.bridge(1))
        self.assertGreater(np.max(np.abs(actual-expected)),1e-5)

    def test_invalid_rho_bool_zero_two_fail_closed(self):
        for bad in (True,0,2,-2,0.5):
            with self.subTest(bad=bad), self.assertRaises(d.ConventionBlocked):
                d.uv_couplings_in_declared_basis(
                    self.y1,self.y2,self.mass,self.alpha,self.beta,
                    self.sigma1,self.sigma2,rho=bad)

    def test_axial_mass_phase_is_excluded(self):
        # F_L and F_R must have the SAME beta. A relative phase would
        # change -M Fbar F into a complex chiral mass, outside this gate.
        difference=0.83
        self.assertGreater(abs(np.exp(1j*difference)-1),0.1)
        # The supported UV transport accepts exactly ONE beta per heavy mode.
        with self.assertRaises(d.ConventionBlocked):
            d.uv_couplings_in_declared_basis(
                self.y1,self.y2,self.mass,self.alpha,np.zeros((3,2)),
                self.sigma1,self.sigma2,rho=1)

    def test_malformed_mass_or_phase_fails_closed(self):
        with self.assertRaises(d.ConventionBlocked):
            d.uv_couplings_in_declared_basis(
                self.y1,self.y2,np.array([1.,-1.,2.]),self.alpha,self.beta,
                self.sigma1,self.sigma2,rho=1)
        with self.assertRaises(d.ConventionBlocked):
            d.uv_couplings_in_declared_basis(
                self.y1,self.y2,self.mass,self.alpha,self.beta,
                np.nan,self.sigma2,rho=1)

    def test_triplet_map_nonunitary_fails_closed(self):
        bad=self.n.copy()
        bad[0,0]*=2
        with self.assertRaises(d.ConventionBlocked):
            d.scalar_native_from_master(self.phi,bad,self.sigma1,self.sigma2)

    def test_symmetric_coordinate_orbit_and_rho_sign(self):
        qmaster=(self.cg.reshape(4,3)@self.n.conj().T).reshape(2,2,3)
        w_native=d.complex_symmetric_w(self.y1,self.y2,self.mass)
        independent,_=b.symmetrise(b.raw_components(w_native,qmaster))
        rows=[{'i':k[0],'j':k[1],'a':k[2],'b':k[3],
               'coefficient_real_if_rho_1':v.real,
               'coefficient_imag_if_rho_1':v.imag,
               'distinct_orbit_size':len(b.orbit(*k))}
              for k,v in sorted(independent.items())]
        pos=d.encode_rephased_coordinate(rows,rho=1)
        neg=d.encode_rephased_coordinate(rows,rho=-1)
        self.assertEqual(len(pos),len(neg))
        self.assertEqual(len(pos),84)
        for a,z in zip(pos,neg):
            self.assertAlmostEqual(a['value_re_gev_inv'],-z['value_re_gev_inv'])
            self.assertAlmostEqual(a['value_im_gev_inv'],-z['value_im_gev_inv'])
        rows[0]['distinct_orbit_size']+=1
        with self.assertRaises(d.ConventionBlocked):
            d.encode_rephased_coordinate(rows,rho=1)

    def test_declared_fullsum_not_extra_half(self):
        qmaster=(self.cg.reshape(4,3)@self.n.conj().T).reshape(2,2,3)
        w_master=d.complex_symmetric_w(*self.bridge(1))
        _,mono=b.symmetrise(b.raw_components(w_master,qmaster))
        correct=w.master_fullsum_operator(mono,self.phi)
        wrong=w.master_fullsum_operator(mono,self.phi,coefficient_factor=0.5)
        self.assertGreater(w.diff(correct,wrong)[1],0.25)


if __name__=='__main__':
    unittest.main()
