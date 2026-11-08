"""Batch 1: source-CG -> master SU2 component and phase-covariance gates."""
from __future__ import annotations

import unittest

import numpy as np

from Numerical.lo_certification.spinor import a_minus4_master_basis_gate as m
from Numerical.lo_certification.spinor import a_minus4_boundary_verify as b


def native_source_CG():
    g=3.**(-.25)
    q=np.zeros((2,2,3),complex)
    q[0,0,0]=g
    q[1,1,0]=-g
    q[0,1,1]=q[1,0,1]=-g
    q[0,0,2]=q[1,1,2]=1j*g
    return q


class TestMasterTripletBasis(unittest.TestCase):
    def setUp(self):
        self.q=native_source_CG()

    def test_registered_source_triplet_map_expected(self):
        r=m.triplet_intertwiner(self.q)
        wanted=np.array([[-1/np.sqrt(2),0,1j/np.sqrt(2)],
                         [0,1,0],
                         [1/np.sqrt(2),0,1j/np.sqrt(2)]],complex)
        np.testing.assert_allclose(r['triplet_N'],wanted,atol=1e-14,rtol=0)
        self.assertLess(r['generator_residual'],3e-13)
        self.assertLess(r['cg_covariance_residual'],3e-13)
        self.assertLess(r['cg_standard_residual'],3e-13)

    def test_independently_recorded_native_matchete_cartan(self):
        native_cartan=np.array([[0,0,-1j],[0,0,0],[1j,0,0]],complex)
        m.triplet_intertwiner(self.q,native_cartan=native_cartan)

    def test_false_native_generator_rejected(self):
        with self.assertRaises(m.BasisBlocked):
            m.triplet_intertwiner(self.q,native_cartan=-np.array([[0,0,-1j],[0,0,0],[1j,0,0]],complex))

    def test_wrong_master_generators_rejected(self):
        with self.assertRaises(m.BasisBlocked):
            m.triplet_intertwiner(self.q,master_generators=[np.eye(3)]*3)

    def test_tampered_native_CG_phase_rejected(self):
        changed=self.q.copy();changed[1,0,1]*=-1
        with self.assertRaises(m.BasisBlocked):m.triplet_intertwiner(changed)

    def test_complex_to_interleaved_real_map(self):
        r=m.triplet_intertwiner(self.q)
        n=r['triplet_N']; t=r['real_R']
        z=np.array([1.2+3.1j,-1.7+.5j,2.2-2.6j])
        re=np.array([v for a in z for v in (a.real,a.imag)])
        mapped=n@z
        expected=np.array([v for a in mapped for v in (a.real,a.imag)])
        np.testing.assert_allclose(t@re,expected,atol=1e-14,rtol=0)
        np.testing.assert_allclose(t.T@t,np.eye(6),atol=1e-14,rtol=0)

    def test_complex_offdiagonal_flavour_tensor_transport(self):
        rng=np.random.default_rng(13082)
        y1=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
        y2=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
        masses=np.array([121.2,178.3,213.4])
        w=b.dense_flavour(y1,y2,masses)
        original=b.raw_components(w,self.q)
        native_indep,old=b.symmetrise(original)
        n=m.triplet_intertwiner(self.q)
        e,j=m.normalized_pair_embedding()
        q_master=(np.sqrt(n['native_CG_norm_sq'])*(e@j)).reshape(2,2,3)
        newindep,new=b.symmetrise(b.raw_components(w,q_master))
        rotated=m.transport_monomials(old,n['real_R'])
        keys=set(rotated)|set(new)
        self.assertLess(max(abs(rotated.get(k,0)-new.get(k,0)) for k in keys)/max(map(abs,new.values())),2e-13)
        metrics=m.verify_transport(old,new,n['real_R'])
        self.assertLess(metrics['commuting_polynomial_relative_max_residual'],2e-12)
        self.assertFalse(metrics['physical_spinor_operator_identity_established'])
        self.assertGreater(len(new),0)

    def test_invalid_independent_orbit_rejected(self):
        entries=[{'i':1,'j':1,'a':5,'b':7,'distinct_orbit_size':4,
                  'coefficient_real_if_rho_1':2.,'coefficient_imag_if_rho_1':-3.}]
        with self.assertRaises(m.BasisBlocked):m.independent_to_monomials(entries)

    def test_wrong_triplet_map_does_not_preserve_polynomial(self):
        rng=np.random.default_rng(909)
        w=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
        _,old=b.symmetrise(b.raw_components(w,self.q))
        e,j=m.normalized_pair_embedding()
        good=m.triplet_intertwiner(self.q)
        q_master=np.sqrt(good['native_CG_norm_sq'])*(e@j)
        _,new=b.symmetrise(b.raw_components(w,q_master.reshape(2,2,3)))
        with self.assertRaises(m.BasisBlocked):m.verify_transport(old,new,np.eye(6))

    def test_dirac_heavy_and_lepton_phase_covariance(self):
        proof=m.phase_covariance_probe()
        self.assertLess(proof['arbitrary_complex_dirac_mass_rephasing_relative_residual'],5e-14)
        self.assertLess(proof['three_independent_lepton_phase_covariance_relative_residual'],5e-14)
        self.assertFalse(proof['native_Matchete_to_master_rephasing_dictionary_certified'])

    def test_three_family_basis_constructed_via_real_master_api(self):
        result=m.three_family_api_basis(b.model_basis_report())
        self.assertEqual(result['full_tensor_shape'],[45,45,12,12])
        self.assertEqual(result['production_build_sm_eft_fermion_dimension'],15)
        self.assertFalse(result['production_build_sm_eft_is_three_family'])


if __name__=='__main__':unittest.main()
