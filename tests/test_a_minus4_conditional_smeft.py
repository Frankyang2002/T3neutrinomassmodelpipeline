"""Independent A(-4) conditional scalar matching invariants and API checks.

Tests using existing project matchers run when installed in the local checkout;
all pure tensor/convention tests also run outside the repository.
"""
from __future__ import annotations

import importlib.util
import math
from pathlib import Path
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/'Numerical/joint_scalar_matching/Aminus4ConditionalSMEFT.py'
spec=importlib.util.spec_from_file_location('a4joint',SRC)
a4=importlib.util.module_from_spec(spec)
spec.loader.exec_module(a4)
HAS_MATCHER=(ROOT/'Numerical/joint_scalar_matching/core.py').is_file() and (ROOT/'RGE/running/intermediate/FullSU2SourceResponse.py').is_file()


def local_positions():
    return np.array([[36,37],[38,39],[40,41]])


def c_from_c5(flav):
    c=np.zeros((45,45,12,12),complex)
    t=a4.master_weinberg_template()
    pos=local_positions()
    for p in range(3):
        for q in range(3):
            for i in range(2):
                for j in range(2):
                    c[pos[p,i],pos[q,j],:4,:4]=flav[p,q]*t[i,j]
    return c


class PureConditionalScalarTests(unittest.TestCase):
    def test_master_weinberg_norm_and_full_pair_symmetries(self):
        t=a4.master_weinberg_template()
        self.assertEqual(t.shape,(2,2,4,4))
        self.assertAlmostEqual(t[0,0,2,2],.25)
        self.assertTrue(np.array_equal(t,t.transpose(1,0,2,3)))
        self.assertTrue(np.array_equal(t,t.transpose(0,1,3,2)))
        self.assertAlmostEqual(t[1,1,0,0],.25)
        self.assertAlmostEqual(t[0,1,2,0],-.125)

    def test_lepton_embedding_keeps_all_45_Weyl_slots(self):
        f=np.diag([1,2,3]).astype(complex)*1e-15
        c=c_from_c5(f)
        self.assertEqual(c.shape,(45,45,12,12))
        self.assertFalse(np.any(c[:36]))
        d=a4.lepton_tensor(c,local_positions(),heavy=False)
        self.assertEqual(d.shape,(3,3,2,2,4,4))
        self.assertTrue(np.array_equal(d[1,1],2e-15*a4.master_weinberg_template()))
        self.assertFalse(np.any(a4.lepton_tensor(c,local_positions(),heavy=True)))

    def test_fails_on_duplicate_lepton_indices_and_wrong_C(self):
        c=np.zeros((45,45,12,12),complex)
        pos=local_positions();pos[0,0]=pos[1,0]
        with self.assertRaisesRegex(a4.ConditionalScalarBlocked,'positions'):
            a4.lepton_tensor(c,pos,heavy=False)
        with self.assertRaisesRegex(a4.ConditionalScalarBlocked,'45x45x12x12'):
            a4.lepton_tensor(np.zeros((6,6,12,12)),local_positions(),heavy=True)

    def test_unbroken_mass2_has_two_and_six_degenerate_eigenvalues(self):
        m=a4.gauge_invariant_mass2(1.2e6,2.3e6)
        self.assertEqual(m.shape,(8,8))
        self.assertTrue(np.array_equal(m,np.diag(np.r_[np.repeat(1.2e6,2),np.repeat(2.3e6,6)])))
        self.assertEqual(np.linalg.eigvalsh(m).shape,(8,))

    def test_no_fake_masses_or_tachyons(self):
        for x,y in ((0,1),(-1,1),(1,float('nan')),(1,float('inf'))):
            with self.subTest(x=x,y=y),self.assertRaisesRegex(a4.ConditionalScalarBlocked,'positive'):
                a4.gauge_invariant_mass2(x,y)

    def test_neutrino_proxy_takagi_singular_values_and_sign(self):
        diag=np.diag([1e-15,2e-15,4e-15]).astype(complex)
        diag[0,1]=diag[1,0]=(0.3+0.1j)*1e-15
        p=a4.neutrino_mass_proxy(diag,246.22)
        minus=a4.neutrino_mass_proxy(-diag,246.22)
        self.assertTrue(np.allclose(p['takagi_singular_values_eV'],minus['takagi_singular_values_eV']))
        self.assertTrue(np.allclose(p['takagi_singular_values_eV'],np.sort(np.linalg.svd(diag,compute_uv=False))*246.22**2/2*1e9))
        self.assertFalse(p['physical_neutrino_masses_certified'])
        with self.assertRaisesRegex(a4.ConditionalScalarBlocked,'not symmetric'):
            a4.neutrino_mass_proxy(np.triu(np.ones((3,3),complex)),246.22)

    def test_full_su2_projection_requires_correct_half_factor(self):
        c5=np.array([[1,2j,3],[2j,2,1j],[3,1j,4]],complex)*1e-15
        c=c_from_c5(c5)
        amplitude=a4.lepton_tensor(c,local_positions(),heavy=False)
        t=a4.master_weinberg_template()
        proj=np.einsum('ijab,pqijab->pq',t.conj(),amplitude)/np.vdot(t,t)
        self.assertTrue(np.allclose(proj,c5,rtol=1e-15,atol=1e-30))
        self.assertTrue(np.allclose(4*amplitude[:,:,0,0,2,2],c5,rtol=1e-15,atol=1e-30))

    def test_reject_wrong_master_array_decode(self):
        with self.assertRaisesRegex(a4.ConditionalScalarBlocked,'real/imag'):
            a4.dec({'value':1},(3,3))
        with self.assertRaisesRegex(a4.ConditionalScalarBlocked,'wrong shape'):
            a4.dec({'real':[1,2],'imag':[0,0]},(3,3))


@unittest.skipUnless(HAS_MATCHER,'Full repository joint scalar matching modules unavailable')
class ExistingRepositoryMatcherTests(unittest.TestCase):
    def test_existing_joint_weinberg_template_is_twice_master(self):
        from Numerical.joint_scalar_matching.core import weinberg_tensor
        self.assertTrue(np.allclose(weinberg_tensor(),2*a4.master_weinberg_template(),atol=1e-15))

    def test_joint_bubble_off_diagonal_scalar_pair_multiplicity(self):
        from Numerical.joint_scalar_matching.core import scalar_bubble_matching,b0bar
        c=np.zeros((3,3,2,2,8,8),complex)
        q=np.zeros((4,4,8,8),float)
        c[0,0,0,0,0,2]=c[0,0,0,0,2,0]=2e-12
        q[2,2,0,2]=q[2,2,2,0]=.4
        m=np.diag(np.r_[np.repeat(1.2e6,2),np.repeat(2.3e6,6)])
        out=scalar_bubble_matching(c,q,m,1e3)
        self.assertEqual(len(out['pair_data']),36)
        # k=0 (S1) and k=2 (S2) give the two ordered contractions.
        expected=-2*(2e-12*.4)*b0bar(1.2e6,2.3e6,1e3)/(2*16*math.pi**2)
        self.assertAlmostEqual(out['amplitude'][0,0,0,0,2,2].real/expected,1,places=12)

    def test_nonzero_joint_scalar_threshold_preserves_full_weinberg(self):
        from Numerical.joint_scalar_matching.core import b0bar
        c=np.zeros((45,45,12,12),complex)
        q=np.zeros((12,12,12,12),float)
        W=np.array([[1,1j,0.4],[1j,2,0.1j],[0.4,0.1j,3]],complex)*1e-12
        T=a4.master_weinberg_template()
        pos=local_positions()
        # Build six *real* HHSS channels that together have the complex
        # Weinberg tensor. This synthetic input is a tensor-algebra witness,
        # not a claim about the A(-4) UV scalar quartic dictionary.
        channel=0
        from itertools import permutations
        for i,j in ((0,0),(0,1),(1,1)):
            for multiply,part in ((1.,T[i,j].real),(1j,T[i,j].imag)):
                a=4+channel
                for p in range(3):
                    for r in range(3):
                        for ii,jj in ((i,j),(j,i)):
                            c[pos[p,ii],pos[r,jj],a,a]=W[p,r]*multiply
                for A in range(4):
                    for B in range(4):
                        for inds in set(permutations((A,B,a,a))):
                            q[inds]=part[A,B]
                channel+=1
        self.assertTrue(np.allclose(q,q.swapaxes(0,1)))
        obj=a4.compute_from_arrays(c,q,1.6e6,1.6e6,pos,1e3)
        expected=-W*b0bar(1.6e6,1.6e6,1e3)/(2*16*math.pi**2)
        self.assertTrue(np.allclose(
            a4.dec(obj['finite_scalar_bubble_C5_gev_inv'],(3,3)),expected,
            rtol=3e-11,atol=1e-26))
        self.assertLess(obj['scalar_bubble_weinberg_tensor_residual'],1e-11)
        self.assertLess(obj['matching_scale_log_derivative_relative_residual'],2e-8)

    def test_zero_heavy_C_retains_direct_Weinberg_and_masses(self):
        c5=np.diag([1.,2.,3.]).astype(complex)*1e-15
        c=c_from_c5(c5)
        q=np.zeros((12,12,12,12),float)
        obj=a4.compute_from_arrays(c,q,1.2e6,2.3e6,local_positions(),1e3)
        self.assertAlmostEqual(obj['max_abs_scalar_bubble_C5_gev_inv'],0)
        self.assertTrue(np.allclose(a4.dec(obj['conditional_SMEFT_C5_gev_inv'],(3,3)),c5,rtol=1e-13,atol=1e-29))
        self.assertEqual(obj['scalar_eigenpair_decomposition']['pair_count'],36)
        self.assertTrue(obj['matching_scale_log_derivative_relative_residual']<1e-9)


if __name__=='__main__':unittest.main()
