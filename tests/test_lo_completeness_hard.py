"""One-loop matching algebra, all-heavy-fermions/common scalar-threshold tests."""
import math
import unittest
import numpy as np
from Numerical.lo_certification.completeness.hard_matching import (
    b0bar,triangle_ratio,weinberg_tensor,evaluate_hard_scalar_channel)


def example(n=8, masses=(300.,500.), scalar=20.):
    k=len(masses)
    c=np.zeros((k,3,3,2,2,n,n),complex)
    q=np.zeros((4,4,n,n),float)
    t=weinberg_tensor()
    flavour=np.array([[1+2j,.8-.7j,.1j],[.8-.7j,2-1j,-1+.3j],
                      [.1j,-1+.3j,.2+1.1j]],complex)*1e-7
    for a in range(4):
        for b in range(4):
            q[a,b,a,b]+=0.5
            q[b,a,a,b]+=0.5
            for j in range(k):
                c[j,:,:,:,:,a,b]=flavour[:,:,None,None]*t[:,:,a,b][None,None,:,:]/(j+1)
    # Source has two separate species, never substitute C_total for C_k.
    s=np.eye(n)*scalar**2
    return c,np.array(masses),q,s


class HardMatchingTests(unittest.TestCase):
    def test_b0_equal_and_near_degenerate(self):
        self.assertAlmostEqual(b0bar(400.,400.,20.),0.,places=13)
        self.assertAlmostEqual(b0bar(400*(1+1e-8),400.,20.),-0.5e-8,places=11)

    def test_triangle_degenerate_three_masses(self):
        self.assertAlmostEqual(triangle_ratio(100.,100.,100.),-0.5,places=13)
        self.assertAlmostEqual(triangle_ratio(100.,100.*(1+1e-9),100.),-0.5,places=7)

    def test_small_scalar_expansion_is_one_loop_not_two(self):
        m2=400.**2
        for a,b in [(1.,1.),(1.,1.1),(10.,9.5)]:
            full=triangle_ratio(m2,a,b)
            leading=(1.+math.log(m2/400.**2))-b0bar(a,b,400.)
            self.assertLess(abs(full-leading),0.003)
        self.assertGreater(abs(triangle_ratio(m2,100.,100.)-
                               (1.-b0bar(100.,100.,400.))),1e-4)

    def test_complex_flavour_and_multi_fermion_consistency(self):
        c,m,q,s=example()
        x=evaluate_hard_scalar_channel(c,m,q,s,300,20,C_total_muF=c.sum(0))
        self.assertTrue(x['projection_passed'],x['residuals'])
        self.assertEqual(len(x['per_pair_terms']),len(m)*8*9//2)
        for a in x['matrices'].values():
            self.assertEqual(a.shape,(3,3))
            self.assertLess(np.linalg.norm(a-a.T),1e-18)
        self.assertGreater(abs(x['matrices']['dim5_sum'][0,1].imag),0)
        self.assertFalse(x['physical_C5_certified'])
        self.assertEqual(x['higher_dimension_remainder_loop_order'],1)

    def test_common_muF_muS_log_partition(self):
        c,m,q,s=example()
        x=evaluate_hard_scalar_channel(c,m,q,s,300,20)
        y=evaluate_hard_scalar_channel(c,m,q,s,700,8)
        np.testing.assert_allclose(x['matrices']['dim5_sum'],y['matrices']['dim5_sum'],rtol=1e-12,atol=1e-19)
        self.assertGreater(np.linalg.norm(x['matrices']['hard_muF']-y['matrices']['hard_muF']),0)
        self.assertLess(np.linalg.norm(x['matrices']['higher_dimension_remainder']),
                        0.05*np.linalg.norm(x['matrices']['full_triangle']))

    def test_phase_covariance_no_hermitian_conjugation(self):
        c,m,q,s=example()
        phase=np.exp(1j*0.73)
        x=evaluate_hard_scalar_channel(c,m,q,s,300,20)
        y=evaluate_hard_scalar_channel(c*phase,m,q,s,300,20)
        for key in x['matrices']:
            np.testing.assert_allclose(y['matrices'][key],phase*x['matrices'][key],rtol=1e-11,atol=1e-18)

    def test_per_species_sum_and_input_rejection(self):
        c,m,q,s=example()
        with self.assertRaisesRegex(ValueError,'do not sum'):
            evaluate_hard_scalar_channel(c,m,q,s,300,20,C_total_muF=2*c.sum(0))
        with self.assertRaisesRegex(ValueError,'Tachyonic'):
            bad=s.copy();bad[0,0]=-1
            evaluate_hard_scalar_channel(c,m,q,bad,300,20)
        with self.assertRaisesRegex(ValueError,'HHSS tensor symmetry'):
            bad=q.copy();bad[0,1,0,0]+=1
            evaluate_hard_scalar_channel(c,m,bad,s,300,20)

    def test_scalar_pair_label_not_sequential_threshold(self):
        c,m,q,s=example()
        # All eight real fields belong to a single diagonalisation. An arbitrary
        # cross-multiplet orthogonal rotation is allowed if every tensor rotates.
        u=np.eye(8)
        u[0,4]=1/math.sqrt(2);u[0,0]=1/math.sqrt(2)
        u[4,0]=-1/math.sqrt(2);u[4,4]=1/math.sqrt(2)
        cr=np.einsum('kfgijab,ac,bd->kfgijcd',c,u,u,optimize=True)
        qr=np.einsum('ABab,ac,bd->ABcd',q,u,u,optimize=True)
        sr=u.T@s@u
        x=evaluate_hard_scalar_channel(c,m,q,s,300,20)
        y=evaluate_hard_scalar_channel(cr,m,qr,sr,300,20)
        np.testing.assert_allclose(x['matrices']['dim5_sum'],y['matrices']['dim5_sum'],rtol=1e-11,atol=1e-18)


if __name__=='__main__': unittest.main()

class ChannelCoverageTests(unittest.TestCase):
    def test_split_quartic_channels_preserve_full_complex_sum(self):
        c,m,q,s=example(n=12)
        src=np.stack([q*.003,q*.997])
        r=evaluate_hard_scalar_channel(c,m,q,s,300,20,
                    lambda_channels=src,channel_names=('suppressed','dominant'))
        self.assertEqual(len(r['quartic_channels']),2)
        for key,full in r['matrices'].items():
            np.testing.assert_allclose(full,sum(v['matrices'][key] for v in
                 r['quartic_channels'].values()),rtol=1e-10,atol=1e-18)
        self.assertGreater(abs(r['quartic_channels']['suppressed']['matrices']['hard_muF'][0,1]),0)
        with self.assertRaisesRegex(ValueError,'exactly sum'):
            evaluate_hard_scalar_channel(c,m,q,s,300,20,
                     lambda_channels=src*2,channel_names=('x','y'))
