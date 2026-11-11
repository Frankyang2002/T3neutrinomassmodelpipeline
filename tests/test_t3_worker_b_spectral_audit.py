import unittest
import numpy as np
from scripts.RunT3WorkerBSpectralAudit import (
    check_hhss,compare_fingerprints,rotate_hhss,scan_source_phase,
    species_sector_fingerprints,diagnose,
)


def example(seed=3):
    rng=np.random.default_rng(seed)
    q=rng.normal(size=(4,4,8,8))
    return (q+q.swapaxes(0,1)+q.swapaxes(2,3)+q.transpose(1,0,3,2))/4


class Tests(unittest.TestCase):
    def test_basis_invariance(self):
        q=example()
        r=np.linalg.qr(np.random.default_rng(5).normal(size=(8,8)))[0]
        h=np.linalg.qr(np.random.default_rng(6).normal(size=(4,4)))[0]
        t=rotate_hhss(q,r,h)
        result=compare_fingerprints(q,t)
        self.assertLess(result['normalized_spectral_lower_bound'],1e-13)
        self.assertFalse(result['exact_orthogonal_rotation_plus_scale_excluded'])

    def test_normalization_invariance(self):
        q=example()
        self.assertLess(compare_fingerprints(q,-2*q)['normalized_spectral_lower_bound'],1e-14)

    def test_real_spectral_obstruction(self):
        q=np.zeros((4,4,8,8));q[0,0,0,0]=1.
        w=np.zeros_like(q);w[0,0,0,0]=1.;w[1,1,1,1]=1.
        report=compare_fingerprints(q,w)
        self.assertTrue(report['exact_orthogonal_rotation_plus_scale_excluded'])
        self.assertGreater(report['normalized_spectral_lower_bound'],0.5)

    def test_reject_bad_symmetry(self):
        q=np.zeros((4,4,8,8));q[0,1,0,0]=1.
        with self.assertRaises(ValueError):check_hhss(q)

    def test_phase_scan(self):
        r=example(4);i=example(5)
        c=2*(np.cos(np.pi/4)*r+np.sin(np.pi/4)*i)
        result=scan_source_phase(c,r,i,samples=32)
        self.assertLess(result['best_sampled_phase']['spectral_distance'],1e-13)

    def test_species_blocks(self):
        q=np.zeros((4,4,8,8));q[0,0,3,3]=2.
        a=species_sector_fingerprints(q,'A_-4')
        d=species_sector_fingerprints(q,'D_-2')
        self.assertEqual(a['S2S2']['norm'],2)
        self.assertEqual(d['S1S1']['norm'],2)

    def test_missing_direction_fails(self):
        q=np.zeros((1,4,4,8,8))
        with self.assertRaises(ValueError):diagnose('A_-4',[('lambdaT3','real')],q,q)

    def test_diagnose(self):
        r=example();i=example(4);h=example(5)
        columns=np.stack([r,i,h])
        x=diagnose('D_-2',[('lambdaT3','real'),('lambdaT3','imag'),('lambdaH2','real')],-columns,columns,samples=16)
        self.assertLess(x['directions']['lambdaH2.real']['signflip_rel_to_candidate'],1e-13)
        self.assertEqual(x['direction_count'],3)


if __name__=='__main__':unittest.main()
