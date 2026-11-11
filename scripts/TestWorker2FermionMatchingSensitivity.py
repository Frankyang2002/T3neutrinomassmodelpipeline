"""Unit tests independent of optional Matchete/Wolfram installations."""
import sys, unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
from Worker2FermionMatchingSensitivity import symmetric,directions,observables

class Worker2Tests(unittest.TestCase):
    def setUp(self):
        rng=np.random.default_rng(12345)
        self.C=symmetric(rng.normal(size=(3,3,2,2,4,4))+1j*rng.normal(size=(3,3,2,2,4,4)))
    def test_exchange_and_orthogonality(self):
        C=self.C
        for name,D in directions(C).items():
            self.assertLess(np.linalg.norm(D-D.transpose(1,0,3,2,4,5)),1e-10)
            self.assertLess(np.linalg.norm(D-D.swapaxes(4,5)),1e-10)
            self.assertAlmostEqual(np.linalg.norm(D)/np.linalg.norm(C),1.,places=12)
            if name!='aligned':self.assertLess(abs(np.vdot(C,D)/np.vdot(C,C)),1e-10)
    def test_signed_relative_norm(self):
        C=self.C
        for D in directions(C).values():
            for p in (.01,.05,.10):
                self.assertAlmostEqual(np.linalg.norm(p*D)/np.linalg.norm(C),p,places=12)
    def test_neutrino_observables(self):
        M=np.diag([.01,.02,.05]).astype(complex)
        o=observables(M)
        self.assertAlmostEqual(o['mass_ratio_m2_m3'],.4)
        self.assertAlmostEqual(o['theta13_deg'],0)
        self.assertAlmostEqual(o['theta12_deg'],0)
    def test_reject_nonsymmetric(self):
        M=np.eye(3,dtype=complex);M[0,1]=.3
        with self.assertRaises(ValueError):observables(M)

if __name__=='__main__':unittest.main()
