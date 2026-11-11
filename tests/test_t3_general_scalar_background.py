import unittest
import numpy as np
from Numerical.running.T3GeneralScalarBackground import (
    MODEL_DIMS, AffineTensor, construct_affine, build_scalar_background,
    replace_source_hhss,
)

class Tests(unittest.TestCase):
    def setUp(self):
        self.state={'q':0.3,'m':1000.0,'k':{'real':3.,'imag':4.}}
        self.Q=np.zeros((12,)*4)
        self.Q[4,4,4,4]=24.
        self.M=np.zeros((12,12))
        self.M[4:,4:]=np.eye(8)*1e6
        self.M[4,4]=1e6+6.; self.M[5,5]=1e6-6.; self.M[4,5]=self.M[5,4]=-8.
        self.q=AffineTensor(self.state,(('q','real'),),self.Q,np.zeros((1,)+(12,)*4),4)
        self.m=AffineTensor(self.state,(('m','real'),('k','real'),('k','imag')),self.M,np.zeros((3,12,12)),2)

    def test_16_model_count(self): self.assertEqual(len(MODEL_DIMS),16)
    def test_mass_eigenbasis(self):
        b=build_scalar_background('D_-2',1000,self.state,self.q,self.m,
                                  expected_quartic_names={'q'},expected_quadratic_names={'m','k'})
        self.assertEqual(b.q_fixed.shape,(12,)*4)
        self.assertTrue(np.allclose(b.scalar_rotation.T@b.scalar_mass_matrix@b.scalar_rotation,
                                    np.diag(b.eigenvalues_gev2)))
        self.assertTrue(np.allclose(b.eigenvalues_gev2[:2],[1e6-10.,1e6+10.]))
    def test_mixed_requires_basis(self):
        with self.assertRaisesRegex(ValueError,'dedicated basis_rotation'):
            build_scalar_background('B_-1',1000,self.state,self.q,self.m,
                 expected_quartic_names={'q'},expected_quadratic_names={'m','k'})
        build_scalar_background('B_-1',1000,self.state,self.q,self.m,
             expected_quartic_names={'q'},expected_quadratic_names={'m','k'},basis_rotation=np.eye(8))
    def test_missing_inventory(self):
        with self.assertRaisesRegex(ValueError,'inventory mismatch'):
            build_scalar_background('A_-4',1000,self.state,self.q,self.m,
                 expected_quartic_names={'q','additional'},expected_quadratic_names={'m','k'})
    def test_source_hhss(self):
        h=np.ones((4,4,8,8))*.4
        q=replace_source_hhss(self.Q,h)
        self.assertEqual(q[1,4,3,5],.4)
        self.assertEqual(q[4,4,4,4],24.)
    def test_affine_complex(self):
        ref={'z':{'real':0.,'imag':0.}}
        def direct(x):
            t=np.zeros((12,12)); t[4,4]=2*x['z']['real']-3*x['z']['imag']; return t
        response=construct_affine(ref,('z',),frozenset({'z'}),direct,2)
        self.assertAlmostEqual(response.evaluate({'z':{'real':2.,'imag':1.}})[4,4],1.)

if __name__=='__main__':unittest.main()
