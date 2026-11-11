import unittest
from types import SimpleNamespace
import numpy as np
from scripts.RunT3WorkerBSourceHHSSAudit import fixed_native_hhss, _pairs
from scripts.RunT3WorkerBQuarticAudit import audit_tensors

class SourceHHSSAuditTests(unittest.TestCase):
    def test_identity(self):
        h=np.zeros((4,4,8,8));h[0,0,2,2]=3
        self.assertTrue(np.array_equal(fixed_native_hhss(h),h))
    def test_rotation_convention(self):
        r=np.eye(8);r[0,0]=r[1,1]=0;r[0,1]=1;r[1,0]=-1
        h=np.zeros((4,4,8,8));h[0,0,0,0]=3
        n=fixed_native_hhss(h,r)
        self.assertEqual(n[0,0,1,1],3)
        back=np.einsum('ia,ABij,jb->ABab',r,n,r,optimize=True)
        self.assertTrue(np.allclose(back,h))
    def test_nonorthogonal_rejected(self):
        with self.assertRaisesRegex(ValueError,'orthogonal'):
            fixed_native_hhss(np.zeros((4,4,8,8)),2*np.eye(8))
    def test_wrong_shape_rejected(self):
        with self.assertRaisesRegex(ValueError,'finite'):
            fixed_native_hhss(np.zeros((4,4,7,8)))
    def test_complex_rotation_permutation(self):
        r=np.eye(8);r[0,0]=r[1,1]=0;r[0,1]=1;r[1,0]=-1
        h=np.zeros((4,4,8,8));h[1,2,0,1]=2;h[2,1,0,1]=2;h[1,2,1,0]=2;h[2,1,1,0]=2
        n=fixed_native_hhss(h,r)
        self.assertTrue(np.allclose(n,n.transpose(1,0,2,3)))
        self.assertTrue(np.allclose(n,n.transpose(0,1,3,2)))
    def test_bad_pairs(self):
        with self.assertRaises(ValueError):_pairs(['D_-2=a','D_-2=b'],('A_-4','D_-2'))
    def test_nonhhss_preserved_and_discrepancy_reported(self):
        q=np.zeros((12,)*4);q[5,5,5,5]=1
        context=SimpleNamespace(model='A_-4',reference_tensor=q,columns=np.zeros((1,)+(12,)*4),
                                directions=(('lambdaH','real'),),full_reference=lambda state:q.copy())
        h=np.zeros((4,4,8,8));h[0,0,0,0]=2
        record=audit_tensors(context,{'lambdaH':0.1},hhss=h)
        self.assertEqual(record['hybrid_nonhhss_preservation_max_abs'],0.)
        self.assertEqual(record['hybrid_hhss_exact_max_abs'],0.)
        self.assertGreater(record['source_hhss_vs_candidate_relative'],0.)
        self.assertEqual(record['status'],'PassedAvailableTensorChecks')
if __name__=='__main__':unittest.main()
