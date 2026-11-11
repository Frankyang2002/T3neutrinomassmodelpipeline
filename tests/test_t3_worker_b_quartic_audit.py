from __future__ import annotations
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from scripts.RunT3WorkerBQuarticAudit import audit_tensors, load_npz_array, parse_pairs
from Numerical.running.T3GeneralScalarBackground import replace_source_hhss


def sym_direction():
    q=np.zeros((12,)*4)
    q[0,0,0,0]=1
    q[4,4,4,4]=2
    q[0,0,4,4]=q[0,4,0,4]=q[0,4,4,0]=q[4,0,0,4]=q[4,0,4,0]=q[4,4,0,0]=3
    return q

class Fake:
    model='A_-4'
    reference={'lambda':1.0}
    directions=(('lambda','real'),)
    reference_tensor=sym_direction()
    columns=np.array([sym_direction()])
    def full_reference(self,state):
        return self.reference_tensor + (float(state['lambda'])-1)*self.columns[0]

class Tests(unittest.TestCase):
    def test_affine_endpoint(self):
        result=audit_tensors(Fake(),{'lambda':2.})
        self.assertEqual(result['status'],'PassedAvailableTensorChecks')
        self.assertFalse(result['source_hhss_evaluated'])
        self.assertGreater(result['candidate_sector_norms']['H2S2'],0)
    def test_source_hhss_six_orbits(self):
        h=np.zeros((4,4,8,8));h[0,0,0,0]=7
        r=audit_tensors(Fake(),{'lambda':2.},hhss=h)
        self.assertEqual(r['hybrid_hhss_exact_max_abs'],0)
        self.assertEqual(r['hybrid_nonhhss_preservation_max_abs'],0)
        self.assertGreater(r['source_hhss_vs_candidate_relative'],0)
    def test_independent_direct(self):
        r=audit_tensors(Fake(),{'lambda':2.},direct=Fake().full_reference({'lambda':2.}))
        self.assertEqual(r['independent_direct_vs_cached_relative'],0.)
    def test_independent_direct_mismatch_not_silenced(self):
        r=audit_tensors(Fake(),{'lambda':2.},direct=np.zeros((12,)*4))
        self.assertGreater(r['independent_direct_vs_cached_relative'],0.)
        self.assertEqual(r['status'],'FailedTensorChecks')
    def test_unsymmetric_hhss_fails(self):
        h=np.zeros((4,4,8,8));h[0,1,0,0]=1
        with self.assertRaises(ValueError):audit_tensors(Fake(),{'lambda':2.},hhss=h)
    def test_npz_shape_and_complex_reject(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'a.npz';np.savez(p,hhss=np.zeros((4,4,8,7)))
            with self.assertRaises(ValueError):load_npz_array(p,'hhss',(4,4,8,8))
    def test_duplicate_selector_fails(self):
        with self.assertRaises(ValueError):parse_pairs(['A_-4=a','A_-4=b'],'source',{'A_-4'})

if __name__=='__main__':unittest.main()
