"""Offline tests for A_-4 native SU(2) contraction diagnostic."""
import importlib.util
from pathlib import Path
import unittest

P = Path(__file__).with_name('check_native_cg_a_m4.py')
spec = importlib.util.spec_from_file_location('cgprobe', P)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


class NativeCGTests(unittest.TestCase):
    def test_delta_internal_sum(self):
        y1 = [1,0,0,1]
        y2 = [complex(i) for i in range(12)]
        expected = [y2[(j*2+i)*3+u] for i in range(2) for j in range(2) for u in range(3)]
        self.assertEqual(m.product_yukawa_cgs(y1,y2),expected)

    def test_bar_yukawa_cgs(self):
        y1=[1j,0,0,1]
        y2=[complex(0,1)]*12
        p=m.product_yukawa_cgs(y1,y2)
        self.assertEqual(p[0],-1)
        self.assertEqual(p[6],-1j)

    def test_exact_unit_identity(self):
        x=[complex(i) for i in range(12)]
        r=m.proportionality(x,x)
        self.assertTrue(r['unity_ratio'])
        self.assertLess(r['residual_relative'],1e-12)

    def test_reversed_sign_not_certified(self):
        x=[complex(i) for i in range(12)]
        r=m.proportionality([-a for a in x],x)
        self.assertFalse(r['unity_ratio'])
        self.assertTrue(r['minus_unity_ratio'])

    def test_unrelated_not_proportional(self):
        r=m.proportionality([1,0,0],[0,1,0])
        self.assertFalse(r['proportional'])

    def test_wrong_dim_fail(self):
        with self.assertRaises(ValueError):
            m.parse_cg({'name':'T3Y1CG','dimensions':[2,3]},'T3Y1CG',[2,2])

    def test_source_claim_fails_closed(self):
        with self.assertRaises(ValueError):
            m.analyze({'format':'t3-native-a-m4-yukawa-mixing-cg-evidence-v1','model':'A_-4',
                       'status':'NativeCGTensorsExtracted_NotPhaseCertified',
                       'previous_release_changed':False,'cg_product_identity_computed':True})


if __name__=='__main__':
    unittest.main(verbosity=2)
