import unittest
import numpy as np

from scripts.RunT3WorkerBHHSSDirectionAudit import (
    metrics, number, put, check_symmetry, audit_responses,
)


def basis(k, value=1.):
    q=np.zeros((4,4,8,8),float)
    q[0,0,k,k]=value
    return q


class DirectionTests(unittest.TestCase):
    def test_complex_conversion(self):
        self.assertEqual(number({'real':2,'imag':3}),2+3j)
        self.assertEqual(put({'x':0},'x',1+2j,True)['x'],{'real':1.,'imag':2.})

    def test_identical(self):
        m=metrics(basis(0),basis(0))
        self.assertAlmostEqual(m['candidate_over_source_scale'],1)
        self.assertAlmostEqual(m['shape_residual_relative_to_source'],0)

    def test_global_sign(self):
        m=metrics(-basis(0),basis(0))
        self.assertAlmostEqual(m['candidate_over_source_scale'],-1)
        self.assertAlmostEqual(m['shape_residual_relative_to_source'],0)

    def test_shape_mismatch(self):
        m=metrics(basis(1),basis(0))
        self.assertAlmostEqual(m['candidate_over_source_scale'],0)
        self.assertAlmostEqual(m['shape_residual_relative_to_source'],1)

    def test_source_only_and_mixed_directions(self):
        state={'lambdaA':0.3,'lambdaB':{'real':0.2,'imag':-0.1},'lambdaC':0.4}
        def direct(x):
            a=number(x['lambdaA']);b=number(x['lambdaB']);c=number(x['lambdaC'])
            return (1+float(a.real))*basis(0)+float(b.real)*basis(1)+float(b.imag)*basis(2)+float(c.real)*basis(3)
        candidate=direct(state)+2*basis(1)
        direction=(('lambdaA','real'),('lambdaB','real'),('lambdaB','imag'))
        cols=np.stack((basis(0),-basis(1),basis(2)))
        out=audit_responses(state,candidate,direction,cols,direct)
        self.assertEqual(out['source_only_directions'],[['lambdaC','real']])
        self.assertEqual(out['source_evaluated_direction_count'],4)
        records={r['coupling']+'.'+r['axis']:r for r in out['directions']}
        self.assertAlmostEqual(records['lambdaB.real']['candidate_over_source_scale'],-1.)
        self.assertTrue(records['lambdaC.real']['candidate_direction_absent'])
        self.assertLess(out['source_affine_interpolation_relative'],1e-12)

    def test_reject_nonlinear(self):
        state={'lambdaA':0.4}
        def direct(x):return float(number(x['lambdaA']).real)**2*basis(0)
        with self.assertRaisesRegex(ValueError,'not affine'):
            audit_responses(state,direct(state),(('lambdaA','real'),),np.stack((basis(0),)),direct)

    def test_reject_nonsymmetric(self):
        q=basis(0);q[1,0,1,1]=1
        with self.assertRaisesRegex(ValueError,'Bose'):
            check_symmetry(q,'test')


if __name__=='__main__':unittest.main()
