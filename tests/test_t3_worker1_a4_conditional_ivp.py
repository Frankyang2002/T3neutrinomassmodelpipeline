import unittest
import numpy as np
from scripts.RunT3Worker1A4ConditionalIVP import compare_matrices, report_pair, matrix_stage_fractions


def enc(q):
    return {'real':np.asarray(q).real.tolist(),'imag':np.asarray(q).imag.tolist()}


def case(q):
    return {'c5_at_low_scale':enc(q),'c5_at_scalar_threshold':enc(q),
            'masses_ev_ascending':[0.01,0.02,0.05],
            'prediction':{'sin2_theta12':0.3,'sin2_theta13':0.02,'sin2_theta23':0.5,
                          'delta_m21_sq_ev2':7e-5,'delta_m3l_sq_ev2':2.5e-3}}

class Tests(unittest.TestCase):
    def test_complex_proportional(self):
        q=np.diag([1,2,3]).astype(complex)
        x=compare_matrices((2+3j)*q,q)
        self.assertLess(x['residual_after_complex_rescaling_relative_to_alternative'],1e-12)
        self.assertGreater(x['residual_after_real_rescaling_relative_to_alternative'],0.5)
    def test_shape_change(self):
        q=np.diag([1,2,3]).astype(complex);r=np.diag([1,4,3]).astype(complex)
        self.assertGreater(compare_matrices(r,q)['residual_after_complex_rescaling_relative_to_alternative'],0.1)
    def test_stage_fractions(self):
        q=np.eye(3)
        cases={key:case(q.copy()) for key in ('1111','0111','1001','1110','1011','1101')}
        self.assertEqual(matrix_stage_fractions(cases)['UV_off'],0)
        alternative={'model':'A_-4','status':'Computed','cases':cases}
        result=report_pair(alternative,alternative)
        self.assertEqual(result['comparison']['1111']['low_energy_c5']['relative_difference_to_baseline'],0)
        self.assertFalse(result['physical_alternative_certified'])
    def test_reject_noncompleted(self):
        with self.assertRaises(ValueError):report_pair({'model':'A_-4','status':'Blocked'}, {'model':'A_-4','status':'Computed'})

if __name__=='__main__':unittest.main()
