import unittest
import numpy as np
from scripts.RunT3WorkerBBasisSearch import (
 real_transform,candidate_rotations,rotate_hhss,fit_metrics,search)

class WorkerBBasisTests(unittest.TestCase):
 def test_rotations_count(self):
  for model in ('A_-4','D_-2'):
   rows=list(candidate_rotations(model))
   self.assertEqual(len(rows),384)
   self.assertTrue(all(np.array_equal(R.T@R,np.eye(8)) for _,R in rows))
 def test_real_complex_conjugation(self):
  R=real_transform(1,(0,),1,True)
  self.assertTrue(np.array_equal(R,np.array([[0,1],[1,0]])))
 def test_shape_rotation(self):
  x=np.zeros((4,4,8,8)); x[0,0,0,0]=1
  _,R=next(candidate_rotations('A_-4'))
  self.assertEqual(rotate_hhss(x,R).shape,x.shape)
 def test_invalid_rotation(self):
  with self.assertRaises(ValueError):real_transform(2,(0,0),0,False)
 def test_fit_identity(self):
  x=np.eye(2)
  self.assertLess(fit_metrics(x,x)['after_2x2_coupling_mix_relative_to_candidate'],1e-12)
 def test_fit_signflip(self):
  x=np.eye(2);m=fit_metrics(x,-x)
  self.assertLess(m['signflip_relative_to_candidate'],1e-12)
 def test_find_permutation(self):
  source=np.zeros((2,4,4,8,8))
  source[0,0,0,0,2]=source[0,0,0,2,0]=1
  source[1,0,0,1,3]=source[1,0,0,3,1]=1
  rotations=list(candidate_rotations('A_-4'))
  chosen=next(R for desc,R in rotations if desc['S2_permutation']==[1,0,2] and
    desc['S1_phase_i_power']==0 and desc['S2_phase_i_power']==0 and
    not desc['S1_conjugate'] and not desc['S2_conjugate'])
  candidate=np.stack([rotate_hhss(q,chosen) for q in source])
  report=search('A_-4',candidate,source,[('lambdaT3','real'),('lambdaT3','imag')],top=1)
  self.assertLess(report['best_lambdaT3']['after_2x2_coupling_mix_relative_to_candidate'],1e-10)
 def test_missing_direction(self):
  x=np.zeros((1,4,4,8,8))
  with self.assertRaises(ValueError):search('A_-4',x,x,[('lambdaT3','real')])

if __name__=='__main__':unittest.main()
