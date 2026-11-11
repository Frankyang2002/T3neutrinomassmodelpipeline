import unittest
import numpy as np
from scripts.RunT3WorkerBHHSSDiscrepancy import diagnose, _metrics, _pair_args


class Tests(unittest.TestCase):
    def test_identical(self):
        s=np.ones((4,4,8,8))
        a=diagnose('A_-4',s,s,top_n=1)
        self.assertAlmostEqual(a['overall']['relative_to_source'],0.)
        self.assertAlmostEqual(a['overall']['least_squares_candidate_over_source'],1.)
        self.assertEqual(a['largest_absolute_entry_differences'][0]['absolute_difference'],0.)

    def test_negative(self):
        s=np.ones((4,4,8,8)); a=diagnose('D_-2',-s,s,top_n=0)
        self.assertAlmostEqual(a['overall']['relative_to_source'],2.)
        self.assertAlmostEqual(a['overall']['least_squares_candidate_over_source'],-1.)
        self.assertAlmostEqual(a['overall']['residual_if_candidate_sign_flipped_relative_to_source'],0.)
        self.assertEqual(a['largest_absolute_entry_differences'],[])

    def test_scale(self):
        s=np.zeros((4,4,8,8));s[0,0,0,0]=5.
        a=diagnose('A_-4',3*s,s)
        self.assertAlmostEqual(a['overall']['least_squares_candidate_over_source'],3.)
        self.assertAlmostEqual(a['overall']['residual_after_best_real_scale_relative_to_source'],0.)

    def test_shape_mismatch_residual(self):
        s=np.zeros((4,4,8,8));s[0,0,0,0]=1.
        c=s.copy();c[0,0,7,7]=1.
        a=diagnose('A_-4',c,s)
        self.assertAlmostEqual(a['overall']['residual_after_best_real_scale_relative_to_source'],1.)
        self.assertEqual(a['largest_absolute_entry_differences'][0]['heavy_species'],['S2','S2'])

    def test_zero_source_block(self):
        c=np.ones((2,2));s=np.zeros((2,2))
        self.assertEqual(_metrics(c,s)['status'],'ZeroSourceNorm')
        self.assertIsNone(_metrics(c,s)['relative_to_source'])

    def test_block_partition(self):
        s=np.ones((4,4,8,8))
        a=diagnose('D_-2',s,s)
        self.assertEqual(set(a['heavy_species_blocks']),{'S1S1','S1S2','S2S1','S2S2'})
        self.assertEqual(len(a['higgs_component_blocks']),10)

    def test_invalid(self):
        with self.assertRaises(ValueError): diagnose('B_-1',np.zeros((4,4,8,8)),np.zeros((4,4,8,8)))
        with self.assertRaises(ValueError): diagnose('A_-4',np.zeros((4,4,8,8)),np.zeros((4,4,7,8)))
        with self.assertRaises(ValueError): _pair_args(['A_-4=foo','A_-4=bar'],('A_-4',),'input')


if __name__=='__main__':unittest.main()
