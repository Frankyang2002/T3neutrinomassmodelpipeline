import json
import unittest
from pathlib import Path

from Numerical.diagnostics.AnalyseMatcheteGenerators import analyse, analyse_representation


class MatcheteGeneratorTests(unittest.TestCase):
    def test_doublet(self):
        from sympy import I
        g1 = [['0','1/2'],['1/2','0']]
        g2 = [['1/2','0'],['0','-1/2']]
        g3 = [['0','-I/2'],['I/2','0']]
        result = analyse_representation({'dimension':2,'generators_input_form':[g1,g2,g3]})
        self.assertEqual(result['status'], 'generator_algebra_verified')
        self.assertEqual(result['cartan_spectrum_descending'], ['1/2','-1/2'])
        self.assertFalse(result['CG_tensor_covariance_verified'])

    def test_triplet(self):
        raw = [
            [['0','0','0'],['0','0','I'],['0','-I','0']],
            [['0','0','-I'],['0','0','0'],['I','0','0']],
            [['0','I','0'],['-I','0','0'],['0','0','0']],
        ]
        result = analyse_representation({'dimension':3,'generators_input_form':raw})
        self.assertEqual(result['status'], 'generator_algebra_verified')
        self.assertEqual(result['cartan_spectrum_descending'], ['1','0','-1'])

    def test_invalid_basis_is_rejected(self):
        record={'dimension':2,'generators_input_form':[[['0','0'],['0','0']]]*3}
        self.assertEqual(analyse_representation(record)['status'],'validation_failed')

if __name__ == '__main__':
    unittest.main()
