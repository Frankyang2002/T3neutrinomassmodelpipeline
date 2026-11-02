import unittest
from Numerical.direct_one_loop.comparison import compare_same_point
from Numerical.direct_one_loop.core import encoded

class ComparisonTests(unittest.TestCase):
    def test_same_scale_and_mismatch_guard(self):
        e=encoded([[1e-11,0,0],[0,2e-11,0],[0,0,4e-11]])
        direct={'source_provenance':{'benchmark_sha256':'hash'},
                'renormalization_scale_gev':1e7,
                'total':{'complex_symmetric_matrix_gev':e},
                'mass_hierarchy':{'max_scalar_squared_to_fermion_squared_ratio':1e-4}}
        ref={'physical_dimension':'mass_matrix_gev','scheme':'synthetic_test',
             'benchmark_sha256':'hash','renormalization_scale_gev':1e7,
             'same_renormalized_uv_parameters':True,'matrix_gev':e}
        r=compare_same_point(direct,ref)
        self.assertAlmostEqual(r['relative_mass_matrix_frobenius_difference'],0)
        ref['renormalization_scale_gev']=100
        self.assertEqual(compare_same_point(direct,ref)['status'],'not_compared_different_scales')
        ref['benchmark_sha256']='wrong'
        with self.assertRaises(ValueError):compare_same_point(direct,ref)

if __name__=='__main__':unittest.main()
