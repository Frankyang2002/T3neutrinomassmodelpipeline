"""A_-4 polynomial quotient tests with ordered 3-flavour real scalar fields."""
from __future__ import annotations
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np

from Numerical.running.GeneralT3WilsonPolynomialBridge import (
    swap_fermion, swap_scalar, symmetrise_polynomial,
    polynomial_residual, monomial_coefficients, bridge_a4,
    master_weinberg_polynomial_check, a4_lepton_dummy_alignment,
)


class ProjectionTests(unittest.TestCase):
    def test_complex_cross_block_source_preserved(self):
        inp={(1,1,2,3,5,7):3+5j,(1,2,2,1,6,12):1-3j}
        out=symmetrise_polynomial(inp)
        self.assertTrue(polynomial_residual(inp,out)['preserved_to_numerical_precision'])
        self.assertEqual(out[(1,1,2,3,5,7)],(3+5j)/4)
        self.assertEqual(out[swap_scalar((1,1,2,3,5,7))],(3+5j)/4)
        self.assertEqual(out[swap_fermion((1,1,2,3,5,7))],(3+5j)/4)

    def test_fully_symmetric_source_unchanged(self):
        k=(1,1,2,3,5,7)
        source={x:2-1j for x in (k,swap_scalar(k),swap_fermion(k),swap_scalar(swap_fermion(k)))}
        self.assertEqual(symmetrise_polynomial(source),source)

    def test_scalar_antisymmetric_null_polynomial(self):
        k=(1,1,2,2,5,7)
        inp={k:5j,swap_scalar(k):-5j}
        self.assertEqual(symmetrise_polynomial(inp),{})
        self.assertEqual(monomial_coefficients(inp),{})

    def test_different_flavours_not_collapsed(self):
        inp={(1,1,2,3,5,7):1,(1,1,2,2,5,7):2}
        self.assertEqual(len(monomial_coefficients(inp)),2)
        out=symmetrise_polynomial(inp)
        self.assertEqual(len(monomial_coefficients(out)),2)
        self.assertTrue(polynomial_residual(inp,out)['preserved_to_numerical_precision'])

    def test_same_external_fermion_factor_two_orbit(self):
        k=(1,1,1,1,5,7)
        inp={k:3+2j}
        out=symmetrise_polynomial(inp)
        self.assertEqual(out[k],(3+2j)/2)
        self.assertEqual(out[swap_scalar(k)],(3+2j)/2)
        self.assertTrue(polynomial_residual(inp,out)['preserved_to_numerical_precision'])

    def test_wrong_scalar_pair_rejected(self):
        with self.assertRaisesRegex(ValueError,'wrong scalar-field pair'):
            symmetrise_polynomial({(1,1,2,2,5,6):1})

    def test_nonfinite_rejected(self):
        with self.assertRaisesRegex(ValueError,'Nonfinite'):
            symmetrise_polynomial({(1,1,2,2,5,7):complex('nan')})

    def test_random_ordered_tensor_monomial_identity(self):
        rng=np.random.default_rng(2026)
        t={}
        for i in (1,2):
            for j in (1,2):
                for p in (1,2,3):
                    for q in (1,2,3):
                        for a in (5,6):
                            for b in range(7,13):
                                t[(i,p,j,q,a,b)] = rng.normal()+1j*rng.normal()
        s=symmetrise_polynomial(t)
        self.assertTrue(polynomial_residual(t,s)['preserved_to_numerical_precision'])
        raw=np.linalg.norm(list(t.values()))
        after=np.linalg.norm(list(s.values()))
        self.assertGreater(after/raw,0.0)
        self.assertLess(after/raw,1.0)

    def test_fermion_symmetric_cross_scalar_norm_ratio(self):
        # The observed A_-4 audit has F-exchange residual exactly zero.
        # For support on only one scalar-block ordering the projection norm
        # falls by sqrt(2), without changing the polynomial.
        source={}
        for p in (1,2,3):
            for q in (1,2,3):
                for a in (5,6):
                    for b in (7,8,9,10,11,12):
                        k=(1,p,2,q,a,b)
                        val=complex(p+q+a, b-p)
                        source[k]=val
                        source[swap_fermion(k)]=val
        projected=symmetrise_polynomial(source)
        self.assertTrue(polynomial_residual(source,projected)['preserved_to_numerical_precision'])
        self.assertAlmostEqual(
            np.linalg.norm(list(projected.values()))/np.linalg.norm(list(source.values())),
            1/(2**0.5), places=12)

    def test_master_weinberg_tensor_convention_when_available(self):
        try:
            import RGE.matching.WeinbergTensorAdapter
        except ModuleNotFoundError:
            self.skipTest('Master RGE source unavailable outside project repository')
        evidence=master_weinberg_polynomial_check()
        self.assertTrue(evidence['master_tensor_reconstructs_polynomial_without_extra_factor'])
        self.assertFalse(evidence['matchete_spinor_normalization_certified'])


def mock_source_boundary():
    return {
        'format':'t3-source-exact-F-tree-flavour-operator-handoff-v1',
        'model':'A_-4',
        'status':'SourceFactoredOrderedFlavourMatricesAtF_NotCanonicalWilsonTensor',
        'operator_terms':[{
            'source_index':1,
            'flavour_matrix_axes_external_indices':['d1','d2'],
            'factors':[{'external_index':'d1'},{'external_index':'d2'}],
            'operator_template_unit_flavor_inputform':(
                'NCM[Transp[Field[l, Fermion, '
                '{Index[d2, SU2L[fund]], Index[d2, Flavor]}, {}]],'
                'DiracProduct[GammaCC,Proj[-1]],'
                'Field[l, Fermion, '
                '{Index[d1, SU2L[fund]], Index[d1, Flavor]}, {}]]')
        }]
    }


class RunnerTests(unittest.TestCase):
    def test_bridge_fixture_preserves_source(self):
        from Numerical.running.GeneralT3OrderedTensor import scalar_basis
        sb=scalar_basis('A_-4')
        coeff=np.eye(3)
        payload={'model':'A_-4',
            'format':'t3-ordered-source-flavour-and-real-scalar-tensor-v1',
            'status':'OrderedNumericalRealComponentTensorFromMatchete',
            'mu_F_gev':1e5,'scalar_basis_explicit':True,
            'matchete_cg_contractions_retained':True,
            'external_three_family_flavour_order_not_symmetrized':True,
            'canonical_C_ijab_normalization_certified':False,
            'tree_operator_count_PL':1,'total_ordered_real_component_count':1,
            'operators':[{'source_index':1,
                'basis_scalar_fields':dict(sb.dimensions),
                'real_scalar_dimension_including_higgs':sb.total_dimension,
                'ordered_not_canonical_symmetrized':True,
                'flavour_matrix_axes_external_indices':['d1','d2'],
                'component_count':1,
                'components':[{'lepton_su2':[1,2],'scalar_real':[5,7],
                    'ordered_flavour_coefficient_gev_inv':{
                        'real':coeff.tolist(),'imag':np.zeros((3,3)).tolist()}}]}]}
        snapshot=json.dumps(payload,sort_keys=True)
        report=bridge_a4(payload,mock_source_boundary(),check_master=False)
        self.assertEqual(snapshot,json.dumps(payload,sort_keys=True))
        self.assertTrue(report['polynomial_identity']['preserved_to_numerical_precision'])
        self.assertEqual(report['source_ordered_coordinate_count'],9)
        self.assertFalse(report['matchete_to_master_overall_matching_factor_certified'])
        self.assertEqual(report['lepton_dummy_alignment']['flavour_axis_for_lepton_field_order_zero_based'],[1,0])
        # The original matrix has only diagonal flavour entries, which mask
        # the flipped SU2-vs-flavour ordering in a numerical benchmark.
        expected=(1,1,2,1,5,7)
        self.assertIn(expected,[tuple(v['index']) for v in report['symmetric_polynomial_coordinates']])

    def test_offdiagonal_flavour_axes_are_swapped(self):
        a=mock_source_boundary()
        operator={'source_index':1,'flavour_matrix_axes_external_indices':['d1','d2']}
        parsed=a4_lepton_dummy_alignment({'operators':[operator]},a)
        self.assertEqual(parsed['flavour_axis_for_lepton_field_order_zero_based'],[1,0])
        self.assertTrue(parsed['tensor_su2_and_flavour_orders_differ'])

    def test_bridge_rekeys_offdiagonal_flavour_matrix(self):
        from Numerical.running.GeneralT3OrderedTensor import scalar_basis
        sb=scalar_basis('A_-4')
        matrix=np.zeros((3,3));matrix[0,1]=4.0
        payload={
            'model':'A_-4','format':'t3-ordered-source-flavour-and-real-scalar-tensor-v1',
            'status':'OrderedNumericalRealComponentTensorFromMatchete','mu_F_gev':1e5,
            'scalar_basis_explicit':True,'matchete_cg_contractions_retained':True,
            'external_three_family_flavour_order_not_symmetrized':True,
            'canonical_C_ijab_normalization_certified':False,
            'tree_operator_count_PL':1,'total_ordered_real_component_count':1,
            'operators':[{'source_index':1,'basis_scalar_fields':dict(sb.dimensions),
                          'real_scalar_dimension_including_higgs':sb.total_dimension,
                          'ordered_not_canonical_symmetrized':True,
                          'flavour_matrix_axes_external_indices':['d1','d2'],
                          'component_count':1,
                          'components':[{'lepton_su2':[1,2],'scalar_real':[5,7],
                             'ordered_flavour_coefficient_gev_inv':{
                                  'real':matrix.tolist(),'imag':np.zeros((3,3)).tolist()}}]}]}
        out=bridge_a4(payload,mock_source_boundary(),check_master=False)
        entries={tuple(row['index']):complex(row['real'],row['imag'])
                 for row in out['symmetric_polynomial_coordinates']}
        # SU2 index 1 belongs to the first Matchete lepton d2, whose
        # flavour comes from source matrix COLUMN 2; the other is ROW 1.
        self.assertEqual(entries[(1,2,2,1,5,7)],1.0+0j)
        self.assertNotIn((1,1,2,2,5,7),entries)
        self.assertTrue(out['polynomial_identity']['preserved_to_numerical_precision'])

    def test_reject_mismatched_flavour_yukawa_axis(self):
        a=mock_source_boundary()
        a['operator_terms'][0]['factors'][0]['external_index']='d2'
        with self.assertRaisesRegex(ValueError,'Yukawa factors'):
            a4_lepton_dummy_alignment({'operators':[{'source_index':1,'flavour_matrix_axes_external_indices':['d1','d2']}]},a)


if __name__=='__main__':unittest.main()
