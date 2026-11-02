"""Standalone exact CG/source tensor tests. No local T3 repository needed."""
from __future__ import annotations

from pathlib import Path
import sys
import unittest
import sympy as sp

from Numerical.lo_certification.couplings.source_cg_v5 import (
    MatcheteCGError,ExactCG,_exact_number,parse_exact_cg,parse_all_cgs,
    inspect_quartic_term,inspect_all_terms,candidate_raw_tensor)
from Numerical.lo_certification.couplings.matchete_rgbeta_dictionary_v5 import (
    exact_basis_map,candidate_cg_reference_maps,native_seed_sector_matches,
    summarize_all_extra_registry,per_named_HHSS_evidence)
from Numerical.lo_certification.couplings.named_dictionary_batch_v5 import run_v5_batch
from Numerical.lo_certification.couplings.extra_hhss_geometry import reference_invariants

EXAMPLE={
 'Name':'T3HiggsPortal1CGInv2',
 'RepsInputForm':'{SU2L[fund], Bar[SU2L[fund]], T3BSMd2, Bar[T3BSMd2]}',
 'TensorInputForm':'SparseArray[Automatic, {2, 2, 2, 2}, 0, {1, {{0, 3, 6}, {{1, 1, 1}, {1, 2, 2}, {2, 2, 1}, {1, 1, 2}, {2, 1, 1}, {2, 2, 2}}}, {1/Sqrt[3], 2/Sqrt[3], -(1/Sqrt[3]), -(1/Sqrt[3]), 2/Sqrt[3], 1/Sqrt[3]}}]'
}


class TestExactNativeMatcheteCG(unittest.TestCase):
    def test_parse_real_source_sparse_portal(self):
        cg=parse_exact_cg(EXAMPLE)
        self.assertEqual(cg.dims,(2,2,2,2))
        self.assertEqual(len(cg.entries),6)
        self.assertEqual(cg.component((0,0,0,0)),1/sp.sqrt(3))
        self.assertEqual(cg.component((0,0,1,1)),2/sp.sqrt(3))
        self.assertEqual(cg.component((0,1,1,1)),0)
        self.assertEqual(cg.invariant_summary()['norm_squared_exact'],'4')

    def test_real_source_complex_phase_values(self):
        for expr,want in (('I/Sqrt[2]',sp.I/sp.sqrt(2)),
                          ('(-I)/Sqrt[2]',-sp.I/sp.sqrt(2)),
                          ('3^(-1/4)',sp.Pow(3,sp.Rational(-1,4))),
                          ('Sqrt[3/11]',sp.sqrt(sp.Rational(3,11))),
                          ('1/2',sp.Rational(1,2))):
            with self.subTest(expr=expr):
                self.assertEqual(sp.simplify(_exact_number(expr)-want),0)

    def test_bars_are_not_dropped_in_rep_metadata(self):
        cg=parse_exact_cg(EXAMPLE)
        self.assertIn('Bar[SU2L[fund]]',cg.reps)

    def test_invalid_or_unsupported_expression_fails_closed(self):
        for value in ('a','Gamma[1]','foo()','2;Delete[]','1/0'):
            with self.subTest(value=value):
                with self.assertRaises(MatcheteCGError):_exact_number(value)

    def test_CSR_offsets_fail_closed(self):
        data=dict(EXAMPLE)
        data['TensorInputForm']=data['TensorInputForm'].replace('{0, 3, 6}','{0, 7, 6}')
        with self.assertRaises(MatcheteCGError):parse_exact_cg(data)

    def test_CSR_coordinate_bounds_fail_closed(self):
        data=dict(EXAMPLE)
        data['TensorInputForm']=data['TensorInputForm'].replace('{{1, 1, 1},','{{3, 1, 1},')
        with self.assertRaises(MatcheteCGError):parse_exact_cg(data)

    def test_representation_dimension_mismatch_is_rejected(self):
        data=dict(EXAMPLE)
        data['RepsInputForm']=data['RepsInputForm'].replace('T3BSMd2','T3BSMd3')
        with self.assertRaises(MatcheteCGError):parse_exact_cg(data)

    def test_exact_sparse_components_in_report(self):
        cg=parse_exact_cg(EXAMPLE)
        report=cg.invariant_summary()
        self.assertEqual(len(report['tensor_exact_nonzero_components']),6)
        self.assertEqual(report['tensor_exact_nonzero_components'][0]['coordinates_zero_based'],[0,0,0,0])
        self.assertTrue(report['tensor_exact_component_digest_sha256'])

    def test_duplicate_CG_registry_symbol_rejected(self):
        with self.assertRaises(MatcheteCGError):parse_all_cgs({'CGRegistry':[EXAMPLE,EXAMPLE]})

    def test_native_term_counts_include_bars(self):
        term={'Index':1,'HiggsFieldCount':2,'Scalar1FieldCount':2,'Scalar2FieldCount':0,
             'CGNames':[],'TermInputForm':(
                 'Bar[Field[H, Scalar, {}, {}]]*Field[H, Scalar, {}, {}]*'
                 'Bar[Field[NewScalar1, Scalar, {}, {}]]*Field[NewScalar1, Scalar, {}, {}]*'
                 'Coupling[lambdaH1Inv1, {}, 0]')}
        doc=inspect_quartic_term(term,{})
        self.assertEqual(doc['counts'],[1,1,1,1,0,0])
        self.assertEqual(doc['native_coupling_symbols'],['lambdaH1Inv1'])

    def test_quartic_term_declared_count_mismatch_rejected(self):
        term={'Index':1,'HiggsFieldCount':4,'Scalar1FieldCount':0,'Scalar2FieldCount':0,
             'CGNames':[],'TermInputForm':'Field[H, Scalar, {}, {}]*Field[H, Scalar, {}, {}]*Field[NewScalar1, Scalar, {}, {}]*Field[NewScalar1, Scalar, {}, {}]'}
        with self.assertRaises(MatcheteCGError):inspect_quartic_term(term,{})

    def test_unknown_CG_is_rejected(self):
        term={'Index':1,'HiggsFieldCount':2,'Scalar1FieldCount':2,'Scalar2FieldCount':0,
             'CGNames':['T3HiggsPortal1CGInv2'],'TermInputForm':(
                 'Field[H, Scalar, {}, {}]*Field[H, Scalar, {}, {}]*'
                 'Field[NewScalar1, Scalar, {}, {}]*Field[NewScalar1, Scalar, {}, {}]*'
                 'CG[T3HiggsPortal1CGInv2, {}]')}
        with self.assertRaises(MatcheteCGError):inspect_quartic_term(term,{})

    def test_matchete_CG_conservation_metadata(self):
        rows=candidate_raw_tensor(parse_exact_cg(EXAMPLE),(0,1,2,3),(2,2,2,2,2,2))
        self.assertEqual(len(rows),24)
        self.assertEqual(rows[0]['raw_tensor'].rows,16)


class TestExactDictionaryMatrix(unittest.TestCase):
    def test_rational_matrix_map_exact(self):
        e1=sp.Matrix([1,0,1]);e2=sp.Matrix([0,1,2])
        got=exact_basis_map([e1,e2],[2*e1+3*e2,sp.Rational(3,5)*e1-e2])
        self.assertTrue(got['full_residual_is_zero'])
        self.assertEqual(got['coefficients'],[['2','3/5'],['3','-1']])
        self.assertTrue(got['reversible'])

    def test_phase_preserved_not_absorbed(self):
        v=sp.Matrix([1,sp.I])
        got=exact_basis_map([v],[sp.I*v])
        self.assertEqual(got['coefficients'],[['I']])

    def test_outside_span_is_not_approximated(self):
        a=exact_basis_map([sp.Matrix([1,0])],[sp.Matrix([0,1])])
        self.assertEqual(a['status'],'target_outside_source_span')
        self.assertIsNone(a['coefficients'])

    def test_dependent_source_rejected(self):
        a=exact_basis_map([sp.Matrix([1,0]),sp.Matrix([2,0])],[sp.Matrix([1,0])])
        self.assertEqual(a['status'],'dependent_source_tensors')

    def test_invalid_vectors_rejected(self):
        with self.assertRaises(MatcheteCGError):
            exact_basis_map([sp.Matrix([1,0])],[sp.Matrix([1,0,0])])

    def test_real_A_triplet_CG_candidate_orientation_is_explicit(self):
        from itertools import product
        counts=(0,2,0,1,0,1)
        types,dims,ref=reference_invariants(1,3,counts)
        self.assertEqual(dims,(2,2,1,3))
        entries=[]
        for i,j,k in product(range(2),range(2),range(3)):
            value=sp.simplify(ref[0][i*6+j*3+k])
            if value!=0:entries.append(((i,j,k),value))
        cg=ExactCG('SyntheticTripletCG',('H','BarH','S2'),(2,2,3),
                   tuple(entries),'unit-test')
        result=candidate_cg_reference_maps('A_-2',counts,cg)
        self.assertGreaterEqual(len(result),2)
        self.assertTrue(all(x['component_residual_zero'] for x in result))
        self.assertTrue(any(x['exact_reference_coefficients']==['1'] for x in result))
        # The existence of multiple Ward-invariant embeddings is *not* a
        # physical Matchete-to-RGBeta phase certification.

    def test_no_invented_native_extra_match(self):
        term={'term_index':3,'counts':[0,2,0,1,0,1],
              'native_coupling_symbols':['lambdaT3'], 'CGNames':['T3MixCG'],
              'term_source_sha256':'hello'}
        hits=native_seed_sector_matches([term],(0,2,0,1,0,1),'lambdaT3ExtraQ020101I1')
        self.assertEqual(len(hits),1)
        self.assertFalse(hits[0]['is_exact_named_coupling_match'])

    def test_extra_quartic_50_inventory_includes_SSSS(self):
        rows=summarize_all_extra_registry('A_-2',{'new_rgbeta_symbols':[
             {'name':'lambdaT3GeneralQ020101I1','counts':[0,2,0,1,0,1]},
             {'name':'lambdaT3GeneralQ000004I1','counts':[0,0,0,0,0,4]}],
             'preexisting_rgbeta_symbols':[]},[])
        self.assertEqual(rows['additional_quartic_registry_entries'],2)
        self.assertEqual(rows['HHSS_entries'],1)

    def test_run_owner_checked_before_sources(self):
        import tempfile
        with tempfile.TemporaryDirectory() as t:
            root=Path(t)
            with self.assertRaises(ValueError):run_v5_batch(root,root/'RGE/foreign')
            self.assertFalse((root/'RGE').exists())

    def test_failed_naming_is_explicit(self):
        prev={'entries':[{'rgbeta_symbol':'lambdaT3GeneralQ020101I1',
                          'counts':[0,2,0,1,0,1],
                          'valid_canonical_invariant_orientation_count':1,
                          'candidate_orientation_projection_maps':[{}]}]}
        entry={'category':'new','counts':(0,2,0,1,0,1),
               'declaration':{'symbol':'lambdaT3GeneralQ020101I1'}}
        r=per_named_HHSS_evidence('A_-2',entry,[],{},prev)
        self.assertEqual(r['matching_named_native_terms'],0)
        self.assertIn('matching_named_Matchete_extra_quartic_absent_in_native_EFT1_seed',r['blocking_dependencies'])
        self.assertIsNone(r['certified_Matchete_RGBeta_coupling_matrix'])


if __name__=='__main__':unittest.main()
