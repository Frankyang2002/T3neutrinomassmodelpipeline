"""All-family polynomial bridge regressions; no live Wolfram needed."""
from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from Numerical.running.GeneralT3OrderedTensor import scalar_basis
from Numerical.running.GeneralT3WilsonEvolution import MODELS, tag, digest
from Numerical.running.GeneralT3WilsonPolynomialBatch import (
    canonical_key, distinct_orbit, canonical_polynomial, polynomial_identity,
    source_term_alignment, ordered_aligned_components, build_polynomial_bridge,
    master_weinberg_polynomial_check,
)
from scripts.RunT3GeneralWilsonPolynomialBatch import main


def mat(where: dict[tuple[int,int], complex]):
    z = np.zeros((3,3),complex)
    for (p,q), val in where.items():
        z[p-1,q-1] = val
    return {'real':z.real.tolist(), 'imag':z.imag.tolist()}


def fixture(model: str, *, reversed_flavour: bool=True,
            source_indices: tuple[int,...]=(2,), values: tuple[complex,...]=(1+2j,)):
    b = scalar_basis(model)
    names = [n for n in b.dimensions if n != 'H']
    assert len(names)>=2
    a, s = b.offsets[names[0]], b.offsets[names[1]]
    ops, terms = [], []
    for idx, val in zip(source_indices, values):
        field_order=(('d$$2','d$$2'),('d$$1','d$$1')) if reversed_flavour else (
            ('d$$1','d$$1'),('d$$2','d$$2'))
        tpl='NCM['+' , '.join(
            f'Field[l, Fermion, {{Index[{u}, SU2L[fund]], Index[{f}, Flavor]}}, {{}}]'
            for u,f in field_order)+',GammaCC]'
        terms.append({'source_index':idx,
                      'chirality':'PL',
                      'flavour_matrix_axes_external_indices':['d$$1','d$$2'],
                      'factors':[{'external_index':'d$$1'},{'external_index':'d$$2'}],
                      'operator_template_unit_flavor_inputform':tpl})
        ops.append({'source_index':idx,
                    'basis_scalar_fields':dict(b.dimensions),
                    'real_scalar_dimension_including_higgs':b.total_dimension,
                    'flavour_matrix_axes_external_indices':['d$$1','d$$2'],
                    'ordered_not_canonical_symmetrized':True,
                    'component_count':1,
                    'components':[{'lepton_su2':[1,2], 'scalar_real':[a,s],
                                   'ordered_flavour_coefficient_gev_inv':mat({(1,2):val})}]})
    payload={
        'format':'t3-mixed-native-matchete-ordered-tensor-v1' if model in ('B_-1','C_-1')
                   else 't3-ordered-source-flavour-and-real-scalar-tensor-v1',
        'status':'OrderedNumericalRealComponentTensorFromMatchete',
        'model':model,'scalar_basis_explicit':True,
        'matchete_cg_contractions_retained':True,
        'external_three_family_flavour_order_not_symmetrized':True,
        'canonical_C_ijab_normalization_certified':False,
        'mu_F_gev':1e5,'tree_operator_count_PL':len(ops),
        'total_ordered_real_component_count':len(ops),
        'operators':ops,
    }
    if model in ('B_-1','C_-1'):
        source={'format':'t3-mixed-factorised-LLSS-source-v1',
                'mu_F_gev':1e5,
                'source_status_sha256':'a'*64,
                'source_mass_basis_sha256':'b'*64,
                'llss_sha256':'c'*64,
                'status':'NumericalOrderedSourceFlavourCGExported',
                'model':model,'operator_terms':terms+[
                    {**t, 'source_index':t['source_index']+100, 'chirality':'PR'}
                    for t in terms]}
        payload.update({'source_status_sha256':'a'*64,
                        'source_mass_basis_sha256':'b'*64,
                        'source_native_fresh_llss_sha256':'c'*64})
    else:
        source={'format':'t3-source-exact-F-tree-flavour-operator-handoff-v1',
                'mu_F_gev':1e5,
                'status':'SourceFactoredOrderedFlavourMatricesAtF_NotCanonicalWilsonTensor',
                'model':model,'operator_terms':terms}
    return payload, source


class PolynomialTests(unittest.TestCase):
    def test_fully_general_complex_monomials(self):
        source={(1,1,2,2,5,7):1+2j, (2,2,1,1,7,5):3-4j,
                (1,3,1,3,6,6):2j, (2,1,2,1,8,9):-5}
        independent, _=canonical_polynomial(source)
        report=polynomial_identity(source,independent)
        self.assertTrue(report['preserved_to_numerical_precision'])
        self.assertEqual(len(distinct_orbit((1,3,1,3,6,6))),1)
        self.assertEqual(len(distinct_orbit((1,1,2,2,5,7))),4)
        self.assertEqual(independent[(1,1,2,2,5,7)], (4-2j)/4)

    def test_all_16_scalar_bases(self):
        for model in MODELS:
            with self.subTest(model=model):
                payload, source=fixture(model)
                result=build_polynomial_bridge(model,payload,source)
                self.assertEqual(result['status'],'PolynomialRepresentativeVerified_NotSpinorNormalisation')
                self.assertTrue(result['polynomial_identity']['preserved_to_numerical_precision'])
                self.assertEqual(result['flavour_permutation_histogram']['[1,0]'],1)
                self.assertFalse(result['overall_matchete_spinor_normalization_certified'])

    def test_mixed_field_order_and_flavour_no_swap(self):
        for model in ('B_-1','C_-1'):
            payload,source=fixture(model,reversed_flavour=False)
            result=build_polynomial_bridge(model,payload,source)
            self.assertEqual(result['flavour_permutation_histogram']['[0,1]'],1)
            row=result['independent_symmetric_coordinates'][0]
            self.assertEqual(row['coordinate'][1],1)
            self.assertEqual(row['coordinate'][3],2)

    def test_flavour_swap_offdiagonal(self):
        payload,source=fixture('A_-4')
        aligned, perms,_=ordered_aligned_components('A_-4',payload,source)
        a,b=[scalar_basis('A_-4').offsets[n] for n in ('NewScalar1','NewScalar2')]
        self.assertEqual(aligned[(1,2,2,1,a,b)],1+2j)
        self.assertNotIn((1,1,2,2,a,b),aligned)
        self.assertEqual(perms[0]['flavour_axis_permutation_for_lepton_order'],[1,0])

    def test_different_source_terms_permuted_before_sum(self):
        payload,source=fixture('A_-4',source_indices=(2,3),values=(1+0j,2+0j))
        # For term 3, first lepton field follows axis 0 rather than axis 1.
        source['operator_terms'][1]['operator_template_unit_flavor_inputform'] = (
            'NCM[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}],'
            'Field[l, Fermion, {Index[d$$2, SU2L[fund]], Index[d$$2, Flavor]}, {}],GammaCC]'
        )
        aligned,perms,_=ordered_aligned_components('A_-4',payload,source)
        a,b=[scalar_basis('A_-4').offsets[n] for n in ('NewScalar1','NewScalar2')]
        self.assertEqual(aligned[(1,2,2,1,a,b)],1)
        self.assertEqual(aligned[(1,1,2,2,a,b)],2)
        self.assertEqual([x['flavour_axis_permutation_for_lepton_order'] for x in perms],
                         [[1,0],[0,1]])

    def test_wrong_flavour_axis_fail_closed(self):
        payload,source=fixture('A_-4')
        source['operator_terms'][0]['factors'][0]['external_index']='wrong'
        with self.assertRaisesRegex(ValueError,'Yukawa'):
            build_polynomial_bridge('A_-4',payload,source)

    def test_undistinguished_field_fail_closed(self):
        payload,source=fixture('A_-4')
        source['operator_terms'][0]['operator_template_unit_flavor_inputform'] = (
            'Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}]*GammaCC'
        )
        with self.assertRaisesRegex(ValueError,'exactly two'):
            build_polynomial_bridge('A_-4',payload,source)

    def test_mixed_stale_matchete_provenance_rejected(self):
        payload, source = fixture('B_-1', reversed_flavour=False)
        source['llss_sha256'] = '0' * 64
        with self.assertRaisesRegex(ValueError,'provenance'):
            build_polynomial_bridge('B_-1',payload,source)

    def test_original_scale_mismatch_rejected(self):
        payload, source = fixture('D_+0')
        source['mu_F_gev'] = 1e4
        with self.assertRaisesRegex(ValueError,'threshold scales'):
            build_polynomial_bridge('D_+0',payload,source)

    def test_mixed_missing_pr_fails_closed(self):
        payload,source=fixture('C_-1',reversed_flavour=False)
        source['operator_terms']=source['operator_terms'][:1]
        with self.assertRaisesRegex(ValueError,'PL/PR'):
            build_polynomial_bridge('C_-1',payload,source)

    def test_repeated_index_same_scalar_orbit_size(self):
        self.assertEqual(len(distinct_orbit((1,1,1,1,5,5))),1)
        self.assertEqual(len(distinct_orbit((1,1,1,1,5,6))),2)
        self.assertEqual(len(distinct_orbit((1,1,1,2,5,5))),2)

    def test_master_reconstruction_if_available(self):
        try:
            result=master_weinberg_polynomial_check()
        except ImportError:
            self.skipTest('Full RGE master package not mounted in isolated test environment')
        self.assertEqual(result['master_weinberg_polynomial_residual'],'0')


class RunnerTests(unittest.TestCase):
    def test_full_16_model_fixture_batch(self):
        with tempfile.TemporaryDirectory() as t:
            base=Path(t)
            root=base/'v3'; out=base/'batch'; root.mkdir()
            rows=[]
            for model in MODELS:
                payload, source=fixture(model, reversed_flavour=(model not in ('B_-1','C_-1')))
                folder=root/tag(model); folder.mkdir()
                sfile=folder/('source_flavour.json' if model in ('B_-1','C_-1') else 'original.json')
                sfile.write_text(json.dumps(source),encoding='utf-8')
                if model not in ('B_-1','C_-1'):
                    payload['source_boundary']=str(sfile)
                    payload['source_boundary_sha256']=digest(sfile)
                pfile=folder/'ordered_tensor.json'
                pfile.write_text(json.dumps(payload),encoding='utf-8')
                rows.append({'model':model,'status':payload['status'],
                             'tree_operator_count_PL':1,'real_component_count':1,
                             'tensor_sha256':digest(pfile)})
            (root/'summary.json').write_text(json.dumps({
                'format':'t3-ordered-source-flavour-and-real-scalar-tensor-v1',
                'ordered_real_tensor_count':16,'blocked_count':0,'results':rows}),encoding='utf-8')
            self.assertEqual(main(['--all','--skip-master-check',
                                   '--boundary-root',str(root),'--output-root',str(out)]),0)
            report=json.loads((out/'summary.json').read_text())
            self.assertEqual(report['models_verified'],16)
            self.assertEqual(report['models_blocked'],0)
            self.assertFalse(report['matchete_to_master_spinor_normalization_certified'])
            self.assertFalse(report['physical_eft1_wilson_rge_computed'])
            a=json.loads((out/'A_m4'/'polynomial_bridge.json').read_text())
            self.assertEqual(a['lepton_flavour_alignments'][0]
                             ['flavour_axis_permutation_for_lepton_order'],[1,0])
            b=json.loads((out/'B_m1'/'polynomial_bridge.json').read_text())
            self.assertEqual(b['lepton_flavour_alignments'][0]
                             ['flavour_axis_permutation_for_lepton_order'],[0,1])


if __name__=='__main__':
    unittest.main()
