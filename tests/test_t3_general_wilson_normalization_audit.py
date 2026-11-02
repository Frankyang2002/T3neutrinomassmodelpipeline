"""General LLSS normalization gate: structural facts, never guessed spinor factor."""
from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from Numerical.running.GeneralT3WilsonNormalizationAudit import (
    STATUS, identify_chirality, verify_native_convention,
    verify_source_chiralities, verify_saved_polynomial, audit_model,
)
from Numerical.running.GeneralT3WilsonPolynomialBatch import build_polynomial_bridge
from Numerical.running.GeneralT3WilsonEvolution import MODELS, tag, digest
from scripts.RunT3GeneralWilsonNormalizationAudit import main
from tests.test_t3_general_wilson_polynomial_batch import fixture


def pl_raw():
    return ('Bar[Coupling[y1, {Index[d$$1, Flavor], Index[d$$9, NFlavor]}, 0]]'
            '*Bar[Coupling[y2, {Index[d$$2, Flavor], Index[d$$9, NFlavor]}, 0]]'
            '*NCM[Transp[Field[l, Fermion, {}, {}]],'
            'DiracProduct[GammaCC, Proj[-1]],Field[l, Fermion, {}, {}]]'
            '/Coupling[MF, {}, 0]')


def pr_raw():
    return ('Coupling[y1, {Index[d$$1, Flavor], Index[d$$9, NFlavor]}, 0]'
            '*Coupling[y2, {Index[d$$2, Flavor], Index[d$$9, NFlavor]}, 0]'
            '*NCM[Bar[Field[l, Fermion, {}, {}]],'
            'DiracProduct[GammaCC, Proj[1]],Transp[Bar[Field[l, Fermion, {}, {}]]]]'
            '/Coupling[MF, {}, 0]')


def source_fixture(model):
    ordered, source = fixture(model, reversed_flavour=model not in ('B_-1','C_-1'))
    for t in source['operator_terms']:
        t['source_term_inputform'] = pr_raw() if t.get('chirality')=='PR' else pl_raw()
    polynomial = build_polynomial_bridge(model, ordered, source)
    return ordered, source, polynomial


class NativeStructureTests(unittest.TestCase):
    def test_native_projector_and_bilinear(self):
        self.assertEqual(identify_chirality(pl_raw()), 'PL')
        self.assertEqual(identify_chirality(pr_raw()), 'PR')

    def test_swapped_projector_rejected(self):
        with self.assertRaisesRegex(ValueError,r'Proj\[1\] term'):
            identify_chirality(pl_raw().replace('Proj[-1]', 'Proj[1]'))

    def test_ambiguous_projector_rejected(self):
        with self.assertRaisesRegex(ValueError,'ambiguous'):
            identify_chirality(pl_raw()+'Proj[1]')

    def test_native_convention_requires_reported_weinberg_factor(self):
        sample=('StringContainsQ[termText, "Proj[-1]"], "PL", '
                'StringContainsQ[termText, "Proj[1]"], "PR", '
                '"L_EFT superset (1/2) C5^(ij) O_ij + h.c., '
                'C5^(ij) = A5^(ij) + A5^(ji)"')
        result=verify_native_convention(sample)
        self.assertFalse(result['GammaCC_to_two_component_spinor_phase_certified'])
        with self.assertRaises(ValueError):
            verify_native_convention(sample.replace('"PR"','"PL"'))


class ModelBatchTests(unittest.TestCase):
    def test_all_16_structural_fixtures(self):
        for model in MODELS:
            with self.subTest(model=model):
                ordered,source,polynomial=source_fixture(model)
                result=audit_model(model,ordered,source,polynomial)
                self.assertEqual(result['status'],STATUS)
                self.assertTrue(result['native_matchete_chirality']['selected_tensor_is_holomorphic_PL'])
                self.assertFalse(result['normalization_findings']['matchete_to_master_spinor_normalization_certified'])

    def test_previous_mixed_wrong_chirality_is_blocked(self):
        ordered,source,polynomial=source_fixture('B_-1')
        for term in source['operator_terms']:
            term['chirality']='PR' if term['chirality']=='PL' else 'PL'
        with self.assertRaisesRegex(ValueError,'chirality selection mismatch'):
            audit_model('B_-1',ordered,source,polynomial)

    def test_missing_matchete_raw_term_rejected(self):
        ordered,source,polynomial=source_fixture('A_-4')
        del source['operator_terms'][0]['source_term_inputform']
        with self.assertRaisesRegex(ValueError,'lacks exact InputForm'):
            audit_model('A_-4',ordered,source,polynomial)

    def test_changed_symmetric_polynomial_rejected(self):
        ordered,source,polynomial=source_fixture('A_-2')
        polynomial['independent_symmetric_coordinates'][0]['coefficient_real']+=1
        with self.assertRaisesRegex(ValueError,'normalized polynomial coefficient differs'):
            audit_model('A_-2',ordered,source,polynomial)

    def test_wrong_orbit_multiplicity_rejected(self):
        ordered,source,polynomial=source_fixture('D_+0')
        polynomial['independent_symmetric_coordinates'][0]['orbit_size']=100
        with self.assertRaisesRegex(ValueError,'orbit'):
            audit_model('D_+0',ordered,source,polynomial)

    def test_full_16_batch_runner_fixture(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            boundary=root/'ordered_v4'; polys=root/'polynomial_v2'; out=root/'audit'
            rows=[]
            for model in MODELS:
                ordered, source, poly=source_fixture(model)
                folder=boundary/tag(model)
                folder.mkdir(parents=True)
                source_path=(folder/'source_flavour.json' if model in ('B_-1','C_-1')
                             else folder/'original_source.json')
                source_path.write_text(json.dumps(source),encoding='utf-8')
                if model not in ('B_-1','C_-1'):
                    ordered['source_boundary']=str(source_path)
                    ordered['source_boundary_sha256']=digest(source_path)
                tensor_path=folder/'ordered_tensor.json'
                tensor_path.write_text(json.dumps(ordered),encoding='utf-8')
                poly['ordered_boundary_sha256']=digest(tensor_path)
                poly['operator_source_sha256']=digest(source_path)
                ppath=polys/tag(model)/'polynomial_bridge.json'
                ppath.parent.mkdir(parents=True)
                ppath.write_text(json.dumps(poly),encoding='utf-8')
                rows.append({'model':model,
                             'status':'OrderedNumericalRealComponentTensorFromMatchete',
                             'tree_operator_count_PL':1,
                             'real_component_count':1,
                             'tensor_sha256':digest(tensor_path)})
            (boundary/'summary.json').write_text(json.dumps({
                'format':'t3-ordered-source-flavour-and-real-scalar-tensor-v1',
                'ordered_real_tensor_count':16,'blocked_count':0,'results':rows
            }),encoding='utf-8')
            with patch('scripts.RunT3GeneralWilsonNormalizationAudit.master_evidence',
                       return_value={'test':True}):
                status=main(['--all','--boundary-root',str(boundary),
                             '--polynomial-root',str(polys),'--output-root',str(out)])
            self.assertEqual(status,0)
            summary=json.loads((out/'summary.json').read_text())
            self.assertEqual(summary['models_structure_verified'],16)
            self.assertEqual(summary['models_blocked'],0)
            self.assertFalse(summary['absolute_spinor_normalization_certified'])

    def test_mixed_conjugate_sector_does_not_double_pl(self):
        ordered,source,polynomial=source_fixture('C_-1')
        result=verify_source_chiralities('C_-1',ordered,source)
        self.assertEqual((result['native_PL_matched_term_count'],
                          result['native_PR_conjugate_term_count']), (1,1))
        self.assertTrue(result['Hermitian_conjugate_not_summed_into_PL_twice'])


if __name__=='__main__':
    unittest.main()
