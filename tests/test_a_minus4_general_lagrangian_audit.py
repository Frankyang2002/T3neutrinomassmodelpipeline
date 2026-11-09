"""Synthetic source-completeness checks; no physical calculation is asserted."""
from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

from scripts.AuditAminus4GeneralLagrangian import (
    EXPECTED_QUARTIC_SYMBOLS, EXPECTED_EFT1_KEYS, audit_source,
)


class GeneralLagrangianAuditTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.source = self.root/'A_m4'
        (self.source/'data').mkdir(parents=True)
        self.beta = self.root/'eft1_rgbeta.json'
        (self.source/'native_status.json').write_text(json.dumps({
            'model':'A_-4','status':'FThresholdTreeWilsonExported',
            'operator_basis_is_full_general_requested':True,
            'full_general_extra_yukawa_count':0,
            'full_general_quartic_count':0,
            'full_general_quadratic_count':0,
        }))
        fields = {'IntegratedFields':['F'],'ActiveHeavyFields':['S1','S2']}
        for filename in ('eft1_after_F_wilson_seed.json','eft1_after_F_scalar_quartic_seed.json'):
            (self.source/'data'/filename).write_text(json.dumps(fields))
        (self.source/'native.log').write_text(
            'lambdaH2: raw SU(2) invariants=2, physical invariants=2\n'
            'Checking candidate: HiggsPortal2Inv1\n  result: InputForm[True]\n'
            'Checking candidate: HiggsPortal2Inv2\n  result: InputForm[True]\n'
        )
        (self.source/'data/eft1_after_F_tree_inputform.txt').write_text(
            ' + '.join(f'Coupling[{x}, {{}}, 0]' for x in sorted(EXPECTED_QUARTIC_SYMBOLS))
        )
        (self.beta).write_text(json.dumps({
            'status':'SuccessCandidateMassBeta',
            'metadata': {'GeneralScalarModel':'A_-4','GeneralScalarStage':'EFT1',
                         'SharedScalar':False, 'dS1':1,'dS2':3,'dF':2,
                         'alpha':-4,'GeneralScalarAdditionalQuartics':0,
                         'GeneralScalarNewBetaKeys':[], 'QuadraticKind':'none'},
            'report_betas':{key:'0' for key in EXPECTED_EFT1_KEYS}
        }))

    def test_complete_inventory_is_not_physical_certificate(self):
        j=audit_source(self.source,self.beta)
        self.assertTrue(j['full_general_renormalizable_eft1_source_ready'])
        self.assertFalse(j['physical_matchete_to_rgbeta_basis_map_certified'])
        self.assertFalse(j['full_dimension_five_operator_closure_verified'])

    def test_missing_higgs_portal_is_blocked(self):
        p=self.source/'data/eft1_after_F_tree_inputform.txt'
        p.write_text(p.read_text().replace('Coupling[lambdaH2Inv1, {}, 0]',''))
        j=audit_source(self.source,self.beta)
        self.assertFalse(j['full_general_renormalizable_eft1_source_ready'])
        self.assertEqual(j['missing_matched_tree_quartic_symbols'],['lambdaH2Inv1'])

    def test_missing_beta_coupling_is_blocked(self):
        j=json.loads(self.beta.read_text());del j['report_betas']['lambdaH2Adj']
        self.beta.write_text(json.dumps(j))
        self.assertFalse(audit_source(self.source,self.beta)['full_general_renormalizable_eft1_source_ready'])

    def test_lagrangian_request_not_equivalent_to_matched_completion(self):
        p=self.source/'native.log'
        p.write_text(p.read_text().replace('physical invariants=2','physical invariants=1'))
        self.assertFalse(audit_source(self.source,self.beta)['full_general_renormalizable_eft1_source_ready'])


if __name__ == '__main__':
    unittest.main()
