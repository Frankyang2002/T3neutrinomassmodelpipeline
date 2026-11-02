"""Run V5 against real native source cache only when checkout is mounted."""
from __future__ import annotations

from pathlib import Path
import json
import unittest

from Numerical.lo_certification.couplings.source_cg_v5 import parse_all_cgs,inspect_all_terms
from Numerical.lo_certification.couplings.named_dictionary_batch_v5 import _source_50_declarations

ROOT=Path(__file__).resolve().parents[1]
HAS_REPO=(ROOT/'RGE/general/generated/general_rgbeta_quartic_registry.json').exists()


@unittest.skipUnless(HAS_REPO,'Local full T3 repository unavailable')
class TestSourceEFT1MatcheteCGs(unittest.TestCase):
    def test_all16_seed_sparse_CGs_and_field_counts(self):
        reg=json.loads((ROOT/'RGE/general/generated/general_rgbeta_quartic_registry.json').read_text(encoding='utf8'))
        for model in reg['models']:
            family,alpha=model.split('_')
            tag=('m'+alpha[1:] if alpha.startswith('-') else 'p'+alpha[1:])
            seed=ROOT/f'output/full16_symbolic/T3_{family}_alpha_{tag}/data/eft1_after_F_scalar_quartic_seed.json'
            with self.subTest(model=model):
                self.assertTrue(seed.is_file())
                data=json.loads(seed.read_text(encoding='utf8'))
                cgs=parse_all_cgs(data)
                terms=inspect_all_terms(data,cgs)
                self.assertGreater(len(terms),0)
                self.assertGreater(len(cgs),0)

    def test_all50_Wolfram_declarations_match_registry(self):
        registry=json.loads((ROOT/'RGE/general/generated/general_rgbeta_quartic_registry.json').read_text(encoding='utf8'))
        total=0;hhss=0
        for model,row in registry['models'].items():
            with self.subTest(model=model):
                data=_source_50_declarations(ROOT,model,row)
                total+=data['total_source_declarations_verified']
                hhss+=sum(sum(x['counts'][:2])==2 for x in data['entries'])
        self.assertEqual(total,registry['totals']['matchete_extra_quartics'])
        self.assertEqual(hhss,24)

    def test_V4_source_report_readonly(self):
        p=ROOT/'output/general_observable/lo_certification_couplings/ws2_named_hhss_v1/summary.json'
        if not p.exists():self.skipTest('V4 local output unavailable')
        prior=json.loads(p.read_text(encoding='utf8'))
        self.assertTrue(prior['all16_source_projections_passed'])
        self.assertEqual(prior['extra_HHSS_source_symbols_examined'],24)

if __name__=='__main__':unittest.main()
