"""Local checkout regressions: require a source-populated T3 repository."""
from __future__ import annotations
import json
from pathlib import Path
import unittest

from Numerical.lo_certification.couplings.named_hhss_projection import named_hhss_sources

ROOT=Path(__file__).resolve().parents[1]
HAS_REPO=(ROOT/'RGE/general/generated/general_rgbeta_quartic_registry.json').exists()


@unittest.skipUnless(HAS_REPO,'Full local source repository unavailable')
class TestLocalNamedHHSSSourceConcordance(unittest.TestCase):
    def test_all16_named_declarations_exact_registry_concordance(self):
        reg=json.loads((ROOT/'RGE/general/generated/general_rgbeta_quartic_registry.json').read_text(encoding='utf8'))
        self.assertEqual(reg['status'],'GeneratedCandidateBasisNotWolframVerified')
        total=0
        for model,values in reg['models'].items():
            with self.subTest(model=model):
                doc=named_hhss_sources(ROOT,model,values)
                n=len(doc['sectors'])
                self.assertEqual(n, sum(sum(x['counts'][:2])==2 and sum(x['counts'][2:])==2
                                        for x in values.get('new_rgbeta_symbols',[])+values.get('preexisting_rgbeta_symbols',[])))
                total+=n
        self.assertEqual(total,24)

    def test_v3_successful_run_is_kept_untouched(self):
        p=ROOT/'output/general_observable/lo_certification_couplings/ws2_quartic_completion_v1/summary.json'
        if not p.exists():self.skipTest('The user has not generated the V3 output yet')
        doc=json.loads(p.read_text(encoding='utf8'))
        self.assertEqual(doc['models_processed'],16)
        self.assertTrue(doc['all16_model_audits_passed'])


if __name__=='__main__':unittest.main()
