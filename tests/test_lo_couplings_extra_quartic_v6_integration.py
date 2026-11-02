"""Real checkout tests: generator declarations and V5 success evidence."""
from __future__ import annotations
from pathlib import Path
import json
import unittest
from Numerical.lo_certification.couplings.extra_quartic_source_v6 import (
    ALL_MODELS,generator_contract,model_registry_quartics)
from Numerical.lo_certification.couplings.extra_quartic_batch_v6 import (_validate_registry,_load_v5)

ROOT=Path(__file__).resolve().parents[1]
AVAILABLE=(ROOT/'Lagrangian/interactions/GeneralScalarOperators.wl').exists()

@unittest.skipUnless(AVAILABLE,'Requires full local repository')
class TestRealV6Generator(unittest.TestCase):
    def test_live_matchete_operator_source_contract(self):
        r=generator_contract(ROOT)
        self.assertTrue(r['Matchete_plus_HC_source_pattern_verified'])
    def test_full50_generator_name_inventory(self):
        registry=json.loads((ROOT/'RGE/general/generated/general_rgbeta_quartic_registry.json').read_text(encoding='utf8'))
        _validate_registry(registry)
        total=hhss=0
        for model in ALL_MODELS:
            rows=model_registry_quartics(model,registry)
            total+=len(rows);hhss+=sum(x['sector_HHSS'] for x in rows)
        self.assertEqual(total,50);self.assertEqual(hhss,24)
    def test_prior_V5_success_report(self):
        s=_load_v5(ROOT,'ws2_matchete_rgbeta_dictionary_v1')
        self.assertEqual(s['extra_HHSS_source_symbols_crosswalked'],24)

if __name__=='__main__':unittest.main()
