"""Local source regressions for D(-2) CG export (no Wolfram required).

These checks exercise the serialized-name grammar and the exporter contract;
they are not a substitute for running Matchete against real tensor data.
"""
from __future__ import annotations

import re
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
WOLFRAM = ROOT / 'Lagrangian' / 'RunGeneralDMinus2Stage1.wl'
CG_PATTERN = re.compile(r'CG\[\s*(?:Bar\[\s*)?([A-Za-z$][A-Za-z0-9$]*)')


def cg_names(inputform: str) -> list[str]:
    return list(dict.fromkeys(CG_PATTERN.findall(inputform)))


class CGRegistryExportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = WOLFRAM.read_text(encoding='utf-8')

    def test_wilson_and_quartic_invariant_references(self):
        term = (
            'Bar[CG[T3MixCG, {Index[d$$1, SU2L[fund]]}]]*'
            'CG[T3Y1CG, {Index[d$$2, T3BSMd3]}]*'
            'CG[Bar[T3MixCG], {Index[d$$3, T3BSMd3]}]'
        )
        self.assertEqual(cg_names(term), ['T3MixCG', 'T3Y1CG'])

    def test_symbolic_cg_extraction_does_not_use_literal_head_cases(self):
        section = self.source.split('ClearAll[T3D2CGNames, T3D2CGRegistry];', 1)[1]
        section = section.split('realWilsonRecords = MapIndexed', 1)[0]
        self.assertIn('StringCases[', section)
        self.assertIn('ToString[term, InputForm]', section)
        self.assertIn('ToExpression[nameText, InputForm]', section)
        self.assertNotIn('HoldComplete @@ terms', section)

    def test_cg_names_are_checked_before_writing_seeds(self):
        src = self.source
        self.assertIn('eft1_after_F_realbasis_cg_registry_audit.json', src)
        self.assertIn('"MissingQuarticCGNames" -> realCGMissingQuartic', src)
        self.assertIn('"MissingWilsonCGNames" -> realCGMissingWilson', src)
        self.assertLess(
            src.index('realCGMissingQuartic = Complement['),
            src.index('eft1_after_F_realbasis_quartic_tensor_seed.json'),
        )
        self.assertIn('!StringStartsQ[tensorText, "SparseArray["]', src)


if __name__ == '__main__':
    unittest.main()
