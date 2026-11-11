"""Offline catalogue/export validation tests for Batch 13 all-16 UV runner."""
from __future__ import annotations
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / 'scripts/RunBatch13UVExportsAll16.py'
spec = importlib.util.spec_from_file_location('batch13_uv', SCRIPT)
module = importlib.util.module_from_spec(spec)
assert spec and spec.loader
spec.loader.exec_module(module)


class Batch13UVTests(unittest.TestCase):
    def test_16_unique_models(self):
        cases = module.expected_models()
        self.assertEqual(len(cases), 16)
        self.assertEqual(len({c[0] for c in cases}), 16)
        self.assertEqual(sum(name in module.EXPECTED_QUADRATIC for name, _, _ in cases), 8)

    def test_negative_alpha_wolfram_tokens(self):
        self.assertEqual(module.integer_token(-4), 'm4')
        self.assertEqual(module.integer_token(0), '0')
        self.assertEqual(module.integer_token(2), '2')

    def test_registry_catalogue_agreement(self):
        self.assertEqual(len(module.load_registry()), 16)

    def test_validation_rejects_conjugate_yukawa(self):
        name, dims, alpha = module.expected_models()[0]
        entry = module.load_registry()[name]
        names = ['gY','g2','g3','yu','yd','ye','y1','y2','MF','mS1Sq','mS2Sq',
                 'lambdaH','lambdaT3', module.EXPECTED_QUADRATIC[name],
                 *[q['name'] for q in entry['new_rgbeta_symbols']], 'y1Conj']
        doc = {'status':'SuccessCandidateMassBeta',
               'metadata':{'GeneralScalarModel':name,'GeneralScalarStage':'UV',
                           'dS1':dims[0],'dS2':dims[1],'dF':dims[2],'alpha':alpha,
                           'QuadraticMassSymbol':module.EXPECTED_QUADRATIC[name],
                           'GeneralScalarAdditionalQuartics':entry['total_additional_quartics'],
                           'PhysicalMatcheteRGBetaMassBasisCertified':False},
               'beta_symbols':names,'report_betas':{k:'0' for k in names},
               'physical_mapping_certified':False}
        with tempfile.TemporaryDirectory() as folder:
            file = Path(folder) / 'fake.json'
            file.write_text(json.dumps(doc), encoding='utf-8')
            errors = module.validate_export(name,dims,alpha,entry,file)
        self.assertTrue(any('Restricted Yukawa violation' in e for e in errors))

    def test_dry_run_writes_16_record_summary(self):
        with tempfile.TemporaryDirectory() as folder:
            rc, data = module.execute(executable='wolframscript',
                output_dir=Path(folder)/'out',registry_path=module.REGISTRY,
                runner_path=module.RUNNER,dry_run=True)
            self.assertEqual(rc,0)
            self.assertEqual(len(data['models']),16)
            self.assertEqual(data['passed'],0)
            for name, dims, alpha in module.expected_models():
                command = data['models'][name]['command']
                self.assertEqual(command[3], 'UV')
                self.assertEqual(command[4:7], list(map(str, dims)))
                self.assertEqual(command[7], module.integer_token(alpha))
                self.assertEqual(command[8], name)
                self.assertEqual(len(command), 10)
            self.assertTrue((Path(folder)/'out/summary.json').is_file())


if __name__ == '__main__':
    unittest.main()
