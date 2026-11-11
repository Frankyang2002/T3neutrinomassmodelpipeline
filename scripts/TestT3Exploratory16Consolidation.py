"""Self-contained fixtures for the 16-model exploratory result consolidator."""
from __future__ import annotations
import json
from pathlib import Path
import tempfile
import unittest
from scripts.ConsolidateT3Exploratory16 import (
    EXPECTED, CASES, MIXED, RECOVERED,
    ConsolidationError, collect, write_outputs,
)


class TestConsolidation(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        root=Path(self.temp.name)
        self.folders={name:root/name for name in
                      ('mixed','nonmixed_initial','nonmixed_recovered')}
        for folder in self.folders.values():
            folder.mkdir()
        self.make_files()

    def tearDown(self):
        self.temp.cleanup()

    def make_files(self):
        for group,folder in self.folders.items():
            models = (MIXED if group=='mixed' else
                      RECOVERED if group=='nonmixed_recovered' else
                      set(EXPECTED)-MIXED-RECOVERED)
            entries=[]
            for model in sorted(models):
                entries.append({'model':model,'status':'Computed'})
                sub=folder/model
                sub.mkdir()
                cases={c:dict(case=c,uv_on=bool(int(c[0])),
                              eft1_on=bool(int(c[1])),
                              wilson_log_on=bool(int(c[2])),
                              smeft_on=bool(int(c[3])),
                              eft1_nfev=170,
                              simultaneous_vs_fixed_log_rel=0.02)
                       for c in CASES}
                (sub/'cases.json').write_text(json.dumps({
                    'status':'Computed','model':model,'cases':cases}))
            (folder/'summary.json').write_text(json.dumps({'model_rows':entries}))

    def test_full_16_models(self):
        data=collect(self.folders)
        self.assertEqual(data['report']['model_count'],16)
        self.assertEqual(data['report']['case_count'],256)
        self.assertEqual(len([r for r in data['rows'] if r['case']=='1111']),16)
        dest=Path(self.temp.name)/'result'
        write_outputs(data,dest)
        self.assertTrue((dest/'all16_cases.csv').exists())
        self.assertIn('B_-1',(dest/'all16_overview.md').read_text())
        with self.assertRaises(ConsolidationError):
            write_outputs(data,dest)

    def test_missing_case_rejected(self):
        path=self.folders['mixed']/'B_-1'/'cases.json'
        data=json.loads(path.read_text())
        del data['cases']['1111']
        path.write_text(json.dumps(data))
        with self.assertRaisesRegex(ConsolidationError,'expected 16'):
            collect(self.folders)

    def test_blocked_rejected(self):
        path=self.folders['nonmixed_initial']/'summary.json'
        data=json.loads(path.read_text())
        data['model_rows'][0]['status']='Blocked'
        path.write_text(json.dumps(data))
        with self.assertRaisesRegex(ConsolidationError,'not Computed'):
            collect(self.folders)

    def test_case_switch_rejected(self):
        path=self.folders['mixed']/'B_-1'/'cases.json'
        data=json.loads(path.read_text())
        data['cases']['1111']['eft1_on']=False
        path.write_text(json.dumps(data))
        with self.assertRaisesRegex(ConsolidationError,'wrong eft1_on'):
            collect(self.folders)


if __name__=='__main__':
    unittest.main()
