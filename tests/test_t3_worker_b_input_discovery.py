from __future__ import annotations
import json
import tempfile
import unittest
from pathlib import Path
from scripts.RunT3WorkerBEndpointAudit import _discover, _extract_endpoint, _explicit_sources

class SourceDiscoveryTests(unittest.TestCase):
    def test_find_valid_source(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            folder=root/'output'/'run1'/'A_-4';folder.mkdir(parents=True)
            p=folder/'running.json'
            p.write_text(json.dumps({'model':'A_-4','scales_gev':{'S':1000},
                         'eft1_at_scalar_threshold':{'mS1Sq':1.,'mS2Sq':2.}}))
            self.assertEqual(_discover(root,'A_-4'),[p])
            mu,state,digest=_extract_endpoint(p,'A_-4')
            self.assertEqual(mu,1000.)
            self.assertEqual(len(digest),64)
            self.assertEqual(_discover(root,'D_-2'),[])
    def test_reject_missing_bilinear(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'running.json'
            p.write_text(json.dumps({'model':'D_-2','scales_gev':{'S':1000},
                           'eft1_at_scalar_threshold':{'mS1Sq':1,'mS2Sq':2}}))
            with self.assertRaises(ValueError): _extract_endpoint(p,'D_-2')
    def test_explicit_source(self):
        result=_explicit_sources(['A_-4=output/run/running.json'],Path('/tmp/project'))
        self.assertEqual(result['A_-4'],Path('/tmp/project/output/run/running.json'))
    def test_ambiguity_not_silently_selected(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            for folder in ['one','two']:
                p=root/'output'/folder/'running.json';p.parent.mkdir(parents=True)
                p.write_text(json.dumps({'model':'A_-4','scales_gev':{'S':1000},
                    'eft1_at_scalar_threshold':{'mS1Sq':1,'mS2Sq':2}}))
            self.assertEqual(len(_discover(root,'A_-4')),2)
if __name__=='__main__': unittest.main()
