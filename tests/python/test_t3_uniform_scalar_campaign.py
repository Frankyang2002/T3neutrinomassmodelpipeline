"""Isolated campaign/regression tests. No physical running or matching."""
from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

from scripts import RunT3UniformScalarCampaign as campaign
from scripts import RunT3UniformScalarBenchmark as uniform


class CampaignLogicTests(unittest.TestCase):
    @staticmethod
    def original():
        ys = np.full((3,3), 1e-4).tolist()
        zeros = np.zeros((3,3)).tolist()
        return {'representation':{'model':'dummy'},
                'scales':{'mu_fermion_threshold_gev':1e5,
                          'mu_scalar_threshold_gev':1e3,'mu_low_gev':1e2},
                'base_state':{'mu_gev':1e7,
                    'sm':{'lambdaH':0.25,'g2':0.59493},
                    'ordinary':{'lambdaT3':{'real':0.1,'imag':0},
                                'lambdaS1':0.1,'lambdaS2':0.12,
                                'lambdaH1':0.02,'lambdaH2':0.03,'lambda12':0.01,
                                'mS1Sq':1306122.448979592,'mS2Sq':1e6,
                                'y1':{'real':ys,'imag':zeros},
                                'y2':{'real':ys,'imag':zeros},
                                'y1Conj':{'real':ys,'imag':zeros}}}}

    def test_reference_must_change_quartics_and_zero_conjugates(self):
        old=self.original()
        new,_=uniform.uniform_config(old,0.1)
        campaign.check_uniform_reference(new,old)
        self.assertEqual(new['base_state']['ordinary']['y1Conj']['real'][0][0],0)
        self.assertEqual(new['base_state']['ordinary']['lambdaS2'],0.01)
        self.assertEqual(new['base_state']['sm']['lambdaH'],0.25)
        self.assertEqual(old['base_state']['ordinary']['lambdaS2'],0.12)

    def test_reference_rejects_old_quartic_or_mass_or_nonzero_conjugate(self):
        old=self.original()
        for name,value in [('lambdaS2',0.12),('mS1Sq',999),
                           ('y1Conj',old['base_state']['ordinary']['y1Conj'])]:
            changed,_=uniform.uniform_config(old,0.1)
            changed['base_state']['ordinary'][name]=deepcopy(value)
            with self.subTest(name=name),self.assertRaisesRegex(ValueError,'boundary'):
                campaign.check_uniform_reference(changed,old)

    def test_reference_rejects_other_lambda_and_sm_change(self):
        old=self.original()
        new,_=uniform.uniform_config(old,0.005)
        with self.assertRaisesRegex(ValueError,'lambdaT3'):
            campaign.check_uniform_reference(new,old)
        new,_=uniform.uniform_config(old,0.1)
        new['base_state']['sm']['lambdaH']=0.01
        with self.assertRaisesRegex(ValueError,'SM inputs'):
            campaign.check_uniform_reference(new,old)

    def test_select_model_inventory(self):
        self.assertEqual(campaign.select_models('pilot'),('A_-4','B_-1'))
        self.assertEqual(len(campaign.select_models('all')),16)
        with self.assertRaisesRegex(ValueError,'Invalid'):
            campaign.select_models('A_-4,A_-4')
        with self.assertRaisesRegex(ValueError,'Invalid'):
            campaign.select_models('not_a_model')

    def test_persisted_reference_checks_hash_and_fit_acceptance(self):
        old=self.original()
        new,_=uniform.uniform_config(old,0.1)
        new['base_state']['ordinary']['y1']['real'][0][0]=0.003
        new['thesis_uniform_scalar_benchmark']={'schema':uniform.SCHEMA}
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            path=root/'source.json';path.write_text('{}',encoding='utf-8')
            manifest=root/'merged.json';manifest.write_text('{}',encoding='utf-8')
            paths={'config':path}
            with patch.object(campaign,'ROOT',root):
                prov=campaign.provenance(paths,manifest)
                row={'schema':campaign.SCHEMA,'stage':'reference','model':'A_-4',
                     'status':'AcceptedApproxFit','chi2':0.01,'lambdaT3':0.1,
                     'provenance':prov,'fitted_yukawas':{
                         k:new['base_state']['ordinary'][k]
                         for k in ('y1','y2')}}
                folder=root/'A_m4'
                campaign.base.json_save(folder/'result.json',row)
                campaign.base.json_save(folder/'fitted_config.json',new)
                readrow,readcfg=campaign.validated_reference(folder,model='A_-4',
                    source_paths=paths,manifest_path=manifest,
                    original_config=old,chi2_cutoff=5.86)
                self.assertEqual(readrow['chi2'],0.01)
                self.assertAlmostEqual(readcfg['base_state']['ordinary']['y1']['real'][0][0],0.003)
                row['chi2']=8
                campaign.base.json_save(folder/'result.json',row)
                with self.assertRaisesRegex(ValueError,'chi2 cutoff'):
                    campaign.validated_reference(folder,model='A_-4',
                        source_paths=paths,manifest_path=manifest,
                        original_config=old,chi2_cutoff=5.86)
                row['chi2']=0.01
                campaign.base.json_save(folder/'result.json',row)
                manifest.write_text('{"changed":true}',encoding='utf-8')
                with self.assertRaisesRegex(ValueError,'fingerprints'):
                    campaign.validated_reference(folder,model='A_-4',
                        source_paths=paths,manifest_path=manifest,
                        original_config=old,chi2_cutoff=5.86)

    def test_reference_stage_is_at_original_lambda_and_sends_original_seed(self):
        old=self.original()
        orig=campaign.base.fitted_vector(old)
        fit={'chi2_final':0.1}
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            source=root/'source.json';source.write_text('{}',encoding='utf-8')
            manifest=root/'merged.json';manifest.write_text('{}',encoding='utf-8')
            zee=SimpleNamespace(MERGED_FILE=manifest,sources=lambda *args:(old,fit,{'config':source}))
            args=SimpleNamespace(canonical_run_id='full16_general_v1',
                                 max_nfev=5,accept_seed=False)
            recorded=[]
            class FakeFactory:
                def __init__(self,model,*args):
                    self.config=None
            def fake_run(factory,lam,args_in,seed):
                recorded.append((factory.config,lam,seed.copy()))
                row={'model':'A_-4','lambdaT3':lam,'status':'AcceptedApproxFit',
                     'chi2':0.3,'fitted_yukawas':campaign.base.yukawa_matrices(seed)}
                return row,None,None
            with patch.object(campaign,'ROOT',root),\
                 patch.object(campaign.uniform,'ModelFactory',FakeFactory),\
                 patch.object(campaign.uniform,'run_point',fake_run),\
                 patch.object(campaign,'record_point',lambda *args:None):
                rows,summary=campaign.process_model('A_-4','reference',root,args,
                                                     {'A_-4':{}},None,zee)
            self.assertEqual(summary['status'],'Accepted')
            self.assertEqual(rows[0]['stage'],'reference')
            self.assertEqual(len(recorded),1)
            self.assertAlmostEqual(recorded[0][1],0.1)
            np.testing.assert_array_equal(recorded[0][2],orig)

    def test_scan_reuses_new_reference_not_original_yukawas(self):
        old=self.original()
        new,_=uniform.uniform_config(old,0.1)
        new['base_state']['ordinary']['y1']['real'][0][0]=0.007
        new['thesis_uniform_scalar_benchmark']={'schema':uniform.SCHEMA}
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            source=root/'source.json';source.write_text('{}',encoding='utf-8')
            manifest=root/'merged.json';manifest.write_text('{}',encoding='utf-8')
            paths={'config':source}
            fit={'chi2_final':0.1}
            zee=SimpleNamespace(MERGED_FILE=manifest,sources=lambda *args:(old,fit,paths))
            args=SimpleNamespace(canonical_run_id='full16_general_v1',
                                 max_nfev=5,accept_seed=False,chi2_cutoff=5.86)
            with patch.object(campaign,'ROOT',root):
                row={'schema':campaign.SCHEMA,'stage':'reference','model':'A_-4',
                     'status':'AcceptedApproxFit','chi2':0.01,'lambdaT3':0.1,
                     'provenance':campaign.provenance(paths,manifest),
                     'fitted_yukawas':{k:deepcopy(new['base_state']['ordinary'][k])
                                         for k in ('y1','y2')}}
                folder=root/'reference'/'A_m4'
                campaign.base.json_save(folder/'result.json',row)
                campaign.base.json_save(folder/'fitted_config.json',new)
                recorded=[]
                class FakeFactory:
                    def __init__(self,model,*arguments):
                        self.config=None
                def fake_run(factory,lam,arguments,seed):
                    recorded.append((lam,float(seed[0])))
                    return ({'model':'A_-4','lambdaT3':lam,'status':'PoorFit',
                             'chi2':10.0},None,None)
                with patch.object(campaign.uniform,'ModelFactory',FakeFactory),\
                     patch.object(campaign.uniform,'run_point',fake_run),\
                     patch.object(campaign,'record_point',lambda *args:None):
                    rows,_=campaign.process_model('A_-4','scan',root,args,
                                                  {'A_-4':{}},None,zee)
            self.assertEqual([r[0] for r in recorded],list(campaign.base.BENCHMARKS))
            self.assertTrue(all(abs(r[1]-0.007)<1e-14 for r in recorded))
            self.assertTrue(all(row['new_reference_chi2']==0.01 for row in rows))


if __name__=='__main__':
    unittest.main()
