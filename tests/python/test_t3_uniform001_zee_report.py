"""Isolated tests for NEW input/report adapters; no RGE/matching is simulated."""
import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from scripts import RunT3Uniform001ZeeSensitivity as sensitivity
from scripts import BuildT3Uniform001ZeeReport as report

MODELS=tuple(f'{fam}_{a:+d}' for fam, aa in {
    'A':(-4,-2,0),'B':(-3,-1,1),'C':(-3,-1,1),
    'D':(-2,0,2),'E':(-4,-2,0,2)}.items() for a in aa)


class Uniform001ZeeStudyTests(unittest.TestCase):
    def test_model_selection_and_order(self):
        self.assertEqual(sensitivity.selected_models('B_-1,A_-4',MODELS),['A_-4','B_-1'])
        self.assertEqual(len(sensitivity.selected_models('all',MODELS)),16)
        for bad in ('A_-4,A_-4','H_-99','',','):
            with self.assertRaises(ValueError):
                sensitivity.selected_models(bad,MODELS)

    def test_resume_provenance_must_match(self):
        models=['A_-4']; gamma=[.99,1.,1.01]
        record=sensitivity.report_payload(models,gamma,'ref_v1','sha','full16_general_v1',[])
        sensitivity.validate_resume(record,models=models,gamma=gamma,
                                    reference_run_id='ref_v1',reference_sha='sha',
                                    canonical_run='full16_general_v1')
        for field,value in [('uniform_quartic',.05),('conjugate_yukawas_uv_zero',False),
                            ('reference_summary_sha256','bad'),('gamma_grid',[.95,1.,1.05])]:
            corrupt={**record,field:value}
            with self.assertRaises(ValueError):
                sensitivity.validate_resume(corrupt,models=models,gamma=gamma,
                                            reference_run_id='ref_v1',reference_sha='sha',
                                            canonical_run='full16_general_v1')

    def test_checkpoint_csv_roundtrip_and_reject_duplicates(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'points.csv'
            rows=[{'model':'A_-4','direction':'T1','gamma':1.,'status':'Success',
                   'chi2':0.,'delta_m21_sq_ev2':7.537e-5,
                   'delta_m3l_sq_ev2':.002511,'sin2_theta12':.3088,
                   'sin2_theta13':.02248,'sin2_theta23':.47,'error':''}]
            sensitivity.atomic_csv(path,rows)
            self.assertEqual(len(sensitivity.stored_points(path,['A_-4'],[1.])),1)
            sensitivity.atomic_csv(path,rows*2)
            with self.assertRaisesRegex(ValueError,'duplicate'):
                sensitivity.stored_points(path,['A_-4'],[1.])

    def test_model_complete_requires_three_directions(self):
        r={'model':'A_-4','direction':'T1','gamma':1.,'status':'Success'}
        ref={'A_-4':(None,None,None,{'chi2':1e-8,'status':'AcceptedApproxFit'},None)}
        result=sensitivity.summarize(['A_-4'],[1.],{sensitivity.row_key(r):r},ref)
        self.assertEqual(result[0]['status'],'Pending')
        all_points={(r['model'],d,1.):{**r,'direction':d} for d in sensitivity.DIRECTIONS}
        self.assertEqual(sensitivity.summarize(['A_-4'],[1.],all_points,ref)[0]['status'],'Complete')

    def test_reference_table_accepted_and_no_foreign_models(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            summary={'schema':'t3-thesis-uniform-001-reference-and-scan-v1',
                     'stage':'reference','uniform_scalar_quartic_uv':.01,
                     'conjugate_yukawas_uv_zero':True,
                     'model_summaries':[{'model':m,'status':'Accepted'} for m in MODELS]}
            (root/'reference_summary.json').write_text(json.dumps(summary))
            table=root/'reference_comparison.csv'
            rows=[{'model':m,'stage':'reference',
                   'lambdaT3':.05 if m.startswith('E_') else .1,
                   'status':'AcceptedApproxFit','chi2':1e-9,
                   'sin2_theta12':.3088,'delta_m21_sq_ev2':7.537e-5} for m in MODELS]
            def save():
                with table.open('w',newline='') as f:
                    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
            save()
            self.assertEqual(len(report.load_reference_table(table,root/'reference_summary.json',MODELS)),16)
            rows[0]['lambdaT3']=.5;save()
            with self.assertRaisesRegex(ValueError,'wrong reference'):
                report.load_reference_table(table,root/'reference_summary.json',MODELS)
            rows[0]['lambdaT3']=.1;rows[0]['status']='PoorFit';save()
            with self.assertRaises(ValueError):
                report.load_reference_table(table,root/'reference_summary.json',MODELS)

    def test_full_report_refuses_wrong_quartic(self):
        with tempfile.TemporaryDirectory() as tmp:
            src=Path(tmp)/'output';src.mkdir()
            (src/'points.csv').write_text('model,direction,gamma,status\n')
            meta={'schema':sensitivity.SCHEMA,'uniform_quartic':.05,
                  'conjugate_yukawas_uv_zero':True,
                  'models':[{'model':m,'status':'Complete'} for m in MODELS]}
            (src/'summary.json').write_text(json.dumps(meta))
            with patch.object(report,'SENSITIVITY',Path(tmp)):
                with self.assertRaises(SystemExit) as cm:
                    report.main(['--mode','sensitivity','--sensitivity-run-id','output'])
            self.assertEqual(cm.exception.code,2)

    def test_checkpoint_writes_utf8_json(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'summary.json';data=sensitivity.report_payload(['A_-4'],[1.],
                    'ref','sha','full16_general_v1',[])
            sensitivity.atomic_json(path,data)
            self.assertEqual(json.loads(path.read_text()),data)


if __name__=='__main__':
    unittest.main()
