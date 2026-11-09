"""v17 synthetic 150-slot census, negative status and source guards."""
from __future__ import annotations
import copy
import tempfile
from pathlib import Path
import unittest
from unittest.mock import patch
from Numerical.parallel_t3.wilson.nonportal_native_dictionary_v16 import (
    audit,V16Blocked,EXPECTED_MODELS,DIMENSIONS,
)
from Numerical.parallel_t3.wilson.full_scalar_census_v17 import audit_from_paths,STATUS
from Numerical.parallel_t3.wilson.portal_real_action_v9 import STATUS as V9_STATUS, FORMAT as V9_FORMAT
from Numerical.parallel_t3.wilson.native_context_quartic_v14 import STATUS as V14_STATUS
from test_t3_parallel_a_nonportal_dictionary_v16 import fake_all

class TestV17(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        m,p=fake_all();cls.v16=audit(m,p,samples=1)
        cls.v9={'format':V9_FORMAT,'status':V9_STATUS,
                'native_Matchete_portal_coupling_slots':58,'models_checked':16,
                'certification':{'physical_Wilson_evolution_authorized':False},'models':{}}
        for t in EXPECTED_MODELS:
            d1,d2,_=DIMENSIONS[t[0]]
            cls.v9['models'][t]=[
                {'rgbeta_registered_couplings':[f'lambdaH{i}']+([f'lambdaH{i}Adj'] if d>1 else [])}
                for i,d in ((1,d1),(2,d2))]
        cls.v14={'status':V14_STATUS,'input_sha256':{'a':'b'},
                 'source_census':{'real_registered':150,'real_diagonal_unit':150,
                                  'portal_response_numeric':110},
                 'real_projector_conditionally_verified':True,
                 'physical_Wilson_evolution_authorized':False,
                 'full_scalar_potential_dictionary_certified':False,
                 'native_numeric_S2_frame_verified':False}
    def run_audit(self,override=None):
        v16=copy.deepcopy(self.v16);v9=copy.deepcopy(self.v9);v14=copy.deepcopy(self.v14)
        if override:override(v16,v9,v14)
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            for s in ('native','portals','rgbeta','v16','v9','v14'):(root/s).write_text('{}',encoding='utf8')
            def reader(p):return {'v16':v16,'v9':v9,'v14':v14}[Path(p).name]
            with patch('Numerical.parallel_t3.wilson.full_scalar_census_v17.audit_v16',return_value=self.v16), \
                 patch('Numerical.parallel_t3.wilson.full_scalar_census_v17.audit_v14',return_value=self.v14), \
                 patch('Numerical.parallel_t3.wilson.full_scalar_census_v17.read',side_effect=reader):
                return audit_from_paths(root,*(root/s for s in ('native','portals','rgbeta','v16','v9','v14')))

    def test_full_real_scalar_150_census(self):
        r=self.run_audit()
        self.assertEqual((r['Matchete_Higgs_quartic_slots'],r['Matchete_Higgs_portal_slots'],
                          r['Matchete_nonportal_scalar_slots']), (16,58,76))
        self.assertEqual(sum(x['total_real_quartic_slots'] for x in r['models'].values()),150)
        self.assertEqual(r['status'],STATUS)
        self.assertFalse(r['certification']['physical_Wilson_evolution_authorized'])

    def test_nonphysical_wilson_claim_in_saved_v16_rejected(self):
        def change(v16,v9,v14):v16['certification']['physical_Wilson_evolution_authorized']=True
        with self.assertRaises(V16Blocked):self.run_audit(change)

    def test_false_physics_in_saved_v14_rejected(self):
        def change(v16,v9,v14):v14['physical_Wilson_evolution_authorized']=True
        with self.assertRaises(V16Blocked):self.run_audit(change)

    def test_tampered_saved_nonportal_matrix_rejected(self):
        def change(v16,v9,v14):v16['models']['E_-2']['cross']['conditional_rgbeta_from_matchete_real_matrix'][0][1]+=1
        with self.assertRaises(V16Blocked):self.run_audit(change)

    def test_tampered_v9_registration_census_rejected(self):
        def change(v16,v9,v14):v9['native_Matchete_portal_coupling_slots']=57
        with self.assertRaises(V16Blocked):self.run_audit(change)

    def test_fake_real_projector_150_claim_rejected(self):
        def change(v16,v9,v14):v14['source_census']['real_diagonal_unit']=149
        # Current revalidated audit is unaffected by a saved change, so it must be caught.
        with self.assertRaises(V16Blocked):self.run_audit(change)

    def test_false_physics_claim_rejected(self):
        def change(v16,v9,v14):v9['certification']['physical_Wilson_evolution_authorized']=True
        with self.assertRaises(V16Blocked):self.run_audit(change)

    def test_stale_native_summary_rejected(self):
        def change(v16,v9,v14):v14['input_sha256']['a']='changed'
        with self.assertRaises(V16Blocked):self.run_audit(change)

if __name__=='__main__':unittest.main()

class TestV16V17ResumeCLIs(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        m,p=fake_all()
        cls.report16=audit(m,p,samples=1)
        cls.report16['input_sha256']={'source':'synthetic'}
        cls.report17={
            'format':'t3-parallel-a-full-real-scalar-census-conditional-v17',
            'status':STATUS,'input_sha256':{'source':'synthetic'},
            'certification':{'physical_Wilson_evolution_authorized':False},
            'models':{'synthetic':{'total_real_quartic_slots':150}},
            'native_RGBeta_real_projector_unit_diagonals':150,
        }

    def test_v16_resume_verifies_existing_file_unchanged(self):
        import json
        from scripts.RunT3ParallelANonportalDictionaryV16 import main
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);out=root/'out.json'
            out.write_text(json.dumps(self.report16),encoding='utf-8')
            old=out.read_bytes()
            with patch('scripts.RunT3ParallelANonportalDictionaryV16.audit_from_paths',return_value=self.report16):
                status=main(['--project-root',str(root),'--output',str(out),'--resume'])
            self.assertEqual(status,0);self.assertEqual(old,out.read_bytes())

    def test_v16_resume_modified_matrix_fails_closed(self):
        import json
        from scripts.RunT3ParallelANonportalDictionaryV16 import main
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);out=root/'out.json'
            bad=copy.deepcopy(self.report16)
            bad['models']['E_-2']['cross']['conditional_rgbeta_from_matchete_real_matrix'][0][1]+=1.
            out.write_text(json.dumps(bad),encoding='utf-8')
            with patch('scripts.RunT3ParallelANonportalDictionaryV16.audit_from_paths',return_value=self.report16):
                self.assertEqual(main(['--project-root',str(root),'--output',str(out),'--resume']),7)

    def test_v17_resume_revalidates_existing_census(self):
        import json
        from scripts.RunT3ParallelAFullScalarCensusV17 import main
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);out=root/'census.json'
            out.write_text(json.dumps(self.report17),encoding='utf-8')
            old=out.read_bytes()
            with patch('scripts.RunT3ParallelAFullScalarCensusV17.audit_from_paths',return_value=self.report17):
                self.assertEqual(main(['--project-root',str(root),'--output',str(out),'--resume']),0)
            self.assertEqual(old,out.read_bytes())

    def test_v17_resume_rejects_physics_authorization(self):
        import json
        from scripts.RunT3ParallelAFullScalarCensusV17 import main
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);out=root/'census.json'
            bad=copy.deepcopy(self.report17)
            bad['certification']['physical_Wilson_evolution_authorized']=True
            out.write_text(json.dumps(bad),encoding='utf-8')
            with patch('scripts.RunT3ParallelAFullScalarCensusV17.audit_from_paths',return_value=self.report17):
                self.assertEqual(main(['--project-root',str(root),'--output',str(out),'--resume']),7)
