"""Source coverage, all-16 status, strict alias and exclusive-write checks."""
import sys
from pathlib import Path
import tempfile
import unittest
import json
from unittest import mock
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from Numerical.lo_certification.couplings import batch


class TestBatchContracts(unittest.TestCase):
    def test_all_16_models_are_present_in_universal_inventory(self):
        models={f'{fam}_{a:+d}' for fam,vals in {
            'A':(-4,-2,0),'B':(-3,-1,1),'C':(-3,-1,1),
            'D':(-2,0,2),'E':(-4,-2,0,2)}.items() for a in vals}
        self.assertEqual(len(models),16)
        self.assertEqual(sum(2*(dims[0]+dims[1]) for dims in batch.MODEL_DIMS.values()),44)
        self.assertEqual(batch.MODEL_DIMS['C'],(2,2,3))

    def test_generated_registry_extra_HHSS_is_not_claimed_normalized(self):
        registry={'new_rgbeta_symbols':[
          {'name':'Q', 'counts':[0,2,0,1,0,1],
           'matchete_basis_name':'Qnative', 'mapping_status':'NotNormalised'},
          {'name':'SSSS', 'counts':[0,0,2,2,0,0],
           'matchete_basis_name':'SSSSnative', 'mapping_status':'NotNormalised'}],
          'preexisting_rgbeta_symbols':[{'counts':[1,1,0,0,1,1],
                                         'existing_name':'lambdaOther'}]}
        result=batch.nonzero_sources(registry, {'lambdaT3'})
        self.assertEqual([x['rgbeta'] for x in result],['Q','lambdaOther'])
        self.assertTrue(all(x['represented_as_certified_native_HHSS'] is False for x in result))

    def test_aliases_are_candidates_not_physical_certificates(self):
        self.assertEqual(batch.CANDIDATE_ALIAS['lambdaH2Inv2'],'lambdaH2Adj')
        self.assertNotIn('lambdaMissing',batch.CANDIDATE_ALIAS)
        class Channel:
            by_name={'lambdaMissing':{},'lambdaT3':{}}
        a,_,_,missing,_=batch.calculate_candidate_quartics(Channel(),{'lambdaT3':0.2})
        self.assertIsNone(a)
        self.assertEqual(missing,['lambdaMissing'])

    def test_output_writes_exclusively(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'status.json'
            batch.write_json(p,{'status':'blocked'})
            with self.assertRaises(FileExistsError):
                batch.write_json(p,{'status':'verified'})
            self.assertEqual(batch.read_json(p)['status'],'blocked')

    def test_changed_registry_missing_prior_HHSS_blocks(self):
        rows={'A_-2': {'new_rgbeta_symbols': [], 'preexisting_rgbeta_symbols': []}}
        d={'status':'GeneratedCandidateBasisNotWolframVerified','models':rows}
        versions={batch.REGISTRY_PATH: {'matches_inspected_main':False}}
        with self.assertRaisesRegex(ValueError, 'prior HHSS registry invariants disappeared'):
            batch.local_source_recheck(Path('.'),('A_-2',),versions,Path('.'),Path('.'),d)

    def test_new_registry_candidate_remains_explicitly_unverified(self):
        rows={'A_-4': {'new_rgbeta_symbols': [
             {'name':'newCandidate','counts':[1,1,1,1,0,0]}],
             'preexisting_rgbeta_symbols': []}}
        d={'status':'GeneratedCandidateBasisNotWolframVerified','models':rows}
        versions={batch.REGISTRY_PATH: {'matches_inspected_main':False}}
        result=batch.local_source_recheck(Path('.'),('A_-4',),versions,Path('.'),Path('.'),d)
        self.assertEqual(result['registry_new_HHSS_candidates']['A_-4'],['new:newCandidate'])
        self.assertFalse(result['physical_source_certified'])

    def test_unreviewed_other_source_change_is_blocked(self):
        versions={'RGE/general/WilsonTensorRGE.py':{'matches_inspected_main':False}}
        with self.assertRaisesRegex(ValueError,'Unreviewed source changes'):
            batch.local_source_recheck(Path('.'),('A_-4',),versions,Path('.'),Path('.'),{})

    def test_git_blob_sha_is_true_git_blob(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'test.txt';p.write_bytes(b'abc\n')
            self.assertEqual(batch.git_blob_id(p),'8baef1b4abc478178b004d62031cf7fe6db6f903')

    def test_running_missing_quartic_channel_is_hard_failure(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'run.json'
            example={'model':'A_-4','status':'NumericalRenormalisableRunningOnly',
              'illustrative_benchmark':True,'physical_mass_parameters_available':False,
              'eft1_wilson_running_performed':False,
              'scales_gev':{'UV':1e7,'F':1e5,'S':1e3},
              'eft1_at_fermion_threshold':{'mS1Sq':1.,'mS2Sq':2.,
                 'lambdaT3':1.,'gY':.3,'g2':.6,'g3':1.,'ye':[], 'yu':[], 'yd':[]},
              'eft1_at_scalar_threshold':{'mS1Sq':1.,'mS2Sq':2.,
                 'lambdaT3':1.,'gY':.3,'g2':.6,'g3':1.,'ye':[], 'yu':[]}}
            p.write_text(json.dumps(example))
            with self.assertRaisesRegex(ValueError,'missing EFT1 running field yd'):
                batch.observed_running_states(p,'A_-4')

    def test_all_16_requires_actual_source_and_does_not_create_output(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            output=root/'output/general_observable/lo_certification_couplings/test'
            with self.assertRaises((ImportError,FileNotFoundError,ValueError)):
                batch.run_batch(root,output)
            self.assertFalse(output.exists())

    def test_wrong_output_owner_is_rejected_first(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            with self.assertRaisesRegex(ValueError,'owned'):
                batch.run_batch(root,root/'RGE/general/unowned')

    def test_no_automatic_physical_spinor_factor(self):
        self.assertNotIn('GammaCC_to_master_Weyl_factor',batch.CANDIDATE_ALIAS)


if __name__=='__main__':
    unittest.main()
