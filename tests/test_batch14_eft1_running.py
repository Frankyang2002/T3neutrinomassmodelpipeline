"""Batch14 F-threshold / EFT1 candidate regression checks.

Project path: tests/test_batch14_eft1_running.py
Run: python -m unittest discover -s tests -p test_batch14_eft1_running.py -v
"""
import json
import tempfile
import unittest
from pathlib import Path
import numpy as np

from Numerical.experimental.Batch13GeneralUVRunning import (
    EXPECTED_CONVENTION, EXPECTED_FORMAT, EXPECTED_STATUS,
    illustrative_initial_state, serialise_state,
)
from Numerical.experimental.Batch14EFT1Running import (
    load_saved_uv_endpoint, run_eft1, validate_eft1_export, serialise_eft1_state,
)
from scripts.RunBatch14EFT1All16 import model_command


def make_payloads():
    # Minimal synthetic inventory; all numerical entries zero-beta.
    names=('gY','g2','g3','yu','yd','ye','lambdaH','y1','y2','MF',
           'mS1Sq','mS2Sq','lambdaS1','lambdaS2','lambdaH1','lambdaH2',
           'lambda12','lambdaT3','lambdaH2Adj','lambdaS2Adj',
           'muT3GeneralS1Sq','lambdaT3GeneralQ000400I1')
    base={'GeneralScalarModel':'A_+0',
          'GeneralScalarBasisStatus':'CandidateNativeRGBetaUnmappedToMatchete',
          'ReportBetaConvention':EXPECTED_CONVENTION,
          'dS1':1,'dS2':3,'dF':2,'alpha':0,
          'GeneralScalarNewBetaKeys':['lambdaT3GeneralQ000400I1'],
          'QuadraticKind':'self_s1','QuadraticMassSymbol':'muT3GeneralS1Sq'}
    def beta(n):
        return '0*'+n if n in ('yu','yd','ye','y1','y2','MF') else '0'
    uv={'format':EXPECTED_FORMAT,'status':EXPECTED_STATUS,
        'physical_mapping_certified':False,
        'metadata':{**base,'GeneralScalarStage':'UV'},
        'report_betas':{n:beta(n) for n in names}}
    eft={'format':EXPECTED_FORMAT,'status':EXPECTED_STATUS,
         'physical_mapping_certified':False,
         'metadata':{**base,'GeneralScalarStage':'EFT1'},
         'report_betas':{n:beta(n) for n in names if n not in ('y1','y2','MF')}}
    return uv,eft


class Batch14Tests(unittest.TestCase):
    def test_command_stage_before_representation(self):
        c=model_command('B_-1',Path('out.json'))
        self.assertEqual(c[3:9],['EFT1','2','2','1','m1','B_-1'])

    def test_eft1_inventory_removes_all_fermions(self):
        uv,eft=make_payloads()
        names=set(validate_eft1_export(eft,uv,'A_+0'))
        self.assertNotIn('MF',names)
        self.assertIn('muT3GeneralS1Sq',names)
        self.assertIn('lambdaT3GeneralQ000400I1',names)

    def test_reject_eft1_heavy_yukawa(self):
        uv,eft=make_payloads()
        eft['report_betas']['y1']='0*y1'
        with self.assertRaisesRegex(ValueError,'heavy fermion'):
            validate_eft1_export(eft,uv,'A_+0')

    def test_reject_missing_quadratic(self):
        uv,eft=make_payloads()
        del eft['report_betas']['muT3GeneralS1Sq']
        with self.assertRaisesRegex(ValueError,'inventory mismatch'):
            validate_eft1_export(eft,uv,'A_+0')

    def test_saved_endpoint_roundtrip_and_zero_beta_eft1(self):
        uv,eft=make_payloads()
        initial=illustrative_initial_state(uv,'A_+0',mu_gev=1e5)
        document={'model':'A_+0','status':'CandidateUVIntegrated',
                  'restricted_yukawa':True,
                  'matchete_rgbeta_normalisation_certified':False,
                  'mu_final_gev':1e5,'final':serialise_state(initial)}
        # Nonzero general scalar quadratic ensures it survives threshold handoff.
        document['final']['muT3GeneralS1Sq']={'real':240.0,'imag':10.0}
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'A_+0_UV_running.json'
            path.write_text(json.dumps(document))
            recovered=load_saved_uv_endpoint(path,uv,'A_+0')
        self.assertEqual(recovered.general_couplings['muT3GeneralS1Sq'],240+10j)
        result=run_eft1(recovered,uv,eft,'A_+0',mu_scalar_gev=1e3)
        self.assertEqual(result['status'],'CandidateFThresholdAndEFT1Integrated')
        self.assertEqual(result['mu_s_gev'],1e3)
        self.assertTrue(result['both_scalars_active_in_eft1'])
        self.assertFalse(result['fermion_threshold_wilson_matching_performed'])
        self.assertEqual(result['state'].general_couplings['muT3GeneralS1Sq'],240+10j)
        self.assertEqual(set(serialise_eft1_state(result['state'])),set(eft['report_betas']))

    def test_endpoint_rejects_missing_complex_coupling(self):
        uv,_=make_payloads()
        initial=illustrative_initial_state(uv,'A_+0',mu_gev=1e5)
        document={'model':'A_+0','status':'CandidateUVIntegrated',
                  'restricted_yukawa':True,
                  'matchete_rgbeta_normalisation_certified':False,
                  'mu_final_gev':1e5,'final':serialise_state(initial)}
        del document['final']['muT3GeneralS1Sq']
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'run.json'
            path.write_text(json.dumps(document))
            with self.assertRaisesRegex(ValueError,'inventory mismatch'):
                load_saved_uv_endpoint(path,uv,'A_+0')


if __name__=='__main__': unittest.main()
