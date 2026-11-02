"""Final Workstream3 closeout: actual saved JSON-contract regression tests.

Synthetic source reports mimic the previously inspected WS3 writer schema; this
is not a substitute for the real Windows 16-model run.
"""
from __future__ import annotations

import copy
from contextlib import redirect_stdout
from io import StringIO
import json
import math
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from Numerical.lo_certification.completeness.inventory import ALL_MODELS
from Numerical.lo_certification.completeness.response_closeout import (
    CloseoutError, complex_matrix, inspect_model, inspect_run,
    source_gate, unit_response_norm)
from scripts import CloseLOPortalWorkstream as cli

COUNTS = {
    'A_-4':0, 'A_-2':2, 'A_+0':2, 'B_-3':1,
    'B_-1':4, 'B_+1':1, 'C_-3':1, 'C_-1':4, 'C_+1':1,
    'D_-2':2, 'D_+0':2, 'D_+2':0,
    'E_-4':0, 'E_-2':2, 'E_+0':2, 'E_+2':0,
}


def matrix(a):
    m=np.asarray(a,complex)
    return {'real':m.real.tolist(),'imag':m.imag.tolist()}


def base_data(model):
    a=np.diag([1e-15,2e-15,3e-15]).astype(complex)
    return {'model':model,'status':'CONDITIONAL_COMPLETED_BASELINE_PORTAL',
            'physical_C5_certified':False,'source_normalization_certified':False,
            'all_independent_HHSS_included':False,
            'response':{'projection_passed':True,'physical_C5_certified':False,
              'dim5_norm':float(np.linalg.norm(a)),
              'matrices':{'dim5_sum':matrix(a)}}}


def extra_data(model,count):
    a=np.diag([1e-15,2e-15,3e-15]).astype(complex)
    units=[]
    for i in range(count):
        unitparts={}
        for name,m in (('Re',a),('Im',1j*a)):
            unitparts[name]={'status':'CONDITIONAL_UNIT_REFERENCE_RESPONSE',
              'physical_C5_certified':False,
              'dim5_matrix_projected_diagnostic_only':matrix(m),
              'projection_residuals':{'hard_muF':2e-14,'dim5_sum':3e-14}}
        units.append({'reference_id':f'reference_invariant_{i+1}',
           'unassigned_source_symbols':[f'lambda_reference_{i+1}'],
           'sector_counts':[2,0,0,0,0,i+2],
           'physical_coefficient_assigned':False,
           'reference_basis_to_named_RGBeta_map_certified':False,
           'unit_parts':unitparts})
    return {'model':model,'status':'REFERENCE_UNIT_SENSITIVITIES_NOT_PHYSICAL',
      'named_Matchete_RGBeta_maps_certified':0,
      'physical_coefficients_inserted':False,'physical_C5_certified':False,
      'linear_response_basis':units}


def summary():
    return {'format':'t3-ws3-conditional-portal-completion-and-reference-response-v1',
      'model_count':16,'validated_input_models':16,'blocked_source_models':0,
      'baseline_directions_verified':14,'extra_named_symbols_verified':24,
      'new_reference_unit_response_count':24,
      'conditional_baseline_ward_passed':16,'conditional_baseline_ward_rejected':0,
      'physical_C5_certified':False,'complete_one_loop_C5_certified':False,
      'physical_named_HHSS_coupling_maps_certified':0,
      'source_GammaCC_to_Weyl_signed_factor':None,
      'extra_reference_coefficients_assigned':False,
      'threshold_sequence':'UV -> ALL BSM fermions @ ONE muF -> EFT1(H,S1,S2,SM fermions) -> BOTH scalars @ ONE muS -> SMEFT',
      'models':[{'model':m,'status':'COMPLETED_BASELINE_CONDITIONAL_C5',
        'portal_ward_projection_passed':True,'physical_C5_certified':False,
        'conditional_sensitivity_only':True,
        'extra_named_inventory':{'extra_named_symbols':COUNTS[m],
           'source_named_symbol_evidence':[{'rgbeta_symbol':f'lambda_reference_{i+1}'}
              for i in range(COUNTS[m])]},
        'baseline_completion':{'added_baseline_directions':(2 if m.startswith('E_')
           else 1 if m.startswith(('A_','D_')) else 0)},
        'extra_reference_basis_count':COUNTS[m],
        'unit_reference_responses_computed':COUNTS[m]}
          for m in ALL_MODELS]}


def fixture(directory):
    root=Path(directory)
    (root/'conditional_portal_response').mkdir(parents=True)
    (root/'extra_reference_unit_responses').mkdir()
    (root/'summary.json').write_text(json.dumps(summary()),encoding='utf-8')
    for m in ALL_MODELS:
        slug=m.replace('-','m').replace('+','p')
        (root/'conditional_portal_response'/(slug+'.json')).write_text(
            json.dumps(base_data(m)),encoding='utf-8')
        (root/'extra_reference_unit_responses'/(slug+'.json')).write_text(
            json.dumps(extra_data(m,COUNTS[m])),encoding='utf-8')


class ResponseCloseoutMathTests(unittest.TestCase):
    def test_exact_Frobenius_operator_norm_orthogonal_i(self):
        a=np.diag([1,2,3]).astype(complex)
        self.assertAlmostEqual(unit_response_norm(a,1j*a),float(np.linalg.norm(a)))

    def test_same_phase_two_parts_max_sqrt2(self):
        a=np.diag([1,2,3]).astype(complex)
        self.assertAlmostEqual(unit_response_norm(a,a),math.sqrt(2)*np.linalg.norm(a))

    def test_canceling_phase_gives_same_sqrt2(self):
        a=np.diag([1,2,3]).astype(complex)
        self.assertAlmostEqual(unit_response_norm(a,-a),math.sqrt(2)*np.linalg.norm(a))

    def test_non_symmetric_matrix_rejected(self):
        m=np.array([[0,1,0],[0,0,0],[0,0,0]],complex)
        with self.assertRaises(CloseoutError):complex_matrix(matrix(m),'invalid')

    def test_nan_matrix_rejected(self):
        obj=matrix(np.eye(3))
        obj['real'][0][0]=float('nan')
        with self.assertRaises(CloseoutError):complex_matrix(obj,'invalid')

    def test_status_rejection_does_not_produce_numeric_envelope(self):
        model='A_-2';meta=summary()['models'][1]
        e=extra_data(model,2)
        e['linear_response_basis'][0]['unit_parts']['Re']={
            'status':'WARD_REJECTED_DO_NOT_USE_AS_C5','physical_C5_certified':False}
        out=inspect_model(model,meta,base_data(model),e,(.01,))
        self.assertEqual(out['ward_rejected_reference_directions'],1)
        self.assertIsNone(out['illustrative_max_relative_frobenius_shift_at_caps']['0.01'])
        self.assertEqual(out['references'][0]['status'],'WARD_REJECTED_NO_C5')

    def test_reference_blocked_does_not_get_zero_imputed(self):
        model='B_-3';meta=summary()['models'][3]
        e=extra_data(model,1)
        e['linear_response_basis'][0]['unit_parts']['Im']={
            'status':'REFERENCE_NUMERICAL_BLOCKED','error':'synthetic',
            'physical_C5_certified':False}
        out=inspect_model(model,meta,base_data(model),e,(.01,))
        self.assertEqual(out['other_blocked_reference_directions'],1)
        self.assertIsNone(out['references'][0]['max_dim5_frobenius_per_unit_complex_coefficient'])

    def test_zero_reference_tensor_valid_zero_sensitivity(self):
        model='B_-3';meta=summary()['models'][3]
        e=extra_data(model,1)
        for label in ('Re','Im'):
            e['linear_response_basis'][0]['unit_parts'][label]={
                'status':'IDENTICALLY_ZERO_REFERENCE_TENSOR','physical_C5_certified':False}
        out=inspect_model(model,meta,base_data(model),e,(.01,))
        self.assertEqual(out['valid_reference_directions'],1)
        self.assertEqual(out['max_single_direction_relative_response_per_unit'],0.0)
        self.assertEqual(out['illustrative_max_relative_frobenius_shift_at_caps']['0.01'],0.0)

    def test_one_unit_response_relative_ratio_and_threshold(self):
        model='B_-3';meta=summary()['models'][3]
        out=inspect_model(model,meta,base_data(model),extra_data(model,1),(.001,.01,.1))
        ref=out['references'][0]
        self.assertAlmostEqual(ref['max_relative_to_baseline_per_unit_coefficient'],1)
        self.assertAlmostEqual(ref['single_reference_coefficient_for_10pct_baseline'],.1)
        self.assertAlmostEqual(out['illustrative_max_relative_frobenius_shift_at_caps']['0.01'],.01)

    def test_no_extra_invariant_not_missing(self):
        model='A_-4';meta=summary()['models'][0]
        out=inspect_model(model,meta,base_data(model),extra_data(model,0),(.01,))
        self.assertEqual(out['status'],'NO_EXTRA_INVARIANT_DIRECTIONS')
        self.assertEqual(out['illustrative_max_relative_frobenius_shift_at_caps']['0.01'],0)

    def test_tampered_flags_not_promoted(self):
        model='B_-3';meta=summary()['models'][3]
        ex=extra_data(model,1)
        ex['physical_coefficients_inserted']=True
        with self.assertRaises(CloseoutError):inspect_model(model,meta,base_data(model),ex,(.01,))
        ex=extra_data(model,1)
        ex['linear_response_basis'][0]['physical_coefficient_assigned']=True
        with self.assertRaises(CloseoutError):inspect_model(model,meta,base_data(model),ex,(.01,))
        base=base_data(model);base['physical_C5_certified']=True
        with self.assertRaises(CloseoutError):inspect_model(model,meta,base,extra_data(model,1),(.01,))

    def test_duplicate_reference_sector_disallowed(self):
        model='A_-2';meta=summary()['models'][1]
        e=extra_data(model,2)
        e['linear_response_basis'][1]['reference_id']=e['linear_response_basis'][0]['reference_id']
        e['linear_response_basis'][1]['sector_counts']=e['linear_response_basis'][0]['sector_counts']
        with self.assertRaises(CloseoutError):inspect_model(model,meta,base_data(model),e,(.01,))

    def test_ward_residual_claim_is_reevaluated(self):
        model='B_-3';meta=summary()['models'][3]
        e=extra_data(model,1)
        e['linear_response_basis'][0]['unit_parts']['Re']['projection_residuals']['dim5_sum']=1e-2
        with self.assertRaises(CloseoutError):inspect_model(model,meta,base_data(model),e,(.01,))

    def test_summary_gate_rejects_wrong_counts_sequence_or_physical(self):
        good=summary();source_gate(good)
        for key,val in (('new_reference_unit_response_count',23),
                        ('baseline_directions_verified',13),('physical_C5_certified',True),
                        ('source_GammaCC_to_Weyl_signed_factor',-1),
                        ('threshold_sequence','UV->SMEFT')):
            bad=copy.deepcopy(good);bad[key]=val
            with self.subTest(key=key), self.assertRaises(CloseoutError):source_gate(bad)


class ResponseCloseoutIntegrationTests(unittest.TestCase):
    def test_real_file_schema_all16_24_and_manifest(self):
        self.assertEqual(set(COUNTS),set(ALL_MODELS))
        self.assertEqual(sum(COUNTS.values()),24)
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)/'portal_completion_02';fixture(root)
            report=inspect_run(root)
            self.assertEqual(report['baseline_completed_models'],16)
            self.assertEqual(report['extra_reference_directions'],24)
            self.assertEqual(report['valid_unit_reference_directions'],24)
            self.assertEqual(report['ward_rejected_reference_directions'],0)
            self.assertEqual(report['models_with_extra_references'],12)
            self.assertEqual(len(report['reference_sensitivity_ranking']),24)
            self.assertFalse(report['physical_C5_certified'])
            for row in report['models']:
                self.assertEqual(row['reference_directions'],COUNTS[row['model']])
                self.assertEqual(row['baseline_ward_passed'],True)
                self.assertIn('baseline_response',row['input_sha256'])

    def test_integrated_one_bad_status_yields_partial_not_false_completion(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)/'portal_completion_02';fixture(root)
            path=root/'extra_reference_unit_responses'/'B_m3.json'
            j=json.loads(path.read_text())
            j['linear_response_basis'][0]['unit_parts']['Re']['status']='WARD_REJECTED_DO_NOT_USE_AS_C5'
            path.write_text(json.dumps(j))
            report=inspect_run(root)
            self.assertEqual(report['ward_rejected_reference_directions'],1)
            self.assertEqual(report['valid_unit_reference_directions'],23)
            self.assertEqual(report['status'],'ALL_16_BASELINES_AUDITED_SOME_UNIT_RESPONSES_BLOCKED')

    def test_wrong_expected_sensitivity_count_blocks_report(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)/'portal_completion_02';fixture(root)
            path=root/'extra_reference_unit_responses'/'B_m3.json'
            obj=json.loads(path.read_text());obj['linear_response_basis'].clear()
            path.write_text(json.dumps(obj))
            with self.assertRaises(CloseoutError):inspect_run(root)

    def test_runner_creates_new_complete_closeout_artifacts(self):
        with tempfile.TemporaryDirectory() as td:
            base=Path(td)/'owned';src=base/'portal_completion_02';dst=base/'final_sensitivity_closeout_01'
            fixture(src)
            with patch.object(cli,'OWNED',base), redirect_stdout(StringIO()):
                self.assertEqual(cli.main(['--input-dir',str(src),'--output-dir',str(dst)]),0)
            self.assertTrue((dst/'WORKSTREAM3_CLOSEOUT.md').is_file())
            self.assertEqual(len((dst/'reference_sensitivity.csv').read_text().splitlines()),25)
            self.assertEqual(len((dst/'model_sensitivity.csv').read_text().splitlines()),17)
            self.assertFalse(json.loads((dst/'summary.json').read_text())['physical_C5_certified'])
            with patch.object(cli,'OWNED',base), redirect_stdout(StringIO()):
                self.assertEqual(cli.main(['--input-dir',str(src),'--output-dir',str(dst)]),2)

    def test_check_only_no_writes_and_no_nonowned_output(self):
        with tempfile.TemporaryDirectory() as td:
            base=Path(td)/'owned';src=base/'portal_completion_02';dst=base/'out'
            fixture(src)
            with patch.object(cli,'OWNED',base), redirect_stdout(StringIO()):
                self.assertEqual(cli.main(['--input-dir',str(src),'--output-dir',str(dst),'--check-only']),0)
                self.assertFalse(dst.exists())
                self.assertEqual(cli.main(['--input-dir',str(src),'--output-dir',str(Path(td)/'escape')]),2)


if __name__=='__main__':unittest.main()
