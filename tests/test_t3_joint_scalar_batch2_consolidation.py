"""Focused synthetic contract tests for the read-only Batch 2 consolidator.

These tests do NOT certify the underlying T3 physics or its source parameters.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from scripts.ConsolidateT3JointScalarBatch2 import (
    MODELS, SUCCESS, LABEL, VEV_TO_EV, V6_FORMAT, V6_RESULT_FORMAT, CSV_COLS,
    V6_INNER_FORMAT, consolidate, csv_row, check_v6_input
)


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, allow_nan=False, indent=2)+'\n', encoding='utf-8')


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def cm(z):
    return {'real': np.asarray(z).real.tolist(), 'imag': np.asarray(z).imag.tolist()}


def fixture(root, failing=None):
    """Purely synthetic 16-entry v6-shaped fixture, no inferred couplings."""
    rows = []
    for k, model in enumerate(MODELS):
        name = model.replace('-', 'm').replace('+', 'p')
        category = ('mixed' if model in ('B_-1', 'C_-1') else
                    'self_s1' if model in ('A_+0', 'D_+0', 'E_+0') else
                    'self_s2' if model in ('A_-2', 'D_-2', 'E_-2') else 'ordinary')
        boundary = 'a'*64
        bg_file = root / 'numeric_backgrounds' / (name+'.json')
        npz_file = root / 'numeric_backgrounds' / (name+'.npz')
        bg = {'model':model, 'physical_C5_certified':False,
              'no_separate_scalar_threshold':True,
              'same_muS_quartic_and_mass':True,
              'source_boundary_sha256':boundary, 'mu_s_gev':1000.,
              'tensor_npz':name+'.npz'}
        write_json(bg_file, bg)
        npz_file.write_bytes(b'synthetic background tensor bytes, not a physics input')
        m = np.diag([1e-15*(k+1), 0.8e-15*(k+1), 0.7e-15*(k+1)]).astype(complex)
        # Exercise a symmetric complex off-diagonal entry (not Hermitian).
        m[0,1] = m[1,0] = 0.01e-15j*(k+1)
        lead = 1.2*m
        bubble = -0.2*m
        evals = np.linalg.svd(m, compute_uv=False)
        metrics = {
            'conditional_C5_singular_values_per_N_gev_inverse':evals.tolist(),
            'conditional_neutrino_mass_singular_values_eV_per_abs_N':(evals*VEV_TO_EV).tolist(),
            'conditional_C5_frobenius_per_N_gev_inverse':float(np.linalg.norm(m)),
            'log_C5_frobenius_per_N_gev_inverse':float(np.linalg.norm(lead)),
            'finite_bubble_C5_frobenius_per_N_gev_inverse':float(np.linalg.norm(bubble)),
            'projection_residual':1e-14, 'matching_scale_cancellation_residual':1e-14,
        }
        parity = {'full_SU2_unordered_pair_beta_parity':True}
        inner = {'format':V6_INNER_FORMAT, 'model':model,
                 'model_source_boundary_sha256':boundary,
                 'background_json_sha256':digest(bg_file),
                 'background_npz_sha256':digest(npz_file),
                 'mu_f_gev':1e5, 'mu_s_gev':1000.,
                 'scalar_eigenstate_masses_gev':[1400.,1600.],
                 'physical_C5_certified':False,
                 'weinberg_beta_projection_residual':1e-14,
                 'weinberg_bubble_projection_residual':1e-14,
                 'relative_matching_scale_cancellation_residual':1e-14,
                 'C5_LLSS_total_one_loop_per_N_gev_inverse':cm(m),
                 'C5_EFT1_induced_leading_log_per_N_gev_inverse':cm(lead),
                 'C5_joint_scalar_finite_bubble_per_N_gev_inverse':cm(bubble)}
        response = {'format':V6_RESULT_FORMAT,'model':model,'response':inner,
                    'physical_C5_certified':False,
                    'complete_one_loop_C5_certified':False,
                    'native_quartic_alias_mapping_certified':False,
                    'muF_muS_mass_eigenfield_alignment_certified':False,
                    'full_SU2_beta_source_parity':parity,
                    'numerical_metrics':metrics}
        result_file = root/'C5_per_unknown_NGamma_exploratory'/(name+'.json')
        write_json(result_file,response)
        extra = {'muS_mass_in_fixed_muF_source_basis':True,
                 'source_mu2':{'real':2300.,'imag':200.}} if category!='ordinary' else {}
        rows.append({'model':model,'status':SUCCESS,
            'fermion_boundary_sha256':boundary,
            'running_snapshot_sha256':'b'*64,
            'source_quartic_seed_sha256':'c'*64,
            'candidate_coupling_map':{'lambdaT3':'lambdaT3'},
            'scalar_mass_parameters_gev2':{'mS1Sq':2e6,'mS2Sq':3e6},
            'scalar_case':category,
            'mass_source':'synthetic only',
            'mu_f_gev':1e5,'mu_s_gev':1000.,
            'v2_pivot_verified':True,'v3_full_SU2_verified':True,
            'full_su2_beta_source_parity':parity,
            'numerical_metrics':metrics,'background':extra})
    if failing is not None:
        rows[1]['status'] = 'ERROR_unverified_source_or_physics_consistency'
        rows[1]['error'] = failing
    write_json(root/'summary.json',{
        'format':V6_FORMAT,'status':'PartialConditionalC5Exploratory' if failing else 'All16ConditionalC5Exploratory',
        'models_checked':16,'models_with_conditional_C5':15 if failing else 16,
        'errors':1 if failing else 0,
        'physical_C5_certified':False,'full_one_loop_neutrino_mass_calculation':False,
        'absolute_GammaCC_to_Weyl_normalization':None,
        'allow_exploratory_map':True,'allow_exploratory_special_bases':True,
        'all_scalars_integrated_at_one_muS':True,'models':rows})


class Batch2ConsolidationTests(unittest.TestCase):
    def test_all16_matrices_and_complex_singular_values(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            fixture(root)
            doc = consolidate(root)
            self.assertEqual(doc['models_with_conditional_C5'],16)
            self.assertEqual(doc['models_blocked'],0)
            self.assertEqual(sum(v['valid_models'] for v in doc['scalar_case_comparison'].values()),16)
            self.assertEqual(doc['label'],LABEL)
            self.assertTrue(all(m['physical_C5_certified'] is False for m in doc['models']))
            self.assertTrue(all(len(m['mass_singular_values_eV_per_abs_N'])==3
                                for m in doc['models']))
            self.assertEqual(len(csv_row(doc['models'][0])),len(CSV_COLS))
            self.assertEqual(doc['models'][0]['conditional_C5_per_N_GeV_inv']['imag'][0][1],1e-17)

    def test_missing_successful_file_stays_error_not_zero_mass(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            fixture(root)
            (root/'C5_per_unknown_NGamma_exploratory'/'A_m2.json').unlink()
            doc=consolidate(root)
            self.assertEqual(doc['models_with_conditional_C5'],15)
            failed=doc['models'][1]
            self.assertEqual(failed['status'],'ERROR_batch2_input_validation')
            self.assertNotIn('conditional_C5_per_N_GeV_inv',failed)
            self.assertIn('required v6 input absent',failed['blocker'])

    def test_preexisting_v6_failure_remains_explicit_blocker(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            fixture(root, failing='demonstrated Ward blocker')
            doc=consolidate(root)
            self.assertEqual(doc['models_with_conditional_C5'],15)
            self.assertEqual(doc['models'][1]['blocker'],'demonstrated Ward blocker')
            self.assertNotIn('mass_singular_values_eV_per_abs_N',doc['models'][1])

    def test_tampered_C5_does_not_enter_final_summary(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            fixture(root)
            path=root/'C5_per_unknown_NGamma_exploratory'/'E_m2.json'
            record=json.loads(path.read_text())
            record['response']['C5_LLSS_total_one_loop_per_N_gev_inverse']['real'][0][0]*=2
            write_json(path,record)
            doc=consolidate(root)
            failed=next(m for m in doc['models'] if m['model']=='E_-2')
            self.assertEqual(failed['status'],'ERROR_batch2_input_validation')
            self.assertIn('C5 != log',failed['blocker'])

if __name__=='__main__':
    unittest.main()
