"""Four focused Batch 3 final-release integrity and physics-scope tests."""
from __future__ import annotations
import copy
import math
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))

# Packaging-only isolated test container lacks the full repository. On the
# user's Windows checkout the actual Numerical.direct_one_loop.core is used.
# Stub only the simple B0/encoding interfaces for isolated release tests.
try:
    from Numerical.direct_one_loop import core as _real_core
except ImportError:
    pack=types.ModuleType('Numerical.direct_one_loop')
    core=types.ModuleType('Numerical.direct_one_loop.core')
    def b0_finite(mf2,ms2,mu2):
        r=ms2/mf2
        return 1-math.log(mf2/mu2)+r*math.log(r)/(1-r)
    def encoded(x):
        a=np.asarray(x,complex)
        return {'real':a.real.tolist(),'imag':a.imag.tolist()}
    core.b0_finite=b0_finite
    core.encoded=encoded
    pack.core=core
    sys.modules['Numerical.direct_one_loop']=pack
    sys.modules['Numerical.direct_one_loop.core']=core

from scripts import BuildT3Batch3FinalRelease as b3


def fake_batch2() -> dict:
    c=np.diag([1e-15,2e-15,3e-15]).astype(complex)
    m=c*b3.VEV_TO_EV
    rows=[]
    for model in b3.MODELS:
        rows.append({'model':model,'status':b3.SUCCESS,
                     'physical_C5_certified':False,
                     'scalar_case':'mixed' if model=='B_-1' else 'ordinary',
                     'conditional_C5_per_N_GeV_inv':b3.direct_core.encoded(c),
                     'conditional_mnu_per_N_eV':b3.direct_core.encoded(m),
                     'mass_singular_values_eV_per_abs_N':np.linalg.svd(m,compute_uv=False).tolist(),
                     'mu_f_gev':1e5,'mu_s_gev':1e3,
                     'scalar_mass_inputs_gev2':{'mS1Sq':1e6,'mS2Sq':1.8e6},
                     'scalar_eigenstate_masses_gev':[1000,1300],
                     'diagnostics':{'all_entries_finite':True,
                                    'SU2_source_beta_parity':True,
                                    'weinberg_projection_residual':1e-15,
                                    'matching_scale_cancellation_residual':1e-15}})
    return {'format':b3.BATCH2_FORMAT,'label':b3.LABEL,
            'models_checked':16,'models_with_conditional_C5':16,
            'models_blocked':0,'physical_C5_certified':False,
            'absolute_GammaCC_to_Weyl_factor':None,
            'complete_one_loop_neutrino_mass_prediction':False,
            'mu_f_gev':1e5,'mu_s_gev':1e3,'models':rows}


def fake_uv(model):
    z=np.diag([3e-11,2e-11,1e-11]).astype(complex)
    return {'model':model,'status':'uv_mass_evaluated_unverified_source_mapping',
            'mass_matrix_gev':b3.direct_core.encoded(z),
            'direct_scale_gev':1e7,
            'mass_singular_values_ev':sorted(np.linalg.svd(z,compute_uv=False)*1e9),
            'max_scalar_squared_over_fermion_squared':1e-4,
            'relative_uv_pole_coefficient':0,
            'leading_heavy_mass_B0_tail':{
                'tail_frobenius_ev':0.0001,
                'tail_over_exact_direct_mass_frobenius':0.0003}}


class Batch3ReleaseTests(unittest.TestCase):
    def test_batch2_16_model_matrices_and_ward_are_revalidated(self):
        b=fake_batch2()
        self.assertEqual(len(b3.validate_batch2(b)),16)
        bad=copy.deepcopy(b)
        bad['models'][0]['conditional_mnu_per_N_eV']['real'][0][0]*=2
        with self.assertRaisesRegex(ValueError,'vev conversion'):
            b3.validate_batch2(bad)
        bad2=copy.deepcopy(b)
        bad2['models'][3]['diagnostics']['weinberg_projection_residual']=1e-2
        with self.assertRaisesRegex(ValueError,'diagnostics not passed'):
            b3.validate_batch2(bad2)

    def test_direct_uv_exact_minus_heavy_expansion_b0_tail(self):
        mf,ms,mu=100000.,1000.,10000000.
        r=(ms/mf)**2
        exact=b3.direct_core.b0_finite(mf*mf,ms*ms,mu*mu)
        amp=np.diag([1e-12,2e-12,3e-12]).astype(complex)
        direct={'renormalization_scale_gev':mu,
                'total':{'complex_symmetric_matrix_gev':b3.direct_core.encoded(amp)},
                'individual_terms':[{'fermion_mass_gev':mf,'scalar_mass_gev':ms,
                                     'matrix_gev':b3.direct_core.encoded(amp)}]}
        tail=b3.tail_of_uv_B0_expansion(direct)
        predicted=amp*(r*r*math.log(r)/(1-r))/exact
        np.testing.assert_allclose(b3.matrix(tail['B0_difference_exact_minus_leading_mass_matrix_gev'],'tail'),
                                   predicted,rtol=1e-12,atol=1e-32)
        self.assertEqual(tail['B0_terms'],1)
        self.assertLess(tail['expansion_parameter_max_ms2_over_mf2'],0.001)

    def test_comparison_refuses_absolute_normalization_or_cross_point_ratios(self):
        eft=b3.validate_batch2(fake_batch2())['B_-1']
        response=b3.compare('B_-1',eft,fake_uv('B_-1'))
        self.assertEqual(response['cross_method_comparison_status'],
                         'not_comparable_unmatched_parameter_point_and_conventions')
        self.assertIsNone(response['relative_uv_vs_eft_matrix_difference'])
        self.assertIsNone(response['ratio_uv_vs_eft_mass_singular_values'])
        self.assertIsNone(response['matrix_difference_gev'])
        self.assertEqual(len(response['blockers_to_direct_physical_comparison']),8)

    def test_final_release_preserves_full_16_and_records_two_or_partial_uv(self):
        import json
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            b2=root/'batch2';b2.mkdir()
            (b2/'summary.json').write_text(json.dumps(fake_batch2()))
            (b2/'all16_conditional_mass_summary.csv').write_text('model,status\n')
            with patch.object(b3,'evaluate_uv',side_effect=lambda model,*a,**kw:fake_uv(model)):
                out=b3.build_release(b2,root,root,root/'out',root)
            self.assertEqual(out['direct_uv_evaluated'],2)
            self.assertEqual(out['physically_same_point_UV_EFT_comparisons'],0)
            self.assertEqual(len(out['all16_conditional_results']),16)
            self.assertEqual(len(out['direct_vs_eft_comparisons']),2)
            with patch.object(b3,'evaluate_uv',side_effect=ValueError('source vertex blocked')):
                blocked=b3.build_release(b2,root,root,root/'out',root)
            self.assertEqual(blocked['direct_uv_evaluated'],0)
            self.assertEqual(blocked['direct_uv_blocked'],2)
            self.assertTrue(all('source vertex blocked' in x['blocker']
                                for x in blocked['direct_vs_eft_comparisons']))

if __name__=='__main__':unittest.main()
