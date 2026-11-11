"""Pure numerical/contract tests; not a substitute for Matchete certificates.

Intended path: tests/test_batch17_conditional_wilson.py
"""
import unittest
import numpy as np
from Numerical.experimental.Batch17ConditionalWilson import (
    MODELS,check_saved_pipeline,mass_parity,
)


def inputs(model='A_-4'):
    handoff={
        'model':model,'status':'CandidateScalarThresholdHandoffAudited',
        'physical_neutrino_mass_prediction':False,
        'both_scalars_integrated_out_in_renormalisable_sm_projection':True,
        'mu_f_gev':1e5,'mu_s_gev':1e3,
        'scalar_threshold':{'additional_quadratic_symbols':[]},
    }
    mass={
        'model':model,'status':'CandidateNativeQuadraticHessianComputed',
        'physical_neutrino_mass_prediction':False,
        'quadratic_kind':'none','quadratic_parameter':None,
        'scalar_CG_basis_mapping_certified':False,
        'conventions':{'rgbeta_potential_plus':{
            'tachyon_free_candidate':True,
            'eigenmass_threshold_hierarchy_candidate':True,
            'mass_squared_eigenvalues_gev2':[4e6]*8}},
    }
    encoded=lambda a:{'real':np.asarray(a,dtype=float).tolist(),
                      'imag':np.zeros_like(np.asarray(a,dtype=float)).tolist()}
    final={n:encoded(x) for n,x in {
        'gY':0.35,'g2':0.64,'g3':1.1,'lambdaH':0.13,
        'yu':np.diag([1e-5,.005,.8]),'yd':np.eye(3)*.01,
        'ye':np.eye(3)*.0001}.items()}
    endpoint={'model':model,'status':'CandidateFThresholdAndEFT1Integrated',
              'physical_neutrino_mass_prediction':False,
              'both_scalars_active_in_eft1':True,
              'mu_f_gev':1e5,'mu_s_gev':1e3,'final':final}
    return handoff,mass,endpoint


class Batch17Checks(unittest.TestCase):
    def test_16_model_inventory(self):
        self.assertEqual(len(MODELS),16)
        self.assertEqual(len(set(MODELS)),16)

    def test_valid_staged_audits(self):
        f,s,sm=check_saved_pipeline('A_-4',*inputs())
        self.assertEqual((f,s),(1e5,1e3))
        self.assertEqual(sm['yu'].shape,(3,3))
        self.assertAlmostEqual(sm['gY'],0.35)

    def test_crossed_model_rejected(self):
        h,m,e=inputs();m['model']='B_-1'
        with self.assertRaisesRegex(ValueError,'model/status'):
            check_saved_pipeline('A_-4',h,m,e)

    def test_missing_quadratic_inventory_rejected(self):
        h,m,e=inputs();h['scalar_threshold']['additional_quadratic_symbols']=['muT3GeneralS1Sq']
        with self.assertRaisesRegex(ValueError,'Quadratic inventory'):
            check_saved_pipeline('A_-4',h,m,e)

    def test_false_physics_certificate_rejected(self):
        h,m,e=inputs();e['physical_neutrino_mass_prediction']=True
        with self.assertRaisesRegex(ValueError,'physical neutrino'):
            check_saved_pipeline('A_-4',h,m,e)

    def test_both_scalars_required(self):
        h,m,e=inputs();e['both_scalars_active_in_eft1']=False
        with self.assertRaisesRegex(ValueError,'both scalar'):
            check_saved_pipeline('A_-4',h,m,e)

    def test_mass_spectrum_parity(self):
        _,m,_=inputs()
        x=mass_parity(m,np.eye(8)*4e6)
        self.assertAlmostEqual(x['spectrum_relative_max_difference'],0)
        self.assertFalse(x['mass_basis_eigenvectors_identified'])

    def test_mismatched_old_background_rejected(self):
        _,m,_=inputs()
        with self.assertRaisesRegex(ValueError,'inconsistent with Batch16'):
            mass_parity(m,np.eye(8)*5e6)

    def test_nonpositive_mass_rejected(self):
        _,m,_=inputs()
        with self.assertRaisesRegex(ValueError,'tachyonic'):
            mass_parity(m,np.eye(8)*-1)

    def test_nonreal_gauge_coupling_rejected(self):
        h,m,e=inputs();e['final']['gY']['imag']=.1
        with self.assertRaisesRegex(ValueError,'nonreal'):
            check_saved_pipeline('A_-4',h,m,e)

if __name__=='__main__':unittest.main()
