"""Offline tests of the A(-4) full-quartic-projection verification gate."""
import json
from pathlib import Path
import unittest
from tempfile import TemporaryDirectory
from unittest.mock import patch

from scripts import VerifyAminus4GlobalQuarticProjection as mod
from Numerical.running import Aminus4EFT1Couplings as eft


class TestGlobalProjectionGate(unittest.TestCase):
    def test_55_probes_cover_all_real_quartic_monomials(self):
        probes = mod.quadratic_probes()
        self.assertEqual(len(probes),55)
        self.assertEqual(len(set(probes)),55)
        self.assertEqual(sum(len(x)==1 for x in probes),10)
        self.assertEqual(sum(len(x)==2 for x in probes),45)
        for p in probes:
            state=mod.unit_point(p)
            self.assertEqual(set(state),eft.BETA_NAMES)
            self.assertEqual(sum(bool(state[k]) for k in eft.SCALAR_COUPLINGS if k!='lambdaT3') + int(bool(state['lambdaT3'].real)) + int(bool(state['lambdaT3'].imag)), len(p))

    def test_probes_require_all_55_comparisons(self):
        calls=[]
        def mock(point,payload):
            calls.append(point)
            return {'reference_and_native_scalar_loop_match':True,
                    'scalar_loop_projection_closure_relative':1.e-14,
                    'mismatched_coupling_betas':[],
                    'coupling_beta_comparison':{}}
        r=mod.verify_scalar_polynomial({},comparison=mock)
        self.assertTrue(r['native_and_reference_scalar_quadratic_polynomials_agree'])
        self.assertEqual(len(calls),55)
        self.assertFalse(r['failures'])

    def test_mixed_cross_term_failure_is_not_hidden_by_unit_probes(self):
        def mock(point,payload):
            a=point['lambdaH']
            b=point['lambdaS2Adj']
            bad=bool(a and b)
            return {'reference_and_native_scalar_loop_match':not bad,
                    'scalar_loop_projection_closure_relative':1.e-14,
                    'mismatched_coupling_betas':['lambdaS2'] if bad else [],
                    'coupling_beta_comparison':{'lambdaS2':{'reference_scalar_loop':0.,
                                                           'native_rgbeta_scalar_loop':1.}}}
        r=mod.verify_scalar_polynomial({},comparison=mock)
        self.assertFalse(r['native_and_reference_scalar_quadratic_polynomials_agree'])
        self.assertEqual(len(r['failures']),1)
        self.assertEqual(r['failures'][0]['nonzero_unit_couplings'],
                         ['lambdaH','lambdaS2Adj'])

    def test_incomplete_global_quartic_inventory_rejected(self):
        original={'metadata':{'SharedScalar':False,'dS1':1,'dS2':3,'dF':2,
                              'alpha':-4,'YS1':'-2','YS2':'-1',
                              'ReportBetaConvention':eft.BETA_CONVENTION},
                  'status':'Success','report_betas':{k:'0' for k in eft.BETA_NAMES}}
        original['metadata'].update(IntegratedField='F',ActiveBSMFields=['S1','S2'])
        projected={'metadata':{**original['metadata'],
                        'GeneralScalarStage':'EFT1','GeneralScalarModel':'A_-4',
                        'GeneralScalarQuarticBetaBasis':'RGBetaFullQuarticProjection',
                        'GeneralScalarAdditionalQuartics':0,
                        'GeneralScalarExistingStructuralCount':0,
                        'GeneralScalarNewBetaKeys':[],
                        'GeneralScalarExistingQuarticBetasProjected':[]},
                   'status':'Success','report_betas':dict(original['report_betas'])}
        with self.assertRaisesRegex(mod.ProjectionError,'all nine'):
            mod.validate_native(original,projected)
        projected['metadata']['GeneralScalarExistingQuarticBetasProjected']=list(mod.QUARTICS)
        self.assertEqual(len(mod.validate_native(original,projected)['quartic_symbols']),9)
        projected['report_betas'].pop('lambdaS2Adj')
        with self.assertRaisesRegex(mod.ProjectionError,'17 active'):
            mod.validate_native(original,projected)

    def test_coupling_adapter_extracts_native_higgs_beta_dynamically(self):
        """Changing a projected coefficient must not trigger a frozen check."""
        import sys, types
        from unittest.mock import patch
        pkg=types.ModuleType('Numerical.core')
        pkg.__path__=[]
        parser=types.ModuleType('Numerical.core.RGBetaEvaluator')
        parser.evaluate_inputform_expression=lambda expression,p: (
            12*p['lambdaH']**2 + 2*p['lambdaH1']**2 +
            6*p['lambdaH2']**2 + float(expression)*p['lambdaH2Adj']**2 +
            4*abs(p['lambdaT3'])**2)
        meta={'SharedScalar':False,'dS1':1,'dS2':3,'dF':2,'alpha':-4,
              'IntegratedField':'F','ActiveBSMFields':['S1','S2'],
              'ReportBetaConvention':eft.BETA_CONVENTION}
        b={'status':'Success','metadata':meta,
           'report_betas':{k:'0' for k in eft.BETA_NAMES}}
        b['report_betas'].update({'gY':'(55*gY^3)/6',
            'g2':'(-5*g2^3)/2','g3':'-7*g3^3'})
        with patch.dict(sys.modules, {'Numerical.core':pkg,
                         'Numerical.core.RGBetaEvaluator':parser}):
            b['report_betas']['lambdaH']='3'
            legacy=eft.source_normalization_audit(b)
            b['report_betas']['lambdaH']='1'
            corrected=eft.source_normalization_audit(b)
        self.assertEqual(legacy['unresolved_normalization_mismatches'],
                         {'lambdaH2Adj':{'reference':1.,'native':3.}})
        self.assertEqual(corrected['unresolved_normalization_mismatches'],{})
        self.assertTrue(corrected['reference_scalar_beta_lambdaH_agrees'])
        self.assertFalse(corrected['scalar_basis_physically_certified'])

    def test_stage_adapter_is_metadata_only_not_a_physical_certificate(self):
        p={'status':'Success','metadata':{'SharedScalar':False,'dS1':1,'dS2':3,
               'dF':2,'alpha':-4,'ReportBetaConvention':eft.BETA_CONVENTION},
           'report_betas':{k:'0' for k in eft.BETA_NAMES}}
        c=mod.compatibility_payload(p)
        self.assertNotIn('IntegratedField',p['metadata'])
        self.assertEqual(c['metadata']['IntegratedField'],'F')
        self.assertEqual(c['metadata']['ActiveBSMFields'],['S1','S2'])

if __name__=='__main__':unittest.main()
