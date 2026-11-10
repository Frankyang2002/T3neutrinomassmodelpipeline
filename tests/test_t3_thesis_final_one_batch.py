import importlib.util
from pathlib import Path
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/'scripts/RunT3ThesisFinalOneBatch.py'
NATIVE=ROOT/'Lagrangian/RunT3FullGeneralTwoThresholdNative.wl'
C5=ROOT/'Lagrangian/RunT3FullGeneralNativeC5Numerics.wl'
FIT=ROOT/'Numerical/parallel_t3/complete/final_one_batch_fit.py'

spec=importlib.util.spec_from_file_location('one_batch_fit',FIT)
fit=importlib.util.module_from_spec(spec)
spec.loader.exec_module(fit)

class TestOneBatch(unittest.TestCase):
    def test_native_two_thresholds_exact_group_order(self):
        src=NATIVE.read_text(encoding='utf-8')
        self.assertIn('RunT3EFT1WilsonFreshThresholdStage[',src)
        self.assertIn('T3ModelFromClass["A",-4,3]',src)
        self.assertIn('RunSequentialT3Matching[',src)
        self.assertIn('t3InverseIndex',src)
        self.assertIn('A_full_F_hard_one_loop_preserved.wxf',src)
        self.assertIn('Blocked_D2NeedsPhysicalScalarContinuation',src)
    def test_native_requires_full_general_builder(self):
        src=NATIVE.read_text(encoding='utf-8')
        for key in ('AllScalarExtensionsRequested','AdditionalQuarticsIncluded',
                    'RealSelfBilinearSplitLightFields','t3InverseIndex'):
            self.assertIn(key,src)
    def test_no_restricted_hard_triangle_import(self):
        src=SCRIPT.read_text(encoding='utf-8')
        for forbidden in ('RunT3FHardTriangleBatch1','RunT3StrictOneLoopBatch2',
                          'RunT3ThesisTwoModelPilot'):
            self.assertNotIn('import '+forbidden,src)
    def test_fit_requires_eighteen(self):
        target={'central_values':{x:0. for x in fit.OBS},
                'one_sigma_errors':{x:1. for x in fit.OBS}}
        with self.assertRaises(ValueError):
            fit.fit_18_real(np.zeros(17), lambda x:{},target,max_nfev=1)
    def test_residuals_five(self):
        target={'central_values':{x:0. for x in fit.OBS},
                'one_sigma_errors':{x:2. for x in fit.OBS}}
        source={x:1. for x in fit.OBS}
        np.testing.assert_allclose(fit.residual5(source,target),np.full(5,.5))
    def test_fitter_mocked_valid_objective(self):
        target={'central_values':{x:0. for x in fit.OBS},
                'one_sigma_errors':{x:1. for x in fit.OBS}}
        def fake(x):
            return {k:float(x[i]) for i,k in enumerate(fit.OBS)}
        result=fit.fit_18_real(np.full(18,0.01),fake,target,max_nfev=2)
        self.assertEqual(result['status'],'SourceDrivenFitCandidate')
        self.assertEqual(len(result['best']['parameters']),18)
    def test_full_fermion_contraction_not_single_generation(self):
        src=C5.read_text(encoding='utf-8')
        self.assertIn('nF=Lookup[nData,"NF",3]',src)
        self.assertIn('If[nF=!=3',src)
        self.assertIn('{v,1,3}',src)
    def test_pole_cancellation_mandatory(self):
        src=C5.read_text(encoding='utf-8')
        self.assertIn('Blocked_UncancelledDimensionalPole',src)
        self.assertIn('nPoleRel<10^-9',src)
    def test_fitter_source_not_surrogate(self):
        src=FIT.read_text(encoding='utf-8')
        self.assertNotIn('hard_triangle(',src)
        self.assertNotIn('scalar_bubble(',src)
    def test_wolfram_runner_has_no_success_without_physical_inputs(self):
        src=SCRIPT.read_text(encoding='utf-8')
        self.assertIn('physical_full_general_neutrino_predictions_complete',src)
        self.assertIn("'complete_numerical_C5_available': False",src)

if __name__=='__main__': unittest.main()
