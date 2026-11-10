"""Offline checks for full-general consolidated matching invariants.

These tests do NOT certify Matchete matching, renormalisation, Wilson ADM,
finite C5, or neutrino observables. Those require an actual Wolfram run.
"""
from pathlib import Path
import importlib.util
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/RunT3ThesisFinalOneBatch.py"
NATIVE = ROOT / "Lagrangian/RunT3FullGeneralTwoThresholdNative.wl"
PHYSICAL_D = ROOT / "Lagrangian/RunT3DMinus2FullGeneralPhysicalScalar.wl"
C5 = ROOT / "Lagrangian/RunT3FullGeneralNativeC5Numerics.wl"
FIT = ROOT / "Numerical/parallel_t3/complete/final_one_batch_fit.py"

spec = importlib.util.spec_from_file_location("one_batch_fit", FIT)
fit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fit)


class TestOneBatch(unittest.TestCase):
    def test_native_fermions_and_joint_scalar_threshold(self):
        src = NATIVE.read_text(encoding="utf-8")
        self.assertIn('T3ModelFromClass["A", -4, 3]', src)
        self.assertIn('RunSequentialT3Matching[', src)
        self.assertIn('RunT3EFT1WilsonFreshThresholdStage[', src)
        self.assertIn('{"S1", "S2"}, {"S1", "S2"}', src)
        self.assertIn('A_full_F_hard_one_loop_preserved.wxf', src)
        self.assertIn('ScalarSoftOnly_HardInheritanceAndWilsonRGPending', src)

    def test_a_parent_does_not_prematurely_register_wilson_couplings(self):
        src = NATIVE.read_text(encoding="utf-8")
        self.assertIn('t3InverseIndex', src)
        self.assertIn('t3YukawaFactor', src)
        self.assertIn('t3WCount =!= 2', src)
        # DefineCoupling must happen only in the fresh stage child, not
        # the parent whose registry has not been initialised for EFT1.
        self.assertNotIn('DefineCoupling[T3EFT1W12', src)
        self.assertIn('Coupling[T3EFT1W12', src)

    def test_d_routes_to_split_physical_scalar_continuation(self):
        src = NATIVE.read_text(encoding="utf-8")
        self.assertIn('RunT3DMinus2FullGeneralPhysicalScalar.wl', src)
        d = PHYSICAL_D.read_text(encoding="utf-8")
        self.assertIn('"RealSelfBilinearSplitLightFields"', d)
        self.assertIn('t3DTreePath', d)
        self.assertIn('t3DHardPath', d)
        self.assertIn('"heavy_flavour_count", 0] =!= 3', d)
        self.assertIn('"baseline" -> 4, "extra" -> 4', d)
        self.assertIn('T3EFT1WExtra', d)
        self.assertIn('T3Scalar2Minus', d)
        self.assertIn('T3Scalar2Plus', d)
        self.assertIn('T32MinusMass', d)
        self.assertIn('T32PlusMass', d)
        self.assertIn('Match[t3DMatchInput, EFTOrder -> 5, LoopOrder -> 1]', d)
        self.assertNotIn('RunSequentialT3Matching[', d)
        self.assertNotIn('RunGeneralDMinus2C5RoughBenchmark', d)
        self.assertIn('D_PhysicalScalarSoftOnly_HardInheritanceAndWilsonRGPending', d)

    def test_full_general_builder_and_coupling_inventory(self):
        a = NATIVE.read_text(encoding="utf-8")
        d = PHYSICAL_D.read_text(encoding="utf-8")
        self.assertIn('AllScalarExtensionsRequested', a)
        for name in ('AdditionalQuarticsIncluded', 'AdditionalQuadraticsIncluded',
                     'AdditionalYukawasIncluded', 'BuildGeneralScalarPotentialCandidates'):
            self.assertIn(name, d)
        self.assertIn('t3DCGUnregistered' if 't3DCGUnregistered' in d else
                      't3DUnregistered', d)

    def test_no_single_fermion_fallback(self):
        d = PHYSICAL_D.read_text(encoding="utf-8")
        self.assertNotIn('NF = 1', d)
        self.assertNotIn('NF=1', d)
        self.assertIn('Sum[y1[p,r]*y2[q,r]/MF[r],{r,1,3}]', d)
        self.assertIn('Sum[y1[p,r]*y2Conj[q,r]/MF[r],{r,1,3}]', d)

    def test_no_false_physical_completion(self):
        for source in (NATIVE, PHYSICAL_D):
            s = source.read_text(encoding="utf-8")
            self.assertIn('NOT', s)
        self.assertIn('"full_general_numerical_c5_verified" -> False',
                      NATIVE.read_text(encoding="utf-8"))
        self.assertIn("'complete_numerical_C5_available': False",
                      SCRIPT.read_text(encoding="utf-8"))

    def test_no_restricted_hard_triangle_import(self):
        src = SCRIPT.read_text(encoding="utf-8")
        for forbidden in ("RunT3FHardTriangleBatch1", "RunT3StrictOneLoopBatch2",
                          "RunT3ThesisTwoModelPilot"):
            self.assertNotIn("import " + forbidden, src)

    def test_fit_requires_18(self):
        target = {"central_values": {x: 0. for x in fit.OBS},
                  "one_sigma_errors": {x: 1. for x in fit.OBS}}
        with self.assertRaises(ValueError):
            fit.fit_18_real(np.zeros(17), lambda x: {}, target, max_nfev=1)

    def test_residuals_five(self):
        target = {"central_values": {x: 0. for x in fit.OBS},
                  "one_sigma_errors": {x: 2. for x in fit.OBS}}
        source = {x: 1. for x in fit.OBS}
        np.testing.assert_allclose(fit.residual5(source, target), np.full(5, .5))

    def test_fit_mocked_objective(self):
        target = {"central_values": {x: 0. for x in fit.OBS},
                  "one_sigma_errors": {x: 1. for x in fit.OBS}}
        def fake(x):
            return {k: float(x[i]) for i, k in enumerate(fit.OBS)}
        result = fit.fit_18_real(np.full(18, .01), fake, target, max_nfev=2)
        self.assertEqual(result["status"], "SourceDrivenFitCandidate")
        self.assertEqual(len(result["best"]["parameters"]), 18)

    def test_numeric_evaluator_cannot_ignore_poles(self):
        src = C5.read_text(encoding="utf-8")
        self.assertIn('nF=Lookup[nData,"NF",3]', src)
        self.assertIn('If[nF=!=3', src)
        self.assertIn('{v,1,3}', src)
        self.assertIn('Blocked_UncancelledDimensionalPole', src)
        self.assertIn('nPoleRel<10^-9', src)

    def test_fitter_not_restricted_surrogate(self):
        src = FIT.read_text(encoding="utf-8")
        self.assertNotIn('hard_triangle(', src)
        self.assertNotIn('scalar_bubble(', src)


if __name__ == '__main__':
    unittest.main()
