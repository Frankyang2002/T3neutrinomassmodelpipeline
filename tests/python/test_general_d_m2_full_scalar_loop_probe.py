"""Contract tests for the full-EFT1 hard one-loop processing probe."""
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import types
import unittest

ROOT = Path(__file__).resolve().parents[2]
RUNNER = ROOT / 'scripts/RunGeneralDMinus2FullScalarLoopProbe.py'
WL = ROOT / 'Lagrangian/RunGeneralDMinus2FullScalarLoopProbe.wl'


def load_runner():
    # A test of classification needs no local Matchete or source EFT1 WXF.
    mod_name = 'scripts.RunGeneralDMinus2ScalarThreshold'
    stub = types.ModuleType(mod_name)
    stub.preflight = lambda *a, **k: {'insertion_file': '/dummy/running.wl'}
    old = sys.modules.get(mod_name)
    sys.modules[mod_name] = stub
    try:
        spec = importlib.util.spec_from_file_location('D2LoopProbeRunner', RUNNER)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    finally:
        if old is None:
            sys.modules.pop(mod_name, None)
        else:
            sys.modules[mod_name] = old


class FullScalarLoopProbeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.runner = load_runner()
        cls.source = WL.read_text(encoding='utf-8')

    def good_audit(self):
        return {
            'Status': 'Passed', 'SMAlphaComparisonCertified': True,
            'MassPowerCountingPromotionOnly': True,
            'AllPhysicalScalarsRegisteredHeavy': True,
            'PureHeavyQuadraticAlphaEqualToFreeLag': True,
            'CheckLagrangianPassed': True, 'FullPhysicalSMIncluded': True,
            'NoLinearHeavyScalarSources': True,
        }

    def good_result(self):
        return {
            'Status': 'FullEFT1OneLoopPostprocessed',
            'OneLoopRawMatchReturned': True,
            'FullPhysicalSMIncluded': True, 'LoopOrder': 1,
            'EFTOrder': 5, 'RemainingHeavyFields': [],
            'MatchResultContainsHbar': True,
            'ProcessedEFTHasHbar': True,
            'HardHbarSectorExtracted': True,
            'SimplificationPerformed': True,
            'ProcessingStages': [
                {'Stage': n, 'Status': 'Returned'}
                for n in ('GreensSimplify', 'EOMSimplify',
                          'EvaluateLoopFunctions', 'ReplaceEffectiveCouplings')
            ],
        }

    def test_success_contract(self):
        self.assertEqual(self.runner.classify(code=0, timed_out=False,
            audit=self.good_audit(), result=self.good_result()),
            'FullEFT1OneLoopPostprocessed')

    def test_timeout(self):
        self.assertEqual(self.runner.classify(code=0, timed_out=True,
            audit=self.good_audit(), result=self.good_result()),
            'ExternalTimeout')

    def test_nonzero_exit(self):
        self.assertEqual(self.runner.classify(code=26, timed_out=False,
            audit=self.good_audit(), result=self.good_result()), 'WolframFailed')

    def test_missing_mass_validation(self):
        bad = self.good_audit()
        bad['PureHeavyQuadraticAlphaEqualToFreeLag'] = False
        self.assertEqual(self.runner.classify(code=0, timed_out=False,
            audit=bad, result=self.good_result()), 'MissingMassOrFullEFT1Validation')

    def test_no_full_sm_rejected(self):
        bad = self.good_audit()
        bad['FullPhysicalSMIncluded'] = False
        self.assertNotEqual(self.runner.classify(code=0, timed_out=False,
            audit=bad, result=self.good_result()), 'FullEFT1OneLoopPostprocessed')

    def test_no_hbar_rejected(self):
        bad = self.good_result()
        bad['MatchResultContainsHbar'] = False
        self.assertNotEqual(self.runner.classify(code=0, timed_out=False,
            audit=self.good_audit(), result=bad), 'FullEFT1OneLoopPostprocessed')

    def test_stale_result_rejected(self):
        self.assertEqual(self.runner.classify(code=0, timed_out=False,
            audit=self.good_audit(), result=None), 'InvalidOneLoopMatchResult')

    def test_full_source_is_used(self):
        self.assertIn('matchInput = D2PromoteScalarMasses[exactTree];', self.source)
        self.assertIn('Match[matchInput, EFTOrder -> 5, LoopOrder -> 1]', self.source)
        self.assertNotIn('matchInput = RelabelIndices[coreTree]', self.source)

    def test_mass_promotion_strictly_only_three_names(self):
        self.assertIn('D2PromoteScalarMasses[expr_]', self.source)
        self.assertIn('D2DemoteScalarMasses[expr_]', self.source)
        for mass in ('MS1', 'T3S2RMass', 'T3S2IMass'):
            self.assertIn(f'Matchete`Coupling[{mass}, {{}}, 1]', self.source)
            self.assertIn(f'Matchete`Coupling[{mass}, {{}}, 0]', self.source)
        self.assertIn('ExactReversibleMassPromotion', self.source)
        self.assertIn('SourceLightMassCounts', self.source)
        self.assertIn('StaleLightMassCounts', self.source)
        self.assertIn('OriginalQuadraticAlphaEqualToFreeLag', self.source)
        self.assertIn('With[{m = massSymbols[[i]]}', self.source)
        self.assertIn('Count[exactTree,', self.source)
        self.assertIn('Count[matchInput,', self.source)
        self.assertNotIn('Count[HoldComplete[exactTree]', self.source)

    def test_mass_promotion_required_by_classification(self):
        bad = self.good_audit()
        bad['MassPowerCountingPromotionOnly'] = False
        self.assertEqual(self.runner.classify(code=0, timed_out=False,
            audit=bad, result=self.good_result()), 'MissingMassOrFullEFT1Validation')

    def test_mass_guard_precedes_match(self):
        self.assertIn('FreeLag[NewScalar1, T3S2R, T3S2I]', self.source)
        self.assertIn('PureHeavyQuadraticAlphaEqualToFreeLag', self.source)
        self.assertLess(self.source.index('If[!massSafe'), self.source.index('Match[matchInput'))
        self.assertLess(self.source.index('CheckLagrangian[matchInput]'),
                        self.source.index('Match[matchInput'))

    def test_stage_separation(self):
        self.assertIn('"InheritedOneLoopAdded" -> False', self.source)
        self.assertIn('"DirectRunningAdded" -> False', self.source)
        self.assertNotIn('Get[insertionPath]', self.source)
        self.assertIn('"FullGeneralScalarThresholdCertified" -> False', self.source)

    def test_bounded_external_run(self):
        self.assertIn('timeout_seconds * 7 + startup_allowance_seconds', RUNNER.read_text())
        self.assertIn('taskkill', RUNNER.read_text())

    def test_missing_postprocessing_rejected(self):
        bad = self.good_result()
        bad['ProcessingStages'][2]['Status'] = 'Aborted'
        self.assertEqual(self.runner.classify(code=0, timed_out=False,
            audit=self.good_audit(), result=bad), 'InvalidOneLoopMatchResult')

    def test_no_candidate_is_not_falsely_certified(self):
        self.assertIn('"FullFlavorHardC5Verified" -> False', self.source)
        self.assertIn('"FixedOrderSMEFTAssembled" -> False', self.source)
        self.assertIn('"InheritedOneLoopAdded" -> False', self.source)
        self.assertIn('"DirectRunningAdded" -> False', self.source)
        self.assertIn('ExtractWeinbergCoefficient[loopPiece]', self.source)

    def test_matching_pipeline_order(self):
        match = self.source.index('Match[matchInput, EFTOrder -> 5, LoopOrder -> 1]')
        green = self.source.index('processStages = {')
        self.assertLess(match, green)
        for stage in ('GreensSimplify', 'EOMSimplify',
                      'EvaluateLoopFunctions', 'ReplaceEffectiveCouplings'):
            self.assertIn(stage, self.source)

    def test_raw_json_object(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / 'x.json'
            path.write_text(json.dumps({'Status': 'Passed'}))
            self.assertEqual(self.runner._json(path), {'Status': 'Passed'})


if __name__ == '__main__':
    unittest.main()
