"""Contract tests for the integrated local general D(-2) sequential runner.

No Mathematica installation or live matching required. Uses small synthetic
reports solely to test control flow and fail-closed validation.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

SCRIPT = Path(__file__).resolve().parents[2] / "scripts/RunGeneralDMinus2SequentialMatching.py"
SPEC = importlib.util.spec_from_file_location("general_d2_seq", SCRIPT)
assert SPEC and SPEC.loader
m = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(m)


def reports():
    return {
        "scalar_hard_one_loop": {"Status": "FullEFT1OneLoopPostprocessed",
            "OneLoopRawMatchReturned": True, "MatchResultContainsHbar": True,
            "ProcessedEFTHasHbar": True, "RemainingHeavyFields": [],
            "WeinbergExtractionStatus": "Success", "WeinbergTermCount": 68,
            "ProcessingStages": [{"Status": "Returned"}] * 4},
        "components": {"Status": "ComponentsAudited", "HardC5LinearInHbar": True,
            "DirectRunningHasExplicitHbar": False, "HardWeinbergTermCount": 68,
            "InheritedPureSMWeinbergTermCount": 12},
        "operator_structure": {"Status": "WeinbergStructureAudited",
            "AllCandidateTermsHaveTwoLeptonsTwoHiggs": True,
            "AllCandidateTermsHaveUnambiguousChirality": True,
            "HardC5IsLinearInHbar": True, "InheritedC5IsLinearInHbar": True},
        "running_symmetry": {"Status": "RunningSymmetryAudited",
            "SymmetricModuloSummationDummies": True,
            "RunningFlavourSymmetryCertifiedModuloDummyRenaming": True,
            "AlphaNormalizedResidualTermCount": 0},
        "candidate_assembly": {"Status": "C5FlavourCandidateExported",
            "FixedOrderC5CandidateAssembled": True,
            "EachHolomorphicTermHasTwoDifferentLeptonFlavourLabels": True,
            "HardSymmetricByConstruction": True, "InheritedSymmetricByConstruction": True,
            "RunningSymmetricModuloSummationDummies": True,
            "TotalCandidateLinearInHbar": True},
        "warsaw_reference": {"Status": "WarsawNormalizationReferenceAudited",
            "ProjectExtractorRecognizesWarsawWeinbergOperator": True,
            "ProjectExtractorLinearOnWarsawReference": True,
            "WarsawOrderedCoefficientWithReferenceCouplingsSetToOne": "1",
            "WarsawOrderedCoefficientUnitWilsonFactorNumeric": True,
            "WarsawSymmetricUnitWilsonFactorAssumingEqualPQandQP": "2",
            "WarsawOperatorNeutrinoVEVComponentIsPlusVSquaredOverTwo": True},
        "benchmark": {"Status": "ProvisionalFinitePartBenchmarkEvaluated",
            "NeutrinoMasses_eV": [0.0, 0.00056384, 0.001938],
            "MassSquaredDifferences_eV2": [3.2e-7, 3.75e-6],
            "FinitePartOnly": True, "RenormalizationCertified": False,
            "C5RelativeSymmetryResidual": 0},
        "exact_pole": {"Status": "ExactNF1PoleCancellationVerified", "NF": 1,
            "ExactPoleCancellationVerifiedForNF1": True,
            "ExactZeroByTogether": True, "ReducedResidueTermCount": 0},
    }


class ReportContractTests(unittest.TestCase):
    def test_every_nominal_report_passes(self):
        for stage, report in reports().items():
            with self.subTest(stage=stage):
                m.validate_report(stage, report)

    def test_warsaw_exporter_actual_json_string_factors(self):
        d = reports()["warsaw_reference"]
        self.assertIsNone(m.validate_report("warsaw_reference", d))
        numeric = {**d,
                   "WarsawOrderedCoefficientWithReferenceCouplingsSetToOne": 1,
                   "WarsawSymmetricUnitWilsonFactorAssumingEqualPQandQP": 2}
        self.assertIsNone(m.validate_report("warsaw_reference", numeric))

    def test_warsaw_report_rejects_wrong_factors_or_false_numeric_flag(self):
        for key, value in (
            ("WarsawOrderedCoefficientWithReferenceCouplingsSetToOne", "0"),
            ("WarsawOrderedCoefficientWithReferenceCouplingsSetToOne", "1.0"),
            ("WarsawOrderedCoefficientWithReferenceCouplingsSetToOne", True),
            ("WarsawSymmetricUnitWilsonFactorAssumingEqualPQandQP", "1"),
            ("WarsawSymmetricUnitWilsonFactorAssumingEqualPQandQP", "2.0"),
            ("WarsawOrderedCoefficientUnitWilsonFactorNumeric", False),
        ):
            with self.subTest(field=key, value=value):
                d = {**reports()["warsaw_reference"], key: value}
                with self.assertRaises(m.PipelineFailure):
                    m.validate_report("warsaw_reference", d)

    def test_bad_status_rejected(self):
        for stage, report in reports().items():
            with self.subTest(stage=stage):
                report = {**report, "Status": "Failed"}
                with self.assertRaises(m.PipelineFailure):
                    m.validate_report(stage, report)

    def test_pole_must_cancel_exactly(self):
        d = reports()["exact_pole"]
        d["ReducedResidueTermCount"] = 1
        with self.assertRaises(m.PipelineFailure):
            m.validate_report("exact_pole", d)

    def test_running_symmetry_modulo_bound_indices(self):
        d = reports()["running_symmetry"]
        d["SymmetricModuloSummationDummies"] = False
        with self.assertRaises(m.PipelineFailure):
            m.validate_report("running_symmetry", d)

    def test_hard_match_must_complete_all_processing_stages(self):
        d = reports()["scalar_hard_one_loop"]
        d["ProcessingStages"][2] = {"Status": "Aborted"}
        with self.assertRaises(m.PipelineFailure):
            m.validate_report("scalar_hard_one_loop", d)

    def test_complex_symmetric_matrix_requires_finite_takagi_masses(self):
        d = reports()["benchmark"]
        d["NeutrinoMasses_eV"] = [0, float("nan"), 0.001]
        with self.assertRaises(m.PipelineFailure):
            m.validate_report("benchmark", d)

    def test_mass_audit_rejects_stale_light_mass_counting(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "profile_full_scalar_one_loop_audit.json"
            d = {k: True for k in (
                "AllPhysicalScalarsRegisteredHeavy", "MassPowerCountingPromotionOnly",
                "PureHeavyQuadraticAlphaEqualToFreeLag", "ExactReversibleMassPromotion",
                "CheckLagrangianPassed", "SMAlphaComparisonCertified",
                "NoLinearHeavyScalarSources", "FullPhysicalSMIncluded")}
            d["Status"] = "Passed"
            p.write_text(json.dumps(d))
            m.validate_mass_audit(Path(tmp))
            d["MassPowerCountingPromotionOnly"] = False
            p.write_text(json.dumps(d))
            with self.assertRaises(m.PipelineFailure):
                m.validate_mass_audit(Path(tmp))

    def test_missing_json_fails_closed(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(m.PipelineFailure):
                m.read_object(Path(d) / "absent.json")

    def test_invalid_json_fails_closed(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "bad.json"
            p.write_text('')
            with self.assertRaises(m.PipelineFailure):
                m.read_object(p)

    def test_boundary_hashes_change_on_modification(self):
        with tempfile.TemporaryDirectory() as d:
            out = Path(d)
            data = out / "data"
            data.mkdir()
            for name in m.BOUNDARY_FILES:
                (data / name).write_bytes(name.encode())
            original = m.boundary_hashes(out)
            self.assertEqual(len(original), 4)
            (data / m.BOUNDARY_FILES[2]).write_bytes(b"modified")
            self.assertNotEqual(original, m.boundary_hashes(out))

    def test_no_legacy_scalar_driver_overwrite(self):
        contents = SCRIPT.read_text(encoding="utf-8")
        self.assertIn("RunGeneralDMinus2SequentialScalarThreshold.wl", contents)
        self.assertNotIn('"Lagrangian/RunGeneralDMinus2ScalarThreshold.wl",', contents)

    def test_order_hard_before_candidate_before_pole(self):
        source = SCRIPT.read_text(encoding="utf-8")
        self.assertLess(source.index('wolfram("scalar_hard_one_loop"'),
                        source.index('wolfram("candidate_assembly"'))
        self.assertLess(source.index('wolfram("candidate_assembly"'),
                        source.index('wolfram("exact_pole"'))

    def test_strict_mass_order_promotion_in_consolidated_scalar_driver(self):
        wl = SCRIPT.parents[1] / "Lagrangian/RunGeneralDMinus2SequentialScalarThreshold.wl"
        t = wl.read_text(encoding="utf-8")
        for symbol in ("MS1", "T3S2RMass", "T3S2IMass"):
            self.assertIn(f"Matchete`Coupling[{symbol}, {{}}, 1]", t)
            self.assertIn(f"Matchete`Coupling[{symbol}, {{}}, 0]", t)
        self.assertIn('Match[matchInput, EFTOrder -> 5, LoopOrder -> 1]', t)
        self.assertIn('PureHeavyQuadraticAlphaEqualToFreeLag', t)
        self.assertIn('NoLinearHeavyScalarSources', t)
        self.assertNotIn('RemoveField[y2Conj]', t)

    def test_scalar_launcher_does_not_overwrite_script_name_with_artifact_name(self):
        # Regression for an actual run: cleanup names must not shadow the
        # "RunGeneralDMinus2SequentialScalarThreshold.wl" launch target.
        source = SCRIPT.read_text(encoding="utf-8")
        self.assertIn("for artifact_filename in (", source)
        self.assertIn("/ script_filename", source)
        self.assertNotIn("ROOT / \"Lagrangian\" / filename", source)

    def test_run_command_rejects_nonexistent_wolfram_source_before_launch(self):
        with tempfile.TemporaryDirectory() as td:
            log_dir = Path(td)
            for bad in (log_dir / "missing.wl", log_dir / "missing.txt"):
                with self.subTest(path=bad), self.assertRaises(m.PipelineFailure):
                    m.run_command("probe", ["wolframscript", "-file", str(bad)],
                                  cap_seconds=30, log_dir=log_dir)

    def test_failed_run_exports_failure_manifest(self):
        with tempfile.TemporaryDirectory() as td:
            out = Path(td) / "out"
            conf = Path(td) / "benchmark.json"
            conf.write_text("{}")
            with patch.object(m, "validate_sources"), \
                 patch.object(m, "run_command", side_effect=m.PipelineFailure("mock process crash")):
                with self.assertRaises(m.PipelineFailure):
                    m.run_pipeline(out, conf, reuse_boundary=False,
                        wolframscript="wolframscript", stage1_timeout=100,
                        python_stage_timeout=100, scalar_stage_timeout=90,
                        wolfram_timeout=180, startup_allowance=120)
            manifest = m.read_object(out / m.REL_RUN / "sequential_matching_manifest.json")
            self.assertEqual(manifest["status"], "Failed")
            self.assertEqual(manifest["failed_step"], "fermion_threshold")
            self.assertIn("mock process crash", manifest["error"])

    def test_resume_after_warsaw_does_not_repeat_scalar_matching(self):
        with tempfile.TemporaryDirectory() as td:
            out = Path(td) / "out"
            profile = out / m.REL_PROFILE
            log_dir = out / m.REL_RUN
            profile.mkdir(parents=True)
            log_dir.mkdir(parents=True)
            for name in m.BOUNDARY_FILES:
                (out / "data" / name).write_bytes(name.encode())
            config = Path(td) / "benchmark.json"
            config.write_text("{}")
            all_reports = reports()
            stages = ("scalar_hard_one_loop", "components", "operator_structure",
                      "running_symmetry", "candidate_assembly", "warsaw_reference")
            prior_steps = []
            for stage in stages:
                (profile / m.STAGE_REPORTS[stage]).write_text(json.dumps(all_reports[stage]))
                rec = {"stage": stage, "exit_code": 0,
                       "timed_out": False, "elapsed_seconds": 0.5, "log": "fake"}
                if stage != "warsaw_reference":
                    rec["reported_status"] = m.EXPECTED_STATUSES[stage]
                prior_steps.append(rec)
            mass = {k: True for k in (
                "AllPhysicalScalarsRegisteredHeavy", "MassPowerCountingPromotionOnly",
                "PureHeavyQuadraticAlphaEqualToFreeLag", "ExactReversibleMassPromotion",
                "CheckLagrangianPassed", "SMAlphaComparisonCertified",
                "NoLinearHeavyScalarSources", "FullPhysicalSMIncluded")}
            mass["Status"] = "Passed"
            (profile / "profile_full_scalar_one_loop_audit.json").write_text(json.dumps(mass))
            for name in ("profile_full_scalar_one_loop_processed.wxf",
                         "profile_full_scalar_one_loop_hard_hbar_sector.wxf",
                         "profile_c5_flavour_candidate_total_provisional.wxf"):
                (profile / name).write_bytes(b"test")
            prior = {"status": "Failed", "failed_step": "warsaw_reference",
                     "mode": "ReuseVerifiedEFT1", "model": "D_-2",
                     "generated_utc": "2026-11-07T00:00:00Z",
                     "benchmark_config": str(config.resolve()),
                     "eft1_boundary_sha256": m.boundary_hashes(out),
                     "steps": prior_steps}
            manifest_path = log_dir / "sequential_matching_manifest.json"
            manifest_path.write_text(json.dumps(prior))
            commands = []

            def fake_run(stage, cmd, **kwargs):
                commands.append(stage)
                (profile / m.STAGE_REPORTS[stage]).write_text(json.dumps(all_reports[stage]))
                return {"stage": stage, "exit_code": 0,
                        "timed_out": False, "log": "fake", "elapsed_seconds": 0.1}

            with patch.object(m, "validate_sources"), \
                 patch.object(m, "initial_preflight", return_value=config), \
                 patch.object(m, "run_command", side_effect=fake_run):
                result = m.run_pipeline(
                    out, config, reuse_boundary=True,
                    wolframscript="wolframscript", stage1_timeout=120,
                    python_stage_timeout=120, scalar_stage_timeout=90,
                    wolfram_timeout=120, startup_allowance=120,
                    resume_from_benchmark=True)
            self.assertEqual(commands, ["benchmark", "exact_pole"])
            self.assertEqual(result["status"],
                             "SequentialC5BenchmarkReproduced_ProvisionalNormalization")
            self.assertEqual(result["mode"], "ResumeValidatedWarsawReference")
            self.assertTrue(result["resumed_from_previous_failed_run"])
            self.assertFalse(result["scalar_matching_executed_this_run"])
            self.assertEqual(len(result["steps"]), 8)
            self.assertTrue(result["boundary_hashes_unchanged"])

            # Invalidated boundary prevents resume BEFORE overwriting manifest.
            manifest_path.write_text(json.dumps(prior))
            (out / "data" / m.BOUNDARY_FILES[0]).write_bytes(b"changed")
            with patch.object(m, "validate_sources"):
                with self.assertRaises(m.PipelineFailure):
                    m.run_pipeline(
                        out, config, reuse_boundary=True,
                        wolframscript="wolframscript", stage1_timeout=120,
                        python_stage_timeout=120, scalar_stage_timeout=90,
                        wolfram_timeout=120, startup_allowance=120,
                        resume_from_benchmark=True)
            self.assertEqual(json.loads(manifest_path.read_text()), prior)

    def test_integrated_schematic_test_from_uv_and_reuse(self):
        # Stub only OS execution; exercise actual manifest, hashes and checks.
        for reuse in (True, False):
            with self.subTest(reuse=reuse), tempfile.TemporaryDirectory() as td:
                out = Path(td) / "out"
                profile = out / m.REL_PROFILE
                profile.mkdir(parents=True)
                for name in m.BOUNDARY_FILES:
                    (out / "data" / name).write_bytes(name.encode())
                conf = Path(td) / "benchmark.json"
                conf.write_text("{}")
                expected = reports()
                commands: list[str] = []
                expected_scripts = {
                    "scalar_hard_one_loop": "RunGeneralDMinus2SequentialScalarThreshold.wl",
                    "components": "RunGeneralDMinus2FixedOrderComponents.wl",
                    "operator_structure": "RunGeneralDMinus2C5FlavourStructureAudit.wl",
                    "running_symmetry": "RunGeneralDMinus2RunningSymmetryAudit.wl",
                    "candidate_assembly": "RunGeneralDMinus2C5FlavourAssemblyCandidate.wl",
                    "warsaw_reference": "RunGeneralDMinus2C5NormalizationReference.wl",
                    "benchmark": "RunGeneralDMinus2C5RoughBenchmark.wl",
                    "exact_pole": "RunGeneralDMinus2C5PoleCancellationAudit.wl",
                }
                def fake_run(stage, cmd, **kwargs):
                    commands.append(stage)
                    if stage in expected_scripts:
                        self.assertEqual(cmd[1], "-file")
                        self.assertEqual(Path(cmd[2]).name, expected_scripts[stage])
                        self.assertTrue(Path(cmd[2]).is_file())
                    if stage in expected:
                        (profile / m.STAGE_REPORTS[stage]).write_text(json.dumps(expected[stage]))
                        if stage == "scalar_hard_one_loop":
                            mass = {k: True for k in (
                                "AllPhysicalScalarsRegisteredHeavy", "MassPowerCountingPromotionOnly",
                                "PureHeavyQuadraticAlphaEqualToFreeLag", "ExactReversibleMassPromotion",
                                "CheckLagrangianPassed", "SMAlphaComparisonCertified",
                                "NoLinearHeavyScalarSources", "FullPhysicalSMIncluded")}
                            mass["Status"] = "Passed"
                            (profile / "profile_full_scalar_one_loop_audit.json").write_text(json.dumps(mass))
                            (profile / "profile_full_scalar_one_loop_processed.wxf").write_bytes(b"p")
                            (profile / "profile_full_scalar_one_loop_hard_hbar_sector.wxf").write_bytes(b"h")
                        elif stage == "candidate_assembly":
                            (profile / "profile_c5_flavour_candidate_total_provisional.wxf").write_bytes(b"c")
                    return {"stage": stage, "exit_code": 0, "elapsed_seconds": 0.01, "log": "fake"}
                with patch.object(m, "validate_sources"), \
                     patch.object(m, "initial_preflight", return_value=conf), \
                     patch.object(m, "run_command", side_effect=fake_run):
                    result = m.run_pipeline(out, conf, reuse_boundary=reuse,
                        wolframscript="wolframscript", stage1_timeout=100,
                        python_stage_timeout=100, scalar_stage_timeout=90,
                        wolfram_timeout=180, startup_allowance=120)
                self.assertEqual(result["status"], "SequentialC5BenchmarkReproduced_ProvisionalNormalization")
                self.assertTrue(result["boundary_hashes_unchanged"])
                self.assertEqual(commands[-1], "exact_pole")
                self.assertEqual(commands[0], "scalar_hard_one_loop" if reuse else "fermion_threshold")
                self.assertEqual(len(commands), 8 if reuse else 11)
                self.assertTrue((out / m.REL_RUN / "sequential_matching_manifest.json").is_file())


if __name__ == "__main__":
    unittest.main()
