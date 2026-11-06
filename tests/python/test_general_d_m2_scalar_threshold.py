"""Tests for the D(-2) real-scalar second-threshold *contract*.

No local Mathematica kernel is assumed.  The live Matchete calculation is
validated on the user's Wolfram installation, never claimed by these tests.
"""
from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from scripts.RunGeneralDMinus2ScalarThreshold import (
    WXF_FILES, make_command, make_profile_command, preflight, run,
    run_tree_profile, PROFILE_PROBES, PROFILE_SMOKE_PROBES, PROFILE_BSM_PROBES,
)


class TestD2ScalarThreshold(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.matching = Path(self.temporary.name)
        self.data = self.matching / "data"
        self.data.mkdir()
        self.wl = self.matching / "test.wl"
        self.wl.write_text("ClearAll[x];", encoding="utf-8")
        self._json(self.data / "general_d_m2_stage1_summary.json", {
            "status": "GeneralDMinus2RealMassBasisFermionThresholdVerified",
            "general_uv": {"extra_quartics": 5, "neutral_scalar_bilinears": 1,
                           "extra_yukawas": 1},
            "first_matching": {"tree_wilson_terms": 8,
                "remaining_physical_BSM_scalars": ["NewScalar1", "T3S2R", "T3S2I"]},
        })
        self._json(self.data / "general_d_m2_eft1_rge/pilot_summary.json", {
            "status": "GeneralDMinus2RealBasisWeinbergBetaVerifiedOneGeneration"})
        self._json(self.data / "general_d_m2_full_flavor/pilot_summary.json", {
            "status": "GeneralDMinus2FullFlavorDirectRunningPrepared"})
        self._json(self.data / "general_d_m2_full_flavor/full_flavor_direct_running.json", {
            "status": "Success", "one_generation_reduction_matches": True,
            "symmetric_under_external_flavor_exchange": True,
            "equal_scale_running_vanishes": True,
            "real_scalar_hard_threshold_performed": False,
        })
        self._json(self.data / "general_d_m2_full_flavor/scalar_threshold_preflight.json", {
            "status": "FullFlavorDirectInsertionPrepared_ScalarHardThresholdNotRun",
            "scalar_matching_integrate_out": ["NewScalar1", "T3S2R", "T3S2I"],
            "scalar_masses_matchete": ["MS1", "T3S2RMass", "T3S2IMass"],
            "physical_mass_relations_imposed": False,
            "insertion_is_pure_SM_and_is_not_to_be_one_loop_matched": True,
            "insertion_coefficient_omits_hbar": True,
        })
        self._json(self.matching / "scalar_mass_basis.json", {
            "RealFields": ["T3S2R", "T3S2I"],
            "RealMassSquared": "MS2^2-2*Abs[kappa]",
            "ImagMassSquared": "MS2^2+2*Abs[kappa]",
        })
        for filename in WXF_FILES:
            (self.data / filename).write_bytes(b"wxf-placeholder")
        for filename in ("eft1_after_F_realbasis_wilson_tensor_seed.json",
                         "eft1_after_F_realbasis_quartic_tensor_seed.json"):
            self._json(self.data / filename, {"Status": "Success"})
        insert = self.data / "general_d_m2_full_flavor/full_flavor_weinberg_insertion.wl"
        insert.write_text("EFT1D2WilsonRunningHeavyInsertion = 0;\n"
                          "EFT1D2WilsonRunningWeinbergCoefficient[p_,q_] := "
                          "EFT1D2C12[p,q]+EFT1D2C1Conj[p,q];\n",
                          encoding="utf-8")

    @staticmethod
    def _json(path: Path, obj: dict):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(obj), encoding="utf-8")

    def test_valid_boundary_passes(self):
        self.assertEqual(preflight(self.matching, self.wl)["status"],
            "GeneralDMinus2ScalarThresholdPreflightPassed")

    def test_preflight_only_never_runs_wolfram(self):
        result = run(self.matching, preflight_only=True, script=self.wl)
        self.assertEqual(result["status"],
            "GeneralDMinus2ScalarThresholdPreflightPassed")
        self.assertFalse((self.data / "general_d_m2_scalar_threshold/wolfram_scalar_threshold.log").exists())

    def test_missing_transition_full_rejected(self):
        (self.data / WXF_FILES[-1]).unlink()
        with self.assertRaisesRegex(FileNotFoundError, "transition_full"):
            preflight(self.matching, self.wl)

    def test_missing_extra_quartic_rejected(self):
        path = self.data / "general_d_m2_stage1_summary.json"
        p = json.loads(path.read_text()); p["general_uv"]["extra_quartics"] = 4
        self._json(path, p)
        with self.assertRaisesRegex(ValueError, "full general"):
            preflight(self.matching, self.wl)

    def test_restricted_boundary_rejected(self):
        path = self.data / "general_d_m2_stage1_summary.json"
        p = json.loads(path.read_text()); p["status"] = "RestrictedProduction"
        self._json(path, p)
        with self.assertRaisesRegex(ValueError, "general real-basis"):
            preflight(self.matching, self.wl)

    def test_direct_running_must_not_be_one_loop_matched(self):
        path = self.data / "general_d_m2_full_flavor/scalar_threshold_preflight.json"
        p = json.loads(path.read_text())
        p["insertion_is_pure_SM_and_is_not_to_be_one_loop_matched"] = False
        self._json(path, p)
        with self.assertRaisesRegex(ValueError, "power convention"):
            preflight(self.matching, self.wl)

    def test_three_distinct_mass_symbols(self):
        self.assertEqual(preflight(self.matching, self.wl)["physical_scalar_fields"],
            ["NewScalar1", "T3S2R", "T3S2I"])

    def test_wrong_mass_sign_fails(self):
        path = self.matching / "scalar_mass_basis.json"
        p = json.loads(path.read_text())
        p["RealMassSquared"] = "MS2^2+2*Abs[kappa]"
        self._json(path, p)
        with self.assertRaisesRegex(ValueError, "mass convention"):
            preflight(self.matching, self.wl)

    def test_missing_insertion_fails(self):
        (self.data / "general_d_m2_full_flavor/full_flavor_weinberg_insertion.wl").unlink()
        with self.assertRaisesRegex(FileNotFoundError, "insertion missing"):
            preflight(self.matching, self.wl)

    def test_missing_new_driver_fails(self):
        with self.assertRaisesRegex(FileNotFoundError, "Wolfram scalar-threshold driver"):
            preflight(self.matching, self.matching / "absent.wl")

    def test_running_insertion_cannot_omit_heavy_zero(self):
        path = self.data / "general_d_m2_full_flavor/full_flavor_weinberg_insertion.wl"
        path.write_text("EFT1D2WilsonRunningWeinbergCoefficient=1; EFT1D2C12=1; EFT1D2C1Conj=1;")
        with self.assertRaisesRegex(ValueError, "HeavyInsertion"):
            preflight(self.matching, self.wl)

    def test_correct_command_target(self):
        cmd = make_command(self.matching, self.wl, "fake-math", self.wl)
        self.assertEqual(cmd[0:2], ["fake-math", "-file"])
        self.assertEqual(cmd[-2], str(self.matching.resolve()))

    def test_failed_preflight_emits_summary_not_success(self):
        (self.data / WXF_FILES[0]).unlink()
        state = run(self.matching, preflight_only=True, script=self.wl)
        self.assertEqual(state["status"], "Failed")
        summary = json.loads((self.data / "general_d_m2_scalar_threshold/pilot_summary.json").read_text())
        self.assertEqual(summary["status"], "Failed")


    def test_profile_command_targets_same_wolfram_script_with_mode(self):
        cmd = make_profile_command(
            self.matching, self.wl, "with_y2Conj", 75,
            wolframscript="fake-math", script=self.wl,
        )
        self.assertEqual(cmd[:2], ["fake-math", "-file"])
        self.assertEqual(cmd[-2:], ["profile:with_y2Conj", "75"])
        self.assertEqual(cmd[-4], str(self.matching.resolve()))

    def test_profile_rejects_invalid_probe_and_timeout(self):
        with self.assertRaisesRegex(ValueError, "Unknown"):
            make_profile_command(self.matching, self.wl, "restricted", script=self.wl)
        with self.assertRaisesRegex(ValueError, "probe_timeout"):
            make_profile_command(self.matching, self.wl, "inspect", 0, script=self.wl)

    def test_tree_profile_does_not_claim_scalar_matching_success(self):
        def fake_probe(matching_dir, insertion, probe, *args):
            return {"probe": probe,
                    "status": "Inspected" if probe == "inspect" else "TimedOut"}
        with patch("scripts.RunGeneralDMinus2ScalarThreshold._execute_profile_probe",
                   side_effect=fake_probe) as mocked:
            result = run_tree_profile(self.matching, script=self.wl)
        self.assertEqual(result["status"], "TreeProfileRecorded")
        self.assertFalse(result["full_general_threshold_matched"])
        self.assertEqual(mocked.call_count, 2)
        self.assertEqual([p["probe"] for p in result["probes"]], ["inspect", "free"])
        self.assertEqual(result["baseline_bottleneck_probe"], "free")
        self.assertEqual(result["baseline_bottleneck_status"], "TimedOut")
        self.assertEqual(result["skipped_probes"], list(PROFILE_PROBES)[2:])
        self.assertFalse((self.data / "general_d_m2_scalar_threshold/pilot_summary.json").exists())
        self.assertTrue((self.data / "general_d_m2_scalar_threshold/profile_tree/profile_tree_summary.json").exists())

    def test_failed_tree_inspection_blocks_other_probe_modes(self):
        with patch("scripts.RunGeneralDMinus2ScalarThreshold._execute_profile_probe",
                   return_value={"probe": "inspect", "status": "WolframFailed"}) as mocked:
            result = run_tree_profile(self.matching, script=self.wl)
        self.assertEqual(result["status"], "ProfileInspectionFailed")
        self.assertEqual(mocked.call_count, 1)

    def test_full_tree_probe_is_explicit_opt_in(self):
        def fake_probe(matching_dir, insertion, probe, *args):
            return {"probe": probe,
                    "status": "Inspected" if probe == "inspect" else "MatchReturned"}
        with patch("scripts.RunGeneralDMinus2ScalarThreshold._execute_profile_probe",
                   side_effect=fake_probe):
            result = run_tree_profile(self.matching, script=self.wl,
                                      include_full_tree=True)
        self.assertEqual(result["probes"][-1]["probe"], "full-tree")

    def test_profile_fails_fast_on_invalid_timeouts(self):
        result = run_tree_profile(self.matching, script=self.wl,
                                  probe_timeout_seconds=0)
        self.assertEqual(result["status"], "Failed")
        self.assertIn("timeout", result["error"])

    def test_wolfram_profile_uses_original_eft1_and_no_one_loop_probe(self):
        actual = (Path(__file__).resolve().parents[2] /
                  "Lagrangian/RunGeneralDMinus2ScalarThreshold.wl")
        code = actual.read_text(encoding="utf-8")
        self.assertIn("profileTree = Expand[currentTree]", code)
        self.assertIn("ReplaceEffectiveCouplings[transitionTree]", code)
        self.assertIn('"ExplicitMatchesExpandedCanonical" -> profileExactRepresentation', code)
        self.assertIn('"full-tree", transitionTree', code)
        self.assertIn("profileTerms = If[Head[profileTree] === Plus", code)
        self.assertIn("TimeConstrained[", code)
        self.assertIn("Match[profileInput, EFTOrder -> 5, LoopOrder -> 0]", code)
        self.assertIn("profile:full-tree", code)
        self.assertLess(code.index('If[profileMode,\n  Print["[profile] Preparing'),
                        code.index('treeResult = UsingFrontEnd[RunT3Matching[transitionTree, 5, 0]]'))

    def test_profile_does_not_split_canonical_eft1_as_one_giant_term(self):
        actual = (Path(__file__).resolve().parents[2] /
                  "Lagrangian/RunGeneralDMinus2ScalarThreshold.wl")
        code = actual.read_text(encoding="utf-8")
        self.assertIn("profileCompactTree = Expand[transitionTree]", code)
        self.assertIn("profileTree = Expand[currentTree]", code)
        self.assertNotIn("profileTree = Expand[transitionTree];", code)
        self.assertIn("profileExactRepresentation", code)
        self.assertLess(code.index('If[!profileExactRepresentation || !profileExact ||'),
                        code.index('If[Length[profileQuarticPositions] === 0'))

    def test_wolfram_script_separates_ABCs(self):
        # Static guard: B takes transitionTree, C takes transitionFullCanonical;
        # directC5 is exported but must not appear in RunT3Matching arguments.
        actual = (Path(__file__).resolve().parents[2] /
                  "Lagrangian/RunGeneralDMinus2ScalarThreshold.wl")
        if not actual.is_file():
            self.skipTest("Run from the generated project overlay")
        code = actual.read_text()
        self.assertIn("treeResult = UsingFrontEnd[RunT3Matching[transitionTree, 5, 0]]", code)
        self.assertIn("loopResult = UsingFrontEnd[RunT3Matching[transitionTree, 5, 1]]", code)
        self.assertIn("probeResult = UsingFrontEnd[RunT3Matching[transitionFullCanonical, 5, 0]]", code)
        self.assertIn("EFT1D2WilsonRunningHeavyInsertion === 0", code)
        self.assertNotIn("RunT3Matching[hardEFT + directC5", code)


    def test_scalar_matching_uses_service_frontend_as_successful_stage1(self):
        # Stage 1 uses UsingFrontEnd[RunT3Matching[...]]; the second
        # threshold must not run Match[] in a kernel without that wrapper.
        actual = (Path(__file__).resolve().parents[2] /
                  "Lagrangian/RunGeneralDMinus2ScalarThreshold.wl")
        code = actual.read_text(encoding="utf-8")
        self.assertIn("UsingFrontEnd[\n        Match[profileInput, EFTOrder -> 5, LoopOrder -> 0]", code)
        self.assertIn("treeResult = UsingFrontEnd[RunT3Matching[transitionTree, 5, 0]]", code)
        self.assertIn("loopResult = UsingFrontEnd[RunT3Matching[transitionTree, 5, 1]]", code)
        self.assertIn("probeResult = UsingFrontEnd[RunT3Matching[transitionFullCanonical, 5, 0]]", code)
        self.assertIn("profileLimit, $TimedOut", code)
        self.assertIn('"DirectWeinbergInsertedIntoMatch" -> False', code)

    def test_bsm_probe_sequence_and_commands(self):
        self.assertEqual(PROFILE_BSM_PROBES, (
            "inspect", "bsm-free", "bsm-standard", "bsm-general",
            "bsm-with-y2", "bsm-with-y2Conj", "bsm-with-both",
        ))
        for probe in PROFILE_BSM_PROBES:
            command = make_profile_command(self.matching, self.wl, probe, 45,
                                           wolframscript="fake-math", script=self.wl)
            self.assertEqual(command[-2:], [f"profile:{probe}", "45"])

    def test_bsm_mode_runs_matching_in_order_when_baselines_return(self):
        def fake_probe(matching_dir, insertion, probe, *args):
            return {"probe": probe,
                    "status": "Inspected" if probe == "inspect" else "MatchReturned"}
        with patch("scripts.RunGeneralDMinus2ScalarThreshold._execute_profile_probe",
                   side_effect=fake_probe) as mocked:
            result = run_tree_profile(self.matching, script=self.wl, bsm_only=True)
        self.assertEqual(result["status"], "BSMProfileRecorded")
        self.assertEqual([entry["probe"] for entry in result["probes"]],
                         list(PROFILE_BSM_PROBES))
        self.assertEqual(mocked.call_count, len(PROFILE_BSM_PROBES))
        self.assertFalse(result["full_general_threshold_matched"])
        self.assertTrue(result["bsm_only"])
        self.assertTrue((self.data / "general_d_m2_scalar_threshold/profile_tree/profile_bsm_summary.json").exists())

    def test_bsm_profile_skips_expensive_wilson_probes_after_baseline_timeout(self):
        def fake_probe(matching_dir, insertion, probe, *args):
            return {"probe": probe, "status": "Inspected" if probe == "inspect"
                    else "TimedOut" if probe == "bsm-general" else "MatchReturned"}
        with patch("scripts.RunGeneralDMinus2ScalarThreshold._execute_profile_probe",
                   side_effect=fake_probe) as mocked:
            result = run_tree_profile(self.matching, script=self.wl, bsm_only=True)
        self.assertEqual(result["status"], "BSMProfileRecorded")
        self.assertEqual(result["baseline_bottleneck_probe"], "bsm-general")
        self.assertEqual(mocked.call_count, 4)
        self.assertEqual(result["skipped_probes"], list(PROFILE_BSM_PROBES)[4:])

    def test_bsm_wolfram_inputs_constructed_without_standalone_lsm(self):
        src = (Path(__file__).resolve().parents[2] /
               "Lagrangian/RunGeneralDMinus2ScalarThreshold.wl").read_text()
        self.assertIn('profileBSMFree = profileHeavyOnly;', src)
        self.assertIn('profileBSMStandard = RelabelIndices[Expand[', src)
        self.assertIn('profileBSMGeneral = RelabelIndices[Expand[', src)
        self.assertIn('profileBSMWithBoth = RelabelIndices[Expand[', src)
        self.assertIn('"bsm-general", profileBSMGeneral', src)
        self.assertIn('"bsm-with-both", profileBSMWithBoth', src)
        self.assertIn('"full-tree", transitionTree', src)
        self.assertIn('treeResult = UsingFrontEnd[RunT3Matching[transitionTree, 5, 0]]', src)
        self.assertIn('loopResult = UsingFrontEnd[RunT3Matching[transitionTree, 5, 1]]', src)

    def test_smoke_probe_order_and_commands(self):
        self.assertEqual(PROFILE_SMOKE_PROBES,
            ("inspect", "sm-only", "s2r-only", "s1-only", "heavy-only", "free"))
        for probe in PROFILE_SMOKE_PROBES:
            cmd = make_profile_command(self.matching, self.wl, probe, 30,
                                       wolframscript="fake-math", script=self.wl)
            self.assertEqual(cmd[-2:], [f"profile:{probe}", "30"])

    def test_smoke_runs_all_independent_probes_even_if_one_times_out(self):
        def fake_probe(matching_dir, insertion, probe, *args):
            return {"probe": probe, "status": "Inspected" if probe == "inspect"
                    else "TimedOut" if probe == "sm-only" else "MatchReturned"}
        with patch("scripts.RunGeneralDMinus2ScalarThreshold._execute_profile_probe",
                   side_effect=fake_probe) as mocked:
            result = run_tree_profile(self.matching, script=self.wl,
                                      smoke_only=True)
        self.assertEqual(result["status"], "SmokeProfileRecorded")
        self.assertEqual(mocked.call_count, len(PROFILE_SMOKE_PROBES))
        self.assertFalse(result["full_general_threshold_matched"])
        self.assertTrue(result["smoke_only"])
        self.assertTrue((self.data / "general_d_m2_scalar_threshold/profile_tree/profile_smoke_summary.json").exists())

    def test_smoke_wolfram_sectors_do_not_modify_full_match(self):
        actual = (Path(__file__).resolve().parents[2] /
                  "Lagrangian/RunGeneralDMinus2ScalarThreshold.wl")
        code = actual.read_text(encoding="utf-8")
        for probe in ("sm-only", "s2r-only", "s1-only", "heavy-only"):
            self.assertIn(f'"{probe}",', code)
        self.assertIn('profileSMOnly = CheckAbort[', code)
        self.assertIn('profileS2ROnly = CheckAbort[', code)
        self.assertIn('profileS1Only = CheckAbort[', code)
        self.assertIn('profileHeavyOnly = CheckAbort[', code)
        self.assertIn('profileResultPath = FileNameJoin[', code)
        self.assertIn('treeResult = UsingFrontEnd[RunT3Matching[transitionTree, 5, 0]]', code)

    def test_incremental_default_probe_order(self):
        self.assertEqual(PROFILE_PROBES, (
            "inspect", "free", "scalar-standard", "scalar-general",
            "scalar-with-y2", "scalar-with-y2Conj",
        ))

    def test_incremental_probe_command_preserves_flags(self):
        for probe in PROFILE_PROBES:
            cmd = make_profile_command(self.matching, self.wl, probe, 60,
                                       wolframscript="fake-math", script=self.wl)
            self.assertEqual(cmd[-2:], [f"profile:{probe}", "60"])

    def test_profile_stops_if_standard_quartics_timeout(self):
        def fake_probe(matching_dir, insertion, probe, *args):
            return {"probe": probe,
                    "status": "Inspected" if probe == "inspect" else (
                        "TimedOut" if probe == "scalar-standard" else "MatchReturned")}
        with patch("scripts.RunGeneralDMinus2ScalarThreshold._execute_profile_probe",
                   side_effect=fake_probe) as mocked:
            result = run_tree_profile(self.matching, script=self.wl)
        self.assertEqual(result["status"], "TreeProfileRecorded")
        self.assertEqual(result["baseline_bottleneck_probe"], "scalar-standard")
        self.assertEqual(mocked.call_count, 3)
        self.assertEqual(result["skipped_probes"], list(PROFILE_PROBES)[3:])

    def test_profile_reaches_both_wilson_sectors_only_if_scalar_general_returns(self):
        def fake_probe(matching_dir, insertion, probe, *args):
            return {"probe": probe,
                    "status": "Inspected" if probe == "inspect" else "MatchReturned"}
        with patch("scripts.RunGeneralDMinus2ScalarThreshold._execute_profile_probe",
                   side_effect=fake_probe) as mocked:
            result = run_tree_profile(self.matching, script=self.wl)
        self.assertEqual(result["status"], "TreeProfileRecorded")
        self.assertEqual(mocked.call_count, len(PROFILE_PROBES))
        self.assertNotIn("baseline_bottleneck_probe", result)

    def test_profile_strict_physical_scalar_sector_is_constructed_from_lsm_free(self):
        actual = (Path(__file__).resolve().parents[2] /
                  "Lagrangian/RunGeneralDMinus2ScalarThreshold.wl")
        code = actual.read_text(encoding="utf-8")
        self.assertIn('uvBuild["LSM"] + FreeLag[NewScalar1, T3S2R, T3S2I]', code)
        self.assertIn('D2ProfileScalarDegree[term_]', code)
        self.assertIn('Power[Field[_, Scalar, ___], n_Integer?Positive]', code)
        self.assertIn('D2ProfileScalarDegree[profileTerms[[#]]] === 4', code)
        self.assertIn('"MissingGeneralQuarticCouplings" -> profileMissingExtra', code)
        self.assertIn('If[Length[profileQuarticPositions] === 0', code)
        self.assertIn('"free", profileFree', code)
        self.assertIn('"scalar-standard", profileStandard', code)
        self.assertIn('"scalar-general", profileGeneral', code)
        self.assertIn('"scalar-with-y2", RelabelIndices[Expand[profileGeneral + profileY2]]', code)
        self.assertIn('"scalar-with-y2Conj", RelabelIndices[Expand[profileGeneral + profileY2Conj]]', code)
        self.assertIn('"full-tree", transitionTree', code)
        self.assertIn('"DirectWeinbergInsertedIntoMatch" -> False', code)


if __name__ == "__main__":
    unittest.main()
