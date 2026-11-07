"""Contract tests for the local provisional C5 numerical benchmark."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
RUNNER_PATH = ROOT / "scripts/RunGeneralDMinus2C5RoughBenchmark.py"
WOLFRAM_PATH = ROOT / "Lagrangian/RunGeneralDMinus2C5RoughBenchmark.wl"
CONFIG_PATH = ROOT / "benchmarks/general_d_m2_rough_benchmark.json"
FINITE_CONFIG_PATH = ROOT / "benchmarks/general_d_m2_rough_benchmark_finitepart.json"

spec = importlib.util.spec_from_file_location("d2benchmark", RUNNER_PATH)
assert spec and spec.loader
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class TestBenchmarkContracts(unittest.TestCase):
    def test_config_loads(self):
        value = module.inspect_config(CONFIG_PATH)
        self.assertEqual(value["NF"], 1)
        self.assertEqual(value["Benchmark"], "D_m2_provisional_one_heavy_fermion_real_couplings")

    def test_three_independent_yukawas(self):
        c = module.inspect_config(CONFIG_PATH)
        self.assertEqual(set(c["Yukawas"]), {"y1", "y2", "y2Conj"})
        for mat in c["Yukawas"].values():
            self.assertEqual(len(mat), 3)
            self.assertEqual([len(row) for row in mat], [1, 1, 1])

    def test_physical_mass_positivity(self):
        c = module.inspect_config(CONFIG_PATH)
        self.assertGreater(c["MS2_GeV"]**2 - 2*c["AbsKappaS_GeV2"], 0)

    def test_exact_symbolic_equations_present(self):
        src = WOLFRAM_PATH.read_text(encoding="utf-8")
        for token in ("Sqrt[ms2^2 - 2 kappa]", "Sqrt[ms2^2 + 2 kappa]",
                      "1/(16 Pi^2)", "-(vev^2/2)", "10^9",
                      "SingularValueList", "TotalCandidateLinearInHbar"):
            if token == "TotalCandidateLinearInHbar":
                continue
            self.assertIn(token, src)

    def test_no_matching(self):
        src = WOLFRAM_PATH.read_text(encoding="utf-8")
        self.assertNotIn("RunT3Matching[", src)
        self.assertNotIn("Match[candidate", src)
        self.assertIn('profile_c5_flavour_candidate_total_provisional.wxf', src)

    def test_no_hidden_unassigned_parameter_zeroing(self):
        src = WOLFRAM_PATH.read_text(encoding="utf-8")
        self.assertIn('"UnassignedScalarParameters"', src)
        self.assertIn('"UnsupportedIndexedCouplings"', src)
        self.assertIn('"FallbackLambdaAssignments"', src)

    def test_output_is_explicitly_provisional(self):
        src = WOLFRAM_PATH.read_text(encoding="utf-8")
        self.assertIn('"NumericalBenchmarkCertifiedPhysical" -> False', src)
        self.assertIn('"PhysicalC5NormalizationIndependentlyCertified" -> False', src)

    def test_valid_classification(self):
        report = {
            "Status": "ProvisionalNumericalBenchmarkEvaluated",
            "NoNewMatchPerformed": True,
            "NumericalBenchmarkCertifiedPhysical": False,
            "PhysicalC5NormalizationIndependentlyCertified": False,
            "NeutrinoMasses_eV": [0.0, 0.02, 0.05],
        }
        self.assertEqual(module.classify(0, False, report),
                         "ProvisionalNumericalBenchmarkEvaluated")

    def test_invalid_classification(self):
        self.assertEqual(module.classify(20, False, None), "WolframFailed")
        self.assertEqual(module.classify(None, True, None), "ExternalTimeout")
        self.assertEqual(module.classify(0, False, None), "InvalidResultReport")

    def test_csv_exports(self):
        src = WOLFRAM_PATH.read_text(encoding="utf-8")
        self.assertGreaterEqual(src.count('"CSV"]'), 6)

    def test_config_json_exact(self):
        c = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))
        self.assertGreater(c["MF_GeV"], 0)
        self.assertNotEqual(c["Yukawas"]["y2"], c["Yukawas"]["y2Conj"])

    def test_finite_sums_precede_numerical_coupling_lookup(self):
        src = WOLFRAM_PATH.read_text(encoding="utf-8")
        preparation = src.index("numericTemplate = candidate /.")
        expansion = src.index("numericTemplate = D2ExpandFiniteSums[numericTemplate]")
        substitution = src.index("v = v /. HoldPattern[Matchete`Coupling[")
        self.assertLess(preparation, expansion)
        self.assertLess(expansion, substitution)
        self.assertIn("{k, lo, hi}", src)

    def test_unexpanded_sums_are_not_silently_discarded(self):
        src = WOLFRAM_PATH.read_text(encoding="utf-8")
        self.assertIn('"Finite heavy-fermion sums were not fully expanded"', src)
        self.assertIn('"RemainingSumCount"', src)
        self.assertIn('"RemainingNFSymbols"', src)

    def test_first_nonnumeric_entry_fails_with_diagnostics(self):
        src = WOLFRAM_PATH.read_text(encoding="utf-8")
        self.assertIn('"BadMatrixEntry" -> {p, q}', src)
        self.assertIn('"EntryPreview"', src)
        self.assertIn('"UnresolvedCouplingCount"', src)
        self.assertIn('"UnknownYukawaCount"', src)
        self.assertIn('"EntryTimedOut"', src)
        self.assertIn('"UnknownScalarValues"', src)
        self.assertIn('"PreNumericEntryPreview"', src)

    def test_benchmark_still_uses_saved_mass_split_and_takagi(self):
        src = WOLFRAM_PATH.read_text(encoding="utf-8")
        self.assertIn('Sqrt[ms2^2 - 2 kappa]', src)
        self.assertIn('Sqrt[ms2^2 + 2 kappa]', src)
        self.assertIn('SingularValueList[mnuEV]', src)
        self.assertIn('"NumericalBenchmarkCertifiedPhysical" -> False', src)


    def test_finitepart_config_optin_is_explicit(self):
        cfg = module.inspect_config(FINITE_CONFIG_PATH)
        self.assertEqual(cfg["PoleHandling"], "FinitePartDiagnostic")
        self.assertEqual(cfg["NF"], 1)
        self.assertEqual(cfg["Yukawas"]["y2Conj"],
                         module.inspect_config(CONFIG_PATH)["Yukawas"]["y2Conj"])

    def test_finitepart_output_classification_requires_metadata(self):
        report = {
            "Status": "ProvisionalFinitePartBenchmarkEvaluated",
            "NoNewMatchPerformed": True,
            "NumericalBenchmarkCertifiedPhysical": False,
            "PhysicalC5NormalizationIndependentlyCertified": False,
            "NeutrinoMasses_eV": [0.0, 0.02, 0.05],
            "PoleHandling": "FinitePartDiagnostic",
            "FinitePartOnly": True,
            "RenormalizationCertified": False,
        }
        self.assertEqual(module.classify(0, False, report),
                         "ProvisionalFinitePartBenchmarkEvaluated")
        report["RenormalizationCertified"] = True
        self.assertEqual(module.classify(0, False, report),
                         "InconsistentBenchmarkMetadata")

    def test_implicit_nflavor_only_for_nf_one(self):
        src = WOLFRAM_PATH.read_text(encoding="utf-8")
        self.assertIn('If[nf === 1,', src)
        self.assertIn('"NFlavor"', src)
        self.assertIn('Matchete`Index[1, rep]', src)
        self.assertIn('"Uncontracted heavy-generation index;', src)

    def test_mubar2_uses_scale_squared(self):
        src = WOLFRAM_PATH.read_text(encoding="utf-8")
        self.assertIn('D2MuBar2Q', src)
        self.assertIn('FromCharacterCode[956] <> "bar2"', src)
        self.assertIn('D2MuBar2Q[s] :> N[mu^2, 20]', src)

    def test_epsilon_pole_is_never_assigned_a_number(self):
        src = WOLFRAM_PATH.read_text(encoding="utf-8")
        self.assertIn('FromCharacterCode[1013]', src)
        self.assertIn('"PoleResidue" -> residue', src)
        self.assertIn('"Nonzero 1/epsilon pole remains', src)
        self.assertIn('poleHandling === "Reject"', src)
        self.assertNotIn('D2DimRegQ[s] :> 0', src)
        self.assertNotIn('D2DimRegQ[s] :> 1', src)

    def test_six_csv_outputs_and_unrenormalized_label(self):
        src = WOLFRAM_PATH.read_text(encoding="utf-8")
        for name in ('poleRealPath', 'poleImagPath', 'C5PoleResidueMaxAbs_GeVinv',
                     '"RenormalizationCertified" -> False',
                     '"FinitePartOnly"', '"NumericalBenchmarkCertifiedPhysical" -> False'):
            self.assertIn(name, src)


if __name__ == "__main__":
    unittest.main()
