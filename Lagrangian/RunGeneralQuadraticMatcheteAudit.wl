(* Experimental 16-model general-T3 scalar operator diagnostic.

   wolframscript -file Lagrangian/RunGeneralQuadraticMatcheteAudit.wl \
       OUTPUT.json CLASS alpha

   Matchete's CheckLagrangian returns False for the B/C(-1) gauge-allowed
   holomorphic scalar mixing because HeavyMassBasisQ requires a diagonal
   heavy mass basis.  Simply declaring S1/S2 Light did NOT suppress that
   check in live Matchete (2026-10-27).  Do not attempt reclassification.

   Instead request DetailedOutput->True and require explicit core checks.
   Matchete-version-dependent opaque diagnostic entries must match those of
   an independently passing ordinary CheckLagrangian baseline.  For B/C(-1), additionally require
   an explicit non-passing HeavyMassBasis diagnostic.  The result is an
   operator-consistency certificate only, NOT a matching certificate.

   No Match[], no changes to production field definitions.
*)
ClearAll["Global`*"];
scriptDir = DirectoryName @ ExpandFileName[$InputFileName];
args = Rest[$ScriptCommandLine];
If[Length[args] =!= 3,
  Print["Usage: OUTPUT.json CLASS alpha"]; Exit[2]];
out = ExpandFileName[args[[1]]];
class = ToUpperCase[args[[2]]];
T3ParseAlpha[s_String] := Which[
  StringMatchQ[s, "m" ~~ DigitCharacter ..], -ToExpression[StringDrop[s, 1]],
  StringMatchQ[s, "p" ~~ DigitCharacter ..], ToExpression[StringDrop[s, 1]],
  StringMatchQ[s, ("-" | "") ~~ DigitCharacter ..], ToExpression[s],
  True, $Failed
];
alpha = T3ParseAlpha[args[[3]]];
If[alpha === $Failed, Print["Invalid alpha"]; Exit[2]];
load = CheckAbort[Check[UsingFrontEnd[Needs["Matchete`"]; True], $Failed], $Aborted];
If[load =!= True, Print["Matchete unavailable"]; Exit[3]];
Get[FileNameJoin[{scriptDir, "model", "T3ModelCatalog.wl"}]];
Get[FileNameJoin[{scriptDir, "model", "LagrangianBuilder.wl"}]];
model = T3ModelFromClass[class, alpha];
If[model === $Failed, Print["Unsupported model"]; Exit[2]];

(* CheckLagrangian[DetailedOutput->True] typically uses STRING keys:
   <|"Hermiticity"->True, ...|>.  ToString["Hermiticity", InputForm]
   adds quotes, so the previous parser silently changed the key into
   "\"Hermiticity\"" and reported a spurious BaselineDetailedCheckFailed.

   Prefer the exact string key.  For symbols from another package version,
   remove the context prefix; trim InputForm quotes defensively. *)
T3DiagnosticKeyName[key_] := Module[{name},
  If[StringQ[key], Return[key]];
  name = Last @ StringSplit[ToString[Unevaluated[key], InputForm], "`"];
  StringReplace[name, {
    StartOfString ~~ "\"" -> "",
    "\"" ~~ EndOfString -> ""
  }]
];

(* Pure string-key regression runs inside the live Wolfram kernel before
   any Matchete operator validation.  Also test Matchete's slash-containing
   compound diagnostic label. *)
If[T3DiagnosticKeyName["Hermiticity"] =!= "Hermiticity" ||
   T3DiagnosticKeyName["Open/Complex/InconsistentSpinChains"] =!=
     "Open/Complex/InconsistentSpinChains",
  Print["Internal Matchete diagnostic key normalization FAILED"];
  Exit[4]
];

(* Baseline-calibrated interpreter for the INSTALLED Matchete version.
   Live evidence from the user's Matchete installation:
     ChargeNeutral is absent, but Charges is present;
     GaugeAnomalies and ZeroModeMixing have non-Boolean values.
   Those are NOT necessarily failed diagnostics. We do not invent their
   semantics; an ordinary CheckLagrangian[baseline] == True certifies the
   original values, and subsequent operators must reproduce them exactly.
   We additionally insist on explicit Boolean/empty values for the known
   structural tests. *)
T3CheckNonMassDiagnostics[lag_, opaqueBaseline_: None] := Module[
  {ordinary, raw, d, required, expected, opaque, opaqueStrings, sameOpaque,
   gaugeOK, chargesOK, failures, summary, basis, genuinePlainResult},
  (* Matchete emits a warning for the expected mixed heavy mass term.
     Check[..., $Failed] would discard the actual False result. *)
  ordinary = CheckAbort[CheckLagrangian[lag], $Aborted];
  raw = CheckAbort[CheckLagrangian[lag, DetailedOutput -> True], $Aborted];
  If[!AssociationQ[raw],
    Return[<|"diagnostic_returned" -> False,
      "non_mass_basis_checks_passed" -> False,
      "standard_check_passed" -> False,
      "heavy_mass_basis_passed" -> False,
      "heavy_mass_basis_reported" -> False,
      "opaque_metadata_baseline_match" -> False,
      "failed_checks" -> {"DetailedOutputNotAssociation"},
      "diagnostic_type" -> ToString[raw, InputForm]|>]
  ];
  d = Association @ KeyValueMap[(T3DiagnosticKeyName[#1] -> #2) &, raw];

  expected = <|
    "Open/Complex/InconsistentSpinChains" -> {},
    "Hermiticity" -> True,
    "UncontractedIndices" -> {},
    "CanonicallyNormalized" -> True,
    "HeavyTadpoles" -> {},
    "AllObjectsDefined" -> True
  |>;
  required = Join[Keys[expected], {"HeavyMassBasis"}];
  failures = Select[Join[required, {"GaugeAnomalies"}], !KeyExistsQ[d, #] &];
  failures = Join[failures,
    Select[Keys[expected],
      KeyExistsQ[d, #] && !SameQ[d[#], expected[#]] &]
  ];
  gaugeOK = (KeyExistsQ[d, "GaugeInvariant"] &&
      TrueQ[d["GaugeInvariant"]]) ||
    (KeyExistsQ[d, "FreeOfGaugeFields"] &&
      TrueQ[d["FreeOfGaugeFields"]]);
  If[!gaugeOK, AppendTo[failures, "GaugeInvariantOrFreeOfGaugeFields"]];
  chargesOK = KeyExistsQ[d, "ChargeNeutral"] || KeyExistsQ[d, "Charges"];
  If[!chargesOK, AppendTo[failures, "ChargeNeutralOrChargesMissing"]];
  If[KeyExistsQ[d, "ChargeNeutral"] && !TrueQ[d["ChargeNeutral"]],
    AppendTo[failures, "ChargeNeutral"]];

  (* The three observed non-Boolean keys are NOT interpreted as failures
     merely because their values are neither True nor {}. Instead demand
     exact agreement with a baseline independently accepted by Matchete's
     ordinary full-Lagrangian checker. This rule applies to ALL unknown keys. *)
  opaque = KeyDrop[d, Join[required, {"GaugeInvariant",
    "FreeOfGaugeFields", "ChargeNeutral"}]];
  opaqueStrings = Association @ KeyValueMap[
    (#1 -> ToString[#2, InputForm]) &, opaque];
  sameOpaque = opaqueBaseline === None ||
    (AssociationQ[opaqueBaseline] && SameQ[opaqueStrings, opaqueBaseline]);
  If[!sameOpaque, AppendTo[failures, "OpaqueMetadataDiffersFromValidBaseline"]];
  genuinePlainResult = ordinary === True || ordinary === False;
  If[!genuinePlainResult, AppendTo[failures, "OrdinaryCheckNotBoolean"]];
  failures = DeleteDuplicates[failures];

  basis = Lookup[d, "HeavyMassBasis", Missing["NotReported"]];
  summary = Association @ KeyValueMap[
    Function[{key, value}, key -> Which[
      KeyExistsQ[opaque, key], If[sameOpaque, "OpaqueMatchesValidBaseline",
        "OpaqueDiffersFromValidBaseline"],
      value === True, "PassTrue", value === {}, "PassEmpty",
      value === False, "FailFalse", True, "UnexpectedValue"
    ]], d];
  <|
    "diagnostic_returned" -> True,
    "standard_check_passed" -> TrueQ[ordinary],
    "standard_check_boolean" -> genuinePlainResult,
    "non_mass_basis_checks_passed" -> (failures === {}),
    "heavy_mass_basis_passed" -> TrueQ[basis],
    "heavy_mass_basis_reported" -> KeyExistsQ[d, "HeavyMassBasis"],
    "opaque_metadata_baseline_match" -> sameOpaque,
    "opaque_metadata_values_inputform" -> opaqueStrings,
    "failed_checks" -> failures,
    "diagnostics" -> summary,
    "diagnostic_key_inputforms" -> (ToString[#, InputForm] & /@ Keys[raw]),
    "diagnostic_key_names" -> Keys[d],
    "comparison_method" -> "ExplicitCoreChecksPlusOrdinaryBaselineAnchoredOpaqueMetadata"
  |>
];

(* Emit the actual failed diagnostic *in the Wolfram log* as well as JSON.
   This is diagnostic only and does not affect the validity criteria. *)
T3PrintDiagnosticFailure[stage_String, diagnostic_] := Module[{},
  Print[stage, ": diagnostic audit FAILED"];
  Print["Actual diagnostic report: ", ToString[diagnostic, InputForm]];
];

T3AuditDetailedValid[diag_Association, mixed_] := TrueQ[
  diag["diagnostic_returned"] &&
  diag["non_mass_basis_checks_passed"] &&
  diag["opaque_metadata_baseline_match"] &&
  diag["standard_check_boolean"] &&
  (TrueQ[diag["standard_check_passed"]] === !TrueQ[mixed]) &&
  diag["heavy_mass_basis_reported"] &&
  If[TrueQ[mixed], !TrueQ[diag["heavy_mass_basis_passed"]],
    TrueQ[diag["heavy_mass_basis_passed"]]]
];

ClearAll[T3GeneralMassAudit];
T3GeneralMassAudit[] := Module[
  {built, extraQ, structures, extraM, expected, isMixed, base,
   baseline, individualQ, individualM, combined, opaqueBaseline,
   status, combinedLag, quadraticNames, quarticNames, failure, checksOK, matrixBasis, extraQNames},
  If[!MatchQ[$FrontEnd, _FrontEndObject],
    Return[<|"status" -> "FrontEndUnavailable"|>]];

  (* Baseline heavy-field assembly includes the ordinary interaction basis.
     The general quadratic coupling is generated only AFTER that passes. *)
  built = BuildT3Lagrangian[model];
  If[!AssociationQ[built] || built["Status"] =!= "Success",
    Return[<|"status" -> "BaselineHeavyBuildFailed"|>]];
  base = built["LUV"];
  extraQ = T3ExtendedBuildExtraQuarticCandidates[model];
  If[!ListQ[extraQ],
    Return[<|"status" -> "QuarticBasisGenerationFailed"|>]];
  structures = T3ExtendedExtraQuadraticStructures[model];
  isMixed = MemberQ[{{"B", -1}, {"C", -1}}, {class, alpha}];
  expected = If[MemberQ[{{"A", 0}, {"A", -2}, {"B", -1}, {"C", -1},
    {"D", 0}, {"D", -2}, {"E", 0}, {"E", -2}}, {class, alpha}], 1, 0];
  If[Length[structures] =!= expected,
    Return[<|"status" -> "QuadraticStructureInventoryMismatch",
      "found" -> Length[structures], "expected" -> expected|>]];
  extraM = T3ExtendedBuildCandidates[model, 2, structures];
  If[!ListQ[extraM] || Length[extraM] =!= expected,
    Return[<|"status" -> "QuadraticOperatorGenerationFailed"|>]];
  quadraticNames = Lookup[extraM, "Name", {}];
  quarticNames = Lookup[extraQ, "Name", {}];

  baseline = T3CheckNonMassDiagnostics[base];
  If[!T3AuditDetailedValid[baseline, False],
    T3PrintDiagnosticFailure["Baseline", baseline];
    Return[<|"status" -> "BaselineDetailedCheckFailed",
      "baseline_diagnostic" -> baseline|>]];
  opaqueBaseline = baseline["opaque_metadata_values_inputform"];

  (* Quartics alone must remain in a valid heavy mass basis. *)
  individualQ = Table[
    With[{diag = T3CheckNonMassDiagnostics[
      (base + extraQ[[i, "Expression"]]) // RelabelIndices, opaqueBaseline]},
      <|"name" -> extraQ[[i, "Name"]],
        "valid" -> T3AuditDetailedValid[diag, False],
        "diagnostic" -> diag|>], {i, Length[extraQ]}];
  If[!And @@ Lookup[individualQ, "valid", {}],
    T3PrintDiagnosticFailure["Individual quartic", Select[individualQ, !TrueQ[#["valid"]] &]];
    Return[<|"status" -> "IndividualQuarticDetailedCheckFailed",
      "quartic_individual_checks" -> individualQ|>]];

  (* For B/C(-1), the one expected failure is precisely HeavyMassBasis.
     No field light/heavy reclassification, no ignored gauge warnings. *)
  If[isMixed, Print["B/C(-1) mass mixing: checking every Matchete ",
    "DetailedOutput component except HeavyMassBasis."]];
  individualM = Table[
    With[{diag = T3CheckNonMassDiagnostics[
      (base + extraM[[i, "Expression"]]) // RelabelIndices, opaqueBaseline]},
      <|"name" -> extraM[[i, "Name"]],
        "valid" -> T3AuditDetailedValid[diag, isMixed],
        "diagnostic" -> diag|>], {i, Length[extraM]}];
  If[!And @@ Lookup[individualM, "valid", {}],
    T3PrintDiagnosticFailure["Individual quadratic", Select[individualM, !TrueQ[#["valid"]] &]];
    Return[<|"status" -> "IndividualQuadraticDetailedCheckFailed",
      "quadratic_individual_checks" -> individualM|>]];

  combinedLag = (base + Total[Lookup[extraQ, "Expression", {}]] +
    Total[Lookup[extraM, "Expression", {}]]) // Expand // RelabelIndices;
  combined = T3CheckNonMassDiagnostics[combinedLag, opaqueBaseline];
  If[!T3AuditDetailedValid[combined, isMixed],
    T3PrintDiagnosticFailure["Combined Lagrangian", combined];
    Return[<|"status" -> "CombinedDetailedCheckFailed",
      "combined_diagnostic" -> combined|>]];

  status = If[isMixed,
    "MatcheteNonMassBasisDetailedChecksPassedHeavyMassBasisBlocked",
    "MatcheteFullDetailedChecksPassed"];
  <|
    "status" -> status,
    "class" -> class, "alpha" -> alpha,
    "baseline_heavy_lagrangian_valid" -> True,
    "additional_quartics" -> Length[extraQ],
    "quadratic_slots" -> Length[extraM],
    "quartic_operators" -> quarticNames,
    "quadratic_operators" -> quadraticNames,
    "quartic_individual_checks" -> individualQ,
    "quadratic_individual_checks" -> individualM,
    "baseline_diagnostic" -> baseline,
    "combined_diagnostic" -> combined,
    "combined_non_mass_basis_checks_passed" -> True,
    "baseline_calibrated_diagnostic_parser" -> True,
    "opaque_diagnostics_uninterpreted_but_baseline_identical" -> True,
    "combined_heavy_mass_basis_passed" -> !isMixed,
    "checked_in_original_heavy_field_definitions" -> True,
    "operator_consistency_validated_except_mass_basis" -> isMixed,
    "matchete_full_lagrangian_valid" -> !isMixed,
    "mass_basis_diagonalized" -> False,
    "heavy_scalar_matching_mass_basis_usable" -> False,
    "mixed_heavy_mass_term_requires_diagonalization" -> isMixed,
    "operator_normalization_certified" -> False,
    "rgbeta_physical_mapping_certified" -> False,
    "mass_eigenbasis_and_matching_done" -> False,
    "matching_executed" -> False
  |>
];
result = CheckAbort[UsingFrontEnd[T3GeneralMassAudit[]], $Aborted];
If[!AssociationQ[result], result = <|"status" -> "Aborted"|>];
If[!DirectoryQ[DirectoryName[out]],
  CreateDirectory[DirectoryName[out], CreateIntermediateDirectories -> True]];
Export[out, result, "RawJSON"];
Print[class, "_", alpha, ": ", result["status"], " -> ", out];
If[!MemberQ[{"MatcheteFullDetailedChecksPassed",
  "MatcheteNonMassBasisDetailedChecksPassedHeavyMassBasisBlocked"},
  result["status"]], Exit[5]];
