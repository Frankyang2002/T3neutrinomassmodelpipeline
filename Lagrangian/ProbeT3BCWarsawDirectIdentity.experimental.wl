(*
  Workstream C: direct native/Warsaw LLHH operator identity, v11.2.

  Full replacement:
    Lagrangian/ProbeT3BCWarsawDirectIdentity.experimental.wl

  The installed Matchete v0.6.0 source Models/SMEFT_Warsaw.m still contains
  a legacy NonCommutativeMultiply (**) expression for its dimension-five
  spinor bilinear. Matchete v0.4+ uses NCM instead of ** (release notes).
  Evaluating the old syntax gives NCM[...] ** NCM[...], not one spinor line.

  This probe constructs the EXACT SAME ordered operator with NCM rather
  than ** and tests whether Matchete reduces it to the saved native PL chain.
  Both signs are tested. A missing equality does not become a certificate.

  This probe does not perform EFT matching, EFT counterterms or EFT running.
*)
ClearAll["Global`*"];
args = Rest[$ScriptCommandLine];
If[Length[args] =!= 1,
  Print["Usage: wolframscript -file Lagrangian/ProbeT3BCWarsawDirectIdentity.experimental.wl OUTPUT.json"];
  Exit[2]
];
outPath = ExpandFileName[First[args]];
If[!DirectoryQ[DirectoryName[outPath]],
  CreateDirectory[DirectoryName[outPath], CreateIntermediateDirectories -> True]
];
format = "t3-bc-warsaw-direct-operator-identity-probe-v11";
writeReport[status_, fields_, code_] := Module[{report, exportResult},
  report = Join[
    <|
      "format" -> format,
      "probe_revision" -> "11.2",
      "direct_equivalence_status" -> status,
      "matched_EFT" -> False,
      "physical_C5" -> Null
    |>,
    fields
  ];
  exportResult = Quiet @ Check[Export[outPath, report, "RawJSON"], $Failed];
  If[exportResult === $Failed,
    Print["T3_WARSAW_DIRECT: REPORT EXPORT FAILED -> ", outPath];
    Exit[9]
  ];
  Print["T3_WARSAW_DIRECT: ", status, " -> ", outPath];
  Exit[code]
];

Print["T3_WARSAW_DIRECT: frontend-aware-Matchete-load"];
load = CheckAbort[
  Quiet @ Check[UsingFrontEnd[Needs["Matchete`"]; True], $Failed],
  $Aborted
];
If[load =!= True,
  writeReport["MatcheteLoadFailed", <|
    "initialization" -> "UsingFrontEnd[Needs[Matchete`]]",
    "loader_result" -> ToString[load, InputForm]
  |>, 3]
];

Print["T3_WARSAW_DIRECT: ResetAll-before-SM-load"];
reset = CheckAbort[ResetAll[]; True, $Aborted];
If[reset =!= True,
  writeReport["ResetAllFailed", <|
    "initialization" -> "ResetAll[]",
    "reset_result" -> ToString[reset, InputForm]
  |>, 4]
];
Print["T3_WARSAW_DIRECT: LoadModel[SM], retaining diagnostics"];
sm = CheckAbort[LoadModel["SM"], $Aborted];
loadMode = "ResetAllThenLoadModel";
If[sm === $Failed || sm === $Aborted,
  Print["T3_WARSAW_DIRECT: SM load failed; retry after ResetAll in frontend"];
  retryReset = CheckAbort[ResetAll[]; True, $Aborted];
  If[retryReset === True,
    sm = CheckAbort[UsingFrontEnd[LoadModel["SM"]], $Aborted];
    loadMode = "ResetAllThenFrontendLoadModel"
  ];
];
If[sm === $Failed || sm === $Aborted || !FreeQ[sm, $Failed, {0, Infinity}],
  writeReport["SMModelLoadFailed", <|
    "initialization" -> loadMode,
    "model_load_result" -> StringTake[ToString[sm, InputForm], UpTo[1200]],
    "Mathematica_version" -> ToString[$Version]
  |>, 4]
];
Print["T3_WARSAW_DIRECT: built-in SM model loaded; head=", Head[sm]];

(* The legacy Warsaw code is recorded ONLY to expose the root cause.
   It cannot be compared directly with NCM: ** is now a different head.
   No model file is edited or reloaded. *)
Print["T3_WARSAW_DIRECT: Warsaw NCM construction"];
gauge = Bar@eps[SU2L][i,j] Bar@eps[SU2L][k,m] H[i] H[k];
warsawLegacy = gauge (Bar@CConj@l[j,p] ** l[m,r]);
warsawNCM = gauge NCM[Bar@CConj@l[j,p], l[m,r]];
native = gauge NCM[
  Transp@l[j,p],
  DiracProduct[GammaCC, Proj[-1]],
  l[m,r]
];

(* Fail closed if migration to NCM did not actually eliminate the old head.
   Confirm that the legacy expression still contains the old head; otherwise
   the installed Matchete conventions differ and this experiment needs review.
*)
legacyHasOldHead = !FreeQ[warsawLegacy, NonCommutativeMultiply, {0, Infinity}];
canonicalHasOldHead = !FreeQ[warsawNCM, NonCommutativeMultiply, {0, Infinity}];
If[!legacyHasOldHead || canonicalHasOldHead,
  writeReport["UnexpectedNoncommutativeHead", <|
    "legacy_contains_noncommutative_multiply" -> legacyHasOldHead,
    "NCM_reference_contains_noncommutative_multiply" -> canonicalHasOldHead,
    "reference_legacy_inputform" -> StringTake[ToString[warsawLegacy, InputForm], UpTo[4500]],
    "reference_ncm_inputform" -> StringTake[ToString[warsawNCM, InputForm], UpTo[4500]]
  |>, 5]
];

(* A zero expression or wrong-sign match cannot be silently accepted. *)
safeReduce[expr_] := Module[{x},
  x = Quiet @ Check[Expand[expr], expr];
  x = Quiet @ Check[RefineDiracProducts[x], x];
  x = Quiet @ Check[Expand[x], x];
  x = Quiet @ Check[RelabelIndices[x], x];
  x
];
refReduced = safeReduce[warsawNCM];
nativeReduced = safeReduce[native];
plus = safeReduce[refReduced - nativeReduced];
minus = safeReduce[refReduced + nativeReduced];
nontrivial = !TrueQ[refReduced === 0] && !TrueQ[nativeReduced === 0];
plusZero = TrueQ[plus === 0];
minusZero = TrueQ[minus === 0];
status = Which[
  !nontrivial, "DegenerateZeroOperator",
  plusZero && !minusZero, "ExactPlusEquality",
  minusZero && !plusZero, "ExactMinusEquality",
  plusZero && minusZero, "BothZeroOrDegenerate",
  True, "NotProvenByNativeSimplifier"
];
limit[s_] := StringTake[s, UpTo[5000]];
writeReport[status, <|
  "Mathematica_version" -> ToString[$Version],
  "source_model" -> "SM",
  "source_model_load_mode" -> loadMode,
  "source_model_load_head" -> ToString[Head[sm], InputForm],
  "reference_warsaw_construction" ->
    "Installed SMEFT_Warsaw.m Lagd5, same field/gauge ordering; legacy ** replaced by Matchete 0.6 NCM",
  "native_construction" -> "Transp[l] GammaCC Proj[-1] l with identical epsilon and Higgs slots",
  "legacy_contains_noncommutative_multiply" -> legacyHasOldHead,
  "NCM_reference_contains_noncommutative_multiply" -> canonicalHasOldHead,
  "nontrivial_reference_and_native" -> nontrivial,
  "positive_sign_equality" -> plusZero,
  "negative_sign_equality" -> minusZero,
  "legacy_internal_inputform" -> limit @ ToString[warsawLegacy, InputForm],
  "reference_ncm_inputform" -> limit @ ToString[refReduced, InputForm],
  "native_internal_inputform" -> limit @ ToString[nativeReduced, InputForm],
  "positive_difference_inputform" -> limit @ ToString[plus, InputForm],
  "negative_difference_inputform" -> limit @ ToString[minus, InputForm],
  "does_not_certify_c5" -> True,
  "production_modified" -> False
|>, If[status === "ExactPlusEquality", 0, 5]];
