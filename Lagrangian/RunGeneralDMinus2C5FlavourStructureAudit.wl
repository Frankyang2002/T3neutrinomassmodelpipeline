(*
  General D(-2): read-only Weinberg operator structure and flavour-index audit.

  Source data: the successful saved full-EFT1 scalar one-loop processed hbar
  sector and the first-threshold EFT1 loop correction.  No Match[] call.
  This is NOT a flavour-index mapping or a certified C5 assembly.

  Usage: wolframscript -file Lagrangian/RunGeneralDMinus2C5FlavourStructureAudit.wl FIRST_THRESHOLD_OUTPUT_DIR
*)
ClearAll["Global`*"];
args = Rest[$ScriptCommandLine];
If[Length[args] =!= 1,
  Print["Usage: OUTPUT_DIR"]; Exit[2]
];
outputDir = ExpandFileName[args[[1]]];
dataDir = FileNameJoin[{outputDir, "data"}];
profileDir = FileNameJoin[{dataDir, "general_d_m2_scalar_threshold", "profile_tree"}];
CreateDirectory[profileDir, CreateIntermediateDirectories -> True];
reportPath = FileNameJoin[{profileDir, "profile_c5_flavour_structure_result.json"}];
hardTermsPath = FileNameJoin[{profileDir, "profile_c5_flavour_structure_hard_terms.json"}];
inheritedTermsPath = FileNameJoin[{profileDir, "profile_c5_flavour_structure_inherited_terms.json"}];
scriptDir = DirectoryName @ ExpandFileName[$InputFileName];

ClearAll[D2Text, D2Fail, D2Terms, D2HasHeavy, D2FieldQ, D2FieldDegree,
  D2FlavourNames, D2Records, D2SourceSummary];
D2Text[x_] := ToString[x, InputForm];
D2Fail[code_Integer, message_String, extra_: <||>] := (
  Print["[C5 structure] ERROR: ", message];
  Export[reportPath, Join[<|"Status" -> "Failed", "Reason" -> message,
    "ExitCode" -> code, "FixedOrderSMEFTAssembled" -> False|>, extra], "RawJSON"];
  Exit[code]
);
D2Terms[x_] := Module[{e = Expand[x]},
  If[e === 0, {}, If[Head[e] === Plus, List @@ e, {e}]]
];
D2HasHeavy[x_] := Module[{s = D2Text[x]},
  AnyTrue[{"NewScalar1", "T3S2R", "T3S2I"},
    StringContainsQ[s, "Field[" <> # <> ","] &]
];
(* Counts Field[..., ...]^n with its true multiplicity, including fields
   inside Bar[..] or NCM[..].  No dependence on string occurrence counts. *)
D2FieldQ[Matchete`Field[name_Symbol, ___], target_String] :=
  SymbolName[name] === target;
D2FieldQ[_, _String] := False;
D2FieldDegree[x_ /; D2FieldQ[x, "H"], "H"] := 1;
D2FieldDegree[x_ /; D2FieldQ[x, "l"], "l"] := 1;
D2FieldDegree[Power[x_, n_Integer?Positive], target_String] :=
  n D2FieldDegree[x, target];
D2FieldDegree[x_?AtomQ, _String] := 0;
D2FieldDegree[x_, target_String] :=
  Total[D2FieldDegree[#, target] & /@ (List @@ x)];
(* Flavour-index labels are diagnostic, not a mapping of flavour dummies
   into the external physical labels p,q. *)
D2FlavourNames[x_] := DeleteDuplicates @ StringCases[D2Text[x],
  RegularExpression["Index\\[([^,\\]\\[]+),\\s*(?:Matchete`)?Flavor\\]"] -> "$1"];
D2Records[terms_List, source_String] := MapIndexed[Function[{t, pos},
  Module[{s = D2Text[t], nH, nL, hol, hc},
    nH = D2FieldDegree[t, "H"];
    nL = D2FieldDegree[t, "l"];
    hol = StringContainsQ[s, "Proj[-1]"];
    hc = StringContainsQ[s, "Proj[1]"];
    <|
      "Index" -> First[pos], "Source" -> source,
      "HiggsLegs" -> nH, "LeptonDoubletLegs" -> nL,
      "LLHHFieldContent" -> (nH === 2 && nL === 2),
      "HasPL" -> hol, "HasPR" -> hc,
      "Chirality" -> Which[hol && !hc, "PL", hc && !hol, "PR",
        hol && hc, "Mixed", True, "Unknown"],
      "FlavourIndexLabels" -> D2FlavourNames[t],
      "HasHeavyScalar" -> D2HasHeavy[t],
      "HasGammaCC" -> !FreeQ[t, GammaCC],
      "HasHbar" -> !FreeQ[t, hbar],
      "InputForm" -> s
    |>
  ]], terms];
D2SourceSummary[records_List] := <|
  "TermCount" -> Length[records],
  "PLTermCount" -> Count[Lookup[records, "Chirality", {}], "PL"],
  "PRTermCount" -> Count[Lookup[records, "Chirality", {}], "PR"],
  "MixedOrUnknownChiralityCount" -> Count[Lookup[records, "Chirality", {}],
    "Mixed" | "Unknown"],
  "TwoHiggsTwoLeptonCount" -> Count[Lookup[records, "LLHHFieldContent", {}], True],
  "NonLLHHCount" -> Count[Lookup[records, "LLHHFieldContent", {}], False],
  "AnyHeavyScalarInOperator" -> MemberQ[Lookup[records, "HasHeavyScalar", {}], True],
  "DistinctFlavourIndexLabels" -> DeleteDuplicates @ Flatten[
    Lookup[records, "FlavourIndexLabels", {}]]
|>;

paths = <|
  "PriorReport" -> FileNameJoin[{profileDir, "profile_fixed_order_components_result.json"}],
  "HardHbar" -> FileNameJoin[{profileDir, "profile_full_scalar_one_loop_hard_hbar_sector.wxf"}],
  "FirstLoop" -> FileNameJoin[{dataDir, "eft1_after_F_realbasis_loop.wxf"}]
|>;
missing = Select[Keys[paths], !FileExistsQ[paths[#]] &];
If[missing =!= {}, D2Fail[3, "Missing saved inputs", <|"MissingKeys" -> missing|>]];
prior = Quiet @ Check[Import[paths["PriorReport"], "RawJSON"], $Failed];
If[!AssociationQ[prior] || Lookup[prior, "Status", ""] =!= "ComponentsAudited" ||
   Lookup[prior, "HardC5LinearInHbar", False] =!= True ||
   Lookup[prior, "InheritedPureSMWeinbergExtractionStatus", ""] =!= "Success",
  D2Fail[4, "Prior fixed-order component audit not certified"]
];

Print["[C5 structure] Loading Matchete and project Weinberg extractor; no matching..."];
If[CheckAbort[UsingFrontEnd[Needs["Matchete`"]; True], $Aborted] =!= True,
  D2Fail[5, "Could not load Matchete"]
];
Get[FileNameJoin[{scriptDir, "PhysicsNotation.wl"}]];
Get[FileNameJoin[{scriptDir, "RunMatching.wl"}]];
hard = Quiet @ Check[Import[paths["HardHbar"], "WXF"], $Failed];
firstLoop = Quiet @ Check[Import[paths["FirstLoop"], "WXF"], $Failed];
If[MemberQ[{hard, firstLoop}, $Failed], D2Fail[6, "Could not import saved WXF"]];

(* Probe the degree counter to avoid the previously diagnosed Field^n bug. *)
hSample = Matchete`Field[H, Matchete`Scalar, {}, {}];
lSample = Matchete`Field[l, Matchete`Fermion, {}, {}];
selfTest = {D2FieldDegree[hSample, "H"], D2FieldDegree[hSample^2, "H"],
  D2FieldDegree[Matchete`Bar[hSample]^2, "H"],
  D2FieldDegree[lSample^2, "l"]};
If[selfTest =!= {1, 2, 2, 2},
  D2Fail[7, "Field-degree self-test failed", <|"SelfTest" -> selfTest|>]
];

(* No scalar threshold matching is applied to the inherited O(hbar) action.
   Zero linear sources in L0 imply that its O(hbar) contribution is the
   pure-SM (heavy-field-zero) projection at this order. *)
inheritedTerms = Select[D2Terms[firstLoop], !D2HasHeavy[#] &];
inherited = Expand[Total[inheritedTerms]];
inheritedHbar = Expand[hbar Coefficient[inherited, hbar, 1]];
If[!TrueQ[Expand[inherited - inheritedHbar] === 0],
  D2Fail[8, "Inherited projected action has terms outside O(hbar)"]
];

hardExt = CheckAbort[ExtractWeinbergCoefficient[hard], $Aborted];
inheritedExt = CheckAbort[ExtractWeinbergCoefficient[inheritedHbar], $Aborted];
If[!AssociationQ[hardExt] || !AssociationQ[inheritedExt] ||
   Lookup[hardExt, "Status", ""] =!= "Success" ||
   Lookup[inheritedExt, "Status", ""] =!= "Success",
  D2Fail[9, "Weinberg extractor returned non-success"]
];
(* Important: these terms are from the project's existing GammaCC extractor.
   Inspect their actual light-field content before claiming they are LLHH. *)
hardRecords = D2Records[Lookup[hardExt, "Terms", {}], "HardScalar1L"];
inheritedRecords = D2Records[Lookup[inheritedExt, "Terms", {}], "InheritedFermion1L"];
hardSummary = D2SourceSummary[hardRecords];
inheritedSummary = D2SourceSummary[inheritedRecords];
If[Length[hardRecords] =!= Lookup[hardExt, "TermCount", -1] ||
   Length[inheritedRecords] =!= Lookup[inheritedExt, "TermCount", -1] ||
   Length[hardRecords] =!= Lookup[prior, "HardWeinbergTermCount", -1] ||
   Length[inheritedRecords] =!= Lookup[prior, "InheritedPureSMWeinbergTermCount", -1],
  D2Fail[10, "Weinberg inventory changed from prior audit"]
];
Export[hardTermsPath, hardRecords, "RawJSON"];
Export[inheritedTermsPath, inheritedRecords, "RawJSON"];

hardC5 = hardExt["Coefficient"];
inheritedC5 = inheritedExt["Coefficient"];
allLLHH = hardSummary["NonLLHHCount"] === 0 &&
  inheritedSummary["NonLLHHCount"] === 0;
noMixed = hardSummary["MixedOrUnknownChiralityCount"] === 0 &&
  inheritedSummary["MixedOrUnknownChiralityCount"] === 0;
result = <|
  "Status" -> If[allLLHH && noMixed,
    "WeinbergStructureAudited", "WeinbergStructureReviewRequired"],
  "NoNewMatching" -> True,
  "HardWeinberg" -> hardSummary,
  "InheritedWeinberg" -> inheritedSummary,
  "AllCandidateTermsHaveTwoLeptonsTwoHiggs" -> allLLHH,
  "AllCandidateTermsHaveUnambiguousChirality" -> noMixed,
  "HardC5HasY2" -> StringContainsQ[D2Text[hardC5], "y2,"],
  "HardC5HasY2Conj" -> StringContainsQ[D2Text[hardC5], "y2Conj,"],
  "HardC5IsLinearInHbar" ->
    TrueQ[Expand[hardC5 - hbar Coefficient[Expand[hardC5], hbar, 1]] === 0],
  "InheritedC5IsLinearInHbar" ->
    TrueQ[Expand[inheritedC5 - hbar Coefficient[Expand[inheritedC5], hbar, 1]] === 0],
  "HardOperatorRecordsFile" -> hardTermsPath,
  "InheritedOperatorRecordsFile" -> inheritedTermsPath,
  "FlavourIndexMappingCertified" -> False,
  "C5SymmetricNormalizationCertified" -> False,
  "FixedOrderSMEFTAssembled" -> False
|>;
Export[reportPath, result, "RawJSON"];
Print["[C5 structure] ", result["Status"],
  "; hard LLHH: ", hardSummary["TwoHiggsTwoLeptonCount"], "/",
  hardSummary["TermCount"], "; inherited LLHH: ",
  inheritedSummary["TwoHiggsTwoLeptonCount"], "/",
  inheritedSummary["TermCount"]];
Exit[If[allLLHH && noMixed, 0, 20]];
