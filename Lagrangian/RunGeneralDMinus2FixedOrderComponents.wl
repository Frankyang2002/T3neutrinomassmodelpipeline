(*
 D(-2) read-only fixed-order components diagnostic.

 Inputs are the ALREADY GENERATED physical real-S2 EFT1 tree and first
 threshold loop correction, plus the PROCESSED full-EFT1 heavy-scalar
 one-loop match from RunGeneralDMinus2FullScalarLoopProbe v3.

 No new Match[] call. No changes to UV, stage-one EFT1, or threshold driver.
 Does not declare C5 fully matched: flavour symmetrization, normalization,
 scalar mass relations and SM canonicalization of the inherited sector
 remain to be checked.

 Fixed-order decomposition:
   C5_SMEFT = [C5_hard_scalar_1L] + [C5_inherited_EFT1_1L]
             + hbar*[C5_direct_running].
 All hard/inherited expressions are stored WITH their native hbar.
 The running source deliberately has NO hbar (hbar applied at assembly).

 Since L0 has no degree-one heavy-scalar terms, its classical stationary
 solution about the zero-heavy-field vacuum is S^(0)=0. Thus the O(hbar)
 inherited action contributes at tree level through its pure-SM terms.
 This diagnostic isolates those terms WITHOUT independently attempting
 EOM canonicalization of a kinetic-term-free difference.
*)
ClearAll["Global`*"];
args = Rest[$ScriptCommandLine];
If[Length[args] =!= 2,
  Print["Usage: OUTPUT_DIR RUNNING_INSERTION.wl"];
  Exit[2]
];
outputDir = ExpandFileName[args[[1]]];
insertionPath = ExpandFileName[args[[2]]];
dataDir = FileNameJoin[{outputDir, "data"}];
profileDir = FileNameJoin[{dataDir, "general_d_m2_scalar_threshold", "profile_tree"}];
CreateDirectory[profileDir, CreateIntermediateDirectories -> True];
resultPath = FileNameJoin[{profileDir, "profile_fixed_order_components_result.json"}];
scriptDir = DirectoryName @ ExpandFileName[$InputFileName];

ClearAll[D2Text, D2Save, D2Fail, D2Terms, D2ContainsHeavy,
  D2HeavyScalarQ, D2HeavyDegree, D2CoefficientStatus];
D2Text[x_] := ToString[x, InputForm];
D2Save[content_Association] := Export[resultPath, content, "RawJSON"];
D2Fail[code_Integer, reason_String, extra_: <||>] := (
  Print["[fixed-order] ERROR: ", reason];
  D2Save[Join[<|"Status" -> "Failed", "Reason" -> reason,
    "ExitCode" -> code, "FixedOrderSMEFTAssembled" -> False|>, extra]];
  Exit[code]
);
D2Terms[x_] := Module[{e = Expand[x]},
  If[e === 0, {}, If[Head[e] === Plus, List @@ e, {e}]]
];
heavyNames = {"NewScalar1", "T3S2R", "T3S2I"};
D2ContainsHeavy[x_] := Module[{s = D2Text[x]},
  AnyTrue[heavyNames, StringContainsQ[s, "Field[" <> # <> ","] &]
];
D2HeavyScalarQ[Matchete`Field[name_Symbol, Matchete`Scalar, ___]] :=
  MemberQ[heavyNames, SymbolName[name]];
D2HeavyScalarQ[_] := False;
D2HeavyDegree[x_ /; D2HeavyScalarQ[x]] := 1;
D2HeavyDegree[Power[x_, n_Integer?Positive]] := n D2HeavyDegree[x];
D2HeavyDegree[x_?AtomQ] := 0;
D2HeavyDegree[x_] := Total[D2HeavyDegree /@ (List @@ x)];

If[!FileExistsQ[insertionPath],
  D2Fail[3, "Direct running insertion file missing"]
];
required = <|
  "Tree" -> FileNameJoin[{dataDir, "eft1_after_F_realbasis_tree.wxf"}],
  "FirstLoop" -> FileNameJoin[{dataDir, "eft1_after_F_realbasis_loop.wxf"}],
  "HardProcessed" -> FileNameJoin[{profileDir,
    "profile_full_scalar_one_loop_processed.wxf"}],
  "HardHbar" -> FileNameJoin[{profileDir,
    "profile_full_scalar_one_loop_hard_hbar_sector.wxf"}],
  "HardCandidate" -> FileNameJoin[{profileDir,
    "profile_full_scalar_one_loop_hard_c5_candidate.txt"}],
  "HardResult" -> FileNameJoin[{profileDir,
    "profile_full_scalar_one_loop_result.json"}],
  "MassAudit" -> FileNameJoin[{profileDir,
    "profile_full_scalar_one_loop_audit.json"}]
|>;
missing = Select[Keys[required], !FileExistsQ[required[#]] &];
If[missing =!= {}, D2Fail[4, "Required saved files missing",
  <|"MissingKeys" -> missing|>]];

Print["[fixed-order] Loading Matchete and the existing extractor; no matching..."];
loaded = CheckAbort[UsingFrontEnd[Needs["Matchete`"]; True], $Aborted];
If[loaded =!= True, D2Fail[5, "Matchete did not load"]];
Get[FileNameJoin[{scriptDir, "PhysicsNotation.wl"}]];
Get[FileNameJoin[{scriptDir, "RunMatching.wl"}]];

massAudit = Quiet @ Check[Import[required["MassAudit"], "RawJSON"], $Failed];
hardResult = Quiet @ Check[Import[required["HardResult"], "RawJSON"], $Failed];
If[!AssociationQ[massAudit] || !AssociationQ[hardResult] ||
   Lookup[massAudit, "Status", ""] =!= "Passed" ||
   Lookup[massAudit, "MassPowerCountingPromotionOnly", False] =!= True ||
   Lookup[massAudit, "CheckLagrangianPassed", False] =!= True ||
   Lookup[hardResult, "Status", ""] =!= "FullEFT1OneLoopPostprocessed" ||
   Lookup[hardResult, "WeinbergExtractionStatus", ""] =!= "Success" ||
   Lookup[hardResult, "WeinbergCandidateFileWritten", False] =!= True,
  D2Fail[6, "Saved v3 raw match/mass/extraction certification missing"]
];

paths = Lookup[required, {"Tree", "FirstLoop", "HardProcessed", "HardHbar"}];
imported = Table[
  Quiet @ Check[Import[paths[[i]], "WXF"], $Failed], {i, Length[paths]}];
If[MemberQ[imported, $Failed], D2Fail[7, "Saved WXF import failed"]];
{tree, firstLoop, hardProcessed, hardHbarSaved} = imported;

(* Check the same multiplicity counting that was independently tested on
   the 79-term saved EFT1 tree. *)
selfField = Matchete`Field[T3S2I, Matchete`Scalar, {}, {}];
selfTest = D2HeavyDegree /@ {
  selfField, selfField^2, selfField^4, Matchete`Bar[selfField]^2};
If[selfTest =!= {1, 2, 4, 2},
  D2Fail[8, "Heavy-degree self-test failed", <|"SelfTest" -> selfTest|>]
];
treeTerms = D2Terms[tree];
treeDegrees = D2HeavyDegree /@ treeTerms;
noLinearTree = Count[treeDegrees, 1] === 0;
If[!noLinearTree || Length[treeTerms] =!= 79,
  D2Fail[9, "Saved EFT1 tree changed or has linear heavy sources",
    <|"TreeTermCount" -> Length[treeTerms],
      "HeavyDegreeDistribution" -> (List @@@ SortBy[Tally[treeDegrees], First])|>]
];

(* Verify the previously saved hard hbar piece against the complete
   processed one-loop expression. Never simply assume the file is current. *)
hardExtracted = Expand[hbar Coefficient[Expand[hardProcessed], hbar, 1]];
hardHbarEqual = TrueQ[Expand[hardHbarSaved - hardExtracted] === 0];
hardHigherOrder = !TrueQ[Coefficient[Expand[hardProcessed], hbar, 2] === 0];
noHeavyInProcessed = !D2ContainsHeavy[hardProcessed];
If[!hardHbarEqual || hardHigherOrder || !noHeavyInProcessed ||
   hardHbarSaved === 0,
  D2Fail[10, "Saved hard hbar provenance or heavy-field audit failed",
    <|"HardSavedMatchesProcessedHbar" -> hardHbarEqual,
      "HardHasExplicitHbarSquared" -> hardHigherOrder,
      "HardProcessedContainsHeavy" -> !noHeavyInProcessed|>]
];

Print["[fixed-order] Re-extracting candidate hard C5..."];
hardExt = CheckAbort[ExtractWeinbergCoefficient[hardHbarSaved], $Aborted];
If[!AssociationQ[hardExt] ||
   Lookup[hardExt, "Status", ""] =!= "Success",
  D2Fail[11, "Existing project extractor failed on saved hard hbar sector"]
];
hardC5 = hardExt["Coefficient"];
hardC5Linear = TrueQ[
  Expand[hardC5 - hbar Coefficient[Expand[hardC5], hbar, 1]] === 0
];
hardC5Text = D2Text[hardC5];
oldC5Text = StringTrim[Import[required["HardCandidate"], "Text"]];
matchingCandidateText = (StringTrim[hardC5Text] === oldC5Text);
If[!hardC5Linear || !matchingCandidateText,
  D2Fail[12, "Hard candidate has changed or has non-linear hbar structure",
    <|"HardC5LinearInHbar" -> hardC5Linear,
      "HardCandidateMatchesSavedText" -> matchingCandidateText|>]
];

(* At O(hbar), S^(0)=0 and the inherited first-threshold correction
   contributes through its S=0 projection. This is NOT a separately
   EOM-reduced SMEFT and is NOT to be one-loop matched again. *)
loopTerms = D2Terms[firstLoop];
inheritedPureSMTerms = Select[loopTerms, !D2ContainsHeavy[#] &];
inheritedWithHeavyTerms = Select[loopTerms, D2ContainsHeavy];
inheritedPureSM = Expand[Total[inheritedPureSMTerms]];
inheritedNonHbar = Coefficient[inheritedPureSM, hbar, 0];
inheritedHbarSquared = Coefficient[inheritedPureSM, hbar, 2];
inheritedHbar = Expand[hbar Coefficient[inheritedPureSM, hbar, 1]];
Export[FileNameJoin[{profileDir,
  "profile_fixed_order_components_inherited_pure_sm.wxf"}],
  inheritedPureSM, "WXF"];

Print["[fixed-order] Inspecting inherited pure-SM Weinberg sector..."];
inheritedExt = CheckAbort[ExtractWeinbergCoefficient[inheritedHbar], $Aborted];
inheritedStatus = If[AssociationQ[inheritedExt],
  Lookup[inheritedExt, "Status", "Unknown"], D2Text[inheritedExt]];
inheritedCount = If[AssociationQ[inheritedExt],
  Lookup[inheritedExt, "TermCount", 0], 0];
inheritedC5 = If[inheritedStatus === "Success",
  Lookup[inheritedExt, "Coefficient", Missing["CoefficientMissing"]],
  If[inheritedStatus === "NotFound", 0, Missing["NotCertified"]]];
If[inheritedStatus === "Success" && !MissingQ[inheritedC5],
  Export[FileNameJoin[{profileDir,
    "profile_fixed_order_components_inherited_c5_candidate.txt"}],
    D2Text[inheritedC5], "Text"]
];

(* Direct one-loop running source has no native hbar. Keep it separate,
   and do not combine it with a hard coefficient whose flavour labels
   and symmetric LLHH normalization have not been verified. *)
Print["[fixed-order] Checking separate direct running insertion..."];
Get[insertionPath];
If[DownValues[EFT1D2WilsonRunningWeinbergCoefficient] === {} ||
   !ValueQ[EFT1D2WilsonRunningHeavyInsertion] ||
   EFT1D2WilsonRunningHeavyInsertion =!= 0,
  D2Fail[13, "Direct running insertion has incompatible contract"]
];
p = Index[T3D2ExternalFlavorP, Flavor];
q = Index[T3D2ExternalFlavorQ, Flavor];
directC5 = EFT1D2WilsonRunningWeinbergCoefficient[p, q];
directContainsLog = !FreeQ[directC5, Log];
directHasHbar = !FreeQ[directC5, hbar];
If[!directContainsLog || directHasHbar,
  D2Fail[14, "Direct running has unexpected Log/hbar bookkeeping",
    <|"RunningHasLog" -> directContainsLog,
      "RunningHasExplicitHbar" -> directHasHbar|>]
];
Export[FileNameJoin[{profileDir,
  "profile_fixed_order_components_running_c5_no_hbar.txt"}],
  D2Text[directC5], "Text"];

inheritedClean = TrueQ[inheritedNonHbar === 0] &&
  TrueQ[inheritedHbarSquared === 0] &&
  MemberQ[{"Success", "NotFound"}, inheritedStatus] &&
  !MissingQ[inheritedC5];
result = <|
  "Status" -> If[inheritedClean, "ComponentsAudited", "InheritedSectorReviewRequired"],
  "DiagnosticOnly" -> True,
  "NoNewMatchPerformed" -> True,
  "TreeTermCount" -> Length[treeTerms],
  "TreeHeavyDegreeDistribution" -> (List @@@ SortBy[Tally[treeDegrees], First]),
  "TreeHasLinearHeavySource" -> !noLinearTree,
  "HardProcessedByteCount" -> ByteCount[hardProcessed],
  "SavedHardSectorMatchesProcessedO1" -> hardHbarEqual,
  "SavedHardSectorContainsHeavy" -> D2ContainsHeavy[hardHbarSaved],
  "HardC5LinearInHbar" -> hardC5Linear,
  "HardCandidateMatchesSavedText" -> matchingCandidateText,
  "HardWeinbergExtractionStatus" -> Lookup[hardExt, "Status", "Unknown"],
  "HardWeinbergTermCount" -> Lookup[hardExt, "TermCount", 0],
  "HardWeinbergHolomorphicTermCount" -> Lookup[hardExt, "HolomorphicTermCount", 0],
  "HardWeinbergConjugateTermCount" -> Lookup[hardExt, "ConjugateTermCount", 0],
  "HardC5HasY2" -> StringContainsQ[hardC5Text, "y2,"],
  "HardC5HasY2Conj" -> StringContainsQ[hardC5Text, "y2Conj,"],
  "FirstThresholdLoopTermCount" -> Length[loopTerms],
  "InheritedPureSMTermCount" -> Length[inheritedPureSMTerms],
  "InheritedHeavyTermCount" -> Length[inheritedWithHeavyTerms],
  "InheritedPureSMHasZeroHbarPart" -> TrueQ[inheritedNonHbar === 0],
  "InheritedPureSMHasHbarSquared" -> !TrueQ[inheritedHbarSquared === 0],
  "InheritedPureSMWeinbergExtractionStatus" -> inheritedStatus,
  "InheritedPureSMWeinbergTermCount" -> inheritedCount,
  "InheritedPureSMC5CandidateWritten" ->
    (inheritedStatus === "Success" && !MissingQ[inheritedC5]),
  "DirectRunningHasLog" -> directContainsLog,
  "DirectRunningHasExplicitHbar" -> directHasHbar,
  "FlavourAssemblyCertified" -> False,
  "HardC5NormalizationCertified" -> False,
  "InheritedPureSMCanonicallyReducedAtSecondThreshold" -> False,
  "PhysicalScalarMassRelationsImposed" -> False,
  "FixedOrderSMEFTAssembled" -> False,
  "FullGeneralScalarThresholdCertified" -> False,
  "ComponentConvention" ->
    "Hard and inherited already contain native hbar; running has no hbar: at O(hbar) C5=hard+inherited+hbar*running, after flavour/normalization checks"
|>;
D2Save[result];
Print["[fixed-order] ", result["Status"],
  "; hard Weinberg terms: ", result["HardWeinbergTermCount"],
  "; inherited pure-SM: ", Length[inheritedPureSMTerms],
  "; inherited Weinberg: ", inheritedStatus,
  "; no new Match[] called."];
Exit[If[inheritedClean, 0, 20]];
