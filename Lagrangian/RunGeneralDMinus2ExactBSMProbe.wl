(* D(-2) diagnostic only: exact saved EFT1 tree minus standalone SM.
   This file DOES NOT modify the authoritative general scalar threshold.
   Usage:
     wolframscript -file Lagrangian/RunGeneralDMinus2ExactBSMProbe.wl \
       STAGE1_OUTPUT_DIR DIRECT_RUNNING_INSERTION.wl MATCH_TIMEOUT_SECONDS

   Uses eft1_after_F_realbasis_tree.wxf as the exact, explicit physical EFT1
   input, rather than incrementally reconstructing selected operators.
   Every operation from Match to ReplaceEffectiveCouplings is timed.
   Successful output is NOT a physical one-loop SMEFT certification.
*)
ClearAll["Global`*"];
args = Rest[$ScriptCommandLine];
If[Length[args] =!= 3,
  Print["Usage: OUTPUT_DIR RUNNING_INSERTION.wl MATCH_TIMEOUT_SECONDS"];
  Exit[2]
];
outputDir = ExpandFileName[args[[1]]];
insertionPath = ExpandFileName[args[[2]]];
limit = Quiet @ Check[ToExpression[args[[3]], InputForm], $Failed];
If[!IntegerQ[limit] || limit < 5 || limit > 3600,
  Print["ERROR: timeout must be an integer between 5 and 3600"];
  Exit[2]
];
dataDir = FileNameJoin[{outputDir, "data"}];
profileDir = FileNameJoin[{dataDir, "general_d_m2_scalar_threshold", "profile_tree"}];
CreateDirectory[profileDir, CreateIntermediateDirectories -> True];
auditPath = FileNameJoin[{profileDir, "profile_exact_bsm_audit.json"}];
resultPath = FileNameJoin[{profileDir, "profile_exact_bsm_result.json"}];
scriptDir = DirectoryName @ ExpandFileName[$InputFileName];

ClearAll[D2Text, D2Save, D2Fail, D2Parts, D2HasHeavy];
D2Text[x_] := ToString[x, InputForm];
D2Save[path_String, content_Association] := Export[path, content, "RawJSON"];
D2Fail[code_Integer, reason_String, extras_: <||>] := (
  Print["[exact BSM] ERROR: ", reason];
  D2Save[resultPath, Join[<|
    "Status" -> "Failed", "Reason" -> reason,
    "ExitCode" -> code, "DiagnosticOnly" -> True,
    "FullGeneralScalarThresholdCertified" -> False
  |>, extras]];
  Exit[code]
);
D2Parts[expr_] := Module[{e = Expand[expr]},
  If[e === 0, {}, If[Head[e] === Plus, List @@ e, {e}]]
];
heavyNames = {"NewScalar1", "T3S2R", "T3S2I"};
D2HasHeavy[expr_] := Module[{t = D2Text[expr]},
  AnyTrue[heavyNames, StringContainsQ[t, "Field[" <> # <> ","] &]
];
If[Environment["T3_EXPERIMENTAL_GENERAL_QUADRATICS"] =!= "1" ||
   Environment["T3_EXPERIMENTAL_GENERAL_FERMIONS"] =!= "1",
  D2Fail[3, "Full-general interaction environment flags are required"]
];
If[!FileExistsQ[insertionPath],
  D2Fail[4, "Preflight direct-running file missing (never added to Match)"]
];
required = FileNameJoin[{dataDir, #}] & /@ {
  "eft1_after_F_realbasis_tree.wxf",
  "eft1_after_F_realbasis_transition_tree.wxf",
  "eft1_after_F_realbasis_wilson_tensor_seed.json",
  "eft1_after_F_realbasis_quartic_tensor_seed.json"
};
If[!And @@ (FileExistsQ /@ required),
  D2Fail[5, "Required first-threshold continuation/CG seed is missing",
    <|"Missing" -> Pick[required, Not /@ (FileExistsQ /@ required)]|>]
];

Print["[exact BSM] Loading Matchete and restoring physical field registry..."];
loaded = CheckAbort[UsingFrontEnd[Needs["Matchete`"]; True], $Aborted];
If[loaded =!= True, D2Fail[6, "Matchete did not load"]];
Get[FileNameJoin[{scriptDir, "PhysicsNotation.wl"}]];
Get[FileNameJoin[{scriptDir, "model", "T3ModelCatalog.wl"}]];
Get[FileNameJoin[{scriptDir, "model", "LagrangianBuilder.wl"}]];
Get[FileNameJoin[{scriptDir, "RunMatching.wl"}]];
model = T3ModelFromClass["D", -2];
If[!AssociationQ[model], D2Fail[7, "D(-2) catalog missing"]];
If[SetT3HeavyFields[{"F"}] === $Failed,
  D2Fail[8, "Cannot declare fermion for UV registry rebuild"]];
uvBuild = CheckAbort[UsingFrontEnd[BuildT3Lagrangian[model]], $Aborted];
If[!AssociationQ[uvBuild] ||
   Lookup[uvBuild, "Status", ""] =!= "Success" ||
   Lookup[uvBuild, "AdditionalQuarticsIncluded", -1] =!= 5 ||
   Lookup[uvBuild, "AdditionalQuadraticsIncluded", -1] =!= 1 ||
   Lookup[uvBuild, "AdditionalYukawasIncluded", -1] =!= 1,
  D2Fail[9, "Unrestricted UV registry could not be rebuilt"]
];
(* Rebuild field registry only; no fermion matching is performed. *)
RemoveField[NewFermion];
RemoveField[NewScalar1];
RemoveField[NewScalar2];
If[SetT3HeavyFields[{"S1"}] === $Failed ||
   DefineT3FieldObjects[model, False, {"S1"}] === $Failed,
  D2Fail[10, "S1 heavy-field registry restoration failed"]
];
DefineField[T3S2R, Scalar, SelfConjugate -> True,
  Mass -> {Heavy, T3S2RMass}];
DefineField[T3S2I, Scalar, SelfConjugate -> True,
  Mass -> {Heavy, T3S2IMass}];
PreserveIntegratedT3MassCouplings[{"S1"}];
If[!KeyExistsQ[GetCouplings[], T3S2Phase],
  DefineCoupling[T3S2Phase, SelfConjugate -> False]
];
fieldRegistry = GetFields[];
If[!AllTrue[{NewScalar1, T3S2R, T3S2I},
    TrueQ[Lookup[Lookup[fieldRegistry, #, <||>], Heavy, False]] &] ||
   KeyExistsQ[fieldRegistry, NewFermion] ||
   KeyExistsQ[fieldRegistry, NewScalar2],
  D2Fail[11, "Wrong physical heavy-field registry"]
];

(* Rehydrate exactly the saved physical EFT1 objects and original CG tensors. *)
{tree, transitionTree} = Quiet @ Check[Import[#, "WXF"], $Failed] & /@
  required[[1 ;; 2]];
If[MemberQ[{tree, transitionTree}, $Failed],
  D2Fail[12, "Saved physical EFT1 WXF import failed"]
];
ClearAll[D2Rehydrate];
D2Rehydrate[x_] := Quiet @ Check[
  ToExpression[D2Text[Unevaluated[x]], InputForm], $Failed
];
transitionTree = D2Rehydrate[transitionTree];
If[transitionTree === $Failed,
  D2Fail[12, "Canonical transition EFT1 rehydration failed"]
];
Scan[Function[sym,
  If[!KeyExistsQ[GetCouplings[], sym],
    Apply[DefineCoupling, {sym, SelfConjugate -> True}]
  ]
], {CNewScalar12, CT3S2R2, CT3S2I2}];
seeds = Import[#, "RawJSON"] & /@ required[[3 ;; 4]];
If[!AllTrue[seeds, AssociationQ] ||
   !AllTrue[seeds, Lookup[#, "Status", ""] === "Success" &],
  D2Fail[13, "Saved EFT1 CG seed audit invalid"]
];
registry = DeleteDuplicatesBy[
  Flatten[Lookup[seeds, "CGRegistry", {}], 1], Lookup[#, "Name", ""] &
];
restoredCG = 0;
Do[
  item = registry[[k]];
  cgName = ToExpression[item["Name"], InputForm];
  current = Quiet @ Check[
    Matchete`PackageScope`$CGproperties[cgName, Indices], $Failed];
  If[current === $Failed || MissingQ[current],
    cgReps = ToExpression[item["RepsInputForm"], InputForm];
    cgTensor = ToExpression[item["TensorInputForm"], InputForm];
    If[CheckAbort[DefineCG[cgName, cgReps, cgTensor]; True,
       $Aborted] =!= True,
      D2Fail[14, "CG restoration aborted", <|"Name" -> item["Name"]|>]
    ];
    restoredCG++
  ],
  {k, Length[registry]}
];
(* Reject any unregistered CG used by the exact saved EFT, not just the
   subset expected from the two seed files. *)
usedCGNames = DeleteDuplicates @ Cases[
  HoldComplete[tree, transitionTree],
  CG[name_, ___] :> Replace[name, Bar[x_] :> x], Infinity];
missingCGNames = Select[usedCGNames, Function[name,
  Module[{r = Quiet @ Check[
      Matchete`PackageScope`$CGproperties[name, Indices], $Failed]},
    r === $Failed || MissingQ[r]
  ]
]];
If[missingCGNames =!= {},
  D2Fail[15, "Unregistered CG tensor present in exact saved EFT1",
    <|"Names" -> (D2Text /@ missingCGNames)|>]
];
(* Same source-to-expanded consistency criterion as the existing profiler. *)
transitionExplicit = Quiet @ CheckAbort[
  ReplaceEffectiveCouplings[transitionTree], $Aborted];
transitionEqual = !MemberQ[{$Failed, $Aborted}, transitionExplicit] &&
  TrueQ[Expand[tree - Expand[transitionExplicit]] === 0];

(* Extract the standalone SM sector from the EXACT saved EFT1 tree.
   Why not subtract uvBuild["LSM"] literally? The first-stage Matchete
   canonicalisation can rename dummy indices independently in each term.
   Thus L_SM[original] and L_SM[saved EFT1] are physically the same
   expression but not necessarily Mathematica-SameQ under Expand[].

   Safety condition: BEFORE removing any pure-SM term, check that the
   complete pure-SM term multiset in saved EFT1 equals the independently
   rebuilt UV L_SM modulo capture-avoiding dummy-index alpha renaming.
   This is NOT permission to drop arbitrary light-only EFT1 operators.
   Any mismatch stops the diagnostic before Match[].

   A canonical term key replaces each distinct Index[name, type] by a
   fixed placeholder, in the order of first occurrence in InputForm.
   The complete Index object (including its representation) is the key:
   identical dummy symbols can legitimately be reused for different
   index representations. This preserves all within-term contractions.
   No actual EFT expression passed to Matchete is renamed by this audit.
*)
ClearAll[D2AlphaIndexKey, D2AlphaMultiset];
D2AlphaIndexKey[term_] := Module[{indices, rules, canon},
  indices = DeleteDuplicates @ Cases[
    HoldComplete[term], HoldPattern[Index[_, _]], Infinity];
  rules = Table[
    With[{original = indices[[i]], rep = indices[[i, 2]],
          slot = i},
      original -> Index[
        Symbol["Global`D2CanonicalDummy" <> ToString[slot]], rep]
    ], {i, Length[indices]}];
  canon = term /. rules;
  ToString[canon, InputForm]
];
D2AlphaMultiset[terms_List] := Sort[D2AlphaIndexKey /@ terms];

exactTree = Expand[tree];
allTerms = D2Parts[exactTree];
heavySource = Select[allTerms, D2HasHeavy];
sourceSM = Select[allTerms, !D2HasHeavy[#] &];
smOriginal = Expand[uvBuild["LSM"]];
smRelabelled = Expand[RelabelIndices[uvBuild["LSM"]]];
smOriginalTerms = D2Parts[smOriginal];
smRelabelledTerms = D2Parts[smRelabelled];

sourceSMKeys = D2AlphaMultiset[sourceSM];
originalSMKeys = D2AlphaMultiset[smOriginalTerms];
relabeledSMKeys = D2AlphaMultiset[smRelabelledTerms];
smAlphaOriginal = TrueQ[sourceSMKeys === originalSMKeys];
smAlphaRelabelled = TrueQ[sourceSMKeys === relabeledSMKeys];
smAlphaCertified = smAlphaOriginal || smAlphaRelabelled;
smLabel = Which[
  smAlphaOriginal, "SavedEFT1-AlphaEquivalentToOriginalLSM",
  smAlphaRelabelled, "SavedEFT1-AlphaEquivalentToRelabelledLSM",
  True, "NotCertified"
];
(* Export bounded comparison previews in case a different local source
   version makes the alpha comparison fail. Multiset comparison uses all
   term multiplicities, not mere set membership. *)
comparisonKeys = If[smAlphaOriginal, originalSMKeys, relabeledSMKeys];
smMissingKeys = Fold[DeleteCases[#1, #2, 1, 1] &,
  comparisonKeys, sourceSMKeys];
smUnexpectedKeys = Fold[DeleteCases[#1, #2, 1, 1] &,
  sourceSMKeys, comparisonKeys];

(* Use the exact saved EFT1 representatives for subtraction, but ONLY after
   the independent alpha-equivalence certificate above is available.
   Otherwise the audit will reject the input before Match[]. *)
smStandalone = Total[sourceSM];
exactBSM = Expand[exactTree - smStandalone];
reconstructionResidual = Expand[exactBSM + smStandalone - exactTree];
bsmTerms = D2Parts[exactBSM];
heavyRemainder = Select[bsmTerms, D2HasHeavy];
remainingSM = Select[bsmTerms, !D2HasHeavy[#] &];
(* Check that no operator containing any physical heavy scalar was altered. *)
heavyPreserved = SameQ[Sort[heavySource], Sort[heavyRemainder]];
extraQNames = {
  "lambdaT3ExtraQ000004I1", "lambdaT3ExtraQ000013I1",
  "lambdaT3ExtraQ001102I1", "lambdaT3ExtraQ020101I1",
  "lambdaT3ExtraQ110002I1"
};
bsmText = D2Text[exactBSM];
sourceText = D2Text[exactTree];
missingExtra = Select[extraQNames,
  !StringContainsQ[bsmText, #] &];
wilsonTerms = Select[bsmTerms,
  StringContainsQ[D2Text[#], "GammaCC"] && D2HasHeavy[#] &];
y2Terms = Select[wilsonTerms,
  StringContainsQ[D2Text[#], "Coupling[y2,"] &];
y2ConjTerms = Select[wilsonTerms,
  StringContainsQ[D2Text[#], "Coupling[y2Conj,"] &];
wilsonUnclassified = Select[wilsonTerms,
  !MemberQ[Join[y2Terms, y2ConjTerms], #] &];
wilsonOverlap = Intersection[y2Terms, y2ConjTerms];
missingScalarNames = Select[heavyNames,
  !StringContainsQ[bsmText, "Field[" <> # <> ","] &];
oldFieldsRemain = AnyTrue[{"NewFermion", "NewScalar2"},
  StringContainsQ[bsmText, "Field[" <> # <> ","] &];
inputSafe = transitionEqual && smAlphaCertified &&
  reconstructionResidual === 0 && heavyPreserved && remainingSM === {} && missingExtra === {} &&
  missingScalarNames === {} && !oldFieldsRemain &&
  Length[wilsonTerms] > 0 && Length[y2Terms] > 0 &&
  Length[y2ConjTerms] > 0 && wilsonUnclassified === {} &&
  wilsonOverlap === {};
audit = <|
  "Status" -> If[inputSafe, "Passed", "Failed"],
  "DiagnosticOnly" -> True,
  "ExactSource" -> "eft1_after_F_realbasis_tree.wxf",
  "SMSource" -> "saved EFT1 pure-SM sector; independently compared with UV LSM",
  "SMOriginalSource" -> "full general BuildT3Lagrangian[model][LSM]",
  "SMAlphaComparisonCertified" -> smAlphaCertified,
  "SavedSMAlphaEquivalentToOriginalLSM" -> smAlphaOriginal,
  "SavedSMAlphaEquivalentToRelabelledLSM" -> smAlphaRelabelled,
  "SavedEFT1PureSMTermCount" -> Length[sourceSM],
  "OriginalLSMTermCount" -> Length[smOriginalTerms],
  "RelabelledLSMTermCount" -> Length[smRelabelledTerms],
  "SMAlphaMissingTermCount" -> Length[smMissingKeys],
  "SMAlphaUnexpectedTermCount" -> Length[smUnexpectedKeys],
  "SMAlphaMissingTermPreview" -> (StringTake[#, UpTo[500]] & /@
    Take[smMissingKeys, UpTo[4]]),
  "SMAlphaUnexpectedTermPreview" -> (StringTake[#, UpTo[500]] & /@
    Take[smUnexpectedKeys, UpTo[4]]),
  "SMIndexConvention" -> smLabel,
  "SMSubtractandContainsHeavyFields" -> D2HasHeavy[smStandalone],
  "CanonicalTransitionExpansionEqualToSavedTree" -> transitionEqual,
  "ReconstructionResidualExactlyZero" -> (reconstructionResidual === 0),
  "AllHeavyTermsPreservedExactly" -> heavyPreserved,
  "StandaloneSMTermsRemaining" -> Length[remainingSM],
  "SMSubtractandExactlyEqualsSavedPureSM" ->
    TrueQ[Expand[smStandalone - Total[sourceSM]] === 0],
  "StandaloneSMResidualPreview" ->
    (StringTake[D2Text[#], UpTo[350]] & /@ Take[remainingSM, UpTo[8]]),
  "ExactTreeTermCount" -> Length[allTerms],
  "ExactBSMTermCount" -> Length[bsmTerms],
  "SourceHeavyTermCount" -> Length[heavySource],
  "BSMHeavyTermCount" -> Length[heavyRemainder],
  "WilsonTermCount" -> Length[wilsonTerms],
  "Y1Y2WilsonTermCount" -> Length[y2Terms],
  "Y1Y2ConjWilsonTermCount" -> Length[y2ConjTerms],
  "UnclassifiedWilsonTerms" -> Length[wilsonUnclassified],
  "OverlappingWilsonTerms" -> Length[wilsonOverlap],
  "MissingExtraQuartics" -> missingExtra,
  "MissingPhysicalScalars" -> missingScalarNames,
  "OldFieldsRemain" -> oldFieldsRemain,
  "CGsRestored" -> restoredCG,
  "UnregisteredCGs" -> (D2Text /@ missingCGNames),
  "ExactTreeBytes" -> ByteCount[exactTree],
  "ExactBSMBytes" -> ByteCount[exactBSM],
  "MassRelationsImposed" -> False,
  "RunningInsertedInMatch" -> False,
  "FullGeneralScalarThresholdCertified" -> False
|>;
D2Save[auditPath, audit];
Print["[exact BSM] Inventory: ", Length[allTerms], " EFT1 terms; ",
  Length[bsmTerms], " BSM terms; ", Length[wilsonTerms], " LLSS terms."];
Print["[exact BSM] SM alpha equivalence: ", smAlphaCertified,
  " (original: ", smAlphaOriginal, "; relabelled: ",
  smAlphaRelabelled, "); saved SM terms: ", Length[sourceSM]];
Print["[exact BSM] Audit status: ", audit["Status"],
  "; remaining pure-SM terms: ", Length[remainingSM]];
If[!inputSafe,
  D2Fail[20, "Exact EFT1-minus-SM audit failed; Match was NOT attempted",
    <|"AuditPath" -> auditPath|>]
];

(* Each processing step is time-bounded. A Python OS-process cap is a
   separate safeguard if a Matchete/frontend native call ignores aborts. *)
ClearAll[D2Step];
SetAttributes[D2Step, HoldRest];
steps = {};
D2Step[label_String, operation_] := Module[{started, result, status},
  Print["[exact BSM] ", label, "..."];
  started = AbsoluteTime[];
  result = TimeConstrained[
    CheckAbort[UsingFrontEnd[operation], $Aborted], limit, $TimedOut
  ];
  status = Which[
    result === $TimedOut, "TimedOut",
    result === $Aborted, "Aborted",
    result === $Failed, "Failed",
    Head[result] === Match, "UnevaluatedMatch",
    True, "Returned"
  ];
  AppendTo[steps, <|"Stage" -> label, "Status" -> status,
    "Seconds" -> AbsoluteTime[] - started,
    "Bytes" -> If[status === "Returned", ByteCount[result], 0],
    "Head" -> D2Text[Head[result]]|>];
  If[status =!= "Returned",
    D2Fail[30, "Timed-out/failed Matchete stage: " <> label,
      <|"Stages" -> steps, "AuditPath" -> auditPath|>]
  ];
  result
];
(* This relabelling uses exactly the same convention as existing profiles. *)
matchInput = RelabelIndices[exactBSM];
raw = D2Step["Match", Match[matchInput, EFTOrder -> 5, LoopOrder -> 0]];
green = D2Step["GreensSimplify", GreensSimplify[raw]];
eom = D2Step["EOMSimplify", EOMSimplify[green]];
loopEval = D2Step["EvaluateLoopFunctions", EvaluateLoopFunctions[eom]];
matched = D2Step["ReplaceEffectiveCouplings",
  ReplaceEffectiveCouplings[loopEval]];

(* Tree-level Weinberg may be zero.  Its absence is an allowed diagnostic
   result and must not be converted to a spurious one-loop C5. *)
finalText = D2Text[matched];
heavyLeft = Select[heavyNames,
  StringContainsQ[finalText, "Field[" <> # <> ","] &];
weinberg = TimeConstrained[
  Quiet @ CheckAbort[ExtractWeinbergCoefficient[matched], $Aborted],
  limit, $TimedOut];
weinbergStatus = If[AssociationQ[weinberg],
  Lookup[weinberg, "Status", "Unknown"], D2Text[weinberg]];
result = <|
  "Status" -> If[heavyLeft === {}, "CompleteTreeDiagnostic", "HeavyFieldsRemain"],
  "DiagnosticOnly" -> True,
  "FullGeneralScalarThresholdCertified" -> False,
  "OneLoopScalarHardMatchingPerformed" -> False,
  "RunningInsertedInMatch" -> False,
  "ExactInputAuditPath" -> auditPath,
  "Stages" -> steps,
  "FinalEFTByteCount" -> ByteCount[matched],
  "FinalEFTLeafCount" -> LeafCount[matched],
  "FinalEFTTermCount" -> Length[D2Parts[matched]],
  "RemainingHeavyFields" -> heavyLeft,
  "WeinbergExtractionStatus" -> weinbergStatus,
  "WeinbergProjectionVerified" -> False,
  "SMStandaloneRestored" -> False,
  "MassRelationsImposed" -> False
|>;
D2Save[resultPath, result];
Export[FileNameJoin[{profileDir, "profile_exact_bsm_matched.wxf"}],
  matched, "WXF"];
Export[FileNameJoin[{profileDir, "profile_exact_bsm_matched_inputform.txt"}],
  finalText, "Text"];
Print["[exact BSM] Tree processing: ", result["Status"],
  "; remaining heavy fields: ", heavyLeft];
Exit[If[heavyLeft === {}, 0, 32]];
