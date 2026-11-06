(* D(-2) FULL PHYSICAL EFT1 scalar one-loop first-pass, diagnostic only.
   Rebuilds the general UV CG/coupling registry, registers the three
   physical heavy scalars, imports the saved physical EFT1 tree, proves
   the SM/BSM partition and the two independent LLSS Yukawa sectors,
   audits the heavy quadratic mass/kinetic terms against FreeLag,
   promoting the three physical scalar masses from the light-mass EFT
   counting used at the fermion threshold to heavy-scale EFT counting.
   Calls CheckLagrangian and then raw Match[full EFT1,LoopOrder->1]
   ONLY if those gates pass. After the raw match it performs the
   established Greens/EOM/loop/effective-coupling reduction chain,
   reports a tentative hard one-loop Weinberg coefficient, and exports
   intermediate WXF data. Does NOT assemble physical SMEFT C5.
   No saved stage-one WXF or authoritative scalar driver is changed.
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
auditPath = FileNameJoin[{profileDir, "profile_full_scalar_one_loop_audit.json"}];
resultPath = FileNameJoin[{profileDir, "profile_full_scalar_one_loop_result.json"}];
scriptDir = DirectoryName @ ExpandFileName[$InputFileName];

ClearAll[D2Text, D2Save, D2Fail, D2Parts, D2HasHeavy,
  D2ElapsedSeconds];
D2Text[x_] := ToString[x, InputForm];
(* Quantize to milliseconds to avoid RawJSON scientific zero tokens. *)
D2ElapsedSeconds[seconds_?NumericQ] := Module[
  {milliseconds = Max[0, Round[1000 seconds]]},
  If[milliseconds === 0, 0, N[milliseconds/1000]]
];
D2Save[path_String, content_Association] := Export[path, content, "RawJSON"];
D2Fail[code_Integer, reason_String, extras_: <||>] := (
  Print["[full scalar 1L] ERROR: ", reason];
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

Print["[full scalar 1L] Loading Matchete and restoring physical field registry..."];
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
Print["[full scalar 1L] Inventory: ", Length[allTerms], " EFT1 terms; ",
  Length[bsmTerms], " BSM terms; ", Length[wilsonTerms], " LLSS terms."];
Print["[full scalar 1L] SM alpha equivalence: ", smAlphaCertified,
  " (original: ", smAlphaOriginal, "; relabelled: ",
  smAlphaRelabelled, "); saved SM terms: ", Length[sourceSM]];
Print["[full scalar 1L] Audit status: ", audit["Status"],
  "; remaining pure-SM terms: ", Length[remainingSM]];
If[!inputSafe,
  D2Fail[20, "Exact EFT1-minus-SM audit failed; Match was NOT attempted",
    <|"AuditPath" -> auditPath|>]
];

(* Registration and mass-term audit, without changing the saved EFT1.
   At the F threshold S1, S2R and S2I were declared LIGHT, so their mass
   parameters were serialized as Coupling[M, {}, 1] (EFTOrder 1).
   Here all three scalar fields are HEAVY and FreeLag generates
   Coupling[M, {}, 0]. The final argument is Matchete EFT power counting,
   NOT a loop order. Promote ONLY the three exact mass labels in a copy
   of the saved EFT1; preserve all numerical mass parameters and the
   original stage-one WXF unchanged.
   Both quadratic and interaction occurrences of these mass labels must
   be promoted to avoid mixed mass counting inside the input theory. *)
ClearAll[D2PromoteScalarMasses, D2DemoteScalarMasses];
D2PromoteScalarMasses[expr_] := expr /. {
  HoldPattern[Matchete`Coupling[MS1, {}, 1]] :>
    Matchete`Coupling[MS1, {}, 0],
  HoldPattern[Matchete`Coupling[T3S2RMass, {}, 1]] :>
    Matchete`Coupling[T3S2RMass, {}, 0],
  HoldPattern[Matchete`Coupling[T3S2IMass, {}, 1]] :>
    Matchete`Coupling[T3S2IMass, {}, 0]
};
D2DemoteScalarMasses[expr_] := expr /. {
  HoldPattern[Matchete`Coupling[MS1, {}, 0]] :>
    Matchete`Coupling[MS1, {}, 1],
  HoldPattern[Matchete`Coupling[T3S2RMass, {}, 0]] :>
    Matchete`Coupling[T3S2RMass, {}, 1],
  HoldPattern[Matchete`Coupling[T3S2IMass, {}, 0]] :>
    Matchete`Coupling[T3S2IMass, {}, 1]
};
massSymbols = {MS1, T3S2RMass, T3S2IMass};
sourceLightMassCounts = Table[
  With[{m = massSymbols[[i]]},
    Count[exactTree,
      HoldPattern[Matchete`Coupling[m, {}, 1]], Infinity]],
  {i, 3}];
matchInput = D2PromoteScalarMasses[exactTree];
(* The saved tree is fully explicit. Match it directly rather than rely
   on possible mass-order metadata in hidden effective coupling definitions
   in the compact transition-tree representation. *)
promotedHeavyMassCounts = Table[
  With[{m = massSymbols[[i]]},
    Count[matchInput,
      HoldPattern[Matchete`Coupling[m, {}, 0]], Infinity]],
  {i, 3}];
staleLightMassCounts = Table[
  With[{m = massSymbols[[i]]},
    Count[matchInput,
      HoldPattern[Matchete`Coupling[m, {}, 1]], Infinity]],
  {i, 3}];
exactReversiblePromotion = TrueQ[
  D2DemoteScalarMasses[matchInput] === exactTree];
promotionSafe = AllTrue[sourceLightMassCounts, # >= 1 &] &&
  sourceLightMassCounts === promotedHeavyMassCounts &&
  staleLightMassCounts === {0, 0, 0} && exactReversiblePromotion;
If[!promotionSafe,
  D2Fail[25, "Mass EFT-order promotion failed strict preservation audit",
    <|"SourceLightMassCounts" -> sourceLightMassCounts,
      "PromotedHeavyMassCounts" -> promotedHeavyMassCounts,
      "StaleLightMassCounts" -> staleLightMassCounts,
      "ExactReversiblePromotion" -> exactReversiblePromotion|>]
];
ClearAll[D2HeavyScalarQ, D2HeavyDegree, D2AllFieldNames];
D2HeavyScalarQ[Matchete`Field[name_Symbol, Matchete`Scalar, ___]] :=
  MemberQ[heavyNames, SymbolName[name]];
D2HeavyScalarQ[_] := False;
D2HeavyDegree[x_ /; D2HeavyScalarQ[x]] := 1;
D2HeavyDegree[Power[x_, n_Integer?Positive]] := n D2HeavyDegree[x];
D2HeavyDegree[x_?AtomQ] := 0;
D2HeavyDegree[x_] := Total[D2HeavyDegree /@ (List @@ x)];
D2AllFieldNames[term_] := DeleteDuplicates @ Cases[
  HoldComplete[term],
  Matchete`Field[name_Symbol, _, ___] :> SymbolName[name],
  Infinity];

selfField = Matchete`Field[T3S2I, Matchete`Scalar, {}, {}];
selfTest = D2HeavyDegree /@ {
  selfField, selfField^2, selfField^4, Matchete`Bar[selfField]^2};
If[selfTest =!= {1, 2, 4, 2},
  D2Fail[24, "Heavy field degree self-test failed",
    <|"SelfTest" -> selfTest|>]
];

registered = GetFields[];
heavies = {NewScalar1, T3S2R, T3S2I};
massLabels = {"MS1", "T3S2RMass", "T3S2IMass"};
heavyFlags = Table[
  TrueQ[Lookup[Lookup[registered, heavies[[i]], <||>], Heavy, False]],
  {i, 3}];
sourceTerms = D2Parts[exactTree];
degrees = D2HeavyDegree /@ sourceTerms;
pureHeavyQuadratics = Pick[sourceTerms,
  MapThread[And, {(# === 2 & /@ degrees),
    (AllTrue[D2AllFieldNames[#], MemberQ[heavyNames, #] &] & /@ sourceTerms)}],
  True];

(* Pure quadratic terms should match exactly the canonically registered
   FreeLag, up to dummy-index alpha renaming.  Preserve both mismatch
   previews, rather than changing the mass basis or guessing a fix. *)
free = CheckAbort[
  UsingFrontEnd[RelabelIndices[Expand[FreeLag[NewScalar1, T3S2R, T3S2I]]]],
  $Aborted];
If[free === $Aborted || free === $Failed,
  D2Fail[25, "FreeLag for the registered heavy scalars failed"]
];
freeTerms = D2Parts[free];
freePureQuadratics = Select[freeTerms,
  D2HeavyDegree[#] === 2 &&
    AllTrue[D2AllFieldNames[#], MemberQ[heavyNames, #] &] &];
originalQuadraticKeys = D2AlphaMultiset[pureHeavyQuadratics];
promotedQuadratics = D2PromoteScalarMasses /@ pureHeavyQuadratics;
sourceKeys = D2AlphaMultiset[promotedQuadratics];
freeKeys = D2AlphaMultiset[freePureQuadratics];
originalQuadraticAlphaEqual = originalQuadraticKeys === freeKeys;
quadraticAlphaEqual = sourceKeys === freeKeys;
missingKeys = Complement[freeKeys, sourceKeys];
unexpectedKeys = Complement[sourceKeys, freeKeys];
massCounts = Table[
  Count[pureHeavyQuadratics, term_ /;
    StringContainsQ[D2Text[term], "Coupling[" <> massLabels[[i]] <> ","] &&
    StringContainsQ[D2Text[term], "Field[" <> heavyNames[[i]] <> ","]],
  {i, 3}];
freeMassCounts = Table[
  Count[freePureQuadratics, term_ /;
    StringContainsQ[D2Text[term], "Coupling[" <> massLabels[[i]] <> ","] &&
    StringContainsQ[D2Text[term], "Field[" <> heavyNames[[i]] <> ","]],
  {i, 3}];

(* Quadratic alpha comparison is useful but not sufficient for the
   literal complete-EFT1 matching. Fail closed on a mismatch. *)
massSafe = promotionSafe && And @@ heavyFlags && quadraticAlphaEqual &&
  AllTrue[massCounts, # > 0 &] && AllTrue[freeMassCounts, # > 0 &] &&
  Count[degrees, 1] === 0 && Count[degrees, 3] === 0;
massAudit = Join[audit, <|
  "Status" -> If[massSafe, "Passed", "Failed"],
  "Mode" -> "FullSavedEFT1_OneLoopMassAudit",
  "FullPhysicalSMIncluded" -> True,
  "AllPhysicalScalarsRegisteredHeavy" -> (And @@ heavyFlags),
  "HeavyRegistrationFlags" -> AssociationThread[heavyNames, heavyFlags],
  "RegistryRecordsInputForm" -> AssociationThread[heavyNames,
    (D2Text[Lookup[registered, #, Missing["NotRegistered"]]] & /@ heavies)],
  "MassLabels" -> AssociationThread[heavyNames, massLabels],
  "MassPowerCountingPromotionOnly" -> promotionSafe,
  "SavedLightMassEFTOrder" -> 1,
  "ThresholdHeavyMassEFTOrder" -> 0,
  "SourceLightMassCounts" -> AssociationThread[heavyNames, sourceLightMassCounts],
  "PromotedHeavyMassCounts" -> AssociationThread[heavyNames, promotedHeavyMassCounts],
  "StaleLightMassCounts" -> AssociationThread[heavyNames, staleLightMassCounts],
  "ExactReversibleMassPromotion" -> exactReversiblePromotion,
  "OriginalQuadraticAlphaEqualToFreeLag" -> originalQuadraticAlphaEqual,
  "MassTermCounts" -> AssociationThread[heavyNames, massCounts],
  "FreeLagMassTermCounts" -> AssociationThread[heavyNames, freeMassCounts],
  "HeavyDegreeDistribution" -> (List @@@ SortBy[Tally[degrees], First]),
  "PureHeavyQuadraticTermCount" -> Length[pureHeavyQuadratics],
  "FreeLagPureQuadraticTermCount" -> Length[freePureQuadratics],
  "PureHeavyQuadraticAlphaEqualToFreeLag" -> quadraticAlphaEqual,
  "FreeLagMissingTermCount" -> Length[missingKeys],
  "FreeLagUnexpectedTermCount" -> Length[unexpectedKeys],
  "FreeLagMissingPreview" -> (StringTake[#, UpTo[450]] & /@ Take[missingKeys, UpTo[5]]),
  "FreeLagUnexpectedPreview" -> (StringTake[#, UpTo[450]] & /@ Take[unexpectedKeys, UpTo[5]]),
  "NoLinearHeavyScalarSources" -> (Count[degrees, 1] === 0),
  "PhysicalMassRelationsImposed" -> False,
  "OneLoopMatchingVerified" -> False
|>];
D2Save[auditPath, massAudit];
Print["[full scalar 1L] Heavy degree distribution: ",
  massAudit["HeavyDegreeDistribution"]];
Print["[full scalar 1L] Registered heavy: ", heavyFlags,
  "; mass terms: ", massCounts,
  "; original quadratic equality: ", originalQuadraticAlphaEqual,
  "; promoted quadratic equality: ", quadraticAlphaEqual,
  "; mass EFT-order promotion only: ", promotionSafe];
If[!massSafe,
  D2Fail[26, "Promoted heavy quadratic sector differs from registered FreeLag; no one-loop Match attempted",
    <|"AuditPath" -> auditPath|>]
];

(* Match the COMPLETE saved EFT1, with ONLY the three heavy-mass
   EFT-order labels promoted to the new threshold.  Nothing is subtracted
   and no inherited hbar or direct running is inserted here. *)
Print["[full scalar 1L] Full EFT1 CheckLagrangian..."];
checkStarted = AbsoluteTime[];
check = TimeConstrained[
  CheckAbort[UsingFrontEnd[CheckLagrangian[matchInput]], $Aborted],
  limit, $TimedOut];
checkElapsed = D2ElapsedSeconds[AbsoluteTime[] - checkStarted];
massAudit = Join[massAudit, <|
  "CheckLagrangianResult" -> D2Text[check],
  "CheckLagrangianSeconds" -> checkElapsed,
  "CheckLagrangianPassed" -> TrueQ[check]|>];
D2Save[auditPath, massAudit];
If[!TrueQ[check],
  D2Fail[27, "Full EFT1 CheckLagrangian failed; Match not attempted",
    <|"AuditPath" -> auditPath, "CheckResult" -> D2Text[check]|>]
];

Print["[full scalar 1L] Starting full physical EFT1 raw scalar Match, LoopOrder->1..."];
matchStarted = AbsoluteTime[];
matched = TimeConstrained[
  CheckAbort[UsingFrontEnd[
    Match[matchInput, EFTOrder -> 5, LoopOrder -> 1]], $Aborted],
  limit, $TimedOut];
matchElapsed = D2ElapsedSeconds[AbsoluteTime[] - matchStarted];
matchedHead = D2Text[Head[matched]];
matchedValid = !MemberQ[{$Failed, $Aborted, $TimedOut}, matched] &&
  !AtomQ[matched] &&
  (StringContainsQ[D2Text[matched], "Field["] ||
   StringContainsQ[D2Text[matched], "FieldStrength["]);
remainingHeavy = If[matchedValid,
  Select[heavyNames, StringContainsQ[D2Text[matched], "Field[" <> # <> ","] &],
  heavyNames];
rawResult = <|
  "Status" -> If[matchedValid && remainingHeavy === {} && !FreeQ[matched, hbar],
    "FullEFT1RawOneLoopMatchReturned", "Failed"],
  "DiagnosticOnly" -> True,
  "FullGeneralScalarThresholdCertified" -> False,
  "FullPhysicalSMIncluded" -> True,
  "MassPowerCountingPromotionOnly" -> promotionSafe,
  "OneLoopRawMatchAttempted" -> True,
  "OneLoopRawMatchReturned" -> matchedValid,
  "LoopOrder" -> 1, "EFTOrder" -> 5,
  "MatchSeconds" -> matchElapsed,
  "MatchResultHead" -> matchedHead,
  "MatchResultBytes" -> If[matchedValid, ByteCount[matched], 0],
  "MatchResultContainsHbar" -> If[matchedValid, !FreeQ[matched, hbar], False],
  "RemainingHeavyFields" -> remainingHeavy,
  "MassAuditPath" -> auditPath,
  "SimplificationPerformed" -> False,
  "HardC5Extracted" -> False,
  "InheritedOneLoopAdded" -> False,
  "DirectRunningAdded" -> False,
  "PhysicalMassRelationsImposed" -> False
|>;
D2Save[resultPath, rawResult];
If[!matchedValid || remainingHeavy =!= {} || FreeQ[matched, hbar],
  D2Fail[30, "Full physical EFT1 one-loop Match invalid, no hbar, timed out, or left heavy fields",
    <|"RawMatchResult" -> rawResult|>]
];
Export[FileNameJoin[{profileDir, "profile_full_scalar_one_loop_raw.wxf"}],
  matched, "WXF"];
Print["[full scalar 1L] Full physical EFT1 raw one-loop Match returned; bytes: ",
  ByteCount[matched], "; hbar: ", !FreeQ[matched, hbar]];


(* V3: process the ALREADY VALIDATED raw one-loop EFT in the same kernel.
   The input remains complete EFT1[L0], with ONLY heavy-mass EFT-order
   promotion. The independently inherited EFT1[L1] and direct running
   are never added here, avoiding double counting at this step.

   This is an exploratory hard one-loop threshold. In particular,
   ExtractWeinbergCoefficient is the project's existing GammaCC-based
   extractor; its normalization and flavour mapping are NOT certified. *)
ClearAll[D2ProcessStage];
processReports = {};
currentEFT = matched;
processStages = {
  "GreensSimplify", "EOMSimplify", "EvaluateLoopFunctions",
  "ReplaceEffectiveCouplings"
};
processWxfNames = {
  "profile_full_scalar_one_loop_green.wxf",
  "profile_full_scalar_one_loop_eom.wxf",
  "profile_full_scalar_one_loop_evaluated.wxf",
  "profile_full_scalar_one_loop_processed.wxf"
};

Do[
  stageName = processStages[[i]];
  Print["[full scalar 1L] ", stageName, "..."];
  stageStart = AbsoluteTime[];
  stageValue = TimeConstrained[
    CheckAbort[
      UsingFrontEnd[Switch[stageName,
        "GreensSimplify", GreensSimplify[currentEFT],
        "EOMSimplify", EOMSimplify[currentEFT],
        "EvaluateLoopFunctions", EvaluateLoopFunctions[currentEFT],
        "ReplaceEffectiveCouplings", ReplaceEffectiveCouplings[currentEFT],
        _, $Failed
      ]],
      $Aborted
    ],
    limit, $TimedOut
  ];
  stageElapsed = D2ElapsedSeconds[AbsoluteTime[] - stageStart];
  stageGood = !MemberQ[{$Failed, $Aborted, $TimedOut}, stageValue] &&
    !AtomQ[stageValue] &&
    (StringContainsQ[D2Text[stageValue], "Field["] ||
     StringContainsQ[D2Text[stageValue], "FieldStrength["]);
  AppendTo[processReports, <|
    "Stage" -> stageName,
    "Status" -> If[stageGood, "Returned", "FailedOrFieldless"],
    "Seconds" -> stageElapsed,
    "Head" -> D2Text[Head[stageValue]],
    "Bytes" -> If[stageGood, ByteCount[stageValue], 0]
  |>];
  D2Save[resultPath, Join[rawResult, <|
    "Status" -> "FullEFT1OneLoopPostprocessingInProgress",
    "ProcessingStages" -> processReports,
    "SimplificationPerformed" -> False,
    "HardC5Extracted" -> False
  |>]];
  If[!stageGood,
    D2Fail[31, "One-loop postprocessing stage invalid/timed out: " <> stageName,
      <|"RawMatchResult" -> rawResult,
        "ProcessingStages" -> processReports|>]
  ];
  Export[FileNameJoin[{profileDir, processWxfNames[[i]]}],
    stageValue, "WXF"];
  currentEFT = stageValue;
, {i, Length[processStages]}];

(* This is the coefficient of explicit hbar in the full matched EFT,
   NOT the complete fixed-order SMEFT result. *)
processed = currentEFT;
If[FreeQ[processed, hbar],
  D2Fail[32, "Processed one-loop EFT lost all explicit hbar terms",
    <|"RawMatchResult" -> rawResult,
      "ProcessingStages" -> processReports|>]
];
loopPiece = Expand[hbar Coefficient[Expand[processed], hbar, 1]];
If[loopPiece === 0,
  D2Fail[32, "Processed one-loop EFT has no isolated O(hbar) sector",
    <|"RawMatchResult" -> rawResult,
      "ProcessingStages" -> processReports|>]
];
remainingHeavyAfter = Select[heavyNames,
  StringContainsQ[D2Text[processed], "Field[" <> # <> ","] &];
If[remainingHeavyAfter =!= {},
  D2Fail[32, "Processed EFT still contains heavy scalar fields",
    <|"RemainingHeavyFields" -> remainingHeavyAfter,
      "ProcessingStages" -> processReports|>]
];
Export[FileNameJoin[{profileDir,
  "profile_full_scalar_one_loop_hard_hbar_sector.wxf"}],
  loopPiece, "WXF"];

Print["[full scalar 1L] Extracting Weinberg sector from hard O(hbar) EFT..."];
extractStart = AbsoluteTime[];
extraction = TimeConstrained[
  CheckAbort[ExtractWeinbergCoefficient[loopPiece], $Aborted],
  limit, $TimedOut
];
extractSeconds = D2ElapsedSeconds[AbsoluteTime[] - extractStart];
extractOK = AssociationQ[extraction];
extractStatus = If[extractOK,
  Lookup[extraction, "Status", "Unknown"], D2Text[extraction]];
extractCount = If[extractOK,
  Lookup[extraction, "TermCount", 0], 0];
extractCoefficient = If[extractOK,
  Lookup[extraction, "Coefficient", Missing["NotFound"]],
  Missing["Failed"]];
If[extractOK && extractStatus === "Success" &&
   !MissingQ[extractCoefficient],
  Export[FileNameJoin[{profileDir,
    "profile_full_scalar_one_loop_hard_c5_candidate.txt"}],
    D2Text[extractCoefficient], "Text"]
];
D2Save[resultPath, Join[rawResult, <|
  "Status" -> "FullEFT1OneLoopPostprocessed",
  "ProcessingStages" -> processReports,
  "SimplificationPerformed" -> True,
  "HardHbarSectorExtracted" -> True,
  "ProcessedEFTBytes" -> ByteCount[processed],
  "ProcessedEFTHasHbar" -> !FreeQ[processed, hbar],
  "RemainingHeavyFields" -> remainingHeavyAfter,
  "WeinbergExtractionStatus" -> extractStatus,
  "WeinbergExtractionSeconds" -> extractSeconds,
  "WeinbergTermCount" -> extractCount,
  "WeinbergCandidateFileWritten" ->
    (extractOK && extractStatus === "Success" &&
     !MissingQ[extractCoefficient]),
  "HardC5Extracted" ->
    (extractOK && extractStatus === "Success" &&
     !MissingQ[extractCoefficient]),
  "FullGeneralScalarThresholdCertified" -> False,
  "FullFlavorHardC5Verified" -> False,
  "InheritedOneLoopAdded" -> False,
  "DirectRunningAdded" -> False,
  "FixedOrderSMEFTAssembled" -> False
|>]];
Print["[full scalar 1L] Full EFT1 hard one-loop EFT postprocessed; ",
  "Weinberg extraction: ", extractStatus,
  "; candidate terms: ", extractCount,
  ". No full C5 assembly has been performed."];
Exit[0];
