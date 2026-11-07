(* Batch 5: independent full-EFT1 scalar-only hard one-loop matching for six ordinary T3 models.
   Consume the SAVED native F-only tree; do not integrate out F again.
   Since Matchete rejects the triple-contracted diagonal NFlavor index in
   y1[i,r] y2[j,r] / MF[r], match a per-mass-eigenvalue *template*: keep the
   exact Yukawa tensors and external flavour labels but use a dummy, unindexed
   MF for the scalar one-loop kernel. Numerical full-flavour assembly MUST
   evaluate this template separately at each heavy mass and sum over r.
   Source: output/general_observable/fermion_first_all16_v1/<slug>/data
   This stage does not add inherited one-loop operators or EFT1 running;
   those remain separate components of the final fixed-order C5.
   This file is a new generalization, NOT yet validated with Matchete in Wolfram.
*)
ClearAll["Global`*"];
args = Rest[$ScriptCommandLine];
If[Length[args] =!= 5,
  Print["Usage: MODEL_CLASS ALPHA F_BOUNDARY_DIR OUTPUT_DIR TIMEOUT_SECONDS"];
  Exit[2]
];
modelClass = args[[1]];
(* Avoid a leading '-' argv token, which wolframscript treats as an option.
   Python passes m4/p2, not -4/+2. Parse digits only. *)
alphaToken = args[[2]];
alpha = Which[
  StringMatchQ[alphaToken, "m" ~~ DigitCharacter ..],
    -ToExpression[StringDrop[alphaToken, 1], InputForm],
  StringMatchQ[alphaToken, "p" ~~ DigitCharacter ..],
     ToExpression[StringDrop[alphaToken, 1], InputForm],
  True, $Failed
];
sourceDir = ExpandFileName[args[[3]]];
outputDir = ExpandFileName[args[[4]]];
limit = Quiet @ Check[ToExpression[args[[5]], InputForm], $Failed];
If[!MemberQ[{"A", "B", "C", "D", "E"}, modelClass] ||
   !IntegerQ[alpha] || !IntegerQ[limit] || limit < 30 || limit > 7200,
  Print["Invalid model/timeout"]; Exit[2]
];
If[!DirectoryQ[outputDir],
  CreateDirectory[outputDir, CreateIntermediateDirectories -> True]];
reportPath = FileNameJoin[{outputDir, "hard_scalar_report.json"}];
auditPath = FileNameJoin[{outputDir, "hard_scalar_input_audit.json"}];
sourceFiles = FileNameJoin[{sourceDir, #}] & /@ {
  "eft1_after_F_tree_inputform.txt",
  "eft1_after_F_wilson_seed.json",
  "eft1_after_F_scalar_quartic_seed.json"
};
ClearAll[T3Text, T3Save, T3Fail, T3Parts, T3HasHeavy,
  T3AllFields, T3Degree, T3AlphaKey];
T3Text[x_] := ToString[x, InputForm];
T3Save[path_, a_Association] := Export[path, a, "RawJSON"];
T3Fail[code_Integer, reason_String, extras_: <||>] := (
  T3Save[reportPath, Join[<|
    "Status" -> "Failed", "Reason" -> reason, "ExitCode" -> code,
    "HardOneLoopMatched" -> False, "FullSequentialC5Certified" -> False
  |>, extras]];
  Print["[ordinary-hard] ERROR: ", reason]; Exit[code]
);
T3Parts[expr_] := Module[{x = Expand[expr]},
  If[x === 0, {}, If[Head[x] === Plus, List @@ x, {x}]]
];
T3HasHeavy[expr_] := AnyTrue[{"NewScalar1", "NewScalar2"},
  StringContainsQ[T3Text[expr], "Field[" <> # <> ","] &];
T3AllFields[expr_] := DeleteDuplicates @ Cases[HoldComplete[expr],
  Matchete`Field[n_Symbol, _, ___] :> SymbolName[n], Infinity];
T3Degree[Matchete`Field[n_Symbol, Matchete`Scalar, ___]] /;
  MemberQ[{NewScalar1, NewScalar2}, n] := 1;
T3Degree[Power[x_, k_Integer?Positive]] := k T3Degree[x];
T3Degree[x_?AtomQ] := 0;
T3Degree[x_] := Total[T3Degree /@ (List @@ x)];
T3AlphaKey[term_] := Module[{ids, rules},
  ids = DeleteDuplicates @ Cases[HoldComplete[term], Index[_, _], Infinity];
  rules = Table[With[{rep = ids[[i, 2]], slot = i},
    ids[[i]] -> Index[Symbol["Global`T3GenericDummy" <> ToString[slot]], rep]],
    {i, Length[ids]}];
  T3Text[term /. rules]
];
If[Environment["T3_EXPERIMENTAL_GENERAL_QUADRATICS"] =!= "1" ||
   Environment["T3_EXPERIMENTAL_GENERAL_FERMIONS"] =!= "1",
  T3Fail[3, "Full-general fermion/quadratic environment flags are required"]
];
If[!And @@ (FileExistsQ /@ sourceFiles),
  T3Fail[4, "Missing saved fermion-first EFT1 tree or tensor seeds"]
];
scriptDir = DirectoryName @ ExpandFileName[$InputFileName];
loaded = CheckAbort[UsingFrontEnd[Needs["Matchete`"] ; True], $Aborted];
If[loaded =!= True, T3Fail[5, "Matchete failed to load"]];
Get[FileNameJoin[{scriptDir, "PhysicsNotation.wl"}]];
Get[FileNameJoin[{scriptDir, "model", "T3ModelCatalog.wl"}]];
Get[FileNameJoin[{scriptDir, "model", "LagrangianBuilder.wl"}]];
Get[FileNameJoin[{scriptDir, "RunMatching.wl"}]];
model = T3ModelFromClass[modelClass, alpha, 3];
If[!AssociationQ[model], T3Fail[6, "Unknown model in T3 catalogue"]];
label = modelClass <> "_" <> If[alpha >= 0, "+", ""] <> ToString[alpha];
(* This first physics backend handles the ordinary mass class ONLY. *)
If[MemberQ[{"A_-2", "A_+0", "D_-2", "D_+0", "E_-2", "E_+0",
    "B_-1", "C_-1"}, label],
  T3Fail[7, "Special scalar mass model requires physical mass-basis matching"]
];
If[SetT3HeavyFields[{"F"}] === $Failed,
  T3Fail[8, "Cannot reconstruct the full-general F-only UV registry"]
];
built = CheckAbort[UsingFrontEnd[BuildT3Lagrangian[model]], $Aborted];
If[!AssociationQ[built] || Lookup[built, "Status", ""] =!= "Success" ||
   !TrueQ[Lookup[built, "AllScalarExtensionsRequested", False]] ||
   Lookup[built, "AdditionalQuadraticsIncluded", -1] =!= 0 ||
   Lookup[built, "AdditionalYukawasIncluded", -1] =!= 0 ||
   Lookup[built, "AdditionalQuarticsIncluded", -1] =!=
     Lookup[built, "AdditionalQuarticsExpected", -2] ||
   Lookup[built, "AdditionalQuarticsIncluded", -1] =!=
     If[MemberQ[{"B", "C"}, modelClass], 1, 0],
  T3Fail[9, "Full-general ordinary inventory incomplete or unexpected: expected 1 additional B/C quartic, 0 E quartics, and no extra mass/Yukawa operators"]
];
(* Use source-derived actual first-threshold EFT1, not a new symbolic estimate. *)
treeString = Import[sourceFiles[[1]], "Text"];
If[!StringQ[treeString] || StringLength[treeString] < 1000 ||
    StringContainsQ[treeString, "Field[NewFermion,"] ||
    !StringContainsQ[treeString, "Field[NewScalar1,"] ||
    !StringContainsQ[treeString, "Field[NewScalar2,"],
  T3Fail[10, "InputForm does not encode an F-free EFT1 containing both heavy scalars"]
];
tree = Quiet @ Check[ToExpression[treeString, InputForm], $Failed];
If[tree === $Failed || AtomQ[tree], T3Fail[11, "Invalid saved tree InputForm"]];
seeds = Quiet @ Check[Import[#, "RawJSON"] & /@ sourceFiles[[2 ;; 3]], $Failed];
If[!ListQ[seeds] || !AllTrue[seeds, AssociationQ] ||
   Lookup[seeds[[1]], "IntegratedFields", {}] =!= {"F"} ||
   Lookup[seeds[[1]], "ActiveHeavyFields", {}] =!= {"S1", "S2"} ||
   Lookup[seeds[[1]], "TreeWilsonTermCount", 0] < 2,
  T3Fail[12, "Saved native F-only LLSS seed incomplete"]
];
(* Full general Lagrangian re-creates the native CG registry; additionally
   import exactly the CG tensors used by first-threshold EFT1. *)
registry = DeleteDuplicatesBy[
  Flatten[Lookup[seeds, "CGRegistry", {}], 1], Lookup[#, "Name", ""] &];
restored = 0;
Do[
  item = registry[[k]];
  name = ToExpression[item["Name"], InputForm];
  props = Quiet @ Check[Matchete`PackageScope`$CGproperties[name, Indices], $Failed];
  If[props === $Failed || MissingQ[props],
    reps = ToExpression[item["RepsInputForm"], InputForm];
    tensor = ToExpression[item["TensorInputForm"], InputForm];
    result = CheckAbort[DefineCG[name, reps, tensor]; True, $Aborted];
    If[result =!= True, T3Fail[13, "Native CG restoration aborted"]];
    restored++
  ],
  {k, Length[registry]}
];
(* Important: retain all exact EFT1 operators and the pure SM. *)
If[ReclassifyT3Fields[model, {"S1", "S2"}, {"S1", "S2"}] === $Failed,
  T3Fail[14, "Cannot register both scalars as heavy and remove F"]
];
fields = GetFields[];
heavyFlags = TrueQ /@ (Lookup[Lookup[fields, #, <||>], Heavy, False] & /@
  {NewScalar1, NewScalar2});
If[heavyFlags =!= {True, True} || KeyExistsQ[fields, NewFermion],
  T3Fail[15, "Invalid scalar-only heavy field registry"]
];
(* Matchete tree was produced with both scalars LIGHT. Promote exactly their
   two physical mass labels to EFTOrder 0 for the scalar threshold. *)
promote[x_] := x /. {
  HoldPattern[Matchete`Coupling[MS1, {}, 1]] :>
    Matchete`Coupling[MS1, {}, 0],
  HoldPattern[Matchete`Coupling[MS2, {}, 1]] :>
    Matchete`Coupling[MS2, {}, 0]
};
demote[x_] := x /. {
  HoldPattern[Matchete`Coupling[MS1, {}, 0]] :>
    Matchete`Coupling[MS1, {}, 1],
  HoldPattern[Matchete`Coupling[MS2, {}, 0]] :>
    Matchete`Coupling[MS2, {}, 1]
};
sourceCounts = (Count[tree, Matchete`Coupling[#, {}, 1], Infinity] & /@
   {MS1, MS2});
(* D(+2) source error, demonstrated by the failed local CheckLagrangian:
     y1[p,r] y2[q,r] / MF[r]  -> NFlavor[r] occurs THREE times.
   This is a diagonal heavy-mass *coefficient*, not a valid Matchete
   contraction. Replacing only MF[r] by symbolic MF[] produces a one-heavy-
   eigenvalue scalar-matching TEMPLATE while retaining both Yukawa matrices
   and the full lepton flavour structure. The original tree is unchanged.
   It is incorrect to evaluate this template once with an average MF and
   claim a full three-mass result. A downstream assembler must evaluate at
   each MF_r using y1[:,r], y2[:,r] and sum the resulting C5 matrices. *)
indexedMFPattern = HoldPattern[Matchete`Coupling[MF, {_}, 0]];
indexedMF = Cases[tree, indexedMFPattern, Infinity];
indexedMFCount = Length[indexedMF];
expectedWilsonTerms = If[MemberQ[{"B"}, modelClass], 2, 4];
indexedMFNFlavor = AllTrue[indexedMF,
   StringContainsQ[T3Text[#], "NFlavor"] &];
If[indexedMFCount =!= seeds[[1]]["TreeWilsonTermCount"] ||
   !indexedMFNFlavor || indexedMFCount =!= expectedWilsonTerms,
  T3Fail[17, "Unexpected diagonal NFlavor mass pattern or operator count: cannot construct a one-eigenvalue template"]
];
Print["[ordinary-hard] Source: ", indexedMFCount,
  " indexed MF[r] factors; constructing one-heavy-eigenvalue matching template"];
perEigenvalueTree = tree /. indexedMFPattern :> Matchete`Coupling[MF, {}, 0];
If[Count[perEigenvalueTree, indexedMFPattern, Infinity] =!= 0 ||
   Count[perEigenvalueTree, Matchete`Coupling[MF, {}, 0], Infinity] =!= indexedMFCount,
  T3Fail[17, "F-mass template substitution did not preserve both LLSS sources"]
];
(* Capture-safe renaming is needed because a textual InputForm import can
   reuse identical dummy names for Flavor and NFlavor indices. *)
renamedTree = CheckAbort[RelabelIndices[perEigenvalueTree], $Aborted];
If[renamedTree === $Failed || renamedTree === $Aborted,
  T3Fail[17, "Per-eigenvalue flavour-index relabelling failed"]
];
(* Reversible mass promotion is audited BEFORE dummy relabelling. *)
promotedUnrelabeled = promote[perEigenvalueTree];
matchInput = promote[renamedTree];
promotedCounts = (Count[matchInput, Matchete`Coupling[#, {}, 0], Infinity] & /@
   {MS1, MS2});
staleCounts = (Count[matchInput, Matchete`Coupling[#, {}, 1], Infinity] & /@
   {MS1, MS2});
reversible = TrueQ[demote[promotedUnrelabeled] === perEigenvalueTree];
promotionOK = AllTrue[sourceCounts, # > 0 &] &&
  sourceCounts === promotedCounts && staleCounts === {0, 0} && reversible;
terms = T3Parts[matchInput];
heavyDegrees = T3Degree /@ terms;
quadratics = Select[terms, T3Degree[#] === 2 &&
    AllTrue[T3AllFields[#], MemberQ[{"NewScalar1", "NewScalar2"}, #] &] &];
free = CheckAbort[UsingFrontEnd[
  RelabelIndices[Expand[FreeLag[NewScalar1, NewScalar2]]]], $Aborted];
If[free === $Aborted || free === $Failed,
  T3Fail[16, "Heavy scalar FreeLag could not be generated"]
];
freeQuadratics = Select[T3Parts[free], T3Degree[#] === 2 &&
    AllTrue[T3AllFields[#], MemberQ[{"NewScalar1", "NewScalar2"}, #] &] &];
quadraticOK = Sort[T3AlphaKey /@ quadratics] ===
   Sort[T3AlphaKey /@ freeQuadratics];
noOddHeavy = Count[heavyDegrees, 1 | 3] === 0;
inputOK = promotionOK && quadraticOK && noOddHeavy;
audit = <|
  "Status" -> If[inputOK, "Passed", "Failed"],
  "Model" -> label,
  "InputSource" -> sourceFiles[[1]],
  "TreeSourceBytes" -> StringLength[treeString],
  "FullGeneralOrdinaryInventoryVerified" -> True,
  "AdditionalQuarticsRetained" -> Lookup[built, "AdditionalQuarticsIncluded", -1],
  "HeavyScalarNames" -> {"NewScalar1", "NewScalar2"},
  "IntegratedFermionStillRegistered" -> KeyExistsQ[fields, NewFermion],
  "RestoredNativeCGCount" -> restored,
  "TreeWilsonTermCount" -> seeds[[1]]["TreeWilsonTermCount"],
  "ModelSpecificExpectedTreeWilsonTermCount" -> expectedWilsonTerms,
  "OriginalLightMassCounts" -> sourceCounts,
  "IndexedFermionMassFactorsInOriginalTree" -> indexedMFCount,
  "DiagonalNFlavorMassIndexEliminatedForTemplate" -> True,
  "ExactFermionMassSpectrumEncodedInTemplate" -> False,
  "PerEigenvalueHeavyMassSumRequired" -> True,
  "OriginalFermionTreeUnmodified" -> True,
  "PromotedHeavyMassCounts" -> promotedCounts,
  "StaleLightMassCounts" -> staleCounts,
  "MassPromotionReversible" -> reversible,
  "PureHeavyQuadraticEqualToFreeLag" -> quadraticOK,
  "PureHeavyQuadraticCount" -> Length[quadratics],
  "OddHeavyScalarDegreeAbsent" -> noOddHeavy,
  "RawOneLoopMatchAttempted" -> False,
  "FullSequentialC5Certified" -> False
|>;
T3Save[auditPath, audit];
If[!inputOK, T3Fail[17, "Saved EFT1 mass/FreeLag audit failed; Match NOT attempted"]];
Print["[ordinary-hard] Checking complete scalar EFT1..."];
check = TimeConstrained[CheckAbort[
  UsingFrontEnd[CheckLagrangian[matchInput]], $Aborted], limit, $TimedOut];
audit = Join[audit, <|"CheckLagrangianResult" -> T3Text[check],
  "CheckLagrangianPassed" -> TrueQ[check]|>];
T3Save[auditPath, audit];
If[check =!= True, T3Fail[18, "Complete EFT1 CheckLagrangian failed; Match NOT attempted"]];
Print["[ordinary-hard] Running independent scalar-only Match at one loop..."];
raw = TimeConstrained[CheckAbort[UsingFrontEnd[
  Match[matchInput, EFTOrder -> 5, LoopOrder -> 1]], $Aborted],
  limit, $TimedOut];
If[MemberQ[{$Failed, $Aborted, $TimedOut}, raw] || AtomQ[raw] ||
   FreeQ[raw, hbar] || T3HasHeavy[raw],
  T3Fail[19, "Scalar-only Match failed, timed out, lacks hbar or retained heavy fields"]
];
Export[FileNameJoin[{outputDir, "hard_raw.wxf"}], raw, "WXF"];
(* Reuse the exact D(-2) reduction order, not a one-scale C5 shortcut. *)
current = raw;
steps = {"GreensSimplify", "EOMSimplify", "EvaluateLoopFunctions",
  "ReplaceEffectiveCouplings"};
Do[
  next = TimeConstrained[CheckAbort[UsingFrontEnd[Switch[steps[[i]],
     "GreensSimplify", GreensSimplify[current],
     "EOMSimplify", EOMSimplify[current],
     "EvaluateLoopFunctions", EvaluateLoopFunctions[current],
     "ReplaceEffectiveCouplings", ReplaceEffectiveCouplings[current],
     _, $Failed]], $Aborted], limit, $TimedOut];
  If[MemberQ[{$Failed, $Aborted, $TimedOut}, next] || AtomQ[next],
    T3Fail[20, "Post-processing failed at " <> steps[[i]]]];
  current = next;
  Export[FileNameJoin[{outputDir, "hard_step_" <> ToString[i] <> ".wxf"}],
    current, "WXF"],
  {i, Length[steps]}
];
hard = Expand[hbar Coefficient[Expand[current], hbar, 1]];
If[hard === 0 || T3HasHeavy[hard],
  T3Fail[21, "Processed hard O(hbar) sector empty or retains heavy fields"]
];
Export[FileNameJoin[{outputDir, "hard_hbar_sector.wxf"}], hard, "WXF"];
extract = TimeConstrained[CheckAbort[
  ExtractWeinbergCoefficient[hard], $Aborted], limit, $TimedOut];
If[!AssociationQ[extract] || Lookup[extract, "Status", ""] =!= "Success" ||
   MissingQ[Lookup[extract, "Coefficient", Missing["NotFound"]]],
  T3Fail[22, "Independent hard Weinberg extraction unavailable"]
];
Export[FileNameJoin[{outputDir, "hard_c5_candidate.txt"}],
  T3Text[extract["Coefficient"]], "Text"];
T3Save[reportPath, <|
  "Status" -> "HardScalarCandidateExported", "Model" -> label,
  "HardOneLoopMatched" -> True,
  "HardHbarSectorExtracted" -> True,
  "PerHeavyMassEigenvalueTemplate" -> True,
  "ModelSpecificExpectedTreeWilsonTermCount" -> expectedWilsonTerms,
  "AdditionalQuarticsRetained" -> Lookup[built, "AdditionalQuarticsIncluded", -1],
  "ThreeHeavyMassNumericalSumImplemented" -> False,
  "OriginalFermionFirstTreeModified" -> False,
  "ScalarMassPromotionOnly" -> True,
  "FermionMatchingRecomputed" -> False,
  "DirectEFT1RunningAdded" -> False,
  "InheritedOneLoopAdded" -> False,
  "WeinbergExtractionStatus" -> extract["Status"],
  "HardWeinbergTermCount" -> Lookup[extract, "TermCount", 0],
  "CompleteC5Assembled" -> False,
  "FullSequentialC5Certified" -> False,
  "SU2NormalizationIndependentlyCertified" -> False
|>];
Print["[ordinary-hard] Independent one-eigenvalue scalar hard template exported; heavy-flavour numerical lift and full C5 still pending."];
Exit[0];
