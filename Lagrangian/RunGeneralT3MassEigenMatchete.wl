(* Experimental native Matchete heavy scalar eigenbasis check (B_-1/C_-1).

   wolframscript -file Lagrangian/RunGeneralT3MassEigenMatchete.wl \
     output/general_t3_mass_eigen_matchete/B_m1_final_parameters.json \
     output/general_t3_mass_eigen_matchete/B_m1_final_matchete.json

   The numerical mass mixing is diagonalised upstream in Python.  Matchete
   is rebuilt with TWO HEAVY Y=-1/2 doublets X,Y and DIAGONAL FreeLag masses.
   All original Yukawa and scalar-potential interaction terms are rewritten
   through the epsilon pseudoreal intertwiner. There is no explicit original
   scalar mixing term in the heavy-mass Lagrangian. The status is not a
   physical CG, Matchete/RGBeta coupling-map or Wilson matching certificate.
*)
ClearAll["Global`*"];
args = Rest[$ScriptCommandLine];
If[Length[args] =!= 2,
  Print["Usage: RunGeneralT3MassEigenMatchete.wl INPUT.json OUTPUT.json"];
  Exit[2]
];
inPath = ExpandFileName[args[[1]]];
outPath = ExpandFileName[args[[2]]];
scriptDir = DirectoryName @ ExpandFileName[$InputFileName];

(* Durable last-completed stage. This file is intentionally distinct from
   the result JSON and survives an abort inside UsingFrontEnd or Matchete. *)
progressPath = FileNameJoin[{DirectoryName[outPath], FileBaseName[outPath] <> ".progress.json"}];
$T3MELastStage = "script-initialized";
$T3MEStageHistory = {};
T3MECheckpoint[stage_String] := Module[{report},
  $T3MELastStage = stage;
  AppendTo[$T3MEStageHistory, stage];
  Print["T3ME_STAGE: ", stage];
  report = <|
    "format" -> "t3-matchete-eigenfield-progress-v1",
    "model" -> If[ValueQ[modelName], modelName, "Unknown"],
    "last_stage" -> stage,
    "stage_history" -> $T3MEStageHistory,
    "physical_matching_performed" -> False
  |>;
  Quiet @ Check[Export[progressPath, report, "RawJSON"], $Failed];
  stage
];
T3MECheckpoint["script-initialized"];

T3MEOutput[data_Association] := Module[{},
  If[!DirectoryQ[DirectoryName[outPath]],
    CreateDirectory[DirectoryName[outPath], CreateIntermediateDirectories -> True]];
  Export[outPath, data, "RawJSON"];
  Print["Mass-eigenstate Matchete result: ", Lookup[data,"status","Unknown"],
        " -> ",outPath];
  data
];
T3MEFail[stage_String, note_: ""] := (
  T3MEOutput[<|"status" -> stage, "diagnostic" -> note,
   "model" -> Lookup[input,"model","Unknown"],
   "mass_basis_ordinary_check_passed" -> False,
   "last_stage" -> $T3MELastStage,
   "stage_history" -> $T3MEStageHistory,
   "all_old_scalar_fields_removed" -> False,
   "physical_matching_performed" -> False,
   "component_basis_map_certified" -> False|>];
  Throw[$Failed,"T3MassEigenStop"]
);

input = Quiet @ Check[Import[inPath,"RawJSON"],$Failed];
If[!AssociationQ[input] ||
   Lookup[input,"format",""] =!= "t3-experimental-mass-eigen-matchete-parameters-v1" ||
   Lookup[input,"status",""] =!= "NumericMassEigenstateBoundaryPreparedMatcheteNotRun" ||
   !MemberQ[{"B_-1","C_-1"},Lookup[input,"model",""]] ||
   Lookup[input,"quartic_source_count",0] =!= 17 ||
   Lookup[input,"yukawa_source_count",0] =!= 2,
  Print["Invalid Matchete mass-basis parameter schema"]; Exit[2]
];
modelName = input["model"];
class = StringTake[modelName,1];
expectedDims = If[class === "B",{2,2,1,-1},{2,2,3,-1}];
If[Lookup[input,"dimensions_alpha",{}] =!= expectedDims,
  Print["Model representation mismatch"];Exit[2]
];

(* Access scalar flavour-space U and squared masses via explicit numeric JSON.
   Never ToExpression a string from an external file.
   IMPORTANT: do not assign to global symbol u: Matchete uses u as the
   SM up-quark field label during LoadModel["SM"]. Giving it matrix OwnValues
   corrupts DefineField::FieldLabel. The $T3ME prefix avoids field labels. *)
complex2[obj_] := Table[
  N[obj["real"][[i,j]]] + I N[obj["imag"][[i,j]]],
  {i,2},{j,2}
];
$T3MEUnitaryMatrix = Quiet @ Check[complex2[input["unitary_U"]],$Failed];
aMatrix = Quiet @ Check[complex2[input["original_mass_squared_A_gev2"]],$Failed];
eigenSq = Lookup[input,"heavy_mass_squared_gev2",{}];
If[!MatrixQ[$T3MEUnitaryMatrix,NumericQ] || Dimensions[$T3MEUnitaryMatrix] =!= {2,2} ||
   !MatrixQ[aMatrix,NumericQ] || Dimensions[aMatrix] =!= {2,2} ||
   Length[eigenSq] =!= 2 || !And@@(NumericQ /@ eigenSq) ||
   !And@@(Positive /@ eigenSq),
  Print["Malformed or tachyonic mass eigenstate input"];Exit[2]
];
matScale = Max[1, Max[Abs[Flatten[aMatrix]]]];
unitarityError = Max[Abs[Flatten[ConjugateTranspose[$T3MEUnitaryMatrix].$T3MEUnitaryMatrix - IdentityMatrix[2]]]];
massError = Max[Abs[Flatten[ConjugateTranspose[$T3MEUnitaryMatrix].aMatrix.$T3MEUnitaryMatrix - DiagonalMatrix[eigenSq]]]];
If[unitarityError > 2 10^-10 || massError > 2 10^-9 matScale,
  Print["Input U/mass diagonalization failed",{unitarityError,massError}];Exit[2]
];

T3MECheckpoint["matchete-load-start"];
loadStatus=CheckAbort[Check[UsingFrontEnd[Needs["Matchete`"];True],$Failed],$Aborted];
T3MECheckpoint["matchete-load-returned"];
If[loadStatus =!= True,Print["Matchete package unavailable"];Exit[3]];
Get[FileNameJoin[{scriptDir,"model","T3ModelCatalog.wl"}]];
Get[FileNameJoin[{scriptDir,"model","LagrangianBuilder.wl"}]];
model=T3ModelFromClass[class,-1];
T3MECheckpoint["source-files-loaded"];
If[model === $Failed, Print["Cannot construct original model"];Exit[3]];

(* Construct native epsilon-intertwiner CGs from the SAME Matchete SU(2)
   representation definition used for the original T3 model. The first slot
   is a virtual original-field index, and the second is a mass-eigenfield
   index. The CC and UU tensor orientation must be checked in a future
   component convention map; do not infer the numerical component phase here. *)
T3MEBridge[conjugated_, name_Symbol] := Module[
  {info,t,canonicalEpsilon,phase},
  info = PhysicalSU2InvariantBasis[{2,2},{conjugated,conjugated}];
  If[!AssociationQ[info] || Lookup[info,"PhysicalCount",0] =!= 1,
    Return[$Failed]];
  t = Normal[info["Basis"][[1]]];
  If[Dimensions[t] =!= {2,2},Return[$Failed]];
  canonicalEpsilon = {{0,1},{-1,0}};
  (* Two fundamental or two antifundamental SU(2) indices have a unique
     antisymmetric invariant. Align its arbitrary native complex phase
     with our numerical epsilon component convention, exactly. *)
  phase = FullSimplify[Total[Flatten[Conjugate[canonicalEpsilon] t]]/2];
  If[TrueQ[PossibleZeroQ[phase]] || !FreeQ[phase,Indeterminate|ComplexInfinity],
    Return[$Failed]];
  If[!TrueQ[FullSimplify[t == phase canonicalEpsilon]],
    Print["Native doublet tensor differs from epsilon up to phase: ",
      ToString[t,InputForm]];
    Return[$Failed]];
  t = FullSimplify[t/phase];
  If[!TrueQ[FullSimplify[t == canonicalEpsilon]] ||
     !TrueQ[FullSimplify[t.ConjugateTranspose[t] == IdentityMatrix[2]]],
    Return[$Failed]];
  If[Check[DefineCG[name,info["CGReps"],t];True,$Failed] =!= True,
    Return[$Failed]];
  <|"CG" -> SymbolName[name],"CGReps"->ToString[info["CGReps"],InputForm],
    "native_intertwiner_inputform"->ToString[t,InputForm],
    "native_phase_before_alignment_inputform"->ToString[phase,InputForm],
    "phase_relative_to_matchete_component_epsilon_verified"->True,
    "physical_RGBeta_component_map_certified"->False|>
];

(* Matchete's public shortcuts are evaluated into internal objects, e.g.
     lambda[]  -> Coupling[lambda, {}, 0]
     S[i]      -> Field[S, Scalar, {i}, ...]
   Hence searching ONLY applied heads such as lambda[___] or S[___] is wrong.
   Inspect the exact label atom at any level of the expanded expression.
   Verbatim prevents a dynamically supplied name being treated as a pattern.
   False-positive control is included; the raw operator must actually contain
   this label, not merely appear in the expected registry. *)
T3MECouplingPresentQ[expr_, symbol_Symbol] := With[{label = symbol},
  !FreeQ[expr, Verbatim[label], {0, Infinity}]
];

T3MECouplingPresenceSelfTest[] := Module[
  {present = Unique["T3MEPresent"], absent = Unique["T3MEAbsent"], test, internal},
  (* Check both public shortcuts AND the Matchete internal Coupling wrapper.
     Using a different head checks that no head-specific assumption survives. *)
  test = 2 present[] + 3 Coupling[present, {}, 0] + 4 present[]^2;
  internal = Coupling[present, {}, 0];
  TrueQ[T3MECouplingPresentQ[test, present]] &&
  TrueQ[T3MECouplingPresentQ[internal, present]] &&
  !TrueQ[T3MECouplingPresentQ[test, absent]] &&
  !TrueQ[T3MECouplingPresentQ[Coupling[absent, {}, 0], present]]
];
If[!TrueQ[T3MECouplingPresenceSelfTest[]],
  Print["T3ME: Matchete Coupling[label,...] detection self-test failed"];
  Exit[4]
];

(* Field labels are likewise embedded inside Matchete Field[...], not
   necessarily used as expression heads. Missing old labels must be checked
   on the actual internalized objects after replacement. *)
T3MEOldFieldsAbsentQ[expr_] :=
  FreeQ[expr, Verbatim[NewScalar1] | Verbatim[NewScalar2], {0, Infinity}];

(* Report actual internal object signatures if the Matchete version differs
   from the expected Field[label,Scalar,{index},...] syntax. This diagnostic
   is bounded and never marks an unsupported representation successful. *)
T3MEScalarSignatureReport[expr_] := <|
  "raw_source_field_examples" -> Take[
    DeleteDuplicates[ToString[#, InputForm] & /@
      Cases[expr,
        f : HoldPattern[Field[NewScalar1 | NewScalar2, ___]] :> f,
        {0, Infinity}]], UpTo[4]],
  "raw_source_coupling_example" -> Take[
    DeleteDuplicates[ToString[#, InputForm] & /@
      Cases[expr,
        c : HoldPattern[Coupling[_, ___]] :> c,
        {0, Infinity}]], UpTo[2]]
|>;

(* A source Field[... ] already has a Matchete-typed SU(2) index
   Index[label,T3BSMd2]. The newly introduced bridge contraction must use
   the SAME typed index form. A bare Module-local label (j$...) makes
   Matchete interpret it as an undefined representation when validating CG.

   The new scalar shortcut NewScalarX[j] / NewScalarY[j] is deliberately
   given the bare label j; DefineField converts it to Index[j,T3BSMd2],
   matching the explicitly typed CG slot. *)
T3MEBridgeIndicesTypedQ[expr_] := Module[{bad},
  bad = Cases[expr,
    HoldPattern[CG[T3MEBridgeBB | T3MEBridgeUU, {a_, b_}]] /;
      (!MatchQ[a, Index[_, T3BSMd2]] ||
       !MatchQ[b, Index[_, T3BSMd2]]),
    {0, Infinity}];
  If[bad =!= {},
    Print["Untyped SU(2) bridge CG indices: ",
      ToString[Take[bad, UpTo[4]], InputForm]];
    Return[False]
  ];
  True
];

(* For source fields the index is ALREADY a Matchete Index[label,rep].
   Passing it back through the public shortcut NewScalarX[i] can wrap the
   index twice. Build the same internal Field object directly and verify it.
   This only changes the transport of the existing typed index; it does not
   alter the SU(2) bridge or the interaction coefficients. *)
T3MENewTypedScalar[label_Symbol, index_Index] :=
  Field[label, Scalar, {index}, {}];

T3MENewScalarIndexAudit[expr_] := Module[{fields, bad},
  fields = Cases[expr,
    HoldPattern[Field[(NewScalarX | NewScalarY), Scalar, inds_, rest___]] :> inds,
    {0, Infinity}];
  bad = Select[fields,
    !MatchQ[#, {Index[name_, T3BSMd2]} /; Head[name] =!= Index] &];
  If[fields === {} || bad =!= {},
    Print["Malformed rotated eigen-scalar index list(s): ",
      ToString[Take[bad, UpTo[6]], InputForm],
      "; eigenfield instances found=", Length[fields]];
    Return[False]
  ];
  True
];

(* Keep the numerically checked U in the input, but present its four complex
   entries as NEUTRAL SYMBOLIC Matchete couplings in interactions.  An invariant
   tensor should not acquire machine-precision coefficients during exact CG
   decomposition.  The mapping to the external numerical U is exported below.
   No claim of exact unitary constraints is made for independent spurions. *)
T3MEDefineMixingSpurions[] := Module[{registered},
  registered = {
    Check[DefineCoupling[T3MEU11, SelfConjugate -> False]; True, $Failed],
    Check[DefineCoupling[T3MEU12, SelfConjugate -> False]; True, $Failed],
    Check[DefineCoupling[T3MEU21, SelfConjugate -> False]; True, $Failed],
    Check[DefineCoupling[T3MEU22, SelfConjugate -> False]; True, $Failed]
  };
  If[registered =!= {True, True, True, True}, Return[$Failed]];
  $T3MEInteractionMixingMatrix = {
    {T3MEU11[], T3MEU12[]},
    {T3MEU21[], T3MEU22[]}
  };
  If[!MatrixQ[$T3MEInteractionMixingMatrix] ||
     Dimensions[$T3MEInteractionMixingMatrix] =!= {2,2}, Return[$Failed]];
  True
];
T3MEMixingSpurionNumericalMap[] := Module[{real, imag, labels},
  real = input["unitary_U"]["real"];
  imag = input["unitary_U"]["imag"];
  labels = {{"T3MEU11", "T3MEU12"}, {"T3MEU21", "T3MEU22"}};
  Association @ Flatten @ Table[
    labels[[i,j]] -> <|"real" -> real[[i,j]], "imag" -> imag[[i,j]]|>,
    {i,2}, {j,2}
  ]
];

(* Both old doublets must be eliminated. Replacements use unique dummy
   indices on each insertion, then RelabelIndices resolves clashes. *)
T3MERotateInteraction[expr_] := Module[{transformed, before},
  before = expr;
  (* IMPORTANT: the source is already evaluated to Matchete Field objects.
     Match Field[...] and Bar[Field[...]], NOT the NewScalar1[i] shortcuts.
     Barred replacements precede bare replacements so conjugate fields are
     not accidentally rewritten using the unconjugated coefficients. *)
  transformed = before /. {
    HoldPattern[Bar[Field[NewScalar1, Scalar, {i_}, rest___]]] :>
      (Bar[$T3MEInteractionMixingMatrix[[1,1]]] Bar[T3MENewTypedScalar[NewScalarX,i]] +
       Bar[$T3MEInteractionMixingMatrix[[1,2]]] Bar[T3MENewTypedScalar[NewScalarY,i]]),
    HoldPattern[Bar[Field[NewScalar2, Scalar, {i_}, rest___]]] :>
      Module[{j},
        -$T3MEInteractionMixingMatrix[[2,1]] CG[T3MEBridgeUU,{i,Index[j,T3BSMd2]}] NewScalarX[j] -
         $T3MEInteractionMixingMatrix[[2,2]] CG[T3MEBridgeUU,{i,Index[j,T3BSMd2]}] NewScalarY[j]
      ],
    HoldPattern[Field[NewScalar1, Scalar, {i_}, rest___]] :>
      ($T3MEInteractionMixingMatrix[[1,1]] T3MENewTypedScalar[NewScalarX,i] +
       $T3MEInteractionMixingMatrix[[1,2]] T3MENewTypedScalar[NewScalarY,i]),
    HoldPattern[Field[NewScalar2, Scalar, {i_}, rest___]] :>
      Module[{j},
        -Bar[$T3MEInteractionMixingMatrix[[2,1]]] CG[T3MEBridgeBB,{i,Index[j,T3BSMd2]}] Bar[NewScalarX[j]] -
         Bar[$T3MEInteractionMixingMatrix[[2,2]]] CG[T3MEBridgeBB,{i,Index[j,T3BSMd2]}] Bar[NewScalarY[j]]
      ]
  };
  If[!T3MEOldFieldsAbsentQ[transformed],
    Print["Old scalar label(s) still occur after internal Field replacement."];
    Print["Raw original Matchete signatures: ",
      ToString[T3MEScalarSignatureReport[before], InputForm]];
    Return[$Failed]
  ];
  If[!TrueQ[T3MEBridgeIndicesTypedQ[transformed]],
    Print["Matchete bridge-index preflight failed before CheckLagrangian"];
    Return[$Failed]
  ];
  If[!TrueQ[T3MENewScalarIndexAudit[transformed]],
    Print["Rotated Field eigen-scalar index audit failed"];
    Return[$Failed]
  ];
  Expand[transformed] // RelabelIndices
];
(* Diagnose gauge-tensor failures without recomputing or redefining any of
   the source Clebsch-Gordan tensors.  Expand the saved valid baseline LInt,
   partition its exact terms by coupling and CG labels, then check each
   disjoint piece with the diagonal FreeLag.  Every original term must be
   accounted for. This is a diagnostic, NOT a weaker acceptance path. *)
T3MEBaselineTermKey[term_] := Module[{couplings, cgs},
  couplings = DeleteDuplicates @ Cases[term,
    HoldPattern[Coupling[label_Symbol, ___]] :> SymbolName[label],
    {0, Infinity}];
  cgs = DeleteDuplicates @ Cases[term,
    HoldPattern[CG[label_Symbol, ___]] :> SymbolName[label],
    {0, Infinity}];
  If[couplings === {}, Return[$Failed]];
  StringRiffle[Sort[couplings], "+"] <> "|CG:" <> StringRiffle[Sort[cgs], "+"]
];
T3MEPartitionBaselineInteractions[expr_] := Module[
  {expanded, terms, keys, pairs, groups, parts},
  expanded = Expand[expr];
  terms = If[Head[expanded] === Plus, List @@ expanded, {expanded}];
  keys = T3MEBaselineTermKey /@ terms;
  If[MemberQ[keys, $Failed] || terms === {} ||
     Length[terms] =!= Length[keys], Return[$Failed]];
  pairs = Transpose[{keys, terms}];
  groups = GatherBy[pairs, First];
  parts = MapIndexed[
    <|"Name" -> ("BaselinePart" <> ToString[First[#2]]),
      "Class" -> "BaselinePartition", "Signature" -> #1[[1,1]],
      "Expression" -> Total[#1[[All,2]]]|> &,
    groups];
  If[parts === {} ||
     !TrueQ[Expand[Total[Lookup[parts, "Expression", {}]]-expanded] === 0],
    Return[$Failed]];
  parts
];

T3MEIsolateNonInvariantSource[built_Association, free_, extra_List] := Module[
  {baseline, candidates, failures = {}, entry, converted, check,
   freeCheck, tested = 0},
  T3MECheckpoint["failure-isolation-started"];
  freeCheck = CheckAbort[
    Quiet @ Check[CheckLagrangian[(built["LSM"]+free) // RelabelIndices], $Failed],
    $Aborted];
  Print["T3ME isolated diagonal FreeLag check: ", InputForm[freeCheck]];
  If[freeCheck =!= True,
    Return[<|"status"->"DiagonalFreeOnlyCheckFailed",
      "free_only_result"->ToString[freeCheck,InputForm],
      "tested_interactions"->0|>]
  ];
  baseline = T3MEPartitionBaselineInteractions[built["LInt"]];
  If[!ListQ[baseline],
    Return[<|"status"->"BaselinePartitionFailed",
      "free_only_result"->"True"|>]
  ];
  Print["T3ME disjoint baseline source partitions: ", Length[baseline]];
  candidates = Join[baseline, extra];
  If[Length[extra] =!= 7 || Length[baseline] < 8,
    Return[<|"status"->"CandidateInventoryWrong",
      "baseline_parts"->Length[baseline], "extra"->Length[extra]|>]
  ];
  Do[
    entry = candidates[[k]];
    tested = k;
    T3MECheckpoint["failure-isolation-operator-" <> ToString[k]];
    converted = T3MERotateInteraction[entry["Expression"]];
    check = If[converted === $Failed, $Failed,
      CheckAbort[Quiet @ Check[
        CheckLagrangian[(built["LSM"] + free + converted) // RelabelIndices],
        $Failed], $Aborted]
    ];
    Print["T3ME individual rotated check ", entry["Name"],
      " ", Lookup[entry,"Signature","AdditionalQuartic"],
      " -> ",InputForm[check]];
    If[check =!= True,
      AppendTo[failures, <|
        "name" -> entry["Name"],
        "class" -> Lookup[entry,"Class","Unknown"],
        "source_signature" -> Lookup[entry,"Signature","ExtraQuartic"],
        "check_result" -> ToString[check,InputForm],
        "rotated_has_old_fields" -> If[converted === $Failed,
          Null,!T3MEOldFieldsAbsentQ[converted]],
        "rotated_new_field_indices_well_typed" -> If[converted === $Failed,
          Null,T3MENewScalarIndexAudit[converted]]
      |>];
      Print["T3ME first failing rotated source operator: ",
        ToString[Last[failures],InputForm]];
      Break[]
    ],
    {k,Length[candidates]}
  ];
  <|
    "status" -> If[failures === {}, "IndividualChecksPass", "IndividualCheckFailed"],
    "free_only_result"->"True",
    "baseline_source_partitions"->Length[baseline],
    "tested_interactions" -> tested,
    "first_failed" -> failures
  |>
];

T3MEPerform[] := Module[
 {built, sourceInt, extraQ, originals, oldPresent, status,
  b1,b2,removeA,removeB,xstatus,ystatus,free,rotated,lag,check,
  detailed,heavy,fullCount,baselineNames,oldYukawas,results,
  extraSymbols,extrasPresent,missingExtra,spurionStatus,missingSpurions,
  failureDetail},
 If[!MatchQ[$FrontEnd,_FrontEndObject],
   T3MEFail["FrontEndUnavailable"]];
 If[Environment["T3_EXPERIMENTAL_GENERAL_QUARTICS"] === "1",
   T3MEFail["DuplicateExperimentalQuarticsOptIn"]];
 T3MECheckpoint["baseline-build-started"];
 built=BuildT3Lagrangian[model];
 If[!AssociationQ[built] || built["Status"] =!= "Success",
   T3MEFail["OriginalT3BaselineInvalid"]];
 baselineNames=built["AllowedInteractions"];
 If[Length[baselineNames] =!= 11 ||
    Count[baselineNames,_?(StringStartsQ[#,"Yukawa"]&)] =!= 2,
   T3MEFail["OriginalT3BaselineVertexInventoryChanged",
     ToString[baselineNames,InputForm]]];
 T3MECheckpoint["extra-quartic-generation-started"];
 extraQ = T3ExtendedBuildExtraQuarticCandidates[model];
 If[!ListQ[extraQ] || Length[extraQ] =!= 7 ||
    Length[DeleteDuplicates[Lookup[extraQ,"Name",{}]]] =!= 7,
   T3MEFail["SevenSourceQuarticsNotGenerated"]];
 originals=Join[
   (<|"name"->#,"class"->"Baseline"|>& /@ baselineNames),
   (<|"name"->#["Name"],"class"->"AdditionalQuartic"|>& /@ extraQ)
 ];
 T3MECheckpoint["original-interaction-assembly-started"];
 sourceInt = Expand[built["LInt"] + Total[Lookup[extraQ,"Expression",{}]]] // RelabelIndices;
 If[T3MEOldFieldsAbsentQ[sourceInt],
   T3MEFail["OriginalScalarLabelsAbsentInInternalSource",
     ToString[T3MEScalarSignatureReport[sourceInt], InputForm]]];
 extraSymbols = Table[
   Symbol["lambdaT3ExtraQ" <> StringJoin[ToString /@ extraQ[[k,"Counts"]]] <>
          "I" <> Last[StringSplit[extraQ[[k,"Name"]],"_"]]],
   {k,Length[extraQ]}];
 missingExtra = Select[extraSymbols,
   !TrueQ[T3MECouplingPresentQ[sourceInt, #]] &];
 If[Length[DeleteDuplicates[extraSymbols]] =!= 7 || missingExtra =!= {},
   Print["Missing additional Matchete Coupling labels: ",
     ToString[missingExtra, InputForm]];
   Print["Observed source internal signatures: ",
     ToString[T3MEScalarSignatureReport[sourceInt], InputForm]];
   T3MEFail["AdditionalQuarticCouplingsAbsentInOriginalSource",
     ToString[missingExtra, InputForm]]];
 oldYukawas = T3MECouplingPresentQ[sourceInt,y1] &&
             T3MECouplingPresentQ[sourceInt,y2];
 If[!oldYukawas,T3MEFail["OriginalYukawaVerticesMissing"]];

 (* Define the CG intertwiners BEFORE deregistering original field labels;
    they depend only on the already registered SU(2) representation. *)
 Print["Internalized-source representative objects: ",
   ToString[T3MEScalarSignatureReport[sourceInt], InputForm]];
 T3MECheckpoint["bridge-barred-started"];
 b1=T3MEBridge[True,T3MEBridgeBB];
 T3MECheckpoint["bridge-unbarred-started"];
 b2=T3MEBridge[False,T3MEBridgeUU];
 If[!AssociationQ[b1] || !AssociationQ[b2],
   T3MEFail["NativeSU2PseudorealBridgeFailed"]];

 (* Retain F as a heavy field; replace only S1/S2. The Matchete heavy mass
    basis is now represented by fresh diagonal kinetic/mass terms from FreeLag. *)
 T3MECheckpoint["remove-old-scalar-1-started"];
 removeA=Quiet@Check[RemoveField[NewScalar1];True,$Failed];
 T3MECheckpoint["remove-old-scalar-2-started"];
 removeB=Quiet@Check[RemoveField[NewScalar2];True,$Failed];
 If[removeA =!= True || removeB =!= True,
   T3MEFail["OriginalScalarDeregistrationFailed"]];
 T3MECheckpoint["register-mixing-spurions-started"];
 spurionStatus = T3MEDefineMixingSpurions[];
 If[spurionStatus =!= True,
   T3MEFail["SymbolicMixingSpurionRegistrationFailed"]];
 T3MECheckpoint["register-heavy-X-started"];
 xstatus=Check[
   DefineField[NewScalarX,Scalar,SelfConjugate->False,
     Indices->{T3BSMd2},Charges->{U1Y[-1/2]},Mass->{Heavy,MSX}];True,
   $Failed];
 T3MECheckpoint["register-heavy-Y-started"];
 ystatus=Check[
   DefineField[NewScalarY,Scalar,SelfConjugate->False,
     Indices->{T3BSMd2},Charges->{U1Y[-1/2]},Mass->{Heavy,MSY}];True,
   $Failed];
 If[xstatus =!= True || ystatus =!= True,
   T3MEFail["HeavyMassEigenfieldRegistrationFailed"]];
 T3MECheckpoint["diagonal-FreeLag-started"];
 free=Quiet@Check[FreeLag[NewFermion,NewScalarX,NewScalarY],$Failed];
 If[free === $Failed,T3MEFail["DiagonalFreeLagFailed"]];
 T3MECheckpoint["rotate-original-interactions-started"];
 rotated=T3MERotateInteraction[sourceInt];
 If[rotated === $Failed, T3MEFail["OriginalScalarFieldsRemainAfterRotation"]];
 oldPresent=!T3MEOldFieldsAbsentQ[rotated];
 If[oldPresent,T3MEFail["OriginalScalarFieldsRemain"]];
 If[!T3MECouplingPresentQ[rotated,y1] ||
    !T3MECouplingPresentQ[rotated,y2],
   T3MEFail["TransformedYukawaVertexMissing"]];
 missingSpurions = Select[{T3MEU11,T3MEU12,T3MEU21,T3MEU22},
   !TrueQ[T3MECouplingPresentQ[rotated,#]] &];
 If[missingSpurions =!= {},
   T3MEFail["MixingSpurionsMissingAfterRotation",
     ToString[missingSpurions,InputForm]]];
 extrasPresent = And @@ (T3MECouplingPresentQ[rotated,#] & /@ extraSymbols);
 If[!extrasPresent || !T3MECouplingPresentQ[rotated,lambdaT3],
   T3MEFail["TransformedQuarticVertexMissing"]];
 T3MECheckpoint["assemble-eigenstate-Lagrangian-started"];
 lag=(built["LSM"] + free + rotated) // Expand // RelabelIndices;
 If[!T3MEOldFieldsAbsentQ[lag],
   T3MEFail["OldScalarFieldInFinalLagrangian"]];
 T3MECheckpoint["final-index-preflight-started"];
 If[!TrueQ[T3MENewScalarIndexAudit[lag]] ||
    !TrueQ[T3MEBridgeIndicesTypedQ[lag]],
   T3MEFail["FinalEigenstateTensorIndexPreflightFailed"]];
 Print["T3ME representative final eigenfield internal objects: ",
   ToString[Take[DeleteDuplicates @ Cases[lag,
      f : HoldPattern[Field[NewScalarX | NewScalarY, ___]] :>
        ToString[f, InputForm], {0, Infinity}], UpTo[4]], InputForm]];

 T3MECheckpoint["ordinary-CheckLagrangian-started"];
 check=CheckAbort[CheckLagrangian[lag],$Aborted];
 T3MECheckpoint["detailed-CheckLagrangian-started"];
 detailed=CheckAbort[CheckLagrangian[lag,DetailedOutput->True],$Aborted];
 heavy=AssociationQ[detailed] &&
   TrueQ[Lookup[detailed,"HeavyMassBasis",False]];
 fullCount=Length[extraQ]+Length[baselineNames];
 results=<|
  "format"->"t3-matchete-native-heavy-mass-eigenbasis-v1",
  "model"->modelName,
  "mass_basis_ordinary_check_passed"->TrueQ[check],
  "heavy_mass_basis_detailed_check_passed"->heavy,
  "all_old_scalar_fields_removed"->True,
  "all_17_quartics_present"->True,
  "both_yukawa_interactions_present"->True,
  "source_baseline_interactions"->Length[baselineNames],
  "source_extra_quartics"->Length[extraQ],
  "source_new_interaction_count"->fullCount,
  "matchete_rotation_uses_exact_symbolic_spurions"->True,
  "rotation_spurion_values_from_numerical_U"->T3MEMixingSpurionNumericalMap[],
  "source_interactions"->originals,
  "diagonal_heavy_fields"->{"NewScalarX","NewScalarY"},
  "epsilon_bridge_BB"->b1,
  "epsilon_bridge_UU"->b2,
  "input_unitarity_error"->unitarityError,
  "input_diagonalization_error_gev2"->massError,
  "matching_to_two_scalar_thresholds_performed"->False,
  "physical_matching_performed"->False,
  "component_basis_map_certified"->False,
  "matchete_rgbeta_quartic_normalization_certified"->False
 |>;
 T3MECheckpoint["both-Matchete-checks-returned"];
 If[!TrueQ[check] || !heavy,
   Print["Native mass eigenbasis CheckLagrangian failed: ",InputForm[check]];
   Print["Native detailed diagnostics: ",ToString[detailed,InputForm]];
   failureDetail = T3MEIsolateNonInvariantSource[
     built, free, extraQ];
   Return[Join[results,<|"status"->"NativeHeavyMassBasisCheckFailed",
       "diagnostic"->ToString[detailed,InputForm],
       "isolated_source_failure"->failureDetail|>]]];
 T3MECheckpoint["complete-valid-eigenbasis"];
 Join[results,<|"status"->
   "MatcheteMassEigenHeavyBasisPassedPhysicalMatchingNotPerformed"|>]
];

run = CheckAbort[
  UsingFrontEnd[Catch[T3MEPerform[],"T3MassEigenStop"]],$Aborted];
If[run === $Aborted || run === $Failed || !AssociationQ[run],
  If[!FileExistsQ[outPath],
    T3MEOutput[<|"status"->"UncaughtMatcheteError","model"->modelName,
      "last_stage"->$T3MELastStage,
      "stage_history"->$T3MEStageHistory,
      "kernel_result_inputform"->ToString[run,InputForm],
      "message_list_inputform"->ToString[$MessageList,InputForm],
      "physical_matching_performed"->False,
      "component_basis_map_certified"->False|>]];
  Exit[5]
];
T3MEOutput[run];
If[run["status"] =!= "MatcheteMassEigenHeavyBasisPassedPhysicalMatchingNotPerformed",
  Exit[5]
];
Exit[0];
