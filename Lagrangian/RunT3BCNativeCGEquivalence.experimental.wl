(* Workstream C complete EFT1 residual exporter, derived only from archived
   RunT3BCXFromSavedY.experimental.wl (Git blob 7647eea9d05f05d70377bf11c237c358257bec14).
   Experimental standalone FORK; no changes to production source, no Match executed
   when T3_BC_FULL_RESIDUAL_ONLY=1.  A passing check is not a physical C5 certificate. *)
(* V3 independent fresh-kernel X-only attempt from the exact Y-only EOMEFT artifact.
   The F-first EOMEFT equivalence and physical one-loop matching remain UNPROVEN.
   Original Y source and all scalar quartics are retained in the independent canonical reassembly.
   Y matching is NOT performed in this invocation. Do not use for physical C5. *)
(* Generated from the exact user-uploaded RunT3BCYOnlyMatching.experimental.wl.
   Do not use as a source of physical Wilson matching certificates.
   Independent candidate is accepted only following full Matchete checks,
   and algebraic equivalence to saved EOMEFT is NOT established. *)
(* GENERATED EXPERIMENTAL Y-ONLY MATCHING FROM INSTALLED NATIVE PROBE.
   This source is NOT a certified finite threshold coefficient.
   Source fidelity and exact edits are checked by the Python builder. *)
(* Experimental second scalar-threshold matching from saved canonical tree EOMEFT (B_-1/C_-1).
   Fork of the independently live-validated mass-eigenfield Lagrangian builder.
   In a fresh kernel F is removed and X,Y are registered HEAVY.
   Match at EFTOrder -> 5, LoopOrder -> 0 then 1 on the *saved canonical tree EOMEFT*.
   Subtract Matchete-matched EFTs to isolate scalar-loop LLHH candidates.
   No inherited hard one-loop piece, EFT1 running, counterterm certification,
   complete scalar threshold C5 or final neutrino mass is claimed.

   wolframscript -file Lagrangian/RunGeneralT3EFT1ScalarThresholdProbe.wl \
     output/general_t3_mass_eigen_matchete/B_m1_final_parameters.json \
     output/general_t3_mass_eigen_matchete/eft1_fermion_first_oneloop/B_m1_final_eft1_fermion_first_oneloop_eft1_tree.inputform.txt \
     output/general_t3_mass_eigen_matchete/eft1_scalar_threshold_probe/B_m1_final_eft1_scalar_threshold_probe.json

   The numerical mass mixing is diagonalised upstream in Python. Matchete
   is rebuilt with F HEAVY and TWO HEAVY Y=-1/2 eigen-doublets X,Y.  The
   scalar kinetic/mass terms are diagonal, and all native SU(2) interaction
   tensors are the independently checked mass-basis tensors. Both scalars are integrated out, while F has already been removed
   by the first-threshold tree-level EFT1 calculation. The EOMEFT
   continuation is deliberately NOT the ReplaceEffectiveCouplings-expanded MatchedEFT.
*)
ClearAll["Global`*"];
args = Rest[$ScriptCommandLine];
If[Length[args] =!= 3,
  Print["Usage: RunGeneralT3EFT1ScalarThresholdProbe.wl INPUT.json EFT1_TREE.inputform.txt OUTPUT.json"];
  Exit[2]
];
inPath = ExpandFileName[args[[1]]];
treePath = ExpandFileName[args[[2]]];
outPath = ExpandFileName[args[[3]]];
scriptDir = DirectoryName @ ExpandFileName[$InputFileName];
(* Opt-in leading Weinberg-sector projection.  This does NOT fix Hermiticity
   of the FULL EFT and cannot be represented as complete scalar EFT matching.
   Its 1-loop LLHH subtheory keeps every LLSS and Higgs-bearing quartic.
   Four-heavy-scalar interactions cannot supply the required two external H
   lines together with LLSS at one loop.  Never run this mode silently. *)
$T3MELeadingC5Only = (Environment["T3_EFT1_LEADING_C5_ONLY"] === "1");
(* Opt-in independent tree construction. Never pretend EOMEFT equality. *)
$T3MEIndependentCanonicalReassembly =
  (Environment["T3_BC_Y_INDEPENDENT_CANONICAL"] === "1");
$T3MEIndependentCheckOnly =
  (Environment["T3_BC_Y_INDEPENDENT_CHECK_ONLY"] === "1");
(* Workstream C: independent, no-Match native CG component condition probe.
   Original native CG source preserved; only an exact component check is added. *)
(* Workstream C: independent, no-Match residual export.  This fork is based
   on the verified archived native source.  Never replace the production
   Y-only or X-only scripts with this diagnostic. *)
$T3BCFullResidualOnly = (Environment["T3_BC_FULL_RESIDUAL_ONLY"] === "1");
If[TrueQ[$T3BCFullResidualOnly],
  If[TrueQ[$T3MELeadingC5Only],
    Print["Full EFT1 residual requires T3_EFT1_LEADING_C5_ONLY unset."];Exit[2]];
  $T3MEIndependentCanonicalReassembly = True;
  $T3MEIndependentCheckOnly = True;
];


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
(* The original T3 topology has H H S1 Bar[S2]. After the pseudoreal
   rotation, it becomes H H (X,Y) (X,Y).  Do NOT leave the product of the
   old four-leg CG and a two-leg epsilon CG in the internal expression:
   that product is not a registered invariant for the four NEW legs on the
   installed Matchete version.  Construct the one-dimensional physical HHXX
   tensor directly, using the existing basis generator.

   Exact tensor relation (in Matchete's own SU(2) basis):
     T_old[i,j,k,l] eps[l,m] = alphaCG T_HHXX[i,j,k,m].
   alphaCG is computed, not guessed. Its validity is tested over ALL 16
   tensor components, and HHXX symmetry in the last pair is mandatory.
   The two X/Y orderings are then collected using bosonic commutativity.
   This is a Matchete-group-tensor identity only, not an RGBeta map. *)
T3MEBuildNativeTopology[] := Module[
  {oldData, newData, oldT, newT, eps, projected, denominator, alphaCG,
   i,j,k,m,l,sourceTensorAgreement, cgStatus, u11,u12,u21,u22,
   common, topology, symCondition},
  oldData = PhysicalSU2InvariantBasis[
    {2,2,2,2}, {False,False,False,True}, {{1,2}}, {1,2}];
  newData = PhysicalSU2InvariantBasis[
    {2,2,2,2}, {False,False,False,False}, {{1,2},{3,4}}, {1,2}];
  If[!AssociationQ[oldData] || !AssociationQ[newData] ||
     oldData["PhysicalCount"] =!= 1 || newData["PhysicalCount"] =!= 1,
    Print["T3ME topology CG invariant-basis multiplicity mismatch"];
    Return[$Failed]
  ];
  oldT = Normal[oldData["Basis"][[1]]];
  newT = Normal[newData["Basis"][[1]]];
  If[Dimensions[oldT] =!= {2,2,2,2} ||
     Dimensions[newT] =!= {2,2,2,2}, Return[$Failed]];
  eps = {{0,1},{-1,0}};
  projected = Table[
    Sum[oldT[[i,j,k,l]] eps[[l,m]], {l,2}],
    {i,2},{j,2},{k,2},{m,2}];
  denominator = Total[Flatten[Conjugate[newT]] Flatten[newT]];
  If[TrueQ[PossibleZeroQ[denominator]], Return[$Failed]];
  alphaCG = FullSimplify[
    Total[Flatten[Conjugate[newT]] Flatten[projected]] / denominator];
  sourceTensorAgreement = TrueQ[FullSimplify[projected == alphaCG newT]];
  symCondition = TrueQ[
    FullSimplify[newT == Transpose[newT,{1,2,4,3}]]];
  If[!sourceTensorAgreement || !symCondition ||
     TrueQ[PossibleZeroQ[alphaCG]],
    Print["Original HH S1 BarS2 -> native HHXY tensor identity FAILED; ",
      "proportionality=", InputForm[alphaCG],
      ", equality=", sourceTensorAgreement,
      ", scalar-pair symmetric=", symCondition];
    Return[$Failed]
  ];
  cgStatus = CheckAbort[
    Check[DefineCG[T3MENativeHHMassCG, newData["CGReps"],newT]; True,
      $Failed],$Aborted];
  If[cgStatus =!= True, Return[$Failed]];
  {u11,u12,u21,u22} = {
    $T3MEInteractionMixingMatrix[[1,1]],
    $T3MEInteractionMixingMatrix[[1,2]],
    $T3MEInteractionMixingMatrix[[2,1]],
    $T3MEInteractionMixingMatrix[[2,2]]};
  common = lambdaT3[] H[i] H[j] CG[T3MENativeHHMassCG,{i,j,k,m}];
  (* The minus sign follows Bar[S2] = -eps (eps Bar[S2]);
     the cross coefficient has BOTH distinguishable X-Y contractions. *)
  topology = PlusHc[-alphaCG common (
    u11 u21 NewScalarX[k] NewScalarX[m] +
    (u11 u22 + u12 u21) NewScalarX[k] NewScalarY[m] +
    u12 u22 NewScalarY[k] NewScalarY[m])];
  If[!T3MEOldFieldsAbsentQ[topology] ||
     !T3MECouplingPresentQ[topology,lambdaT3] ||
     !T3MECouplingPresentQ[topology,T3MEU11] ||
     !T3MECouplingPresentQ[topology,T3MEU12] ||
     !T3MECouplingPresentQ[topology,T3MEU21] ||
     !T3MECouplingPresentQ[topology,T3MEU22],
    Return[$Failed]];
  <|"Expression" -> (topology // Expand // RelabelIndices),
    "InvariantCG" -> "T3MENativeHHMassCG",
    "OriginalTimesEpsilonCGFactor" -> ToString[alphaCG,InputForm],
    "OriginalTimesEpsilonTensorExact" -> sourceTensorAgreement,
    "SymmetricMassScalarCG" -> symCondition,
    "OriginalHHScalarContractions" -> 1,
    "NativeMassScalarContractions" -> 1,
    "MappedRGBetaQuarticNormalizationCertified" -> False|>
];

(* Remove only the old lambdaT3 sector BEFORE doing generic field transport.
   No other baseline or extra quartic source is removed. Strictly reconstruct
   the original expression so that the baseline inventory cannot be weakened. *)
T3MESeparateOriginalTopology[expr_] := Module[
  {expanded,terms,selected,remainder},
  expanded = Expand[expr];
  terms = If[Head[expanded] === Plus,List@@expanded,{expanded}];
  selected = Select[terms,TrueQ[T3MECouplingPresentQ[#,lambdaT3]] &];
  If[selected === {} || Length[selected] >= Length[terms],Return[$Failed]];
  (* Multiset-safe: reconstruct by SELECTING all non-topology terms in
     order, never by using Complement on possibly repeated terms. *)
  remainder = Total[Select[terms,!TrueQ[T3MECouplingPresentQ[#,lambdaT3]] &]];
  If[!TrueQ[Expand[Total[selected]+remainder-expanded]===0] ||
     T3MECouplingPresentQ[remainder,lambdaT3] ||
     !T3MECouplingPresentQ[Total[selected],lambdaT3],Return[$Failed]];
  <|"OriginalTopology"->Total[selected],"RemainingInteractions"->remainder,
    "OriginalTopologyExpandedTermCount"->Length[selected]|>
];

(* Portal2 H^dag H S2^dag S2 contains TWO S2 legs. Transporting each
   independently leaves T3HiggsPortal2CGInvI * epsilon * epsilon in the
   Matchete expression; it is invariant as a product, but Matchete does not
   interpret it as a single invariant tensor for the four external fields.

   Reconstruct BOTH independent portals as native four-field invariants.
   Old slot orientations: {Hbar,H,Bar[S2],S2} -> {T,F,T,F}.
   New slot orientations: {Hbar,H,X,Bar[X]} -> {T,F,F,T}.
   Each old tensor is contracted with epsilon on the last two slots, and
   the resulting tensor is registered directly against the new CGReps.
   No tensor normalization or phase is silently guessed. Projected rank,
   span and four-index orientation are verified exactly. *)
T3MEBuildNativeH2Portals[] := Module[
  {oldData,newData,eps,oldTensors,newBasis,projected,cgNames,
   originalRank,newRank,unionRank,register,coefficients,
   i,j,k,l,m,n,fields,terms,inv,expression},
  oldData=PhysicalSU2InvariantBasis[
    {2,2,2,2},{True,False,True,False},{},{1,2}];
  newData=PhysicalSU2InvariantBasis[
    {2,2,2,2},{True,False,False,True},{},{1,2}];
  If[!AssociationQ[oldData] || !AssociationQ[newData] ||
     oldData["PhysicalCount"] =!= 2 ||
     newData["PhysicalCount"] =!= 2,
    Print["T3ME native H2 portal physical invariant count mismatch"];
    Return[$Failed]];
  eps={{0,1},{-1,0}};
  oldTensors=Normal /@ oldData["Basis"];
  newBasis=Normal /@ newData["Basis"];
  If[!And@@(Dimensions[#]==={2,2,2,2}& /@ Join[oldTensors,newBasis]),
    Return[$Failed]];
  projected=Table[
    Table[Sum[oldTensors[[inv,i,j,m,n]] eps[[m,k]] eps[[n,l]],
      {m,2},{n,2}],{i,2},{j,2},{k,2},{l,2}],
    {inv,2}];
  originalRank=MatrixRank[Flatten /@ projected];
  newRank=MatrixRank[Flatten /@ newBasis];
  unionRank=MatrixRank[Flatten /@ Join[newBasis,projected]];
  If[{originalRank,newRank,unionRank}=!={2,2,2},
    Print["Native portal2 projected tensors do not span native basis: ",
      InputForm[{originalRank,newRank,unionRank}]];
    Return[$Failed]];
  cgNames={T3MENativeH2CGInv1,T3MENativeH2CGInv2};
  register=Table[CheckAbort[
    Check[DefineCG[cgNames[[inv]],newData["CGReps"],projected[[inv]]];True,
      $Failed],$Aborted],{inv,2}];
  If[register =!= {True,True},
    Print["Native H2 portal CG tensor registration failed: ",InputForm[register]];
    Return[$Failed]];
  coefficients={
    $T3MEInteractionMixingMatrix[[2,1]],
    $T3MEInteractionMixingMatrix[[2,2]]};
  fields={NewScalarX,NewScalarY};
  terms=Table[
    Sum[coefficients[[a]] Bar[coefficients[[b]]] *
      fields[[a]][k] Bar[fields[[b]][l]],
      {a,2},{b,2}],{inv,2}];
  expression=Association@Table[
    With[{cg=cgNames[[inv]],coupling=If[inv===1,lambdaH2Inv1,lambdaH2Inv2]},
      SymbolName[coupling] -> (
        PlusHc[coupling[] Bar[H[i]] H[j] CG[cg,{i,j,k,l}] terms[[inv]]]/2
          // Expand // RelabelIndices)],
    {inv,2}];
  If[Length[expression]=!=2 ||
     !And@@(T3MECouplingPresentQ[expression[#],Symbol[#]]& /@
       {"lambdaH2Inv1","lambdaH2Inv2"}) ||
     !And@@(T3MEOldFieldsAbsentQ /@ Values[expression]),
    Print["Native H2 portal expression inventory failed"];
    Return[$Failed]];
  <|"Expressions"->expression,
    "OldInvariantCount"->oldData["PhysicalCount"],
    "NewInvariantCount"->newData["PhysicalCount"],
    "TransportedRank"->originalRank,
    "CombinedBasisRank"->unionRank,
    "ExactTensorProjectionVerified"->True,
    "PhysicalCouplingMapCertified"->False|>
];

(* Separate these real portal sectors without dropping any conjugate terms.
   Strict algebraic reconstruction is checked for the entire source. *)
T3MESeparateOriginalH2Portals[expr_] := Module[
  {expanded,terms,portalTerms,rest,required},
  expanded=Expand[expr];
  terms=If[Head[expanded]===Plus,List@@expanded,{expanded}];
  required={lambdaH2Inv1,lambdaH2Inv2};
  (* Use an explicit nested Function: unlike pattern replacement, this does
     not capture the source term or dynamically bound coupling symbol. *)
  portalTerms=Select[terms,
    Function[term,AnyTrue[required,T3MECouplingPresentQ[term,#]&]]];
  rest=Total[Select[terms,
    Function[term,!AnyTrue[required,T3MECouplingPresentQ[term,#]&]]]];
  If[portalTerms==={} ||
     !And@@(T3MECouplingPresentQ[Total[portalTerms],#]& /@ required) ||
     AnyTrue[required,T3MECouplingPresentQ[rest,#]&] ||
     !TrueQ[Expand[Total[portalTerms]+rest-expanded]===0],
    Return[$Failed]];
  <|"OriginalH2Portals"->Total[portalTerms],
    "RemainingInteractions"->rest,
    "SourceTermCount"->Length[portalTerms]|>
];

(* Native scalar-quartic recoupling for all three baseline sources that
   contain TWO or FOUR S2 legs. The mixing makes S2 an epsilon-conjugate
   linear combination of X,Y; leaving the source CG multiplied by bridge
   CGs can fail native Matchete's individual invariant-tensor checks.

   (Bar S1) S1 (Bar S2) S2 -> (Bar mass) mass mass (Bar mass)
   (Bar S2) S2 (Bar S2) S2 -> mass (Bar mass) mass (Bar mass).
   Contract each S2 index with the same canonical epsilon used in the
   numerical mass rotation, and register one new native CG for each source
   contraction. Check both new invariant spans exactly, never equate bare
   tensor normalization to a numerical RGBeta coupling. *)
T3MEBuildNativeS2ScalarQuartics[] := Module[
  {eps, oldCross, newCross, oldSelf, newSelf,
   crossOld, crossNew, crossProjected, crossRank, crossNewRank, crossUnionRank,
   selfOld, selfNew, selfProjected, selfRank, selfNewRank, selfUnionRank,
   crossCG, selfCG, registers, fields, u1, u2, crossPolynomial,
   selfPolynomial, crossTerms, selfTerm, terms, q,i,j,k,l,m,n,o,p,
   a,b,c,d,expected},
  eps={{0,1},{-1,0}};
  oldCross=PhysicalSU2InvariantBasis[
    {2,2,2,2},{True,False,True,False}];
  newCross=PhysicalSU2InvariantBasis[
    {2,2,2,2},{True,False,False,True}];
  oldSelf=PhysicalSU2InvariantBasis[
    {2,2,2,2},{True,False,True,False},{{1,3},{2,4}}];
  newSelf=PhysicalSU2InvariantBasis[
    {2,2,2,2},{False,True,False,True},{{1,3},{2,4}}];
  If[!And@@(AssociationQ /@ {oldCross,newCross,oldSelf,newSelf}) ||
     {oldCross["PhysicalCount"],newCross["PhysicalCount"],
      oldSelf["PhysicalCount"],newSelf["PhysicalCount"]}=!={2,2,1,1},
    Print["Native S2 scalar quartic invariant multiplicities differ from source"];
    Return[$Failed]];
  crossOld=Normal /@ oldCross["Basis"];
  crossNew=Normal /@ newCross["Basis"];
  selfOld=Normal /@ oldSelf["Basis"];
  selfNew=Normal /@ newSelf["Basis"];
  If[!And@@(Dimensions[#]==={2,2,2,2}& /@
        Join[crossOld,crossNew,selfOld,selfNew]), Return[$Failed]];
  (* Two epsilon contractions on the two source S2 legs. *)
  crossProjected=Table[
    Table[Sum[crossOld[[q,i,j,m,n]] eps[[m,k]] eps[[n,l]],
       {m,2},{n,2}],{i,2},{j,2},{k,2},{l,2}],
    {q,2}];
  (* Four epsilon contractions for S2 self, preserving both boson pairs. *)
  selfProjected=Table[
    Sum[selfOld[[1,m,n,o,p]] eps[[m,i]] eps[[n,j]]
      eps[[o,k]] eps[[p,l]],{m,2},{n,2},{o,2},{p,2}],
    {i,2},{j,2},{k,2},{l,2}];
  crossRank=MatrixRank[Flatten /@ crossProjected];
  crossNewRank=MatrixRank[Flatten /@ crossNew];
  crossUnionRank=MatrixRank[Flatten /@ Join[crossNew,crossProjected]];
  selfRank=MatrixRank[{Flatten[selfProjected]}];
  selfNewRank=MatrixRank[Flatten /@ selfNew];
  selfUnionRank=MatrixRank[Flatten /@ Append[selfNew,selfProjected]];
  If[{crossRank,crossNewRank,crossUnionRank,selfRank,selfNewRank,
      selfUnionRank}=!={2,2,2,1,1,1},
    Print["Native S2 scalar quartic epsilon transport rank failure: ",
      InputForm[{crossRank,crossNewRank,crossUnionRank,selfRank,selfNewRank,
        selfUnionRank}]];
    Return[$Failed]];
  If[TrueQ[$T3BCFullResidualOnly],
    $T3BCNativeCrossProjected=crossProjected];
  crossCG={T3MENativeCrossCGInv1,T3MENativeCrossCGInv2};
  selfCG=T3MENativeS2SelfCGInv1;
  registers=Join[
    Table[CheckAbort[
       Check[DefineCG[crossCG[[q]],newCross["CGReps"],
         crossProjected[[q]]];True,$Failed],$Aborted],{q,2}],
    {CheckAbort[Check[DefineCG[selfCG,newSelf["CGReps"],
      selfProjected];True,$Failed],$Aborted]}];
  If[registers =!= {True,True,True},
    Print["Native S2 scalar quartic CG registration failure: ",
      InputForm[registers]];Return[$Failed]];
  fields={NewScalarX,NewScalarY};
  u1=$T3MEInteractionMixingMatrix[[1]];
  u2=$T3MEInteractionMixingMatrix[[2]];
  (* Each old S2 or Bar[S2] contributes minus epsilon, so the paired
     signs cancel in cross and in S2-self operators. *)
  crossPolynomial=Sum[
    Bar[u1[[a]]] u1[[b]] u2[[c]] Bar[u2[[d]]] *
    Bar[fields[[a]][i]] fields[[b]][j] *
    fields[[c]][k] Bar[fields[[d]][l]],
    {a,2},{b,2},{c,2},{d,2}];
  selfPolynomial=Sum[
    u2[[a]] Bar[u2[[b]]] u2[[c]] Bar[u2[[d]]] *
    fields[[a]][i] Bar[fields[[b]][j]] *
    fields[[c]][k] Bar[fields[[d]][l]],
    {a,2},{b,2},{c,2},{d,2}];
  crossTerms=Table[
    With[{coupling=If[q===1,lambda12Inv1,lambda12Inv2],
          cg=crossCG[[q]]},
      SymbolName[coupling] -> (coupling[]
        CG[cg,{i,j,k,l}] crossPolynomial // Expand // RelabelIndices)],
    {q,2}];
  selfTerm=(lambdaS2Inv1[]/2 CG[selfCG,{i,j,k,l}] selfPolynomial)
    // Expand // RelabelIndices;
  terms=Append[Association[crossTerms],"lambdaS2Inv1"->selfTerm];
  expected={"lambda12Inv1","lambda12Inv2","lambdaS2Inv1"};
  If[Sort[Keys[terms]]=!=Sort[expected] ||
     !And@@(T3MECouplingPresentQ[terms[#],Symbol[#]]& /@ expected) ||
     !And@@(T3MEOldFieldsAbsentQ /@ Values[terms]),
    Print["Native S2 scalar quartic coupling or field inventory failed"];
    Return[$Failed]];
  <|"Expressions"->terms,
    "CrossScalar"-><|"OldInvariantCount"->2,"NewInvariantCount"->2,
      "TransportedRank"->crossRank,"CombinedBasisRank"->crossUnionRank,
      "ExactTensorProjectionVerified"->True,
      "PhysicalCouplingMapCertified"->False|>,
    "S2Self"-><|"OldInvariantCount"->1,"NewInvariantCount"->1,
      "TransportedRank"->selfRank,"CombinedBasisRank"->selfUnionRank,
      "ExactTensorProjectionVerified"->True,
      "PhysicalCouplingMapCertified"->False|>|>
];

T3MESeparateOriginalS2ScalarQuartics[expr_] := Module[
  {required,expanded,parts,selected,rest},
  required={lambda12Inv1,lambda12Inv2,lambdaS2Inv1};
  expanded=Expand[expr];
  parts=If[Head[expanded]===Plus,List@@expanded,{expanded}];
  selected=Select[parts,
    Function[term,AnyTrue[required,T3MECouplingPresentQ[term,#]&]]];
  rest=Total[Select[parts,
    Function[term,!AnyTrue[required,T3MECouplingPresentQ[term,#]&]]]];
  If[selected==={} ||
     !And@@(T3MECouplingPresentQ[Total[selected],#]& /@ required) ||
     AnyTrue[required,T3MECouplingPresentQ[rest,#]&] ||
     !TrueQ[Expand[Total[selected]+rest-expanded]===0],
    Print["Native S2 scalar-quartic source separation inventory failed"];
    Return[$Failed]];
  <|"OriginalQuartics"->Total[selected],
    "RemainingInteractions"->rest,"SourceTermCount"->Length[selected]|>
];


(* The second T3 Yukawa y2 couples Bar[l] F Bar[S2]. The Bar[S2]
   SU(2) index is pseudoreal-transported through epsilon to an UNBARRED
   mass-eigen doublet. An old T3Y2CG * bridge contraction cannot be handed
   unchanged to CheckLagrangian: construct its native 2- or 3-leg invariant.

   B_-1: dF=1, original {CRep[1],CRep[1]} -> {CRep[1],1}.
   C_-1: dF=3, original {CRep[1],2,CRep[1]} -> {CRep[1],2,1}.
   The proportionality factor alpha is derived from the native tensor
   components and checked entrywise. No common sign is assumed across B/C.
   The explicit minus sign below follows Bar[S2] = -eps * Phi2. *)
T3MEBuildNativeYukawa2[model_Association] := Module[
  {dF,oldData,newData,oldT,newT,eps,projected,denominator,
   alphaCG,exact,registered,i,j,k,m,p,r,fermion,labels,term,
   u21,u22,rank,sourceOrientation},
  dF=model["Fermion","SU2"];
  If[!MemberQ[{1,3},dF],Return[$Failed]];
  oldData=PhysicalSU2InvariantBasis[
    {2,dF,2},{True,False,True},{},{1}];
  newData=PhysicalSU2InvariantBasis[
    {2,dF,2},{True,False,False},{},{1}];
  If[!AssociationQ[oldData] || !AssociationQ[newData] ||
     oldData["PhysicalCount"] =!= 1 ||
     newData["PhysicalCount"] =!= 1,
    Print["Native y2 SU(2) invariant multiplicity failed for dF=",dF];
    Return[$Failed]];
  oldT=Normal[oldData["Basis"][[1]]];
  newT=Normal[newData["Basis"][[1]]];
  eps={{0,1},{-1,0}};
  projected=If[dF===1,
    Table[Sum[oldT[[i,m]] eps[[m,k]],{m,2}],{i,2},{k,2}],
    Table[Sum[oldT[[i,j,m]] eps[[m,k]],{m,2}],
      {i,2},{j,3},{k,2}]];
  If[Dimensions[projected] =!= Dimensions[newT],Return[$Failed]];
  denominator=Total[Flatten[Conjugate[newT]] Flatten[newT]];
  If[TrueQ[PossibleZeroQ[denominator]],Return[$Failed]];
  alphaCG=FullSimplify[
    Total[Flatten[Conjugate[newT]] Flatten[projected]]/denominator];
  exact=TrueQ[FullSimplify[projected==alphaCG newT]];
  If[!exact || TrueQ[PossibleZeroQ[alphaCG]],
    Print["Native Yukawa2 epsilon tensor identity failed: ",
      InputForm[alphaCG]];
    Return[$Failed]];
  registered=CheckAbort[Check[
    DefineCG[T3MENativeY2CG,newData["CGReps"],newT];True,$Failed],$Aborted];
  If[registered =!= True,Return[$Failed]];
  fermion=If[dF===1,NewFermion[r],NewFermion[j,r]];
  labels=If[dF===1,{i,k},{i,j,k}];
  u21=$T3MEInteractionMixingMatrix[[2,1]];
  u22=$T3MEInteractionMixingMatrix[[2,2]];
  term=PlusHc[-alphaCG y2[p,r] NCM[Bar[l[i,p]],fermion]
      CG[T3MENativeY2CG,labels]
      (u21 NewScalarX[k]+u22 NewScalarY[k])];
  term=Expand[term]//RelabelIndices;
  If[!T3MEOldFieldsAbsentQ[term] ||
     !T3MECouplingPresentQ[term,y2] ||
     !T3MECouplingPresentQ[term,T3MEU21] ||
     !T3MECouplingPresentQ[term,T3MEU22],
    Print["Native y2 reconstructed interaction inventory failed"];
    Return[$Failed]];
  <|"Expression"->term,
    "OldInvariantCount"->1,"NewInvariantCount"->1,
    "ExactTensorProjectionVerified"->exact,
    "OldTimesEpsilonCGFactor"->ToString[alphaCG,InputForm],
    "FermionSU2Dimension"->dF,
    "NativeCG"->"T3MENativeY2CG",
    "PhysicalCouplingMapCertified"->False|>
];

(* Account for every source y2 term, including its hermitian conjugate,
   before substituting the native invariant. Nothing else is discarded. *)
T3MESeparateOriginalYukawa2[expr_] := Module[
  {expanded,terms,selected,rest},
  expanded=Expand[expr];
  terms=If[Head[expanded]===Plus,List@@expanded,{expanded}];
  selected=Select[terms,TrueQ[T3MECouplingPresentQ[#,y2]]&];
  rest=Total[Select[terms,!TrueQ[T3MECouplingPresentQ[#,y2]]&]];
  If[selected==={} || Length[selected]>=Length[terms] ||
     !T3MECouplingPresentQ[Total[selected],y2] ||
     T3MECouplingPresentQ[rest,y2] ||
     !TrueQ[Expand[Total[selected]+rest-expanded]===0],
    Print["Original y2 source interaction separation incomplete"];
    Return[$Failed]];
  <|"OriginalYukawa2"->Total[selected],
    "RemainingInteractions"->rest,"SourceTermCount"->Length[selected]|>
];

(* Native mass-eigenstate reconstruction of the FIRST Yukawa.
   Unlike y2, y1 contains the original unbarred S1 doublet and therefore
   requires NO epsilon on its scalar leg. Reconstruct the full PlusHc vertex
   instead of rewriting a previously internalised Field inside a fermion
   expression; exact original/new three-leg tensors must agree entrywise. *)
T3MEBuildNativeYukawa1[model_Association] := Module[
  {dF,oldData,newData,oldT,newT,denominator,alphaCG,exact,registered,
   i,j,k,p,r,fermion,labels,u11,u12,term},
  dF=model["Fermion","SU2"];
  If[!MemberQ[{1,3},dF],Return[$Failed]];
  oldData=PhysicalSU2InvariantBasis[
    {2,dF,2},{True,True,False},{},{1}];
  newData=PhysicalSU2InvariantBasis[
    {2,dF,2},{True,True,False},{},{1}];
  If[!AssociationQ[oldData] || !AssociationQ[newData] ||
     oldData["PhysicalCount"] =!= 1 ||
     newData["PhysicalCount"] =!= 1,
    Print["Native y1 invariant multiplicity failed"];Return[$Failed]];
  oldT=Normal[oldData["Basis"][[1]]];
  newT=Normal[newData["Basis"][[1]]];
  If[Dimensions[oldT] =!= Dimensions[newT],Return[$Failed]];
  denominator=Total[Flatten[Conjugate[newT]] Flatten[newT]];
  If[TrueQ[PossibleZeroQ[denominator]],Return[$Failed]];
  alphaCG=FullSimplify[
    Total[Flatten[Conjugate[newT]] Flatten[oldT]]/denominator];
  exact=TrueQ[FullSimplify[oldT==alphaCG newT]];
  If[!exact || TrueQ[PossibleZeroQ[alphaCG]],
    Print["Native Yukawa1 tensor identity failed"];
    Return[$Failed]];
  registered=CheckAbort[Check[
    DefineCG[T3MENativeY1CG,newData["CGReps"],newT];True,$Failed],
    $Aborted];
  If[registered =!= True,Return[$Failed]];
  fermion=If[dF===1,CConj[NewFermion[r]],CConj[NewFermion[j,r]]];
  labels=If[dF===1,{i,k},{i,j,k}];
  u11=$T3MEInteractionMixingMatrix[[1,1]];
  u12=$T3MEInteractionMixingMatrix[[1,2]];
  term=PlusHc[alphaCG y1[p,r] NCM[Bar[l[i,p]],fermion]
    CG[T3MENativeY1CG,labels]
    (u11 NewScalarX[k]+u12 NewScalarY[k])];
  term=Expand[term]//RelabelIndices;
  If[!T3MEOldFieldsAbsentQ[term] ||
     !T3MECouplingPresentQ[term,y1] ||
     !T3MECouplingPresentQ[term,T3MEU11] ||
     !T3MECouplingPresentQ[term,T3MEU12],
    Print["Native y1 reconstructed vertex inventory failed"];
    Return[$Failed]];
  <|"Expression"->term,
    "OldInvariantCount"->1,"NewInvariantCount"->1,
    "ExactTensorProjectionVerified"->exact,
    "OldToNewCGFactor"->ToString[alphaCG,InputForm],
    "FermionSU2Dimension"->dF,
    "NativeCG"->"T3MENativeY1CG",
    "PhysicalCouplingMapCertified"->False|>
];

(* Keep every y1 source term, including h.c., accounted for before
   reconstruction. The source term inventory and all other couplings stay
   untouched. *)
T3MESeparateOriginalYukawa1[expr_] := Module[
  {expanded,terms,selected,rest},
  expanded=Expand[expr];
  terms=If[Head[expanded]===Plus,List@@expanded,{expanded}];
  selected=Select[terms,TrueQ[T3MECouplingPresentQ[#,y1]]&];
  rest=Total[Select[terms,!TrueQ[T3MECouplingPresentQ[#,y1]]&]];
  If[selected==={} || Length[selected]>=Length[terms] ||
     !T3MECouplingPresentQ[Total[selected],y1] ||
     T3MECouplingPresentQ[rest,y1] ||
     !TrueQ[Expand[Total[selected]+rest-expanded]===0],
    Print["Original y1 source interaction separation incomplete"];
    Return[$Failed]];
  <|"OriginalYukawa1"->Total[selected],
    "RemainingInteractions"->rest,"SourceTermCount"->Length[selected]|>
];

(* General direct-CG transformation of all SEVEN additional scalar quartics.
   Species are H,Hdag,S1,S1dag,S2,S2dag; the original source tensors
   were fully symmetrised in identical fields by
   T3ExtendedBosonicTensorBasis. The old S2,Bar[S2] legs are transported
   through the numerical convention epsilon={{0,1},{-1,0}}:
       S2    -> -conj(U2a) epsilon Bar[Ma]
       BarS2 -> -U2a epsilon Ma
       S1    -> U1a Ma
       BarS1 -> conj(U1a) Bar[Ma]
   with a= X,Y. Every source tensor is projected index by index and
   registered as ONE native Matchete four-leg CG. No composite source-CG
   times epsilon-CG survives in these quartics. No RG coefficient map is
   inferred from the Matchete tensor normalization. *)
T3MEExtraProjectedTensorSlot[tensor_,pos_Integer] := Module[
  {eps={{0,1},{-1,0}},i,j,k,l,old},
  Table[With[{idx={i,j,k,l}},
      Sum[Extract[tensor,ReplacePart[idx,pos->old]] *
        eps[[old,idx[[pos]]]],{old,2}]],
    {i,2},{j,2},{k,2},{l,2}]
];
T3MEExtraMassField[species_String,idx_] := Module[
  {u11=$T3MEInteractionMixingMatrix[[1,1]],
   u12=$T3MEInteractionMixingMatrix[[1,2]],
   u21=$T3MEInteractionMixingMatrix[[2,1]],
   u22=$T3MEInteractionMixingMatrix[[2,2]]},
  Switch[species,
    "H",H[idx],
    "Hdag",Bar[H[idx]],
    "S1",u11 NewScalarX[idx]+u12 NewScalarY[idx],
    "S1dag",Bar[u11] Bar[NewScalarX[idx]]+Bar[u12] Bar[NewScalarY[idx]],
    "S2",-Bar[u21] Bar[NewScalarX[idx]]-Bar[u22] Bar[NewScalarY[idx]],
    "S2dag",-u21 NewScalarX[idx]-u22 NewScalarY[idx],
    _, $Failed]
];
T3MEBuildNativeExtraQuartics[model_Association,extra_List] := Module[
  {expressions=<||>,certs=<||>,couplings={},names,counts,indexNumber,
   digits,coupling,cgName,source,oldData,newData,oldT,newBasis,
   newConjugated,positions,projected,rankNew,rankUnion,exact,
   cgReps,dims,smPositions,register,indices,product,term,entry},
  If[Length[extra] =!= 7 ||
     Length[DeleteDuplicates[Lookup[extra,"Name",{}]]] =!= 7,
    Return[$Failed]];
  Do[
    entry=extra[[q]];
    counts=entry["Counts"];
    If[!ListQ[counts] || Length[counts] =!= 6 || Total[counts] =!= 4,
      Return[$Failed]];
    names=T3ExtendedSpeciesList[counts];
    If[Length[names] =!= 4,Return[$Failed]];
    digits=StringJoin[ToString/@counts];
    indexNumber=Last[StringSplit[entry["Name"],"_"]];
    If[!StringMatchQ[indexNumber,DigitCharacter..],Return[$Failed]];
    coupling=Symbol["lambdaT3ExtraQ"<>digits<>"I"<>indexNumber];
    cgName=Symbol["T3MENativeExtraQ"<>digits<>"I"<>indexNumber<>"CG"];
    If[!T3MECouplingPresentQ[entry["Expression"],coupling],
      Print["Additional quartic source/coupling mismatch: ",entry["Name"]];
      Return[$Failed]];
    oldData=T3ExtendedBosonicTensorBasis[model,names];
    If[!AssociationQ[oldData],Return[$Failed]];
    source=FromDigits[indexNumber];
    If[!IntegerQ[source] || source<1 || source>Length[oldData["Basis"]],
      Return[$Failed]];
    oldT=Normal[oldData["Basis"][[source]]];
    dims=T3ExtendedFieldDimensions[model,names];
    newConjugated=MemberQ[{"Hdag","S1dag","S2"},#]& /@ names;
    smPositions=Flatten@Position[names,"H"|"Hdag"];
    newData=PhysicalSU2InvariantBasis[
      dims,newConjugated,{},smPositions];
    If[!AssociationQ[newData],Return[$Failed]];
    newBasis=Normal /@ newData["Basis"];
    If[Dimensions[oldT] =!= {2,2,2,2} ||
       !And@@(Dimensions[#]==={2,2,2,2}& /@ newBasis),
      Return[$Failed]];
    positions=Flatten@Position[names,"S2"|"S2dag"];
    projected=Fold[T3MEExtraProjectedTensorSlot[#1,#2]&,oldT,positions];
    rankNew=MatrixRank[Flatten/@newBasis];
    rankUnion=MatrixRank[Flatten/@Join[newBasis,{projected}]];
    exact=TrueQ[rankUnion===rankNew] &&
      !TrueQ[projected==ConstantArray[0,{2,2,2,2}]];
    If[!exact,
      Print["Additional quartic projected tensor is not in native SU(2) ",
        "invariant span: ",entry["Name"]];
      Return[$Failed]];
    cgReps=newData["CGReps"];
    register=CheckAbort[Check[
      DefineCG[cgName,cgReps,projected];True,$Failed],$Aborted];
    If[register =!= True,
      Print["Failed native additional-quartic CG registration: ",
        entry["Name"]];Return[$Failed]];
    indices=Array[Unique["t3meextra"]&,4];
    product=Times@@MapThread[T3MEExtraMassField,{names,indices}];
    If[FreeQ[product,Verbatim[NewScalarX]|Verbatim[NewScalarY],
         {0,Infinity}],Return[$Failed]];
    term=PlusHc[coupling[] CG[cgName,indices] product];
    term=Expand[term]//RelabelIndices;
    If[!T3MEOldFieldsAbsentQ[term] ||
       !T3MECouplingPresentQ[term,coupling],Return[$Failed]];
    expressions[entry["Name"]]=term;
    certs[entry["Name"]]=<|
      "SourceCoupling"->SymbolName[coupling],
      "SourceCounts"->counts,
      "SourceTensorIndex"->source,
      "NativeCG"->SymbolName[cgName],
      "S2EpsilonProjectedLegs"->Length[positions],
      "NativeInvariantSpanRank"->rankNew,
      "CombinedBasisRank"->rankUnion,
      "ExactTensorProjectionVerified"->exact|>;
    AppendTo[couplings,coupling],
    {q,Length[extra]}];
  If[Length[expressions] =!= 7 || Length[certs] =!= 7 ||
     Length[DeleteDuplicates[couplings]] =!= 7 ||
     !And@@(T3MECouplingPresentQ[Total[Values[expressions]],#]& /@
       couplings),Return[$Failed]];
  <|"Expressions"->expressions,
    "Certificates"->certs,
    "SourceOperatorCount"->7,"RecoupledOperatorCount"->7,
    "AllProjectedTensorsInvariant"->True,
    "PhysicalCouplingMapCertified"->False|>
];

(* All seven quartic source terms are retained, including h.c.; a single
   missing source symbol or any unexplained remainder is fatal. *)
T3MESeparateOriginalExtraQuartics[expr_,symbols_List] := Module[
  {expanded,terms,selected,rest},
  If[Length[symbols] =!= 7 ||
     Length[DeleteDuplicates[symbols]] =!= 7,Return[$Failed]];
  expanded=Expand[expr];
  terms=If[Head[expanded]===Plus,List@@expanded,{expanded}];
  selected=Select[terms,Function[part,
    AnyTrue[symbols,T3MECouplingPresentQ[part,#]&]]];
  rest=Total[Select[terms,Function[part,
    !AnyTrue[symbols,T3MECouplingPresentQ[part,#]&]]]];
  If[selected==={} ||
     !And@@(T3MECouplingPresentQ[Total[selected],#]& /@ symbols) ||
     AnyTrue[symbols,T3MECouplingPresentQ[rest,#]&] ||
     !TrueQ[Expand[Total[selected]+rest-expanded]===0],
    Print["Additional quartic source expression separation incomplete"];
    Return[$Failed]];
  <|"OriginalExtraQuartics"->Total[selected],
    "RemainingInteractions"->rest,"SourceTermCount"->Length[selected]|>
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

T3MEIsolateNonInvariantSource[built_Association, free_, extra_List, nativeTopology_, nativeH2_Association, nativeS2_Association, nativeY2_Association, nativeY1_Association, nativeExtra_Association] := Module[
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
    converted = Which[
      TrueQ[T3MECouplingPresentQ[entry["Expression"],lambdaT3]],
        nativeTopology,
      TrueQ[T3MECouplingPresentQ[entry["Expression"],lambdaH2Inv1]],
        nativeH2["Expressions"]["lambdaH2Inv1"],
      TrueQ[T3MECouplingPresentQ[entry["Expression"],lambdaH2Inv2]],
        nativeH2["Expressions"]["lambdaH2Inv2"],
      TrueQ[T3MECouplingPresentQ[entry["Expression"],lambda12Inv1]],
        nativeS2["Expressions"]["lambda12Inv1"],
      TrueQ[T3MECouplingPresentQ[entry["Expression"],lambda12Inv2]],
        nativeS2["Expressions"]["lambda12Inv2"],
      TrueQ[T3MECouplingPresentQ[entry["Expression"],lambdaS2Inv1]],
        nativeS2["Expressions"]["lambdaS2Inv1"],
      TrueQ[T3MECouplingPresentQ[entry["Expression"],y2]],
        nativeY2["Expression"],
      TrueQ[T3MECouplingPresentQ[entry["Expression"],y1]],
        nativeY1["Expression"],
      KeyExistsQ[nativeExtra["Expressions"],Lookup[entry,"Name",""]],
        nativeExtra["Expressions"][entry["Name"]],
      True,T3MERotateInteraction[entry["Expression"]]
    ];
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

(* Fail-closed partition audit of the *saved* Matchete EOMEFT.  This does
   not modify the Lagrangian used for matching.  All partitions are defined
   on exact expanded terms, so their sum reconstructs the original EFT. *)
ClearAll[T3MEHermiticitySectorAudit, T3MEHermiticitySectorCheck];
T3MEHermiticitySectorCheck[expr_] := Module[{result, keys},
  result=CheckAbort[
    Quiet[CheckLagrangian[expr,DetailedOutput->True],
      {CheckLagrangian::Hermiticity}],$Aborted];
  keys={"Hermiticity","CanonicallyNormalized","HeavyMassBasis",
        "GaugeInvariant","AllObjectsDefined"};
  If[!AssociationQ[result],
    Return[<|"result"->ToString[Head[result],InputForm],
      "hermiticity"->"Unavailable"|>]];
  <|"result"->"DetailedCheckReturned", "hermiticity"->Lookup[result,"Hermiticity","Unavailable"],
    "canonically_normalized"->Lookup[result,"CanonicallyNormalized","Unavailable"],
    "heavy_mass_basis"->Lookup[result,"HeavyMassBasis","Unavailable"],
    "gauge_invariant"->Lookup[result,"GaugeInvariant","Unavailable"],
    "all_objects_defined"->Lookup[result,"AllObjectsDefined","Unavailable"]|>
];
T3MEHermiticitySectorAudit[eft_, dest_String] := Module[
 {expanded, terms, indexed, llss, quartics, rest, xx, yy, xy,
  qSelf, qPortal, groups, results=<||>, sector, testExpr, idx, text,
  audit, name, counts, originalCount, notHermitian, preview, failGroups,
  qLabels, qCouplingResults, qSelected,qOther,label},
 T3MECheckpoint["scalar-threshold-hermiticity-sector-audit-started"];
 expanded=Expand[eft];
 terms=If[Head[expanded]===Plus,List@@expanded,{expanded}];
 indexed=Map[Function[t,<|"term"->t,"text"->ToString[t,InputForm]|>],terms];
 llss=Select[indexed,StringContainsQ[#["text"],"GammaCC"]&];
 quartics=Select[indexed,Function[item,
   With[{st=item["text"]},
     !StringContainsQ[st,"GammaCC"] &&
     (StringContainsQ[st,"Field[NewScalarX, Scalar"] ||
      StringContainsQ[st,"Field[NewScalarY, Scalar"]) &&
     (StringCount[st,"Field[NewScalarX, Scalar"]+
      StringCount[st,"Field[NewScalarY, Scalar"]+
      StringCount[st,"Field[H, Scalar"])===4]]];
 rest=Select[indexed,Function[item,
   !MemberQ[llss,item] && !MemberQ[quartics,item]]];
 (* The XX, XY and YY sectors include each corresponding conjugate; this
    avoids treating one chirality in isolation as a Hermitian Lagrangian. *)
 xx=Select[llss,StringContainsQ[#["text"],"Field[NewScalarX, Scalar"] &&
   !StringContainsQ[#["text"],"Field[NewScalarY, Scalar"]&];
 yy=Select[llss,StringContainsQ[#["text"],"Field[NewScalarY, Scalar"] &&
   !StringContainsQ[#["text"],"Field[NewScalarX, Scalar"]&];
 xy=Select[llss,StringContainsQ[#["text"],"Field[NewScalarX, Scalar"] &&
   StringContainsQ[#["text"],"Field[NewScalarY, Scalar"]&];
 qSelf=Select[quartics,StringCount[#["text"],"Field[H, Scalar"]===0&];
 qPortal=Select[quartics,StringCount[#["text"],"Field[H, Scalar"]===2&];
 groups=<|
   "full_expanded"->indexed,
   "background_no_LLSS_or_BSM_quartics"->rest,
   "background_plus_all_LLSS"->Join[rest,llss],
   "background_plus_all_BSM_quartics"->Join[rest,quartics],
   "without_LLSS"->Join[rest,quartics],
   "without_BSM_quartics"->Join[rest,llss],
   "background_plus_XX_LLSS"->Join[rest,xx],
   "background_plus_XY_LLSS"->Join[rest,xy],
   "background_plus_YY_LLSS"->Join[rest,yy],
   "background_plus_scalar_self_quartics"->Join[rest,qSelf],
   "background_plus_Higgs_portal_quartics"->Join[rest,qPortal]
 |>;
 counts=<|"all_terms"->Length[indexed],"LLSS"->Length[llss],
   "BSM_quartics"->Length[quartics],"background"->Length[rest],
   "XX_LLSS"->Length[xx],"XY_LLSS"->Length[xy],
   "YY_LLSS"->Length[yy],"scalar_self_quartics"->Length[qSelf],
   "Higgs_portal_quartics"->Length[qPortal]|>;
 (* Do not silently classify away unrecognized LLSS/quartic sectors. *)
 If[Length[Join[llss,quartics,rest]] =!= Length[indexed],
   Return[<|"status"->"PartitionCountMismatch","counts"->counts|>]];
 Do[
   sector=groups[name];
   testExpr=If[sector==={},0,Total[sector[[All,"term"]]]];
   results[name]=Join[
     <|"term_count"->Length[sector]|>,
     T3MEHermiticitySectorCheck[testExpr]],
   {name,Keys[groups]}];
 (* Coupling-resolved test isolates which original scalar quartic label
    actually controls the failed Hermiticity check. These are pure diagnostic
    deletions, NEVER used as a matching Lagrangian. *)
 qLabels=Sort@DeleteDuplicates@Flatten[
   StringCases[#["text"],
     RegularExpression["Coupling\\[(lambda[A-Za-z0-9]+),"]->"$1"]& /@ qSelf];
 qCouplingResults=<||>;
 Do[
   qSelected=Select[qSelf,
     StringContainsQ[#["text"],"Coupling["<>label<>","]&];
   qOther=Select[qSelf,
     !StringContainsQ[#["text"],"Coupling["<>label<>","]&];
   If[qSelected =!= {},
     qCouplingResults[label]=<|
       "term_count"->Length[qSelected],
       "background_plus_coupling"->T3MEHermiticitySectorCheck[
         Total[Join[rest,qSelected][[All,"term"]]]]["hermiticity"],
       "background_plus_self_quartics_without_coupling"->
         T3MEHermiticitySectorCheck[
           Total[Join[rest,qOther][[All,"term"]]]]["hermiticity"]|>],
   {label,qLabels}];
 failGroups=Select[Keys[results],results[#]["hermiticity"]===False&];
 (* Only compact, bounded source excerpts. Never export the huge formatted
    Mathematica CheckLagrangian result or silently alter any operator. *)
 preview=AssociationMap[
   Function[k,Take[Map[StringTake[#["text"],UpTo[250]]&,
      If[k==="LLSS",llss,If[k==="BSM_quartics",quartics,rest]]],UpTo[3]]],
   {"LLSS","BSM_quartics","background"}];
 audit=<|
   "status"->"HermiticitySectorAuditCompletedNotMatching",
   "partition_conserves_exact_terms"->True,
   "counts"->counts,
   "sector_checks"->results,
   "scalar_self_quartic_coupling_isolation"->qCouplingResults,
   "sectors_failing_hermiticity"->failGroups,
   "sample_inputform_prefixes"->preview,
   "scalar_threshold_matching_performed"->False,
   "physical_sequential_C5_computed"->False|>;
 If[Quiet@Check[Export[dest,audit,"RawJSON"],$Failed]===$Failed,
   Return[<|"status"->"AuditExportFailed"|>]];
 T3MECheckpoint["scalar-threshold-hermiticity-sector-audit-returned"];
 Print["Scalar-threshold Hermiticity sector audit: ",
   ToString[<|"counts"->counts,"failing_sectors"->failGroups|>,InputForm]];
 Print["Scalar-threshold self-quartic coupling isolation: ",
   ToString[qCouplingResults,InputForm]];
 audit
];


(* Strictly identify the quartic-only spectator sector at this specific
   1-loop LLHH boundary, without redefining couplings or editing the
   original saved canonical tree. Every original term must belong to either
   the retained or excluded multiset. Exactly 122 spectator scalar self
   quartics are expected for B_-1 and C_-1; unexpected structures abort. *)
ClearAll[T3MELeadingC5OneLoopProjection];
T3MELeadingC5OneLoopProjection[eft_] := Module[
 {expanded, terms, records, spectator, kept, majors, higgsQuartics,
  background, rest, expr, total, retCount, expCount, sourceHashes},
 expanded=Expand[eft];
 terms=If[Head[expanded]===Plus,List@@expanded,{expanded}];
 records=Map[<|"term"->#,"text"->ToString[#,InputForm]|>&,terms];
 spectator=Select[records, Function[item,
   With[{text=item["text"]},
    !StringContainsQ[text,"GammaCC"] &&
    StringCount[text,"Field[H, Scalar"]===0 &&
    (StringCount[text,"Field[NewScalarX, Scalar"]+
     StringCount[text,"Field[NewScalarY, Scalar"])===4]]];
 kept=Select[records,!MemberQ[spectator,#]&];
 majors=Select[kept,StringContainsQ[#["text"],"GammaCC"]&];
 higgsQuartics=Select[kept,Function[item,
   With[{text=item["text"]},
    !StringContainsQ[text,"GammaCC"] &&
    StringCount[text,"Field[H, Scalar"]===2 &&
    (StringCount[text,"Field[NewScalarX, Scalar"]+
     StringCount[text,"Field[NewScalarY, Scalar"])===2]]];
 background=Select[kept,Function[item,
   !MemberQ[majors,item] && !MemberQ[higgsQuartics,item]]];
 expCount=If[modelName==="B_-1",20,24];
 retCount=Length[kept];
 expr=Total[kept[[All,"term"]]];
 If[Length[spectator]=!=122 || Length[majors]=!=expCount ||
    Length[higgsQuartics]=!=52 || Length[background]=!=21 ||
    retCount=!=21+52+expCount ||
    Length[spectator]+retCount=!=Length[terms] ||
    !TrueQ[Expand[Total[terms]-Total[spectator[[All,"term"]]]-expr]===0] ||
    FreeQ[expr,GammaCC] || FreeQ[expr,lambdaT3],
   Return[<|"status"->"RejectedUnexpectedOperatorInventory",
     "saved_term_count"->Length[terms],
     "excluded_count"->Length[spectator],"kept_count"->retCount,
     "LLSS_count"->Length[majors],"Higgs_quartic_count"->Length[higgsQuartics],
     "background_count"->Length[background]|>]];
 <|"status"->"ExactLeadingC5OperatorPartition",
   "original_saved_tree_unchanged"->True,
   "saved_term_count"->Length[terms],
   "excluded_count"->Length[spectator],
   "kept_count"->retCount,
   "LLSS_count"->Length[majors],
   "Higgs_quartic_count"->Length[higgsQuartics],
   "background_count"->Length[background],
   "candidate"->expr|>
];

(* Verify a possible loss of Hermitian representation at the F-first EOMEFT
   handoff.  At tree level through dimension five, F exchange creates LLSS;
   the renormalizable scalar interactions themselves are already present in
   the independently validated UV mass-eigenfield Lagrangian. Rebuild those
   interactions from the same native SU(2) source, but only permit replacing
   EOMEFT if their sum with the EXACT saved LLSS equals EOMEFT algebraically.
   This is an equality gate, not a projection onto a Hermitian subspace. *)
(* Workstream C: full, unmodified symbolic residual and coupling partition.
   This is a representation comparison, NOT tensor identity proof, EFT loop
   subtraction, renormalization or a physical C5 certificate.  Export the
   original and the constructed candidate in full, so no terms are inferred
   from bounded diagnostic prefixes.  The names of the coupling groups are
   taken from native ASTs, not from text search and never silently dropped. *)
ClearAll[T3BCExportFullEFT1Residual];
T3BCExportFullEFT1Residual[saved_, composite_, residual_, dest_String] := Module[
  {prefix, raw, termList, groups, category, part,
   groupRecords={}, sumGroups, files=<||>, file, key, ix=0,
   termsWithGamma, write, result},
  prefix=FileNameJoin[{DirectoryName[dest],FileBaseName[dest]}];
  raw=Expand[residual];
  termList=If[raw===0,{},If[Head[raw]===Plus,List@@raw,{raw}]];
  termsWithGamma=Count[termList,t_ /; !FreeQ[t,GammaCC]];
  If[termsWithGamma=!=0,
    Print["Full EFT1 residual unexpectedly contains GammaCC terms."];
    Return[$Failed]];
  write[label_String,expr_] := Module[{path},
    path=prefix<>"_"<>label<>".inputform.txt";
    If[Quiet@Check[Export[path,ToString[InputForm[expr]],"Text"],$Failed]===$Failed,
      Return[$Failed]];
    FileNameTake[path]];
  files["saved_eomeft"] = write["saved_eomeft_full",saved];
  files["reconstructed_candidate"] = write["reconstructed_candidate_full",composite];
  files["symbolic_residual"] = write["full_symbolic_residual",raw];
  If[MemberQ[Values[files],$Failed],Return[$Failed]];
  category[t_] := Module[{all,labels,fields,fieldCounts,fieldKey,couplingKey},
    all=DeleteDuplicates@Cases[t,
      HoldPattern[Coupling[s_Symbol, ___]] :> SymbolName[s],{0,Infinity}];
    labels=Sort@Select[all,StringStartsQ[#,"lambda"]&];
    couplingKey=Which[Length[labels]===0,"NO_QUARTIC_COUPLING",
      Length[labels]===1,First[labels],True,"MULTIPLE_QUARTIC_COUPLINGS"];
    fields=Cases[t,HoldPattern[Field[s_Symbol,Scalar,___]] :>
      SymbolName[s],{0,Infinity}];
    fieldCounts=Counts[fields];
    fieldKey=If[fieldCounts===<||>,"NO_SCALARS",StringRiffle[
      (#+ToString[fieldCounts[#]]& /@ Sort@Keys[fieldCounts]),"_"]];
    couplingKey<>"__"<>fieldKey];
  groups=GroupBy[termList,category];
  sumGroups=If[groups===<||>,0,Total[Total /@ Values[groups]]];
  If[!TrueQ[Expand[sumGroups-raw]===0],
    Print["Full residual coupling partition does not reconstitute source."];
    Return[$Failed]];
  Do[
    ix++;
    part=Total[groups[key]];
    file=write["residual_group_"<>IntegerString[ix,10,2],part];
    If[file===$Failed,Return[$Failed]];
    AppendTo[groupRecords,<|
      "structural_sector"->key,"additive_terms"->Length[groups[key]],
      "file"->file,
      "native_CG_labels"->Sort@DeleteDuplicates@Cases[part,
        HoldPattern[CG[s_Symbol,___]]:>SymbolName[s],{0,Infinity}],
      "field_H"->StringContainsQ[ToString[InputForm[part]],"Field[H, Scalar"],
      "field_X"->StringContainsQ[ToString[InputForm[part]],"Field[NewScalarX"],
      "field_Y"->StringContainsQ[ToString[InputForm[part]],"Field[NewScalarY"]|>],
    {key,Sort@Keys[groups]}];
  If[Total[Lookup[groupRecords,"additive_terms",{}]] =!= Length[termList],
    Return[$Failed]];
  result=<|
    "format"->"t3-bc-complete-eft1-reconstruction-residual-v1",
    "status"->"CompleteSymbolicResidualExportedEquivalenceUnproven",
    "model"->modelName,
    "original_vs_independent_candidate"->"saved_minus_reconstruction",
    "matching_scope"->"FullEFT",
    "contains_majorana_residual"->False,
    "residual_zero_in_matchete_expression_algebra"->TrueQ[raw===0],
    "residual_additive_terms"->Length[termList],
    "coupling_group_partition_conserved"->True,
    "coupling_groups"->groupRecords,
    "files"->files,
    "CG_and_index_tensor_equivalence_proven"->False,
    "physical_matching_certified"->False,
    "renormalized_finite_C5_computed"->False|>;
  If[Quiet@Check[Export[prefix<>"_full_residual_manifest.json",
                        result,"RawJSON"],$Failed]===$Failed,Return[$Failed]];
  result
];

(* Exactly compare the two registered S2-cross CG tensors with the component
   constraints derived independently by T3BCSU2ComponentEquivalence.py.
   All tensor entries are the SOURCE arrays passed to DefineCG above.
   The slots are (i,j,k,l) exactly as in native CrossCGInv1/2.
   No Matchete matching is executed in this diagnostic. *)
ClearAll[T3BCNativeCrossTensorCheck];
T3BCNativeCrossTensorCheck[tensors_,dest_String]:=Module[
  {a,b,t1,t2,a1,b1,a2,b2,span1,span2,required1,required2,
   nonzero1,nonzero2,dir,stem,report},
  If[Dimensions[tensors]=!={2,2,2,2,2},Return[$Failed]];
  t1=tensors[[1]];t2=tensors[[2]];
  a=Table[KroneckerDelta[i,j] KroneckerDelta[k,l],
    {i,2},{j,2},{k,2},{l,2}];
  b=Table[KroneckerDelta[i,k] KroneckerDelta[j,l],
    {i,2},{j,2},{k,2},{l,2}];
  a1=FullSimplify[t1[[1,1,2,2]]];
  b1=FullSimplify[t1[[1,1,1,1]]-a1];
  a2=FullSimplify[t2[[1,1,2,2]]];
  b2=FullSimplify[t2[[1,1,1,1]]-a2];
  span1=TrueQ[FullSimplify[t1==a1 a+b1 b]];
  span2=TrueQ[FullSimplify[t2==a2 a+b2 b]];
  required1=span1 && TrueQ[FullSimplify[a1==1]] &&
    TrueQ[FullSimplify[b1==-1]];
  required2=span2 && TrueQ[FullSimplify[a2==b2]];
  nonzero1=!TrueQ[FullSimplify[t1==ConstantArray[0,{2,2,2,2}]]];
  nonzero2=!TrueQ[FullSimplify[t2==ConstantArray[0,{2,2,2,2}]]];
  dir=DirectoryName[dest];stem=FileBaseName[dest];
  Export[FileNameJoin[{dir,stem<>"_native_cross_cg_arrays.inputform.txt"}],
     ToString[InputForm[tensors]],"Text"];
  report=<|
    "format"->"t3-bc-native-cross-cg-exact-component-check-v1",
    "model"->modelName,
    "source"->"Exact arrays passed to DefineCG for native S2 cross quartics",
    "component_order"->"CG[i,j,k,l]",
    "basis_A"->"KroneckerDelta[i,j] KroneckerDelta[k,l]",
    "basis_B"->"KroneckerDelta[i,k] KroneckerDelta[j,l]",
    "CrossInv1_coeff_A_inputform"->ToString[InputForm[a1]],
    "CrossInv1_coeff_B_inputform"->ToString[InputForm[b1]],
    "CrossInv2_coeff_A_inputform"->ToString[InputForm[a2]],
    "CrossInv2_coeff_B_inputform"->ToString[InputForm[b2]],
    "CrossInv1_in_singlet_span"->span1,
    "CrossInv2_in_singlet_span"->span2,
    "CrossInv1_exact_required_A_minus_B"->required1,
    "CrossInv2_required_equal_A_B_coefficients"->required2,
    "CrossInv1_nonzero"->nonzero1,
    "CrossInv2_nonzero"->nonzero2,
    "native_cross_conditions_met"->(required1 && required2 && nonzero1 && nonzero2),
    "physical_matching_certified"->False,
    "renormalized_finite_C5_computed"->False|>;
  If[Quiet@Check[Export[FileNameJoin[{dir,stem<>
     "_native_cross_cg_check.json"}],report,"RawJSON"],$Failed]===$Failed,
    Return[$Failed]];
  Print["T3ME_NATIVE_CG_CHECK: model=",modelName,
    "; CrossInv1=A-B: ",required1,
    "; CrossInv2 in span(A+B): ",required2,
    "; nonzero: ",{nonzero1,nonzero2}];
  report
];

ClearAll[T3MEVerifiedRenormalizableReconstruction];
T3MEVerifiedRenormalizableReconstruction[
  saved_, built_Association, extraParts_Association,
  nativeTopology_Association,nativeH2_Association,
  nativeS2_Association,nativeExtra_Association] := Module[
 {expanded,terms,majorana,rest,rotatedRest,renormFree,renorm,composite,
  residual,exact,relabeled,attempted=False,fullCheck=$Failed,detail=$Failed,
  qCount,clean,report},
 T3MECheckpoint["scalar-threshold-verified-renormalizable-reconstruction-started"];
 expanded=Expand[saved];
 terms=If[Head[expanded]===Plus,List@@expanded,{expanded}];
 majorana=Select[terms,!FreeQ[#,GammaCC]&];
 rest=Select[terms,FreeQ[#,GammaCC]&];
 qCount=If[modelName==="B_-1",20,24];
 If[Length[majorana]=!=qCount ||
    AnyTrue[majorana,FreeQ[#,Verbatim[NewScalarX]|Verbatim[NewScalarY],{0,Infinity}]&],
    Return[<|"status"->"UnexpectedMatchedMajoranaInventory",
      "majorana_terms"->Length[majorana],"expected"->qCount,
      "equality_proven"->False,"accepted"->False|>]];
 rotatedRest=T3MERotateInteraction[extraParts["RemainingInteractions"]];
 If[rotatedRest===$Failed ||
    !T3MEOldFieldsAbsentQ[rotatedRest],
    Return[<|"status"->"OriginalRenormalizableRotationFailed",
      "equality_proven"->False,"accepted"->False|>]];
 renormFree=CheckAbort[FreeLag[NewScalarX,NewScalarY],$Aborted];
 If[renormFree===$Failed || renormFree===$Aborted,
    Return[<|"status"->"RenormalizableFreeLagFailed",
      "equality_proven"->False,"accepted"->False|>]];
 renorm=Expand[built["LSM"]+renormFree+rotatedRest+
    Total[Values[nativeH2["Expressions"]]]+
    Total[Values[nativeS2["Expressions"]]]+
    Total[Values[nativeExtra["Expressions"]]]+
    nativeTopology["Expression"]]//RelabelIndices;
 If[!FreeQ[renorm,Verbatim[NewFermion],{0,Infinity}] ||
    FreeQ[renorm,Verbatim[NewScalarX],{0,Infinity}] ||
    FreeQ[renorm,Verbatim[NewScalarY],{0,Infinity}] ||
    !FreeQ[renorm,GammaCC],
    Return[<|"status"->"RenormalizableFieldInventoryInvalid",
      "equality_proven"->False,"accepted"->False|>]];
 composite=Expand[renorm+Total[majorana]]//RelabelIndices;
 (* First try literal equality, then strictly algebraic index relabelling.
    We intentionally do not accept numerical sample agreement, nor simply
    matching the scalar-quartic coupling count. *)
 residual=Expand[saved-composite];
  If[TrueQ[$T3BCFullResidualOnly],
    $T3BCFullResidualReport=T3BCExportFullEFT1Residual[
      saved,composite,residual,outPath];
    If[!AssociationQ[$T3BCFullResidualReport],
      T3MEFail["FullEFT1ResidualExportFailed"]];
    T3MECheckpoint["full-EFT1-symbolic-residual-exported"]];
 exact=TrueQ[residual===0];
 If[!exact,
   relabeled=CheckAbort[RelabelIndices[residual],$Aborted];
   exact=(relabeled=!= $Aborted &&
     (TrueQ[relabeled===0] ||
       TrueQ[TimeConstrained[Simplify[relabeled==0],15,False]]))];
 (* This is an INDEPENDENT CheckLagrangian test, even when the saved
    EOMEFT and reassembled source differ. It does not certify equality. *)
 attempted=True;
 fullCheck=CheckAbort[
   Quiet[CheckLagrangian[composite],{CheckLagrangian::Hermiticity}],$Aborted];
 detail=CheckAbort[
   Quiet[CheckLagrangian[composite,DetailedOutput->True],
     {CheckLagrangian::Hermiticity}],$Aborted];
 clean=TrueQ[fullCheck] && AssociationQ[detail] &&
   And@@(TrueQ[Lookup[detail,#,False]]& /@
     {"Hermiticity","CanonicallyNormalized","HeavyMassBasis",
      "GaugeInvariant","AllObjectsDefined"});
 report=<|
   "status"->If[clean && exact,"VerifiedEquivalentReconstructionAccepted",
     If[clean,"IndependentCanonicalReassemblyValidatedEOMEquivalenceUnproved",
       If[exact,"EquivalentReconstructionFailedMatcheteChecks",
         "IndependentCanonicalReassemblyFailedMatcheteChecks"]]],
   "matched_majorana_term_count"->Length[majorana],
   "renormalizable_tree_term_count"->Length[rest],
   "algebraic_equality_proven"->exact,
   "independent_matchete_check_attempted"->attempted,
   "independent_matchete_check_passed"->clean,
   "original_saved_EOMEFT_algebraically_equivalent"->exact,
   "independent_canonical_source_verified_only"->(clean && !exact),
   "independent_hermiticity"->If[AssociationQ[detail],
     Lookup[detail,"Hermiticity","Unavailable"],"Unavailable"],
   "residual_nonzero_inputform_prefix"->If[exact,"",
     StringTake[ToString[InputForm[residual]],UpTo[600]]],
   "physical_sequential_C5_computed"->False,
   "accepted"->(clean && exact)|>;
 If[clean,AssociateTo[report,"candidate"->composite]];
 T3MECheckpoint["scalar-threshold-verified-renormalizable-reconstruction-returned"];
 report
];


(* Standalone X-only diagnostic performed after successful native Y-only matching.
   NO running couplings, NO inherited Y-loop in the X source.  The full raw Y
   one-loop expression is separately preserved.  An EOMEFT failure is fatal
   for this X diagnostic: never substitute MatchedEFT without explicit proof. *)
ClearAll[T3BCXFollowupCheck, T3BCXFollowup];
$T3BCXRunMatch=(Environment["T3_BC_X_FOLLOWUP_MATCH"]==="1");
T3BCXFollowupCheck[lag_] := Module[{ordinary, detail},
 ordinary=CheckAbort[Quiet[CheckLagrangian[lag],
     {CheckLagrangian::Hermiticity}],$Aborted];
 detail=CheckAbort[Quiet[CheckLagrangian[lag,DetailedOutput->True],
     {CheckLagrangian::Hermiticity}],$Aborted];
 <|"ordinary_check"->TrueQ[ordinary],
   "hermiticity"->If[AssociationQ[detail],
      TrueQ[Lookup[detail,"Hermiticity",False]],False],
   "canonical"->If[AssociationQ[detail],
      TrueQ[Lookup[detail,"CanonicallyNormalized",False]],False],
   "heavy_mass_basis"->If[AssociationQ[detail],
      TrueQ[Lookup[detail,"HeavyMassBasis",False]],False],
   "gauge_invariant"->If[AssociationQ[detail],
      TrueQ[Lookup[detail,"GaugeInvariant",False]],False],
   "all_objects_defined"->If[AssociationQ[detail],
      TrueQ[Lookup[detail,"AllObjectsDefined",False]],False]|>
];
T3BCXFollowup[yTree_Association,yLoop_Association,dest_String]:= Module[
 {stem,base,dir,rawTree,rawLoop,treeText,loopText,diagnostic,
  freeLight,freeHeavy,freeDelta,freeDeltaZero,removed,defined,rebased,check,
  directHeavyCheck,directHeavyOK,rebasedOK,rebaseAudit,matchingInput,matchingInputLabel,
  xTree,xLoop,loopDifference,loopTextX,fullTextX,candidates,terms,
  xReport,output,sourceCheck,xFullCheck,matchedText,readOnlyTree,
  success,valid,matchReport,loopSectors,report},
 dir=DirectoryName[dest];base=FileBaseName[dest];
 stem=FileNameJoin[{dir,base<>"_x_followup"}];
 report=<|"format"->"t3-bc-native-x-only-followup-v1",
   "status"->"XFollowupNotStarted", "model"->modelName,
   "y_only_native_matchete_completed"->True,
   "source_is_y_tree_EOMEFT"->True,
   "y_loop_inherited_in_x_matching"->False,
   "running_threshold_couplings_applied"->False,
   "saved_fermion_EOMEFT_equivalence_certified"->False,
   "physical_matching_certified"->False,
   "two_threshold_physical_prediction"->False|>;
 T3MECheckpoint["x-followup-input-export-started"];
 If[!KeyExistsQ[yTree,"EOMEFT"]||!KeyExistsQ[yLoop,"EOMEFT"],
   AssociateTo[report,"status"->"XFollowupMissingYOnlyEOMEFT"];
   Export[stem<>".json",report,"RawJSON"];Return[report]];
 rawTree=yTree["EOMEFT"];rawLoop=yLoop["EOMEFT"];
 treeText=ToString[InputForm[rawTree]];
 loopText=ToString[InputForm[rawLoop]];
 Export[stem<>"_y_tree_eomeft.inputform.txt",treeText,"Text"];
 Export[stem<>"_y_full_eomeft.inputform.txt",loopText,"Text"];
 If[!StringContainsQ[treeText,"Field[NewScalarX"] ||
    StringContainsQ[treeText,"Field[NewScalarY"] ||
    StringContainsQ[treeText,"Field[NewFermion"],
   AssociateTo[report,"status"->"XFollowupYTreeEOMFieldInventoryInvalid"];
   Export[stem<>".json",report,"RawJSON"];Return[report]];
 T3MECheckpoint["x-followup-light-EFT2-check-started"];
 sourceCheck=T3BCXFollowupCheck[rawTree];
 AssociateTo[report,"light_x_eft2_detailed_check"->sourceCheck];
 (* The stage-after-EOM / stage-after-coupling-replacement comparison is
    diagnostic only. Never change the canonical EOMEFT continuation source. *)
 If[KeyExistsQ[yTree,"MatchedEFT"],
   AssociateTo[report,"light_x_matched_eft_comparator_check"->
     T3BCXFollowupCheck[yTree["MatchedEFT"]]]];
 If[!And@@(TrueQ/@Values[sourceCheck]),
   AssociateTo[report,"status"->"XFollowupYTreeEOMCheckFailed"];
   Export[stem<>".json",report,"RawJSON"];Return[report]];
 T3MECheckpoint["x-followup-X-free-lagrangian-rebase-started"];
 freeLight=CheckAbort[FreeLag[NewScalarX]//RelabelIndices,$Aborted];
 If[MemberQ[{$Failed,$Aborted},freeLight],
   AssociateTo[report,"status"->"XFollowupLightFreeLagFailed"];
   Export[stem<>".json",report,"RawJSON"];Return[report]];
 removed=Quiet@Check[RemoveField[NewScalarX];True,$Failed];
 defined=If[removed===True,Quiet@Check[
   DefineField[NewScalarX,Scalar,SelfConjugate->False,
     Indices->{T3BSMd2},Charges->{U1Y[-1/2]},Mass->{Heavy,MSX}];True,$Failed],$Failed];
 If[defined=!=True,
   AssociateTo[report,"status"->"XFollowupHeavyXRegistrationFailed"];
   Export[stem<>".json",report,"RawJSON"];Return[report]];
 freeHeavy=CheckAbort[FreeLag[NewScalarX]//RelabelIndices,$Aborted];
 If[MemberQ[{$Failed,$Aborted},freeHeavy],
   report=Join[report,<|"status"->"XFollowupHeavyFreeLagFailed"|>];
   Export[stem<>".json",report,"RawJSON"];Return[report]];
 (* Distinguish an actual heavy-registry failure from a bad free-lagrangian
    subtraction. No approximation and no MatchedEFT fallback are allowed.
    A different free Lagrangian makes the unchanged source INELIGIBLE for
    the X match, even if it passes CheckLagrangian under the heavy registry. *)
 freeDelta=CheckAbort[Expand[freeHeavy-freeLight]//RelabelIndices,$Aborted];
 freeDeltaZero=TrueQ[freeDelta===0];
 T3MECheckpoint["x-followup-direct-heavy-registry-check-started"];
 directHeavyCheck=T3BCXFollowupCheck[rawTree];
 directHeavyOK=And@@(TrueQ/@Values[directHeavyCheck]);
 rebased=CheckAbort[Expand[rawTree-freeLight+freeHeavy]//RelabelIndices,$Aborted];
 If[MemberQ[{$Failed,$Aborted},rebased] ||
    FreeQ[rebased,NewScalarX] || !FreeQ[rebased,NewScalarY],
   report=Join[report,<|"status"->"XFollowupRebasedFieldInventoryInvalid",
     "direct_x_heavy_eft2_detailed_check"->directHeavyCheck,
     "free_lagrangians_algebraically_identical"->freeDeltaZero|>];
   Export[stem<>".json",report,"RawJSON"];Return[report]];
 Export[stem<>"_x_heavy_rebased_eft2.inputform.txt",ToString[InputForm[rebased]],"Text"];
 T3MECheckpoint["x-followup-rebased-check-started"];
 check=T3BCXFollowupCheck[rebased];
 rebasedOK=And@@(TrueQ/@Values[check]);
 (* Export a separate flat audit even if a subsequent X matching stage fails.
    The six individual Matchete flags are never hidden behind a nested
    optional report property. *)
 rebaseAudit=<|
   "status"->"XHeavyRegistryVersusRebaseComparisonCompleted",
   "model"->modelName,
   "source_light_checks"->sourceCheck,
   "direct_heavy_checks"->directHeavyCheck,
   "rebased_heavy_checks"->check,
   "source_light_check_passed"->And@@(TrueQ/@Values[sourceCheck]),
   "direct_heavy_check_passed"->directHeavyOK,
   "rebased_heavy_check_passed"->rebasedOK,
   "free_lagrangians_algebraically_identical"->freeDeltaZero,
   "free_light_and_heavy_structurally_identical"->TrueQ[freeLight===freeHeavy],
   "free_difference_inputform_prefix"->If[freeDelta===$Aborted,"ABORTED",
     StringTake[ToString[InputForm[freeDelta]],UpTo[1200]]],
   "rebase_minus_source_inputform_prefix"->StringTake[
     ToString[InputForm[Expand[rebased-rawTree]]],UpTo[1200]],
   "matched_eft_fallback_used"->False,
   "full_sequential_matching_certified"->False|>;
 If[Quiet@Check[Export[stem<>"_rebase_crosschecks.json",rebaseAudit,"RawJSON"],$Failed]===$Failed,
   report=Join[report,<|"status"->"XFollowupDiagnosticExportFailed"|>];
   Export[stem<>".json",report,"RawJSON"];Return[report]];
 report=Join[report,<|
   "rebased_x_heavy_eft2_detailed_check"->check,
   "direct_x_heavy_eft2_detailed_check"->directHeavyCheck,
   "exact_rebase_power_counting_performed"->True,
   "source_check_passed"->And@@(TrueQ/@Values[sourceCheck]),
   "direct_heavy_check_passed"->directHeavyOK,
   "rebased_check_passed"->rebasedOK,
   "free_lagrangians_algebraically_identical"->freeDeltaZero,
   "rebase_crosschecks_file"->FileNameTake[stem<>"_rebase_crosschecks.json"]|>];
 (* Matchete is allowed to match only a source with a verified heavy-field
    registry AND a justified free Lagrangian. Rebase is preferred; the
    unchanged exact EOMEFT is eligible only if freeLight==freeHeavy proven. *)
 matchingInput=rebased;
 matchingInputLabel="heavy_free_rebased_y_tree_eomeft";
 If[!rebasedOK,
   If[directHeavyOK && freeDeltaZero,
     matchingInput=rawTree;
     matchingInputLabel="unchanged_y_tree_eomeft_proven_equal_free_lagrangians",
     report=Join[report,<|"status"->"XFollowupRebasedEFT2CheckFailed"|>];
     Export[stem<>".json",report,"RawJSON"];
     Print["X followup rebase FAIL: direct heavy check=",directHeavyOK,
       "; rebase check=",rebasedOK,"; exact free difference zero=",freeDeltaZero,
       "; detailed evidence: ",stem<>"_rebase_crosschecks.json"];
     Return[report]]];
 report=Join[report,<|"matching_input_selection"->matchingInputLabel,
   "selected_x_input_full_matchete_check_passed"->True,
   "x_heavy_check_satisfied"->True|>];
 If[!TrueQ[$T3BCXRunMatch],
   report=Join[report,<|"status"->"XFollowupEFT2XHeavyCheckPassedOnly"|>];
   Export[stem<>".json",report,"RawJSON"];Return[report]];
 report=Join[report,<|"status"->"XFollowupTreeMatchingInProgress",
   "x_matching_source_validation_passed"->True|>];
 Export[stem<>".json",report,"RawJSON"];
 T3MECheckpoint["x-followup-tree-matching-started"];
 xTree=CheckAbort[RunT3Matching[matchingInput,5,0],$Aborted];
 T3MECheckpoint["x-followup-tree-matching-returned"];

 If[!AssociationQ[xTree]||Lookup[xTree,"Status",""]=!="Success"||
   !KeyExistsQ[xTree,"MatchedEFT"],
   AssociateTo[report,"status"->"XFollowupTreeMatchingFailed",
     "matchete_tree_stage_status"->If[AssociationQ[xTree],
       Lookup[xTree,"Status","NoStatus"],ToString[Head[xTree],InputForm]]];
   Export[stem<>".json",report,"RawJSON"];Return[report]];
 report=Join[report,<|"status"->"XFollowupOneLoopMatchingInProgress",
   "x_tree_matching_completed"->True|>];
 Export[stem<>".json",report,"RawJSON"];
 T3MECheckpoint["x-followup-one-loop-matching-started"];
 xLoop=CheckAbort[RunT3Matching[matchingInput,5,1],$Aborted];
 T3MECheckpoint["x-followup-one-loop-matching-returned"];

 If[!AssociationQ[xLoop]||Lookup[xLoop,"Status",""]=!="Success"||
   !KeyExistsQ[xLoop,"MatchedEFT"],
   AssociateTo[report,"status"->"XFollowupOneLoopMatchingFailed",
     "matchete_one_loop_stage_status"->If[AssociationQ[xLoop],
       Lookup[xLoop,"Status","NoStatus"],ToString[Head[xLoop],InputForm]]];
   Export[stem<>".json",report,"RawJSON"];Return[report]];
 T3MECheckpoint["x-followup-operator-export-started"];
 matchedText=ToString[InputForm[xTree["MatchedEFT"]]];
 fullTextX=ToString[InputForm[xLoop["MatchedEFT"]]];
 loopDifference=Expand[xLoop["MatchedEFT"]-xTree["MatchedEFT"]];
 loopTextX=ToString[InputForm[loopDifference]];
 Export[stem<>"_x_tree_sm.inputform.txt",matchedText,"Text"];
 Export[stem<>"_x_full_sm.inputform.txt",fullTextX,"Text"];
 Export[stem<>"_x_loop_difference.inputform.txt",loopTextX,"Text"];
 If[StringContainsQ[fullTextX,"Field[NewScalarX"] ||
    StringContainsQ[fullTextX,"Field[NewScalarY"] ||
    StringContainsQ[fullTextX,"Field[NewFermion"] ||
    !StringContainsQ[loopTextX,"hbar"],
   AssociateTo[report,"status"->"XFollowupMatchedFieldInventoryOrLoopInvalid"];
   Export[stem<>".json",report,"RawJSON"];Return[report]];
 terms=With[{e=Expand[loopDifference]},If[Head[e]===Plus,List@@e,{e}]];
 candidates=Select[terms,Function[term,
   With[{txt=ToString[InputForm[term]]},
     StringContainsQ[txt,"GammaCC"] &&
     StringCount[txt,"Field[H, Scalar"]===2 &&
     StringCount[txt,"Field[l, Fermion"]===2 &&
     !StringContainsQ[txt,"Field[NewScalarX"] &&
     !StringContainsQ[txt,"Field[NewScalarY"]]]];
 Export[stem<>"_x_raw_llhh.inputform.txt",
   ToString[InputForm[If[candidates==={},0,Total[candidates]]]],"Text"];
 AssociateTo[report,"status"->"XFollowupYTreeOnlyXRawMatchingNotCertified",
   "x_tree_matching_completed"->True,"x_one_loop_matching_completed"->True,
   "x_raw_llhh_term_count"->Length[candidates],
   "x_one_loop_correction_term_count"->Length[terms]];
 T3MECheckpoint["x-followup-matching-exports-returned"];
 Export[stem<>".json",report,"RawJSON"];
 report
];

T3MEScalarThresholdPerform[] := Module[
 {built, sourceInt, extraQ, originals, oldPresent, status,
  b1,b2,removeA,removeB,xstatus,ystatus,free,rotated,lag,check,
  detailed,heavy,fullCount,baselineNames,oldYukawas,results,
  extraSymbols,extrasPresent,missingExtra,spurionStatus,missingSpurions,
  failureDetail, topoParts, nativeTopology, nonTopologyRotated, topoInventory,
  portalParts, nativeH2, otherRotated, portalsRotated,
  scalarParts, nativeS2, scalarRotated, yukawaParts, nativeY2,
  y1Parts, nativeY1, extraParts, nativeExtra, extraRotated,
  matchingSource, matchingResult, matchedEFT, eftPath,
  terms, llssTerms, llss, llssPath, xPresent, yPresent, fAbsent,
  xLLSS, yLLSS, oneLoopMatchingResult, oneLoopMatchedEFT,
  loopDifference, loopTerms, hardCandidates, hardRaw, hardPath,
  oneLoopEFTPath, oneLoopCorrectionPath, loopText, hardText,
  hasLoopHbar, fLoopAbsent, xLoopPresent, yLoopPresent, hardFAbsent,
  scalarEFT, treeText, treeResult, loopResult, treeMatched, loopMatched,
  correction, correctionText, completeText, deltaTerms, scalarCandidates,
  llhh, llhhText, scalarTreePath, scalarLoopPath, scalarCorrectionPath,
  scalarLLHHPath, coreCheck, detailedCore, noF, noOld,
  massX, massY, checkHeavyX, checkHeavyY,
  firstTreeLength, firstLLSSCount, treeX, treeY, firstStageOriginalCount,
  lightTreeDetail, originalTreeDetail, hermiticityTrial, hermiticityTrialDetail,
  hermiticityRewrite = "None", hermiticityAttempts = {}, treeHermiticityReport,
  hermiticityDiagnosticPath, hermiticitySectorPath, hermiticitySectorResult, stageTag, rewriteName, candidate,
  retainedFermionMass, pureBaseline, baselineFields,
  reconstruction, reconstructionProof, reconstructionPath, reconstructionClean,
  lightEigenFree, heavyEigenFree, rebasedSource, rebaseDelta, rebaseReport,
  matchAuditPath, stageSummary, stageDiagnostic, k, expr, snapText},
 If[!MatchQ[$FrontEnd,_FrontEndObject],
   T3MEFail["FrontEndUnavailable"]];
 If[Environment["T3_EXPERIMENTAL_GENERAL_QUARTICS"] === "1",
   T3MEFail["DuplicateExperimentalQuarticsOptIn"]];

 T3MECheckpoint["saved-EFT1-tree-source-preflight-started"];
 If[!FileExistsQ[treePath],
   T3MEFail["SavedEFT1TreeNotFound", treePath]];
 treeText = Quiet @ Check[Import[treePath,"Text"],$Failed];
 If[!StringQ[treeText] || StringLength[treeText]<1000 ||
    !StringContainsQ[treeText,"Field[NewScalarX"] ||
    !StringContainsQ[treeText,"Field[NewScalarY"] ||
    !StringContainsQ[treeText,"GammaCC"] ||
    StringContainsQ[treeText,"Field[NewFermion"],
   T3MEFail["SavedTreeEFT1ContentInvalid"]];
 (* This is a local Matchete artifact, not an arbitrary user program.
    Restrict obvious Wolfram evaluation/IO syntax before InputForm parsing. *)
 If[AnyTrue[{";","Get[","Import[","Export[","RunProcess[","Run[",
              "ToExpression[","Set[","SetDelayed[","Quit[","Exit[",
              "Remove[","SystemOpen[","URLRead["},StringContainsQ[treeText,#]&],
   T3MEFail["SavedEFT1TreeContainsUnsupportedExecutableSyntax"]];
 firstTreeLength=StringLength[treeText];
 T3MECheckpoint["saved-EFT1-tree-source-preflight-passed"];
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

 (* Retain F as Heavy; replace S1/S2 by Light mass eigen-doublets. Their
    kinetic/mass terms are diagonal under FreeLag; only F is integrated out. *)
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
 T3MECheckpoint["register-light-X-for-fermion-first-started"];
 xstatus=Check[
   DefineField[NewScalarX,Scalar,SelfConjugate->False,
     Indices->{T3BSMd2},Charges->{U1Y[-1/2]},Mass->{Light,MSX}];True,
   $Failed];
 T3MECheckpoint["register-light-Y-for-fermion-first-started"];
 ystatus=Check[
   DefineField[NewScalarY,Scalar,SelfConjugate->False,
     Indices->{T3BSMd2},Charges->{U1Y[-1/2]},Mass->{Light,MSY}];True,
   $Failed];
 If[xstatus =!= True || ystatus =!= True,
   T3MEFail["LightMassEigenfieldRegistrationFailed"]];
 T3MECheckpoint["fermion-heavy-two-scalars-light-FreeLag-started"];
 free=Quiet@Check[FreeLag[NewFermion,NewScalarX,NewScalarY],$Failed];
 If[free === $Failed,T3MEFail["DiagonalFreeLagFailed"]];
 T3MECheckpoint["native-topology-basis-started"];
 topoParts = T3MESeparateOriginalTopology[sourceInt];
 If[!AssociationQ[topoParts],
   T3MEFail["OriginalTopologyOperatorSeparationFailed"]];
 nativeTopology = T3MEBuildNativeTopology[];
 If[!AssociationQ[nativeTopology] ||
    !TrueQ[nativeTopology["OriginalTimesEpsilonTensorExact"]] ||
    !TrueQ[nativeTopology["SymmetricMassScalarCG"]],
   T3MEFail["NativeTopologyTensorConstructionFailed"]];
 T3MECheckpoint["native-H2-portals-started"];
 portalParts=T3MESeparateOriginalH2Portals[topoParts["RemainingInteractions"]];
 If[!AssociationQ[portalParts],
   T3MEFail["OriginalH2PortalSeparationFailed"]];
 nativeH2=T3MEBuildNativeH2Portals[];
 If[!AssociationQ[nativeH2] ||
    !TrueQ[nativeH2["ExactTensorProjectionVerified"]],
   T3MEFail["NativeH2PortalTensorConstructionFailed"]];
 T3MECheckpoint["native-S2-scalar-quartics-started"];
 scalarParts=T3MESeparateOriginalS2ScalarQuartics[
   portalParts["RemainingInteractions"]];
 If[!AssociationQ[scalarParts],
   T3MEFail["OriginalS2ScalarQuarticSeparationFailed"]];
 nativeS2=T3MEBuildNativeS2ScalarQuartics[];
 If[!AssociationQ[nativeS2] ||
    !TrueQ[nativeS2["CrossScalar"]["ExactTensorProjectionVerified"]] ||
    !TrueQ[nativeS2["S2Self"]["ExactTensorProjectionVerified"]],
   T3MEFail["NativeS2ScalarQuarticTensorConstructionFailed"]];
 T3MECheckpoint["native-Yukawa2-basis-started"];
 yukawaParts=T3MESeparateOriginalYukawa2[
   scalarParts["RemainingInteractions"]];
 If[!AssociationQ[yukawaParts],
   T3MEFail["OriginalYukawa2SeparationFailed"]];
 nativeY2=T3MEBuildNativeYukawa2[model];
 If[!AssociationQ[nativeY2] ||
    !TrueQ[nativeY2["ExactTensorProjectionVerified"]] ||
    nativeY2["OldInvariantCount"] =!= 1 ||
    nativeY2["NewInvariantCount"] =!= 1,
   T3MEFail["NativeYukawa2TensorConstructionFailed"]];
 T3MECheckpoint["native-Yukawa1-basis-started"];
 y1Parts=T3MESeparateOriginalYukawa1[yukawaParts["RemainingInteractions"]];
 If[!AssociationQ[y1Parts],
   T3MEFail["OriginalYukawa1SeparationFailed"]];
 nativeY1=T3MEBuildNativeYukawa1[model];
 If[!AssociationQ[nativeY1] ||
    !TrueQ[nativeY1["ExactTensorProjectionVerified"]] ||
    nativeY1["OldInvariantCount"] =!= 1 ||
    nativeY1["NewInvariantCount"] =!= 1,
   T3MEFail["NativeYukawa1TensorConstructionFailed"]];
 T3MECheckpoint["native-extra-quartics-started"];
 extraParts=T3MESeparateOriginalExtraQuartics[
   y1Parts["RemainingInteractions"],extraSymbols];
 If[!AssociationQ[extraParts],
   T3MEFail["OriginalExtraQuarticSeparationFailed"]];
 nativeExtra=T3MEBuildNativeExtraQuartics[model,extraQ];
 If[!AssociationQ[nativeExtra] ||
    !TrueQ[nativeExtra["AllProjectedTensorsInvariant"]] ||
    nativeExtra["SourceOperatorCount"] =!= 7 ||
    nativeExtra["RecoupledOperatorCount"] =!= 7,
   T3MEFail["NativeExtraQuarticTensorConstructionFailed"]];

 (* All source and native SU(2) CGs have now been defined in the new kernel.
    Remove the already-integrated F, retain its mass as a real coupling,
    and re-register both light EFT1 eigen-doublets as Heavy scalars. *)
 T3MECheckpoint["scalar-threshold-remove-integrated-F-started"];
 If[Quiet@Check[RemoveField[NewFermion];True,$Failed] =!= True,
   T3MEFail["IntegratedFermionDeregistrationFailed"]];
 If[Quiet@Check[DefineCoupling[MF,SelfConjugate->True];True,$Failed] =!= True,
   T3MEFail["IntegratedFermionMassCouplingRegistrationFailed"]];
 (* The canonical EOMEFT was serialized in a different Matchete kernel.
    Before changing Light to Heavy, check its Hermiticity with exactly the
    same light-X/Y field assignments used in the first matching stage.
    The absent F has already been removed, so the EFT is the correct
    field inventory for this independent check. *)
 T3MECheckpoint["saved-EFT1-tree-parse-started"];
 scalarEFT=CheckAbort[Quiet@Check[ToExpression[treeText,InputForm],$Failed],$Aborted];
 T3MECheckpoint["saved-EFT1-tree-parse-returned"];
 If[scalarEFT === $Failed || scalarEFT === $Aborted || !FreeQ[scalarEFT,NewFermion] ||
    FreeQ[scalarEFT,NewScalarX] || FreeQ[scalarEFT,NewScalarY] ||
    FreeQ[scalarEFT,GammaCC],
   T3MEFail["SavedEFT1TreeParseOrFieldPreflightFailed"]];
 T3MECheckpoint["scalar-threshold-light-EFT1-hermiticity-preflight-started"];
 lightTreeDetail=CheckAbort[
   Quiet[CheckLagrangian[scalarEFT,DetailedOutput->True],
     {CheckLagrangian::Hermiticity}],$Aborted];
 T3MECheckpoint["scalar-threshold-light-EFT1-hermiticity-preflight-returned"];
 Print["Saved canonical EFT1 before heavy reclassification: ",
   ToString[<|"Hermiticity"->If[AssociationQ[lightTreeDetail],
       Lookup[lightTreeDetail,"Hermiticity",Missing["Absent"]],$Failed],
     "CanonicallyNormalized"->If[AssociationQ[lightTreeDetail],
       Lookup[lightTreeDetail,"CanonicallyNormalized",Missing["Absent"]],$Failed]|>,InputForm]];
 (* Capture the actual scalar free Lagrangian while X and Y are still
    registered LIGHT.  The incoming EFT was serialized with this exact
    mass-power-counting assignment, so it must not be sent unchanged to
    Match after redefining those fields as Heavy. *)
 T3MECheckpoint["scalar-threshold-light-scalar-FreeLag-capture-started"];
 lightEigenFree=CheckAbort[FreeLag[NewScalarX,NewScalarY] // RelabelIndices,$Aborted];
 If[MemberQ[{$Aborted,$Failed},lightEigenFree],
   T3MEFail["ScalarThresholdLightFreeLagCaptureFailed"]];
 T3MECheckpoint["scalar-threshold-light-scalar-FreeLag-capture-returned"];
 T3MECheckpoint["scalar-threshold-reclassify-X-Y-started"];
 If[Quiet@Check[RemoveField[NewScalarX];True,$Failed] =!= True ||
    Quiet@Check[RemoveField[NewScalarY];True,$Failed] =!= True,
   T3MEFail["EFT1LightScalarDeregistrationFailed"]];
 checkHeavyX=Quiet@Check[
   DefineField[NewScalarX,Scalar,SelfConjugate->False,
     Indices->{T3BSMd2},Charges->{U1Y[-1/2]},Mass->{Light,MSX}];True,$Failed];
 checkHeavyY=Quiet@Check[
   DefineField[NewScalarY,Scalar,SelfConjugate->False,
     Indices->{T3BSMd2},Charges->{U1Y[-1/2]},Mass->{Heavy,MSY}];True,$Failed];
 If[checkHeavyX =!= True || checkHeavyY =!= True,
   T3MEFail["ScalarThresholdHeavyFieldRegistrationFailed"]];
 (* Re-registering the heavy fields does not retroactively rebuild free
    terms constructed in a kernel with Mass->{Light,...}.  Obtain the heavy
    free terms under the CURRENT registry, then replace the old light free
    terms only in the explicitly opted-in leading-C5 matching subtheory.
    Any mismatch is exposed by CheckLagrangian; this is NOT an equality claim
    between light and heavy power counting. *)
 T3MECheckpoint["scalar-threshold-heavy-scalar-FreeLag-capture-started"];
 heavyEigenFree=CheckAbort[FreeLag[NewScalarX,NewScalarY] // RelabelIndices,$Aborted];
 If[MemberQ[{$Aborted,$Failed},heavyEigenFree],
   T3MEFail["ScalarThresholdHeavyFreeLagCaptureFailed"]];
 T3MECheckpoint["scalar-threshold-heavy-scalar-FreeLag-capture-returned"];
 (* Save the diagnostic status of the full theory before any optional
    leading-C5 projection. This mode does not certify full-EFT Hermiticity. *)
 fullSavedTreeHermiticity=If[AssociationQ[lightTreeDetail],
   Lookup[lightTreeDetail,"Hermiticity","Unavailable"],"Unavailable"];
 leadingC5Partition=<|"status"->"NotRequested"|>;
 If[TrueQ[$T3MELeadingC5Only],
   T3MECheckpoint["scalar-threshold-leading-C5-sector-projection-started"];
   leadingC5Partition=T3MELeadingC5OneLoopProjection[scalarEFT];
   If[Lookup[leadingC5Partition,"status",""]=!=
       "ExactLeadingC5OperatorPartition",
     T3MEFail["LeadingC5SectorInventoryFailed",
       ToString[KeyDrop[leadingC5Partition,{"candidate"}],InputForm]]];
   scalarEFT=leadingC5Partition["candidate"];
   T3MECheckpoint["scalar-threshold-leading-C5-sector-projection-returned"];
   (* A power-counting-aware replacement of the two scalar free
      Lagrangians, not a replacement of interaction vertices or LLSS.
      Keep this scoped to the explicit restricted matching experiment. *)
   T3MECheckpoint["scalar-threshold-light-to-heavy-free-rebase-started"];
   rebaseDelta=Expand[heavyEigenFree-lightEigenFree] // RelabelIndices;
   rebasedSource=Expand[scalarEFT-lightEigenFree+heavyEigenFree] // RelabelIndices;
   If[rebasedSource===$Aborted || rebasedSource===$Failed ||
      FreeQ[rebasedSource,NewScalarX] || FreeQ[rebasedSource,NewScalarY] ||
      FreeQ[rebasedSource,GammaCC],
      T3MEFail["ScalarThresholdFreeLagRebaseInvalid"]];
   rebaseReport=<|
     "status"->"HeavyScalarFreeLagRebaseAttemptedNotYetVerifiedByMatch",
     "mode"->"LeadingC5OneLoopSectorOnly",
     "light_free_head"->ToString[Head[lightEigenFree],InputForm],
     "heavy_free_head"->ToString[Head[heavyEigenFree],InputForm],
     "free_lagrangians_structurally_identical"->TrueQ[lightEigenFree===heavyEigenFree],
     "power_counting_delta_is_zero"->TrueQ[rebaseDelta===0],
     "delta_contains_X"->!FreeQ[rebaseDelta,NewScalarX],
     "delta_contains_Y"->!FreeQ[rebaseDelta,NewScalarY],
     "delta_inputform_prefix"->StringTake[ToString[InputForm[rebaseDelta]],UpTo[1800]],
     "full_saved_eft_unchanged"->True,
     "physical_C5_computed"->False|>;
   Export[FileNameJoin[{DirectoryName[outPath],
      FileBaseName[outPath]<>"_free_lagrangian_rebase.json"}],
      rebaseReport,"RawJSON"];
   scalarEFT=rebasedSource;
   T3MECheckpoint["scalar-threshold-light-to-heavy-free-rebase-returned"];
   Print["Opt-in leading-C5 one-loop subtheory: excluded ",
     leadingC5Partition["excluded_count"],
     " spectator four-heavy-scalar quartics; retained LLSS=",
     leadingC5Partition["LLSS_count"],
     "; Higgs-bearing quartics=",leadingC5Partition["Higgs_quartic_count"],
     "; FULL EFT Hermiticity is NOT certified"];
 ];
  If[!TrueQ[$T3MELeadingC5Only],
    T3MECheckpoint["Y-only-full-EFT-free-rebase-started"];
    scalarEFT=Expand[scalarEFT-lightEigenFree+heavyEigenFree] // RelabelIndices;
    If[scalarEFT===$Failed || scalarEFT===$Aborted ||
       FreeQ[scalarEFT,NewScalarX] || FreeQ[scalarEFT,NewScalarY] ||
       FreeQ[scalarEFT,GammaCC],
      T3MEFail["YOnlyFullEFTFreeLagRebaseInvalid"]];
    T3MECheckpoint["Y-only-full-EFT-free-rebase-returned"]];
 T3MECheckpoint["scalar-threshold-EFT1-lagrangian-check-started"];
 (* Do not use Check[CheckLagrangian[...,],$Failed] here: Matchete may issue
    nonfatal diagnostic messages even when its returned check is True.
    Match the first threshold and production SafeStage use CheckAbort only. *)
 coreCheck=CheckAbort[CheckLagrangian[scalarEFT],$Aborted];
 T3MECheckpoint["scalar-threshold-ordinary-lagrangian-check-returned"];
 detailedCore=CheckAbort[
   CheckLagrangian[scalarEFT,DetailedOutput->True],$Aborted];
 T3MECheckpoint["scalar-threshold-detailed-lagrangian-check-returned"];
 Print["Scalar-threshold CheckLagrangian: ",InputForm[coreCheck]];
 (* Check whether only Matchete's representation of Hermitian operators
    changed during the EOMEFT/serialization transition. Only attempt
    documented equivalence-preserving simplification functions; NEVER
    insert HC terms or replace the EFT by (L + HC[L])/2.

    Each rewrite must pass the *entire* ordinary CheckLagrangian plus the
    detailed Hermiticity, canonical, heavy-mass and gauge checks. *)
 originalTreeDetail=detailedCore;
 If[!TrueQ[coreCheck] || !AssociationQ[detailedCore] ||
    !TrueQ[Lookup[detailedCore,"HeavyMassBasis",False]],
   If[AssociationQ[detailedCore] &&
      Lookup[detailedCore,"Hermiticity",Missing["Unknown"]] === False &&
      TrueQ[Lookup[detailedCore,"HeavyMassBasis",False]] &&
      TrueQ[Lookup[detailedCore,"CanonicallyNormalized",False]] &&
      TrueQ[Lookup[detailedCore,"GaugeInvariant",False]],
     Do[
       rewriteName=stageTag;
       T3MECheckpoint["scalar-threshold-hermiticity-equivalence-"<>rewriteName<>"-started"];
       candidate=CheckAbort[
         Quiet[Switch[rewriteName,
            "OriginalRelabel", RelabelIndices[scalarEFT],
            "ExpandOnly", Expand[scalarEFT],
            "ExpandRelabel", RelabelIndices[Expand[scalarEFT]],
            "HcSimplify", RelabelIndices[HcSimplify[scalarEFT]],
            "GreensSimplify", RelabelIndices[GreensSimplify[scalarEFT]],
            _, $Failed]],$Aborted];
       If[candidate===$Failed || candidate===$Aborted ||
          !FreeQ[candidate,Verbatim[NewFermion],{0,Infinity}] ||
          FreeQ[candidate,Verbatim[NewScalarX],{0,Infinity}] ||
          FreeQ[candidate,Verbatim[NewScalarY],{0,Infinity}] ||
          FreeQ[candidate,GammaCC],
         AppendTo[hermiticityAttempts,<|"method"->rewriteName,
           "status"->"RewriteFailedOrInvalid"|>];Continue[]];
       hermiticityTrial=CheckAbort[
         Quiet[CheckLagrangian[candidate],{CheckLagrangian::Hermiticity}],$Aborted];
       hermiticityTrialDetail=CheckAbort[
         Quiet[CheckLagrangian[candidate,DetailedOutput->True],
           {CheckLagrangian::Hermiticity}],$Aborted];
       AppendTo[hermiticityAttempts,<|"method"->rewriteName,
         "ordinary"->ToString[hermiticityTrial,InputForm],
         "hermiticity"->If[AssociationQ[hermiticityTrialDetail],
           Lookup[hermiticityTrialDetail,"Hermiticity","Unavailable"],"Unavailable"],
         "canonical"->If[AssociationQ[hermiticityTrialDetail],
           Lookup[hermiticityTrialDetail,"CanonicallyNormalized","Unavailable"],"Unavailable"]|>];
       If[TrueQ[hermiticityTrial] && AssociationQ[hermiticityTrialDetail] &&
          TrueQ[Lookup[hermiticityTrialDetail,"Hermiticity",False]] &&
          TrueQ[Lookup[hermiticityTrialDetail,"CanonicallyNormalized",False]] &&
          TrueQ[Lookup[hermiticityTrialDetail,"HeavyMassBasis",False]] &&
          TrueQ[Lookup[hermiticityTrialDetail,"GaugeInvariant",False]] &&
          TrueQ[Lookup[hermiticityTrialDetail,"AllObjectsDefined",False]],
          scalarEFT=candidate;
          coreCheck=hermiticityTrial;
          detailedCore=hermiticityTrialDetail;
          hermiticityRewrite=rewriteName;
          Break[]],
       {stageTag,{"OriginalRelabel","ExpandOnly","ExpandRelabel","HcSimplify","GreensSimplify"}}]
   ];
 ];
 (* If a genuine Hermiticity failure remains, compare the saved EOMEFT
    against the original UV renormalizable X/Y interactions plus exactly its
    matched LLSS. Only exact equality AND a full new Matchete check permit
    a canonical representation replacement. *)
 reconstruction=<|"status"->"NotAttemptedOriginalPasses","accepted"->False|>;
 (* No fall-through to native Match when the residual branch cannot run. *)
 If[TrueQ[$T3BCFullResidualOnly] &&
    (TrueQ[coreCheck] || !AssociationQ[detailedCore] ||
     Lookup[detailedCore,"Hermiticity",True]=!=False),
   T3MEFail["FullEFT1ResidualOriginalComparisonBranchUnavailable"]];
 reconstructionPath=FileNameJoin[{DirectoryName[outPath],
   FileBaseName[outPath]<>"_reconstruction_audit.json"}];
 If[!TrueQ[coreCheck] && AssociationQ[detailedCore] &&
    Lookup[detailedCore,"Hermiticity",True]===False,
   reconstruction=T3MEVerifiedRenormalizableReconstruction[
     scalarEFT,built,extraParts,nativeTopology,nativeH2,nativeS2,nativeExtra];
   reconstructionProof=KeyDrop[reconstruction,{"candidate"}];
   Quiet@Check[Export[reconstructionPath,reconstructionProof,"RawJSON"],$Failed];
   Print["Scalar-threshold reconstruction: ",
      ToString[KeyTake[reconstructionProof,{"status","algebraic_equality_proven",
        "independent_matchete_check_passed"}],InputForm]];
   If[TrueQ[$T3BCFullResidualOnly],
      If[!AssociationQ[$T3BCFullResidualReport] ||
         Lookup[$T3BCFullResidualReport,"model",""]=!=modelName,
         T3MEFail["FullEFT1ResidualReportMissing"]];
      With[{cgcheck=T3BCNativeCrossTensorCheck[
           $T3BCNativeCrossProjected,outPath]},
         If[!AssociationQ[cgcheck],
           T3MEFail["NativeCrossCGComponentCheckFailed"]]];
      Return[<|
        "status"->"EFT1FullResidualExportedBeforeMatching",
        "model"->modelName,
        "residual_manifest_file"->FileBaseName[outPath]<>
          "_full_residual_manifest.json",
        "original_eomeft_equal_to_independent_candidate"->
          TrueQ[Lookup[reconstruction,"algebraic_equality_proven",False]],
        "reconstruction_check_passed"->
          TrueQ[Lookup[reconstruction,"independent_matchete_check_passed",False]],
        "matching_performed"->False,"physical_matching_certified"->False|>]];
   If[TrueQ[Lookup[reconstruction,"accepted",False]],
     scalarEFT=reconstruction["candidate"];
     coreCheck=True;
     detailedCore=CheckAbort[
       CheckLagrangian[scalarEFT,DetailedOutput->True],$Aborted];
     If[!AssociationQ[detailedCore] ||
        !And@@(TrueQ[Lookup[detailedCore,#,False]]& /@
          {"Hermiticity","CanonicallyNormalized","HeavyMassBasis",
           "GaugeInvariant","AllObjectsDefined"}),
       T3MEFail["ReconstructionValidationNotStable"]];
     hermiticityRewrite="VerifiedRenormalizableReconstruction"]];
 (* If and only if the independently rebuilt full EFT passes EVERY ordinary
    Matchete check, allow an opt-in experimental source-defined candidate.
    It contains EXACT saved LLSS, with every reconstructed scalar quartic.
    This is NOT an algebraically certified replacement of EOMEFT. *)
 If[!TrueQ[coreCheck] && TrueQ[$T3MEIndependentCanonicalReassembly] &&
    AssociationQ[reconstruction] &&
    TrueQ[Lookup[reconstruction,"independent_matchete_check_passed",False]] &&
    KeyExistsQ[reconstruction,"candidate"],
   T3MECheckpoint["independent-canonical-full-EFT-candidate-started"];
   candidate=reconstruction["candidate"];
   If[FreeQ[candidate,GammaCC] ||
      FreeQ[candidate,Verbatim[NewScalarX],{0,Infinity}] ||
      FreeQ[candidate,Verbatim[NewScalarY],{0,Infinity}] ||
      !FreeQ[candidate,Verbatim[NewFermion],{0,Infinity}],
     T3MEFail["IndependentCanonicalReassemblyFieldInventoryFailed"]];
   scalarEFT=candidate;
   coreCheck=CheckAbort[CheckLagrangian[scalarEFT],$Aborted];
   detailedCore=CheckAbort[
     CheckLagrangian[scalarEFT,DetailedOutput->True],$Aborted];
   If[!TrueQ[coreCheck] || !AssociationQ[detailedCore] ||
      !And@@(TrueQ[Lookup[detailedCore,#,False]]& /@
         {"Hermiticity","CanonicallyNormalized","HeavyMassBasis",
          "GaugeInvariant","AllObjectsDefined"}),
     T3MEFail["IndependentCanonicalReassemblyCheckNotStable"]];
   hermiticityRewrite="IndependentCanonicalSourceNotProvenEquivalentToSavedEOMEFT";
   T3MECheckpoint["independent-canonical-full-EFT-candidate-verified"];
 ];
 (* Always export compact diagnostics so a true Hermiticity failure can be
    distinguished from a registration problem without a 100-kB box dump. *)
 hermiticityDiagnosticPath=FileNameJoin[{DirectoryName[outPath],
   FileBaseName[outPath]<>"_hermiticity_diagnostics.json"}];
 treeHermiticityReport=<|
   "model"->modelName,
   "matching_scope"->If[TrueQ[$T3MELeadingC5Only],
     "LeadingC5OneLoopSectorOnly","FullEFT"],
   "full_saved_eft_hermiticity"->fullSavedTreeHermiticity,
   "excluded_four_heavy_scalar_quartic_terms"->Lookup[leadingC5Partition,
     "excluded_count",0],
   "source"->FileNameTake[treePath],
   "light_EFT1_hermiticity"->If[AssociationQ[lightTreeDetail],
     Lookup[lightTreeDetail,"Hermiticity","Unavailable"],"Unavailable"],
   "heavy_EFT1_original_hermiticity"->If[AssociationQ[originalTreeDetail],
     Lookup[originalTreeDetail,"Hermiticity","Unavailable"],"Unavailable"],
   "heavy_EFT1_original_canonical"->If[AssociationQ[originalTreeDetail],
     Lookup[originalTreeDetail,"CanonicallyNormalized","Unavailable"],"Unavailable"],
   "heavy_EFT1_original_mass_basis"->If[AssociationQ[originalTreeDetail],
     Lookup[originalTreeDetail,"HeavyMassBasis","Unavailable"],"Unavailable"],
   "equivalence_rewrite_accepted"->hermiticityRewrite,
   "reconstruction_audit_file"->If[FileExistsQ[reconstructionPath],
      FileNameTake[reconstructionPath],""],
   "equivalence_rewrite_attempts"->hermiticityAttempts,
   "final_ordinary_CheckLagrangian"->ToString[coreCheck,InputForm],
   "final_hermiticity"->If[AssociationQ[detailedCore],
     Lookup[detailedCore,"Hermiticity","Unavailable"],"Unavailable"],
   "physical_sequential_C5_computed"->False|>;
 Quiet@Check[Export[hermiticityDiagnosticPath,treeHermiticityReport,"RawJSON"],$Failed];
 Print["Scalar-threshold Hermiticity audit: ",
   ToString[KeyTake[treeHermiticityReport,{
     "light_EFT1_hermiticity","heavy_EFT1_original_hermiticity",
     "equivalence_rewrite_accepted","final_hermiticity"}],InputForm]];
 If[!TrueQ[coreCheck] || !AssociationQ[detailedCore] ||
    !TrueQ[Lookup[detailedCore,"Hermiticity",False]] ||
    !TrueQ[Lookup[detailedCore,"HeavyMassBasis",False]],
   hermiticitySectorPath=FileNameJoin[{DirectoryName[outPath],
     FileBaseName[outPath]<>"_hermiticity_sector_audit.json"}];
   hermiticitySectorResult=T3MEHermiticitySectorAudit[
     scalarEFT,hermiticitySectorPath];
   Print["Scalar-threshold Hermiticity sector diagnostic: ",hermiticitySectorPath];
   Print["Scalar-threshold diagnostic file: ",hermiticityDiagnosticPath];
   T3MEFail["SavedEFT1TreeHeavyScalarCheckFailed",
     "ordinary="<>ToString[coreCheck,InputForm]<>
     "; Hermiticity="<>ToString[If[AssociationQ[detailedCore],
       Lookup[detailedCore,"Hermiticity",False],"Unavailable"],InputForm]<>
     "; diagnostic="<>FileNameTake[hermiticityDiagnosticPath]]];
 (* In check-only mode do not call Match. A passed independent check
    still cannot prove the EFT1 EOMEFT serialization is equivalent. *)
 If[TrueQ[$T3MEIndependentCheckOnly],
   If[hermiticityRewrite =!=
      "IndependentCanonicalSourceNotProvenEquivalentToSavedEOMEFT",
     T3MEFail["IndependentCanonicalCandidateNotEstablished"]];
   T3MECheckpoint["independent-canonical-check-only-passed"];
   Return[<|
     "format"->"t3-bc-y-only-independent-canonical-probe-v1",
     "status"->"IndependentCanonicalFullEFTCheckPassedEOMEquivalenceUnproved",
     "model"->modelName,
     "saved_eomeft_algebraic_equality_proven"->False,
     "saved_eomeft_hermiticity_passed"->False,
     "independent_canonical_lagrangian_check_passed"->True,
     "matched_majorana_term_count"->Lookup[reconstruction,
       "matched_majorana_term_count",Missing["Missing"]],
     "quartic_self_sectors_preserved"->True,
     "heavy_scalar_fields"->{"NewScalarY"},
     "light_scalar_fields"->{"H","NewScalarX"},
     "hermiticity_diagnostic_file"->FileNameTake[hermiticityDiagnosticPath],
     "reconstruction_audit_file"->FileNameTake[reconstructionPath],
     "matching_performed"->False,
     "physical_matching_performed"->False,
     "physical_sequential_C5_computed"->False|>]];
 If[hermiticityRewrite =!=
    "IndependentCanonicalSourceNotProvenEquivalentToSavedEOMEFT",
   T3MEFail["IndependentCanonicalCandidateNotEstablished"]];
 T3MECheckpoint["scalar-threshold-match-loader-started"];
 matchingSource=FileNameJoin[{scriptDir,"RunMatching.wl"}];
 If[!FileExistsQ[matchingSource],T3MEFail["MatchingFunctionsMissing"]];
 Get[matchingSource];
 If[!TrueQ[Length[DownValues[RunT3Matching]]>0],
   T3MEFail["MatchingFunctionsUnavailable"]];
 (* Store the exact Matchete stage evidence even if the subsequent
    loop-correction guard fails.  In particular, RunT3Matching's Success
    only certifies that its stages returned rather than $Failed; warnings
    about heavy fields being light require independent inspection. *)
 matchAuditPath=FileNameJoin[{DirectoryName[outPath],
   FileBaseName[outPath]<>"_matching_stage_audit.json"}];

 (* A new-kernel continuation prevents the second threshold from inheriting
    one-loop Matchete internal state from prior Y matching.  The Y-EOMEFT
    artifacts are immutable and separately source-hash-verified by the Python
    launcher. This is an experimental raw-operator check only; the saved
    F-first EOMEFT is NOT proved equivalent to the independent reconstruction. *)
 If[Environment["T3_BC_X_FRESH_EOM_TREE"] =!= "",
   Return[Module[{pathTree,pathFull,txtTree,txtFull,parsedTree,parsedFull,
           unsupported,xfresh,decision},
     pathTree=Environment["T3_BC_X_FRESH_EOM_TREE"];
     pathFull=Environment["T3_BC_X_FRESH_EOM_FULL"];
     T3MECheckpoint["x-fresh-load-saved-y-eomeft-started"];
     If[!StringQ[pathTree] || !StringQ[pathFull] ||
        !FileExistsQ[pathTree] || !FileExistsQ[pathFull],
       T3MEFail["XFreshSavedYEOMArtifactsMissing"]];
     txtTree=Quiet@Check[Import[pathTree,"Text"],$Failed];
     txtFull=Quiet@Check[Import[pathFull,"Text"],$Failed];
     unsupported={";","Get[","Import[","Export[","RunProcess[",
       "Run[","ToExpression[","Set[","SetDelayed[","Quit[","Exit[",
       "Remove[","SystemOpen[","URLRead["};
     If[!StringQ[txtTree] || !StringQ[txtFull] ||
        StringLength[txtTree]<500 || StringLength[txtFull]<500 ||
        AnyTrue[unsupported,StringContainsQ[txtTree,#]||StringContainsQ[txtFull,#]&] ||
        !StringContainsQ[txtTree,"Field[NewScalarX"] ||
        StringContainsQ[txtTree,"Field[NewScalarY"] ||
        StringContainsQ[txtTree,"Field[NewFermion"],
       T3MEFail["XFreshSavedYEOMInputInvalid"]];
     parsedTree=CheckAbort[Quiet@Check[ToExpression[txtTree,InputForm],$Failed],$Aborted];
     parsedFull=CheckAbort[Quiet@Check[ToExpression[txtFull,InputForm],$Failed],$Aborted];
     If[MemberQ[{$Failed,$Aborted},parsedTree] ||
        MemberQ[{$Failed,$Aborted},parsedFull],
       T3MEFail["XFreshSavedYEOMParseFailed"]];
     T3MECheckpoint["x-fresh-saved-y-eomeft-parsed"];
     xfresh=T3BCXFollowup[<|"EOMEFT"->parsedTree|>,
                        <|"EOMEFT"->parsedFull|>,outPath];
     decision=Lookup[xfresh,"status","XFollowupReturnedNoStatus"];
     T3MECheckpoint["x-fresh-X-followup-returned"];
     Return[<|
        "format"->"t3-bc-x-fresh-kernel-eomeft-source-v3",
        "status"->"SavedYTreeEOMConsumedForFreshXNoYRematch",
        "model"->modelName,
        "x_followup_status"->decision,
        "fresh_kernel_without_y_matching"->True,
        "x_source_is_saved_y_tree_eomeft"->True,
        "y_one_loop_correction_inherited_in_x_source"->False,
        "saved_fermion_eomeft_equivalence_certified"->False,
        "running_boundary_couplings_applied"->False,
        "physical_sequential_matching_certified"->False|>]
   ]]
 ];

 T3MECheckpoint["scalar-threshold-tree-Match-started"];
 treeResult=CheckAbort[RunT3Matching[scalarEFT,5,0],$Aborted];
 T3MECheckpoint["scalar-threshold-tree-Match-returned"];
 If[!AssociationQ[treeResult] || Lookup[treeResult,"Status",""]=!="Success" ||
    !KeyExistsQ[treeResult,"MatchedEFT"],
   T3MEFail["ScalarThresholdTreeMatchingFailed",
     ToString[If[AssociationQ[treeResult],Lookup[treeResult,"Status","Unknown"],
       treeResult],InputForm]]];
 treeMatched=treeResult["MatchedEFT"];
 scalarTreePath=FileNameJoin[{DirectoryName[outPath],
   FileBaseName[outPath]<>"_sm_tree.inputform.txt"}];
 If[Quiet@Check[Export[scalarTreePath,ToString[InputForm[treeMatched]],"Text"],
    $Failed]===$Failed,T3MEFail["ScalarThresholdTreeExportFailed"]];
 T3MECheckpoint["scalar-threshold-one-loop-Match-started"];
 loopResult=CheckAbort[RunT3Matching[scalarEFT,5,1],$Aborted];
 T3MECheckpoint["scalar-threshold-one-loop-Match-returned"];
 If[!AssociationQ[loopResult] || Lookup[loopResult,"Status",""]=!="Success" ||
    !KeyExistsQ[loopResult,"MatchedEFT"],
   T3MEFail["ScalarThresholdOneLoopMatchingFailed",
     ToString[If[AssociationQ[loopResult],Lookup[loopResult,"Status","Unknown"],
       loopResult],InputForm]]];
 loopMatched=loopResult["MatchedEFT"];
 (* Export stage diagnostics *before* checking that the resulting EFT
    actually contains the scalar-loop correction. *)
 stageSummary=Association@Table[k->Module[{v,textRep},
   v=Lookup[treeResult,k,Missing["StageAbsent"]];
   textRep=ToString[InputForm[v]];
   <|"head"->ToString[Head[v],InputForm],
     "characters"->StringLength[textRep],
     "explicit_hbar"->StringContainsQ[textRep,"hbar"],
     "contains_X"->StringContainsQ[textRep,"Field[NewScalarX"],
     "contains_Y"->StringContainsQ[textRep,"Field[NewScalarY"],
     "contains_GammaCC"->StringContainsQ[textRep,"GammaCC"],
     "prefix"->StringTake[textRep,UpTo[1100]]|>],
   {k,{"RawEFT","EOMEFT","MatchedEFT"}}];
 stageDiagnostic=Association@Table[k->Module[{v,textRep},
   v=Lookup[loopResult,k,Missing["StageAbsent"]];
   textRep=ToString[InputForm[v]];
   <|"head"->ToString[Head[v],InputForm],
     "characters"->StringLength[textRep],
     "explicit_hbar"->StringContainsQ[textRep,"hbar"],
     "contains_X"->StringContainsQ[textRep,"Field[NewScalarX"],
     "contains_Y"->StringContainsQ[textRep,"Field[NewScalarY"],
     "contains_GammaCC"->StringContainsQ[textRep,"GammaCC"],
     "prefix"->StringTake[textRep,UpTo[1100]]|>],
   {k,{"RawEFT","EOMEFT","MatchedEFT"}}];
 Export[matchAuditPath,<|
   "status"->"MatchStageAuditBeforeLoopValidation",
   "model"->modelName,"matching_scope"->If[TrueQ[$T3MELeadingC5Only],
     "LeadingC5OneLoopSectorOnly","FullEFT"],
   "tree_stage"->stageSummary,"one_loop_stage"->stageDiagnostic,
   "full_eft_validated"->!TrueQ[$T3MELeadingC5Only],
   "physical_sequential_C5_computed"->False|>,"RawJSON"];
 correction=Expand[loopMatched-treeMatched];
 correctionText=ToString[InputForm[correction]];
 completeText=ToString[InputForm[loopMatched]];
 If[correction === 0 || !StringContainsQ[correctionText,"hbar"] ||
    !StringContainsQ[completeText,"Field[NewScalarX"] ||
    StringContainsQ[completeText,"Field[NewScalarY"] ||
    StringContainsQ[completeText,"Field[NewFermion"],
   (Print["Scalar-threshold diagnostic: tree/loop stage audit -> ",matchAuditPath];
    Print["  loop correction zero: ",correction===0,
      "; has hbar: ",StringContainsQ[correctionText,"hbar"],
      "; matched EFT retains X/Y/F: ",
      Or@@(StringContainsQ[completeText,#]& /@
        {"Field[NewScalarY","Field[NewFermion"})];
    T3MEFail["ScalarThresholdLoopCorrectionInvalid",
       "See "<>FileNameTake[matchAuditPath]])];
 scalarLoopPath=FileNameJoin[{DirectoryName[outPath],
   FileBaseName[outPath]<>"_sm_tree_plus_loop.inputform.txt"}];
 scalarCorrectionPath=FileNameJoin[{DirectoryName[outPath],
   FileBaseName[outPath]<>"_scalar_one_loop_correction.inputform.txt"}];
 If[Quiet@Check[Export[scalarLoopPath,completeText,"Text"],$Failed]===$Failed ||
    Quiet@Check[Export[scalarCorrectionPath,correctionText,"Text"],$Failed]===$Failed,
   T3MEFail["ScalarThresholdLoopExportFailed"]];
 T3MECheckpoint["scalar-threshold-LLHH-extract-started"];
 deltaTerms=With[{expanded=Expand[correction]},
   If[Head[expanded]===Plus,List@@expanded,{expanded}]];
 scalarCandidates=Select[deltaTerms,Function[term,
   Module[{s=ToString[InputForm[term]]},
     StringContainsQ[s,"GammaCC"] &&
     StringCount[s,"Field[H, Scalar"] === 2 &&
     StringCount[s,"Field[l, Fermion"] === 2 &&
     !StringContainsQ[s,"Field[NewScalarX"] &&
     !StringContainsQ[s,"Field[NewScalarY"] &&
     !StringContainsQ[s,"Field[NewFermion"]]]];
 llhh=If[scalarCandidates==={},0,Total[scalarCandidates]];
 llhhText=ToString[InputForm[llhh]];
 scalarLLHHPath=FileNameJoin[{DirectoryName[outPath],
   FileBaseName[outPath]<>"_raw_scalar_llhh.inputform.txt"}];
 If[Quiet@Check[Export[scalarLLHHPath,llhhText,"Text"],$Failed]===$Failed,
   T3MEFail["ScalarThresholdLLHHExportFailed"]];
 T3MECheckpoint["scalar-threshold-LLHH-extract-returned"];
 Print["T3 scalar-threshold tree+loop matching: total loop terms=",
   Length[deltaTerms],"; raw LLHH candidates=",Length[scalarCandidates]];
  (* Additional EFT2 operator inventory after native Y-only matching.
     Only terms containing two light leptons and exactly two H/X scalars
     are exported; the complete matched EFT2 is separately preserved above.
     These are raw one-loop operators and are not renormalized Wilson seeds. *)
  $T3BCYOnlySectorCounts=<||>;
  Module[{baseTerms, labels, eligible, part, nH, nX, filePath, label},
    baseTerms=Select[deltaTerms,Function[term,
      With[{s=ToString[InputForm[term]]},
        StringContainsQ[s,"GammaCC"] &&
        StringCount[s,"Field[l, Fermion"]===2 &&
        !StringContainsQ[s,"Field[NewScalarY"] &&
        !StringContainsQ[s,"Field[NewFermion"]]]];
    labels={{"LLHH",2,0},{"LLHX",1,1},{"LLXX",0,2}};
    Do[
      label=entry[[1]];nH=entry[[2]];nX=entry[[3]];
      eligible=Select[baseTerms,Function[term,
        With[{s=ToString[InputForm[term]]},
          StringCount[s,"Field[H, Scalar"]===nH &&
          StringCount[s,"Field[NewScalarX, Scalar"]===nX]]];
      part=If[eligible==={},0,Total[eligible]];
      filePath=FileNameJoin[{DirectoryName[outPath],
        FileBaseName[outPath]<>"_raw_"<>ToLowerCase[label]<>".inputform.txt"}];
      If[Quiet@Check[Export[filePath,ToString[InputForm[part]],"Text"],$Failed]===$Failed,
        T3MEFail["YOnlyEFT2OperatorSectorExportFailed",label]];
      $T3BCYOnlySectorCounts[label]=Length[eligible],
      {entry,labels}];
    Print["Y-only EFT2 raw Wilson sector counts: ",$T3BCYOnlySectorCounts]
  ];

 (* Experimental third threshold: only Y-tree EOMEFT is used, so any X
    one-loop result omits inherited Y-loop hard matching by construction. *)
 $T3BCXFollowupReport=T3BCXFollowup[treeResult,loopResult,outPath];
 Print["Native X followup status: ",Lookup[$T3BCXFollowupReport,"status","Unknown"]];

 <|
   "format"->"t3-bc-y-only-matchete-eigenfield-probe-v1",
   "status"->"YOnlyIndependentCanonicalRawEFT2NotCertified",
   "model"->modelName,
   "saved_eomeft_algebraic_equality_proven"->False,
   "independent_canonical_lagrangian_check_passed"->True,
   "independent_canonical_quartics_preserved"->True,
   "matching_scope"->If[TrueQ[$T3MELeadingC5Only],
     "LeadingC5OneLoopSectorOnly","FullEFT"],
   "full_saved_eft_hermiticity"->fullSavedTreeHermiticity,
   "excluded_four_heavy_scalar_quartic_terms"->Lookup[leadingC5Partition,
      "excluded_count",0],
   "full_eft_validated"->False,
   "native_canonical_candidate_checked_only"->True,
   "full_saved_eft_unchanged_on_disk"->True,
   "input_eft1_tree"->FileNameTake[treePath],
   "input_eft1_tree_character_count"->firstTreeLength,
   "eigenfield_input_yukawas_and_quartics_reconstructed"->True,
   "fermion_removed_before_scalar_matching"->True,
   "heavy_scalar_fields"->{"NewScalarY"},
    "light_scalar_fields"->{"H","NewScalarX"},
   "scalar_heavy_mass_basis_detailed_check_passed"->True,
   "hermiticity_check_passed"->True,
   "hermiticity_applies_to_matching_subtheory_only"->TrueQ[$T3MELeadingC5Only],
   "canonical_EOMEFT_hermiticity_rewrite"->hermiticityRewrite,
   "hermiticity_diagnostic_file"->FileNameTake[hermiticityDiagnosticPath],
   "reconstruction_audit_file"->FileNameTake[reconstructionPath],
   "tree_matching_completed"->True,"one_loop_matching_completed"->True,
   "scalar_matching_stage_audit_file"->FileNameTake[matchAuditPath],
   "scalar_free_lagrangian_rebased"->TrueQ[$T3MELeadingC5Only],
   "scalar_one_loop_correction_hbar_present"->True,
   "raw_scalar_llhh_term_count"->Length[scalarCandidates],
    "raw_EFT2_operator_sector_counts"->$T3BCYOnlySectorCounts,
    "x_followup_status"->Lookup[$T3BCXFollowupReport,"status","NotRun"],
   "x_followup_is_physical"->False,
   "Y_only_eigenfield_matching_experimental"->True,
    "renormalized_sequential_matching_certified"->False,
   "scalar_loop_difference_term_count"->Length[deltaTerms],
   "sm_tree_file"->FileNameTake[scalarTreePath],
   "sm_one_loop_file"->FileNameTake[scalarLoopPath],
   "scalar_correction_file"->FileNameTake[scalarCorrectionPath],
   "raw_scalar_llhh_file"->FileNameTake[scalarLLHHPath],
   "heavy_Y_absent_and_light_X_retained_in_matched_eft"->True,
   "eft1_running_performed"->False,
   "hard_pole_counterterms_combined"->False,
   "hard_pole_rge_certified"->False,
   "physical_sequential_C5_computed"->False,
   "physical_matching_performed"->False
 |>
];

run=CheckAbort[
  UsingFrontEnd[Catch[T3MEScalarThresholdPerform[],"T3MassEigenStop"]],$Aborted];
If[run === $Aborted || run === $Failed || !AssociationQ[run],
  If[!FileExistsQ[outPath],
    T3MEOutput[<|"status"->"UncaughtScalarThresholdError",
      "model"->modelName,"last_stage"->$T3MELastStage,
      "stage_history"->$T3MEStageHistory,
      "kernel_result_inputform"->ToString[run,InputForm],
      "message_list_inputform"->ToString[$MessageList,InputForm],
      "physical_matching_performed"->False|>]];
  Exit[5]];
T3MEOutput[run];
If[!MemberQ[{
     "IndependentCanonicalFullEFTCheckPassedEOMEquivalenceUnproved",
     "YOnlyIndependentCanonicalRawEFT2NotCertified",
     "SavedYTreeEOMConsumedForFreshXNoYRematch",
     "EFT1FullResidualExportedBeforeMatching"},
     run["status"]],Exit[5]];
