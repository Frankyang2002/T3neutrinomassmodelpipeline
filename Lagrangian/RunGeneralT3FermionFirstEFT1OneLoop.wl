(* Experimental one-loop F-first matching at the initial EFT1 boundary (B_-1/C_-1).
   Fork of the independently live-validated mass-eigenfield Lagrangian builder.
   At the F threshold only F is HEAVY; X and Y are LIGHT and remain explicit.
   Match at EFTOrder -> 5, LoopOrder -> 0 and then -> 1 on the same
   UV tree Lagrangian. Subtract Matchete-matched EFTs to isolate the
   one-loop correction, and export a RAW hard LLHH candidate sector.
   No pole subtraction, scheme-normalized C5, EFT1 RGEs, or scalar
   threshold is claimed. This is not a final neutrino-mass prediction.

   wolframscript -file Lagrangian/RunGeneralT3FermionFirstEFT1OneLoop.wl \
     output/general_t3_mass_eigen_matchete/B_m1_final_parameters.json \
     output/general_t3_mass_eigen_matchete/eft1_fermion_first_oneloop/B_m1_final_eft1_fermion_first_oneloop.json

   The numerical mass mixing is diagonalised upstream in Python. Matchete
   is rebuilt with F HEAVY and TWO LIGHT Y=-1/2 eigen-doublets X,Y.  The
   scalar kinetic/mass terms are diagonal, and all native SU(2) interaction
   tensors are the independently checked mass-basis tensors. No scalar is
   integrated out. This is a one-loop FIRST-THRESHOLD matching probe.
*)
ClearAll["Global`*"];
args = Rest[$ScriptCommandLine];
If[Length[args] =!= 2,
  Print["Usage: RunGeneralT3FermionFirstEFT1OneLoop.wl INPUT.json OUTPUT.json"];
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

T3MEPerform[] := Module[
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
  hasLoopHbar, fLoopAbsent, xLoopPresent, yLoopPresent, hardFAbsent},
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
 T3MECheckpoint["rotate-original-interactions-started"];
 otherRotated=T3MERotateInteraction[extraParts["RemainingInteractions"]];
 If[otherRotated === $Failed,
   T3MEFail["OriginalScalarFieldsRemainAfterRotation"]];
 portalsRotated=Total[Values[nativeH2["Expressions"]]];
 scalarRotated=Total[Values[nativeS2["Expressions"]]];
 extraRotated=Total[Values[nativeExtra["Expressions"]]];
 rotated=Expand[otherRotated+portalsRotated+scalarRotated+extraRotated+
    nativeTopology["Expression"]+nativeY2["Expression"]+
    nativeY1["Expression"]]
   // RelabelIndices;
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

 T3MECheckpoint["fermion-first-ordinary-CheckLagrangian-started"];
 check=CheckAbort[CheckLagrangian[lag],$Aborted];
 T3MECheckpoint["fermion-first-detailed-CheckLagrangian-started"];
 detailed=CheckAbort[CheckLagrangian[lag,DetailedOutput->True],$Aborted];
 heavy=AssociationQ[detailed] && TrueQ[Lookup[detailed,"HeavyMassBasis",False]];
 T3MECheckpoint["fermion-first-lagrangian-checks-returned"];
 If[!TrueQ[check] || !heavy,
   Print["F-first light-scalar CheckLagrangian: ",InputForm[check]];
   Print["F-first detailed result: ",ToString[detailed,InputForm]];
   T3MEFail["FermionFirstLightScalarLagrangianCheckFailed"]];

 (* The first EFT matching threshold must integrate out ONLY the fermion.
    NewScalarX and NewScalarY are registered Light above, never Heavy.
    Do not integrate out the scalars or insert the previously obtained
    all-heavy one-loop C5 into this intermediate EFT. *)
 T3MECheckpoint["fermion-first-production-Match-loader-started"];
 matchingSource=FileNameJoin[{scriptDir,"RunMatching.wl"}];
 If[!FileExistsQ[matchingSource],
   T3MEFail["MatchingFunctionsMissing",matchingSource]];
 Get[matchingSource];
 If[!TrueQ[Length[DownValues[RunT3Matching]]>0],
   T3MEFail["RunT3MatchingUnavailable"]];
 T3MECheckpoint["fermion-first-tree-Match-started"];
 matchingResult=CheckAbort[RunT3Matching[lag,5,0],$Aborted];
 T3MECheckpoint["fermion-first-tree-Match-returned"];
 If[!AssociationQ[matchingResult] ||
    Lookup[matchingResult,"Status",""] =!= "Success" ||
    !KeyExistsQ[matchingResult,"MatchedEFT"],
   Print["F-first matching status: ",
     ToString[If[AssociationQ[matchingResult],
       KeyDrop[matchingResult,{"RawEFT","GreenEFT","EOMEFT","LoopEFT","MatchedEFT"}],
       matchingResult],InputForm]];
   T3MEFail["FermionFirstTreeMatchingFailed",
     ToString[If[AssociationQ[matchingResult],
       Lookup[matchingResult,"Status","Unknown"],matchingResult],InputForm]]];
 matchedEFT=matchingResult["MatchedEFT"];
 eftPath=FileNameJoin[{DirectoryName[outPath],
   FileBaseName[outPath]<>"_eft1_tree.inputform.txt"}];
 If[Quiet@Check[Export[eftPath,ToString[InputForm[matchedEFT]],"Text"],$Failed]===$Failed,
   T3MEFail["EFT1TreeExportFailed"]];
 (* Textual diagnostics intentionally follow the existing
    RunMatching.wl EFT1 Wilson seed extractor: Matchete packages internal
    GammaCC/Field symbols in contexts not safely matched by simple global
    Wolfram patterns. Select the dimension-five LLSS sector conservatively. *)
 terms=With[{expanded=Expand[matchedEFT]},
   If[Head[expanded]===Plus,List@@expanded,{expanded}]];
 llssTerms=Select[terms,Function[term,
   Module[{s=ToString[InputForm[term]]},
     StringContainsQ[s,"GammaCC"] &&
     (StringContainsQ[s,"Field[NewScalarX"] ||
      StringContainsQ[s,"Field[NewScalarY"])]]];
 llss=If[llssTerms==={},0,Total[llssTerms]];
 llssPath=FileNameJoin[{DirectoryName[outPath],
   FileBaseName[outPath]<>"_llss_tree.inputform.txt"}];
 If[Quiet@Check[Export[llssPath,ToString[InputForm[llss]],"Text"],$Failed]===$Failed,
   T3MEFail["EFT1LLSSExportFailed"]];
 T3MECheckpoint["fermion-first-LLSS-sector-extracted"];
 xPresent=!FreeQ[matchedEFT,Verbatim[NewScalarX],{0,Infinity}];
 yPresent=!FreeQ[matchedEFT,Verbatim[NewScalarY],{0,Infinity}];
 fAbsent=FreeQ[matchedEFT,Verbatim[NewFermion],{0,Infinity}];
 xLLSS=AnyTrue[llssTerms,StringContainsQ[ToString[InputForm[#]],"Field[NewScalarX"]&];
 yLLSS=AnyTrue[llssTerms,StringContainsQ[ToString[InputForm[#]],"Field[NewScalarY"]&];
 If[Length[llssTerms]<1 || !xPresent || !yPresent || !fAbsent || !xLLSS || !yLLSS,
   Print["EFT1 LLSS diagnostics: ",InputForm[<|
       "llss_terms"->Length[llssTerms],"X"->xPresent,
       "Y"->yPresent,"F_absent"->fAbsent,
       "X_LLSS"->xLLSS,"Y_LLSS"->yLLSS|>]];
   Return[<|
     "status"->"FermionFirstTreeMatchedButLLSSInventoryIncomplete",
     "model"->modelName,"tree_matching_completed"->True,
     "light_scalar_fields"->{"NewScalarX","NewScalarY"},
     "matched_eft_file"->FileNameTake[eftPath],
     "llss_inputform_file"->FileNameTake[llssPath],
     "llss_term_count"->Length[llssTerms],
     "fermion_field_absent"->fAbsent,
     "x_scalar_present"->xPresent,"y_scalar_present"->yPresent,
     "physical_C5_computed"->False,"sequential_EFT1_completed"->False|>]];
 T3MECheckpoint["fermion-first-LLSS-tree-success"];

 (* Both calls act on the SAME original tree-level UV Lagrangian. Their
    subtraction removes the O(hbar^0) matched EFT without introducing
    loop-order contamination from an already loop-corrected input. This
    follows RunSequentialT3Matching's first-threshold subtraction. *)
 T3MECheckpoint["fermion-first-one-loop-Match-started"];
 oneLoopMatchingResult=CheckAbort[RunT3Matching[lag,5,1],$Aborted];
 T3MECheckpoint["fermion-first-one-loop-Match-returned"];
 If[!AssociationQ[oneLoopMatchingResult] ||
    Lookup[oneLoopMatchingResult,"Status",""] =!= "Success" ||
    !KeyExistsQ[oneLoopMatchingResult,"MatchedEFT"],
   Print["F-first one-loop matching result: ",
     ToString[If[AssociationQ[oneLoopMatchingResult],
       KeyDrop[oneLoopMatchingResult,{"RawEFT","GreenEFT","EOMEFT","LoopEFT","MatchedEFT"}],
       oneLoopMatchingResult],InputForm]];
   T3MEFail["FermionFirstOneLoopMatchingFailed",
     ToString[If[AssociationQ[oneLoopMatchingResult],
       Lookup[oneLoopMatchingResult,"Status","Unknown"],oneLoopMatchingResult],InputForm]]];
 oneLoopMatchedEFT=oneLoopMatchingResult["MatchedEFT"];
 T3MECheckpoint["fermion-first-one-loop-EFT-difference-started"];
 loopDifference=Expand[oneLoopMatchedEFT-matchedEFT];
 loopText=ToString[InputForm[loopDifference]];
 hasLoopHbar=StringContainsQ[loopText,"hbar"];
 fLoopAbsent=!StringContainsQ[loopText,"Field[NewFermion"];
 xLoopPresent=StringContainsQ[loopText,"Field[NewScalarX"];
 yLoopPresent=StringContainsQ[loopText,"Field[NewScalarY"];
 If[loopDifference === 0 || !hasLoopHbar || !fLoopAbsent,
   Print["F-first loop-difference diagnostics: ",InputForm[<|
     "nonzero"->(loopDifference=!=0),"has_hbar"->hasLoopHbar,
     "no_fermion_field"->fLoopAbsent|>]];
   T3MEFail["FermionFirstOneLoopDifferenceInvalid"]];
 oneLoopEFTPath=FileNameJoin[{DirectoryName[outPath],
   FileBaseName[outPath]<>"_eft1_tree_plus_loop.inputform.txt"}];
 oneLoopCorrectionPath=FileNameJoin[{DirectoryName[outPath],
   FileBaseName[outPath]<>"_eft1_one_loop_correction.inputform.txt"}];
 If[Quiet@Check[Export[oneLoopEFTPath,
      ToString[InputForm[oneLoopMatchedEFT]],"Text"],$Failed]===$Failed ||
    Quiet@Check[Export[oneLoopCorrectionPath,loopText,"Text"],$Failed]===$Failed,
   T3MEFail["FermionFirstOneLoopEFTExportFailed"]];

 (* Only report a raw operator SECTOR. Do not strip fields, assert a
    physical coefficient or drop divergent poles: the hard coefficient
    requires separate renormalization, flavour and MSbar checks. Since
    gamma-CC Matchete expressions may contain additional loop structures,
    insist on two Higgs legs and no external light scalar legs. *)
 loopTerms=With[{expanded=Expand[loopDifference]},
   If[Head[expanded]===Plus,List@@expanded,{expanded}]];
 hardCandidates=Select[loopTerms,Function[term,
   Module[{s=ToString[InputForm[term]]},
     StringContainsQ[s,"GammaCC"] &&
     StringCount[s,"Field[H,"]>=2 &&
     !StringContainsQ[s,"Field[NewScalarX"] &&
     !StringContainsQ[s,"Field[NewScalarY"] &&
     !StringContainsQ[s,"Field[NewFermion"]]]];
 hardRaw=If[hardCandidates==={},0,Total[hardCandidates]];
 hardText=ToString[InputForm[hardRaw]];
 hardFAbsent=!StringContainsQ[hardText,"Field[NewFermion"];
 hardPath=FileNameJoin[{DirectoryName[outPath],
   FileBaseName[outPath]<>"_raw_hard_llhh_sector.inputform.txt"}];
 If[Quiet@Check[Export[hardPath,hardText,"Text"],$Failed]===$Failed || !hardFAbsent,
   T3MEFail["FermionFirstRawHardSectorExportFailed"]];
 T3MECheckpoint["fermion-first-one-loop-sector-extracted"];
 Print["F-first ONE-LOOP EFT1: LLSS tree terms=",Length[llssTerms],
   "; loop-difference terms=",Length[loopTerms],
   "; raw hard LLHH candidates=",Length[hardCandidates],
   "; hbar=",hasLoopHbar];
 <|
   "format"->"t3-mass-eigen-fermion-first-eft1-oneloop-v1",
   "status"->"FermionFirstOneLoopRawSectorsExtractedNotRenormalized",
   "model"->modelName,
   "matching_approach"->"FermionFirstTreeAndOneLoopScalarsRemainLight",
   "matching_EFT_order"->5,"matching_loop_order"->1,
   "mass_basis_ordinary_check_passed"->True,
   "fermion_heavy_mass_basis_detailed_check_passed"->True,
   "source_quartic_count"->17,"source_yukawa_count"->2,
   "scalar_light_fields"->{"NewScalarX","NewScalarY"},
   "integrated_heavy_fields"->{"NewFermion"},
   "tree_matching_completed"->True,
   "fermion_field_absent"->fAbsent,
   "x_scalar_present"->xPresent,"y_scalar_present"->yPresent,
   "llss_contains_x"->xLLSS,"llss_contains_y"->yLLSS,
   "llss_term_count"->Length[llssTerms],
   "matched_eft_file"->FileNameTake[eftPath],
   "llss_inputform_file"->FileNameTake[llssPath],
   "one_loop_matched_eft_file"->FileNameTake[oneLoopEFTPath],
   "one_loop_difference_file"->FileNameTake[oneLoopCorrectionPath],
   "raw_hard_llhh_file"->FileNameTake[hardPath],
   "one_loop_matching_completed"->True,
   "loop_difference_hbar_present"->hasLoopHbar,
   "loop_difference_fermion_absent"->fLoopAbsent,
   "loop_difference_term_count"->Length[loopTerms],
   "raw_hard_llhh_candidate_term_count"->Length[hardCandidates],
   "raw_hard_llhh_no_fermion_field"->hardFAbsent,
   "raw_hard_llhh_scheme_renormalized"->False,
   "raw_hard_llhh_pole_rge_certified"->False,
   "numerical_U_parameters"->T3MEMixingSpurionNumericalMap[],
   "eft1_rge_performed"->False,
   "scalar_threshold_matching_performed"->False,
   "one_loop_hard_fermion_matching_performed"->False,
   "physical_C5_computed"->False,
   "sequential_EFT1_completed"->False,
   "physical_matching_performed"->False,
   "matchete_rgbeta_quartic_normalization_certified"->False
 |>

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
If[run["status"] =!= "FermionFirstOneLoopRawSectorsExtractedNotRenormalized",
  Exit[5]
];
Exit[0];
