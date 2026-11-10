(* Consolidated full-general T3 threshold continuation, local-only.
   One common F threshold followed by one common scalar threshold.
   A(-4) reuses the three-flavour F boundary, reparameterises LLSS
   flavour contractions and promotes both scalar masses in a dedicated
   heavy-scalar Matchete continuation.
   D(-2) dispatches to the physical split-real-scalar continuation.
   This source does NOT claim finite renormalised C5 or EFT1 Wilson running.
   CLI: wolframscript -file SCRIPT OUT_DIR A_-4|D_-2 [A_REUSE_F_DIR]
*)
ClearAll["Global`*"];
t3args = Rest[$ScriptCommandLine];
If[!MemberQ[{2, 3}, Length[t3args]] ||
   !MemberQ[{"A_-4", "D_-2"}, t3args[[2]]],
  Print["Usage: OUT_DIR A_-4|D_-2 [A_REUSE_F_DIR]"]; Exit[2]];
t3out = ExpandFileName[t3args[[1]]];
t3modelName = t3args[[2]];
t3reuse = If[Length[t3args] === 3, ExpandFileName[t3args[[3]]], None];
t3root = DirectoryName @ ExpandFileName[$InputFileName];
t3project = DirectoryName[t3root];
If[!DirectoryQ[t3out], CreateDirectory[t3out, CreateIntermediateDirectories -> True]];
t3status = FileNameJoin[{t3out, "native_two_threshold_status.json"}];
t3stage = "initial";
ClearAll[t3Write, t3Stop, t3Save];
t3Write[status_String, reason_String : "", more_ : <||>] := Module[{rec},
  rec = Join[<|
    "format" -> "t3-full-general-two-threshold-local-continuation-v3",
    "model" -> t3modelName, "status" -> status, "stage" -> t3stage,
    "reason" -> reason, "heavy_fermion_flavours" -> 3,
    "fermion_groups" -> {{"F"}}, "scalar_groups" -> {{"S1", "S2"}},
    "group_count" -> 2, "one_loop_matching_requested" -> True,
    "eft1_wilson_rge_inserted" -> False,
    "renormalized_finite_wilsons_verified" -> False,
    "complete_hard_plus_soft_matching_performed" -> False,
    "full_general_numerical_c5_verified" -> False|>, more];
  Export[t3status, rec, "RawJSON"];
  Print[t3modelName, " ", status, " at ", t3stage,
    If[reason === "", "", ": " <> reason]];
  rec
];
t3Stop[msg_String, more_ : <||>] := (
  t3Write["Blocked_NativeStage", msg, more]; Exit[5]);
t3Save[name_String, obj_, kind_String : "WXF"] := Module[
  {p = FileNameJoin[{t3out, name}]},
  If[Quiet @ Check[Export[p, obj, kind], $Failed] === $Failed ||
     !FileExistsQ[p], t3Stop["Could not export " <> name]];
  p
];
If[!And @@ (Environment[#] === "1" & /@ {
    "T3_EXPERIMENTAL_GENERAL_FERMIONS",
    "T3_EXPERIMENTAL_GENERAL_QUARTICS",
    "T3_EXPERIMENTAL_GENERAL_QUADRATICS",
    "T3_EXPERIMENTAL_FLAVOR_MASSES"}),
  t3Stop["All four native full-general switches must be set to 1"]];

(* This is a separate real-mass-basis kernel continuation, NOT the generic
   NewScalar2 kernel. It consumes the saved F output, never re-matches F. *)
If[t3modelName === "D_-2",
  t3DFile = FileNameJoin[{t3root, "RunT3DMinus2FullGeneralPhysicalScalar.wl"}];
  If[!FileExistsQ[t3DFile], t3Stop["Missing D(-2) physical continuation " <> t3DFile]];
  Get[t3DFile];
  t3Stop["D(-2) physical continuation returned without a status"]
];

(* The A(-4) F stage contains all three F flavours.  If previously saved
   canonical stage-2 inputs exist, use them exactly. *)
t3stage = "A-matchete-load";
If[CheckAbort[UsingFrontEnd[Needs["Matchete`"]; True], $Aborted] =!= True,
  t3Stop["Matchete could not be loaded"]];
Get[FileNameJoin[{t3root, "model", "T3ModelCatalog.wl"}]];
Get[FileNameJoin[{t3root, "model", "LagrangianBuilder.wl"}]];
Get[FileNameJoin[{t3root, "RunMatching.wl"}]];
t3model = T3ModelFromClass["A", -4, 3];
If[!AssociationQ[t3model], t3Stop["Three-flavour A(-4) catalog failed"]];
If[SetT3HeavyFields[{"F"}] === $Failed,
  t3Stop["Common F threshold declaration failed"]];
t3stage = "A-F-boundary-reuse-or-match";
If[StringQ[t3reuse],
  t3treePath = FileNameJoin[{t3reuse, "input_tree.wxf"}];
  t3loopPath = FileNameJoin[{t3reuse, "input_loop.wxf"}];
  If[!And @@ (FileExistsQ /@ {t3treePath, t3loopPath}),
    t3Stop["Requested F reuse lacks tree/loop WXF"]];
  t3tree = Quiet @ Check[Import[t3treePath, "WXF"], $Failed];
  t3hard = Quiet @ Check[Import[t3loopPath, "WXF"], $Failed];
  If[MemberQ[{t3tree, t3hard}, $Failed] || FreeQ[t3tree, GammaCC] ||
     FreeQ[t3hard, hbar],
    t3Stop["Reused F WXF lacks LLSS or one-loop hard source"]];
  t3reused = True,
  t3build = CheckAbort[UsingFrontEnd[BuildT3Lagrangian[t3model]], $Aborted];
  If[!AssociationQ[t3build] || Lookup[t3build, "Status", ""] =!= "Success" ||
     !TrueQ[Lookup[t3build, "FullValidation", False]] ||
     !TrueQ[Lookup[t3build, "AllScalarExtensionsRequested", False]] ||
     Lookup[t3build, "AdditionalQuarticsIncluded", -1] =!= 0,
    t3Stop["Full general A(-4) UV builder failed"]];
  t3f = CheckAbort[UsingFrontEnd[RunSequentialT3Matching[
      t3build["LUV"], t3model, {{"F"}}, 5, 1,
      t3out, t3modelName]], $Aborted];
  If[!AssociationQ[t3f] || Lookup[t3f, "Status", ""] =!= "Success" ||
     Length[Lookup[t3f, "Stages", {}]] =!= 1,
    t3Stop["F-only native matching failed"]];
  t3tree = First[t3f["Stages"]]["TreeLagrangian"];
  t3hard = First[t3f["Stages"]]["OneLoopCorrection"];
  t3reused = False
];
t3Save["A_full_F_tree_original.wxf", t3tree];
t3Save["A_full_F_hard_one_loop_preserved.wxf", t3hard];

(* IMPORTANT FIX: DO NOT invoke DefineCoupling here. This parent kernel
   has no reconstructed EFT1 Matchete field/coupling registry. The fresh
   scalar kernel RunThresholdStageEFT1Wilson.wl performs ResetAll[],
   LoadModel["SM"], DefineT3FieldObjects[], DefineT3Couplings[], then
   registers T3EFT1W12/T3EFT1WExtra. A symbolic Coupling[...] is valid
   for WXF transport without parent-kernel registration. *)
t3stage = "A-tree-flavour-Wilson-contract";
ClearAll[t3InverseIndex, t3YukawaFactor, t3ContractTerm];
t3InverseIndex[x_] := Replace[Unevaluated[x], {
   HoldPattern[Power[Coupling[MF, {r:Index[_, NFlavor]}, 0], -1]] :> r,
   _ :> Missing["NotIndexedInverseMass"]}];
t3YukawaFactor[x_] := Replace[Unevaluated[x], {
   HoldPattern[Coupling[n:(y1|y2),
      {p:Index[_, Flavor], r:Index[_, NFlavor]}, 0]] :>
     <|"name" -> n, "p" -> p, "r" -> r, "barred" -> False|>,
   HoldPattern[Bar[Coupling[n:(y1|y2),
      {p:Index[_, Flavor], r:Index[_, NFlavor]}, 0]]] :>
     <|"name" -> n, "p" -> p, "r" -> r, "barred" -> True|>,
   _ :> Missing["NotYukawa"]}];
t3WCount = 0;
t3ContractTerm[term_] := Module[
  {factors, inds, massPos, heavy, ys, y1pos, y2pos,
   first, second, remaining, w},
  factors = If[Head[term] === Times, List @@ term, {term}];
  inds = t3InverseIndex /@ factors;
  massPos = Select[Range[Length[inds]], !MissingQ[inds[[#]]] &];
  If[massPos === {}, Return[term]];
  If[Length[massPos] =!= 1, t3Stop["Multiple inverse MF in a tree monomial"]];
  heavy = inds[[First[massPos]]];
  ys = t3YukawaFactor /@ factors;
  y1pos = Select[Range[Length[ys]], AssociationQ[ys[[#]]] &&
    Lookup[ys[[#]], "name", None] === y1 &&
    Lookup[ys[[#]], "r", None] === heavy &];
  y2pos = Select[Range[Length[ys]], AssociationQ[ys[[#]]] &&
    Lookup[ys[[#]], "name", None] === y2 &&
    Lookup[ys[[#]], "r", None] === heavy &];
  If[Length[y1pos] =!= 1 || Length[y2pos] =!= 1,
    t3Stop["No unique y1/y2/MF heavy-flavour contraction"]];
  first = ys[[First[y1pos]]]; second = ys[[First[y2pos]]];
  If[first["barred"] =!= second["barred"],
    t3Stop["Mixed conjugation in F boundary"]];
  w = Coupling[T3EFT1W12, {first["p"], second["p"]}, 0];
  If[TrueQ[first["barred"]], w = Bar[w]];
  remaining = Delete[factors, List /@ Sort[{
     First[massPos], First[y1pos], First[y2pos]}]];
  t3WCount++;
  Times @@ Append[remaining, w]
];
t3terms = If[Head[t3tree] === Plus, List @@ t3tree, {t3tree}];
t3treeW = Plus @@ (t3ContractTerm /@ t3terms);
If[t3WCount =!= 2 || FreeQ[t3treeW, T3EFT1W12] ||
   !FreeQ[t3treeW, Power[Coupling[MF, {Index[_, NFlavor]}, 0], -1]],
  t3Stop["Unexpected A(-4) Wilson tree inventory",
    <|"converted_W_terms" -> t3WCount|>]];
t3Save["A_EFT1_tree_contracted.wxf", t3treeW];
t3Save["A_EFT1_tree_contracted.txt", ToString[t3treeW, InputForm], "Text"];
t3Save["A_EFT1_Wilson_definitions.json", <|
   "W12[p,q]" -> "Sum[y1[p,r]*y2[q,r]/MF[r],{r,1,3}]",
   "conjugate_W12" -> "Conjugate[W12[p,q]] for positive real MF[r]",
   "all_three_fermions_integrated_together" -> True,
   "one_loop_hard_reparameterised" -> False,
   "renormalised_Wilson_boundary_verified" -> False|>, "RawJSON"];

(* The dedicated physical A scalar-threshold continuation is analogous to
   the verified D(-2) direct-Match strategy.  It rebuilds the Matchete
   registry and promotes BOTH scalar mass power-counting labels from 1 to 0
   before calling Match with both scalars heavy in one threshold.
   The original F hard source stays untouched and separate. *)
t3AFile = FileNameJoin[{t3root, "RunT3AMinus4FullGeneralPhysicalScalar.wl"}];
If[!FileExistsQ[t3AFile],
  t3Stop["Missing A(-4) physical scalar continuation " <> t3AFile]];
Get[t3AFile];
t3Stop["A(-4) physical scalar continuation returned without a status"];
