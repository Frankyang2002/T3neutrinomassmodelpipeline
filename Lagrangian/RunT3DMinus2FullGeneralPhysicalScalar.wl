(* D(-2) physical scalar-threshold continuation, loaded by
   RunT3FullGeneralTwoThresholdNative.wl in its own Wolfram process.

   Exact source: successful NF=3, one-loop F-first Batch2 EXPORTS.
   No fermion matching is re-run. The neutral S2 self-bilinear has already
   been diagonalised into T3Scalar2Minus/T3Scalar2Plus; both eigenfields and
   NewScalar1 are made HEAVY and removed in one Matchete scalar Match[].

   W12[p,q] = Sum_r y1[p,r] y2[q,r]/MF[r]
   WExtra[p,q] = Sum_r y1[p,r] y2Conj[q,r]/MF[r]
   The tree EFT1 boundary is algebraically reparameterised, not reduced
   to a single fermion flavour. F one-loop hard terms are saved separately.

   LIMITATION: no complete F-hard Wilson finite counterterms or EFT1 Wilson
   running are implemented here. Its soft one-loop Weinberg coefficient is
   explicitly NONRENORMALISED and NOT a physical C5 prediction.
*)
t3stage = "D-verify-NF3-physical-F-source";
t3DSource = FileNameJoin[{t3project, "output", "parallel_t3",
  "full_general_batch2", "thesis_full_general_batch2_raw_v1",
  "D_-2", "native_matching"}];
t3DData = FileNameJoin[{t3DSource, "data"}];
t3DStatus = FileNameJoin[{t3DSource, "native_status.json"}];
t3DBasis = FileNameJoin[{t3DSource, "self_mass_basis.json"}];
t3DTreePath = FileNameJoin[{t3DData, "eft1_after_F_tree_inputform.txt"}];
t3DHardPath = FileNameJoin[{t3DData, "eft1_after_F_one_loop_inputform.txt"}];
t3DWilsonSeed = FileNameJoin[{t3DData, "eft1_after_F_wilson_seed.json"}];
t3DQuarticSeed = FileNameJoin[{t3DData, "eft1_after_F_scalar_quartic_seed.json"}];
t3DNeeded = {t3DStatus, t3DBasis, t3DTreePath, t3DHardPath,
  t3DWilsonSeed, t3DQuarticSeed};
If[!And @@ (FileExistsQ /@ t3DNeeded),
  t3Stop["Missing verified three-flavour D physical F source",
    <|"missing" -> Select[t3DNeeded, !FileExistsQ[#] &]|>]];
t3DMeta = Quiet @ Check[Import[t3DStatus, "RawJSON"], $Failed];
t3DMass = Quiet @ Check[Import[t3DBasis, "RawJSON"], $Failed];
If[!AssociationQ[t3DMeta] || !AssociationQ[t3DMass] ||
   Lookup[t3DMeta, "status", ""] =!= "FThresholdTreeWilsonExported" ||
   Lookup[t3DMeta, "heavy_flavour_count", 0] =!= 3 ||
   !TrueQ[Lookup[t3DMeta, "raw_one_loop_F_matching_included", False]] ||
   Lookup[t3DMass, "mass_basis_kind", ""] =!=
     "RealSelfBilinearSplitLightFields" ||
   Lookup[t3DMass, "split_species", 0] =!= 2 ||
   Lookup[t3DMass, "minus_field", ""] =!= "T3Scalar2Minus" ||
   Lookup[t3DMass, "plus_field", ""] =!= "T3Scalar2Plus" ||
   TrueQ[Lookup[t3DMass, "complex_self_bilinear_phase_supported", True]],
  t3Stop["F source is not the validated NF=3 real-bilinear mass basis"]];

t3stage = "D-rebuild-Matchete-physical-registry";
If[CheckAbort[UsingFrontEnd[Needs["Matchete`"]; True], $Aborted] =!= True,
  t3Stop["Matchete could not be loaded for D(-2)"]];
Get[FileNameJoin[{t3root, "model", "T3ModelCatalog.wl"}]];
Get[FileNameJoin[{t3root, "model", "LagrangianBuilder.wl"}]];
Get[FileNameJoin[{t3root, "RunMatching.wl"}]];
t3DModel = T3ModelFromClass["D", -2, 3];
If[!AssociationQ[t3DModel] || SetT3HeavyFields[{"F"}] === $Failed,
  t3Stop["D(-2) three-flavour model/registry unavailable"]];
t3DBuild = CheckAbort[UsingFrontEnd[BuildT3Lagrangian[t3DModel]], $Aborted];
If[!AssociationQ[t3DBuild] ||
   Lookup[t3DBuild, "Status", ""] =!= "Success" ||
   !TrueQ[Lookup[t3DBuild, "FullValidation", False]] ||
   Lookup[t3DBuild, "AdditionalQuarticsIncluded", -1] =!= 5 ||
   Lookup[t3DBuild, "AdditionalQuadraticsIncluded", -1] =!= 1 ||
   Lookup[t3DBuild, "AdditionalYukawasIncluded", -1] =!= 1,
  t3Stop["D(-2) complete UV coupling inventory failed"]];
(* Rebuild field objects only. DO NOT repeat F matching. *)
RemoveField[NewFermion]; RemoveField[NewScalar1]; RemoveField[NewScalar2];
If[SetT3HeavyFields[{"S1"}] === $Failed ||
   DefineT3FieldObjects[t3DModel, False, {"S1"}] === $Failed,
  t3Stop["Cannot register D scalar1 at scalar threshold"]];
DefineField[T3Scalar2Minus, Scalar, SelfConjugate -> True,
  Mass -> {Heavy, T32MinusMass}];
DefineField[T3Scalar2Plus, Scalar, SelfConjugate -> True,
  Mass -> {Heavy, T32PlusMass}];
PreserveIntegratedT3MassCouplings[{"S1"}];
(* These are new independent EFT1 Wilson tensors, registered AFTER the
   fresh physical field basis exists. Do not Check[] warning messages. *)
Quiet[DefineCoupling[T3EFT1W12, Indices -> {Flavor, Flavor},
   SelfConjugate -> False]];
Quiet[DefineCoupling[T3EFT1WExtra, Indices -> {Flavor, Flavor},
   SelfConjugate -> False]];
t3DCouplings = GetCouplings[];
If[!AssociationQ[t3DCouplings] ||
   !And @@ (KeyExistsQ[t3DCouplings, #] & /@
     {T3EFT1W12, T3EFT1WExtra}),
  t3Stop["Could not register D(-2) physical EFT1 flavour-Wilson tensors"]];
t3DFields = GetFields[];
If[!AssociationQ[t3DFields] ||
   !And @@ (TrueQ[Lookup[Lookup[t3DFields, #, <||>], Heavy, False]] & /@
     {NewScalar1, T3Scalar2Minus, T3Scalar2Plus}) ||
   KeyExistsQ[t3DFields, NewScalar2] ||
   KeyExistsQ[t3DFields, NewFermion],
  t3Stop["D physical scalar heavy-field inventory is inconsistent"]];
If[DefineT3InvariantCGs[t3DModel] === $Failed ||
   !ListQ[BuildGeneralScalarPotentialCandidates[t3DModel]],
  t3Stop["D invariant CG reconstruction failed"]];

t3stage = "D-rehydrate-exact-three-flavour-EFT1";
t3DSeeds = Quiet @ Check[
  Import[#, "RawJSON"] & /@ {t3DWilsonSeed, t3DQuarticSeed}, $Failed];
If[!ListQ[t3DSeeds] || Length[t3DSeeds] =!= 2 ||
   !AllTrue[t3DSeeds, AssociationQ],
  t3Stop["D EFT1 CG seed metadata malformed"]];
t3DCGEntries = DeleteDuplicatesBy[
  Flatten[Lookup[t3DSeeds, "CGRegistry", {}], 1],
  Lookup[#, "Name", ""] &];
Do[
  t3DCGName = ToExpression[t3DCGItem["Name"], InputForm];
  t3DCGKnown = Quiet @ Check[
    Matchete`PackageScope`$CGproperties[t3DCGName, Indices], $Failed];
  If[t3DCGKnown === $Failed || MissingQ[t3DCGKnown],
    t3DCGReps = ToExpression[t3DCGItem["RepsInputForm"], InputForm];
    t3DCGTensor = ToExpression[t3DCGItem["TensorInputForm"], InputForm];
    If[CheckAbort[
       DefineCG[t3DCGName, t3DCGReps, t3DCGTensor]; True,
       $Aborted] =!= True,
      t3Stop["Could not restore D CG " <> t3DCGItem["Name"]]]
  ],
  {t3DCGItem, t3DCGEntries}
];
t3DTreeText = Import[t3DTreePath, "Text"];
t3DHardText = Import[t3DHardPath, "Text"];
If[!StringQ[t3DTreeText] || !StringQ[t3DHardText] ||
   !StringContainsQ[t3DTreeText, "T3Scalar2Minus"] ||
   !StringContainsQ[t3DTreeText, "T3Scalar2Plus"] ||
   !StringContainsQ[t3DHardText, "hbar"],
  t3Stop["D physical F source text inventory invalid"]];
(* Batch2 stores exact Matchete InputForm. This evaluation is done only
   after Matchete model and CG symbols have been registered. *)
t3DTree = Quiet @ Check[ToExpression[t3DTreeText, InputForm], $Failed];
If[t3DTree === $Failed || FreeQ[t3DTree, GammaCC],
  t3Stop["D saved first-threshold tree cannot be rehydrated"]];
t3DUsedCG = DeleteDuplicates @ Cases[HoldComplete[t3DTree],
  CG[name_, ___] :> Replace[name, Bar[x_] :> x], Infinity];
t3DUnregistered = Select[t3DUsedCG, Function[name,
  Module[{r = Quiet @ Check[
      Matchete`PackageScope`$CGproperties[name, Indices], $Failed]},
    r === $Failed || MissingQ[r]]]];
If[t3DUnregistered =!= {},
  t3Stop["D first-threshold tree contains unregistered CG tensors",
    <|"CGs" -> (ToString[#, InputForm] & /@ t3DUnregistered)|>]];
(* Preserve F hard source as an exact InputForm file; it is NOT inserted
   into a second one-loop match or discarded. *)
t3Save["D_full_F_hard_one_loop_preserved.txt", t3DHardText, "Text"];
t3Save["D_full_F_tree_original.wxf", t3DTree];

(* Contract the full indexed MF[r] inverses in all 8 LLSS monomials.
   The two different Yukawa sectors and both Hermitian conjugates remain. *)
ClearAll[t3DMassIndex, t3DYukawa, t3DContract];
t3DMassIndex[x_] := Replace[Unevaluated[x], {
  HoldPattern[Power[Coupling[MF, {r:Index[_, NFlavor]}, 0], -1]] :> r,
  _ :> Missing["NoMF"]}];
t3DYukawa[x_] := Replace[Unevaluated[x], {
  HoldPattern[Coupling[n:(y1|y2|y2Conj),
      {p:Index[_, Flavor], r:Index[_, NFlavor]}, 0]] :>
    <|"name" -> n, "p" -> p, "r" -> r, "barred" -> False|>,
  HoldPattern[Bar[Coupling[n:(y1|y2|y2Conj),
      {p:Index[_, Flavor], r:Index[_, NFlavor]}, 0]]] :>
    <|"name" -> n, "p" -> p, "r" -> r, "barred" -> True|>,
  _ :> Missing["NotYukawa"]}];
t3DWilsonCount = <|"baseline" -> 0, "extra" -> 0|>;
t3DContract[term_] := Module[
  {factors, masses, massPositions, flav, y1pos, y2pos, r,
   left, right, kind, w, remainder},
  factors = If[Head[term] === Times, List @@ term, {term}];
  masses = t3DMassIndex /@ factors;
  massPositions = Select[Range[Length[masses]], !MissingQ[masses[[#]]] &];
  If[massPositions === {}, Return[term]];
  If[Length[massPositions] =!= 1,
    t3Stop["D has ambiguous inverse heavy-mass index"]];
  r = masses[[First[massPositions]]];
  flav = t3DYukawa /@ factors;
  y1pos = Select[Range[Length[flav]], AssociationQ[flav[[#]]] &&
    flav[[#]]["name"] === y1 && flav[[#]]["r"] === r &];
  y2pos = Select[Range[Length[flav]], AssociationQ[flav[[#]]] &&
    MemberQ[{y2, y2Conj}, flav[[#]]["name"]] &&
    flav[[#]]["r"] === r &];
  If[Length[y1pos] =!= 1 || Length[y2pos] =!= 1,
    t3Stop["D inverse-MF term has no unique y1/y2 flavour contraction"]];
  left = flav[[First[y1pos]]]; right = flav[[First[y2pos]]];
  If[left["barred"] =!= right["barred"],
    t3Stop["D tree mixes barred/unbarred Yukawa sectors"]];
  kind = If[right["name"] === y2, "baseline", "extra"];
  w = Coupling[If[kind === "baseline", T3EFT1W12, T3EFT1WExtra],
    {left["p"], right["p"]}, 0];
  If[TrueQ[left["barred"]], w = Bar[w]];
  remainder = Delete[factors, List /@ Sort[{
    First[massPositions], First[y1pos], First[y2pos]}]];
  t3DWilsonCount[kind]++;
  Times @@ Append[remainder, w]
];
t3DTerms = If[Head[t3DTree] === Plus, List @@ t3DTree, {t3DTree}];
t3DContracted = Plus @@ (t3DContract /@ t3DTerms);
If[t3DWilsonCount =!= <|"baseline" -> 4, "extra" -> 4|> ||
   FreeQ[t3DContracted, T3EFT1W12] ||
   FreeQ[t3DContracted, T3EFT1WExtra] ||
   !FreeQ[t3DContracted,
     HoldPattern[Coupling[MF, {Index[_, NFlavor]}, 0]]] ||
   !FreeQ[t3DContracted, NewFermion],
  t3Stop["D tree flavour Wilson contract not complete",
    <|"contract_counts" -> t3DWilsonCount|>]];
t3Save["D_EFT1_tree_contracted.wxf", t3DContracted];
t3Save["D_EFT1_tree_contracted.txt",
  ToString[t3DContracted, InputForm], "Text"];
t3Save["D_EFT1_Wilson_definitions.json", <|
  "W12[p,q]" -> "Sum[y1[p,r]*y2[q,r]/MF[r],{r,1,3}]",
  "WExtra[p,q]" -> "Sum[y1[p,r]*y2Conj[q,r]/MF[r],{r,1,3}]",
  "mass_basis" -> {"NewScalar1", "T3Scalar2Minus", "T3Scalar2Plus"},
  "three_independent_fermion_flavours" -> True,
  "one_loop_hard_reparameterised" -> False,
  "physical_self_bilinear_complex_phase_supported" -> False|>, "RawJSON"];

(* Scalar masses were light at F threshold (EFT order 1) and become
   heavy at scalar threshold (EFT order 0). Only these three mass labels
   change power counting. No numerical scalar spectrum is invented. *)
t3stage = "D-promote-all-physical-scalar-masses";
t3DPhysicalMassNames = {MS1, T32MinusMass, T32PlusMass};
t3DMassCounts = AssociationThread[
  (ToString[#, InputForm] & /@ t3DPhysicalMassNames),
  (With[{m = #}, Count[t3DContracted,
       HoldPattern[Coupling[m, {}, 1]], Infinity]] & /@
    t3DPhysicalMassNames)];
If[!AllTrue[Values[t3DMassCounts], # > 0 &],
  t3Stop["D source is missing one or more scalar mass symbols",
    <|"original_scalar_mass_occurrences" -> t3DMassCounts|>]];
t3DMatchInput = t3DContracted /. {
  HoldPattern[Coupling[MS1, {}, 1]] :> Coupling[MS1, {}, 0],
  HoldPattern[Coupling[T32MinusMass, {}, 1]] :> Coupling[T32MinusMass, {}, 0],
  HoldPattern[Coupling[T32PlusMass, {}, 1]] :> Coupling[T32PlusMass, {}, 0]};
t3DDemoted = t3DMatchInput /. {
  HoldPattern[Coupling[MS1, {}, 0]] :> Coupling[MS1, {}, 1],
  HoldPattern[Coupling[T32MinusMass, {}, 0]] :> Coupling[T32MinusMass, {}, 1],
  HoldPattern[Coupling[T32PlusMass, {}, 0]] :> Coupling[T32PlusMass, {}, 1]};
If[!TrueQ[t3DDemoted === t3DContracted] ||
   AnyTrue[t3DPhysicalMassNames, Function[m,
     !FreeQ[t3DMatchInput, Coupling[m, {}, 1]]]],
  t3Stop["D scalar power-counting mass promotion not reversible"]];
t3Save["D_EFT1_physical_scalar_match_input.wxf", t3DMatchInput];

t3stage = "D-physical-scalar-CheckLagrangian";
t3DChecked = TimeConstrained[
  CheckAbort[UsingFrontEnd[CheckLagrangian[t3DMatchInput]], $Aborted],
  2400, $TimedOut];
If[t3DChecked =!= True,
  t3Stop["D physical scalar EFT1 CheckLagrangian failed",
    <|"check_result" -> ToString[t3DChecked, InputForm],
      "reused_F_boundary" -> t3DTreePath|>]];
t3stage = "D-joint-physical-scalar-soft-one-loop";
t3DMatched = TimeConstrained[
  CheckAbort[UsingFrontEnd[
    Match[t3DMatchInput, EFTOrder -> 5, LoopOrder -> 1]], $Aborted],
  3600, $TimedOut];
If[MemberQ[{$Failed, $Aborted, $TimedOut}, t3DMatched] ||
   AtomQ[t3DMatched] || FreeQ[t3DMatched, hbar] ||
   !FreeQ[t3DMatched,
     Field[(NewScalar1|T3Scalar2Minus|T3Scalar2Plus), ___]],
  t3Stop["D joint physical scalar Match returned invalid EFT or left heavy scalars",
    <|"result_head" -> ToString[Head[t3DMatched], InputForm]|>]];
t3Save["D_scalar_soft_only_raw_matched.wxf", t3DMatched];
(* Standard Matchete postprocessing, each part retained; no hard F source
   or hbar-order running is spliced into this new one-loop Match. *)
t3DProcessed = t3DMatched;
Do[
  t3DProcessed = TimeConstrained[CheckAbort[UsingFrontEnd[
    Switch[step,
      "GreensSimplify", GreensSimplify[t3DProcessed],
      "EOMSimplify", EOMSimplify[t3DProcessed],
      "EvaluateLoopFunctions", EvaluateLoopFunctions[t3DProcessed],
      "ReplaceEffectiveCouplings", ReplaceEffectiveCouplings[t3DProcessed]]],
    $Aborted], 2400, $TimedOut];
  If[MemberQ[{$Failed, $Aborted, $TimedOut}, t3DProcessed],
    t3Stop["D scalar-soft postprocessing failed at " <> step]];
  t3Save["D_scalar_soft_only_" <> step <> ".wxf", t3DProcessed],
  {step, {"GreensSimplify", "EOMSimplify", "EvaluateLoopFunctions",
    "ReplaceEffectiveCouplings"}}
];
t3DOneLoop = hbar Coefficient[Expand[t3DProcessed], hbar, 1];
If[t3DOneLoop === 0,
  t3Stop["D physical scalar result has no O(hbar) soft contribution"]];
t3DWeinberg = CheckAbort[
  ExtractWeinbergCoefficient[t3DOneLoop], $Aborted];
If[!AssociationQ[t3DWeinberg] ||
   Lookup[t3DWeinberg, "Status", ""] =!= "Success" ||
   Lookup[t3DWeinberg, "TermCount", 0] <= 0,
  t3Stop["D physical scalar soft Weinberg extraction failed"]];
t3Save["D_scalar_soft_only_unrenormalized.wxf",
  t3DWeinberg["Coefficient"]];
t3Save["D_scalar_soft_only_unrenormalized.txt",
  ToString[t3DWeinberg["Coefficient"], InputForm], "Text"];
t3stage = "D-physical-soft-success-hard-inheritance-pending";
t3Write["D_PhysicalScalarSoftOnly_HardInheritanceAndWilsonRGPending",
  "Three-flavour physical scalar threshold matched; F hard, EFT1 Wilson RG and finite renormalisation are absent. NOT physical C5.",
  <|"first_threshold_reused" -> True,
    "F_source_three_flavours" -> True,
    "D_physical_mass_basis" -> {"NewScalar1", "T3Scalar2Minus", "T3Scalar2Plus"},
    "joint_scalar_match_performed" -> True,
    "wilson_tree_contractions" -> t3DWilsonCount,
    "physical_scalar_soft_Weinberg_terms" -> t3DWeinberg["TermCount"],
    "D_hard_source_text" -> "D_full_F_hard_one_loop_preserved.txt",
    "D_soft_source_WXF" -> "D_scalar_soft_only_unrenormalized.wxf"|>];
Exit[5];
