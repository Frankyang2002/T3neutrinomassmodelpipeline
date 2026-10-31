(* Workstream A: live, read-only-source Matchete two-threshold experiment.
   Full Lagrangian and extra/zero variants are matched in independent Wolfram
   processes. This script does not run EFT1 Wilson RGEs or renormalize poles.
   Usage: wolframscript -file Lagrangian/RunT3ExtendedHierarchyNative.wl
          OUT_DIRECTORY A|D|E m2|p0 seeded|zero stage1|hierarchy
   All exported coefficients are in native Matchete coordinates. *)
ClearAll["Global`*"];
hArgs = Rest[$ScriptCommandLine];
If[Length[hArgs] =!= 5,
  Print["Usage: OUT_DIR A|D|E m2|p0 seeded|zero stage1|hierarchy"];
  Exit[2]
];
hOut = ExpandFileName[hArgs[[1]]];
hClass = ToUpperCase[hArgs[[2]]];
hTag = hArgs[[3]];
hVariant = hArgs[[4]];
hMode = hArgs[[5]];
hModel = hClass <> "_" <> hTag;
hAllowed = <|"A_m2" -> -2, "A_p0" -> 0,
  "D_m2" -> -2, "D_p0" -> 0,
  "E_m2" -> -2, "E_p0" -> 0|>;
If[!KeyExistsQ[hAllowed, hModel] ||
   !MemberQ[{"seeded", "zero"}, hVariant] ||
   !MemberQ[{"stage1", "hierarchy"}, hMode],
  Print["Invalid model/variant/mode"];
  Exit[2]
];
If[!And @@ (Environment[#] === "1" & /@ {
    "T3_EXPERIMENTAL_GENERAL_FERMIONS",
    "T3_EXPERIMENTAL_FLAVOR_MASSES",
    "T3_EXPERIMENTAL_GENERAL_QUARTICS"}) ||
    Environment["T3_EXPERIMENTAL_GENERAL_QUADRATICS"] === "1",
  Print["Required general fermions/flavour masses/quartics=1; quadratics !=1"];
  Exit[2]
];
If[!DirectoryQ[hOut], CreateDirectory[hOut, CreateIntermediateDirectories -> True]];
hRoot = DirectoryName @ ExpandFileName[$InputFileName];
hStep = "start";
ClearAll[HStep, HText, HExport, HFail, HResult];
HText[x_] := ToString[x, InputForm];
HStep[s_String] := (hStep = s; Print["EXT_HIERARCHY_STAGE: ",s];
  Export[FileNameJoin[{hOut,"progress.json"}],
    <|"model"->hModel,"variant"->hVariant,"mode"->hMode,"stage"->s,
      "hierarchical_c5_physically_certified"->False|>,"RawJSON"];
  s);
HExport[name_String, expression_] := Module[{p},
  p = FileNameJoin[{hOut, name <> ".wxf"}];
  If[Quiet@Check[Export[p, expression, "WXF"], $Failed] === $Failed ||
      !FileExistsQ[p], HFail["WXF export failed: " <> name]];
  p
];
HFail[msg_String] := Throw[<|"status"->"NativeHierarchyBlocked",
    "failed_step"->hStep,"reason"->msg|>,"HierarchyStop"];
HResult[obj_Association] := Module[{r},
  r = Join[<|
    "format"->"t3-extra-yukawa-hierarchical-native-v1",
    "model"->hModel,"alpha"->hAllowed[hModel],
    "variant"->hVariant,"mode"->hMode,"heavy_flavours"->3,
    "last_step"->hStep,
    "fermion_threshold_integrated_out"->False,
    "scalar_threshold_integrated_out"->False,
    "eft1_Wilson_RGE_evolved"->False,
    "uncancelled_poles_excluded"->False,
    "absolute_matchete_rgbeta_dictionary_certified"->False,
    "physical_low_energy_prediction_certified"->False,
    "full_scalar_quadratic_basis_included"->False|>, obj];
  Export[FileNameJoin[{hOut,"native.json"}],r,"RawJSON"];
  Print[hModel,"/",hVariant,"/",hMode,": ",Lookup[r,"status","Unknown"]];
  r
];
HStep["load-Matchete"];
(* Separate completed top-level load: Bar and Coupling below must resolve to
   Matchete` symbols, not symbols created prematurely in Global`. *)
hLoaded = CheckAbort[Check[UsingFrontEnd[Needs["Matchete`"]; True],$Failed],$Aborted];
If[hLoaded =!= True,
  HResult[<|"status"->"NativeHierarchyBlocked",
    "reason"->"Matchete package load failed"|>]; Exit[3]
];
(* The entire following expression is parsed AFTER Needs has finished. *)
hRun = CheckAbort[Catch[Module[
  {model, build, uv, extra, rules, plan, matched, stages, first,
   final, c5Data, c5Raw, out, last, success, poleQ, nativeC5,
   stage1Tree, stage1Loop, stage1Full, finalEft, stageCount,
   freshStage, freshPath, poleDiagnostic, eft1RunningLoaded},
  Get[FileNameJoin[{hRoot,"model","T3ModelCatalog.wl"}]];
  Get[FileNameJoin[{hRoot,"model","LagrangianBuilder.wl"}]];
  Get[FileNameJoin[{hRoot,"RunMatching.wl"}]];
  model = T3ModelFromClass[hClass,hAllowed[hModel],3];
  If[!AssociationQ[model],HFail["Invalid T3 representation"]];
  HStep["configure-F-only-heavy"];
  If[SetT3HeavyFields[{"F"}] === $Failed,
    HFail["Cannot set first heavy threshold to F only"]];
  HStep["build-and-check-native-Lagrangian"];
  build = CheckAbort[UsingFrontEnd[BuildT3Lagrangian[model]],$Aborted];
  If[!AssociationQ[build] || Lookup[build,"Status",""] =!= "Success" ||
     !TrueQ[Lookup[build,"FullValidation",False]] ||
     Lookup[build,"AdditionalYukawasIncluded",0] =!= 1 ||
     !TrueQ[Lookup[build,"MatchingMassBasisReady",False]],
    HFail["Native model did not validate exactly one extra Yukawa with F heavy"]];
  extra = If[hAllowed[hModel] === 0, y1Conj, y2Conj];
  uv = build["LUV"];
  If[FreeQ[uv,extra],HFail["Requested extra Yukawa absent from native UV action"]];
  If[hVariant === "zero",
    HStep["zero-extra-at-UV-action-before-matching"];
    rules = {
      HoldPattern[Bar[Coupling[y1Conj,___]]] :> 0,
      HoldPattern[Bar[Coupling[y2Conj,___]]] :> 0,
      HoldPattern[Coupling[y1Conj,___]] :> 0,
      HoldPattern[Coupling[y2Conj,___]] :> 0};
    uv = Expand[uv /. rules];
    If[!FreeQ[uv,y1Conj] || !FreeQ[uv,y2Conj],
      HFail["Zero-extra UV projection left additional Yukawa symbols"]]
  ];
  HStep["native-fermion-threshold-Matchete"];
  plan = If[hMode === "hierarchy",{{"F"},{"S1","S2"}},{{"F"}}];
  matched = CheckAbort[UsingFrontEnd[
    RunSequentialT3Matching[uv,model,plan,5,1,hOut,hModel]],$Aborted];
  If[!AssociationQ[matched],
    HFail["Sequential matching aborted without an association"]];
  stages = Lookup[matched,"Stages",{}];
  stageCount = Length[stages];
  If[stageCount >= 1,
    first = First[stages];
    If[!And @@ (KeyExistsQ[first,#] & /@ {
        "TreeLagrangian","OneLoopCorrection","CompactLagrangian"}),
      HFail["Missing tree/loop split in stage-one matching"]];
    stage1Tree = HExport["stage1_tree",first["TreeLagrangian"]];
    stage1Loop = HExport["stage1_one_loop",first["OneLoopCorrection"]];
    stage1Full = HExport["stage1_full",first["CompactLagrangian"]]
  ];
  If[Lookup[matched,"Status",""] =!= "Success",
    HFail["Sequential matching returned " <> ToString[Lookup[matched,"Status","Unknown"]] <>
      "; completed stage count=" <> ToString[stageCount] <>
      "; see fresh_kernel_stage_2/stdout.log if present"]];
  If[stageCount =!= Length[plan],HFail["Incomplete threshold stage count"]];
  If[hMode === "stage1",
    HStep["complete-stage1"];
    Return[<|"status"->"FermionEFT1NativeMatched",
      "fermion_threshold_integrated_out"->True,
      "scalar_threshold_integrated_out"->False,
      "stage_count"->stageCount,"stage1_tree_wxf"->stage1Tree,
      "stage1_one_loop_wxf"->stage1Loop,"stage1_full_wxf"->stage1Full,
      "tree_one_loop_split_retained"->True,
      "matching_performed_in_Matchete"->True|>]
  ];
  HStep["extract-stage2-Weinberg-sector"];
  (* The authoritative native stage-2 result is saved by the ORIGINAL fresh
     threshold driver. Preserve its pole/RGE and EFT1-insertion diagnostics;
     a successful Matchete Match alone does not cancel UV poles. *)
  freshPath = FileNameJoin[{hOut,"debug","fresh_kernel_stage_2","result.wxf"}];
  freshStage = Quiet@Check[Import[freshPath,"WXF"],$Failed];
  If[!AssociationQ[freshStage] || Lookup[freshStage,"Status",""] =!= "Success",
    HFail["Sequential second stage missing native fresh-kernel provenance"]];
  eft1RunningLoaded = TrueQ[Lookup[freshStage,"EFT1WilsonRunningLoaded",False]];
  poleDiagnostic = Lookup[freshStage,"C5PoleRGEConsistency",<||>];
  If[!AssociationQ[poleDiagnostic],poleDiagnostic=<||>];
  final = Last[stages];
  If[!KeyExistsQ[final,"Lagrangian"],HFail["Missing final Matchete EFT"]];
  finalEft = HExport["stage2_full",final["Lagrangian"]];
  c5Data = CheckAbort[ExtractWeinbergCoefficient[final["Lagrangian"]],$Aborted];
  If[!AssociationQ[c5Data] || Lookup[c5Data,"Status",""] =!= "Success",
    HFail["Two threshold stages finished but Weinberg extraction is " <>
      ToString[If[AssociationQ[c5Data],Lookup[c5Data,"Status","Unknown"],c5Data]]]];
  c5Raw = c5Data["Coefficient"];
  nativeC5 = HExport["stage2_c5_ordered",c5Raw];
  Export[FileNameJoin[{hOut,"stage2_c5_ordered.inputform.txt"}],
    HText[c5Raw],"Text"];
  poleQ = !FreeQ[c5Raw, \[Epsilon]];
  HStep["complete-hierarchical-native-expression"];
  <|"status"->"NativeHierarchicalC5ExtractedNotRenormalized",
    "fermion_threshold_integrated_out"->True,
    "scalar_threshold_integrated_out"->True,
    "stage_count"->stageCount,"stage1_tree_wxf"->stage1Tree,
    "stage1_one_loop_wxf"->stage1Loop,
    "stage1_full_wxf"->stage1Full,
    "stage2_full_wxf"->finalEft,
    "stage2_c5_ordered_wxf"->nativeC5,
    "stage2_c5_contains_epsilon_pole"->poleQ,
    "stage2_c5_is_finite_renormalized"->False,
    "EFT1_Wilson_running_insertion_added"->eft1RunningLoaded,
    "fresh_kernel_stage2_result_wxf"->freshPath,
    "stage2_native_pole_rge_diagnostic"->poleDiagnostic,
    "stage2_direct_log_equals_twice_hard_pole_residue"->
      Lookup[poleDiagnostic,"DirectLogEqualsTwicePoleResidueOneGeneration",False],
    "stage2_C5_pole_rge_consistency_path"->
      Lookup[freshStage,"C5PoleRGEConsistencyPath","Unavailable"],
    "matching_performed_in_Matchete"->True|>
],"HierarchyStop"],$Aborted];
If[!AssociationQ[hRun], hRun = <|"status"->"NativeHierarchyBlocked",
  "reason"->"Native Matchete expression aborted or returned an invalid result",
  "failed_step"->hStep|>];
hSaved = HResult[hRun];
If[!MemberQ[{"FermionEFT1NativeMatched",
   "NativeHierarchicalC5ExtractedNotRenormalized"},
   Lookup[hSaved,"status",""]], Exit[6]];
Exit[0];
