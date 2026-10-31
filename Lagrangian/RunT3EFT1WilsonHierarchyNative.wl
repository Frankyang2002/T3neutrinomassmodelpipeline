(* Local-only, experimental two-threshold EFT1 flavour-tensor matching.
   Integrates all three F flavours in ONE Matchete matching step, then
   contracts the resulting tree-level y1*y2/MF[r] and extra*y2/MF[r]
   monomials into two EFT1 Wilson couplings W[p,q] before the scalar step.
   Stage-1 one-loop hard terms and EFT1 Wilson RGEs are NOT included in
   the stage-2 soft contribution. No finite physical C5 is claimed.

   wolframscript -file Lagrangian/RunT3EFT1WilsonHierarchyNative.wl
      OUTPUT_DIR A|D|E m2|p0 seeded|zero prepare|tree|soft
*)
ClearAll["Global`*"];
wArgs = Rest[$ScriptCommandLine];
If[Length[wArgs] =!= 5,
  Print["Usage: OUT_DIR A|D|E m2|p0 seeded|zero prepare|tree|soft"];
  Exit[2]];
wOut = ExpandFileName[wArgs[[1]]];
wClass = ToUpperCase[wArgs[[2]]]; wTag = wArgs[[3]];
wVariant = wArgs[[4]]; wMode = wArgs[[5]];
wModel = wClass <> "_" <> wTag;
wAlpha = <|"A_p0"->0,"A_m2"->-2,"D_p0"->0,"D_m2"->-2,
  "E_p0"->0,"E_m2"->-2|>;
If[!KeyExistsQ[wAlpha,wModel] || !MemberQ[{"seeded","zero"},wVariant] ||
    !MemberQ[{"prepare","tree","soft"},wMode],
  Print["Invalid model/variant/mode"];Exit[2]];
If[Environment["T3_EXPERIMENTAL_GENERAL_FERMIONS"] =!= "1" ||
   Environment["T3_EXPERIMENTAL_FLAVOR_MASSES"] =!= "1" ||
   Environment["T3_EXPERIMENTAL_GENERAL_QUARTICS"] =!= "1" ||
   Environment["T3_EXPERIMENTAL_GENERAL_QUADRATICS"] === "1",
   Print["Required experimental flavours/fermions/quartics=1, quadratics=0"];
   Exit[2]];
If[!DirectoryQ[wOut],CreateDirectory[wOut,CreateIntermediateDirectories->True]];
wRoot=DirectoryName @ ExpandFileName[$InputFileName];
wLast="initial";
ClearAll[WCheckpoint,WStop,WExport,WReport];
WCheckpoint[step_String] := (wLast=step;Print["EFT1_WILSON_STAGE: ",step];
  Export[FileNameJoin[{wOut,"progress.json"}],
    <|"model"->wModel,"variant"->wVariant,"mode"->wMode,
      "last_step"->wLast,"full_hierarchical_c5_certified"->False|>,"RawJSON"];
  step);
WStop[reason_String] := Throw[<|"status"->"Blocked",
    "reason"->reason,"failed_step"->wLast|>,"EFT1WStop"];
WExport[file_String,data_,format_String:"WXF"] := Module[{path},
  path=FileNameJoin[{wOut,file}];
  If[Quiet@Check[Export[path,data,format],$Failed] === $Failed || !FileExistsQ[path],
    WStop["Failed export: "<>file]];
  path];
WReport[data_Association] := Module[{record},
  record=Join[<|
    "format"->"t3-eft1-flavour-wilson-soft-v1",
    "model"->wModel,"variant"->wVariant,"mode"->wMode,
    "heavy_generations_integrated_together"->3,
    "independent_heavy_masses_preserved_in_wilson_definition"->True,
    "tree_level_eft1_wilson_contraction_applied"->False,
    "stage1_hard_one_loop_inherited_at_stage2"->False,
    "eft1_wilson_running_included"->False,
    "renormalized_hierarchical_c5_certified"->False,
    "absolute_matchete_rgbeta_dictionary_certified"->False,
    "physical_neutrino_prediction_certified"->False,
    "last_step"->wLast|>,data];
  Export[FileNameJoin[{wOut,"native.json"}],record,"RawJSON"];
  Print[wModel,"/",wVariant,"/",wMode,": ",Lookup[record,"status","Unknown"]];
  record];

WCheckpoint["matchete-load"];
wLoaded=CheckAbort[Check[UsingFrontEnd[Needs["Matchete`"];True],$Failed],$Aborted];
If[wLoaded =!= True,
  WReport[<|"status"->"Blocked","reason"->"Cannot load Matchete"|>];Exit[3]];

(* All Matchete symbols in the following Module are parsed AFTER Needs. *)
wOutput = CheckAbort[Catch[Module[
 {model,build,uv,seedRules,originalStage,steps,stage1,originalTree,
  originalOneLoop,transitionTree,treeContracted,firstCounts,secondCounts,
  f1,f2,mass,heavy,candidates,left,right,lepton1,lepton2,
  matching,softLoop,stage2,c5,outFile,def,definitionPath,
  stage1Hbar,zeroed,stage1Result,record},
 Get[FileNameJoin[{wRoot,"model","T3ModelCatalog.wl"}]];
 Get[FileNameJoin[{wRoot,"model","LagrangianBuilder.wl"}]];
 Get[FileNameJoin[{wRoot,"RunMatching.wl"}]];
 Get[FileNameJoin[{wRoot,"RunT3EFT1WilsonStageLauncher.wl"}]];
 model=T3ModelFromClass[wClass,wAlpha[wModel],3];
 If[!AssociationQ[model],WStop["Model construction failed"]];
 WCheckpoint["build-UV-with-F-heavy"];
 If[SetT3HeavyFields[{"F"}] === $Failed,
    WStop["Cannot select common fermion matching group"]];
 build=CheckAbort[UsingFrontEnd[BuildT3Lagrangian[model]],$Aborted];
 If[!AssociationQ[build] || Lookup[build,"Status",""] =!= "Success" ||
   !TrueQ[Lookup[build,"FullValidation",False]] ||
   Lookup[build,"AdditionalYukawasIncluded",0] =!= 1,
   WStop["Native UV builder did not accept extra Yukawa"]];
 uv=build["LUV"];
 If[wVariant === "zero",
   WCheckpoint["zero-extra-at-UV"];
   seedRules={
     HoldPattern[Bar[Coupling[y1Conj,___]]]:>0,
     HoldPattern[Coupling[y1Conj,___]]:>0,
     HoldPattern[Bar[Coupling[y2Conj,___]]]:>0,
     HoldPattern[Coupling[y2Conj,___]]:>0};
   uv=Expand[uv/.seedRules];
   If[!FreeQ[uv,y1Conj|y2Conj],WStop["Zero-extra UV projection failed"]]
 ];
 WCheckpoint["F-only-Matchete-stage1-three-heavy-flavours"];
 stage1Result=CheckAbort[UsingFrontEnd[RunSequentialT3Matching[
    uv,model,{{"F"}},5,1,wOut,wModel]],$Aborted];
 If[!AssociationQ[stage1Result] || Lookup[stage1Result,"Status",""] =!= "Success" ||
    Length[Lookup[stage1Result,"Stages",{}]] =!= 1,
    WStop["Native F-threshold matching failed"]];
 stage1=First[stage1Result["Stages"]];
 originalTree=stage1["TreeLagrangian"];
 originalOneLoop=stage1["OneLoopCorrection"];
 transitionTree=Lookup[stage1,"TransitionTreeLagrangian",originalTree];
 WExport["stage1_tree_original.wxf",originalTree];
 WExport["stage1_hard_one_loop_omitted_from_stage2.wxf",originalOneLoop];
 WExport["stage1_transition_tree_original.wxf",transitionTree];

 (* W[p,q] is an INDEPENDENT EFT1 coefficient.  It is deliberately not
    replaced with 3 Yukawas and an indexed inverse mass when passed to
    Matchete. The definitions below retain three distinct MF[r]. *)
 WCheckpoint["define-two-flavour-Wilson-couplings"];
 If[Quiet@Check[
     DefineCoupling[T3EFT1W12,Indices->{Flavor,Flavor},SelfConjugate->False];
     DefineCoupling[T3EFT1WExtra,Indices->{Flavor,Flavor},SelfConjugate->False];
     True,False] =!= True,
   WStop["Failed to register EFT1 tensor couplings"]];

 ClearAll[WIndexMass,WParseYukawa,WReplaceTerm,WContractTree];
 WIndexMass[expr_] := Replace[Unevaluated[expr],{
    HoldPattern[Power[Coupling[MF,{r:Index[_,NFlavor]},0],-1]] :> r,
    _ :> Missing["NotIndexedInverseMass"]}];
 WParseYukawa[expr_] := Replace[Unevaluated[expr],{
    HoldPattern[Coupling[n:(y1|y2|y1Conj|y2Conj),
      {p:Index[_,Flavor],r:Index[_,NFlavor]},0]] :>
      <|"name"->n,"p"->p,"r"->r,"barred"->False|>,
    HoldPattern[Bar[Coupling[n:(y1|y2|y1Conj|y2Conj),
      {p:Index[_,Flavor],r:Index[_,NFlavor]},0]]] :>
      <|"name"->n,"p"->p,"r"->r,"barred"->True|>,
    _ :> Missing["NotYukawa"]}];
 firstCounts=<|"baseline"->0,"extra"->0,"unmatched_mass"->0|>;
 WReplaceTerm[term_] := Module[
   {factors,inds,masses,ix,parsed,y1cand,y2cand,positions,
    leftY,rightY,wName,isBar,wCoupling,otherFactors,isExtra},
   factors=If[Head[term] === Times,List@@term,{term}];
   inds=WIndexMass /@ factors;
   masses=Select[Range[Length[inds]],!MissingQ[inds[[#]]] &];
   If[masses==={},Return[term]];
   If[Length[masses] =!= 1,
     firstCounts["unmatched_mass"]++;
     WStop["Unexpected number of MF[r] inverse factors in tree monomial"]];
   ix=First[masses];heavy=inds[[ix]];
   parsed=WParseYukawa /@ factors;
   y1cand=Select[Range[Length[parsed]],
     AssociationQ[parsed[[#]]] && Lookup[parsed[[#]],"r",None]===heavy &&
       MemberQ[{y1,y1Conj},Lookup[parsed[[#]],"name",None]] &];
   y2cand=Select[Range[Length[parsed]],
     AssociationQ[parsed[[#]]] && Lookup[parsed[[#]],"r",None]===heavy &&
       MemberQ[{y2,y2Conj},Lookup[parsed[[#]],"name",None]] &];
   If[Length[y1cand] =!= 1 || Length[y2cand] =!= 1,
     firstCounts["unmatched_mass"]++;
     WStop["MF[r] did not have exactly one y1 and one y2 factor"]];
   leftY=parsed[[First[y1cand]]];rightY=parsed[[First[y2cand]]];
   If[leftY["barred"] =!= rightY["barred"],
     WStop["Mixed Yukawa complex-conjugation channels"]];
   isExtra=MemberQ[{y1Conj,y2Conj},leftY["name"]] ||
           MemberQ[{y1Conj,y2Conj},rightY["name"]];
   If[isExtra && !(If[wTag==="p0",leftY["name"]===y1Conj,
                           rightY["name"]===y2Conj]),
      WStop["Unexpected extra-Yukawa orientation"]];
   If[!isExtra && (leftY["name"] =!= y1 || rightY["name"] =!= y2),
      WStop["Unexpected baseline Yukawa orientation"]];
   wName=If[isExtra,T3EFT1WExtra,T3EFT1W12];
   wCoupling=Coupling[wName,{leftY["p"],rightY["p"]},0];
   If[TrueQ[leftY["barred"]],wCoupling=Bar[wCoupling]];
   positions=Sort[{ix,First[y1cand],First[y2cand]}];
   otherFactors=Delete[factors,List/@positions];
   If[isExtra,firstCounts["extra"]++,firstCounts["baseline"]++];
   Times@@Append[otherFactors,wCoupling]
 ];
 WContractTree[expr_] := Module[{terms,converted},
   terms=If[Head[expr]===Plus,List@@expr,{expr}];
   converted=WReplaceTerm/@terms;
   Plus@@converted];
 WCheckpoint["contract-F-index-into-EFT1-W-pq"];
 treeContracted=WContractTree[transitionTree];
 If[firstCounts["unmatched_mass"] =!= 0 ||
    firstCounts["baseline"] =!= 2 ||
    firstCounts["extra"] =!= If[wVariant==="seeded",2,0],
    WStop["Unexpected tree-level Wilson contraction counts: "<>ToString[firstCounts,InputForm]]];
 If[!FreeQ[treeContracted,Coupling[MF,{Index[_,NFlavor]},0]],
   WStop["Indexed inverse MF remains after tree Wilson contraction"]];
 If[!FreeQ[treeContracted,T3EFT1W12] ||
   (wVariant==="seeded" && FreeQ[treeContracted,T3EFT1WExtra]) ||
   (wVariant==="zero" && !FreeQ[treeContracted,T3EFT1WExtra]),
   WStop["Wilson tensor inventory failed"]];
 WExport["stage1_tree_wilson.wxf",treeContracted];
 WExport["stage1_tree_wilson.inputform.txt",ToString[InputForm[treeContracted]],"Text"];
 def=<|
  "T3EFT1W12[p,q]"->
    "Sum[y1[p,r] y2[q,r] / MF[r], {r,1,3}]",
  "T3EFT1WExtra[p,q]"->If[wTag==="p0",
    "Sum[y1Conj[p,r] y2[q,r] / MF[r], {r,1,3}]",
    "Sum[y1[p,r] y2Conj[q,r] / MF[r], {r,1,3}]"],
  "conjugation"->"Bar[W[p,q]] = Conjugate[W[p,q]] for real positive MF[r]",
  "mass_assumption"->"Three diagonal independent masses, integrated at one common F threshold",
  "normalisation"->"Exact reparameterization of Matchete stage-1 tree terms, not RGBeta mapping"|>;
 definitionPath=WExport["wilson_definitions.json",def,"RawJSON"];
 If[wMode==="prepare",
   WCheckpoint["prepared-Wilson-tree-EFT1"];
   Return[<|"status"->"FermionThresholdWilsonTreeContracted",
      "tree_level_eft1_wilson_contraction_applied"->True,
      "stage1_fermion_threshold_passed"->True,
      "wilson_counts"->firstCounts,
      "wilson_definitions_path"->definitionPath|>]
 ];
 WCheckpoint["fresh-kernel-scalar-threshold-with-Wilson-tensors"];
 (* Only the tree-level F-threshold Wilson EFT is propagated.
    The existing original stage-1 hard O(hbar) correction has been saved,
    but is NOT incorrectly treated as though it had already been
    projected into an EFT1 Wilson basis. *)
 softLoop=If[wMode==="soft",1,0];
 stage2=CheckAbort[RunT3EFT1WilsonFreshThresholdStage[
   treeContracted,0,treeContracted,treeContracted,
   model,{"S1","S2"},{"S1","S2"},2,5,softLoop,wOut],$Aborted];
 If[!AssociationQ[stage2] || Lookup[stage2,"Status",""] =!= "Success",
   WStop["Fresh scalar stage failed; inspect debug/fresh_kernel_stage_2/stdout.log"]];
 WCheckpoint["extract-scalar-threshold-Weinberg-sector"];
 c5=Lookup[stage2,"AuthoritativeThresholdC5",$Failed];
 If[c5===$Failed,WStop["Fresh scalar stage omitted authoritative Weinberg coefficient"]];
 If[wMode==="soft" &&
    (TrueQ[PossibleZeroQ[c5]] || FreeQ[c5,T3EFT1W12]),
    WStop["Scalar one-loop C5 has no verified baseline EFT1 Wilson term"]];
 WExport["stage2_weinberg_unrenormalized.wxf",c5];
 WExport["stage2_weinberg_unrenormalized.inputform.txt",ToString[InputForm[c5]],"Text"];
 WCheckpoint["experimental-soft-scalar-matching-complete"];
 <|"status"->"ScalarThresholdWilsonSoftMatchedNotRenormalized",
    "tree_level_eft1_wilson_contraction_applied"->True,
    "stage1_fermion_threshold_passed"->True,
    "scalar_threshold_native_passed"->True,
    "stage2_loop_order"->softLoop,
    "wilson_counts"->firstCounts,
    "wilson_definitions_path"->definitionPath,
    "stage2_c5_has_epsilon"->(!FreeQ[c5,\[Epsilon]]),
    "stage2_c5_source"->"Tree EFT1 Wilson matching with new scalar threshold only",
    "eft1_wilson_running_insertion_added"->False,
    "stage1_hard_one_loop_inherited_at_stage2"->False,
    "stage2_fresh_result_path"->FileNameJoin[{wOut,"debug","fresh_kernel_stage_2","result.wxf"}]|>
 ],"EFT1WStop"],$Aborted];
If[!AssociationQ[wOutput],wOutput=<|"status"->"Blocked",
  "reason"->"Wolfram abort or invalid native expression","failed_step"->wLast|>];
wSaved=WReport[wOutput];
If[!MemberQ[{"FermionThresholdWilsonTreeContracted",
   "ScalarThresholdWilsonSoftMatchedNotRenormalized"},
  Lookup[wSaved,"status","" ]],Exit[6]];
Exit[0];
