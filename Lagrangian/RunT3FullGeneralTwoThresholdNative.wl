(* Local consolidated-pipeline repair for the verified 2026-11-10 logs.
   Matches three F generations simultaneously at muF, then both scalars at muS.
   In A(-4), the ordinary Matchete second threshold cannot accept three
   repeated NFlavor indices in y1[p,r] y2[q,r] / MF[r].  Reparameterise the
   EXACT tree-level LLSS source into independent W[p,q] before stage 2.

   The one-loop F hard source remains EXPORTED but is NOT passed to scalar
   matching: it has not yet been brought to the independent EFT1 Wilson basis.
   Scalar-soft stage therefore CANNOT be reported as complete fixed-order C5.

   In D(-2), NEVER call the unrotated generic first threshold.  The saved
   first threshold already used split neutral real scalar eigenfields. The
   missing continuation must preserve that mass basis at the scalar stage.

   CLI: wolframscript -file thisfile OUT_DIR A_-4|D_-2 [REUSE_F_INPUT_DIR]
   REUSE_F_INPUT_DIR is the existing A(-4) debug/fresh_kernel_stage_2
   directory with input_tree.wxf and input_loop.wxf.  It avoids rematching F.
*)
ClearAll["Global`*"];
t3args=Rest[$ScriptCommandLine];
If[!MemberQ[{2,3},Length[t3args]] || !MemberQ[{"A_-4","D_-2"},t3args[[2]]],
 Print["Usage: OUT_DIR A_-4|D_-2 [REUSE_F_INPUT_DIR]"];Exit[2]];
t3out=ExpandFileName[t3args[[1]]];
t3modelName=t3args[[2]];
t3reuse=If[Length[t3args]===3,ExpandFileName[t3args[[3]]],None];
t3root=DirectoryName@ExpandFileName[$InputFileName];
t3project=DirectoryName[t3root];
If[!DirectoryQ[t3out],CreateDirectory[t3out,CreateIntermediateDirectories->True]];
t3status=FileNameJoin[{t3out,"native_two_threshold_status.json"}];
t3stage="init";
ClearAll[t3Write,t3Stop,t3Save];
t3Write[status_String,reason_String:"",more_:<||>]:=Module[{rec},
 rec=Join[<|
   "format"->"t3-full-general-two-threshold-native-repair-v2",
   "model"->t3modelName,"status"->status,"stage"->t3stage,
   "reason"->reason,"heavy_fermion_flavours"->3,
   "fermion_groups"->{{"F"}},"scalar_groups"->{{"S1","S2"}},
   "group_count"->2,"one_loop_matching_requested"->True,
   "eft1_wilson_rge_inserted"->False,
   "renormalized_finite_wilsons_verified"->False,
   "complete_hard_plus_soft_matching_performed"->False,
   "full_general_numerical_c5_verified"->False|>,more];
 Export[t3status,rec,"RawJSON"];
 Print[t3modelName," ",status," at ",t3stage,
  If[reason==="","",": "<>reason]];rec];
t3Stop[msg_String,more_:<||>]:=(t3Write["Blocked_NativeStage",msg,more];Exit[5]);
t3Save[name_String,obj_,kind_String:"WXF"]:=Module[{p=FileNameJoin[{t3out,name}]},
 If[Quiet@Check[Export[p,obj,kind],$Failed]===$Failed||!FileExistsQ[p],
  t3Stop["Could not export "<>name]];p];
If[!And@@(Environment[#]==="1"&/@{
 "T3_EXPERIMENTAL_GENERAL_FERMIONS","T3_EXPERIMENTAL_GENERAL_QUARTICS",
 "T3_EXPERIMENTAL_GENERAL_QUADRATICS","T3_EXPERIMENTAL_FLAVOR_MASSES"}),
 t3Stop["All four native full-general switches must be set to 1"]];

(* D(-2) physical F-first is already present.  No new F matching is needed. *)
If[t3modelName==="D_-2",
 t3stage="D-real-scalar-F-source-reuse";
 t3source=FileNameJoin[{t3project,"output","parallel_t3",
   "full_general_batch2","thesis_full_general_batch2_raw_v1","D_-2",
   "native_matching"}];
 t3native=FileNameJoin[{t3source,"native_status.json"}];
 t3basis=FileNameJoin[{t3source,"self_mass_basis.json"}];
 If[!And@@(FileExistsQ/@{t3native,t3basis}),
  t3Stop["Successful physical D(-2) F source/basis metadata not found"]];
 t3meta=Quiet@Check[Import[t3native,"RawJSON"],$Failed];
 t3mass=Quiet@Check[Import[t3basis,"RawJSON"],$Failed];
 If[!AssociationQ[t3meta]||!AssociationQ[t3mass]||
   Lookup[t3meta,"status",""]=!="FThresholdTreeWilsonExported"||
   Lookup[t3meta,"heavy_flavour_count",0]=!=3||
   Lookup[t3mass,"mass_basis_kind",""]=!="RealSelfBilinearSplitLightFields"||
   Lookup[t3mass,"split_species",0]=!=2,
   t3Stop["Physical split-scalar source is not validated"]];
 t3Write["Blocked_D2NeedsPhysicalScalarContinuation",
   "D(-2) successful F matching uses T3Scalar2Minus/T3Scalar2Plus; "<>
   "the generic NewScalar2 second-threshold kernel is incompatible. "<>
   "Do not repeat the unrotated F match or substitute an NF=1 scalar fit.",
   <|"verified_F_source"->t3native,
     "physical_scalar_fields_at_second_threshold"->
       {"NewScalar1","T3Scalar2Minus","T3Scalar2Plus"},
     "real_bilinear_mass_basis_reused"->True|>];
 Exit[5]
];

(* A(-4) has no self bilinear; its F-threshold physical basis is unchanged. *)
t3stage="A-matchete-load";
If[CheckAbort[UsingFrontEnd[Needs["Matchete`"];True],$Aborted]=!=True,
 t3Stop["Matchete could not be loaded"]];
Get[FileNameJoin[{t3root,"model","T3ModelCatalog.wl"}]];
Get[FileNameJoin[{t3root,"model","LagrangianBuilder.wl"}]];
Get[FileNameJoin[{t3root,"RunMatching.wl"}]];
Get[FileNameJoin[{t3root,"RunT3EFT1WilsonStageLauncher.wl"}]];
t3model=T3ModelFromClass["A",-4,3];
If[!AssociationQ[t3model],t3Stop["A(-4) three-flavour catalog failed"]];
If[SetT3HeavyFields[{"F"}]===$Failed,t3Stop["Common F threshold declaration failed"]];
t3stage="A-F-boundary-reuse-or-match";
If[StringQ[t3reuse],
 t3treePath=FileNameJoin[{t3reuse,"input_tree.wxf"}];
 t3loopPath=FileNameJoin[{t3reuse,"input_loop.wxf"}];
 If[!And@@(FileExistsQ/@{t3treePath,t3loopPath}),
  t3Stop["Requested F result reuse lacks tree/loop WXF"]];
 t3tree=Quiet@Check[Import[t3treePath,"WXF"],$Failed];
 t3hard=Quiet@Check[Import[t3loopPath,"WXF"],$Failed];
 If[MemberQ[{t3tree,t3hard},$Failed]||FreeQ[t3tree,GammaCC]||
    FreeQ[t3hard,hbar],t3Stop["Reused F WXF lacks tree LLSS or one-loop hard"]];
 t3reused=True,
 t3build=CheckAbort[UsingFrontEnd[BuildT3Lagrangian[t3model]],$Aborted];
 If[!AssociationQ[t3build]||Lookup[t3build,"Status",""]=!="Success"||
    !TrueQ[Lookup[t3build,"FullValidation",False]]||
    !TrueQ[Lookup[t3build,"AllScalarExtensionsRequested",False]]||
    Lookup[t3build,"AdditionalQuarticsIncluded",-1]=!=0,
    t3Stop["Full general A(-4) UV builder failed"]];
 (* F-only matching, all 3 flavours together; not three separate thresholds. *)
 t3f=CheckAbort[UsingFrontEnd[RunSequentialT3Matching[
    t3build["LUV"],t3model,{{"F"}},5,1,t3out,t3modelName]],$Aborted];
 If[!AssociationQ[t3f]||Lookup[t3f,"Status",""]=!="Success"||
    Length[Lookup[t3f,"Stages",{}]]=!=1,
   t3Stop["F-only native matching failed: "<>
    If[AssociationQ[t3f],ToString[Lookup[t3f,"Status","Unknown"]],
      ToString[InputForm[Head[t3f]]]]]];
 t3tree=First[t3f["Stages"]]["TreeLagrangian"];
 t3hard=First[t3f["Stages"]]["OneLoopCorrection"];
 t3reused=False
];
t3Save["A_full_F_tree_original.wxf",t3tree];
t3Save["A_full_F_hard_one_loop_preserved.wxf",t3hard];

(* Exact algebraic reparameterisation at tree level:
     W[p,q] = Sum_r y1[p,r] y2[q,r] / MF[r] .
   This eliminates the three-fold repeated r in each scalar matching term.
   Its numerical definition is EXPORTED: the 3 F flavours are not discarded.
   It does NOT resolve the independent one-loop hard Wilson operators. *)
t3stage="A-tree-flavour-Wilson-contract";
If[Quiet@Check[
 DefineCoupling[T3EFT1W12,Indices->{Flavor,Flavor},SelfConjugate->False];
 DefineCoupling[T3EFT1WExtra,Indices->{Flavor,Flavor},SelfConjugate->False];
 True,False]=!=True,t3Stop["Could not register flavour Wilson tensors"]];
ClearAll[t3InverseIndex,t3YukawaFactor,t3ContractTerm];
t3InverseIndex[x_]:=Replace[Unevaluated[x],{
 HoldPattern[Power[Coupling[MF,{r:Index[_,NFlavor]},0],-1]]:>r,
 _:>Missing["NotIndexedInverseMass"]}];
t3YukawaFactor[x_]:=Replace[Unevaluated[x],{
 HoldPattern[Coupling[n:(y1|y2),{p:Index[_,Flavor],r:Index[_,NFlavor]},0]]:>
   <|"name"->n,"p"->p,"r"->r,"barred"->False|>,
 HoldPattern[Bar[Coupling[n:(y1|y2),{p:Index[_,Flavor],r:Index[_,NFlavor]},0]]]:>
   <|"name"->n,"p"->p,"r"->r,"barred"->True|>,
 _:>Missing["NotYukawa"]}];
t3WCount=0;
t3ContractTerm[term_]:=Module[
 {factors,inds,massPos,heavy,ys,y1pos,y2pos,first,second,remaining,w},
 factors=If[Head[term]===Times,List@@term,{term}];
 inds=t3InverseIndex/@factors;
 massPos=Select[Range[Length[inds]],!MissingQ[inds[[#]]]&];
 If[massPos==={},Return[term]];
 If[Length[massPos]=!=1,t3Stop["More than one inverse MF in tree term"]];
 heavy=inds[[First[massPos]]];
 ys=t3YukawaFactor/@factors;
 y1pos=Select[Range[Length[ys]],AssociationQ[ys[[#]]]&&
    Lookup[ys[[#]],"name",None]===y1&&Lookup[ys[[#]],"r",None]===heavy&];
 y2pos=Select[Range[Length[ys]],AssociationQ[ys[[#]]]&&
    Lookup[ys[[#]],"name",None]===y2&&Lookup[ys[[#]],"r",None]===heavy&];
 If[Length[y1pos]=!=1||Length[y2pos]=!=1,
  t3Stop["Cannot map inverse MF to unique y1 and y2"]];
 first=ys[[First[y1pos]]];second=ys[[First[y2pos]]];
 If[first["barred"]=!=second["barred"],
  t3Stop["Inconsistent conjugation of heavy-flavour source"]];
 w=Coupling[T3EFT1W12,{first["p"],second["p"]},0];
 If[TrueQ[first["barred"]],w=Bar[w]];
 remaining=Delete[factors,List/@Sort[{
   First[massPos],First[y1pos],First[y2pos]}]];
 t3WCount++;
 Times@@Append[remaining,w]
];
t3terms=If[Head[t3tree]===Plus,List@@t3tree,{t3tree}];
t3treeW=Plus@@(t3ContractTerm/@t3terms);
If[t3WCount=!=2||FreeQ[t3treeW,T3EFT1W12]||
  !FreeQ[t3treeW,Power[Coupling[MF,{Index[_,NFlavor]},0],-1]],
 t3Stop["Wrong A tree flavour-Wilson inventory; expected 2 terms",
  <|"converted_W_terms"->t3WCount|>]];
t3Save["A_EFT1_tree_contracted.wxf",t3treeW];
t3Save["A_EFT1_tree_contracted.txt",ToString[t3treeW,InputForm],"Text"];
t3Save["A_EFT1_Wilson_definitions.json",<|
  "W12[p,q]"->"Sum[y1[p,r]*y2[q,r]/MF[r],{r,1,3}]",
  "conjugate_W12"->"Conjugate[W12[p,q]] for positive real MF[r]",
  "heavy_fermions_integrated_together"->3,
  "physical_RGBeta_normalization_certified"->False,
  "one_loop_hard_reparameterised"->False|>,"RawJSON"];

(* Existing experimental scalar matcher registers the W couplings in the
   fresh kernel.  Only F TREE operators enter; hard 1L remains preserved.
   This is a scalar-soft matching diagnostic, NOT an entire C5. *)
t3stage="A-joint-scalar-soft-one-loop";
t3scalar=CheckAbort[RunT3EFT1WilsonFreshThresholdStage[
 t3treeW,0,t3treeW,t3treeW,t3model,{"S1","S2"},{"S1","S2"},
 2,5,1,t3out],$Aborted];
If[!AssociationQ[t3scalar]||Lookup[t3scalar,"Status",""]=!="Success",
 t3Stop["Joint scalar-soft native matching failed; see "<>
  "debug/fresh_kernel_stage_2/stdout.log",
  <|"F_tree_flavour_Wilsons_contracted"->True,
    "first_threshold_reused"->t3reused,
    "native_scalar_result_head"->ToString[InputForm[Head[t3scalar]]]|>]];
t3soft=Lookup[t3scalar,"AuthoritativeThresholdC5",$Failed];
If[t3soft===$Failed,
 t3Stop["Joint scalar match succeeded but did not return a Weinberg coefficient"]];
t3Save["A_scalar_soft_only_unrenormalized.wxf",t3soft];
t3Save["A_scalar_soft_only_unrenormalized.txt",
 ToString[t3soft,InputForm],"Text"];
t3stage="A-soft-success-hard-inheritance-pending";
t3Write["ScalarSoftOnly_HardInheritanceAndWilsonRGPending",
 "Joint scalar-soft matching succeeded, but F hard one-loop operators "<>
 "and EFT1 Wilson RG were not included. NOT a physical full-general C5.",
 <|"A_scalar_threshold_native_matched"->True,
   "first_threshold_reused"->t3reused,
   "three_flavour_Wilson_definitions_preserved"->True,
   "A_hard_source_WXF"->"A_full_F_hard_one_loop_preserved.wxf",
   "A_soft_source_WXF"->"A_scalar_soft_only_unrenormalized.wxf"|>];
Exit[5];
