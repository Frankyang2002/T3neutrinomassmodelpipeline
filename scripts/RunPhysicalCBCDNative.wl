(* TRACK C ONLY. Fresh installed-Matchete action/CG/portal extraction for ONE of
   B_-3, B_+1, C_-3, C_+1. Run the SAME script four times, in fresh kernels.
   No production engine, matching, RGE or EFT; no overwrite.
   Usage: wolframscript -file scripts/RunPhysicalCBCDNative.wl B_-3 output.json
*)
ClearAll["Global`*"];
$tcArgs = Rest[$ScriptCommandLine];
If[Length[$tcArgs] != 2, Print["Usage: MODEL OUTPUT.json"]; Exit[2]];
$tcName = $tcArgs[[1]];
$tcSpecs = <|"B_-3"->{"B",-3,1},"B_+1"->{"B",1,1},
 "C_-3"->{"C",-3,3},"C_+1"->{"C",1,3}|>;
If[!KeyExistsQ[$tcSpecs,$tcName],Print["Excluded model"];Exit[2]];
$tcOut = ExpandFileName[$tcArgs[[2]]];
If[FileExistsQ[$tcOut],Print["Refusing overwrite: ",$tcOut];Exit[2]];
$tcRoot = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
$tcFlags = {"T3_EXPERIMENTAL_GENERAL_FERMIONS",
 "T3_EXPERIMENTAL_GENERAL_QUADRATICS", "T3_EXPERIMENTAL_GENERAL_QUARTICS",
 "T3_EXPERIMENTAL_FLAVOR_MASSES"};
$tcPrior = AssociationMap[Function[z,With[{v=Environment[z]},If[StringQ[v],v,"UNSET"]]],$tcFlags];
(* A native diagonal-flavour mass action is requested, not a general 3x3 mass. *)
SetEnvironment[Thread[$tcFlags -> {"0","0","0","1"}]];
If[!And@@MapThread[SameQ,{Environment/@$tcFlags,{"0","0","0","1"}}],Exit[3]];
If[CheckAbort[Check[UsingFrontEnd[Needs["Matchete`"]];True,$Failed],$Aborted]=!=True,
 Print["Matchete not available"];Exit[4]];
Get[FileNameJoin[{$tcRoot,"Lagrangian","model","T3ModelCatalog.wl"}]];
Get[FileNameJoin[{$tcRoot,"Lagrangian","model","LagrangianBuilder.wl"}]];

ClearAll[TCFail,TCStr,TCTerms,TCNum,TCSelect,TCTensor,TCArray,TCEq,TCAllEq,TCRun];
TCFail[msg_] := Throw[<|"status"->"blocked","reason"->ToString[msg,InputForm]|>,"TCStop"];
TCStr[x_] := ToString[x,InputForm];
TCTerms[e_,words_List] := Select[
 If[Head[Expand[e]]===Plus,List@@Expand[e],{Expand[e]}],
 Function[t,AllTrue[words,StringContainsQ[TCStr[t],#]&]]];
TCNum[t_] := Times@@Select[If[Head[t]===Times,List@@t,{t}],NumberQ];
TCSelect[e_,words_List,n_Integer,expect_List] := Module[{terms,coeff},
 terms=TCTerms[e,words];coeff=TCNum/@terms;
 If[Length[terms]=!=n || coeff=!=expect,
  TCFail["action coefficient/count "<>TCStr[words]<>" "<>TCStr[coeff]]];
 <|"coefficients"->(TCStr/@coeff),"terms_inputform"->(TCStr/@terms)|>];
TCTensor[name_Symbol,dims_List] := Module[{raw,arr},
 raw=CheckAbort[Quiet@Check[GetCGTensor[name],$Failed],$Aborted];
 If[raw===$Failed || raw===$Aborted,TCFail["GetCGTensor "<>SymbolName[name]]];
 arr=If[Head[raw]===SparseArray,Normal[raw],raw];
 If[Dimensions[arr]=!=dims || !And@@(NumericQ/@(N/@Flatten[arr])),
  TCFail["invalid tensor "<>SymbolName[name]]];arr];
TCArray[a_] := <|"shape"->Dimensions[a],"real_flat"->(N[Re[#],16]&/@Flatten[a]),
 "imag_flat"->(N[Im[#],16]&/@Flatten[a])|>;
TCEq[a_,b_] := TrueQ[FullSimplify[a-b==0]];
TCAllEq[a_,b_] := And@@MapThread[TCEq,{Flatten[a],Flatten[b]}];

TCRun[] := Module[{spec,model,built,lf,li,df,mass,y1,y2,p11,p12,p21,p22,
 mix,portal,kin,sourcepaths,checks,identity4,head,dig,c1,c2,mm},
 spec=$tcSpecs[$tcName];model=T3ModelFromClass[spec[[1]],spec[[2]]];
 If[!AssociationQ[model] || model["Fermion","SU2"]=!=spec[[3]] ||
   model["Fermion","Y"]===0,TCFail["catalog/charged Dirac contract"]];
 built=CheckAbort[Check[UsingFrontEnd[BuildT3Lagrangian[model]],$Failed],$Aborted];
 If[!AssociationQ[built] || Lookup[built,"Status",""]=!="Success" ||
    !TrueQ[Lookup[built,"FullValidation",False]],
   TCFail["CheckLagrangian restricted action failed"]];
 lf=built["LFree"];li=built["LInt"];
 mass=TCSelect[lf,{"Coupling[MF,","NewFermion","NCM[Bar["},1,{-1}];
 If[!StringContainsQ[First[mass["terms_inputform"]],"NFlavor"] ||
    !StringContainsQ[First[mass["terms_inputform"]],"Bar[Field[NewFermion"],
  TCFail["Dirac flavour-indexed mass must pair Bar[F] F"]];
 If[StringContainsQ[First[mass["terms_inputform"]],"GammaCC"],
  TCFail["Majorana mass structure in charged Dirac model"]];
 kin=Join[TCTerms[lf,{"NewScalar1","Lorentz"}],TCTerms[lf,{"NewScalar2","Lorentz"}]];
 If[Length[kin]=!=2 || !AllTrue[kin,TCNum[#]===1&],TCFail["scalar kinetic normalisation"]];
 portal=<|
  "y1_original"->TCSelect[li,{"Coupling[y1,","CG[T3Y1CG,","DiracProduct[GammaCC"},1,{1}],
  "y2_original"->TCSelect[li,{"Coupling[y2,","CG[T3Y2CG,","DiracProduct[Proj[1]]"},1,{1}],
  "y1_hc"->TCSelect[li,{"Bar[Coupling[y1,","CG[Bar[T3Y1CG],"},1,{1}],
  "y2_hc"->TCSelect[li,{"Bar[Coupling[y2,","CG[Bar[T3Y2CG],"},1,{1}],
  "mixing"->TCSelect[li,{"Coupling[lambdaT3,","CG[T3MixCG,"},1,{1}],
  "mixing_hc"->TCSelect[li,{"Bar[Coupling[lambdaT3,","CG[Bar[T3MixCG],"},1,{1}],
  "S1Inv1"->TCSelect[li,{"Coupling[lambdaH1Inv1,","CG[T3HiggsPortal1CGInv1,"},1,{1}],
  "S1Inv2"->TCSelect[li,{"Coupling[lambdaH1Inv2,","CG[T3HiggsPortal1CGInv2,"},1,{1}],
  "S2Inv1"->TCSelect[li,{"Coupling[lambdaH2Inv1,","CG[T3HiggsPortal2CGInv1,"},1,{1}],
  "S2Inv2"->TCSelect[li,{"Coupling[lambdaH2Inv2,","CG[T3HiggsPortal2CGInv2,"},1,{1}]|>;
 y1=TCTensor[T3Y1CG,If[spec[[3]]===1,{2,2},{2,3,2}]];
 y2=TCTensor[T3Y2CG,If[spec[[3]]===1,{2,2},{2,3,2}]];
 mix=TCTensor[T3MixCG,{2,2,2,2}];
 p11=TCTensor[T3HiggsPortal1CGInv1,{2,2,2,2}];
 p12=TCTensor[T3HiggsPortal1CGInv2,{2,2,2,2}];
 p21=TCTensor[T3HiggsPortal2CGInv1,{2,2,2,2}];
 p22=TCTensor[T3HiggsPortal2CGInv2,{2,2,2,2}];
 identity4=Table[KroneckerDelta[i,j] KroneckerDelta[a,b],
 {i,1,2},{j,1,2},{a,1,2},{b,1,2}];
 c1=If[spec[[3]]===1,y1[[1,1]],y1[[1,2,1]]];
 c2=If[spec[[3]]===1,y2[[1,2]],y2[[1,2,2]]];
 mm=mix[[2,2,All,All]];
 checks=<|
  "charged_dirac_field"->(!TrueQ[model["Fermion","Y"]===0]),
  "source_yukawa_orientations_and_Hc"->True,
  "dirac_minus_one_not_majorana_half"->True,
  "mass_one_flavour_index"->True,
  "CG_y1_charge_channel"->TCEq[c1,If[spec[[3]]===1,1,3^(-1/4)]],
  "CG_y2_charge_channel"->TCEq[c2,If[spec[[3]]===1,1,-3^(-1/4)]],
  "H0_S1_1_S2_2_mix"->TCEq[mm[[1,2]],-2/Sqrt[3]],
  "H0_no_other_mixing"->(Count[Flatten[mm],u_/;!TrueQ[PossibleZeroQ[u]]]===1),
  "S1_full_offshell_portal_identity"->TCAllEq[(p11+Sqrt[3] p12)/2,identity4],
  "S2_full_offshell_portal_identity"->TCAllEq[(p21+Sqrt[3] p22)/2,identity4],
  "S1_H0_Inv1_projector"->TCAllEq[p11[[2,2,All,All]],DiagonalMatrix[{0,1}]],
  "S2_H0_Inv1_projector"->TCAllEq[p21[[2,2,All,All]],DiagonalMatrix[{0,1}]],
  "S1_H0_Inv2_projector"->TCAllEq[p12[[2,2,All,All]],DiagonalMatrix[{2,1}]/Sqrt[3]],
  "S2_H0_Inv2_projector"->TCAllEq[p22[[2,2,All,All]],DiagonalMatrix[{2,1}]/Sqrt[3]]|>;
 If[!AllTrue[Values[checks],TrueQ],TCFail["failed native component checks "<>TCStr[checks]]];
 sourcepaths={"Lagrangian/model/T3ModelCatalog.wl",
  "Lagrangian/model/LagrangianBuilder.wl","Lagrangian/model/T3Fields.wl",
  "Lagrangian/model/SU2Invariants.wl","Lagrangian/interactions/T3Topology.wl",
  "Lagrangian/interactions/ScalarPotential.wl"};
 dig=Association@Table[With[{f=FileNameJoin[Prepend[FileNameSplit[z],$tcRoot]]},
  If[!FileExistsQ[f],TCFail["source missing "<>z]];
  z->IntegerString[FileHash[f,"SHA256"],16,64]],{z,sourcepaths}];
 head=StringTrim[Quiet@Check[RunProcess[{"git","-C",$tcRoot,"rev-parse","HEAD"},
  "StandardOutput"],"UNAVAILABLE"]];
 <|"format"->"track-c-bc-native-action-v1",
 "status"->"native_action_and_tensors_extracted_spinor_absolute_pending",
 "model"->$tcName,"local_git_head"->head,"source_sha256"->dig,
 "restricted_flags_effective"->AssociationMap[Environment,$tcFlags],
 "full_CheckLagrangian_valid"->True,"mass"->mass,
 "action"->portal,"kinetic_terms"->(TCStr/@kin),
 "native_tensors"-><|"T3Y1CG"->TCArray[y1],"T3Y2CG"->TCArray[y2],
  "T3MixCG"->TCArray[mix],"T3HiggsPortal1CGInv1"->TCArray[p11],
  "T3HiggsPortal1CGInv2"->TCArray[p12],
  "T3HiggsPortal2CGInv1"->TCArray[p21],
  "T3HiggsPortal2CGInv2"->TCArray[p22]|>,
 "native_checks"->checks,
 "native_numeric_MF_eigenvalues_assigned"->False,
 "real_triplet_cartesian_dirac_mass_metric"->"Bar[F_a] delta_ab F_b (not F^T G F)",
 "benchmark_identity_portal_is_conditional_choice"->True,
 "absolute_installed_GammaCC_phase_verified"->False,
 "physical_certification"->False|>
];
$tcResult=CheckAbort[Catch[TCRun[],"TCStop"],$Aborted];
If[!AssociationQ[$tcResult] || Lookup[$tcResult,"status",""]=!=
 "native_action_and_tensors_extracted_spinor_absolute_pending",
 Print["BLOCKED: ",If[AssociationQ[$tcResult],Lookup[$tcResult,"reason","unknown"],"abort"]];Exit[5]];
If[!DirectoryQ[DirectoryName[$tcOut]],
 CreateDirectory[DirectoryName[$tcOut],CreateIntermediateDirectories->True]];
If[Quiet@Check[Export[$tcOut,$tcResult,"RawJSON"],$Failed]===$Failed,Exit[6]];
Print["TRACK C native action/CG/portal extraction PASS: ",$tcName," -> ",$tcOut];
Exit[0];
