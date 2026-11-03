(* Track B independent restricted A/D/E Matchete source-action + complete CG witness.
   Must be run with ONE model in a FRESH Wolfram kernel. No Match/RG/EFT.
   Usage: wolframscript -file scripts/RunTrackBADENativeAction.wl A_-4 OUTPUT.json
   Never overwrites any existing witness; never edits production source. *)
ClearAll["Global`*"];
$bdArgs=Take[$ScriptCommandLine,-2];
If[Length[$ScriptCommandLine]<2,Print["Usage MODEL OUTPUT.json"];Exit[2]];
$bdLabel=$bdArgs[[1]];$bdOut=ExpandFileName[$bdArgs[[2]]];
$bdCases=<|"A_-4"->{"A",-4},"A_-2"->{"A",-2},
  "D_-2"->{"D",-2},"D_+0"->{"D",0},"D_+2"->{"D",2},
  "E_-4"->{"E",-4},"E_-2"->{"E",-2},"E_+0"->{"E",0},"E_+2"->{"E",2}|>;
If[!KeyExistsQ[$bdCases,$bdLabel],Print["Unknown Track B model: ",$bdLabel];Exit[2]];
If[FileExistsQ[$bdOut],Print["BLOCKED: do not overwrite native evidence: ",$bdOut];Exit[2]];
$bdRoot=DirectoryName@DirectoryName@ExpandFileName[$InputFileName];
$bdFlags={"T3_EXPERIMENTAL_GENERAL_FERMIONS",
 "T3_EXPERIMENTAL_GENERAL_QUADRATICS","T3_EXPERIMENTAL_GENERAL_QUARTICS",
 "T3_EXPERIMENTAL_FLAVOR_MASSES"};
$bdBefore=AssociationMap[Function[f,With[{v=Environment[f]},If[StringQ[v],v,"UNSET"]]],$bdFlags];
SetEnvironment[Thread[$bdFlags->{"0","0","0","1"}]];
If[!And@@MapThread[(Environment[#1]===#2)&,{$bdFlags,{"0","0","0","1"}}],
  Print["BLOCKED: restricted/diagonal-flavour flags not imposed"];Exit[3]];
If[CheckAbort[Check[UsingFrontEnd[Needs["Matchete`"]];True,$Failed],$Aborted]=!=True,
 Print["BLOCKED: installed Matchete unavailable"];Exit[4]];
Get[FileNameJoin[{$bdRoot,"Lagrangian","model","T3ModelCatalog.wl"}]];
Get[FileNameJoin[{$bdRoot,"Lagrangian","model","LagrangianBuilder.wl"}]];

ClearAll[BDFail,BDText,BDTerms,BDSelect,BDNumber,BDCheckTerms,BDGetTensor,
  BDPayload,BDHashes,BDMain];
BDFail[s_String]:=Throw[<|"status"->"BLOCKED","reason"->s|>,"BDAbort"];
BDText[x_]:=ToString[x,InputForm];
BDTerms[x_]:=Module[{e=Expand[x]},If[Head[e]===Plus,List@@e,{e}]];
BDSelect[e_,tokens_List]:=Select[BDTerms[e],With[{s=BDText[#]},
    AllTrue[tokens,StringContainsQ[s,#]&]]&];
BDNumber[t_]:=Module[{parts=If[Head[t]===Times,List@@t,{t}]},
 Times@@Select[parts,NumberQ]];
BDCheckTerms[e_,tokens_List,expected_List]:=Module[{terms,nums},
 terms=BDSelect[e,tokens];nums=BDNumber/@terms;
 If[Sort[nums]=!=Sort[expected],
   BDFail["Native action coefficient/count changed for "<>BDText[tokens]<>
      "; got "<>BDText[nums]<>" instead of "<>BDText[expected]]];
 <|"count"->Length[terms],"prefactors"->(BDText/@nums),
   "terms_inputform"->(BDText/@terms)|>
];
BDGetTensor[name_Symbol,expectedDims_List]:=Module[{a,rr},
 rr=CheckAbort[Quiet@Check[GetCGTensor[name],$Failed],$Aborted];
 If[rr===$Failed||rr===$Aborted,BDFail["GetCGTensor failed: "<>SymbolName[name]]];
 a=If[Head[rr]===SparseArray,Normal[rr],rr];
 If[Dimensions[a]=!=expectedDims||
    !AllTrue[Flatten[a],TrueQ[NumericQ[FullSimplify[#]]]&],
   BDFail["Native tensor shape/numeric failure: "<>SymbolName[name]<>
       " expected "<>BDText[expectedDims]<>" got "<>BDText[Dimensions[a]]]];
 a
];
BDPayload[t_]:=<|"shape"->Dimensions[t],
 "real_flat"->(N[Re[#],17]&/@Flatten[t]),
 "imag_flat"->(N[Im[#],17]&/@Flatten[t])|>;
BDHashes[]:=Association@Table[
 With[{f=FileNameJoin[Prepend[FileNameSplit[p],$bdRoot]]},
  If[!FileExistsQ[f],BDFail["Missing source: "<>p]];
  p->IntegerString[FileHash[f,"SHA256"],16,64]],
 {p,{"Lagrangian/model/T3ModelCatalog.wl",
      "Lagrangian/model/T3Fields.wl",
      "Lagrangian/model/SU2Invariants.wl",
      "Lagrangian/model/LagrangianBuilder.wl",
      "Lagrangian/interactions/T3Topology.wl",
      "Lagrangian/interactions/ScalarPotential.wl"}}];
BDMain[]:=Module[{spec,model,build,lf,lint,d1,d2,df,alp,
  mass,action=<||>,tn=<||>,q,ks,name,tens,expDims,massString,
  y1o,y2o,y1h,y2h,cm,checks,kin1,kin2,repoHead},
 spec=$bdCases[$bdLabel];model=T3ModelFromClass@@spec;
 If[!AssociationQ[model],BDFail["T3 catalog failed for model"]];
 {d1,d2,df,alp}={model["Scalar1","SU2"],model["Scalar2","SU2"],
                   model["Fermion","SU2"],spec[[2]]};
 If[df=!=2||!MemberQ[{1,3},d1]||!MemberQ[{1,3},d2],
  BDFail["Unexpected native Dirac representations"]];
 build=CheckAbort[Check[UsingFrontEnd[BuildT3Lagrangian[model]],$Failed],$Aborted];
 If[!AssociationQ[build]||Lookup[build,"Status",""]=!="Success"||
    !TrueQ[Lookup[build,"FullValidation",False]],
  BDFail["Live UV Matchete CheckLagrangian not successful"]];
 lf=build["LFree"];lint=build["LInt"];
 action["scalar1_mass"]=BDCheckTerms[lf,{"Coupling[MS1,","^2","NewScalar1"},{-1}];
 action["scalar2_mass"]=BDCheckTerms[lf,{"Coupling[MS2,","^2","NewScalar2"},{-1}];
 action["fermion_mass"]=BDCheckTerms[lf,{"Coupling[MF,","NewFermion"},{-1}];
 massString=First[action["fermion_mass"]["terms_inputform"]];
 If[!StringContainsQ[massString,"Coupling[MF, {Index["]||
    !StringContainsQ[massString,"NFlavor"]||
    StringCount[massString,"NFlavor"]<3||
    !StringContainsQ[massString,"NCM[Bar["]||
    StringContainsQ[massString,"GammaCC"],
  BDFail["Native diagonal flavour-indexed DIRAC mass monomial failed"]];
 (* Yukawa: original/H.c. should be one each, with literal +1 coefficient. *)
 Do[
  name=SymbolName[y];
  action[name<>"_original"]=BDCheckTerms[lint,{"Coupling["<>name<>",","CG[T3Y"<>StringTake[name,-1]<>"CG,"},{1}];
  action[name<>"_hc"]=BDCheckTerms[lint,{"Bar[Coupling["<>name<>",","CG[Bar[T3Y"<>StringTake[name,-1]<>"CG],"},{1}];
 ,{y,{y1,y2}}];
 y1o=First[action["y1_original"]["terms_inputform"]];
 y2o=First[action["y2_original"]["terms_inputform"]];
 y1h=First[action["y1_hc"]["terms_inputform"]];
 y2h=First[action["y2_hc"]["terms_inputform"]];
 If[!StringContainsQ[y1o,"GammaCC"]||StringContainsQ[y2o,"GammaCC"]||
    !StringContainsQ[y1h,"GammaCC"]||!StringContainsQ[y2h,"GammaCC"]||
    !StringContainsQ[y1o,"Proj[1]"]||!StringContainsQ[y2o,"Proj[1]"]||
    !StringContainsQ[y1h,"Proj[-1]"]||!StringContainsQ[y2h,"Proj[-1]"],
  BDFail["Native y1/y2 GammaCC, projector, H.c. orientation failed"]];
 action["mix_original"]=BDCheckTerms[lint,{"Coupling[lambdaT3,","CG[T3MixCG,"},{1}];
 action["mix_hc"]=BDCheckTerms[lint,{"Bar[Coupling[lambdaT3,","CG[Bar[T3MixCG],"},{1}];
 (* Verify all portal invariants including their H.c.-expanded action factors. *)
 Do[
  ks=If[If[which===1,d1,d2]===1,1,2];
  Do[
   name="T3HiggsPortal"<>ToString[which]<>"CGInv"<>ToString[k];
   action["portal"<>ToString[which]<>"_Inv"<>ToString[k]]=
     BDCheckTerms[lint,{"Coupling[lambdaH"<>ToString[which]<>"Inv"<>ToString[k]<>",",
                        "CG["<>name<>","},
                      If[ks===1&&If[which===1,d1,d2]===1,{1},{1/2,1/2}]];
  ,{k,ks}];
 ,{which,{1,2}}];
 kin1=Select[BDSelect[lf,{"NewScalar1","Lorentz"}],
    !StringContainsQ[BDText[#],"Coupling[MS1,"]&];
 kin2=Select[BDSelect[lf,{"NewScalar2","Lorentz"}],
    !StringContainsQ[BDText[#],"Coupling[MS2,"]&];
 If[Length[kin1]=!=1||Length[kin2]=!=1||
    BDNumber[First[kin1]]=!=1||BDNumber[First[kin2]]=!=1,
  BDFail["Live scalar kinetic normalisation missing or not +1"]];
 (* Native tensors: dims omit SU(2) singlets; Python only inserts singleton axes. *)
 tn["T3Y1CG"]=BDPayload[BDGetTensor[T3Y1CG,DeleteCases[{2,df,d1},1]]];
 tn["T3Y2CG"]=BDPayload[BDGetTensor[T3Y2CG,DeleteCases[{2,df,d2},1]]];
 tn["T3MixCG"]=BDPayload[BDGetTensor[T3MixCG,DeleteCases[{2,2,d1,d2},1]]];
 Do[
  ks=If[If[which===1,d1,d2]===1,1,2];
  Do[
    name="T3HiggsPortal"<>ToString[which]<>"CGInv"<>ToString[k];
    expDims=DeleteCases[{2,2,If[which===1,d1,d2],If[which===1,d1,d2]},1];
    tn[name]=BDPayload[BDGetTensor[Symbol[name],expDims]];
   ,{k,ks}];
  ,{which,{1,2}}];
 checks=<|
  "restricted_general_flags_zero"->AllTrue[Take[$bdFlags,3],Environment[#]==="0"&],
  "indexed_diagonal_flavour_mass_optin_one"->(Environment[Last[$bdFlags]]==="1"),
  "full_native_Matchete_action_validated"->True,
  "exact_Dirac_MF_mass_and_flavour_indices"->True,
  "both_source_Yukawas_original_HC_chirality"->True,
  "native_source_scalar_Hessian_terms_extracted"->True,
  "full_portal_invariant_family_extracted"->True,
  "canonical_scalar_kinetic_coefficients_plus_one"->True|>;
 If[!AllTrue[Values[checks],TrueQ],BDFail["Internal native checks failed"]];
 repoHead=Quiet@Check[StringTrim[RunProcess[{"git","-C",$bdRoot,"rev-parse","HEAD"},"StandardOutput"]],"UNAVAILABLE"];
 <|"format"->"t3-track-b-ade-native-action-v1",
 "status"->"NativeRestrictedADEActionExtracted","model"->$bdLabel,
 "dimensions_alpha"->{d1,d2,df,alp},
 "hypercharges_inputform"-><|"S1"->BDText[model["Scalar1","Y"]],
   "S2"->BDText[model["Scalar2","Y"]],"F"->BDText[model["Fermion","Y"]]|>,
 "source_git_head"->repoHead,"source_sha256"->BDHashes[],
 "native_matchete_checked_lagrangian"->True,"checks"->checks,
 "native_action"-><|"fermion_mass_coefficient"->"-1",
  "flavour_mass_indexed"->True,"assembled_native_terms"->action,
  "scalar_kinetic_terms"->(BDText/@Join[kin1,kin2])|>,
 "native_tensors"->tn,
 "native_cartesian_triplet_T3_convention"->"generator matrix must be separately reviewed",
 "absolute_installed_Matchete_GammaCC_phase_measured"->False,
 "physical_source_certified"->False,"matching_RGE_EFT_executed"->False,
 "experimental_flags_prior"->$bdBefore,
 "experimental_flags_effective"->AssociationMap[Environment,$bdFlags]|>
];
$bdResult=CheckAbort[Catch[BDMain[],"BDAbort"],$Aborted];
If[!AssociationQ[$bdResult],
 Print["BLOCKED: native Mathematica abort or invalid result"];Exit[7]];
If[!DirectoryQ[DirectoryName[$bdOut]],
 CreateDirectory[DirectoryName[$bdOut],CreateIntermediateDirectories->True]];
If[Quiet@Check[Export[$bdOut,$bdResult,"RawJSON"],$Failed]===$Failed,
 Print["BLOCKED: JSON export failed"];Exit[8]];
Print[$bdLabel,": ",Lookup[$bdResult,"status","BLOCKED"]," -> ",$bdOut];
If[Lookup[$bdResult,"status",""]=!="NativeRestrictedADEActionExtracted",Exit[6]];
Exit[0];
