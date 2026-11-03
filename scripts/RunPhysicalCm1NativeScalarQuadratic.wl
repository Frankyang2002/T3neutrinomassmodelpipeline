(* C_-1 Majorana-triplet: live native Matchete scalar/Yukawa/mass witness.
   The native Inv1 is an exchange tensor. The identity is exactly
   1/2 Inv1 + Sqrt[3]/2 Inv2 for BOTH scalar doublets (all 4 axes).

   DO NOT overwrite prior failed evidence. Use a fresh output directory.

   Source: restricted assembled LUV and GetCGTensor, never provisional tensors.
   No Match[], EFT, RGBeta, RGE or GitHub changes.
   Usage: wolframscript -file scripts/RunPhysicalCm1NativeScalarQuadratic.wl OUTPUT.json
*)
ClearAll["Global`*"];
$bmArgs = Rest[$ScriptCommandLine];
If[Length[$bmArgs] =!= 1, Print["Usage: OUTPUT.json"]; Exit[2]];
$bmOut = ExpandFileName[First[$bmArgs]];
If[FileExistsQ[$bmOut], Print["BLOCKED: witness already exists: ", $bmOut]; Exit[2]];
$bmRoot = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
$bmFlags = {"T3_EXPERIMENTAL_GENERAL_FERMIONS", "T3_EXPERIMENTAL_GENERAL_QUADRATICS",
  "T3_EXPERIMENTAL_GENERAL_QUARTICS", "T3_EXPERIMENTAL_FLAVOR_MASSES"};
$bmPriorFlags = AssociationMap[Function[x, With[{v=Environment[x]},If[StringQ[v],v,"UNSET"]]],$bmFlags];
SetEnvironment[Thread[$bmFlags -> ConstantArray["0",Length[$bmFlags]]]];
If[!AllTrue[$bmFlags, Environment[#] === "0" &],
 Print["BLOCKED: restricted flags not isolated"]; Exit[3]];
If[CheckAbort[Check[UsingFrontEnd[Needs["Matchete`"]];True,$Failed],$Aborted] =!= True,
 Print["BLOCKED: Matchete unavailable"]; Exit[4]];
Get[FileNameJoin[{$bmRoot,"Lagrangian","model","T3ModelCatalog.wl"}]];
Get[FileNameJoin[{$bmRoot,"Lagrangian","model","LagrangianBuilder.wl"}]];

ClearAll[BMFail,BMStr,BMTerms,BMNumber,BMSelect,BMTensor,BMArray,BMEqual,BMChecks,BMRun];
BMFail[msg_String] := Throw[<|"status"->"blocked","reason"->msg|>,"BMStop"];
BMStr[x_] := ToString[x,InputForm];
BMTerms[e_, tokens_List] := Select[
 If[Head[Expand[e]]===Plus,List@@Expand[e],{Expand[e]}],
 Function[t,AllTrue[tokens,StringContainsQ[BMStr[t],#]&]]];
BMNumber[t_] := Times @@ Select[If[Head[t]===Times,List@@t,{t}],NumberQ];
BMSelect[e_, tokens_List, count_Integer, expected_List] := Module[{found,coeff},
 found=BMTerms[e,tokens]; coeff=BMNumber/@found;
 If[Length[found]=!=count || coeff=!=expected,
  BMFail["unexpected native term coefficients for "<>BMStr[tokens]<>
   ": count="<>ToString[Length[found]]<>" coefficients="<>BMStr[coeff]]];
 <|"count"->count,"coefficients"->(BMStr/@coeff),"terms_inputform"->(BMStr/@found)|>
];
BMTensor[name_Symbol,dims_List] := Module[{raw,arr},
 raw=CheckAbort[Quiet@Check[GetCGTensor[name],$Failed],$Aborted];
 If[raw===$Failed || raw===$Aborted,BMFail["GetCGTensor failed: "<>SymbolName[name]]];
 arr=If[Head[raw]===SparseArray,Normal[raw],raw];
 If[Dimensions[arr]=!=dims || !AllTrue[Flatten[arr],TrueQ[NumericQ[FullSimplify[#]]]&],
  BMFail["invalid registered tensor: "<>SymbolName[name]]];
 arr
];
BMArray[arr_] := <|"shape"->Dimensions[arr],
 "real_flat"->(N[Re[#],16]&/@Flatten[arr]),
 "imag_flat"->(N[Im[#],16]&/@Flatten[arr])|>;
BMEqual[x_,y_] := TrueQ[FullSimplify[x-y==0]];
BMChecks[x_,y_] := And@@Flatten[MapThread[BMEqual,{x,y},ArrayDepth[x]]];
BMRun[] := Module[{model,built,lint,lf,cs1,cs2,cm,p11,p12,p21,p22,mass,
  kin1,kin2,portal,mu,proj11,proj12,proj21,proj22,checks,digest,paths,githead,
  unit4,identity1,identity2},
 model=T3ModelFromClass["C",-1];
 If[!AssociationQ[model],BMFail["C_-1 catalog unavailable"]];
 built=CheckAbort[Check[UsingFrontEnd[BuildT3Lagrangian[model]],$Failed],$Aborted];
 If[!AssociationQ[built] || Lookup[built,"Status",""]=!="Success" ||
    !TrueQ[Lookup[built,"FullValidation",False]],
  BMFail["native restricted C_-1 action not CheckLagrangian-valid"]];
 lint=built["LInt"]; lf=built["LFree"];
 mass=BMSelect[lf,{"Coupling[MF,","GammaCC","NewFermion"},1,{-1/2}];
 kin1=BMTerms[lf,{"NewScalar1","Lorentz"}];
 kin2=BMTerms[lf,{"NewScalar2","Lorentz"}];
 If[Length[kin1]=!=1 || Length[kin2]=!=1 ||
    BMNumber[First[kin1]]=!=1 || BMNumber[First[kin2]]=!=1,
  BMFail["canonical scalar kinetic terms absent"]];
 portal=<|
  "S1Inv1"->BMSelect[lint,{"Coupling[lambdaH1Inv1,","CG[T3HiggsPortal1CGInv1,"},1,{1}],
  "S1Inv2"->BMSelect[lint,{"Coupling[lambdaH1Inv2,","CG[T3HiggsPortal1CGInv2,"},1,{1}],
  "S2Inv1"->BMSelect[lint,{"Coupling[lambdaH2Inv1,","CG[T3HiggsPortal2CGInv1,"},1,{1}],
  "S2Inv2"->BMSelect[lint,{"Coupling[lambdaH2Inv2,","CG[T3HiggsPortal2CGInv2,"},1,{1}],
  "mixing"->BMSelect[lint,{"Coupling[lambdaT3,","CG[T3MixCG,"},1,{1}],
  "mixing_hc"->BMSelect[lint,{"Bar[Coupling[lambdaT3,","CG[Bar[T3MixCG],"},1,{1}],
  "y1"->BMSelect[lint,{"Coupling[y1,","CG[T3Y1CG,","DiracProduct[Proj[1]]"},1,{1}],
  "y2"->BMSelect[lint,{"Coupling[y2,","CG[T3Y2CG,","DiracProduct[Proj[1]]"},1,{1}],
  "y1_hc"->BMSelect[lint,{"Bar[Coupling[y1,","CG[Bar[T3Y1CG],"},1,{1}],
  "y2_hc"->BMSelect[lint,{"Bar[Coupling[y2,","CG[Bar[T3Y2CG],"},1,{1}]
 |>; 
 cs1=BMTensor[T3Y1CG,{2,3,2}];cs2=BMTensor[T3Y2CG,{2,3,2}];
 cm=BMTensor[T3MixCG,{2,2,2,2}];
 p11=BMTensor[T3HiggsPortal1CGInv1,{2,2,2,2}];
 p12=BMTensor[T3HiggsPortal1CGInv2,{2,2,2,2}];
 p21=BMTensor[T3HiggsPortal2CGInv1,{2,2,2,2}];
 p22=BMTensor[T3HiggsPortal2CGInv2,{2,2,2,2}];
 proj11=(p11[[2,2,All,All]]+Conjugate[Transpose[p11[[2,2,All,All]]]])/2;
 proj12=(p12[[2,2,All,All]]+Conjugate[Transpose[p12[[2,2,All,All]]]])/2;
 proj21=(p21[[2,2,All,All]]+Conjugate[Transpose[p21[[2,2,All,All]]]])/2;
 proj22=(p22[[2,2,All,All]]+Conjugate[Transpose[p22[[2,2,All,All]]]])/2;
 mu=cm[[2,2,All,All]];
 (* For each complex doublet S, the SOURCE four-index operator basis is
    exchange Inv1 and a mixed Inv2. Do not mistake Inv1's H0 projection,
    DiagonalMatrix[{0,1}], for the identity. The following verifies the
    identity at ALL Higgs and scalar indices, not just after H0 projection.
    Source action has +1/2 PlusHc for real portal coefficients; once both
    terms are assembled, their quadratic coefficients are +1 each. *)
 unit4=Table[KroneckerDelta[i,j] KroneckerDelta[a,b],
   {i,1,2},{j,1,2},{a,1,2},{b,1,2}];
 identity1=(p11+Sqrt[3] p12)/2;
 identity2=(p21+Sqrt[3] p22)/2;
 checks=<|
  "C_minus1_dimensions"->(model["Scalar1","SU2"]===2 &&
    model["Scalar2","SU2"]===2 && model["Fermion","SU2"]===3),
  "majorana_source_mass_minus_half"->True,
  "triplet_majorana_mass_indices"->StringContainsQ[First[mass["terms_inputform"]],"T3BSMd3"],
  "neutral_y1_triplet_component"->BMEqual[cs1[[1,2,1]],3^(-1/4)],
  "neutral_y2_triplet_component"->BMEqual[cs2[[1,2,2]],-3^(-1/4)],
  "neutral_F_cartesian_axis_unique"->(
     Count[cs1[[1,All,1]],x_/;!TrueQ[PossibleZeroQ[x]]]===1 &&
     Count[cs2[[1,All,2]],x_/;!TrueQ[PossibleZeroQ[x]]]===1),
  "native_mix_H0_component"->BMEqual[mu[[1,2]],-2/Sqrt[3]],
  "native_mix_H0_single_charge_channel"->(Count[Flatten[mu],x_/;!TrueQ[PossibleZeroQ[x]]]===1),
  "portal1_exchange_H0"->BMChecks[proj11,DiagonalMatrix[{0,1}]],
  "portal2_exchange_H0"->BMChecks[proj21,DiagonalMatrix[{0,1}]],
  "portal1_second_H0"->BMChecks[proj12,DiagonalMatrix[{2/Sqrt[3],1/Sqrt[3]}]],
  "portal2_second_H0"->BMChecks[proj22,DiagonalMatrix[{2/Sqrt[3],1/Sqrt[3]}]],
  "portal1_full_four_index_identity"->BMChecks[identity1,unit4],
  "portal2_full_four_index_identity"->BMChecks[identity2,unit4],
  "portal1_vev_identity"->BMChecks[(proj11+Sqrt[3]proj12)/2,IdentityMatrix[2]],
  "portal2_vev_identity"->BMChecks[(proj21+Sqrt[3]proj22)/2,IdentityMatrix[2]]
 |>; 
 If[!AllTrue[Values[checks],TrueQ],
  BMFail["native component witness does not satisfy C_-1 charge/mass contract: "<>BMStr[checks]]];
 paths={"Lagrangian/model/T3ModelCatalog.wl","Lagrangian/model/LagrangianBuilder.wl",
   "Lagrangian/model/SU2Invariants.wl","Lagrangian/model/T3Fields.wl",
   "Lagrangian/interactions/ScalarPotential.wl","Lagrangian/interactions/T3Topology.wl"};
 digest=Association@Table[With[{f=FileNameJoin[Prepend[FileNameSplit[path],$bmRoot]]},
  If[!FileExistsQ[f],BMFail["missing source "<>path]];
  path->IntegerString[FileHash[f,"SHA256"],16,64]],{path,paths}];
 githead=Quiet@Check[RunProcess[{"git","-C",$bmRoot,"rev-parse","HEAD"},"StandardOutput"],"UNAVAILABLE"];
 <|
  "format"->"t3-b1-cm1-live-native-scalar-quadratic-v1",
  "status"->"NativeCm1TripletScalarAndMajoranaActionExtracted_ExactPortalIdentityMap_PhysicalDictionaryPending",
  "model"->"C_-1","source_git_head"->StringTrim[githead],
  "source_sha256"->digest,"native_lagrangian_validated"->True,
  "restricted_environment"-><|"prior"->$bmPriorFlags,
     "effective"->AssociationMap[Environment,$bmFlags]|>,
  "source_hypercharges"-><|"S1"->"-1/2","S2"->"+1/2","F"->"0"|>,
  "native_mass_term"->mass,"native_action_terms"->portal,
  "scalar_free_kinetic"->(BMStr/@Join[kin1,kin2]),
  "native_tensors"-><|"T3Y1CG"->BMArray[cs1],"T3Y2CG"->BMArray[cs2],
   "T3MixCG"->BMArray[cm],"T3HiggsPortal1CGInv1"->BMArray[p11],
   "T3HiggsPortal1CGInv2"->BMArray[p12],
   "T3HiggsPortal2CGInv1"->BMArray[p21],
   "T3HiggsPortal2CGInv2"->BMArray[p22]|>,
  "native_H0_projection"-><|
   "mixing_S1_S2bar"->BMArray[mu],
   "S1Inv1"->BMArray[proj11],"S1Inv2"->BMArray[proj12],
   "S2Inv1"->BMArray[proj21],"S2Inv2"->BMArray[proj22]|>,
  "exact_native_identity_portal_map"-><|
    "full_four_index_proof"->True,
    "formula"->"(Inv1 + Sqrt[3]*Inv2)/2 = delta_H delta_S",
    "Inv1_coefficient"->"1/2",
    "Inv2_coefficient"->"Sqrt[3]/2",
    "applies_to"->{"S1", "S2"},
    "frozen_coupling_assignment_is_selected_not_native_mandated"->True
  |>,
  "higgs_vev"->"H[2]=v/sqrt(2)",
  "matrix_convention"->"A_S2,S1=-(lambdaT3*v^2/2)*CG[T3MixCG][H2,H2,S1,S2bar]",
  "native_checks"->checks,"benchmark_portal_identification_certified"->False,
  "triplet_mass_metric_in_weight_basis_verified_by_live_native"->False,
  "absolute_Matchete_GammaCC_matrix_phase_measured"->False,
  "physical_source_certified"->False,"matching_or_running_performed"->False|>
];
$bmResult=CheckAbort[Catch[BMRun[],"BMStop"],$Aborted];
If[!AssociationQ[$bmResult] || Lookup[$bmResult,"status",""]=!="NativeCm1TripletScalarAndMajoranaActionExtracted_ExactPortalIdentityMap_PhysicalDictionaryPending",
 Print["BLOCKED: ",If[AssociationQ[$bmResult],Lookup[$bmResult,"reason","unknown"],"Wolfram abort"]];Exit[5]];
If[!DirectoryQ[DirectoryName[$bmOut]],
 CreateDirectory[DirectoryName[$bmOut],CreateIntermediateDirectories->True]];
If[Quiet@Check[Export[$bmOut,$bmResult,"RawJSON"],$Failed]===$Failed,
 Print["BLOCKED: native report export failed"];Exit[6]];
Print["C_-1 live triplet Majorana mass -1/2, neutral CGs +/-3^(-1/4), portal and mixing: EXTRACTED"];
Print["Native H0 mixing CG: -2/Sqrt[3]; S1,S2 identity portals: 1/2 Inv1 + Sqrt[3]/2 Inv2"];
Print["Witness: ",$bmOut];
Print["Physical source-certified models: 0/16; triplet metric and spinor phase pending"];
Exit[0];
