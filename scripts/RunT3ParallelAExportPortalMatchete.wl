(* A v5: all-16 REGISTERED native Matchete Higgs portal CGs with the
   corrected in-PlusHc phase. No source file or GitHub repository is changed.
   wolframscript -file scripts/RunT3ParallelAExportPortalMatchete.wl OUTPUT.json
*)
ClearAll["Global`*"];
If[Length[$ScriptCommandLine] =!= 2,
  Print["Usage: RunT3ParallelAExportPortalMatchete.wl OUTPUT.json"];Exit[64]];
out = ExpandFileName[$ScriptCommandLine[[2]]];
If[FileExistsQ[out],Print["Refusing to overwrite ",out];Exit[6]];
root = DirectoryName @ DirectoryName @ ExpandFileName[$InputFileName];
sourceFiles = {"Lagrangian/model/T3ModelCatalog.wl",
  "Lagrangian/model/LagrangianBuilder.wl",
  "Lagrangian/model/SU2Invariants.wl",
  "Lagrangian/model/T3Fields.wl",
  "Lagrangian/interactions/T3Topology.wl",
  "Lagrangian/interactions/ScalarPotential.wl",
  "Lagrangian/interactions/GeneralScalarOperators.wl"};
If[!AllTrue[sourceFiles,FileExistsQ[FileNameJoin[{root,#}]]&],
  Print["Required Matchete native source missing"];Exit[2]];
source = Import[FileNameJoin[{root,"Lagrangian","interactions","ScalarPotential.wl"}],"Text"];
If[!StringContainsQ[source,"T3HiggsPortalHermitianPhase"] ||
   !StringContainsQ[source,"phase coupling[] fieldProduct"] ||
   !StringContainsQ[source,"1/2 PlusHc["],
  Print["Corrected Hermitian portal phase inside PlusHc is not in the current source"];
  Exit[3]];
loaded=CheckAbort[Check[UsingFrontEnd[Needs["Matchete`"];True],$Failed],$Aborted];
If[loaded =!= True,Print["Matchete/FrontEnd unavailable"];Exit[3]];
Get[FileNameJoin[{root,"Lagrangian","model","T3ModelCatalog.wl"}]];
Get[FileNameJoin[{root,"Lagrangian","model","LagrangianBuilder.wl"}]];
T3A5Fail[s_] := (Print["A v5 Matchete ERROR: ",s];Throw[$Failed,"T3APortal5"]);
T3A5Text[x_] := ToString[FullSimplify[x],InputForm];
T3A5Portal[which_Integer,d_Integer,k_Integer,expr_] := Module[
  {cg,name,coupling,raw,tensor,flat,dims,expected,phase,hc,norm,herm,anti,kind,form},
  name="T3HiggsPortal"<>ToString[which]<>"CGInv"<>ToString[k];
  coupling="lambdaH"<>ToString[which]<>"Inv"<>ToString[k];
  form=ToString[expr,InputForm];
  If[!StringContainsQ[form,"CG["<>name<>","] ||
     !StringContainsQ[form,"Coupling["<>coupling<>","],
    T3A5Fail["Assembled Lagrangian lacks "<>name<>"/"<>coupling]];
  cg=Symbol[name];
  raw=CheckAbort[Quiet@Check[GetCGTensor[cg],$Failed],$Aborted];
  If[raw===$Failed || raw===$Aborted,
    T3A5Fail["Matchete`GetCGTensor failed for "<>name]];
  tensor=If[Head[raw]===SparseArray,Normal[raw],raw];
  dims=Dimensions[tensor];
  expected=Select[{2,2,d,d},#>1&];
  If[dims=!=expected,T3A5Fail["Native CG dimension mismatch "<>name<>
      " dims="<>ToString[dims,InputForm]]];
  flat=Flatten[{tensor}];
  If[Length[flat]=!=4 d d || !AllTrue[FullSimplify/@flat,NumericQ] ||
     !AnyTrue[flat,!TrueQ[PossibleZeroQ[#]]&],
    T3A5Fail["Incomplete/non-numeric registered CG "<>name]];
  phase=T3HiggsPortalHermitianPhase[cg,d];
  If[!MemberQ[{1,I},phase],T3A5Fail["Invalid real portal phase "<>name]];
  tensor=ArrayReshape[tensor,{2,2,d,d}];
  hc=Conjugate[Transpose[tensor,{2,1,4,3}]];
  norm=Max[Abs[Flatten[N[tensor,24]]]];
  herm=Max[Abs[Flatten[N[tensor-hc,24]]]];
  anti=Max[Abs[Flatten[N[tensor+hc,24]]]];
  If[!TrueQ[norm>10^-12],T3A5Fail["Zero tensor "<>name]];
  kind=Which[TrueQ[herm<10^-10 norm],"Hermitian",
             TrueQ[anti<10^-10 norm],"AntiHermitian",True,"Mixed"];
  If[kind==="Mixed" ||
    (kind==="Hermitian" && phase=!=1) ||
    (kind==="AntiHermitian" && phase=!=I),
    T3A5Fail["CG/PlusHc phase mismatch "<>name]];
  <|"which"->which,"invariant"->k,"coupling"->coupling,"cg"->name,
    "dimensions"->{2,2,d,d},"native_dimensions"->dims,
    "components_inputform"->(T3A5Text/@flat),
    "phase_inside_plus_hc"->If[phase===I,"I","1"],
    "native_hermiticity"->kind,
    "native_tensor_accessor"->"Matchete`GetCGTensor",
    "assembled_lagrangian_contains_cg_and_coupling"->True|>
];
T3A5Cases=Join[
  ( {"A",#}&/@{-4,-2,0}),({"B",#}&/@{-3,-1,1}),
  ({"C",#}&/@{-3,-1,1}),({"D",#}&/@{-2,0,2}),
  ({"E",#}&/@{-4,-2,0,2})];
T3A5Label[cls_,a_] := cls<>"_"<>If[a<0,"-","+"]<>ToString[Abs[a]];
(* The B/C alpha=-1 general quadratic is a mixed heavy-scalar mass.
   BuildT3Lagrangian retains the registered portal CGs and symbolic LInt,
   but intentionally refuses HeavyMassBasis and Match[].  Tensor extraction
   MAY use this symbolic result; physical matching/running MAY NOT.  No
   other failed full validation is accepted. *)
T3A5Model[cls_,a_] := Module[
  {tag,model,built,d1,d2,dF,expr,rows={},which,d,k,n,
   status,fullQ,massReadyQ,massBlocked,mixedQ,allowed,name},
  tag=T3A5Label[cls,a];Print["A v5 Matchete registered portals: ",tag];
  model=T3ModelFromClass[cls,a];
  If[!AssociationQ[model],T3A5Fail["Model catalog failed "<>tag]];
  {d1,d2,dF}={model["Scalar1","SU2"],model["Scalar2","SU2"],model["Fermion","SU2"]};
  (* Do NOT wrap the entire builder in Check: its expected B/C(-1)
     HeavyMassBasis messages would make Check discard the returned symbolic
     Association. Validate that Association and its explicit Status instead. *)
  built=CheckAbort[UsingFrontEnd[BuildT3Lagrangian[model]],$Aborted];
  If[!AssociationQ[built],T3A5Fail["Native builder did not return symbolic data "<>tag]];
  status=Lookup[built,"Status","MissingStatus"];
  fullQ=TrueQ[Lookup[built,"FullValidation",False]];
  massReadyQ=TrueQ[Lookup[built,"MatchingMassBasisReady",False]];
  massBlocked=Lookup[built,"AdditionalQuadraticsHeavyMassBasisBlocked",Missing["Count"]];
  mixedQ=status==="RequiresMixedScalarMassDiagonalization";
  Which[
    status==="Success" && fullQ && massReadyQ && massBlocked===0,
      Null,
    mixedQ && MemberQ[{"B_-1","C_-1"},tag] &&
      !fullQ && !massReadyQ && massBlocked===1 &&
      TrueQ[Lookup[built,"AdditionalQuadraticsRequested",False]] &&
      Lookup[built,"AdditionalQuadraticsExpected",-1]===1 &&
      Lookup[built,"AdditionalQuadraticsIncluded",-1]===1 &&
      TrueQ[Lookup[built,"T3IngredientsPresent",False]],
      Print["A v5: ",tag,
        " has a registered symbolic portal basis, but HEAVY MASS BASIS is blocked; no Match[]."],
    True,
      T3A5Fail["Native build status not valid for portal-only export "<>tag<>
        " status="<>ToString[status,InputForm]<>
        " mixedBlocked="<>ToString[massBlocked,InputForm]]
  ];
  expr=Lookup[built,"LInt",$Failed];
  allowed=Lookup[built,"AllowedInteractions",$Failed];
  If[expr===$Failed || !ListQ[allowed],
     T3A5Fail["Missing assembled LInt / allowed interaction inventory "<>tag]];
  Do[
    d=If[which===1,d1,d2];n=If[d===1,1,2];
    Do[
      name="HiggsPortal"<>ToString[which]<>"Inv"<>ToString[k];
      If[!MemberQ[allowed,name],
        T3A5Fail["Native portal not accepted in symbolic assembly "<>tag<>"/"<>name]];
      AppendTo[rows,T3A5Portal[which,d,k,expr]],
      {k,n}],
    {which,1,2}];
  <|"dimensions"->{d1,d2,dF},"alpha"->a,
    "native_build_status"->status,
    "full_lagrangian_validation"->fullQ,
    "matching_mass_basis_ready"->massReadyQ,
    "mixed_scalar_mass_basis_blocked"->mixedQ,
    "additional_quadratics_heavy_mass_basis_blocked"->massBlocked,
    "portal_cg_operator_inventory_verified"->True,
    "portals"->rows|>
];
(* As above, an outer Check would reject the expected mass-basis warnings
   even after every model passed its explicit symbolic-build status gate. *)
run=CheckAbort[Catch[Association@Table[
  With[{cls=First[c],a=Last[c]},T3A5Label[cls,a]->T3A5Model[cls,a]],
  {c,T3A5Cases}],"T3APortal5"],$Aborted];
If[!AssociationQ[run] || Length[Keys[run]]=!=16 ||
   Total[Length[Lookup[#,"portals",{}]]&/@Values[run]]=!=58,
  Print["A v5 Matchete all-model native portal export BLOCKED"];Exit[5]];
a4=Select[run["A_-4"]["portals"],#["coupling"] === "lambdaH2Inv1" &];
If[Length[a4]=!=1 || Lookup[First[a4],"phase_inside_plus_hc",""]=!="I",
  Print["A_-4 corrected anti-Hermitian portal phase missing"];Exit[5]];
sha=Association@Table[
  s->IntegerString[FileHash[FileNameJoin[{root,s}],"SHA256"],16,64],
  {s,sourceFiles}];
report=<|"format"->"t3-parallel-a-native-matchete-portal-cg-phases-v5",
  "status"->"NativeRegisteredPortalsWithCorrectedRealHermitianPhases",
  "model_count"->16,"portal_coupling_slots"->58,
  "models_requiring_mixed_scalar_mass_diagonalization"->
    Select[Keys[run],TrueQ[run[#]["mixed_scalar_mass_basis_blocked"]]&],
  "mixed_scalar_mass_basis_matching_authorized"->False,
  "source_sha256"->sha,
  "source_potential_convention"->"real lambda /2 PlusHc[phase*CG*fields]",
  "native_RGBeta_generator_components_certified"->False,
  "physical_coupling_map_certified"->False,"models"->run|>;
If[!DirectoryQ[DirectoryName[out]],
  CreateDirectory[DirectoryName[out],CreateIntermediateDirectories->True]];
write=Quiet@Check[Export[out,report,"RawJSON"],$Failed];
If[write===$Failed || !FileExistsQ[out],
  Print["A v5 Matchete JSON write failed"];Exit[7]];
Print["58/58 native Matchete Higgs-portal CGs with checked in-PlusHc phases -> ",out];
