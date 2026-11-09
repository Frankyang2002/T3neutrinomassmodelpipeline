(* Workstream A v16: native Matchete scalar self/cross CGs, all 16 models.
   Read-only source registration. Never Match[], RG-run, or integrate out fields.
   wolframscript -file scripts/RunT3ParallelAExportNonportalMatcheteV16.wl OUT.json
*)
ClearAll["Global`*"];
If[Length[$ScriptCommandLine] =!= 2,
  Print["Usage: RunT3ParallelAExportNonportalMatcheteV16.wl OUTPUT.json"]; Exit[64]];
T3A16Out = ExpandFileName[$ScriptCommandLine[[2]]];
If[FileExistsQ[T3A16Out] || FileExistsQ[T3A16Out <> ".partial"],
  Print["v16 refuses existing output or partial: ", T3A16Out]; Exit[6]];
T3A16Root = DirectoryName @ DirectoryName @ ExpandFileName[$InputFileName];
T3A16Sources = {
  "Lagrangian/model/T3ModelCatalog.wl",
  "Lagrangian/model/LagrangianBuilder.wl",
  "Lagrangian/model/SU2Invariants.wl",
  "Lagrangian/model/T3Fields.wl",
  "Lagrangian/interactions/T3Topology.wl",
  "Lagrangian/interactions/ScalarPotential.wl",
  "Lagrangian/interactions/GeneralScalarOperators.wl",
  "scripts/RunT3ParallelAExportNonportalMatcheteV16.wl"
};
If[!AllTrue[T3A16Sources, FileExistsQ[FileNameJoin[{T3A16Root, #}]] &],
  Print["v16 required native source missing"]; Exit[2]];
T3A16PotentialSource = Import[FileNameJoin[{T3A16Root,
  "Lagrangian","interactions","ScalarPotential.wl"}], "Text"];
If[!And @@ (StringContainsQ[T3A16PotentialSource, #] & /@ {
    "BuildScalarSelfInvariantCandidates", "BuildCrossScalarInvariantCandidates",
    "PhysicalSU2InvariantBasis", "SymmetrizeTensorPairs"}),
  Print["v16 Matchete scalar-source contract changed"]; Exit[3]];
T3A16Loaded = CheckAbort[Check[UsingFrontEnd[Needs["Matchete`"]; True], $Failed], $Aborted];
If[T3A16Loaded =!= True, Print["v16 Matchete/FrontEnd unavailable"]; Exit[3]];
Get[FileNameJoin[{T3A16Root,"Lagrangian","model","T3ModelCatalog.wl"}]];
Get[FileNameJoin[{T3A16Root,"Lagrangian","model","LagrangianBuilder.wl"}]];
T3A16Error[message_String] := (Print["v16 BLOCKED: ", message]; Throw[$Failed, "T3A16"]);
T3A16Cases = Join[( {"A",#}& /@ {-4,-2,0}),
   ({"B",#}& /@ {-3,-1,1}), ({"C",#}& /@ {-3,-1,1}),
   ({"D",#}& /@ {-2,0,2}), ({"E",#}& /@ {-4,-2,0,2})];
T3A16Label[c_,a_] := c <> "_" <> If[a < 0,"-","+"] <> ToString[Abs[a]];
T3A16ExpectedCount[d1_,d2_,family_] := Switch[family,
  "self1", If[d1===3,2,1],
  "self2", If[d2===3,2,1],
  "cross", Which[d1===3 && d2===3,3,d1>1 && d2>1,2,True,1],
  _, 0];
T3A16Flatten[t_] := If[Head[t]===SparseArray,Flatten[Normal[t]],Flatten[{t}]];
T3A16Entry[family_,k_Integer,d1_Integer,d2_Integer,expr_,allowed_] := Module[
  {symbol,cgSymbol,operator,symbolName,cgName,dim,native,flat,expected,form,
   raw,sourceFields,actionFactor},
  Switch[family,
    "self1", dim=d1;symbolName="lambdaS1Inv"<>ToString[k];
      cgName="T3Scalar1SelfCGInv"<>ToString[k];
      operator="Scalar1SelfInv"<>ToString[k];
      expected={d1,d1,d1,d1}; actionFactor="1/2";
      sourceFields={"Bar[S1]","S1","Bar[S1]","S1"},
    "self2", dim=d2;symbolName="lambdaS2Inv"<>ToString[k];
      cgName="T3Scalar2SelfCGInv"<>ToString[k];
      operator="Scalar2SelfInv"<>ToString[k];
      expected={d2,d2,d2,d2};actionFactor="1/2";
      sourceFields={"Bar[S2]","S2","Bar[S2]","S2"},
    "cross",dim=0;symbolName="lambda12Inv"<>ToString[k];
      cgName="T3CrossScalarCGInv"<>ToString[k];
      operator="ScalarCrossPortalInv"<>ToString[k];
      expected={d1,d1,d2,d2};actionFactor="1";
      sourceFields={"Bar[S1]","S1","Bar[S2]","S2"},
    _,T3A16Error["invalid family"]];
  If[!MemberQ[allowed,operator],T3A16Error["missing allowed native operator "<>operator]];
  form=ToString[expr,InputForm];
  If[!StringContainsQ[form,"Coupling["<>symbolName<>","],
    T3A16Error["native LInt missing coupling "<>symbolName]];
  If[Times@@expected===1,
    raw=1;native={};flat={"1"};cgName="None",
    If[!StringContainsQ[form,"CG["<>cgName<>","],
      T3A16Error["native LInt missing registered CG "<>cgName]];
    cgSymbol=Symbol[cgName];
    raw=CheckAbort[Quiet@Check[GetCGTensor[cgSymbol],$Failed],$Aborted];
    If[raw===$Failed || raw===$Aborted,
      T3A16Error["GetCGTensor failed for "<>cgName]];
    native=Dimensions[raw];flat=ToString[FullSimplify[#],InputForm]& /@ T3A16Flatten[raw];
  ];
  If[Length[flat] =!= Times@@expected ||
     !AllTrue[flat,StringQ] ||
     !AnyTrue[flat,#=!= "0" &],
    T3A16Error["incomplete native tensor "<>symbolName]];
  If[native =!= Select[expected,#>1&] && native =!= {},
    T3A16Error["unexpected native CG dimensions for "<>symbolName<>" "<>
      ToString[native,InputForm]]];
  <|"family"->family,"invariant"->k,"coupling"->symbolName,
    "cg"->cgName,"operator"->operator,
    "dimensions"->expected,"native_dimensions"->native,
    "fields"->sourceFields,"matchete_potential_factor"->actionFactor,
    "components_inputform"->flat,"native_tensor_accessor"->
      If[cgName==="None","Scalar singlet 1","Matchete`GetCGTensor"],
    "assembled_lagrangian_contains_cg_and_coupling"->True|>
];
T3A16Model[cl_,a_] := Module[{tag,model,built,status,fullQ,readyQ,blocked,mixedQ,
   d1,d2,df,allowed,expr,rows={},family,n,k},
  tag=T3A16Label[cl,a]; Print["v16 native Matchete nonportal CG: ",tag];
  model=T3ModelFromClass[cl,a];
  If[!AssociationQ[model],T3A16Error["catalog failed "<>tag]];
  {d1,d2,df}={model["Scalar1","SU2"],model["Scalar2","SU2"],model["Fermion","SU2"]};
  built=CheckAbort[UsingFrontEnd[BuildT3Lagrangian[model]],$Aborted];
  If[!AssociationQ[built],T3A16Error["builder returned no symbolic Association "<>tag]];
  status=Lookup[built,"Status",None];
  fullQ=TrueQ[Lookup[built,"FullValidation",False]];
  readyQ=TrueQ[Lookup[built,"MatchingMassBasisReady",False]];
  blocked=Lookup[built,"AdditionalQuadraticsHeavyMassBasisBlocked",Missing["Count"]];
  mixedQ=status==="RequiresMixedScalarMassDiagonalization";
  Which[
    status==="Success" && fullQ && readyQ && blocked===0,Null,
    MemberQ[{"B_-1","C_-1"},tag] && mixedQ && !fullQ && !readyQ &&
      blocked===1 && TrueQ[Lookup[built,"AdditionalQuadraticsRequested",False]] &&
      Lookup[built,"AdditionalQuadraticsExpected",-1]===1 &&
      Lookup[built,"AdditionalQuadraticsIncluded",-1]===1 &&
      TrueQ[Lookup[built,"T3IngredientsPresent",False]],
      Print["v16 ",tag," registered CG ONLY; mixed scalar mass basis blocked"],
    True,T3A16Error["unexpected model build status "<>tag<>" "<>ToString[status,InputForm]]
  ];
  expr=Lookup[built,"LInt",$Failed]; allowed=Lookup[built,"AllowedInteractions",$Failed];
  If[expr===$Failed || !ListQ[allowed],T3A16Error["native LInt or operator inventory missing "<>tag]];
  Do[
    n=T3A16ExpectedCount[d1,d2,family];
    Do[AppendTo[rows,T3A16Entry[family,k,d1,d2,expr,allowed]],{k,n}],
    {family,{"self1","self2","cross"}}];
  <|"dimensions"->{d1,d2,df},"alpha"->a,
    "native_build_status"->status,"matching_mass_basis_ready"->readyQ,
    "mixed_scalar_mass_basis_blocked"->mixedQ,
    "scalar_cg_operator_inventory_verified"->True,"nonportal"->rows|>
];
T3A16Run=CheckAbort[Catch[Association@Table[
   With[{cl=First[c],a=Last[c]},T3A16Label[cl,a]->T3A16Model[cl,a]],
   {c,T3A16Cases}],"T3A16"],$Aborted];
If[!AssociationQ[T3A16Run] || Length[T3A16Run]=!=16 ||
   Total[Length[Lookup[#,"nonportal",{}]]& /@ Values[T3A16Run]]=!=76,
  Print["v16 16-model/76-slot native export BLOCKED"]; Exit[5]];
T3A16Hashes=Association@Table[
  p->IntegerString[FileHash[FileNameJoin[{T3A16Root,p}],"SHA256"],16,64],
  {p,T3A16Sources}];
T3A16Data=<|"format"->"t3-parallel-a-native-matchete-nonportal-cg-v16",
  "status"->"NativeMatchete76NonportalCGComponentsCaptured_MassBasisStillConditional",
  "models_checked"->16,"nonportal_slots"->76,
  "source_sha256"->T3A16Hashes,"models"->T3A16Run,
  "physical_scalar_dictionary_certified"->False,
  "physical_Wilson_evolution_authorized"->False,
  "EFT_order"->"UV -> integrate out ALL fermions -> integrate out ALL scalars -> SMEFT"|>;
If[!DirectoryQ[DirectoryName[T3A16Out]],
  CreateDirectory[DirectoryName[T3A16Out],CreateIntermediateDirectories->True]];
T3A16Tmp=T3A16Out<>".partial";
If[Quiet@Check[Export[T3A16Tmp,T3A16Data,"RawJSON"],$Failed]===$Failed,
  Print["v16 could not write native diagnostic"];Exit[7]];
T3A16Read=Quiet@Check[Import[T3A16Tmp,"RawJSON"],$Failed];
If[!AssociationQ[T3A16Read] || Lookup[T3A16Read,"nonportal_slots",-1]=!=76,
  Print["v16 output readback check failed"];Exit[7]];
If[Quiet@Check[RenameFile[T3A16Tmp,T3A16Out],$Failed]===$Failed,
  Print["v16 output atomic rename failed"];Exit[7]];
Print["v16 native Matchete: 76/76 self/cross scalar CGs -> ",T3A16Out];
Print["B_-1/C_-1 registered CG ONLY; no physical mixed-mass matching/Wilson evolution"];
Exit[0];
