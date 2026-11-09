(* v18: source-only Matchete CG witness, all16 topology or B/C(-1) extra quartics.
   Never Match[], HeavyMassBasis[], RGBeta running, or physical Wilson running.
   wolframscript -file scripts/RunT3ParallelAExportComplexMatcheteV18.wl OUTPUT.json TOPOLOGY|BCEXTRA
*)
ClearAll["Global`*"];
If[Length[$ScriptCommandLine] =!= 3,
 Print["Usage: OUTPUT.json TOPOLOGY|BCEXTRA"]; Exit[64]];
T3A18Out=ExpandFileName[$ScriptCommandLine[[2]]];
T3A18Mode=$ScriptCommandLine[[3]];
If[!MemberQ[{"TOPOLOGY","BCEXTRA"},T3A18Mode],Print["Invalid v18 mode"];Exit[64]];
If[FileExistsQ[T3A18Out] || FileExistsQ[T3A18Out<>".partial"],
 Print["Refusing overwrite of native v18 evidence ",T3A18Out];Exit[6]];
T3A18Flags={"T3_EXPERIMENTAL_GENERAL_FERMIONS","T3_EXPERIMENTAL_GENERAL_QUADRATICS",
 "T3_EXPERIMENTAL_GENERAL_QUARTICS","T3_EXPERIMENTAL_FLAVOR_MASSES"};
T3A18Want=(T3A18Mode==="BCEXTRA");
If[AnyTrue[DeleteCases[T3A18Flags,"T3_EXPERIMENTAL_GENERAL_QUARTICS"],Environment[#]==="1"&] ||
   (Environment["T3_EXPERIMENTAL_GENERAL_QUARTICS"]==="1")=!=T3A18Want,
 Print["v18 BLOCKED: optional flags incorrect for mode ",T3A18Mode];Exit[3]];
T3A18Root=DirectoryName@DirectoryName@ExpandFileName[$InputFileName];
T3A18Sources={
 "Lagrangian/model/T3ModelCatalog.wl",
 "Lagrangian/model/LagrangianBuilder.wl",
 "Lagrangian/model/SU2Invariants.wl",
 "Lagrangian/model/T3Fields.wl",
 "Lagrangian/interactions/T3Topology.wl",
 "Lagrangian/interactions/ScalarPotential.wl",
 "Lagrangian/interactions/GeneralScalarOperators.wl",
 "Lagrangian/interactions/GeneralScalarVerifiedQuartics.wl",
 "scripts/RunT3ParallelAExportComplexMatcheteV18.wl"};
If[!AllTrue[T3A18Sources,FileExistsQ[FileNameJoin[{T3A18Root,#}]]&],
 Print["v18 source files missing"];Exit[2]];
T3A18Loaded=CheckAbort[Check[UsingFrontEnd[Needs["Matchete`"]];True,$Failed],$Aborted];
If[T3A18Loaded=!=True,Print["v18 Matchete/FrontEnd unavailable"];Exit[4]];
Get[FileNameJoin[{T3A18Root,"Lagrangian","model","T3ModelCatalog.wl"}]];
Get[FileNameJoin[{T3A18Root,"Lagrangian","model","LagrangianBuilder.wl"}]];
T3A18Fail[s_String]:=(Print["v18 BLOCKED: ",s];Throw[$Failed,"T3A18"]);
T3A18Cases=If[T3A18Want,{{"B",-1},{"C",-1}},
 Join[({"A",#}&/@{-4,-2,0}),({"B",#}&/@{-3,-1,1}),
 ({"C",#}&/@{-3,-1,1}),({"D",#}&/@{-2,0,2}),
 ({"E",#}&/@{-4,-2,0,2})]];
T3A18Label[c_,a_]:=c<>"_"<>If[a<0,"-","+"]<>ToString[Abs[a]];
T3A18Flatten[raw_]:=Flatten[If[Head[raw]===SparseArray,Normal[raw],{raw}]];
T3A18CG[cg_,dims_List] := Module[{raw,native,flat},
 raw=CheckAbort[Quiet@Check[GetCGTensor[cg],$Failed],$Aborted];
 If[raw===$Failed || raw===$Aborted,T3A18Fail["native GetCGTensor failed for "<>SymbolName[cg]]];
 native=Dimensions[raw];
 flat=ToString[FullSimplify[#],InputForm]& /@ T3A18Flatten[raw];
 If[native=!=Select[dims,#>1&] || Length[flat]=!=Times@@dims ||
    !AllTrue[flat,StringQ] || !AnyTrue[flat,#=!= "0"&],
    T3A18Fail["native CG shape/components invalid: "<>SymbolName[cg]]];
 <|"native_dimensions"->native,"components_inputform"->flat|>
];
T3A18Model[cl_,alpha_] := Module[{tag,model,built,allowed,expr,d1,d2,df,status,
 rows={},cg,cgName,name,key,counts,names,dims,expected,k,n,record,extras},
 tag=T3A18Label[cl,alpha]; Print["v18 native complex CG: ",tag," / ",T3A18Mode];
 model=T3ModelFromClass[cl,alpha];
 If[!AssociationQ[model],T3A18Fail["invalid model catalog "<>tag]];
 {d1,d2,df}={model["Scalar1","SU2"],model["Scalar2","SU2"],model["Fermion","SU2"]};
 built=CheckAbort[UsingFrontEnd[BuildT3Lagrangian[model]],$Aborted];
 If[!AssociationQ[built],T3A18Fail["native builder did not produce Association "<>tag]];
 status=Lookup[built,"Status",None];
 If[status=!= "Success" || !TrueQ[Lookup[built,"T3IngredientsPresent",False]],
   T3A18Fail["native symbolic Lagrangian failed "<>tag<>" "<>ToString[status,InputForm]]];
 If[TrueQ[Lookup[built,"AdditionalQuadraticsRequested",False]] ||
    TrueQ[Lookup[built,"AdditionalYukawasRequested",False]],
   T3A18Fail["unexpected B/C mixed-mass/Yukawa extension "<>tag]];
 allowed=Lookup[built,"AllowedInteractions",$Failed];expr=Lookup[built,"LInt",$Failed];
 If[!ListQ[allowed] || expr===$Failed,T3A18Fail["LInt/operator list missing "<>tag]];
 If[!T3A18Want,
   If[!MemberQ[allowed,"T3ScalarMix_GenericCG"] ||
      !StringContainsQ[ToString[expr,InputForm],"CG[T3MixCG"],
      T3A18Fail["source topology CG not assembled "<>tag]];
   record=T3A18CG[T3MixCG,{2,2,d1,d2}];
   rows={<|"symbol"->"lambdaT3","cg"->"T3MixCG",
      "dimensions"->{2,2,d1,d2},"fields"->{"H","H","S1","Bar[S2]"},
      "source_term"->"PlusHc[lambdaT3[] H H S1 Bar[S2] CG[T3MixCG]]",
      "assembled_lagrangian_contains_operator"->True|>~Join~record},
   If[!TrueQ[Lookup[built,"AdditionalQuarticsIncluded",-1]===7],
     T3A18Fail["expected seven B/C extra quartics not all included "<>tag]];
   extras=T3ExtendedExtraQuarticStructures[model];
   Do[
     counts=extras[[i]];key=StringJoin[ToString/@counts];
     names=T3ExtendedSpeciesList[counts];
     dims=T3ExtendedFieldDimensions[model,names];
     n=Lookup[T3ExtendedStructureExpectedMultiplicities[model],
       T3ExtendedStructureKey[counts],0];
     If[!IntegerQ[n] || n<=0,T3A18Fail["no native B/C structure multiplicity "<>tag]];
     Do[
       name="T3ExtraQ_"<>T3ExtendedStructureKey[counts]<>"_"<>ToString[k];
       cgName="T3ExtraCGQ"<>key<>"I"<>ToString[k];
       If[!MemberQ[allowed,name] ||
          !StringContainsQ[ToString[expr,InputForm],"CG["<>cgName<>""],
         T3A18Fail["native B/C extra operator not assembled "<>name]];
       cg=Symbol[cgName];record=T3A18CG[cg,dims];
       AppendTo[rows,<|"symbol"->"lambdaT3ExtraQ"<>key<>"I"<>ToString[k],
         "cg"->cgName,"operator"->name,"counts"->counts,"fields"->names,
         "dimensions"->dims,"assembled_lagrangian_contains_operator"->True|>~Join~record],
       {k,n}],{i,Length[extras]}];
   If[Length[rows]=!=7,T3A18Fail["B/C source extra CG count !=7 "<>tag]]
 ];
 <|"dimensions"->{d1,d2,df},"alpha"->alpha,
  "native_builder_status"->status,
  "matching_mass_basis_physically_certified"->False,
  "symbolic_quartic_only"->True,"records"->rows|>
];
T3A18Data=CheckAbort[Catch[Association@Table[
 With[{cl=First[c],a=Last[c]},T3A18Label[cl,a]->T3A18Model[cl,a]],
 {c,T3A18Cases}],"T3A18"],$Aborted];
T3A18Expected=If[T3A18Want,14,16];
If[!AssociationQ[T3A18Data] || Length[T3A18Data]=!=Length[T3A18Cases] ||
 Total[Length[Lookup[#,"records",{}]]&/@Values[T3A18Data]]=!=T3A18Expected,
 Print["v18 native complex CG inventory incomplete"];Exit[5]];
T3A18Hashes=Association@Table[p->IntegerString[
 FileHash[FileNameJoin[{T3A18Root,p}],"SHA256"],16,64],{p,T3A18Sources}];
T3A18Report=<|"format"->"t3-parallel-a-matchete-complex-cg-v18",
 "status"->"NativeMatcheteComplexCGSourceOnly_MassAndBasisConditional",
 "mode"->T3A18Mode,"models_checked"->Length[T3A18Data],
 "native_complex_CG_slots"->T3A18Expected,
 "models"->T3A18Data,"source_sha256"->T3A18Hashes,
 "independent_RGBeta_native_component_frame_certified"->False,
 "BC_minus1_mixed_scalar_mass_diagonalized"->False,
 "physical_Wilson_evolution_authorized"->False,
 "EFT_order"->"UV -> integrate out ALL fermions -> integrate out ALL scalars -> SMEFT"|>;
If[!DirectoryQ[DirectoryName[T3A18Out]],
 CreateDirectory[DirectoryName[T3A18Out],CreateIntermediateDirectories->True]];
T3A18Tmp=T3A18Out<>".partial";
If[Quiet@Check[Export[T3A18Tmp,T3A18Report,"RawJSON"],$Failed]===$Failed,
 Print["v18 output write failed"];Exit[7]];
T3A18Recheck=Quiet@Check[Import[T3A18Tmp,"RawJSON"],$Failed];
If[!AssociationQ[T3A18Recheck] ||
 Lookup[T3A18Recheck,"native_complex_CG_slots",-1]=!=T3A18Expected,
 Print["v18 JSON readback failed"];Exit[7]];
If[Quiet@Check[RenameFile[T3A18Tmp,T3A18Out],$Failed]===$Failed,
 Print["v18 output final rename failed"];Exit[7]];
Print["v18 ",T3A18Mode,": ",T3A18Expected," native complex Matchete CG slots -> ",T3A18Out];
Print["v18 CG source only; B/C mixed scalar mass basis still NOT certified."];
Exit[0];
