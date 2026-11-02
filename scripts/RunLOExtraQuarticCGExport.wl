(* V6 live source export of Matchete's ACTUALLY REGISTERED additional quartic CGs.
   Does not call Match, IntegrateOut, RGBeta, or any numerical solver.
   Required: $env:T3_EXPERIMENTAL_GENERAL_QUARTICS='1';
   Disable experimental quadratics and fermions for a pure scalar-source audit.
   Usage: wolframscript -file scripts/RunLOExtraQuarticCGExport.wl OUT.json A m2
   Output is restricted to the Workstream 2 owner directory. *)
ClearAll["Global`*"];
v6Args = Take[$ScriptCommandLine, -3];
If[Length[v6Args] =!= 3, Print["Usage: OUT.json A|B|C|D|E mN|pN"]; Exit[2]];
v6Dir = DirectoryName @ ExpandFileName[$InputFileName];
v6Root = DirectoryName[v6Dir];
v6Out = ExpandFileName[v6Args[[1]]];
v6Owner = ExpandFileName @ FileNameJoin[{v6Root, "output", "general_observable", "lo_certification_couplings"}];
If[!StringStartsQ[v6Out, v6Owner <> $PathnameSeparator],
  Print["Workstream 2 output directory violation"]; Exit[2]];
v6Class = ToUpperCase[v6Args[[2]]];
v6Tag = v6Args[[3]];
v6Alpha = Which[
 StringMatchQ[v6Tag, "m" ~~ DigitCharacter ..], -ToExpression[StringDrop[v6Tag, 1]],
 StringMatchQ[v6Tag, "p" ~~ DigitCharacter ..], ToExpression[StringDrop[v6Tag, 1]],
 True, $Failed];
If[v6Alpha === $Failed || !MemberQ[{"A","B","C","D","E"},v6Class], Exit[2]];
v6ModelName = v6Class <> "_" <> If[v6Alpha < 0,"-","+"] <> ToString[Abs[v6Alpha]];
If[Environment["T3_EXPERIMENTAL_GENERAL_QUARTICS"] =!= "1" ||
   Environment["T3_EXPERIMENTAL_GENERAL_QUADRATICS"] === "1" ||
   Environment["T3_EXPERIMENTAL_GENERAL_FERMIONS"] === "1",
 Print["Require quartics=1; quadratics/fermions !=1 for this source-only audit"]; Exit[2]];

v6Loaded = CheckAbort[Check[UsingFrontEnd[Needs["Matchete`"] ; True], $Failed], $Aborted];
If[v6Loaded =!= True, Print["Matchete unavailable"]; Exit[3]];
Get[FileNameJoin[{v6Root,"Lagrangian","model","T3ModelCatalog.wl"}]];
Get[FileNameJoin[{v6Root,"Lagrangian","model","LagrangianBuilder.wl"}]];
v6Model = T3ModelFromClass[v6Class,v6Alpha];
If[!AssociationQ[v6Model], Print["Unsupported model"]; Exit[2]];
ClearAll[V6Fail,V6CGRead,V6ExportOne,V6ExportAll];
V6Fail[s_String] := Throw[<|"status"->"LiveExportBlocked", "reason"->s|>,"V6Stop"];
V6CGRead[name_String] := Module[{raw,arr,dims,flat},
 raw=CheckAbort[Quiet@Check[GetCGTensor[Symbol[name]], $Failed], $Aborted];
 If[raw === $Failed || raw === $Aborted,V6Fail["GetCGTensor failed: " <> name]];
 arr=If[Head[raw]===SparseArray,Normal[raw],raw];
 dims=Dimensions[arr];flat=Flatten[arr];
 If[dims==={} || Times@@dims=!=Length[flat] ||
    !AllTrue[flat,TrueQ[NumericQ[FullSimplify[#]]] &] ||
    !AnyTrue[flat,!TrueQ[PossibleZeroQ[#]] &],
   V6Fail["Malformed native CG components: " <> name]];
 <|"CG_dimensions"->dims,
   "CG_components_inputform"->(ToString[FullSimplify[#],InputForm]& /@ flat),
   "implicit_singlet_CG"->False|>
];
V6ExportOne[counts_List,idx_Integer,allowed_List] := Module[
 {key,opName,coupling,cgName,names,dims,info},
 key=StringJoin[ToString /@ counts];
 opName="T3ExtraQ_"<>StringRiffle[ToString /@ counts,"_"]<>"_"<>ToString[idx];
 coupling="lambdaT3ExtraQ"<>key<>"I"<>ToString[idx];
 cgName="T3ExtraCGQ"<>key<>"I"<>ToString[idx];
 If[!MemberQ[allowed,opName],V6Fail["Validated interaction missing: "<>opName]];
 names=T3ExtendedSpeciesList[counts];
 dims=T3ExtendedFieldDimensions[v6Model,names];
 info=If[AllTrue[dims,#===1&],
   <|"CG_dimensions"->{},"CG_components_inputform"->{"1"},"implicit_singlet_CG"->True|>,
   V6CGRead[cgName]];
 Join[<|"canonical_field_counts"->counts,"invariant_index"->idx,
    "operator_name"->opName,"coupling_symbol"->coupling,"CG_symbol"->cgName,
    "present_in_validated_Lagrangian"->True,
    "source_interaction_normalization"->"PlusHc[coupling[] product CG] (source literal)",
    "physical_RGBeta_mapping_certified"->False|>,info]
];
V6ExportAll[] := Module[{built,verified,expected,records},
 built=CheckAbort[Check[BuildT3Lagrangian[v6Model],$Failed],$Aborted];
 If[!AssociationQ[built] || Lookup[built,"Status",""]=!="Success" ||
    !TrueQ[Lookup[built,"FullValidation",False]],
   V6Fail["Full Matchete Lagrangian build/validation was not successful"]];
 verified=T3ExtendedVerifiedEntries[v6Model];
 If[!ListQ[verified],V6Fail["Verified 16-model quartic table unavailable"]];
 expected=Total[verified[[All,2]]];
 If[Lookup[built,"AdditionalQuarticsIncluded",-1]=!=expected ||
    Lookup[built,"AdditionalQuarticsExpected",-1]=!=expected,
   V6Fail["Generator and assembled extra quartic counts disagree"]];
 records=Flatten[Table[
    V6ExportOne[verified[[k,1]],i,built["AllowedInteractions"]],
    {k,Length[verified]},{i,verified[[k,2]]}],1];
 If[Length[records]=!=expected,V6Fail["Incomplete native tensor export"]];
 <|"format"->"t3-ws2-live-Matchete-extra-quartic-CGs-v6",
   "model"->v6ModelName,"status"->"LiveExtraQuarticCGExportVerified",
   "source_level"->"experimental_UV_Lagrangian_not_EFT1_matching_seed",
   "full_matchete_lagrangian_validation_passed"->True,
   "extra_quartic_count"->expected,
   "native_GetCGTensor_used"->True,
   "operator_normalization_source"->"T3ExtendedBuildCandidates: PlusHc[lambda[] product CG]",
   "operators"->records,
   "physical_Matchete_RGBeta_conversion_certified"->False,
   "EFT1_running_performed"->False,
   "matching_performed"->False|>
];
v6Result=CheckAbort[Catch[UsingFrontEnd[V6ExportAll[]],"V6Stop"],
   <|"status"->"LiveExportBlocked","reason"->"Aborted during Matchete live extraction"|>];
If[!AssociationQ[v6Result],v6Result=<|"status"->"LiveExportBlocked","reason"->"Unexpected Wolfram result"|>];
If[!DirectoryQ[DirectoryName[v6Out]],CreateDirectory[DirectoryName[v6Out],CreateIntermediateDirectories->True]];
If[Quiet@Check[Export[v6Out,v6Result,"RawJSON"],$Failed]===$Failed,
 Print["Could not write V6 live output"]; Exit[7]];
Print[v6ModelName,": ",v6Result["status"]," -> ",v6Out];
If[v6Result["status"] =!= "LiveExtraQuarticCGExportVerified",Exit[6]];
Exit[0];
