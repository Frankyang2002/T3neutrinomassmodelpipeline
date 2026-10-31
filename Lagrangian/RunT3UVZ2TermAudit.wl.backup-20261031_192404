(* Workstream A: 16-model read-only Z2 term audit.
   wolframscript -file Lagrangian/RunT3UVZ2TermAudit.wl OUTPUT.json A|B|C|D|E m4|m3|m2|m1|p0|p1|p2

   Rebuilds the COMPLETE BSM free+interaction Lagrangian with general fermions,
   indexed MF and quartics, but deliberately WITHOUT extra scalar bilinears;
   those require the separate 8-model eigenbasis drivers.  Counts odd-parity
   NEW field occurrences in each expanded term.  Never calls Match[].

   Evenness of the assembled expression does NOT mean Z2 is registered as a
   discrete symmetry with Matchete, or prove all possible fermion operators.
*)
ClearAll["Global`*"];
$T3ZAArgs=Rest[$ScriptCommandLine];
If[Length[$T3ZAArgs]=!=3,Print["Usage: OUTPUT.json FAMILY alphaTag"];Exit[2]];
$T3ZAOut=ExpandFileName[$T3ZAArgs[[1]]];
$T3ZAFamily=ToUpperCase[$T3ZAArgs[[2]]];
$T3ZAAlphaTag=$T3ZAArgs[[3]];
If[!MemberQ[{"A","B","C","D","E"},$T3ZAFamily] ||
   !StringMatchQ[$T3ZAAlphaTag, ("m"|"p") ~~ DigitCharacter ..],
   Print["Invalid model arguments"];Exit[2]];
$T3ZAAlpha=If[StringStartsQ[$T3ZAAlphaTag,"m"],-1,1] *
  ToExpression[StringDrop[$T3ZAAlphaTag,1]];
$T3ZAModel=$T3ZAFamily<>"_"<>$T3ZAAlphaTag;
$T3ZADir=DirectoryName@ExpandFileName[$InputFileName];
$T3ZAFlags={"T3_EXPERIMENTAL_FLAVOR_MASSES",
  "T3_EXPERIMENTAL_GENERAL_FERMIONS", "T3_EXPERIMENTAL_GENERAL_QUARTICS"};
If[!AllTrue[$T3ZAFlags,Environment[#]==="1"&] ||
    Environment["T3_EXPERIMENTAL_GENERAL_QUADRATICS"]==="1",
  Print["Require flavour/fermion/quartic flags=1; quadratics=0"];Exit[2]];
(* Load Matchete BEFORE parsing any unqualified Matchete public symbol. *)
$T3ZALoad=CheckAbort[UsingFrontEnd[Needs["Matchete`"] ; True],$Aborted];
If[$T3ZALoad=!=True,Print["Matchete could not initialize"];Exit[3]];
Get[FileNameJoin[{$T3ZADir,"model","T3ModelCatalog.wl"}]];
Get[FileNameJoin[{$T3ZADir,"model","LagrangianBuilder.wl"}]];

T3ZATerms[expr_] := Module[{expanded=Expand[expr]},
  If[Head[expanded]===Plus,List@@expanded,{expanded}]
];
T3ZABSMFields[term_] := Cases[term,
  HoldPattern[Field[name_Symbol,___]] /;
    MemberQ[{NewFermion,NewScalar1,NewScalar2},name] :> SymbolName[name],
  {0,Infinity}];
T3ZAInspect[expr_] := Module[{terms,fields,count,odd,zero,distribution},
  terms=T3ZATerms[expr];
  fields=T3ZABSMFields /@ terms;
  count=Length /@ fields;
  odd=Flatten@Position[Mod[count,2],1];
  zero=Flatten@Position[count,0];
  distribution=Counts[count];
  <|"expanded_term_count"->Length[terms],
    "odd_bsm_term_count"->Length[odd],
    "no_bsm_field_term_count"->Length[zero],
    "bsm_field_count_distribution"->Association@KeyValueMap[ToString[#1]->#2&,distribution],
    "odd_term_examples_inputform"->(ToString[terms[[#]],InputForm]& /@ Take[odd,UpTo[4]]),
    "zero_field_examples_inputform"->(ToString[terms[[#]],InputForm]& /@ Take[zero,UpTo[4]])|>
];
T3ZAExecute[] := Module[{model,built,int,free,ok,status,rep,mfIndices,indexedMF},
  model=T3ModelFromClass[$T3ZAFamily,$T3ZAAlpha];
  If[model===$Failed,Return[<|"status"->"UnsupportedModel"|>]];
  built=BuildT3Lagrangian[model];
  If[!AssociationQ[built],Return[<|"status"->"BuilderFailed"|>]];
  int=T3ZAInspect[built["LInt"]];
  free=T3ZAInspect[built["LFree"]];
  mfIndices=Cases[built["LFree"],
    HoldPattern[Coupling[name_,ids_List,___]] /; name===MF :> ids,
    {0,Infinity}];
  indexedMF=mfIndices=!={} &&
    AllTrue[mfIndices,MatchQ[#,{Index[_,NFlavor]}]&];
  ok=indexedMF && built["Status"]==="Success" && TrueQ[built["FullValidation"]] &&
    int["expanded_term_count"]>0 && free["expanded_term_count"]>0 &&
    int["odd_bsm_term_count"]===0 && free["odd_bsm_term_count"]===0 &&
    int["no_bsm_field_term_count"]===0 && free["no_bsm_field_term_count"]===0;
  status=If[ok,"PassedActualBSMTermParityAndMatchete","TermParityOrBuilderFailed"];
  <|"status"->status,"model"->$T3ZAModel,
    "builder_status"->Lookup[built,"Status","Missing"],
    "matchete_full_check_passed"->TrueQ[built["FullValidation"]],
    "interaction_terms"->int,"free_terms"->free,
    "indexed_mf_in_original_free_lagrangian"->indexedMF,
    "mf_mass_indices_inputform"->(ToString[#,InputForm]& /@ DeleteDuplicates[mfIndices]),
    "extra_yukawas_included"->Lookup[built,"AdditionalYukawasIncluded",-1],
    "extra_quartics_included"->Lookup[built,"AdditionalQuarticsIncluded",-1],
    "extra_quadratics_included"->Lookup[built,"AdditionalQuadraticsIncluded",-1],
    "discrete_z2_explicitly_registered_in_matchete"->False,
    "physical_matching_performed"->False|>
];
$T3ZAResult=CheckAbort[UsingFrontEnd[T3ZAExecute[]],$Aborted];
If[!AssociationQ[$T3ZAResult],$T3ZAResult=<|"status"->"Aborted"|>];
If[!DirectoryQ[DirectoryName[$T3ZAOut]],
  CreateDirectory[DirectoryName[$T3ZAOut],CreateIntermediateDirectories->True]];
Export[$T3ZAOut,Join[<|"format"->"t3-workstream-a-term-parity-v1",
  "model"->$T3ZAModel|>,$T3ZAResult],"RawJSON"];
Print[$T3ZAModel,": ",Lookup[$T3ZAResult,"status","Missing"]];
If[Lookup[$T3ZAResult,"status",""]=!="PassedActualBSMTermParityAndMatchete",
  Exit[6]];
Exit[0];
