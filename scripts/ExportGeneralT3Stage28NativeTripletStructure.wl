(* Stage 28: live RGBeta v1.2.0 SU(2) S2 structural/conjugation evidence.
   This deliberately DOES NOT pretend tGen[SU2L[S2],A,i,j] is a numeric
   component matrix.  S2 is an abstract symmetric-square representation.
   Run from root:
     wolframscript -file scripts/ExportGeneralT3Stage28NativeTripletStructure.wl \
       output/general_t3_stage22/stage28_live_rgbeta_triplet_structure.json
*)
ClearAll["Global`*"];
If[Length[$ScriptCommandLine] =!= 2,
  Print["Usage: ExportGeneralT3Stage28NativeTripletStructure.wl OUTPUT.json"];Exit[64]];
outPath = ExpandFileName[$ScriptCommandLine[[2]]];
root = DirectoryName @ DirectoryName @ ExpandFileName[$InputFileName];
If[FileExistsQ[outPath],DeleteFile[outPath]];
reg = Quiet@Check[Import[FileNameJoin[{root,"RGE","general","generated",
  "general_rgbeta_quartic_registry.json"}],"RawJSON"],$Failed];
If[!AssociationQ[reg] || Length[Keys[Lookup[reg,"models",<||>]]]=!=16,
  Print["Stage28: Incomplete model registry"];Exit[2]];
loaded=CheckAbort[Check[Needs["RGBeta`"];True,$Failed],$Aborted];
If[loaded=!=True,Print["Stage28: RGBeta unavailable"];Exit[3]];
Get[FileNameJoin[{root,"RGE","running","rgbeta","T3RGBetaModel.wl"}]];
Get[FileNameJoin[{root,"scripts","Stage22RGBetaNativeAccess.wl"}]];

T3S28Fail[msg_String]:=(Print["Stage28 ERROR: ",msg];Throw[$Failed,"T3Stage28"]);
T3S28Resolve[suffix_String]:=Module[{names},
  names=DeleteDuplicates@Flatten[Names/@{"RGBeta`"<>suffix,
    "RGBeta`Private`"<>suffix,"RGBeta`*`"<>suffix,"*`"<>suffix}];
  names=Select[names,StringStartsQ[#,"RGBeta`"]&&
    Last[StringSplit[#,"`"]]===suffix&];
  If[Length[names]=!=1,T3S28Fail["Need one RGBeta symbol "<>suffix<>
    "; found: "<>ToString[names,InputForm]]];
  Symbol[First[names]]
];
T3S28Models={"A_+0","A_-2","D_+0","D_-2","E_+0","E_-2"};
T3S28Inspect[model_String]:=Module[
  {dims,meta,build,triplets={},rep,dim,cas,tr,gen,bar,orth,
   expected,im,metric,scal,barScal,pos,barPos},
  Print["Stage28 native SU2 triplet structure: ",model];
  meta=reg["models",model];dims=meta["dimensions"];
  If[!ListQ[dims]||Length[dims]=!=3,T3S28Fail["Bad dimensions "<>model]];
  expected=Switch[model,"A_+0"|"A_-2",{1,3},
    "D_+0"|"D_-2",{3,1},"E_+0"|"E_-2",{3,3},_,{}];
  If[Take[dims,2]=!=expected,T3S28Fail["Triplet species dimension mismatch "<>model]];
  build=CheckAbort[Check[Apply[T3RGBetaBuild,Join[dims,{meta["alpha"]}]],$Failed],$Aborted];
  If[!AssociationQ[build],T3S28Fail["RGBeta model build failed "<>model]];
  rep=SU2L[S2];dim=Dim[rep];cas=Casimir2[rep];tr=TraceNormalization[rep];
  If[!TrueQ[dim===3&&cas===2&&tr===2],
    T3S28Fail["Native SU2 symmetric-square dimension/Casimir/trace mismatch "<>model<>
      ": "<>ToString[{dim,cas,tr},InputForm]]];
  gen=tGen[rep,t28A,t28x,t28y];
  If[Head[gen]=!=tGen,T3S28Fail["Unexpected numeric or rewritten S2 generator "<>model]];
  bar=tGen[Bar[rep],t28A,t28x,t28y];
  If[!TrueQ[bar===-tGen[rep,t28A,t28y,t28x]],
    T3S28Fail["Native conjugate generator convention changed "<>model]];
  orth=delS2[SU2L,t28r,t28i,t28j]*delS2[SU2L,t28s,t28i,t28j];
  If[!TrueQ[orth===del[rep,t28r,t28s]],
    T3S28Fail["Native S2 symmetric-doublet contraction changed "<>model<>
      ": "<>ToString[orth,InputForm]]];
  im=T3S28Resolve["$fieldIndexMap"];
  metric=T3S28Resolve["$scalarContraction"];
  Do[
    If[dims[[k]]===3,
      scal=If[k===1,S1,S2];barScal=Bar[scal];
      If[!KeyExistsQ[$scalars,scal]||
         !MemberQ[T3Stage22NativeGet[$scalars[scal],"GaugeRep"],rep],
        T3S28Fail["Triplet field does not use SU2L[S2] "<>model]];
      If[!KeyExistsQ[im["Scalars"],scal]||!KeyExistsQ[im["Scalars"],barScal],
        T3S28Fail["Missing complex triplet/conjugate index "<>model]];
      pos=im["Scalars",scal];barPos=im["Scalars",barScal];
      If[!TrueQ[metric[[pos,barPos]]===2]||
         !TrueQ[metric[[barPos,pos]]===2],
         T3S28Fail["Native scalar conjugate metric factor changed "<>model]];
      AppendTo[triplets,SymbolName[scal]]]
    ,{k,2}];
  If[Length[triplets]=!=Count[expected,3],T3S28Fail["Triplet count differs "<>model]];
  <|"model"->model,"scalar_su2_dimensions"->expected,
    "triplet_fields"->triplets,
    "native_S2_representation"->"SU2L[S2]",
    "dimension"->dim,"quadratic_casimir"->cas,"trace_normalisation"->tr,
    "native_generator_inputform"->ToString[gen,InputForm],
    "native_generator_is_symbolic"->True,
    "native_conjugate_generator_minus_transpose"->True,
    "native_delS2_orthogonality"->True,
    "native_complex_conjugate_metric_two"->True,
    "native_numerical_generator_components_exported"->False|>
];
T3S28Run[]:=Module[{rows},
  rows=T3S28Inspect/@T3S28Models;
  <|"format"->"stage28-live-rgbeta-s2-structure-v1",
    "source"->"RGBeta runtime group representation and scalar index internals",
    "models_checked"->Length[rows],"models"->rows,
    "all_group_representation_checks_passed"->True,
    "all_native_barred_generator_rules_passed"->True,
    "all_native_delS2_contractions_passed"->True,
    "native_numerical_triplet_generators_available"->False,
    "physical_cross_package_component_basis_certified"->False,
    "physical_coupling_mapping_certified"->False,
    "production_modified"->False|>
];
result=CheckAbort[Check[Catch[T3S28Run[],"T3Stage28"],$Failed],$Aborted];
If[!AssociationQ[result]||Lookup[result,"models_checked",0]=!=6,
  Print["Stage28 native structure export FAILED"];Exit[5]];
If[!DirectoryQ[DirectoryName[outPath]],
  CreateDirectory[DirectoryName[outPath],CreateIntermediateDirectories->True]];
write=Quiet@Check[Export[outPath,result,"RawJSON"],$Failed];
If[write===$Failed||!FileExistsQ[outPath],
  Print["Stage28 JSON export FAILED"];Exit[6]];
Print["Stage28 native symbolic triplet structure: 6/6 models"]; 
Print["RGBeta S2 numerical generator components: NOT SPECIFIED by abstract tGen"]; 
Print["Physical cross-package component basis and coupling map: UNCERTIFIED"]; 
Print["Report: ",outPath];
