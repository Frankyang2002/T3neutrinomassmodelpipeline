(* Stage 24 read-only audit of the FOUR remaining B/C reused quartic
   GroupInvariant functions, independently evaluated in live RGBeta.

   wolframscript -file scripts/ExportGeneralT3BCQuarticNativeCovariance.wl \
     output/general_t3_stage22/stage24_bc_wolfram_native_covariance.json

   Assumes local Stage22 110101 repair has been applied. Does not alter
   production, write beta functions, or certify physical coupling maps.
*)
ClearAll["Global`*"];
If[Length[$ScriptCommandLine] =!= 2,
  Print["Usage: ExportGeneralT3BCQuarticNativeCovariance.wl OUTPUT.json"];Exit[64]];
outPath=ExpandFileName[$ScriptCommandLine[[2]]];
root=DirectoryName@DirectoryName@ExpandFileName[$InputFileName];
If[FileExistsQ[outPath],DeleteFile[outPath]];
load=CheckAbort[Check[Needs["RGBeta`"];True,$Failed],$Aborted];
If[load=!=True,Print["RGBeta package unavailable"];Exit[3]];
Get[FileNameJoin[{root,"RGE","running","rgbeta","T3RGBetaModel.wl"}]];
Get[FileNameJoin[{root,"scripts","RGBetaNativeAccess.wl"}]];
T3S24Fail[msg_] := (Print["Stage24: ",msg];Throw[$Failed,"T3S24"]);

(* RGBeta's del[rep,a_,a_] equals Dim[rep] for a repeated DUMMY
   index. Replace tensor heads BEFORE assigning numeric component labels,
   so a component KroneckerDelta[i,i] correctly evaluates to one. *)
T3S24Delta[a_Integer,b_Integer]:=KroneckerDelta[a,b];
T3S24Epsilon[a_Integer,b_Integer]:=Signature[{a,b}];
T3S24Component[expr_,inds_List]:=Module[{v},
  v=expr/.{
    HoldPattern[del[SU2L[fund],a_,b_]]:>T3S24Delta[a,b],
    HoldPattern[eps[SU2L[fund],a_,b_]]:>T3S24Epsilon[a,b]
  };
  v=FullSimplify[v/.Thread[{t3i1,t3i2,t3i3,t3i4}->inds]];
  If[!FreeQ[v,del|eps|delS2|_Missing|Indeterminate|$Failed]||
       !NumericQ[v]||!TrueQ[FullSimplify[Im[v]==0]],
    T3S24Fail["Unresolved SU2 component: "<>ToString[v,InputForm]]];
  v
];
T3S24Tensor[expr_]:=Table[T3S24Component[expr,{i,j,k,l}],
  {i,1,2},{j,1,2},{k,1,2},{l,1,2}];

T3S24Generators[fields_List,n_Integer]:=Module[{g=PauliMatrix[n]/2},
  If[Length[fields]=!=4,T3S24Fail["Wrong field rank"]];
  (If[StringStartsQ[#,"Bar["],-Conjugate[g],g]&) /@ fields
];
T3S24Audit[expr_,fields_List]:=Module[{tensor,res,nz,gens,idx},
  tensor=T3S24Tensor[expr];
  res=Table[
    gens=T3S24Generators[fields,n];
    Flatten@Table[
      idx={i,j,k,l};
      FullSimplify[Sum[
        Sum[gens[[axis,b,idx[[axis]]]] *
          Extract[tensor,ReplacePart[idx,axis->b]],{b,1,2}],
        {axis,1,4}]],
      {i,1,2},{j,1,2},{k,1,2},{l,1,2}],
    {n,1,3}];
  nz=(Count[#,Except[0]]&) /@ res;
  <|"su2_covariant"->(nz==={0,0,0}),
    "nonzero_residual_components_by_generator"->nz,
    "component_tensor_flat_inputform"->(ToString[#,InputForm]& /@ Flatten[tensor])|>
];

symbols={"lambdaS1barS2S2bar2","lambdaS1bar2S2bar2",
 "lambdaS1S1bar2S2bar","lambdaHHdagS1barS1bar"};
fieldSets=<|
 "lambdaS1barS2S2bar2"->{"Bar[S1]","S2","Bar[S2]","Bar[S2]"},
 "lambdaS1bar2S2bar2"->{"Bar[S1]","Bar[S1]","Bar[S2]","Bar[S2]"},
 "lambdaS1S1bar2S2bar"->{"S1","Bar[S1]","Bar[S1]","Bar[S2]"},
 "lambdaHHdagS1barS1bar"->{"Bar[H]","Bar[H]","Bar[S1]","Bar[S1]"}
|>;
proposed=<|
 "lambdaS1barS2S2bar2"->((eps[SU2L@fund,#1,#4]*del[SU2L@fund,#2,#3]+
                              eps[SU2L@fund,#1,#3]*del[SU2L@fund,#2,#4])/2 &),
 "lambdaS1bar2S2bar2"->((eps[SU2L@fund,#1,#4]*eps[SU2L@fund,#2,#3]+
                           eps[SU2L@fund,#1,#3]*eps[SU2L@fund,#2,#4])/2 &),
 "lambdaS1S1bar2S2bar"->((del[SU2L@fund,#1,#3]*eps[SU2L@fund,#2,#4]+
                           del[SU2L@fund,#1,#2]*eps[SU2L@fund,#3,#4])/2 &),
 "lambdaHHdagS1barS1bar"->((eps[SU2L@fund,#1,#4]*eps[SU2L@fund,#2,#3]+
                             eps[SU2L@fund,#1,#3]*eps[SU2L@fund,#2,#4])/2 &)
|>;

T3S24Run[]:=Module[{result=<||>,model,spec,build,ext,name,record,fs,
                    original,alternative,oldAudit,newAudit,rows},
  Do[
    model=spec[[1]];
    Print["Stage24 live native RGBeta B/C covariance: ",model];
    build=CheckAbort[Check[Apply[T3RGBetaBuild,spec[[2]]],$Failed],$Aborted];
    If[!AssociationQ[build],T3S24Fail["Build failed: "<>model]];
    T3GeneralQuarticAddedSymbols={};
    ext=FileNameJoin[{root,"RGE","general","generated",model<>"_extra_quartics.wl"}];
    If[CheckAbort[Check[Get[ext];True,$Failed],$Aborted]=!=True,
      T3S24Fail["Generated quartic extension failed: "<>model]];
    rows=<||>;
    Do[
      record=T3Stage22NativeRecord[Symbol[name]];
      fs=fieldSets[name];
      If[StringReplace[record["fields_inputform"],WhitespaceCharacter..->""]=!=
         "{"<>StringRiffle[fs,","]<>"}",
        T3S24Fail["Unexpected registered field order: "<>model<>"/"<>name]];
      original=record["_native_invariant"][t3i1,t3i2,t3i3,t3i4];
      alternative=proposed[name][t3i1,t3i2,t3i3,t3i4];
      oldAudit=T3S24Audit[original,fs];
      newAudit=T3S24Audit[alternative,fs];
      If[TrueQ[oldAudit["su2_covariant"]]||!TrueQ[newAudit["su2_covariant"]],
        T3S24Fail["Native/proposed status changed: "<>model<>"/"<>name]];
      AssociateTo[rows,name-><|
        "native_expression_inputform"->ToString[original,InputForm],
        "proposed_expression_inputform"->ToString[alternative,InputForm],
        "native"->oldAudit,"proposed"->newAudit
      |>],
      {name,symbols}];
    AssociateTo[result,model->rows],
    {spec,{{"B_-1",{2,2,1,-1}},{"C_-1",{2,2,3,-1}}}}];
  <|"format"->"stage24-live-rgbeta-bc-other-quartic-covariance-v1",
    "source"->"RGBeta native registered GroupInvariant functions",
    "convention"->"fundamental:U, Bar fields:Conjugate[U]",
    "models"->result,
    "native_non_covariant_slots"->8,
    "proposed_covariant_slots"->8,
    "physical_coupling_mapping_certified"->False,
    "production_modified"->False|>
];
result=CheckAbort[Check[Catch[T3S24Run[],"T3S24"],$Failed],$Aborted];
If[!AssociationQ[result]||Length[Keys[Lookup[result,"models",<||>]]]=!=2,
  Print["Stage24 live native covariance FAILED"];Exit[5]];
If[!DirectoryQ[DirectoryName[outPath]],
  CreateDirectory[DirectoryName[outPath],CreateIntermediateDirectories->True]];
write=Quiet@Check[Export[outPath,result,"RawJSON"],$Failed];
If[write===$Failed||!FileExistsQ[outPath],
  Print["Stage24 live JSON export FAILED"];Exit[6]];
Print["Stage24 live native B/C covariance evidence: ",outPath];
Print["This is NOT a physical Matchete-RGBeta coupling certificate"];
