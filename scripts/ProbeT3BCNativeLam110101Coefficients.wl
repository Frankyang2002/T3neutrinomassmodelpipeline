(* Narrow independent live RGBeta Lam component test, B_-1/C_-1.
   Uses the installed RGBeta`PackageScope`Lam, NOT reconstructed tensors.
   Checks 16 components per coupling for direct and conjugate orientations.
   Each direct/conjugate field orientation is visited by the outer Do
   iterator; this must be a COMMA before {slot,1,2}, not a semicolon.
   The separate Matchete side has not been run live by this script.
   Usage: wolframscript -file scripts/ProbeT3BCNativeLam110101Coefficients.wl NEW.json
*)
ClearAll["Global`*"];
If[Length[$ScriptCommandLine] =!= 2,
  Print["Usage: ProbeT3BCNativeLam110101Coefficients.wl NEW.json"]; Exit[64]];
outFile = ExpandFileName[$ScriptCommandLine[[2]]];
If[FileExistsQ[outFile], Print["Refusing existing output: ",outFile]; Exit[64]];
root=DirectoryName @ DirectoryName @ ExpandFileName[$InputFileName];
loaded=CheckAbort[Check[Needs["RGBeta`"]; True,$Failed],$Aborted];
If[loaded =!= True, Print["RGBeta not available"];Exit[3]];
Get[FileNameJoin[{root,"RGE","running","rgbeta","T3RGBetaModel.wl"}]];
Get[FileNameJoin[{root,"scripts","RGBetaNativeAccess.wl"}]];

lamName="RGBeta`PackageScope`Lam";
mapName="RGBeta`PackageScope`$fieldIndexMap";
If[Names[lamName] =!= {lamName} || Names[mapName] =!= {mapName},
  Print["Required native package symbol not found"]; Exit[4]];
lamFunction=Symbol[lamName];
If[With[{f=lamFunction},Length[DownValues[f]]] < 1,
  Print["Native Lam has no definition"];Exit[4]];

ClearAll[T3BCFail,T3BCMatrixQ,T3BCMatrixIdentity,T3BCSelectCoupling,
  T3BCDelta,T3BCEpsilon,T3BCComponent,T3BCProbeModel,
  T3BCLookupDiagnostic,T3BCTruncate,T3BCLastDiagnostic];
T3BCLastDiagnostic=<||>;
(* Fail closed but preserve a JSON diagnostic rather than losing the
   nonzero sparse-array entry when the probe's extraction is wrong. *)
T3BCFail[s_String] := (
  Print["ERROR: ",s];
  Throw[<|"status"->"BlockedNativeLamComponentExtraction",
      "error"->s,"diagnostic"->T3BCLastDiagnostic,
      "physical_coupling_mapping_certified"->False,
      "physical_matching_certified"->False|>,"T3BCCoeff"]
);
T3BCTruncate[expr_, n_Integer:1600] := Module[{v=ToString[expr,InputForm]},
  StringTake[v,Min[n,StringLength[v]]]
];
T3BCLookupDiagnostic[array_SparseArray,fieldMap_,positions_,entry_,model_String] := Module[
  {rules,mentions,targets={"lambdaT3GeneralQ110101I1",
    "lambdaHHdagS1barS2barCross"}},
  rules=Most@ArrayRules[array];
  mentions=Select[rules,Function[rule,
    AnyTrue[targets,Function[target,StringContainsQ[
      T3BCTruncate[Last[rule],100000],target]]]]];
  (* Do not silently infer field positions from the first matching term. *)
  <|"model"->model,
    "sparse_dimensions"->Dimensions[array],
    "scalar_field_map_inputform"->T3BCTruncate[fieldMap,3000],
    "requested_positions"->positions,
    "requested_entry_inputform"->T3BCTruncate[entry,4000],
    "matrix_symbol_names"->Names["*`Matrix"],
    "coupling_bearing_rule_count"->Length[mentions],
    "coupling_bearing_sparse_rules_excerpt"->(T3BCTruncate[#,2400]&/@Take[mentions,UpTo[6]]),
    "extraction_success"->False|>
];
(* The installed v1.2.0 native Lam *evaluates* its Matrix[coupling][]
   declaration into a bare coupling symbol in the sparse array.  The live
   2026-10-31 diagnostic actually printed
       (lambdaFirst*del*eps + lambdaSecond*del*eps)/4.
   Hence extracting only Matrix[...] heads is incorrect.  Use exact, linear
   symbolic coefficients from the evaluated Lam entry.  This is not a
   reconstructed interaction: the entry is Extract[realNativeLam,position].

   For h.c. RGBeta may retain Bar[coupling] or Conjugate[coupling], depending
   on installed package simplifications.  Accept *exactly one* of these
   representations, recording which one occurred.  Never allow an absent
   coefficient to pass as a zero tensor. *)
T3BCSelectCoupling[expression_, target_Symbol, hc_:False] := Module[
  {forms, coefficients, which, coefficient, form},
  forms=DeleteDuplicates[If[TrueQ[hc],
    {Bar[target],Conjugate[target]}, {target}]];
  coefficients=Coefficient[Expand[expression],#,1]&/@forms;
  which=Flatten@Position[
    ((!TrueQ[FullSimplify[# == 0]])&/@coefficients), True];
  If[Length[which] =!= 1,
    T3BCFail["Cannot uniquely extract "<>SymbolName[target]<>
      If[TrueQ[hc]," conjugate"," direct"]<>
      " from native evaluated Lam: nonzero coefficient forms="<>
      ToString[Length[which]]<>
      "; candidate forms="<>ToString[forms,InputForm]]];
  form=forms[[First[which]]];
  coefficient=coefficients[[First[which]]];
  If[!FreeQ[coefficient,target],
    T3BCFail["Nonlinear self-dependence of native coupling "<>SymbolName[target]]];
  {coefficient,form}
];
T3BCDelta[i_Integer,j_Integer] := KroneckerDelta[i,j];
T3BCEpsilon[i_Integer,j_Integer] := Signature[{i,j}];
T3BCComponent[expression_,indices_List] := Module[{v},
  (* Convert RGBeta's del/eps into inert symbols BEFORE setting any numeric
     indices: native del[rep, a, a] would otherwise evaluate to Dim[rep]. *)
  v=expression /. {
    HoldPattern[del[SU2L[fund],a_,b_]] :> T3BCDelta[a,b],
    HoldPattern[eps[SU2L[fund],a_,b_]] :> T3BCEpsilon[a,b]
  };
  v=v /. Thread[{t3a,t3b,t3c,t3d}->indices];
  v=FullSimplify[v];
  If[!NumericQ[v] || !TrueQ[FullSimplify[Im[v] == 0]] ||
     !FreeQ[v,_Missing|$Failed|$Aborted|Indeterminate|del|eps|Matrix],
    T3BCFail["Unresolved native component "<>ToString[indices,InputForm]<>" = "<>ToString[v,InputForm]]];
  v
];
T3BCProbeModel[model_String,dF_Integer] := Module[
  {built,extension,symbols,records,fieldMap,fieldOrders,positions,
   expr,arrays,array,idx,terms,checks={},failures={},fieldOrder,entry,
   coeff,inv,observed,expect,comp,conjFlag,slot,k,chosen,term,
   selectedTerms,orientations={},maxError=0.,count=0,ok},
  T3BCLastDiagnostic=<||>;
  Print["Checking live Lam source coefficients for ",model," ..."];
  built=CheckAbort[Check[T3RGBetaBuild[2,2,dF,-1],$Failed],$Aborted];
  If[!AssociationQ[built],T3BCFail[model<>" build failed"]];
  extension=FileNameJoin[{root,"RGE","general","generated",model<>"_extra_quartics.wl"}];
  If[!FileExistsQ[extension],T3BCFail[model<>" missing extension"]];
  If[CheckAbort[Check[Get[extension];True,$Failed],$Aborted] =!= True,
    T3BCFail[model<>" extension failed"]];
  symbols=Symbol/@{"lambdaT3GeneralQ110101I1","lambdaHHdagS1barS2barCross"};
  If[!And@@(KeyExistsQ[$quartics,#]&/@symbols),
    T3BCFail[model<>" missing registered couplings"]];
  records=T3Stage22NativeRecord/@symbols;
  If[!And@@(StringReplace[Lookup[#,"fields_inputform",""],WhitespaceCharacter..->""]===
      "{H,Bar[H],Bar[S1],Bar[S2]}"&/@records),
    T3BCFail[model<>" wrong source field order"]];
  If[!And@@(StringMatchQ[StringTrim[Lookup[#,"coupling_factor_inputform",""]],
      "6*Matrix["~~__~~"][##1] &"]&/@records),
    T3BCFail[model<>" nonstandard native normalization"]];
  If[Lookup[records,"self_conjugate"] =!= {False,False},
    T3BCFail[model<>" nonstandard self-conjugation"]];
  If[!And@@(Length[ToExpression[Lookup[#,"unique_arrangements_inputform"]]]==24&/@records),
    T3BCFail[model<>" expected 24 distinct permutations"]];
  fieldMap=Symbol[mapName]["Scalars"];
  fieldOrders={{H,Bar[H],Bar[S1],Bar[S2]},
               {Bar[H],H,S1,S2}};
  positions=(fieldMap/@#)&/@fieldOrders;
  If[!And@@(AllTrue[#,IntegerQ[#]&&#>0&]&/@positions),
    T3BCFail[model<>" invalid native scalar index positions"]];
  expr=CheckAbort[Check[lamFunction[t3a,t3b,t3c,t3d],$Failed],$Aborted];
  arrays=Cases[expr,_SparseArray,Infinity];
  If[Length[arrays] =!= 1,T3BCFail[model<>" expected exactly one Lam sparse tensor"]];
  array=First[arrays];
  (* For four distinct fields, Lam has one of the 24 arrangements in each
     ordered field slot.  RGBeta's internal factor 6 gives 6/24=1/4.
     The action contains 24 orderings and is normalized by 1/4!, with the
     additional complex-to-real field conversion treated separately. *)
  Do[
    conjFlag=(slot==2);
    entry=Extract[array,positions[[slot]]];
    T3BCLastDiagnostic=T3BCLookupDiagnostic[
      array,fieldMap,positions[[slot]],entry,model];
    T3BCLastDiagnostic["orientation"]=If[conjFlag,"conjugate","direct"];
    If[!StringContainsQ[T3BCLastDiagnostic["requested_entry_inputform"],
          "lambdaT3GeneralQ110101I1"] &&
       !StringContainsQ[T3BCLastDiagnostic["requested_entry_inputform"],
          "lambdaHHdagS1barS2barCross"],
      Print[model,": requested ",T3BCLastDiagnostic["orientation"],
        " entry at ",positions[[slot]]," contains no 110101 name;",
        " recording coupling-bearing sparse locations"]
    ];
    selectedTerms={};
    Do[
      chosen=T3BCSelectCoupling[entry,symbols[[k]],conjFlag];
      coeff=chosen[[1]];
      term=chosen[[2]];
      AppendTo[selectedTerms,term*coeff];
      AppendTo[orientations,<|
        "orientation"->If[conjFlag,"conjugate","direct"],
        "coupling"->SymbolName[symbols[[k]]],
        "evaluated_coupling_form"->ToString[term,InputForm]|>];
      inv=records[[k]]["_native_invariant"][t3a,t3b,t3c,t3d];
      Do[
        observed=T3BCComponent[coeff,idx];
        expect=T3BCComponent[inv,idx]/4;
        comp=FullSimplify[observed-expect];
        count++;
        If[!TrueQ[comp==0],
          AppendTo[failures,<|"orientation"->If[conjFlag,"conjugate","direct"],
             "coupling"->SymbolName[symbols[[k]]],"indices"->idx,
             "observed"->ToString[observed,InputForm],
             "expected"->ToString[expect,InputForm]|>]];
        maxError=Max[maxError,Abs[N[comp]]],
        {idx,Tuples[{1,2},4]}],
      {k,1,2}];
    (* A source drift/extra coupling in the same ordered channel must not
       silently disappear when projecting only the two expected symbols. *)
    If[!TrueQ[FullSimplify[Expand[entry-Total[selectedTerms]] == 0]],
      T3BCFail["Native Lam selected entry has extra/nonlinear coupling terms for "<>
        model<>"/"<>If[conjFlag,"conjugate","direct"]]],
    {slot,1,2}];
  ok=(count==64 && failures=={});
  Print[model,": evaluated ",count," direct/conjugate native Lam components; mismatches=",Length[failures],
        "; max residual=",ScientificForm[maxError,3]];
  <|"model"->model,
    "status"->If[ok,"LiveLam110101SourceComponentCoefficientsPassed",
                         "LiveLam110101SourceComponentCoefficientsMismatch"],
    "live_native_lam_called"->True,"symbol"->lamName,
    "source_symbols"->SymbolName/@symbols,
    "source_scalar_field_positions"->positions,
    "component_tests"->count,"expected_component_tests"->64,
    "expected_lam_over_registered_invariant"->"1/4",
    "native_evaluated_coupling_forms"->orientations,
    "coefficient_extraction"->"symbolic Coefficient on evaluated native Lam",
    "failed_component_count"->Length[failures],
    "first_mismatches"->Take[failures,UpTo[8]],
    "maximum_absolute_component_residual"->N[maxError],
    "registered_complex_quartic_normalization"->"6*Matrix[...]",
    "unique_arrangements_per_coupling"->24,
    "matchete_runtime_evaluated"->False,
    "physical_coupling_mapping_certified"->False,
    "physical_matching_certified"->False|>
];

result=CheckAbort[Check[
  <|"format"->"t3-bc-native-rgbeta-lam-110101-component-coefficients-v1",
    "status"->"LiveLam110101SourceCoefficientComponentAudit",
    "models"->Association@Table[
       spec[[1]]->Catch[T3BCProbeModel[spec[[1]],spec[[2]]],"T3BCCoeff"],
       {spec,{{"B_-1",1},{"C_-1",3}}}],
    "source"->"installed RGBeta`PackageScope`Lam",
    "action_level_matchete_independently_evaluated"->False,
    "physical_coupling_mapping_certified"->False,
    "physical_matching_certified"->False|>,$Failed],$Aborted];
If[!AssociationQ[result] || Length[Keys[result["models"]]] =!= 2,
  Print["Native Lam component coefficient audit failed to construct report"];Exit[5]];
If[!DirectoryQ[DirectoryName[outFile]],CreateDirectory[DirectoryName[outFile],CreateIntermediateDirectories->True]];
If[Check[Export[outFile,result,"RawJSON"],$Failed] === $Failed || !FileExistsQ[outFile],
   Print["Failed to write component audit report"];Exit[6]];
Print["Live native Lam coefficient report: ",outFile];
If[!AllTrue[Values[result["models"]],
      Lookup[#,"failed_component_count",-1]===0 &&
      Lookup[#,"status",""]==="LiveLam110101SourceComponentCoefficientsPassed"&],
  Print["Native Lam coefficient audit BLOCKED/FAILED; diagnostic JSON saved: ",outFile];
  Do[
    If[Lookup[result["models",label],"status",""]===
         "BlockedNativeLamComponentExtraction",
       Print[label,": ",Lookup[result["models",label],"error","unknown"],
         "; requested entry: ",
         Lookup[Lookup[result["models",label],"diagnostic",<||>],
                "requested_entry_inputform","missing"],
         "; bearing rules: ",
         Lookup[Lookup[result["models",label],"diagnostic",<||>],
                "coupling_bearing_rule_count",-1]]
    ],{label,{"B_-1","C_-1"}}];
  Exit[7]];
Print["Native Lam source coefficient test PASSED for both models."];
Print["Physical Matchete/RGBeta action and matching: NOT CERTIFIED"];
