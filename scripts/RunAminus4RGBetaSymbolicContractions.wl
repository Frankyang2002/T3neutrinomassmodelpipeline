(* A(-4) EFT1: read-only native RGBeta group-index contraction witness.
   wolframscript -file scripts/RunAminus4RGBetaSymbolicContractions.wl OUTPUT.json

   Do not supply user-made tGen matrices, CG components, index-reduction
   rules, or change package definitions.  This tests only which contractions
   the installed RGBeta runtime actually evaluates.  It does not certify
   a native-to-Matchete component-axis or physical-coupling map.
*)
ClearAll["Global`*"];
If[Length[Rest[$ScriptCommandLine]] =!= 1,
  Print["Usage: OUTPUT.json"]; Exit[2]];
$ACOut = ExpandFileName[Last[$ScriptCommandLine]];
If[FileExistsQ[$ACOut] || FileExistsQ[$ACOut <> ".partial"],
  Print["BLOCKED: output already exists: ", $ACOut]; Exit[2]];
$ACFlags = {"T3_EXPERIMENTAL_GENERAL_FERMIONS",
  "T3_EXPERIMENTAL_GENERAL_QUADRATICS",
  "T3_EXPERIMENTAL_GENERAL_QUARTICS",
  "T3_EXPERIMENTAL_FLAVOR_MASSES"};
If[AnyTrue[$ACFlags, Environment[#] === "1" &],
  Print["BLOCKED: unrestricted T3 feature flag enabled"]; Exit[3]];
$ACRoot = DirectoryName @ DirectoryName @ ExpandFileName[$InputFileName];
$ACModel = FileNameJoin[{$ACRoot,"RGE","running","rgbeta","T3RGBetaModel.wl"}];
If[!FileExistsQ[$ACModel], Print["BLOCKED: model source missing"]; Exit[3]];
$ACSHA = Quiet @ Check[IntegerString[FileHash[$ACModel,"SHA256"],16,64],$Failed];
If[!StringQ[$ACSHA] || StringLength[$ACSHA] =!= 64,
  Print["BLOCKED: model source hash unavailable"]; Exit[3]];
If[CheckAbort[Check[Needs["RGBeta`"]; True, $Failed], $Aborted] =!= True,
  Print["BLOCKED: installed RGBeta unavailable"]; Exit[4]];
Get[$ACModel];
$ACBuild = CheckAbort[Quiet @ Check[
  T3RGBetaBuildEFT1[1,3,2,-4], $Failed], $Aborted];
If[!AssociationQ[$ACBuild], Print["BLOCKED: A(-4) EFT1 build failed"]; Exit[4]];

ClearAll[ACString, ACProbe, ACRegistryInvariant, ACName,
  ACRegistryRec, ACNativeDefinitions];
ACString[x_] := ToString[x,InputForm];
ACName[x_] := Last @ StringSplit[ACString[x],"`"];
ACRegistryRec[sym_] := If[AssociationQ[$quartics],
  Lookup[$quartics,sym,<||>],<||>];
ACRegistryInvariant[sym_] := Module[{r=ACRegistryRec[sym], keys},
  If[!AssociationQ[r], Return[$Failed]];
  keys=Select[Keys[r],ACName[#]==="Invariant" &];
  If[Length[keys]=!=1,Return[$Failed]];
  Lookup[r,First[keys]]
];

(* The actual products are passed as HoldComplete, not substituted from an
   assumed SU(2) representation.  RGBeta's own up-values and FullSimplify
   decide whether each symbolic contraction reduces.  Unreduced values are
   recorded as unresolved, never labelled failed physical identities. *)
ACProbe[id_String,expr_HoldComplete,expected_HoldComplete] := Module[
  {v,e,reduced,exact=False,st},
  v = CheckAbort[TimeConstrained[Quiet@Check[
    ReleaseHold[expr],$Failed],12,$TimedOut],$Aborted];
  e = ReleaseHold[expected];
  reduced = If[MemberQ[{$Failed,$TimedOut,$Aborted},v],v,
    CheckAbort[TimeConstrained[Quiet@Check[
      FullSimplify[Expand[v-e]],$Failed],12,$TimedOut],$Aborted]];
  exact = TrueQ[reduced===0];
  st = Which[exact,"EVALUATED_MATCH",v===$TimedOut || reduced===$TimedOut,
    "TIMED_OUT",v===$Failed || reduced===$Failed ||
    v===$Aborted || reduced===$Aborted,"EVALUATION_BLOCKED",
    True,"SYMBOLICALLY_UNRESOLVED_OR_DIFFERENT"];
  <|"id"->id,"status"->st,
    "source_expression_inputform"->ACString[expr],
    "native_evaluated_inputform"->ACString[v],
    "target_inputform"->ACString[e],
    "native_difference_inputform"->ACString[reduced],
    "native_exact_target_match"->TrueQ[exact]|>
];
$ACGroupProbes = {
  ACProbe["fund_casimir_3_over_4",
    HoldComplete[tGen[SU2L[fund],acA,aci,acj]*
      tGen[SU2L[fund],acA,acj,ack]],
    HoldComplete[(3/4)*del[SU2L[fund],aci,ack]]],
  ACProbe["triplet_casimir_2",
    HoldComplete[tGen[SU2L[S2],acA,aci,acj]*
      tGen[SU2L[S2],acA,acj,ack]],
    HoldComplete[2*del[SU2L[S2],aci,ack]]],
  ACProbe["triplet_fully_contracted_generator_norm_6",
    HoldComplete[tGen[SU2L[S2],acA,aci,acj]*
      tGen[SU2L[S2],acA,acj,aci]],
    HoldComplete[6]],
  ACProbe["fund_fully_contracted_generator_norm_3_over_2",
    HoldComplete[tGen[SU2L[fund],acA,aci,acj]*
      tGen[SU2L[fund],acA,acj,aci]],
    HoldComplete[3/2]]
};
(* Diagnostic quartic group contractions: capture the exact native
   invariant function with free symbolic indices.  No assumed values. *)
$ACQuarticNames = {"lambdaS2","lambdaS2Adj","lambdaH2",
  "lambdaH2Adj","lambdaT3"};
$ACQuarticRows = Table[Module[{sym,fn,expression},
  sym=ToExpression[name];fn=ACRegistryInvariant[sym];
  expression=If[fn===$Failed,$Failed,
    CheckAbort[TimeConstrained[Quiet@Check[
      fn[aci,acj,ack,acl],$Failed],12,$TimedOut],$Aborted]];
  <|"name"->name,"field_order_inputform"->"recorded in previous witness",
    "native_invariant_inputform"->ACString[expression],
    "captured"->!MemberQ[{$Failed,$TimedOut,$Aborted},expression]|>
],{name,$ACQuarticNames}];

(* Test dummy-index freshness at the native registry-function boundary.
   The production source uses the literal symbol A for each SU(2) adjoint
   contraction.  Reusing it in two independently instantiated quartic
   kernels could matter if the RGBeta projector does not alpha-rename its
   contracted indices.  This probe establishes RAW reuse only; it does NOT
   infer whether the projector renames them later. *)
ACDoubleCall[id_String,sym1_,sym2_] := Module[{f1,f2,x,y,shared},
  f1=ACRegistryInvariant[sym1];f2=ACRegistryInvariant[sym2];
  x=If[f1===$Failed,$Failed,f1[aci,acj,ack,acl]];
  y=If[f2===$Failed,$Failed,f2[acp,acq,acr,acs]];
  shared=TrueQ[!FreeQ[x,A] && !FreeQ[y,A]];
  <|"id"->id,"first_call_inputform"->ACString[x],
    "second_call_inputform"->ACString[y],
    "same_literal_A_present_in_both_raw_calls"->shared|>
];
$ACDoubleCalls = {
  ACDoubleCall["S2Adj_times_S2Adj",lambdaS2Adj,lambdaS2Adj],
  ACDoubleCall["H2Adj_times_H2Adj",lambdaH2Adj,lambdaH2Adj],
  ACDoubleCall["S2Adj_times_H2Adj",lambdaS2Adj,lambdaH2Adj],
  ACDoubleCall["S2Singlet_times_S2Singlet_control",lambdaS2,lambdaS2]
};

(* Inspect package-defined reductions/projector symbols.  This is strictly
   observational; no internal symbol is assigned a new rule or value. *)
$ACAllNames = DeleteDuplicates @ Join[Names["RGBeta`PackageScope`*"],
  Names["RGBeta`FieldsAndCouplings`PackagePrivate`*"],
  Names["RGBeta`*`*Project*"], Names["RGBeta`*`*Quartic*"]];
$ACProjectorNames = Select[$ACAllNames,
  StringContainsQ[ToLowerCase[#],
    "project"|"quartic"|"potential"|"fierz"|"simplifygroup"|"contract"] &];
$ACProjectorNames = Take[Sort[$ACProjectorNames],UpTo[80]];
ACNativeDefinitions[name_String] := Module[{defs},
  defs = CheckAbort[TimeConstrained[Quiet@Check[
    With[{s=ToExpression[name]},
      ACString[{Attributes[s],DownValues[s],OwnValues[s]}]],
      "unavailable"],8,"timed_out"],"aborted"];
  <|"symbol"->name,"definition_prefix_inputform"->
      StringTake[defs,UpTo[18000]],
    "full_definition_length"->StringLength[defs],
    "truncated"->TrueQ[StringLength[defs]>18000]|>
];
$ACSymbolRows = ACNativeDefinitions /@ Take[$ACProjectorNames,UpTo[25]];
$ACProbeOK = Length[$ACGroupProbes]===4 &&
  Length[$ACQuarticRows]===5 &&
  AllTrue[$ACQuarticRows,TrueQ[Lookup[#,"captured",False]] &];
$ACStatus = If[$ACProbeOK,
  "NATIVE_SYMBOLIC_CONTRACTION_DIAGNOSTICS_CAPTURED",
  "BLOCKED_NATIVE_CONTRACTION_DIAGNOSTICS_INCOMPLETE"];
$ACReport = <|
  "format"->"a-minus4-rgbeta-symbolic-contractions-v1",
  "status"->$ACStatus,"model"->"A_-4","stage"->"EFT1",
  "dimensions_alpha"->{1,3,2,-4},
  "source_model_path"->"RGE/running/rgbeta/T3RGBetaModel.wl",
  "source_model_sha256"->$ACSHA,
  "rgbeta_package_file"->ACString[FindFile["RGBeta`"]],
  "generator_contraction_probes"->$ACGroupProbes,
  "quartic_invariant_samples"->$ACQuarticRows,
  "native_dummy_index_double_call_probes"->$ACDoubleCalls,
  "native_projector_related_symbol_names"->$ACProjectorNames,
  "native_projector_symbol_definitions"->$ACSymbolRows,
  "native_group_reductions_certified"->False,
  "quartic_projection_certified"->False,
  "physical_matchete_coupling_map_certified"->False,
  "native_to_reference_axes_certified"->False,
  "source_modified"->False|>;
$ACJSON = CheckAbort[Check[ExportString[$ACReport,"RawJSON"],$Failed],$Aborted];
If[!StringQ[$ACJSON] || StringLength[$ACJSON]<100,
  Print["BLOCKED: native diagnostic is not JSON-serializable"]; Exit[5]];
If[!DirectoryQ[DirectoryName[$ACOut]],
  CreateDirectory[DirectoryName[$ACOut],CreateIntermediateDirectories->True]];
$ACTmp=$ACOut<>".partial";
If[CheckAbort[Check[Export[$ACTmp,$ACJSON,"Text"],$Failed],$Aborted]===$Failed ||
  !FileExistsQ[$ACTmp] || FileByteCount[$ACTmp]<100,
  Print["BLOCKED: native report temporary write failed"];Exit[5]];
$ACRead = CheckAbort[Check[Import[$ACTmp,"RawJSON"],$Failed],$Aborted];
If[!AssociationQ[$ACRead] || Lookup[$ACRead,"format",""]=!=
  "a-minus4-rgbeta-symbolic-contractions-v1",
  Print["BLOCKED: native report read-back failed"];Exit[5]];
If[CheckAbort[Check[RenameFile[$ACTmp,$ACOut],$Failed],$Aborted]===$Failed,
  Print["BLOCKED: native report finalize failed"];Exit[5]];
Print["Native RGBeta group contractions:",
  Count[$ACGroupProbes,x_/;TrueQ[Lookup[x,"native_exact_target_match",False]]],
  "/4 evaluated to tested targets"];
Print["Native quartet invariant diagnostics: ",
  Count[$ACQuarticRows,x_/;TrueQ[Lookup[x,"captured",False]]],"/5 captured"];
Print["Raw adjoint dummy-index reuse: ",
  Count[$ACDoubleCalls,x_/;TrueQ[Lookup[x,
    "same_literal_A_present_in_both_raw_calls",False]]],"/4 paired probes"];
Print["Quartic projector / group-reduction native symbols found: ",
  Length[$ACProjectorNames]];
Print["Physical basis / beta certified: NO"];
Print["Report: ",$ACOut];
Exit[If[$ACProbeOK,0,6]];
