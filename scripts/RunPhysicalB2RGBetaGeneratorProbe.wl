(* B2 G2 stage 2: fresh-kernel native *symbolic* generator/field-slot probe.
   Do NOT call directly: Python launcher supplies six-source-validated Git
   provenance through four environment variables. No RGE, matching or EFT.
   Captures no numeric eps/delS2 axes or physical action coefficients.
*)
ClearAll["Global`*"];
If[Length[Rest[$ScriptCommandLine]] =!= 1,
  Print["Usage: OUTPUT.json via Python launcher"]; Exit[2]];
$B6Output = ExpandFileName[Last[$ScriptCommandLine]];
If[FileExistsQ[$B6Output] || FileExistsQ[$B6Output <> ".partial"],
  Print["BLOCKED: refusing existing native JSON/partial: ", $B6Output]; Exit[2]];
$B6Flags = {"T3_EXPERIMENTAL_GENERAL_FERMIONS",
 "T3_EXPERIMENTAL_GENERAL_QUADRATICS",
 "T3_EXPERIMENTAL_GENERAL_QUARTICS",
 "T3_EXPERIMENTAL_FLAVOR_MASSES"};
If[AnyTrue[$B6Flags, Environment[#] === "1" &],
 Print["BLOCKED: unrestricted optional T3 flags"]; Exit[3]];
$B6Head = Environment["T3_B2_PROVENANCE_HEAD"];
$B6Blobs = <|
 "RGE/running/rgbeta/T3RGBetaModel.wl" -> Environment["T3_B2_PROVENANCE_RGBETA"],
 "Lagrangian/interactions/T3Topology.wl" -> Environment["T3_B2_PROVENANCE_TOPOLOGY"],
 "Lagrangian/model/SU2Invariants.wl" -> Environment["T3_B2_PROVENANCE_SU2"]|>;
B6SHAQ[x_] := StringQ[x] && StringMatchQ[x, RegularExpression["[0-9a-f]{40}"]];
If[!B6SHAQ[$B6Head] || !AllTrue[Values[$B6Blobs], B6SHAQ],
 Print["BLOCKED: Python source provenance missing"]; Exit[3]];
$B6Root = DirectoryName @ DirectoryName @ ExpandFileName[$InputFileName];
If[CheckAbort[Check[Needs["RGBeta`"]; True, $Failed], $Aborted] =!= True,
 Print["BLOCKED: installed RGBeta unavailable"]; Exit[4]];
Get[FileNameJoin[{$B6Root,"RGE","running","rgbeta","T3RGBetaModel.wl"}]];
ClearAll[B6String, B6Read, B6Gauge, B6Class, B6Index, B6DownValues,
 B6IndexScopeNames, B6NativeIndexMap];
B6String[x_] := ToString[x, InputForm];
B6Read[a_,name_String] := Module[{keys},
 If[!AssociationQ[a],Return[Missing["NotAssociation",name]]];
 keys=Select[Keys[a],Last[StringSplit[B6String[#],"`"]] === name &];
 If[Length[keys] =!= 1,Return[Missing["NotUnique",name]]];
 Lookup[a,First[keys]]];
B6Gauge[registry_,field_] := B6Read[Lookup[registry,field,<||>],"GaugeRep"];
B6Index[assoc_,field_] := Module[{v=If[AssociationQ[assoc],Lookup[assoc,field,-1],-1]},
 If[IntegerQ[v],v,-1]];
B6DownValues[names_List] := If[Length[names] === 1,
 With[{symbol=ToExpression[First[names]]}, B6String[DownValues[symbol]]],
 "unresolved native Tscal symbol"];
(* $fieldIndexMap is PackageScope, NOT PackageExport, in RGBeta's
   FieldsAndCouplings.m. An unqualified symbol in this Global` probe is
   not the live index registry. Resolve one actual RGBeta package symbol
   and read its own value, failing closed on absence or ambiguity. *)
B6IndexScopeNames[] := Names["RGBeta`*`$fieldIndexMap"];
B6NativeIndexMap[] := Module[{names = B6IndexScopeNames[], value},
 If[Length[names] =!= 1, Return[<||>]];
 value = Quiet@Check[ToExpression[First[names]], $Failed];
 If[AssociationQ[value], value, <||>]
];
$B6Classes = <|"A"->{1,3,2,0},"B"->{2,2,1,-1},
 "C"->{2,2,3,-1},"D"->{3,1,2,0},"E"->{3,3,2,0}|>;
B6Class[family_String] := Module[
 {d1,d2,df,alpha,yf,ok,second,repS1,repS2,repF,repY1,repY2,
  expectedS1,expectedS2,expectedF,indexMap,mapping,ids,regOK},
 {d1,d2,df,alpha} = $B6Classes[family];
 yf = (alpha+1)/2;
 ok = CheckAbort[Check[
   ResetModel[];T3RGBetaAddSM[];T3RGBetaAddHeavyFermion[df,yf];
   AddScalar[S1,GaugeRep->T3RGBetaGaugeRep[d1,alpha/2]];
   AddScalar[S2,GaugeRep->T3RGBetaGaugeRep[d2,(alpha+2)/2]];
   T3RGBetaAddT3Yukawas[d1,d2,df,yf];
   T3RGBetaAddMasses[d1,d2,df,yf];
   T3RGBetaAddBasicQuartics[d1,d2];
   T3RGBetaAddMixingQuartic[d1,d2]; True, $Failed],$Aborted];
 If[ok =!= True,Return[<|"family"->family,"dimensions_alpha"->{d1,d2,df,alpha},
    "registered"->False,"physical_action_certified"->False,
    "numeric_SU2_axes_certified"->False|>]];
 second = If[yf===0,F,Fc];
 repS1 = B6Gauge[$scalars,S1];repS2 = B6Gauge[$scalars,S2];
 repF = B6Gauge[$fermions,F];
 expectedS1 = T3RGBetaGaugeRep[d1,alpha/2];
 expectedS2 = T3RGBetaGaugeRep[d2,(alpha+2)/2];
 expectedF = T3RGBetaGaugeRep[df,yf];
 repY1 = B6Read[Lookup[$yukawas,y1,<||>],"Fields"];
 repY2 = B6Read[Lookup[$yukawas,y2,<||>],"Fields"];
 indexMap = B6NativeIndexMap[];
 mapping = Lookup[indexMap,"All",<||>];
 ids = <|"S1"->B6Index[mapping,S1],
  "Bar[S1]"->B6Index[mapping,Bar[S1]],
  "S2"->B6Index[mapping,S2],
  "Bar[S2]"->B6Index[mapping,Bar[S2]]|>;
 regOK = AssociationQ[$quartics] && KeyExistsQ[$quartics,lambdaT3] &&
   AssociationQ[$yukawas] && KeyExistsQ[$yukawas,y1] && KeyExistsQ[$yukawas,y2];
 <|"family"->family,"dimensions_alpha"->{d1,d2,df,alpha},
  "registered"->True,
  "scalar_representations_match_registration"->TrueQ[repS1===expectedS1 && repS2===expectedS2],
  "fermion_representations_match_registration"->TrueQ[repF===expectedF],
  "all_scalar_field_and_conjugate_indices_present"->
     TrueQ[AllTrue[Values[ids],IntegerQ[#] && #>0 &] &&
       Length[DeleteDuplicates[Values[ids]]]===4],
  "yukawa_field_orders_match"->TrueQ[regOK && repY1==={Bar[S1],l,F} &&
       repY2==={S2,l,second}],
  "S1_gauge_reps"->B6String[repS1],
  "S2_gauge_reps"->B6String[repS2],
  "F_gauge_reps"->B6String[repF],
  "field_index_ids"->ids,
  "native_package_field_index_map_available"->TrueQ[
    AssociationQ[indexMap] && KeyExistsQ[indexMap,"All"] &&
    AssociationQ[mapping]],
  "physical_action_certified"->False,"numeric_SU2_axes_certified"->False|>
];
$B6IndexSymbolNames = B6IndexScopeNames[];
$B6IndexSymbolOK = Length[$B6IndexSymbolNames]===1 &&
 StringStartsQ[First[$B6IndexSymbolNames],"RGBeta`"] &&
 StringEndsQ[First[$B6IndexSymbolNames],"`$fieldIndexMap"];
$B6TscalNames = Names["*`Tscal"];
$B6TscalDVs = B6DownValues[$B6TscalNames];
$B6TscalOK = Length[$B6TscalNames]===1 &&
 StringContainsQ[$B6TscalDVs,"tGen["] &&
 StringContainsQ[$B6TscalDVs,"Bar["] &&
 StringContainsQ[$B6TscalDVs,"$fieldIndexMap"];
$B6BarRule = TrueQ[
 tGen[Bar[SU2L[fund]],b6A,b6i,b6j] ===
 -tGen[SU2L[fund],b6A,b6j,b6i]];
$B6EPSAbstract = Head[eps[SU2L[fund],1,2]] === eps;
$B6S2Abstract = Head[delS2[SU2L,1,1,2]] === delS2;
$B6Rows = B6Class /@ Keys[$B6Classes];
$B6OK = $B6IndexSymbolOK && $B6TscalOK && $B6BarRule &&
 $B6EPSAbstract && $B6S2Abstract &&
 AllTrue[$B6Rows,TrueQ[Lookup[#,"registered",False]] &&
  TrueQ[Lookup[#,"scalar_representations_match_registration",False]] &&
  TrueQ[Lookup[#,"fermion_representations_match_registration",False]] &&
  TrueQ[Lookup[#,"all_scalar_field_and_conjugate_indices_present",False]] &&
  TrueQ[Lookup[#,"native_package_field_index_map_available",False]] &&
  TrueQ[Lookup[#,"yukawa_field_orders_match",False]] &];
$B6Data = <|
 "format"->"b2-rgbeta-native-symbolic-generator-orientation-v1",
 "status"->If[$B6OK,"NATIVE_SYMBOLIC_GENERATOR_RULES_CAPTURED_NUMERIC_AXES_PENDING",
                "BLOCKED_NATIVE_SYMBOLIC_GENERATOR_WITNESS"],
 "git_head_at_capture"->$B6Head,"git_source_blobs_at_capture"->$B6Blobs,
 "git_provenance_capture_method"->"Python six-source clean-filter guard",
 "rgbeta_package_file"->B6String[FindFile["RGBeta`"]],
 "field_index_map_symbol_names"->$B6IndexSymbolNames,
 "field_index_map_unique_package_scope"->$B6IndexSymbolOK,
 "Tscal_symbol_names"->$B6TscalNames,
 "Tscal_downvalues_inputform"->$B6TscalDVs,
 "Tscal_rule_has_conjugate_scalar_blocks"->$B6TscalOK,
 "tGen_bar_transpose_identity_native"->$B6BarRule,
 "eps_sample_remains_symbolic"->$B6EPSAbstract,
 "delS2_sample_remains_symbolic"->$B6S2Abstract,
 "classes"->$B6Rows,
 "RGBeta_numeric_eps_components_certified"->False,
 "RGBeta_numeric_delS2_components_certified"->False,
 "RGBeta_to_Matchete_basis_identified"->False,
 "physical_coupling_dictionary_certified"->False,
 "physical_source_certified_models"->0,
 "safe_to_integrate_into_production"->False,"EFT_performed"->False|>;
$B6Encoded=CheckAbort[Check[ExportString[$B6Data,"RawJSON"],$Failed],$Aborted];
If[!StringQ[$B6Encoded] || StringLength[$B6Encoded]<32,
 Print["BLOCKED: native JSON encoding failed"];Exit[5]];
$B6Verify=CheckAbort[Check[ImportString[$B6Encoded,"RawJSON"],$Failed],$Aborted];
If[!AssociationQ[$B6Verify] ||
 Lookup[$B6Verify,"format",None]=!="b2-rgbeta-native-symbolic-generator-orientation-v1",
 Print["BLOCKED: native JSON roundtrip failed"];Exit[5]];
If[!DirectoryQ[DirectoryName[$B6Output]],
 CreateDirectory[DirectoryName[$B6Output],CreateIntermediateDirectories->True]];
$B6Tmp=$B6Output<>".partial";
$B6Wrote=CheckAbort[Check[
 Export[$B6Tmp,$B6Encoded,"Text",CharacterEncoding->"UTF-8"],$Failed],$Aborted];
If[$B6Wrote===$Failed || !FileExistsQ[$B6Tmp] || FileByteCount[$B6Tmp]<32,
 Print["BLOCKED: native diagnostic output write failed"];Exit[6]];
$B6Readback=CheckAbort[Check[Import[$B6Tmp,"RawJSON"],$Failed],$Aborted];
If[!AssociationQ[$B6Readback] ||
 Lookup[$B6Readback,"format",None]=!="b2-rgbeta-native-symbolic-generator-orientation-v1",
 Print["BLOCKED: native diagnostic JSON readback failed"];Exit[6]];
$B6Final=CheckAbort[Check[RenameFile[$B6Tmp,$B6Output],$Failed],$Aborted];
If[$B6Final===$Failed || !FileExistsQ[$B6Output],
 Print["BLOCKED: native diagnostic JSON rename failed"];Exit[6]];
If[!$B6OK,
 Print["BLOCKED: native symbolic generator checks not all measured; inspect ",$B6Output];Exit[7]];
Print["B2 G2 native symbolic generator 5/5; tGen Bar=yes; numeric axes=NO; physical=0/16: ",$B6Output];
Exit[0];
