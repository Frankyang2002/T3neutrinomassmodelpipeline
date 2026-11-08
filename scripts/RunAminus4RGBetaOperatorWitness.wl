(* A(-4) EFT1 read-only native RGBeta quartic-operator witness.

   wolframscript -file scripts/RunAminus4RGBetaOperatorWitness.wl OUTPUT.json

   Loads the installed RGBeta package and the project's unchanged EFT1
   builder, then records the ACTUAL $quartics registry and the package's
   formal Lam[a,b,c,d] tensor. No tGen/delS2 component matrices are injected.
   This is a symbolic native witness, not an identification of RGBeta's
   component basis with Matchete or a physical quartic beta certificate.
   This script does not call BetaTerm, run RGEs, edit source, or overwrite
   existing native evidence.
*)
ClearAll["Global`*"];
If[Length[Rest[$ScriptCommandLine]] =!= 1,
  Print["Usage: OUTPUT.json"]; Exit[2]];
$OPOut = ExpandFileName[Last[$ScriptCommandLine]];
If[FileExistsQ[$OPOut] || FileExistsQ[$OPOut <> ".partial"],
  Print["Refusing to overwrite native witness/partial: ", $OPOut]; Exit[2]];
$OPForbiddenFlags = {"T3_EXPERIMENTAL_GENERAL_FERMIONS",
  "T3_EXPERIMENTAL_GENERAL_QUADRATICS",
  "T3_EXPERIMENTAL_GENERAL_QUARTICS",
  "T3_EXPERIMENTAL_FLAVOR_MASSES"};
If[AnyTrue[$OPForbiddenFlags, Environment[#] === "1" &],
  Print["BLOCKED: unrestricted T3 feature flag enabled"]; Exit[3]];

$OPRoot = DirectoryName @ DirectoryName @ ExpandFileName[$InputFileName];
$OPModel = FileNameJoin[{$OPRoot, "RGE", "running", "rgbeta", "T3RGBetaModel.wl"}];
If[!FileExistsQ[$OPModel], Print["Missing T3RGBetaModel.wl"]; Exit[3]];
$OPSourceSHA = Quiet @ Check[
  IntegerString[FileHash[$OPModel, "SHA256"], 16, 64], $Failed];
If[!StringQ[$OPSourceSHA] || StringLength[$OPSourceSHA] =!= 64,
  Print["Cannot hash local native model source"]; Exit[3]];
If[CheckAbort[Check[Needs["RGBeta`"]; True, $Failed], $Aborted] =!= True,
  Print["Cannot load installed RGBeta package"]; Exit[4]];
Get[$OPModel];
$OPBuild = CheckAbort[Quiet @ Check[
  T3RGBetaBuildEFT1[1, 3, 2, -4], $Failed], $Aborted];
If[!AssociationQ[$OPBuild], Print["A(-4) EFT1 build failed"]; Exit[4]];

ClearAll[OPString, OPRead, OPIndexMap, OPIndex, OPOperator, OPUnitRules,
  OPFormalSample, OPNumericQ];
OPString[x_] := ToString[x, InputForm];
OPRead[rec_, name_String] := Module[{matches},
  If[!AssociationQ[rec], Return[Missing["NotAssociation"]]];
  matches = Select[Keys[rec],
    Last[StringSplit[OPString[#], "`"]] === name &];
  If[Length[matches] =!= 1, Return[Missing["MissingOrAmbiguous", name]]];
  Lookup[rec, First[matches]]
];
OPIndexMap[] := Module[{names, data},
  names = Names["RGBeta`PackageScope`$fieldIndexMap"];
  If[Length[names] =!= 1, Return[<||>]];
  data = Quiet @ Check[ToExpression[First[names]], $Failed];
  If[AssociationQ[data], Lookup[data, "Scalars", <||>], <||>]
];
OPIndex[indices_, phi_] := Module[{index},
  index = If[AssociationQ[indices], Lookup[indices, phi, -1], -1];
  If[IntegerQ[index] && index > 0, index, -1]
];
OPOperator[name_String] := Module[
  {symbol, info, fields, inv, cup, cupBar, arrangements, self,
   sample, sampleCup, sampleCupBar, good},
  symbol = ToExpression[name];
  info = If[AssociationQ[$quartics], Lookup[$quartics, symbol, <||>], <||>];
  fields = OPRead[info, "Fields"];
  inv = OPRead[info, "Invariant"];
  cup = OPRead[info, "Coupling"];
  cupBar = OPRead[info, "CouplingBar"];
  arrangements = OPRead[info, "UniqueArrangements"];
  self = OPRead[info, "SelfConjugate"];
  good = ListQ[fields] && Length[fields] === 4 &&
    Head[inv] =!= Missing && Head[cup] =!= Missing &&
    Head[cupBar] =!= Missing && ListQ[arrangements] &&
    MemberQ[{True, False}, self];
  sample = If[good, Quiet @ Check[inv[opA, opB, opC, opD], $Failed], $Failed];
  sampleCup = If[good, Quiet @ Check[cup[], $Failed], $Failed];
  sampleCupBar = If[good, Quiet @ Check[cupBar[], $Failed], $Failed];
  <|"name" -> name,
    "registration_complete" -> TrueQ[good && sample =!= $Failed &&
      sampleCup =!= $Failed && sampleCupBar =!= $Failed],
    "fields" -> OPString[fields],
    "invariant_inputform" -> OPString[sample],
    "coupling_inputform" -> OPString[sampleCup],
    "coupling_bar_inputform" -> OPString[sampleCupBar],
    "self_conjugate" -> If[MemberQ[{True, False}, self], self, Null],
    "native_unique_arrangements" -> If[ListQ[arrangements], Length[arrangements], -1]
  |>
];
OPNumericQ[v_] := NumericQ[v] && FreeQ[v, _Real | Indeterminate |
  ComplexInfinity | DirectedInfinity | _Missing];

$OPNames = {"lambdaH", "lambdaS1", "lambdaS2", "lambdaH1",
  "lambdaH2", "lambda12", "lambdaT3", "lambdaH2Adj", "lambdaS2Adj"};
$OPNativeNames = If[AssociationQ[$quartics], Sort[SymbolName /@ Keys[$quartics]], {}];
$OPRows = OPOperator /@ $OPNames;
$OPScalars = OPIndexMap[];
$OPScalarIndices = <|
  "H" -> OPIndex[$OPScalars, H],
  "Bar[H]" -> OPIndex[$OPScalars, Bar[H]],
  "S1" -> OPIndex[$OPScalars, S1],
  "Bar[S1]" -> OPIndex[$OPScalars, Bar[S1]],
  "S2" -> OPIndex[$OPScalars, S2],
  "Bar[S2]" -> OPIndex[$OPScalars, Bar[S2]]|>;

(* Probe the installed package's actual tensor, not a user-defined replacement.
   The resulting TStructure/SparseArray may retain symbolic group indices.
   Both the native source code and native output are inspected downstream. *)
$OPLamNames = Names["RGBeta`PackageScope`Lam"];
$OPLam = If[Length[$OPLamNames] === 1,
  CheckAbort[TimeConstrained[Quiet @ Check[
    With[{fn = ToExpression[First[$OPLamNames]]}, fn[opA, opB, opC, opD]],
    $Failed], 120, $TimedOut], $Aborted], $Failed];
$OPLamString = If[MemberQ[{$Failed, $Aborted, $TimedOut}, $OPLam],
  "", OPString[$OPLam]];
$OPLamSymbols = Select[$OPNames, StringContainsQ[$OPLamString, #] &];
$OPLamAvailable = StringLength[$OPLamString] > 50 &&
  StringContainsQ[$OPLamString, "TStructure["] &&
  StringContainsQ[$OPLamString, "SparseArray["] &&
  Sort[$OPLamSymbols] === Sort[$OPNames] &&
  !StringStartsQ[$OPLamString, "Lam["];

(* Do not fake RGBeta native numeric component axes. Ask the installed
   symbolic group primitives with integer indices; they may remain abstract. *)
$OPFundNumerics = Table[tGen[SU2L[fund], a, i, j],
  {a, 1, 3}, {i, 1, 2}, {j, 1, 2}];
$OPTripletNumerics = Table[tGen[SU2L[S2], a, i, j],
  {a, 1, 3}, {i, 1, 3}, {j, 1, 3}];
$OPTripletCGNumerics = Table[delS2[SU2L, a, i, j],
  {a, 1, 3}, {i, 1, 2}, {j, 1, 2}];
$OPFundNumeric = AllTrue[Flatten[$OPFundNumerics], OPNumericQ];
$OPTripletNumeric = AllTrue[Flatten[$OPTripletNumerics], OPNumericQ];
$OPCGNumeric = AllTrue[Flatten[$OPTripletCGNumerics], OPNumericQ];

$OPRegistryOK = Sort[$OPNames] === $OPNativeNames &&
  AllTrue[$OPRows, TrueQ[#["registration_complete"]] &] &&
  AllTrue[Values[$OPScalarIndices], IntegerQ[#] && # > 0 &] &&
  Length[DeleteDuplicates[Values[$OPScalarIndices]]] === 6;
$OPStatus = If[$OPRegistryOK && $OPLamAvailable,
  "NATIVE_SYMBOLIC_OPERATOR_TENSOR_CAPTURED_COMPONENT_MAP_PENDING",
  "BLOCKED_INCOMPLETE_NATIVE_OPERATOR_WITNESS"];
$OPReport = <|
  "format" -> "a-minus4-rgbeta-operator-witness-v1",
  "status" -> $OPStatus,
  "model" -> "A_-4", "stage" -> "EFT1", "dimensions_alpha" -> {1, 3, 2, -4},
  "source_model_path" -> "RGE/running/rgbeta/T3RGBetaModel.wl",
  "source_model_sha256" -> $OPSourceSHA,
  "rgbeta_package_file" -> OPString[FindFile["RGBeta`"]],
  "native_quartic_symbols" -> $OPNativeNames,
  "scalar_species_indices" -> $OPScalarIndices,
  "quartic_operators" -> $OPRows,
  "native_Lam_symbol_names" -> $OPLamNames,
  "native_Lam_formal_tensor_available" -> TrueQ[$OPLamAvailable],
  "native_Lam_visible_quartic_symbols" -> $OPLamSymbols,
  "native_Lam_formal_tensor_inputform" -> If[$OPLamAvailable,
    StringTake[$OPLamString, UpTo[2000000]], ""],
  "native_Lam_string_length" -> StringLength[$OPLamString],
  "native_Lam_truncated" -> TrueQ[StringLength[$OPLamString] > 2000000],
  "native_fund_generator_components_numeric" -> TrueQ[$OPFundNumeric],
  "native_triplet_generator_components_numeric" -> TrueQ[$OPTripletNumeric],
  "native_triplet_CG_components_numeric" -> TrueQ[$OPCGNumeric],
  "native_fund_generator_components_inputform" ->
    If[$OPFundNumeric, OPString[$OPFundNumerics], "unavailable"],
  "native_triplet_generator_components_inputform" ->
    If[$OPTripletNumeric, OPString[$OPTripletNumerics], "unavailable"],
  "native_triplet_CG_components_inputform" ->
    If[$OPCGNumeric, OPString[$OPTripletCGNumerics], "unavailable"],
  "native_operator_registry_complete" -> TrueQ[$OPRegistryOK],
  "matchete_physical_coupling_map_certified" -> False,
  "native_to_reference_component_axes_certified" -> False,
  "rgbeta_quartic_betas_certified" -> False,
  "production_modified" -> False
|>;
$OPEncoded = CheckAbort[Check[ExportString[$OPReport, "RawJSON"], $Failed], $Aborted];
If[!StringQ[$OPEncoded] || StringLength[$OPEncoded] < 100,
  Print["Native RawJSON serialization failed"]; Exit[5]];
If[!DirectoryQ[DirectoryName[$OPOut]],
  CreateDirectory[DirectoryName[$OPOut], CreateIntermediateDirectories -> True]];
$OPPartial = $OPOut <> ".partial";
$OPWrite = CheckAbort[Check[Export[$OPPartial, $OPEncoded, "Text"], $Failed], $Aborted];
If[$OPWrite === $Failed || !FileExistsQ[$OPPartial] ||
  FileByteCount[$OPPartial] < 100,
  Print["Native JSON temporary write failed"]; Exit[5]];
$OPReadback = CheckAbort[Check[Import[$OPPartial, "RawJSON"], $Failed], $Aborted];
If[!AssociationQ[$OPReadback] ||
  Lookup[$OPReadback, "format", ""] =!= "a-minus4-rgbeta-operator-witness-v1",
  Print["Native JSON readback failed"]; Exit[5]];
If[CheckAbort[Check[RenameFile[$OPPartial, $OPOut], $Failed], $Aborted] === $Failed,
  Print["Native JSON finalization failed"]; Exit[5]];
Print["A(-4) EFT1 native operator registry: ", If[$OPRegistryOK, "PASS", "BLOCKED"]];
Print["Installed RGBeta Lam tensor: ", If[$OPLamAvailable, "CAPTURED", "BLOCKED"]];
Print["Native numeric SU(2) triplet axes: ", If[$OPTripletNumeric && $OPCGNumeric,
  "AVAILABLE", "NOT EXPOSED"]];
Print["Matchete/physical tensor map certified: NO"];
Print["Report: ", $OPOut];
Exit[If[$OPStatus === "NATIVE_SYMBOLIC_OPERATOR_TENSOR_CAPTURED_COMPONENT_MAP_PENDING", 0, 6]];
