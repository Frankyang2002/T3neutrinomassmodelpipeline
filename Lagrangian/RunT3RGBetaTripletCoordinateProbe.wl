(* Workstream A: probe whether the INSTALLED RGBeta runtime supplies numeric
   SU(2) symmetric-triplet coordinates. No user-selected component conventions
   are substituted. Does not modify registered source files or call Match[].

   wolframscript -file Lagrangian/RunT3RGBetaTripletCoordinateProbe.wl OUTPUT.json

   Literal integer index probes are a diagnostic only. Failure to resolve
   these indices is NOT a finding that the abstract representation is wrong.
*)
ClearAll["Global`*"];
$T3CPArgs = Rest[$ScriptCommandLine];
If[Length[$T3CPArgs] =!= 1, Print["Usage: OUTPUT.json"]; Exit[2]];
$T3CPOut = ExpandFileName[$T3CPArgs[[1]]];
If[FileExistsQ[$T3CPOut], Print["Refusing to overwrite ", $T3CPOut]; Exit[2]];
$T3CPDir = DirectoryName @ ExpandFileName[$InputFileName];
$T3CPRoot = DirectoryName[$T3CPDir];
$T3CPLoad = CheckAbort[Check[Needs["RGBeta`"]; True, $Failed], $Aborted];
If[$T3CPLoad =!= True, Print["RGBeta load failed"]; Exit[3]];
Get[FileNameJoin[{$T3CPRoot,"RGE","running","rgbeta","T3RGBetaModel.wl"}]];
$T3CPBuild = CheckAbort[Check[T3RGBetaBuild[1,3,2,-2], $Failed], $Aborted];
If[!AssociationQ[$T3CPBuild], Print["RGBeta model build failed"]; Exit[4]];

ClearAll[T3CPString, T3CPExactNumericQ, T3CPRecord];
T3CPString[x_] := ToString[x, InputForm];
T3CPExactNumericQ[x_] := TrueQ[NumericQ[x]] &&
  FreeQ[x, _Real | Indeterminate | ComplexInfinity | DirectedInfinity | _Missing];
T3CPRecord[array_, shape_List] := Module[{flat, exactQ, strings, head, sample},
  flat = Flatten[array];
  exactQ = Dimensions[array] === shape && Length[flat] === Times@@shape &&
    AllTrue[flat, T3CPExactNumericQ];
  strings = If[TrueQ[exactQ], T3CPString /@ flat, {}];
  sample = T3CPString /@ Take[flat, UpTo[5]];
  <|"shape" -> Dimensions[array], "expected_shape" -> shape,
    "exact_numeric_components_exposed" -> TrueQ[exactQ],
    "components_inputform" -> strings,
    "diagnostic_first_values_inputform" -> sample|>
];

(* Evaluate the actually installed native symbols with integer index labels.
   Avoid giving these symbols definitions or injecting guessed matrices. *)
$T3CPDel = Table[delS2[SU2L, k, i, j], {k,1,3}, {i,1,2}, {j,1,2}];
$T3CPGen = Table[tGen[SU2L[S2], a, i, j], {a,1,3}, {i,1,3}, {j,1,3}];
$T3CPEps = Table[eps[SU2L[fund], i, j], {i,1,2}, {j,1,2}];
$T3CPD = T3CPRecord[$T3CPDel,{3,2,2}];
$T3CPG = T3CPRecord[$T3CPGen,{3,3,3}];
$T3CPE = T3CPRecord[$T3CPEps,{2,2}];
$T3CPBoth = TrueQ[$T3CPD["exact_numeric_components_exposed"]] &&
             TrueQ[$T3CPG["exact_numeric_components_exposed"]];
$T3CPResult = <|
  "format" -> "t3-workstream-a-rgbeta-triplet-coordinate-probe-v1",
  "status" -> If[$T3CPBoth,"NativeRGBetaNumericTripletCandidatesExposed",
      "NativeRGBetaTripletCoordinatesNotExposedByIntegerIndexProbe"],
  "source" -> "installed RGBeta runtime after T3RGBetaBuild[1,3,2,-2]",
  "model" -> "A_-2",
  "native_triplet_CG_delS2" -> $T3CPD,
  "native_triplet_generators_tGen" -> $T3CPG,
  "native_fundamental_epsilon" -> $T3CPE,
  "integer_index_probe_only" -> True,
  "numeric_triplet_CG_and_generators_both_exposed" -> $T3CPBoth,
  "cross_package_component_basis_certified" -> False,
  "barred_field_map_certified" -> False,
  "signed_physical_coupling_dictionary_certified" -> False,
  "production_modified" -> False
|>;
If[!DirectoryQ[DirectoryName[$T3CPOut]],
  CreateDirectory[DirectoryName[$T3CPOut],CreateIntermediateDirectories->True]];
$T3CPWrite = Quiet @ Check[Export[$T3CPOut,$T3CPResult,"RawJSON"],$Failed];
If[$T3CPWrite === $Failed || !FileExistsQ[$T3CPOut],
  Print["Failed to write RGBeta probe JSON"]; Exit[5]];
Print["RGBeta native triplet coordinate probe: ",$T3CPResult["status"]];
Print["delS2 numeric: ",$T3CPD["exact_numeric_components_exposed"],
      "; tGen numeric: ",$T3CPG["exact_numeric_components_exposed"]];
Print["Report: ",$T3CPOut];
Exit[0];
