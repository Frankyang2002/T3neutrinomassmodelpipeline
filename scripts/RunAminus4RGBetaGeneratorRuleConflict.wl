(* Read-only A(-4) EFT1 RGBeta generic-generator versus SU(2) Fierz witness.
   Run in fresh kernel:
   wolframscript -file scripts/RunAminus4RGBetaGeneratorRuleConflict.wl OUTPUT.json
   This script does not modify installed RGBeta or T3 production source.
   It tests actual Mathematica rule precedence; source-level disagreement alone
   is NOT promoted to a beta-function diagnosis or Matchete certificate. *)
ClearAll["Global`*"];
If[Length[Rest[$ScriptCommandLine]] =!= 1,
  Print["Usage: OUTPUT.json"]; Exit[2]];
$RGOut = ExpandFileName[Last[$ScriptCommandLine]];
If[FileExistsQ[$RGOut] || FileExistsQ[$RGOut <> ".partial"],
  Print["BLOCKED: output or partial already exists: ", $RGOut]; Exit[2]];
$RGFlags = {"T3_EXPERIMENTAL_GENERAL_FERMIONS",
  "T3_EXPERIMENTAL_GENERAL_QUADRATICS",
  "T3_EXPERIMENTAL_GENERAL_QUARTICS",
  "T3_EXPERIMENTAL_FLAVOR_MASSES"};
If[AnyTrue[$RGFlags, Environment[#] === "1" &],
  Print["BLOCKED: experimental T3 flag enabled"]; Exit[3]];
$RGRoot = DirectoryName @ DirectoryName @ ExpandFileName[$InputFileName];
$RGModel = FileNameJoin[{$RGRoot, "RGE", "running", "rgbeta", "T3RGBetaModel.wl"}];
If[!FileExistsQ[$RGModel], Print["BLOCKED: model source missing"]; Exit[3]];
$RGModelSHA = Quiet@Check[IntegerString[FileHash[$RGModel, "SHA256"], 16, 64], $Failed];
If[!StringQ[$RGModelSHA] || StringLength[$RGModelSHA] =!= 64,
  Print["BLOCKED: model hash unavailable"]; Exit[3]];
If[CheckAbort[Check[Needs["RGBeta`"]; True, $Failed], $Aborted] =!= True,
  Print["BLOCKED: installed RGBeta unavailable"]; Exit[4]];
$RGInit = FindFile["RGBeta`"];
$RGGroupSource = FileNameJoin[{DirectoryName@DirectoryName@$RGInit,
  "Package", "GroupsAndIndices.m"}];
$RGGroupHash = If[FileExistsQ[$RGGroupSource],
  Quiet@Check[IntegerString[FileHash[$RGGroupSource,"SHA256"],16,64],$Failed],
  $Failed];
If[!StringQ[$RGGroupHash] || StringLength[$RGGroupHash] =!= 64,
  Print["BLOCKED: cannot hash installed GroupsAndIndices.m at: ",
    $RGGroupSource]; Exit[4]];
Get[$RGModel];
$RGBuild = CheckAbort[Quiet@Check[
    T3RGBetaBuildEFT1[1,3,2,-4],$Failed],$Aborted];
If[!AssociationQ[$RGBuild],
  Print["BLOCKED: native A(-4) EFT1 build failed"];Exit[4]];

ClearAll[RGString, RGEval, RGCase, RGComparison];
RGString[x_] := ToString[x, InputForm];
RGEval[held_HoldComplete] := CheckAbort[
  TimeConstrained[Quiet@Check[
    RefineGroupStructures[ReleaseHold[held]],$Failed],20,$TimedOut],$Aborted];
RGComparison[x_, y_] := If[
  MemberQ[{$Failed,$TimedOut,$Aborted},x] ||
  MemberQ[{$Failed,$TimedOut,$Aborted},y],$Failed,
  CheckAbort[TimeConstrained[Quiet@Check[
    FullSimplify[Expand[x-y]],$Failed],20,$TimedOut],$Aborted]];
RGCase[id_, directHeld_HoldComplete, genericExpectedHeld_HoldComplete] := Module[
  {direct, expected, diff},
  direct=RGEval[directHeld];
  expected=ReleaseHold[genericExpectedHeld];
  diff=RGComparison[direct,expected];
  <|"id"->id,
    "input_inputform"->RGString[directHeld],
    "native_reduced_inputform"->RGString[direct],
    "target_inputform"->RGString[expected],
    "difference_inputform"->RGString[diff],
    "native_matches_target"->TrueQ[diff===0],
    "native_evaluated"->!MemberQ[{$Failed,$TimedOut,$Aborted},direct]|>
];

(* Same-orientation internal fund index, summed over gauge A and fund j.
   A generic rule in GroupsAndIndices.m:190 predicts -3/4 delta(i,k).
   The independent SU(2) Fierz rule:430-431 predicts +1/4 delta(i,k).
   Crucially, first expand four-DISTINCT-index pair before identifying l=j;
   this avoids the competing generic rule at the initial evaluation stage. *)
$RGDirect = RGEval[HoldComplete[
  tGen[SU2L[fund], rgA, rgi, rgj] *
  tGen[SU2L[fund], rgA, rgk, rgj]]];
$RGFourDistinct = RGEval[HoldComplete[
  tGen[SU2L[fund], rgA, rgi, rgj] *
  tGen[SU2L[fund], rgA, rgk, rgl]]];
$RGViaFierz = If[MemberQ[{$Failed,$TimedOut,$Aborted},$RGFourDistinct],
  $Failed,CheckAbort[TimeConstrained[Quiet@Check[
    RefineGroupStructures[$RGFourDistinct /. rgl -> rgj],$Failed],20,$TimedOut],
  $Aborted]];
$RGGenericTarget = -(3/4)*del[SU2L[fund],rgi,rgk];
$RGFierzTarget = (1/4)*del[SU2L[fund],rgi,rgk];
$RGDirectVsGeneric = RGComparison[$RGDirect,$RGGenericTarget];
$RGDirectVsFierz = RGComparison[$RGDirect,$RGFierzTarget];
$RGViaFierzVsFierz = RGComparison[$RGViaFierz,$RGFierzTarget];
$RGPathDiff = RGComparison[$RGDirect,$RGViaFierz];

(* The reverse-index Casimir control is distinct from the same-orientation
   contraction: the former is a valid SU(2) identity and passed earlier. *)
$RGReverseCase = RGCase["fundamental_reverse_index_casimir_control",
  HoldComplete[tGen[SU2L[fund],rgA,rgi,rgj]*
               tGen[SU2L[fund],rgA,rgj,rgk]],
  HoldComplete[(3/4)*del[SU2L[fund],rgi,rgk]]];
$RGTripletCase = RGCase["triplet_same_orientation_native_only",
  HoldComplete[tGen[SU2L[S2],rgA,rgi,rgj]*
               tGen[SU2L[S2],rgA,rgk,rgj]],
  HoldComplete[-2*del[SU2L[S2],rgi,rgk]]];
$RGConflict = TrueQ[$RGPathDiff =!= 0] &&
  !MemberQ[{$Failed,$TimedOut,$Aborted},$RGPathDiff] &&
  TrueQ[$RGViaFierzVsFierz===0];
$RGStatus = If[$RGConflict,
  "NATIVE_SAME_ORIENTATION_FIERZ_PATH_DISAGREEMENT",
  If[TrueQ[$RGPathDiff===0],
    "NATIVE_SAME_ORIENTATION_PATHS_AGREE",
    "NATIVE_SAME_ORIENTATION_PATH_COMPARISON_UNRESOLVED"]];
$RGReport = <|
  "format"->"a-minus4-rgbeta-generator-rule-conflict-v1",
  "status"->$RGStatus,
  "model"->"A_-4","stage"->"EFT1",
  "model_source_sha256"->$RGModelSHA,
  "installed_groups_source_path"->$RGGroupSource,
  "installed_groups_source_sha256"->$RGGroupHash,
  "installed_rgbeta_init"->RGString[$RGInit],
  "test_kind"->"same-orientation generator contraction; no beta functions",
  "fund_direct_inputform"->RGString[HoldComplete[
    tGen[SU2L[fund],rgA,rgi,rgj]*tGen[SU2L[fund],rgA,rgk,rgj]]],
  "fund_direct_reduction_inputform"->RGString[$RGDirect],
  "fund_four_distinct_reduction_inputform"->RGString[$RGFourDistinct],
  "fund_fierz_then_contract_reduction_inputform"->RGString[$RGViaFierz],
  "fund_expected_generic_inputform"->RGString[$RGGenericTarget],
  "fund_expected_fierz_inputform"->RGString[$RGFierzTarget],
  "fund_direct_minus_generic_inputform"->RGString[$RGDirectVsGeneric],
  "fund_direct_minus_fierz_inputform"->RGString[$RGDirectVsFierz],
  "fund_fierz_minus_expected_inputform"->RGString[$RGViaFierzVsFierz],
  "fund_direct_minus_fierz_path_inputform"->RGString[$RGPathDiff],
  "native_direct_matches_generic"->TrueQ[$RGDirectVsGeneric===0],
  "native_direct_matches_fierz"->TrueQ[$RGDirectVsFierz===0],
  "native_fierz_path_matches_reference"->TrueQ[$RGViaFierzVsFierz===0],
  "native_two_paths_disagree"->TrueQ[$RGConflict],
  "reverse_index_control"->$RGReverseCase,
  "triplet_same_orientation_diagnostic"->$RGTripletCase,
  "physical_scalar_beta_error_proven"->False,
  "physical_matchete_coupling_map_certified"->False,
  "source_modified"->False|>;
$RGJSON = CheckAbort[Check[ExportString[$RGReport,"RawJSON"],$Failed],$Aborted];
If[!StringQ[$RGJSON] || StringLength[$RGJSON] < 150,
  Print["BLOCKED: output is not JSON serializable"];Exit[5]];
If[!DirectoryQ[DirectoryName[$RGOut]],
  CreateDirectory[DirectoryName[$RGOut],CreateIntermediateDirectories->True]];
$RGPartial = $RGOut <> ".partial";
If[CheckAbort[Check[Export[$RGPartial,$RGJSON,"Text"],$Failed],$Aborted] === $Failed ||
  !FileExistsQ[$RGPartial] || FileByteCount[$RGPartial] < 150,
  Print["BLOCKED: native JSON write failed"];Exit[5]];
$RGRead = CheckAbort[Check[Import[$RGPartial,"RawJSON"],$Failed],$Aborted];
If[!AssociationQ[$RGRead] || Lookup[$RGRead,"format",""] =!=
  "a-minus4-rgbeta-generator-rule-conflict-v1",
  Print["BLOCKED: native JSON round trip failed"];Exit[5]];
If[CheckAbort[Check[RenameFile[$RGPartial,$RGOut],$Failed],$Aborted] === $Failed,
  Print["BLOCKED: report finalize failed"];Exit[5]];
Print["Native direct/fierz paths: ",$RGStatus];
Print["Native fundamental direct matches generic: ",
  TrueQ[$RGDirectVsGeneric===0]];
Print["Native fundamental direct matches Fierz: ",
  TrueQ[$RGDirectVsFierz===0]];
Print["Native fundamental Fierz path matches SU(2) reference: ",
  TrueQ[$RGViaFierzVsFierz===0]];
Print["Triplet same-orientation native matches generic target: ",
  TrueQ[$RGTripletCase["native_matches_target"]]];
Print["Scalar beta mismatch explained: NOT YET"];
Print["Report: ",$RGOut];
Exit[0];
