(*
  General D(-2) provisional C5: exact 1/epsilon residue audit.
  Read-only: does not run Match[], change the candidate, renormalize,
  or set a dimensional regulator numerically.

  For the NF=1 benchmark ONLY, finite heavy-flavour sums and NFlavor
  indices are reduced exactly, before any floating-point substitution.
  This tests algebraic cancellation of a simple dimensional pole in the
  already assembled candidate, not full physical renormalization.

  Usage:
    wolframscript -file Lagrangian/RunGeneralDMinus2C5PoleCancellationAudit.wl \
      FIRST_THRESHOLD_OUTPUT_DIR 1
*)
ClearAll["Global`*"];
args = Rest[$ScriptCommandLine];
If[Length[args] =!= 2,
  Print["Usage: FIRST_THRESHOLD_OUTPUT_DIR NF (currently 1)"]; Exit[2]
];
outputDir = ExpandFileName[args[[1]]];
nf = Quiet @ Check[ToExpression[args[[2]], InputForm], $Failed];
If[nf =!= 1,
  Print["[C5 pole] Only the independently validated NF=1 contraction is supported."];
  Exit[3]
];
profileDir = FileNameJoin[{outputDir, "data", "general_d_m2_scalar_threshold", "profile_tree"}];
root = FileNameJoin[{profileDir, "profile_c5_pole_cancellation"}];
reportPath = root <> "_result.json";
residuePath = root <> "_residue_exact.wxf";
residueTextPath = root <> "_residue_exact.txt";
candidatePath = FileNameJoin[{profileDir, "profile_c5_flavour_candidate_total_provisional.wxf"}];
benchmarkReportPath = FileNameJoin[{profileDir, "profile_c5_rough_benchmark_result.json"}];
flavourReportPath = FileNameJoin[{profileDir, "profile_c5_flavour_candidate_result.json"}];

ClearAll[D2Text, D2Name, D2JSONSafe, D2Fail, D2DimRegQ,
  D2PoleResidue, D2ExpandFiniteSums, D2ReduceNF1, D2TermCount];
D2Text[x_] := ToString[x, InputForm];
D2Name[s_Symbol] := SymbolName[Unevaluated[s]];
D2Name[x_] := D2Text[x];
D2JSONSafe[x_Association] := Association @ KeyValueMap[ToString[#1] -> D2JSONSafe[#2] &, x];
D2JSONSafe[x_List] := D2JSONSafe /@ x;
D2JSONSafe[x_String] := x;
D2JSONSafe[x_Integer] := x;
D2JSONSafe[True] := True;
D2JSONSafe[False] := False;
D2JSONSafe[x_Real] := Module[{v = N[x, MachinePrecision]}, If[TrueQ[v == 0], 0, v]];
D2JSONSafe[x_] := D2Text[x];
D2Fail[code_Integer, msg_String, details_: <||>] := (
  Print["[C5 pole] ERROR: ", msg];
  Export[reportPath, D2JSONSafe[Join[<|
    "Status" -> "Failed", "Reason" -> msg, "ExitCode" -> code,
    "NoNewMatchPerformed" -> True,
    "RenormalizationCertified" -> False|>, details]], "RawJSON"];
  Exit[code]
);
D2TermCount[x_] := Module[{v = Expand[x]},
  If[v === 0, 0, If[Head[v] === Plus, Length[List @@ v], 1]]
];
D2DimRegQ[s_Symbol] := MemberQ[
  {FromCharacterCode[1013], FromCharacterCode[949]}, D2Name[s]
];
D2DimRegQ[_] := False;
D2PoleResidue[expression_] := Module[
  {poles, unsupported, tagged, remaining, residue, higher},
  poles = Cases[expression,
    (HoldPattern[Power[s_Symbol, -1]] /; D2DimRegQ[s]) :> D2Name[s], Infinity];
  unsupported = Cases[expression,
    (HoldPattern[Power[s_Symbol, power_Integer]] /;
      D2DimRegQ[s] && power =!= -1) :> D2Text[Power[s, power]], Infinity];
  If[unsupported =!= {}, Return[<|
    "Status" -> "UnsupportedRegulatorPowers", "Preview" -> Take[unsupported, UpTo[3]]|>]];
  tagged = expression /. (
    HoldPattern[Power[s_Symbol, -1]] /; D2DimRegQ[s]) :> D2PoleMarker;
  remaining = DeleteDuplicates @ Cases[tagged,
    (s_Symbol /; D2DimRegQ[s]) :> D2Name[s], Infinity];
  If[remaining =!= {}, Return[<|
    "Status" -> "UnresolvedDimensionalRegulator", "Preview" -> remaining|>]];
  residue = Coefficient[Expand[tagged], D2PoleMarker, 1];
  higher = !FreeQ[residue, D2PoleMarker] ||
    Coefficient[Expand[tagged], D2PoleMarker, 2] =!= 0;
  If[higher, Return[<|"Status" -> "HigherOrderPole"|>]];
  <|"Status" -> "Separated", "PoleOccurrences" -> Length[poles],
    "Residue" -> residue|>
];
D2ExpandFiniteSums[expression_] := FixedPoint[
  Function[x, x /. HoldPattern[
    Sum[body_, {r_Symbol, lo_Integer, hi_Integer}]] :>
      If[hi < lo, 0, Total[Table[body /. r -> k, {k, lo, hi}]]]],
  expression, 12
];
D2ReduceNF1[expression_] := Module[{x},
  x = expression /. (s_Symbol /; D2Name[s] === "NF") :> 1;
  x = D2ExpandFiniteSums[x];
  x /. (HoldPattern[Matchete`Index[label_, rep_]] /;
      D2Name[rep] === "NFlavor") :> Matchete`Index[1, rep]
];

If[!FileExistsQ[candidatePath] || !FileExistsQ[flavourReportPath],
  D2Fail[4, "Saved provisional C5 or flavour audit missing"]
];
flavour = Quiet @ Check[Import[flavourReportPath, "RawJSON"], $Failed];
If[!AssociationQ[flavour] ||
  Lookup[flavour, "Status", ""] =!= "C5FlavourCandidateExported" ||
  !TrueQ[Lookup[flavour, "RunningSymmetricModuloSummationDummies", False]] ||
  !TrueQ[Lookup[flavour, "TotalCandidateLinearInHbar", False]],
  D2Fail[5, "Saved flavour audit is not certified for this provisional candidate"]
];
Print["[C5 pole] Loading Matchete and candidate WXF; no Match..."];
If[CheckAbort[UsingFrontEnd[Needs["Matchete`"]; True], $Aborted] =!= True,
  D2Fail[6, "Matchete did not load"]
];
candidate = Quiet @ Check[Import[candidatePath, "WXF"], $Failed];
If[candidate === $Failed || candidate === $Aborted || candidate === 0,
  D2Fail[7, "Candidate WXF is invalid"]
];
If[FreeQ[candidate, hbar], D2Fail[8, "Candidate has no explicit loop factor"]];
If[!FreeQ[candidate, D2PoleMarker],
  D2Fail[9, "Reserved pole marker already occurs in saved C5"]
];

Print["[C5 pole] Extracting exact pole before numeric substitutions..."];
pole = TimeConstrained[
  CheckAbort[D2PoleResidue[candidate], $Aborted], 90, $Aborted];
If[pole === $Aborted || !AssociationQ[pole] ||
    Lookup[pole, "Status", ""] =!= "Separated",
  D2Fail[10, "Could not isolate exact dimensional pole", <|
    "PoleStatus" -> If[AssociationQ[pole], Lookup[pole, "Status", ""], D2Text[pole]],
    "Details" -> If[AssociationQ[pole], Lookup[pole, "Preview", {}], {}]|>]
];
residue = pole["Residue"];
(* This is a strictly algebraic check with exact coefficients and symbolic
   masses/Yukawas; no Rationalize[] of approximate benchmark values. *)
resNF1 = TimeConstrained[
  CheckAbort[Expand[D2ReduceNF1[residue]], $Aborted], 90, $Aborted];
If[resNF1 === $Aborted, D2Fail[11, "NF=1 exact contraction timed out"]];
remainingSums = Cases[resNF1, HoldPattern[Sum[___]], Infinity];
remainingNF = Cases[resNF1,
  (s_Symbol /; D2Name[s] === "NF") :> D2Text[s], Infinity];
remainingNFlavor = Cases[resNF1,
  (HoldPattern[Matchete`Index[label_, rep_]] /;
    D2Name[rep] === "NFlavor" && !IntegerQ[label]) :> D2Text[label], Infinity];
If[remainingSums =!= {} || remainingNF =!= {} || remainingNFlavor =!= {},
  D2Fail[12, "Heavy-fermion flavour contraction remains incomplete", <|
    "RemainingSumCount" -> Length[remainingSums],
    "RemainingNFCount" -> Length[remainingNF],
    "RemainingNFlavorIndexCount" -> Length[remainingNFlavor]|>]
];
zeroByExpand = TrueQ[resNF1 === 0];
reduced = If[zeroByExpand, 0, TimeConstrained[
  CheckAbort[Quiet[Cancel[Together[resNF1]]], $Aborted], 90, $Aborted]];
reductionAborted = (reduced === $Aborted);
zeroByTogether = TrueQ[reduced === 0];
If[reductionAborted,
  reduced = resNF1;
  zeroByTogether = False
];
(* Do not call an expression with a residual pole renormalized. *)
Export[residuePath, reduced, "WXF"];
Export[residueTextPath, D2Text[reduced], "Text"];
bench = If[FileExistsQ[benchmarkReportPath],
  Quiet @ Check[Import[benchmarkReportPath, "RawJSON"], <||>], <||>];
If[!AssociationQ[bench], bench = <||>];
status = If[zeroByExpand || zeroByTogether,
  "ExactNF1PoleCancellationVerified", "ExactNF1PoleCancellationNotEstablished"];
report = <|
  "Status" -> status,
  "NoNewMatchPerformed" -> True,
  "SavedCandidateUnmodified" -> True,
  "NF" -> nf,
  "OriginalSimplePoleOccurrences" -> Lookup[pole, "PoleOccurrences", 0],
  "RawResidueTermCount" -> D2TermCount[residue],
  "AfterNF1ExpansionTermCount" -> D2TermCount[resNF1],
  "ReducedResidueTermCount" -> D2TermCount[reduced],
  "ExactZeroByExpand" -> zeroByExpand,
  "ExactZeroByTogether" -> zeroByTogether,
  "ExactPoleCancellationVerifiedForNF1" -> (zeroByExpand || zeroByTogether),
  "ReductionTimedOut" -> reductionAborted,
  "ResidualPreview" -> StringTake[D2Text[reduced], UpTo[1600]],
  "ResidualWXF" -> residuePath,
  "ResidualInputForm" -> residueTextPath,
  "PreviousNumericPoleMaxAbs_GeVinv" -> Lookup[bench, "C5PoleResidueMaxAbs_GeVinv", Missing["Unknown"]],
  "RenormalizationCertified" -> False,
  "PhysicalNormalizationIndependentlyCertified" -> False
|>;
Export[reportPath, D2JSONSafe[report], "RawJSON"];
Print["[C5 pole] ", status, "; exact residual terms: ", D2TermCount[reduced],
  "; report: ", reportPath];
Exit[0];
