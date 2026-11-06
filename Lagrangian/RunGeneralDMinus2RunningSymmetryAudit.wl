(*
  Read-only D(-2) direct-running flavour-symmetry audit.

  Purpose: distinguish an actual p<->q asymmetry from an insufficiently
  simplified/representation-dependent comparison. NEVER modify or symmetrize
  the original direct-running coefficient. No new Match[] is performed.

  Usage:
    wolframscript -file Lagrangian/RunGeneralDMinus2RunningSymmetryAudit.wl \
      FIRST_THRESHOLD_OUTPUT_DIR DIRECT_RUNNING_INSERTION.wl
*)
ClearAll["Global`*"];
args = Rest[$ScriptCommandLine];
If[Length[args] =!= 2,
  Print["Usage: FIRST_THRESHOLD_OUTPUT_DIR DIRECT_RUNNING_INSERTION.wl"];
  Exit[2]
];
outputDir = ExpandFileName[args[[1]]];
insertionPath = ExpandFileName[args[[2]]];
profileDir = FileNameJoin[{outputDir, "data", "general_d_m2_scalar_threshold", "profile_tree"}];
If[!DirectoryQ[profileDir],
  CreateDirectory[profileDir, CreateIntermediateDirectories -> True]
];
root = FileNameJoin[{profileDir, "profile_c5_running_symmetry_audit"}];
reportPath = root <> "_result.json";
logPrefix = "[running symmetry] ";
ClearAll[D2Text, D2Fail, D2Terms, D2Export, D2SymbolName];
D2Text[x_] := ToString[x, InputForm];
D2Terms[x_] := Module[{v = Expand[x]},
  If[v === 0, {}, If[Head[v] === Plus, List @@ v, {v}]]
];
D2Export[suffix_String, expression_] := Module[{path = root <> "_" <> suffix},
  Export[path <> ".wxf", expression, "WXF"];
  Export[path <> ".txt", D2Text[expression], "Text"];
  path
];
D2Fail[code_Integer, reason_String] := (
  Print[logPrefix, "ERROR: ", reason];
  Export[reportPath, <|
    "Status" -> "Failed", "Reason" -> reason, "ExitCode" -> code,
    "NoNewMatchPerformed" -> True,
    "DirectRunningSourceUnmodified" -> True
  |>, "RawJSON"];
  Exit[code]
);

If[!FileExistsQ[insertionPath], D2Fail[3, "Direct-running insertion file missing"]];
Print[logPrefix, "Loading Matchete and direct-running definition; no Match..."];
loaded = CheckAbort[UsingFrontEnd[Needs["Matchete`"] ; True], $Aborted];
If[loaded =!= True, D2Fail[4, "Matchete failed to load"]];
getResult = CheckAbort[Get[insertionPath]; True, $Aborted];
If[getResult =!= True, D2Fail[5, "Direct-running insertion failed to load"]];
If[DownValues[EFT1D2WilsonRunningWeinbergCoefficient] === {} ||
   !ValueQ[EFT1D2WilsonRunningHeavyInsertion] ||
   EFT1D2WilsonRunningHeavyInsertion =!= 0,
  D2Fail[6, "Direct-running generator contract not satisfied"]
];

p = Matchete`Index[T3D2ExternalFlavorP, Flavor];
q = Matchete`Index[T3D2ExternalFlavorQ, Flavor];
runPQ = CheckAbort[EFT1D2WilsonRunningWeinbergCoefficient[p, q], $Aborted];
runQP = CheckAbort[EFT1D2WilsonRunningWeinbergCoefficient[q, p], $Aborted];
If[MemberQ[{$Failed, $Aborted}, runPQ] ||
   MemberQ[{$Failed, $Aborted}, runQP],
  D2Fail[7, "Running function failed on external flavour indices"]
];
If[FreeQ[runPQ, Log] || !FreeQ[runPQ, hbar] ||
   FreeQ[runQP, Log] || !FreeQ[runQP, hbar],
  D2Fail[8, "Running expression lacks Log or unexpectedly contains hbar"]
];

(* Fresh calls generate r$... symbols for bound Sum indices.  Equality of
   flavour tensors requires alpha-equivalence of those local binders.
   Canonicalize ONLY a single-index finite Sum binder: do not rename free
   flavour, gauge or Lorentz indices, and do not alter the source generator.
   Inactive[Sum] prevents Mathematica from re-creating fresh binders. *)
ClearAll[D2AlphaDummy, D2AlphaCanon, D2ActiveSums, D2InertSums];
D2AlphaCanon[x_] := x /. {
  HoldPattern[Sum[body_, {r_Symbol, lo_, hi_}]] :>
    Inactive[Sum][body /. r -> D2AlphaDummy,
      {D2AlphaDummy, lo, hi}]
};
D2ActiveSums[x_] := Length @ Cases[{x}, _Sum, Infinity];
D2InertSums[x_] := Length @ Cases[{x},
  HoldPattern[Inactive[Sum][_, {_, _, _}]], Infinity];
If[!FreeQ[runPQ, D2AlphaDummy] || !FreeQ[runQP, D2AlphaDummy],
  D2Fail[9, "Chosen alpha-equivalence dummy is already present in source"]
];
alphaPQ = D2AlphaCanon[runPQ];
alphaQP = D2AlphaCanon[runQP];
originalSums = {D2ActiveSums[runPQ], D2ActiveSums[runQP]};
normalizedSums = {D2InertSums[alphaPQ], D2InertSums[alphaQP]};
If[Min[originalSums] < 1 || originalSums =!= normalizedSums ||
   D2ActiveSums[alphaPQ] =!= 0 || D2ActiveSums[alphaQP] =!= 0,
  D2Fail[10, "Unrecognized Sum structure; refusing partial alpha normalization"]
];
alphaResidual = Expand[alphaPQ - alphaQP];
alphaEqual = TrueQ[alphaResidual === 0];
alphaResidualRoot = D2Export["alpha_normalized_pq_minus_qp", alphaResidual];
Print[logPrefix, "Alpha-normalized p-q symmetry: ", alphaEqual,
  "; finite sums checked: ", originalSums];

residual = Expand[runPQ - runQP];
residualTerms = D2Terms[residual];
(* Keep the simple exact polynomial check distinct from further symbolic
   simplification; neither should silently change the matching result. *)
expandedEqual = TrueQ[residual === 0];
simplified = TimeConstrained[
  CheckAbort[Quiet[Simplify[residual]], $Aborted], 15, $Aborted
];
simplifiedZero = TrueQ[simplified === 0];
simplifyReturned = simplified =!= $Aborted && simplified =!= $Failed;

pqRoot = D2Export["pq_as_defined", runPQ];
qpRoot = D2Export["qp_swapped", runQP];
residualRoot = D2Export["pq_minus_qp", residual];
(* Write both halves ONLY as audit metadata, not as a replacement for the
   supplied definition.  The required normalization depends on whether
   the generator denotes C5 or ordered A5. *)
symmetricPartRoot = D2Export["symmetric_part_NOT_APPLIED", (runPQ + runQP)/2];
antiPartRoot = D2Export["antisymmetric_part_NOT_APPLIED", residual/2];
GetDefinitionPath = root <> "_generator_definition.txt";
Export[GetDefinitionPath,
  D2Text[DownValues[EFT1D2WilsonRunningWeinbergCoefficient]], "Text"];

Print[logPrefix, "Expanded p-q symmetry: ", expandedEqual,
  "; simplified symmetry: ", simplifiedZero,
  "; residual terms: ", Length[residualTerms]];

report = <|
  "Status" -> "RunningSymmetryAudited",
  "NoNewMatchPerformed" -> True,
  "DirectRunningSourceUnmodified" -> True,
  "SourceInsertionFile" -> insertionPath,
  "DefinitionFile" -> GetDefinitionPath,
  "InputHasLog" -> True,
  "InputHasExplicitHbar" -> False,
  "InputHeavyInsertionIsZero" -> True,
  "SymmetricByExpand" -> expandedEqual,
  "SimplifyCompleted" -> simplifyReturned,
  "SymmetricBySimplify" -> simplifiedZero,
  "SymmetricModuloSummationDummies" -> alphaEqual,
  "AlphaNormalizationFullyCoveredSums" -> True,
  "BoundSumCountsPQandQP" -> originalSums,
  "AlphaNormalizedResidualTermCount" -> Length[D2Terms[alphaResidual]],
  "AlphaNormalizedResidualRoot" -> alphaResidualRoot,
  "ResidualTermCount" -> Length[residualTerms],
  "ResidualTermPreviews" ->
    (StringTake[D2Text[#], UpTo[900]] & /@ Take[residualTerms, UpTo[10]]),
  "PQExpressionRoot" -> pqRoot,
  "QPExpressionRoot" -> qpRoot,
  "ResidualExpressionRoot" -> residualRoot,
  "SymmetricPartDiagnosticRoot" -> symmetricPartRoot,
  "AntisymmetricPartDiagnosticRoot" -> antiPartRoot,
  "RunningFlavourSymmetryCertifiedModuloDummyRenaming" -> alphaEqual,
  "RunningNormalizationCertified" -> False,
  "FixedOrderSMEFTCertified" -> False
|>;
Export[reportPath, report, "RawJSON"];
Print[logPrefix, "Audited; report: ", reportPath];
Exit[0];
