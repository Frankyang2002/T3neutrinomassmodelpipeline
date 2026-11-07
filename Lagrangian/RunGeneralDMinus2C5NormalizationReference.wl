(*
  General D(-2): independent normalization reference for the PROVISIONAL C5.

  Read-only.  Loads the Matchete SMEFT_Warsaw model to inspect its Q_llHH
  normalization with the project's existing ExtractWeinbergCoefficient.
  Separately verifies the explicit epsilon/Higgs-vev factor v^2/2.

  Does NOT rerun matching; does NOT replace CG tensors by unity in the
  physical EFT; does NOT change any provisional coefficient or mass relation.

  Usage:
    wolframscript -file Lagrangian/RunGeneralDMinus2C5NormalizationReference.wl \
      FIRST_THRESHOLD_OUTPUT_DIR
*)
ClearAll["Global`*"];
args = Rest[$ScriptCommandLine];
If[Length[args] =!= 1,
  Print["Usage: FIRST_THRESHOLD_OUTPUT_DIR"];
  Exit[2]
];
outputDir = ExpandFileName[args[[1]]];
profileDir = FileNameJoin[{outputDir, "data", "general_d_m2_scalar_threshold", "profile_tree"}];
If[!DirectoryQ[profileDir],
  Print["[C5 normalization] Missing profile directory: ", profileDir];
  Exit[3]
];
reportPath = FileNameJoin[{profileDir, "profile_c5_normalization_reference_result.json"}];
refTermsPath = FileNameJoin[{profileDir, "profile_c5_normalization_warsaw_reference_terms.txt"}];
refCoeffPath = FileNameJoin[{profileDir, "profile_c5_normalization_warsaw_reference_coefficient.txt"}];
scriptDir = DirectoryName[ExpandFileName[$InputFileName]];

ClearAll[D2TextNorm, D2SafeNorm, D2FailNorm, D2TermsNorm, D2DegreeNorm,
  D2FieldQNorm, D2EpsilonNorm, D2HiggsNorm, D2UnitWilsonNorm];
D2TextNorm[x_] := ToString[x, InputForm];
D2SafeNorm[x_Association] := Association @ KeyValueMap[
  ToString[#1] -> D2SafeNorm[#2] &, x];
D2SafeNorm[x_List] := D2SafeNorm /@ x;
D2SafeNorm[x_String] := x;
D2SafeNorm[x_Integer] := x;
D2SafeNorm[True] := True;
D2SafeNorm[False] := False;
D2SafeNorm[Null] := Null;
D2SafeNorm[x_Real] := N[x, MachinePrecision];
D2SafeNorm[x_] := D2TextNorm[x];
D2FailNorm[code_Integer, reason_String, extra_: <||>] := Module[{record},
  Print["[C5 normalization] ERROR: ", reason];
  record = Join[<|"Status" -> "Failed", "Reason" -> reason,
     "ExitCode" -> code, "NoNewMatching" -> True,
     "C5NormalizationIndependentlyCertified" -> False|>, extra];
  Export[reportPath, D2SafeNorm[record], "RawJSON"];
  Exit[code]
];
D2TermsNorm[x_] := Module[{expanded = Expand[x]},
  If[expanded === 0, {},
    If[Head[expanded] === Plus, List @@ expanded, {expanded}]]
];
D2FieldQNorm[Matchete`Field[name_Symbol, ___], target_String] :=
  SymbolName[name] === target;
D2FieldQNorm[_, _String] := False;
D2DegreeNorm[x_ /; D2FieldQNorm[x, "H"], "H"] := 1;
D2DegreeNorm[x_ /; D2FieldQNorm[x, "l"], "l"] := 1;
D2DegreeNorm[Power[x_, n_Integer?Positive], target_String] :=
  n D2DegreeNorm[x, target];
D2DegreeNorm[x_?AtomQ, _String] := 0;
D2DegreeNorm[x_, target_String] := Total[
  D2DegreeNorm[#, target] & /@ (List @@ x)];
D2UnitWilsonNorm[x_] := x /. {
  HoldPattern[Matchete`Bar[Matchete`Coupling[___]]] :> 1,
  HoldPattern[Matchete`Coupling[___]] :> 1
};

priorPath = FileNameJoin[{profileDir, "profile_c5_flavour_candidate_result.json"}];
candidatePath = FileNameJoin[{profileDir, "profile_c5_flavour_candidate_total_provisional.wxf"}];
If[!FileExistsQ[priorPath] || !FileExistsQ[candidatePath],
  D2FailNorm[4, "Provisional flavour candidate or its audit is missing"]
];
prior = Quiet @ Check[Import[priorPath, "RawJSON"], $Failed];
If[!AssociationQ[prior] ||
   Lookup[prior, "Status", ""] =!= "C5FlavourCandidateExported" ||
   Lookup[prior, "RunningSymmetricModuloSummationDummies", False] =!= True ||
   Lookup[prior, "TotalCandidateLinearInHbar", False] =!= True ||
   Lookup[prior, "FixedOrderC5CandidateAssembled", False] =!= True,
  D2FailNorm[5, "Required flavour and fixed-order audits are not passed"]
];

Print["[C5 normalization] Loading official Warsaw reference; no matching invoked..."];
If[CheckAbort[UsingFrontEnd[Needs["Matchete`"]; True], $Aborted] =!= True,
  D2FailNorm[6, "Matchete could not be loaded"]
];
Get[FileNameJoin[{scriptDir, "RunMatching.wl"}]];
If[DownValues[ExtractWeinbergCoefficient] === {},
  D2FailNorm[7, "Project Weinberg extractor is unavailable"]
];

(* Do not evaluate or re-register any D(-2) BSM fields.  The Warsaw basis
   is independently loaded in a fresh Matchete kernel. *)
warsaw = CheckAbort[UsingFrontEnd[LoadModel["SMEFT_Warsaw"]], $Aborted];
If[MemberQ[{$Aborted, $Failed, Null}, warsaw] ||
   !FreeQ[warsaw, LoadModel],
  D2FailNorm[8, "Matchete Warsaw reference model did not load",
    <|"ReturnHead" -> D2TextNorm[Head[warsaw]]|>]
];

(* Strict operator-content selection, not merely an unrestricted GammaCC
   substring match; other SMEFT operators may also contain charge conjugation. *)
warsawTerms = D2TermsNorm[warsaw];
warsawWeinberg = Select[warsawTerms,
  !FreeQ[#, Matchete`GammaCC] &&
    D2DegreeNorm[#, "H"] === 2 && D2DegreeNorm[#, "l"] === 2 &
];
Export[refTermsPath, StringRiffle[D2TextNorm /@ warsawWeinberg, "\n\n"], "Text"];
If[warsawWeinberg === {},
  D2FailNorm[9, "No explicit LLHH GammaCC terms found in Warsaw reference",
    <|"WarsawTotalTermCount" -> Length[warsawTerms],
      "WarsawReferenceTermsPath" -> refTermsPath|>]
];

refExpr = Total[warsawWeinberg];
refExt = CheckAbort[ExtractWeinbergCoefficient[refExpr], $Aborted];
If[!AssociationQ[refExt] || Lookup[refExt, "Status", ""] =!= "Success",
  D2FailNorm[10, "Project extractor did not recognize Warsaw Q_llHH",
    <|"ReferenceExtractionStatus" -> If[AssociationQ[refExt],
      Lookup[refExt, "Status", "Unknown"], D2TextNorm[refExt]],
      "WarsawReferenceTermsPath" -> refTermsPath|>]
];

refHol = Lookup[refExt, "HolomorphicTerms", {}];
refHC = Lookup[refExt, "ConjugateTerms", {}];
refOrdered = Lookup[refExt, "Coefficient", $Failed];
If[refHol === {} || refHC === {} || refOrdered === $Failed ||
   MissingQ[refOrdered],
  D2FailNorm[11, "Warsaw reference cannot be separated into PL/PR coefficients"]
];
refNames = DeleteDuplicates @ Cases[refHol,
  Matchete`Coupling[name_Symbol, ___] :> SymbolName[name], Infinity];
refFactor = Quiet @ Check[Factor[D2UnitWilsonNorm[refOrdered]], $Failed];
Export[refCoeffPath, D2TextNorm[refOrdered], "Text"];

(* Linearity of the project's extractor under an independent coefficient
   rescaling of the reference operator: catches hidden additive terms or
   double symmetrization in the extraction routine. *)
scaledExt = CheckAbort[ExtractWeinbergCoefficient[3 refExpr], $Aborted];
scaledCoefficient = If[AssociationQ[scaledExt],
  Lookup[scaledExt, "Coefficient", $Failed], $Failed];
linearReference = scaledCoefficient =!= $Failed &&
  TrueQ[Expand[scaledCoefficient - 3 refOrdered] === 0];

(* Component-level EWSB verification of the Warsaw operator convention
   Q_llHH = epsilon_{in} epsilon_{jm} H^i H^j (l_p^n)^T C l_q^m.
   With H=(0,v/sqrt(2)) the neutrino-neutrino component is +v^2/2.
   For L \supset +(1/2) C5 Q + h.c. and the conventional Majorana
   mass term -(1/2) nu^T C m_nu nu + h.c., m_nu=-(v^2/2) C5.
   This does not establish the relative normalization of arbitrary CGs
   in the D(-2) matched EFT; that requires explicit tensor projection. *)
D2EpsilonNorm[i_Integer, j_Integer] := Which[
  i === 1 && j === 2, 1,
  i === 2 && j === 1, -1,
  True, 0
];
D2HiggsNorm[i_Integer] := If[i === 2, vD2Reference/Sqrt[2], 0];
ewsbNN = Simplify[Sum[
  D2EpsilonNorm[i, 1] D2EpsilonNorm[j, 1] *
    D2HiggsNorm[i] D2HiggsNorm[j],
  {i, 1, 2}, {j, 1, 2}]];
ewsbCorrect = TrueQ[Simplify[ewsbNN - vD2Reference^2/2] === 0];

refFactorNumeric = NumberQ[refFactor] && refFactor =!= 0;
report = <|
  "Status" -> "WarsawNormalizationReferenceAudited",
  "NoNewMatching" -> True,
  "ProvisionalCandidateUnchanged" -> True,
  "PriorFlavourCandidateStatus" -> prior["Status"],
  "WarsawModel" -> "SMEFT_Warsaw",
  "WarsawTotalTermCount" -> Length[warsawTerms],
  "WarsawWeinbergTermCount" -> Length[warsawWeinberg],
  "WarsawHolomorphicTermCount" -> Length[refHol],
  "WarsawConjugateTermCount" -> Length[refHC],
  "WarsawReferenceCouplingNames" -> refNames,
  "ProjectExtractorRecognizesWarsawWeinbergOperator" -> True,
  "ProjectExtractorLinearOnWarsawReference" -> linearReference,
  "WarsawOrderedCoefficientWithReferenceCouplingsSetToOne" -> D2TextNorm[refFactor],
  "WarsawOrderedCoefficientUnitWilsonFactorNumeric" -> refFactorNumeric,
  "WarsawSymmetricUnitWilsonFactorAssumingEqualPQandQP" ->
    If[refFactorNumeric, D2TextNorm[2 refFactor], "NotEstablished"],
  "WarsawOperatorNeutrinoVEVComponentIsPlusVSquaredOverTwo" -> ewsbCorrect,
  "MajoranaMassConvention" -> "For L+(1/2)C5 Q_llHH+h.c. and Lmass=-(1/2)nu^T C m nu+h.c., m=-v^2 C5/2",
  "WarsawReferenceTermsPath" -> refTermsPath,
  "WarsawReferenceOrderedCoefficientPath" -> refCoeffPath,
  "FullDMinus2SU2CGNormalizationCertified" -> False,
  "FullPhysicalC5NormalizationCertified" -> False,
  "PhysicalMassRelationsImposed" -> False,
  "NumericalBenchmarkPerformed" -> False
|>;
Export[reportPath, D2SafeNorm[report], "RawJSON"];
Print["[C5 normalization] Warsaw reference: ", Length[warsawWeinberg],
  " LLHH terms; holomorphic: ", Length[refHol],
  "; HC: ", Length[refHC],
  "; reference factor: ", D2TextNorm[refFactor],
  "; extractor linear: ", linearReference,
  "; EWSB: ", ewsbCorrect];
If[!TrueQ[linearReference] || !TrueQ[ewsbCorrect],
  Exit[20]
];
Exit[0];
