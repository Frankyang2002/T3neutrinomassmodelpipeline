(*
  Provisional numerical D(-2) C5 benchmark from the SAVED symmetric C5.
  Read-only inputs; NO Matchete Match[] calls.  The reference Warsaw
  normalization has passed, but the full model CG normalization is NOT
  independently certified.  All parameter choices are toy values.

  Invocation:
    wolframscript -file Lagrangian/RunGeneralDMinus2C5RoughBenchmark.wl \
      OUTPUT_DIR benchmarks/general_d_m2_rough_benchmark.json
*)
ClearAll["Global`*"];
args = Rest[$ScriptCommandLine];
If[Length[args] =!= 2,
  Print["Usage: FIRST_THRESHOLD_OUTPUT_DIR BENCHMARK_JSON"];
  Exit[2]
];
outputDir = ExpandFileName[args[[1]]];
configPath = ExpandFileName[args[[2]]];
profileDir = FileNameJoin[{outputDir, "data", "general_d_m2_scalar_threshold", "profile_tree"}];
candidatePath = FileNameJoin[{profileDir, "profile_c5_flavour_candidate_total_provisional.wxf"}];
flavourReportPath = FileNameJoin[{profileDir, "profile_c5_flavour_candidate_result.json"}];
normReportPath = FileNameJoin[{profileDir, "profile_c5_normalization_reference_result.json"}];
reportPath = FileNameJoin[{profileDir, "profile_c5_rough_benchmark_result.json"}];
inventoryPath = FileNameJoin[{profileDir, "profile_c5_rough_benchmark_couplings.json"}];
c5RealPath = FileNameJoin[{profileDir, "profile_c5_rough_benchmark_C5_real_GeVinv.csv"}];
c5ImagPath = FileNameJoin[{profileDir, "profile_c5_rough_benchmark_C5_imag_GeVinv.csv"}];
mnuRealPath = FileNameJoin[{profileDir, "profile_c5_rough_benchmark_mnu_real_eV.csv"}];
mnuImagPath = FileNameJoin[{profileDir, "profile_c5_rough_benchmark_mnu_imag_eV.csv"}];
poleRealPath = FileNameJoin[{profileDir, "profile_c5_rough_benchmark_C5_pole_residue_real.csv"}];
poleImagPath = FileNameJoin[{profileDir, "profile_c5_rough_benchmark_C5_pole_residue_imag.csv"}];

ClearAll[D2Text, D2SafeJSON, D2Fail, D2Name, D2Label, D2NumericalCoupling,
  D2Yukawa, D2Evaluate, D2ExpandFiniteSums, D2EvaluateEntry,
  D2MuBar2Q, D2DimRegQ, D2SplitPole, D2ResidualNames, D2FiniteNumberQ];
D2Text[x_] := ToString[x, InputForm];
D2Name[s_Symbol] := SymbolName[Unevaluated[s]];
D2Name[x_] := D2Text[x];
D2SafeJSON[x_Association] := Association @ KeyValueMap[ToString[#1] -> D2SafeJSON[#2] &, x];
D2SafeJSON[x_List] := D2SafeJSON /@ x;
D2SafeJSON[x_String] := x;
D2SafeJSON[x_Integer] := x;
D2SafeJSON[True] := True;
D2SafeJSON[False] := False;
D2SafeJSON[Null] := Null;
(* Matchete/Wolfram can print zero as 0.e-7, which RawJSON rejects. *)
D2SafeJSON[x_Real] := Module[{v = N[x, MachinePrecision]},
  If[TrueQ[v == 0], 0, v]
];
D2SafeJSON[x_] := D2Text[x];
D2Fail[code_Integer, reason_String, info_: <||>] := Module[{rec},
  Print["[C5 benchmark] ERROR: ", reason];
  rec = Join[<|"Status" -> "Failed", "Reason" -> reason,
    "ExitCode" -> code, "NumericalBenchmarkCertifiedPhysical" -> False,
    "NoNewMatchPerformed" -> True|>, info];
  Export[reportPath, D2SafeJSON[rec], "RawJSON"];
  Exit[code]
];

If[!AllTrue[{candidatePath, flavourReportPath, normReportPath, configPath}, FileExistsQ],
  Print["[C5 benchmark] Missing files: ",
    Select[{candidatePath, flavourReportPath, normReportPath, configPath}, !FileExistsQ[#] &]];
  D2Fail[3, "Required candidate, audits or configuration missing"]
];
flavour = Quiet @ Check[Import[flavourReportPath, "RawJSON"], $Failed];
norm = Quiet @ Check[Import[normReportPath, "RawJSON"], $Failed];
config = Quiet @ Check[Import[configPath, "RawJSON"], $Failed];
If[!AssociationQ[flavour] || !AssociationQ[norm] || !AssociationQ[config],
  D2Fail[4, "Could not read JSON prerequisites"]
];
If[Lookup[flavour, "Status", ""] =!= "C5FlavourCandidateExported" ||
   !TrueQ[Lookup[flavour, "RunningSymmetricModuloSummationDummies", False]] ||
   !TrueQ[Lookup[flavour, "FixedOrderC5CandidateAssembled", False]] ||
   Lookup[norm, "Status", ""] =!= "WarsawNormalizationReferenceAudited" ||
   !TrueQ[Lookup[norm, "ProjectExtractorRecognizesWarsawWeinbergOperator", False]] ||
   !TrueQ[Lookup[norm, "ProjectExtractorLinearOnWarsawReference", False]] ||
   !TrueQ[Lookup[norm, "WarsawOperatorNeutrinoVEVComponentIsPlusVSquaredOverTwo", False]] ||
   ToString[Lookup[norm, "WarsawOrderedCoefficientWithReferenceCouplingsSetToOne", ""]] =!= "1",
  D2Fail[5, "Flavour or normalization-reference prerequisite not passed"]
];

massNames = {"MF", "MS1", "MS2", "T3S2RMass", "T3S2IMass", "MuScalar"};
mf = Lookup[config, "MF_GeV", Missing["MF"]];
ms1 = Lookup[config, "MS1_GeV", Missing["MS1"]];
ms2 = Lookup[config, "MS2_GeV", Missing["MS2"]];
kappa = Lookup[config, "AbsKappaS_GeV2", Missing["Kappa"]];
mu = Lookup[config, "MuScalar_GeV", Missing["Mu"]];
vev = Lookup[config, "Vev_GeV", Missing["Vev"]];
nf = Lookup[config, "NF", Missing["NF"]];
yukawas = Lookup[config, "Yukawas", <||>];
scalarCouplings = Lookup[config, "ScalarCouplings", <||>];
dimless = Lookup[config, "DimensionlessParameters", <||>];
otherLambda = Lookup[config, "OtherLambdaDefault", Missing["Lambda"]];
poleHandling = Lookup[config, "PoleHandling", "Reject"];
If[!MemberQ[{"Reject", "FinitePartDiagnostic"}, poleHandling],
  D2Fail[6, "PoleHandling must be Reject or FinitePartDiagnostic"]
];
If[!And @@ (NumericQ /@ {mf, ms1, ms2, kappa, mu, vev, otherLambda}) ||
   !IntegerQ[nf] || nf < 1 || nf > 6 ||
   Min[N[{mf, ms1, ms2, mu, vev}]] <= 0 || kappa < 0 ||
   ms2^2 <= 2 kappa ||
   !AssociationQ[yukawas] || !AssociationQ[scalarCouplings] ||
   !AssociationQ[dimless],
  D2Fail[6, "Invalid toy parameters or nonpositive scalar mass squared"]
];
If[!AllTrue[{"y1", "y2", "y2Conj"},
  Function[n, KeyExistsQ[yukawas, n] && MatrixQ[yukawas[n], NumericQ] &&
    Dimensions[yukawas[n]] === {3, nf}]],
  D2Fail[7, "Each independent Yukawa must be a numeric 3 x NF matrix"]
];
mr = Sqrt[ms2^2 - 2 kappa];
mi = Sqrt[ms2^2 + 2 kappa];
Print["[C5 benchmark] Masses GeV: MF=", mf, ", MS1=", ms1,
  ", mR=", N[mr], ", mI=", N[mi], "; NF=", nf];

Print["[C5 benchmark] Loading candidate WXF; NO Match called..."];
load = CheckAbort[UsingFrontEnd[Needs["Matchete`"]; True], $Aborted];
If[load =!= True, D2Fail[8, "Matchete could not load"]];
candidate = Quiet @ Check[Import[candidatePath, "WXF"], $Failed];
If[candidate === $Failed || candidate === $Aborted || candidate === 0,
  D2Fail[9, "Cannot import nonzero provisional C5 candidate"]
];
(* Preflight inventory is exported before any attempted numerical evaluation.
   No scalar coupling or unknown constant is silently set to zero. *)
allCouplings = DeleteDuplicates @ Cases[candidate,
  HoldPattern[Matchete`Coupling[n_Symbol, idx_List, ord_]] :>
    <|"Name" -> D2Name[n], "IndexCount" -> Length[idx],
      "PowerCounting" -> D2Text[ord]|>, Infinity];
scalarNames = Sort @ DeleteDuplicates @ Cases[allCouplings,
  a_Association /; Lookup[a, "IndexCount", -1] === 0 :> a["Name"]];
flavNames = Sort @ DeleteDuplicates @ Cases[allCouplings,
  a_Association /; Lookup[a, "IndexCount", -1] > 0 :> a["Name"]];
Export[inventoryPath, <|
  "Couplings" -> allCouplings, "ScalarParameterNames" -> scalarNames,
  "IndexedCouplingNames" -> flavNames,
  "BenchmarkIsIllustrative" -> True, "Configuration" -> configPath|>, "RawJSON"];
Print["[C5 benchmark] Coupling names: ", scalarNames,
  "; indexed: ", flavNames];

(* Mass substitutions are physical second-threshold masses; the first
   threshold WXF is never modified.  The three mass symbols have order 0
   in the saved second-threshold candidate. *)
scalarValues = Join[
  <|"MF" -> N[mf, 20], "MS1" -> N[ms1, 20],
    "MS2" -> N[ms2, 20], "T3S2RMass" -> N[mr, 20],
    "T3S2IMass" -> N[mi, 20], "MuScalar" -> N[mu, 20]|>,
  Association @ KeyValueMap[#1 -> N[#2, 20] &, dimless],
  Association @ KeyValueMap[#1 -> N[#2, 20] &, scalarCouplings]
];
(* Catch an unsupported indexed coupling without fabricating a matrix. *)
unsupportedIndexed = Complement[flavNames, {"y1", "y2", "y2Conj"}];
If[unsupportedIndexed =!= {},
  D2Fail[10, "Other flavour-indexed couplings need explicit benchmark matrices",
    <|"UnsupportedIndexedCouplings" -> unsupportedIndexed,
      "InventoryPath" -> inventoryPath|>]
];
(* Dimensionless quartics not individually specified use the declared
   illustrative OtherLambdaDefault; ALL such assignments are recorded. *)
filledLambdaNames = Select[scalarNames,
  StringStartsQ[#, "lambda"] && !KeyExistsQ[scalarValues, #] &];
Scan[(scalarValues[#] = N[otherLambda, 20]) &, filledLambdaNames];
unknownScalar = Complement[scalarNames, Keys[scalarValues]];
If[unknownScalar =!= {},
  D2Fail[11, "Unassigned scalar parameters; add values to JSON config",
    <|"UnassignedScalarParameters" -> unknownScalar,
      "InventoryPath" -> inventoryPath|>]
];
Print["[C5 benchmark] Illustrative fallback quartics: ", filledLambdaNames];

D2Label[Matchete`Index[label_, _]] := label;
D2Label[Matchete`Bar[x_]] := D2Label[x];
D2Label[x_] := x;
D2Yukawa[name_String, p_Integer, r_Integer] := Module[{mat},
  mat = Lookup[yukawas, name, Missing["NoYukawa"]];
  If[!MatrixQ[mat, NumericQ] || p < 1 || p > 3 || r < 1 || r > nf,
    Indeterminate, N[mat[[p, r]], 20]]
];
D2NumericalCoupling[n_Symbol, idx_List, order_] := Module[{name, labels},
  name = D2Name[n];
  If[idx === {}, Return[Lookup[scalarValues, name, D2MissingValue[name]]]];
  labels = D2Label /@ idx;
  If[MemberQ[{"y1", "y2", "y2Conj"}, name] && Length[labels] === 2,
    Return[D2Yukawa[name, labels[[1]], labels[[2]]]]
  ];
  D2UnknownIndexed[name, labels]
];
(* The saved coefficient contains Sum[..., {r$..., 1, NF}].  In the
   previous benchmark the Yukawa lookup happened FIRST.  That left the
   heavy-fermion flavour index symbolic inside an unevaluated coupling,
   leading to nonnumeric entries or lengthy symbolic simplification.

   Resolve NF and hbar ONCE, then explicitly unroll finite Sum iterators,
   and ONLY THEN substitute Yukawa matrix entries.  Do not remove sums,
   project flavour indices, or assign any unknown coupling to zero. *)
D2ExpandFiniteSums[expr_] := FixedPoint[
  Function[x, x /. HoldPattern[Sum[body_, {r_Symbol, lo_Integer,
        hi_Integer}]] :> If[hi < lo, 0,
      Total[Table[body /. r -> k, {k, lo, hi}]]]],
  expr, 12
];

(* Matchete's renormalization scale is mu_bar^2, not MuScalar itself.
   InputForm can contain either Mubar2 or the Greek-mu name mu-bar2. *)
D2MuBar2Q[s_Symbol] := MemberQ[
  {"Mubar2", FromCharacterCode[956] <> "bar2"}, D2Name[s]];
D2MuBar2Q[_] := False;
(* The Greek lunate epsilon (U+03F5) denotes a dimensional pole.
   This must NEVER be mistaken for a Yukawa index or set numerically. *)
D2DimRegQ[s_Symbol] := MemberQ[
  {FromCharacterCode[1013], FromCharacterCode[949]}, D2Name[s]];
D2DimRegQ[_] := False;
muBarSymbols = DeleteDuplicates @ Cases[candidate,
  s_Symbol /; D2MuBar2Q[s] :> D2Text[s], Infinity];
dimRegSymbols = DeleteDuplicates @ Cases[candidate,
  s_Symbol /; D2DimRegQ[s] :> D2Text[s], Infinity];
Print["[C5 benchmark] Scale-squared symbols: ", muBarSymbols,
  "; dimensional regulator symbols: ", dimRegSymbols];
numericTemplate = candidate /. {
  s_Symbol /; D2Name[s] === "NF" :> nf,
  s_Symbol /; D2Name[s] === "hbar" :> 1/(16 Pi^2),
  s_Symbol /; D2MuBar2Q[s] :> N[mu^2, 20]
};
numericTemplate = D2ExpandFiniteSums[numericTemplate];
remainingSums = Cases[numericTemplate, HoldPattern[Sum[___]], Infinity];
remainingNF = Cases[numericTemplate,
  s_Symbol /; D2Name[s] === "NF" :> D2Text[s], Infinity];
Print["[C5 benchmark] Finite sums prepared: remaining Sum heads=",
  Length[remainingSums], "; remaining NF symbols=", Length[remainingNF]];
If[remainingSums =!= {} || remainingNF =!= {},
  D2Fail[12, "Finite heavy-fermion sums were not fully expanded",
    <|"RemainingSumCount" -> Length[remainingSums],
      "RemainingSumPreview" -> (StringTake[D2Text[#], UpTo[500]] & /@
        Take[remainingSums, UpTo[3]]),
      "RemainingNFSymbols" -> DeleteDuplicates[remainingNF]|>]
];

(* After summing explicit Sum binders, Matchete's inherited EFT can
   retain implicit NFlavor dummy labels (d$$1 etc.). For NF=1 and ONLY
   NF=1, all such typed labels unambiguously denote the sole generation.
   For NF>1 an explicit contraction implementation is required; never
   invent a flavour assignment. Ordinary external Flavor labels remain. *)
implicitNFlavorBefore = DeleteDuplicates @ Cases[numericTemplate,
  HoldPattern[Matchete`Index[label_, rep_]] /;
    D2Name[rep] === "NFlavor" :> D2Text[label], Infinity];
If[nf === 1,
  numericTemplate = numericTemplate /.
    (HoldPattern[Matchete`Index[label_, rep_]] /;
      D2Name[rep] === "NFlavor") :> Matchete`Index[1, rep]
];
implicitNFlavorAfter = DeleteDuplicates @ Cases[numericTemplate,
  HoldPattern[Matchete`Index[label_, rep_]] /;
    D2Name[rep] === "NFlavor" && !IntegerQ[label] :> D2Text[label], Infinity];
Print["[C5 benchmark] Unresolved typed NFlavor indices before/after: ",
  implicitNFlavorBefore, " / ", implicitNFlavorAfter];
If[implicitNFlavorAfter =!= {},
  D2Fail[12, "Uncontracted heavy-generation index; NF>1 needs explicit tensor contraction",
    <|"NF" -> nf, "UnresolvedNFlavorLabels" -> implicitNFlavorAfter|>]
];

(* Explicitly separate Laurent pole and finite part.  In the benchmark
   there are simple 1/epsilon poles from the one-loop coefficient.  A
   nonzero pole is physically UNRENORMALIZED; removing it is only allowed
   in the opt-in FinitePartDiagnostic mode, never silently. *)
D2SplitPole[expression_] := Module[
  {poles, unsupported, tagged, residue, finite, leftovers},
  poles = Cases[expression,
    HoldPattern[Power[s_Symbol, -1]] /; D2DimRegQ[s] :> D2Text[s], Infinity];
  unsupported = Cases[expression,
    HoldPattern[Power[s_Symbol, n_Integer]] /;
      D2DimRegQ[s] && n =!= -1 :> D2Text[Power[s, n]], Infinity];
  If[unsupported =!= {},
    Return[<|"Status" -> "UnsupportedRegulatorPower",
      "Preview" -> Take[unsupported, UpTo[3]]|>]
  ];
  tagged = expression /.
    (HoldPattern[Power[s_Symbol, -1]] /; D2DimRegQ[s]) :>
      D2RegulatorPoleMarker;
  leftovers = DeleteDuplicates @ Cases[tagged,
    s_Symbol /; D2DimRegQ[s] :> D2Text[s], Infinity];
  If[leftovers =!= {},
    Return[<|"Status" -> "RegulatorNotFullyIsolated",
      "Preview" -> leftovers|>]
  ];
  residue = Coefficient[Expand[tagged], D2RegulatorPoleMarker, 1];
  finite = Expand[tagged /. D2RegulatorPoleMarker -> 0];
  If[!FreeQ[residue, D2RegulatorPoleMarker],
    Return[<|"Status" -> "HigherOrderRegulatorPole", "Preview" -> {}|>]
  ];
  <|"Status" -> "Separated", "Finite" -> finite,
    "PoleResidue" -> residue, "SimplePoleOccurrences" -> Length[poles]|>
];

(* Every Yukawa lookup now receives two concrete integer indices. *)
D2Evaluate[expr_, p_Integer, q_Integer] := Module[{v, split},
  v = expr /. {
    HoldPattern[Matchete`Index[lab_Symbol, flav_]] /;
        D2Name[lab] === "T3D2ExternalFlavorP" :> Matchete`Index[p, flav],
    HoldPattern[Matchete`Index[lab_Symbol, flav_]] /;
        D2Name[lab] === "T3D2ExternalFlavorQ" :> Matchete`Index[q, flav]
  };
  (* Keep typed Yukawa indices in diagnostics in case an unsummed
     NFlavor/Flavor label survives a future candidate representation. *)
  $D2LastTypedYukawaIndices = Take[DeleteDuplicates @ Cases[v,
    HoldPattern[Matchete`Coupling[n_Symbol, idx_List, order_]] /;
      MemberQ[{"y1", "y2", "y2Conj"}, D2Name[n]] :>
        D2Text[idx], Infinity], UpTo[12]];
  v = v /. HoldPattern[Matchete`Coupling[n_Symbol, idx_List, order_]] :>
    D2NumericalCoupling[n, idx, order];
  v = v //. {
    HoldPattern[Matchete`Bar[z_?NumericQ]] :> Conjugate[z],
    HoldPattern[Matchete`Bar[Matchete`Bar[z_]]] :> z
  };
  $D2LastPreNumericEntry = v;
  split = D2SplitPole[v];
  If[Lookup[split, "Status", ""] =!= "Separated", Return[split]];
  split = Join[split, <|
    "Finite" -> Quiet @ CheckAbort[
      TimeConstrained[N[split["Finite"], 18], 25, $Aborted], $Aborted],
    "PoleResidue" -> Quiet @ CheckAbort[
      TimeConstrained[N[split["PoleResidue"], 18], 25, $Aborted], $Aborted]
  |>];
  split
];

D2EvaluateEntry[p_Integer, q_Integer] := Module[{entry, fin, residue},
  entry = D2Evaluate[numericTemplate, p, q];
  fin = Lookup[entry, "Finite", $Failed];
  residue = Lookup[entry, "PoleResidue", $Failed];
  Print["[C5 benchmark] C5[", p, ",", q, "] finite numeric=",
    NumericQ[fin], "; pole residue numeric=", NumericQ[residue],
    "; pole occurrences=", Lookup[entry, "SimplePoleOccurrences", -1]];
  If[!TrueQ[NumericQ[fin]] || !TrueQ[NumericQ[residue]] ||
     !FreeQ[{fin, residue},
       Indeterminate | ComplexInfinity | DirectedInfinity[_]],
    D2Fail[12, "Nonnumeric finite or pole C5 entry after coupling substitution",
      <|"BadMatrixEntry" -> {p, q},
        "PoleSplitStatus" -> Lookup[entry, "Status", "Unknown"],
        "PoleSplitPreview" -> Lookup[entry, "Preview", {}],
        "EntryPreview" -> StringTake[D2Text[fin], UpTo[1800]],
        "PoleResiduePreview" -> StringTake[D2Text[residue], UpTo[800]],
        "UnresolvedCouplingCount" -> Length[Cases[fin,
          HoldPattern[Matchete`Coupling[___]], Infinity]],
        "UnresolvedBarCount" -> Length[Cases[fin,
          HoldPattern[Matchete`Bar[___]], Infinity]],
        "UnresolvedSumCount" -> Length[Cases[fin,
          HoldPattern[Sum[___]], Infinity]],
        "UnknownYukawaCount" -> Length[Cases[fin,
          HoldPattern[D2Yukawa[___]], Infinity]],
        "UnknownScalarValues" -> DeleteDuplicates[Cases[fin,
          HoldPattern[D2MissingValue[n_]] :> D2Text[n], Infinity]],
        "EntryTimedOut" -> MemberQ[{fin, residue}, $Aborted],
        "PreNumericEntryPreview" -> StringTake[
          D2Text[$D2LastPreNumericEntry], UpTo[1800]],
        "TypedYukawaIndexPreview" -> $D2LastTypedYukawaIndices,
        "RemainingSymbolsInputForm" -> Take[
          DeleteDuplicates[Cases[fin,
            s_Symbol :> D2Text[s], Infinity]], UpTo[24]],
        "InventoryPath" -> inventoryPath|>]
  ];
  entry
];
entryGrid = Table[D2EvaluateEntry[p, q], {p, 1, 3}, {q, 1, 3}];
matrixC5 = Map[Lookup[#, "Finite"] &, entryGrid, {2}];
matrixPole = Map[Lookup[#, "PoleResidue"] &, entryGrid, {2}];
poleCount = Total[Flatten[Map[Lookup[#, "SimplePoleOccurrences", 0] &,
  entryGrid, {2}]]];
maxPoleAbs = Max[Abs[Flatten[N[matrixPole, 15]]]];
Print["[C5 benchmark] Dimensional-pole occurrences: ", poleCount,
  "; maximal C5 pole residue: ", maxPoleAbs];
If[!MatrixQ[matrixPole, NumericQ] ||
  !FreeQ[matrixPole, Indeterminate | ComplexInfinity | DirectedInfinity[_]],
  D2Fail[12, "Dimensional pole residues did not evaluate numerically"]
];
Export[poleRealPath, N[Re[matrixPole], 12], "CSV"];
Export[poleImagPath, N[Im[matrixPole], 12], "CSV"];
If[poleHandling === "Reject" && TrueQ[maxPoleAbs > 0],
  D2Fail[15, "Nonzero 1/epsilon pole remains in assembled C5; benchmark requires an explicit finite-part diagnostic choice or a renormalized coefficient",
    <|"DimRegPolePresent" -> True,
      "C5PoleResidueMaxAbs_GeVinv" -> N[maxPoleAbs, 15],
      "PoleOccurrenceCount" -> poleCount,
      "PoleHandling" -> poleHandling,
      "PoleResidueRealCSV" -> poleRealPath,
      "PoleResidueImagCSV" -> poleImagPath,
      "RenormalizationCertified" -> False|>]
];
(* Strong numeric gate: catches any unresolved Coupling, Matchete Bar,
   unevaluated Sum, loop function, flavour or hidden free symbol. *)
If[!MatrixQ[matrixC5, NumericQ] ||
  !FreeQ[matrixC5, Indeterminate | ComplexInfinity | DirectedInfinity[_]] ||
  AnyTrue[Flatten[matrixC5],
    (!NumericQ[N[Re[#], MachinePrecision]] ||
      !NumericQ[N[Im[#], MachinePrecision]]) &],
  D2Fail[12, "Candidate did not evaluate to a finite numeric 3x3 matrix",
    <|"EntryHeads" -> (D2Text[Head[#]] & /@ Flatten[matrixC5]),
      "EntryPreviews" -> (StringTake[D2Text[#], UpTo[1000]] & /@
        Take[Flatten[matrixC5], UpTo[4]]),
      "UnresolvedCouplingCount" -> Length[Cases[matrixC5,
        HoldPattern[Matchete`Coupling[___]], Infinity]],
      "UnresolvedBarCount" -> Length[Cases[matrixC5,
        HoldPattern[Matchete`Bar[___]], Infinity]],
      "UnresolvedSumCount" -> Length[Cases[matrixC5,
        HoldPattern[Sum[___]], Infinity]],
      "UnknownYukawaCount" -> Length[Cases[matrixC5,
        HoldPattern[D2Yukawa[___]], Infinity]],
      "UnknownScalarValues" -> DeleteDuplicates[Cases[matrixC5,
        HoldPattern[D2MissingValue[n_]] :> D2Text[n], Infinity]],
      "EntryAbortedCount" -> Count[Flatten[matrixC5], $Aborted],
      "InventoryPath" -> inventoryPath|>]
];
matrixC5 = N[matrixC5, 15];
symResidual = Max[Abs[Flatten[matrixC5 - Transpose[matrixC5]]]];
symScale = Max[10^-100, Max[Abs[Flatten[matrixC5]]]];
If[!TrueQ[N[symResidual/symScale] < 10^-8],
  D2Fail[13, "Numeric C5 is not symmetric",
    <|"RelativeSymmetryResidual" -> D2Text[symResidual/symScale]|>]
];
(* The reference convention is Lsupset+(1/2) C5 Q +h.c., therefore
   mnu=-(v^2/2)C5.  In general complex-symmetric cases the physical
   Takagi masses are singular values, NOT ordinary eigenvalues. *)
mnuEV = N[-(vev^2/2) matrixC5 10^9, 15];
massesEV = Sort[N[SingularValueList[mnuEV], 12]];
If[!VectorQ[massesEV, NumericQ] || Length[massesEV] =!= 3,
  D2Fail[14, "Takagi singular values could not be obtained"]
];
Export[c5RealPath, N[Re[matrixC5], 12], "CSV"];
Export[c5ImagPath, N[Im[matrixC5], 12], "CSV"];
Export[mnuRealPath, N[Re[mnuEV], 12], "CSV"];
Export[mnuImagPath, N[Im[mnuEV], 12], "CSV"];
report = <|
  "Status" -> If[poleHandling === "FinitePartDiagnostic",
    "ProvisionalFinitePartBenchmarkEvaluated", "ProvisionalNumericalBenchmarkEvaluated"],
  "NoNewMatchPerformed" -> True,
  "PhysicalC5NormalizationIndependentlyCertified" -> False,
  "NumericalBenchmarkCertifiedPhysical" -> False,
  "WarsawReferencePassed" -> True,
  "FlavourSymmetryCertified" -> True,
  "BenchmarkConfigPath" -> configPath,
  "CandidatePath" -> candidatePath,
  "NF" -> nf,
  "MS2_GeV" -> ms2,
  "AbsKappaS_GeV2" -> kappa,
  "T3S2RMass_GeV" -> N[mr, 15],
  "T3S2IMass_GeV" -> N[mi, 15],
  "OneOver16PiSquared" -> N[1/(16 Pi^2), 14],
  "PoleHandling" -> poleHandling,
  "DimRegPolePresent" -> TrueQ[maxPoleAbs > 0],
  "C5PoleResidueMaxAbs_GeVinv" -> N[maxPoleAbs, 15],
  "PoleOccurrenceCount" -> poleCount,
  "RenormalizationCertified" -> False,
  "FinitePartOnly" -> (poleHandling === "FinitePartDiagnostic"),
  "ImplicitNFlavorLabelsMappedToOne" -> (nf === 1 && implicitNFlavorBefore =!= {}),
  "MuBar2ParameterNames" -> muBarSymbols,
  "RegulatorNames" -> dimRegSymbols,
  "PoleResidueRealCSV" -> poleRealPath,
  "PoleResidueImagCSV" -> poleImagPath,
  "FallbackLambdaAssignments" -> filledLambdaNames,
  "ScalarCouplingAssignments" -> scalarValues,
  "C5RelativeSymmetryResidual" -> N[symResidual/symScale, 14],
  "NeutrinoMasses_eV" -> massesEV,
  "MassSquaredDifferences_eV2" -> {massesEV[[2]]^2 - massesEV[[1]]^2,
    massesEV[[3]]^2 - massesEV[[1]]^2},
  "C5RealCSV" -> c5RealPath, "C5ImagCSV" -> c5ImagPath,
  "NeutrinoMassRealCSV" -> mnuRealPath,
  "NeutrinoMassImagCSV" -> mnuImagPath,
  "CouplingInventoryPath" -> inventoryPath,
  "BenchmarkPurpose" -> "Illustrative unrenormalized one-loop finite-part diagnostic; not a fitted or scheme-certified neutrino spectrum"
|>;
Export[reportPath, D2SafeJSON[report], "RawJSON"];
Print["[C5 benchmark] Status: ", report["Status"]];
Print["[C5 benchmark] Neutrino masses (Takagi singular values, eV): ", massesEV];
Print["[C5 benchmark] C5 relative symmetry residual: ", symResidual/symScale];
Exit[0];
