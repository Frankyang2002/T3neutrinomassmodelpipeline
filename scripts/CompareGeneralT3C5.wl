(* Experimental *formal* C5 regression for Matchete baseline/extended runs.
   Requires no FrontEnd and no new Lagrangian construction.
   Reads C5 InputForm exactly as exported by Lagrangian/RunModel.wl.

   wolframscript -file scripts/CompareGeneralT3C5.wl \
      output/general_t3_trial output/general_t3_trial/c5_algebraic_report.json

   There is NO substitution of UV boundary values in this script. In particular,
   zeroing extra quartics is a formal algebraic test, not an RGE evolution.
*)

ClearAll[GeneralC5Parse, GeneralC5ExtraCouplingQ, GeneralC5ZeroExtra,
  GeneralC5Describe, GeneralC5OneModel, GeneralC5Report];

GeneralC5Parse[path_String] := Module[{txt, held},
  txt = Quiet@Check[Import[path, "Text"], $Failed];
  If[!StringQ[txt] || StringTrim[txt] === "", Return[$Failed]];
  held = Quiet@Check[ToExpression[txt, InputForm, HoldComplete], $Failed];
  If[held === $Failed, Return[$Failed]];
  Quiet@Check[ReleaseHold[held], $Failed]
];

GeneralC5ExtraCouplingQ[c_Coupling] :=
  StringStartsQ[Last[StringSplit[ToString[First[c], InputForm], "`"]], "lambdaT3Extra"];
GeneralC5ExtraCouplingQ[_] := False;

(* Match Bar[Coupling] BEFORE the interior Coupling to avoid Bar[0]. *)
GeneralC5ZeroExtra[expr_] := expr /. {
  HoldPattern[Bar[c_Coupling]] /; GeneralC5ExtraCouplingQ[c] :> 0,
  HoldPattern[Conjugate[c_Coupling]] /; GeneralC5ExtraCouplingQ[c] :> 0,
  c_Coupling /; GeneralC5ExtraCouplingQ[c] :> 0
};

(* Identify only INDEPENDENT coupling atoms. Cases[..., Coupling[...], Infinity]
   also descends inside Bar[Coupling[...]], making the bare coupling appear
   spuriously even when only its conjugate is physically present. *)
GeneralC5Describe[expr_] := Module[{conjugateAtoms, withoutConjugates, bareAtoms},
  conjugateAtoms = DeleteDuplicates@Join[
    Cases[expr, HoldPattern[Bar[c_Coupling]] /;
      GeneralC5ExtraCouplingQ[c] :> Bar[c], Infinity],
    Cases[expr, HoldPattern[Conjugate[c_Coupling]] /;
      GeneralC5ExtraCouplingQ[c] :> Conjugate[c], Infinity]
  ];
  withoutConjugates = expr /. {
    HoldPattern[Bar[c_Coupling]] /; GeneralC5ExtraCouplingQ[c] :>
      Unique["t3ConjugateMasked"],
    HoldPattern[Conjugate[c_Coupling]] /; GeneralC5ExtraCouplingQ[c] :>
      Unique["t3ConjugateMasked"]
  };
  bareAtoms = Cases[withoutConjugates,
    c_Coupling /; GeneralC5ExtraCouplingQ[c] :> c, Infinity];
  DeleteDuplicates[Join[conjugateAtoms, bareAtoms]]
];

GeneralC5ProvenZero[expr_] := Module[{reduced = Quiet@Check[Expand[expr], $Failed]},
  If[reduced === $Failed, Return[False]];
  If[TrueQ[reduced === 0] || TrueQ[PossibleZeroQ[reduced]],
    Return[True]];
  TrueQ[TimeConstrained[
    Quiet@Check[FullSimplify[reduced == 0], False], 30, False]]
];

GeneralC5OneModel[folder_String] := Module[
  {name, baselinePath, extendedPath, originalText, newText, result,
   oldExpr, newExpr, zeroExpr, zeroDifference, fullDifference, reduced,
   status, names, atoms, probes, probed, polynomialQ, independentQ,
   contributions, outputFolder, physicalDifference, coefficients, residual,
   responseZero, residualProven},

  name = FileNameTake[folder];
  baselinePath = FileNameJoin[{folder, "baseline", "c5_coefficient.txt"}];
  extendedPath = FileNameJoin[{folder, "extended", "c5_coefficient.txt"}];
  result = <|"model" -> name, "status" -> "NotRun", "zero_extra_limit_proven" -> False,
    "extra_c5_couplings" -> {}, "extra_c5_atom_count" -> 0,
    "linear_in_extra_coupling_atoms" -> Null,
    "coefficient_independent_of_other_extra_couplings" -> Null|>;
  If[!FileExistsQ[baselinePath] || !FileExistsQ[extendedPath],
    Return[Association[result, "status" -> "MissingC5Artifact"]]
  ];
  originalText = Import[baselinePath, "Text"];
  newText = Import[extendedPath, "Text"];
  names = Sort@DeleteDuplicates@StringCases[newText,
    RegularExpression["lambdaT3Extra[A-Za-z0-9]+"]];
  result = Association[result, "extra_c5_couplings" -> names];

  (* Identity as exported is a stronger and much cheaper check for unchanged models. *)
  If[StringReplace[originalText, WhitespaceCharacter .. -> ""] ===
     StringReplace[newText, WhitespaceCharacter .. -> ""],
    Return[Association[result, "status" -> "TextIdentical",
      "zero_extra_limit_proven" -> True]]
  ];

  oldExpr = TimeConstrained[GeneralC5Parse[baselinePath], 90, $Aborted];
  newExpr = TimeConstrained[GeneralC5Parse[extendedPath], 90, $Aborted];
  If[oldExpr === $Aborted || newExpr === $Aborted,
    Return[Association[result, "status" -> "ParseTimedOut"]]
  ];
  If[oldExpr === $Failed || newExpr === $Failed,
    Return[Association[result, "status" -> "SymbolicParseFailed"]]
  ];

  (* The comparison here uses the exported FORMAL coefficients. In particular
     it does not assume a numerical value or a normalisation for RGBeta couplings. *)
  zeroExpr = GeneralC5ZeroExtra[newExpr];
  zeroDifference = TimeConstrained[Quiet@Check[
    Expand[zeroExpr - oldExpr], $Failed], 120, $Aborted];
  If[zeroDifference === $Aborted,
    Return[Association[result, "status" -> "ReductionTimedOut"]]
  ];
  If[zeroDifference === $Failed,
    Return[Association[result, "status" -> "ReductionFailed"]]
  ];
  reduced = TrueQ[zeroDifference === 0] || TrueQ[PossibleZeroQ[zeroDifference]];
  If[!reduced,
    reduced = TrueQ[TimeConstrained[Quiet@Check[
      FullSimplify[zeroDifference == 0], False], 45, False]]
  ];
  result = Association[result, "zero_extra_limit_proven" -> reduced];

  (* Retain an exact symbolic delta for inspection. The file is never fed
     automatically to production numerical analysis. *)
  outputFolder = FileNameJoin[{folder, "general_c5_algebraic"}];
  If[!DirectoryQ[outputFolder], CreateDirectory[outputFolder,
    CreateIntermediateDirectories -> True]];
  fullDifference = TimeConstrained[Quiet@Check[
    Expand[newExpr - oldExpr], $Failed], 120, $Aborted];
  If[fullDifference === $Aborted || fullDifference === $Failed,
    Return[Association[result, "status" -> If[reduced,
      "ZeroExtraProvenButDeltaExportFailed", "DeltaExportFailed"]]]
  ];
  Export[FileNameJoin[{outputFolder, "delta_c5_inputform.txt"}],
    ToString[fullDifference, InputForm], "Text"];
  Export[FileNameJoin[{outputFolder, "zero_extra_residual_inputform.txt"}],
    ToString[zeroDifference, InputForm], "Text"];

  (* Compare to the same extended expression with EXTRA couplings switched
     off. This is exactly the coupling response; it cannot be polluted by a
     baseline subtraction that is equal only after nontrivial simplification.
     The independent baseline equality was proved above and is preserved. *)
  physicalDifference = TimeConstrained[Quiet@Check[
    Expand[newExpr - GeneralC5ZeroExtra[newExpr]], $Failed], 120, $Aborted];
  If[physicalDifference === $Aborted || physicalDifference === $Failed,
    Return[Association[result, "status" -> "ResponseExtractionFailed"]]
  ];
  Export[FileNameJoin[{outputFolder, "delta_c5_response_inputform.txt"}],
    ToString[physicalDifference, InputForm], "Text"];

  atoms = GeneralC5Describe[physicalDifference];
  probes = Table[Unique["c5Probe"], {Length[atoms]}];
  probed = Expand[physicalDifference /. Thread[atoms -> probes]];
  polynomialQ = TrueQ[PolynomialQ[probed, probes]];
  coefficients = Coefficient[probed, #] & /@ probes;
  responseZero = If[Length[probes] > 0,
    Expand[probed /. Thread[probes -> ConstantArray[0, Length[probes]]]],
    probed];
  residual = Expand[probed - Total[MapThread[Times,
    {probes, coefficients}]]];
  residualProven = GeneralC5ProvenZero[residual];

  contributions = Table[
    Module[{part = coefficients[[j]]},
      <|"coupling_atom" -> ToString[atoms[[j]], InputForm],
        "linear_coefficient_inputform" -> ToString[part, InputForm],
        "coefficient_formally_zero" -> GeneralC5ProvenZero[part],
        "depends_on_other_extra_atoms" -> !FreeQ[part,
          Alternatives @@ probes]|>
    ],
    {j, Length[atoms]}
  ];
  (* Zero derivative artifacts are not physical couplings. Preserve their
     diagnostic number but do not report them as active contributions. *)
  contributions = Select[contributions,
    !TrueQ[Lookup[#, "coefficient_formally_zero", False]] &];
  independentQ = polynomialQ && residualProven &&
    GeneralC5ProvenZero[responseZero] &&
    Length[contributions] > 0 &&
    AllTrue[contributions,
      !TrueQ[Lookup[#, "depends_on_other_extra_atoms", True]] &];

  Export[FileNameJoin[{outputFolder, "delta_c5_atoms.json"}],
    <|"model" -> name, "atoms" -> contributions,
      "response_reconstruction_residual_inputform" ->
        ToString[residual, InputForm],
      "response_reconstruction_residual_proven_zero" -> residualProven,
      "response_at_zero_inputform" -> ToString[responseZero, InputForm],
      "response_at_zero_proven_zero" -> GeneralC5ProvenZero[responseZero],
      "polynomial_in_extra_atoms" -> polynomialQ|>, "RawJSON"];

  Association[result, "status" -> If[reduced,
    "ZeroExtraLimitProven", "ZeroExtraLimitNotProven"],
    "extra_c5_atom_count" -> Length[contributions],
    "linear_in_extra_coupling_atoms" -> polynomialQ && independentQ,
    "coefficient_independent_of_other_extra_couplings" -> independentQ,
    "response_reconstruction_residual_proven_zero" -> residualProven,
    "response_at_zero_proven_zero" -> GeneralC5ProvenZero[responseZero],
    "delta_inputform_file" -> "general_c5_algebraic/delta_c5_inputform.txt",
    "atom_coefficients_file" -> "general_c5_algebraic/delta_c5_atoms.json"
  ]
];

GeneralC5Report[root_String, target_String] := Module[{folders, results, report},
  folders = Select[FileNames["*", root], DirectoryQ];
  results = GeneralC5OneModel /@ folders;
  report = <|
    "status" -> If[Length[results] === 16 &&
      AllTrue[results, MemberQ[{ "TextIdentical", "ZeroExtraLimitProven" },
        # ["status"]] &], "AllZeroExtraLimitsProven", "ReviewRequired"],
    "models_checked" -> Length[results],
    "formal_only" -> True,
    "remark" -> "This verifies reduction of formal Matchete C5 at zero new quartics; it does not validate RGBeta normalisation, numerical RGE closure or quadratic mixing.",
    "results" -> results
  |>;
  If[!DirectoryQ[DirectoryName[target]], CreateDirectory[DirectoryName[target],
    CreateIntermediateDirectories -> True]];
  Export[target, report, "RawJSON"];
  Scan[Function[row, Print[row["model"], "  ", row["status"],
    "  new=", Length[row["extra_c5_couplings"]],
    "  zeroLimit=", row["zero_extra_limit_proven"]]], results];
  Print["Report: ", target];
  If[report["status"] === "AllZeroExtraLimitsProven", 0, 2]
];

(* Standalone kernel self-test: this uses no Matchete installation or trial files. *)
GeneralC5SelfTest[] := Module[
  {a, b, linearConjugate, linearMixed, quadratic, atoms, probes,
   substituted, coefficients, residual, cases, ok},
  a = Coupling[lambdaT3ExtraQ020020I1, {}, 0];
  b = Coupling[lambdaT3ExtraQ020200I1, {}, 0];
  linearConjugate = q1 Bar[a] + q2 Bar[b];
  linearMixed = q1 Bar[a] + q2 b;
  quadratic = q1 Bar[a] b + q2 a^2;
  cases = {
    <|"name" -> "two_conjugate_atoms", "expression" -> linearConjugate,
      "expected_atom_count" -> 2, "expect_linear" -> True|>,
    <|"name" -> "conjugate_and_direct", "expression" -> linearMixed,
      "expected_atom_count" -> 2, "expect_linear" -> True|>,
    <|"name" -> "nonlinear_rejected", "expression" -> quadratic,
      "expected_atom_count" -> 3, "expect_linear" -> False|>
  };
  ok = True;
  Scan[Function[test,
    atoms = GeneralC5Describe[test["expression"]];
    probes = Table[Unique["selfTestProbe"], {Length[atoms]}];
    substituted = Expand[(test["expression"] -
      GeneralC5ZeroExtra[test["expression"]]) /.
      Thread[atoms -> probes]];
    coefficients = Coefficient[substituted, #] & /@ probes;
    residual = Expand[substituted -
      Total[MapThread[Times, {probes, coefficients}]]];
    If[Length[atoms] =!= test["expected_atom_count"] ||
      GeneralC5ProvenZero[residual] =!= test["expect_linear"],
      ok = False;
      Print["FAILED: ", test["name"], " atoms=", Length[atoms],
        " residual=", ToString[residual, InputForm]],
      Print["PASS: ", test["name"]]
    ];
  ], cases];
  If[ok, Print["General C5 atom self-test: PASS"]; 0,
    Print["General C5 atom self-test: FAIL"]; 2]
];

(* The wolframscript executable may consume arguments beginning with '--'
   as its own options before populating $ScriptCommandLine.  In -file mode
   use a plain positional SELFTEST token, never --selftest. *)
args = Rest[$ScriptCommandLine];
If[args === {"SELFTEST"} || args === {"--selftest"},
  Exit[GeneralC5SelfTest[]]
];
If[Length[args] < 1 || Length[args] > 2,
  Print["Usage: wolframscript -file scripts/CompareGeneralT3C5.wl ROOT [REPORT_JSON]"];
  Print["For the self-test: wolframscript -file scripts/CompareGeneralT3C5.wl SELFTEST"];
  Print["Received script arguments: ", ToString[args, InputForm]];
  Exit[64]
];
root = ExpandFileName[args[[1]]];
If[!DirectoryQ[root], Print["Missing trial output directory: ", root]; Exit[66]];
target = If[Length[args] >= 2, ExpandFileName[args[[2]]],
  FileNameJoin[{root, "c5_algebraic_report.json"}]];
Exit[GeneralC5Report[root, target]];
