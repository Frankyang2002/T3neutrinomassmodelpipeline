(* TestWeinberg.wl
   Consolidated direct normalization regression for the matched Weinberg operator.
   Projects the actual Matchete LLHH terms onto nu nu H0 H0. *)

ClearAll[
  InternalHeadName,
  InternalSymbolName,
  ParseWeinbergRawTerms,
  SU2DummyName,
  SU2FundIndexQ,
  NeutralIndexRules,
  EpsilonValue,
  NeutralProjectTerm,
  SafeRatio
];

InternalHeadName[x_] := Quiet @ Check[
  SymbolName[Head[x]],
  ToString[Head[x], InputForm]
];

InternalSymbolName[x_Symbol] := SymbolName[x];
InternalSymbolName[x_] := ToString[x, InputForm];

ParseWeinbergRawTerms[path_String] := Module[
  {text, chunks},

  If[!FileExistsQ[path],
    Print["ERROR: missing file: ", path];
    Return[$Failed]
  ];

  text = Import[path, "Text"];

  chunks = Select[
    StringTrim /@ StringSplit[
      text,
      RegularExpression["--- Weinberg term [0-9]+ ---\\s*"]
    ],
    StringLength[#] > 0 &
  ];

  Quiet @ Check[
    ToExpression[#, InputForm] & /@ chunks,
    $Failed
  ]
];

SU2DummyName[index_] := Module[{args},
  If[InternalHeadName[index] =!= "Index", Return[Missing["NotIndex"]]];
  args = List @@ index;
  If[Length[args] < 2, Return[Missing["BadIndex"]]];
  InternalSymbolName[args[[1]]]
];

SU2FundIndexQ[index_] := Module[{args},
  If[InternalHeadName[index] =!= "Index", Return[False]];
  args = List @@ index;
  If[Length[args] < 2, Return[False]];
  ToString[args[[2]], InputForm] === "SU2L[fund]"
];

NeutralIndexRules[term_] := Module[
  {hIndices, lIndices, rules},

  (* In Matchete's SU(2) fundamental ordering:
       H = (H+, H0), so H0 is component 2.
       L = (nu, e),  so nu is component 1. *)
  hIndices = Cases[
    Unevaluated[term],
    object_ /;
      InternalHeadName[Unevaluated[object]] === "Field" &&
      Quiet @ Check[
        InternalSymbolName[(List @@ Unevaluated[object])[[1]]] === "H",
        False
      ] :>
      Quiet @ Check[
        Cases[
          (List @@ Unevaluated[object])[[3]],
          idx_ /; TrueQ[SU2FundIndexQ[idx]] :>
            SU2DummyName[idx],
          Infinity
        ],
        {}
      ],
    Infinity
  ];

  lIndices = Cases[
    Unevaluated[term],
    object_ /;
      InternalHeadName[Unevaluated[object]] === "Field" &&
      Quiet @ Check[
        InternalSymbolName[(List @@ Unevaluated[object])[[1]]] === "l",
        False
      ] :>
      Quiet @ Check[
        Cases[
          (List @@ Unevaluated[object])[[3]],
          idx_ /; TrueQ[SU2FundIndexQ[idx]] :>
            SU2DummyName[idx],
          Infinity
        ],
        {}
      ],
    Infinity
  ];

  hIndices = DeleteDuplicates @ Cases[Flatten[hIndices], _String];
  lIndices = DeleteDuplicates @ Cases[Flatten[lIndices], _String];

  rules = Join[
    Thread[hIndices -> 2],
    Thread[lIndices -> 1]
  ];

  Association[rules]
];

EpsilonValue[a_Integer, b_Integer] := Which[
  a == 1 && b == 2, 1,
  a == 2 && b == 1, -1,
  True, 0
];

EpsilonValue[___] := Missing["UnresolvedEpsilon"];

NeutralProjectTerm[term_] := Module[
  {
    indexMap,
    hFields,
    spinorChains,
    cgFactors,
    epsilonPairs,
    epsilonValues,
    denominator,
    projected,
    unresolved
  },

  indexMap = NeutralIndexRules[term];

  If[
    MemberQ[Values[indexMap], _Missing] ||
    !SubsetQ[{1, 2}, DeleteDuplicates[Values[indexMap]]],
    Return[
      <|
        "Projected" -> term,
        "UnresolvedCG" -> {},
        "IndexMap" -> indexMap,
        "EpsilonPairs" -> {},
        "ProjectionFailure" -> "Neutral SU(2) index map is incomplete"
      |>
    ]
  ];

  (* Collect the external structures directly.  The Matchete expressions in
     c5_raw.txt are multiplicative LLHH terms, so dividing out the exact
     factors is more robust than relying on ReplaceAll to descend through
     Matchete's context-qualified/held heads. *)
  hFields = Cases[
    term,
    object_ /;
      InternalHeadName[object] === "Field" &&
      Quiet @ Check[
        InternalSymbolName[(List @@ object)[[1]]] === "H",
        False
      ],
    Infinity
  ];

  spinorChains = Cases[
    term,
    object_ /; InternalHeadName[object] === "NCM",
    Infinity
  ];

  cgFactors = Cases[
    term,
    object_ /; InternalHeadName[object] === "CG",
    Infinity
  ];

  epsilonPairs = (
    Module[{args = List @@ #, indices, names},
      If[Length[args] =!= 2, Return[{}]];
      indices = args[[2]];
      names = Cases[
        indices,
        idx_ /; TrueQ[SU2FundIndexQ[idx]] :> SU2DummyName[idx],
        Infinity
      ];
      Lookup[indexMap, #, Missing["UnknownIndex"]] & /@ names
    ] & /@ cgFactors
  );

  If[
    !AllTrue[epsilonPairs, MatchQ[#, {_Integer, _Integer}] &],
    Return[
      <|
        "Projected" -> term,
        "UnresolvedCG" -> cgFactors,
        "IndexMap" -> indexMap,
        "EpsilonPairs" -> epsilonPairs,
        "ProjectionFailure" -> "Could not resolve SU(2) epsilon component pairs"
      |>
    ]
  ];

  epsilonValues = (EpsilonValue @@ #) & /@ epsilonPairs;

  If[MemberQ[epsilonValues, _Missing],
    Return[
      <|
        "Projected" -> term,
        "UnresolvedCG" -> cgFactors,
        "IndexMap" -> indexMap,
        "EpsilonPairs" -> epsilonPairs,
        "ProjectionFailure" -> "Could not evaluate SU(2) epsilon components"
      |>
    ]
  ];

  denominator = Times @@ Join[hFields, spinorChains, cgFactors];

  projected = Quiet @ Check[
    FactorTerms @ Cancel @ Together[
      term * Times @@ epsilonValues / denominator
    ],
    Simplify[term * Times @@ epsilonValues / denominator]
  ];

  unresolved = Cases[
    projected,
    object_ /; MemberQ[{"Field", "NCM", "CG"}, InternalHeadName[object]],
    Infinity
  ];

  <|
    "Projected" -> projected,
    "UnresolvedCG" -> Select[
      unresolved,
      InternalHeadName[#] === "CG" &
    ],
    "UnresolvedExternal" -> Select[
      unresolved,
      MemberQ[{"Field", "NCM"}, InternalHeadName[#]] &
    ],
    "IndexMap" -> indexMap,
    "EpsilonPairs" -> epsilonPairs,
    "EpsilonValues" -> epsilonValues
  |>
];

SafeRatio[a_, b_] := Quiet @ Check[
  FactorTerms @ Cancel @ Together[a/b],
  Simplify[a/b]
];


(* ------------------------------------------------------------------------- *)
(* Locate output files                                                       *)
(* ------------------------------------------------------------------------- *)

scriptDirectory = DirectoryName @ ExpandFileName[$InputFileName];
projectRoot = ExpandFileName @ FileNameJoin[{scriptDirectory, "..", ".."}];

ResolveWeinbergPaths[] := Module[
  {
    requested,
    requestedParent,
    outputRoot,
    canonicalRun,
    rawCandidates,
    validRuns,
    ranked,
    runDirectory,
    rawPath,
    coefficientPath
  },

  If[Length[$ScriptCommandLine] >= 2,
    requested = ExpandFileName[$ScriptCommandLine[[2]]];

    (* Current pipeline layout:
         <run>/c5_raw.txt
         <run>/data/c5_coefficient.txt *)
    If[
      FileExistsQ[FileNameJoin[{requested, "c5_raw.txt"}]] &&
      FileExistsQ[FileNameJoin[{requested, "data", "c5_coefficient.txt"}]],
      Return[
        <|
          "OutputDirectory" -> requested,
          "RawPath" -> FileNameJoin[{requested, "c5_raw.txt"}],
          "CoefficientPath" ->
            FileNameJoin[{requested, "data", "c5_coefficient.txt"}]
        |>
      ]
    ];

    (* Also accept the data directory itself as the explicit argument. *)
    requestedParent = DirectoryName[requested];
    If[
      FileNameTake[requested] === "data" &&
      FileExistsQ[FileNameJoin[{requestedParent, "c5_raw.txt"}]] &&
      FileExistsQ[FileNameJoin[{requested, "c5_coefficient.txt"}]],
      Return[
        <|
          "OutputDirectory" -> requestedParent,
          "RawPath" -> FileNameJoin[{requestedParent, "c5_raw.txt"}],
          "CoefficientPath" -> FileNameJoin[{requested, "c5_coefficient.txt"}]
        |>
      ]
    ];

    (* Retain compatibility with older layouts that colocated both files. *)
    If[
      FileExistsQ[FileNameJoin[{requested, "c5_raw.txt"}]] &&
      FileExistsQ[FileNameJoin[{requested, "c5_coefficient.txt"}]],
      Return[
        <|
          "OutputDirectory" -> requested,
          "RawPath" -> FileNameJoin[{requested, "c5_raw.txt"}],
          "CoefficientPath" -> FileNameJoin[{requested, "c5_coefficient.txt"}]
        |>
      ]
    ];

    Print[
      "ERROR: explicit Weinberg output path does not contain the required ",
      "c5_raw.txt / c5_coefficient.txt artifacts: ",
      requested
    ];
    Return[$Failed];
  ];

  outputRoot = FileNameJoin[{projectRoot, "output"}];
  If[!DirectoryQ[outputRoot],
    Print["ERROR: output directory does not exist: ", outputRoot];
    Return[$Failed];
  ];

  (* Prefer the canonical direct pipeline output when it exists. *)
  canonicalRun = FileNameJoin[{outputRoot, "T3_B_alpha_m1"}];
  If[
    FileExistsQ[FileNameJoin[{canonicalRun, "c5_raw.txt"}]] &&
    FileExistsQ[
      FileNameJoin[{canonicalRun, "data", "c5_coefficient.txt"}]
    ],
    Return[
      <|
        "OutputDirectory" -> canonicalRun,
        "RawPath" -> FileNameJoin[{canonicalRun, "c5_raw.txt"}],
        "CoefficientPath" ->
          FileNameJoin[{canonicalRun, "data", "c5_coefficient.txt"}]
      |>
    ]
  ];

  (* Fallback for batch/smoke/full outputs. *)
  rawCandidates = FileNames["c5_raw.txt", outputRoot, Infinity];
  validRuns = Select[
    DirectoryName /@ rawCandidates,
    StringContainsQ[#, "T3_B_alpha_m1"] &&
    FileExistsQ[FileNameJoin[{#, "data", "c5_coefficient.txt"}]] &
  ];

  If[validRuns === {},
    Print[
      "ERROR: no T3_B_alpha_m1 output with <run>/c5_raw.txt and ",
      "<run>/data/c5_coefficient.txt was found below ",
      outputRoot
    ];
    Return[$Failed];
  ];

  ranked = Reverse @ SortBy[
    DeleteDuplicates[validRuns],
    FileDate[FileNameJoin[{#, "c5_raw.txt"}]] &
  ];
  runDirectory = First[ranked];
  rawPath = FileNameJoin[{runDirectory, "c5_raw.txt"}];
  coefficientPath =
    FileNameJoin[{runDirectory, "data", "c5_coefficient.txt"}];

  <|
    "OutputDirectory" -> runDirectory,
    "RawPath" -> rawPath,
    "CoefficientPath" -> coefficientPath
  |>
];

weinbergPaths = ResolveWeinbergPaths[];
If[weinbergPaths === $Failed, Exit[1]];

outputDirectory = weinbergPaths["OutputDirectory"];
rawPath = weinbergPaths["RawPath"];
coefficientPath = weinbergPaths["CoefficientPath"];

Print["========================================================================"];
Print["WEINBERG -> NEUTRINO MASS NORMALIZATION TEST"];
Print["========================================================================"];
Print[];
Print["Output directory: ", outputDirectory];
Print[];


(* ------------------------------------------------------------------------- *)
(* Read the actual Matchete output                                           *)
(* ------------------------------------------------------------------------- *)

terms = ParseWeinbergRawTerms[rawPath];

If[terms === $Failed,
  Print["FAIL: could not parse c5_raw.txt"];
  Exit[1]
];

If[!FileExistsQ[coefficientPath],
  Print["FAIL: missing c5_coefficient.txt"];
  Exit[1]
];

c5 = Quiet @ Check[
  ToExpression[Import[coefficientPath, "Text"], InputForm],
  $Failed
];

If[c5 === $Failed,
  Print["FAIL: could not parse c5_coefficient.txt"];
  Exit[1]
];

(* The production extractor identifies the P_L sector as holomorphic.
   Reproduce that selection using the textual Proj[-1] marker. *)
holomorphicTerms = Select[
  terms,
  StringContainsQ[
    ToString[Unevaluated[#], InputForm],
    "Proj[-1]"
  ] &
];

Print["Total raw Weinberg terms: ", Length[terms]];
Print["Holomorphic terms: ", Length[holomorphicTerms]];
Print[];

If[holomorphicTerms === {},
  Print["FAIL: no holomorphic terms found"];
  Exit[1]
];


(* ------------------------------------------------------------------------- *)
(* Project LLHH onto nu nu H0 H0                                             *)
(* ------------------------------------------------------------------------- *)

projectedData = NeutralProjectTerm /@ holomorphicTerms;
projectedTerms = Lookup[projectedData, "Projected"];
unresolved = Flatten @ Lookup[projectedData, "UnresolvedCG"];
projectionFailures = DeleteMissing @ Lookup[
  projectedData,
  "ProjectionFailure",
  Missing["NotAvailable"]
];

If[projectionFailures =!= {},
  Print["FAIL: neutral SU(2) index projection failed."];
  Print[projectionFailures];
  Print["Index maps used:"];
  Print[Lookup[projectedData, "IndexMap"]];
  Print["Epsilon component pairs:"];
  Print[Lookup[projectedData, "EpsilonPairs"]];
  Exit[1]
];

unresolvedExternal = Flatten @ Lookup[
  projectedData,
  "UnresolvedExternal",
  {}
];

If[unresolvedExternal =!= {},
  Print["FAIL: neutral projection left external Field/NCM structure unresolved."];
  Print[InputForm /@ unresolvedExternal];
  Print["Projected terms:"];
  Print[InputForm /@ projectedTerms];
  Print["Index maps used:"];
  Print[Lookup[projectedData, "IndexMap"]];
  Exit[1]
];

If[unresolved =!= {},
  Print["FAIL: some SU(2) CG tensors did not resolve."];
  Print["Unresolved CG tensors:"];
  Print[InputForm /@ unresolved];
  Print[];
  Print["Index maps used:"];
  Print[Lookup[projectedData, "IndexMap"]];
  Print["Epsilon component pairs:"];
  Print[Lookup[projectedData, "EpsilonPairs"]];
  Print["Epsilon values:"];
  Print[Lookup[projectedData, "EpsilonValues", {}]];
  Exit[1]
];

neutralCoefficient = Quiet @ Check[
  FactorTerms @ Cancel @ Together[Total[projectedTerms]],
  Simplify[Total[projectedTerms]]
];

ratio = SafeRatio[neutralCoefficient, c5];

Print["Exported C5:"];
Print[InputForm[c5]];
Print[];

Print["Projected coefficient of nu nu H0 H0:"];
Print[InputForm[neutralCoefficient]];
Print[];

Print["r = C_(nu nu H0 H0) / C5:"];
Print[InputForm[ratio]];
Print[];

Print["Therefore, with H0 -> v/Sqrt[2],"];
Print["  L_EFT -> (r C5 v^2 / 2) nu nu + h.c."];
Print["and comparison with"];
Print["  L_mass = -(1/2) m_nu nu nu + h.c."];
Print["gives"];
Print["  m_nu = -r v^2 C5."];
Print[];

Which[
  TrueQ[Simplify[ratio - 1] === 0],
    Print["RESULT: r = 1"];
    Print["Use m_nu = -v^2 C5."],
  TrueQ[Simplify[ratio - 1/2] === 0],
    Print["RESULT: r = 1/2"];
    Print["Use m_nu = -(v^2/2) C5."],
  True,
    Print["FAIL: non-standard Weinberg normalization ratio: ", InputForm[ratio]];
    Exit[1]
];

Print[];
Print["PASS: direct Weinberg component projection and normalization."];
Exit[0];

