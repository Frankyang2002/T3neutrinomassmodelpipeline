(* Direct normalization test for the matched Matchete Weinberg operator.
   This does not alter any production files. *)

ClearAll[
  InternalHeadName,
  InternalSymbolName,
  ParseWeinbergRawTerms,
  SU2DummyName,
  NeutralIndexRules,
  EpsilonValue,
  NeutralProjectTerm,
  SafeRatio
];

InternalHeadName[x_] := Quiet @ Check[
  SymbolName[Unevaluated[Head[x]]],
  ToString[Unevaluated[Head[x]], InputForm]
];

InternalSymbolName[x_Symbol] := SymbolName[Unevaluated[x]];
InternalSymbolName[x_] := ToString[Unevaluated[x], InputForm];

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
  args = List @@ Unevaluated[index];
  If[Length[args] < 2, Return[Missing["BadIndex"]]];
  InternalSymbolName[args[[1]]]
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
          idx_ /; InternalHeadName[Unevaluated[idx]] === "Index" :>
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
          idx_ /; InternalHeadName[Unevaluated[idx]] === "Index" :>
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
  {indexMap, projected, epsRules, unresolved},

  indexMap = NeutralIndexRules[term];

  projected = Unevaluated[term] /. {
    (* Strip the fermion spinor chain only after its SU(2) index has been read. *)
    object_ /; InternalHeadName[Unevaluated[object]] === "NCM" :> 1,

    (* Remove the two external Higgs fields after recording H0 = component 2. *)
    object_ /;
      InternalHeadName[Unevaluated[object]] === "Field" &&
      Quiet @ Check[
        InternalSymbolName[(List @@ Unevaluated[object])[[1]]] === "H",
        False
      ] :> 1
  };

  (* Evaluate epsilon tensors after assigning every external neutral component. *)
  epsRules = {
    object_ /; InternalHeadName[Unevaluated[object]] === "CG" :>
      Module[{args, tensor, indices, values},
        args = List @@ Unevaluated[object];
        If[Length[args] =!= 2, Return[object]];

        tensor = args[[1]];
        indices = args[[2]];

        If[
          FreeQ[
            ToString[Unevaluated[tensor], InputForm],
            "eps"
          ],
          Return[object]
        ];

        values = indices /. {
          HoldPattern[Bar[idx_]] :> idx,
          idx_ /; InternalHeadName[Unevaluated[idx]] === "Index" :>
            Lookup[indexMap, SU2DummyName[idx], Missing["UnknownIndex"]]
        };

        If[
          MatchQ[values, {_Integer, _Integer}],
          EpsilonValue @@ values,
          object
        ]
      ]
  };

  projected = projected /. epsRules;

  unresolved = Cases[
    projected,
    object_ /; InternalHeadName[Unevaluated[object]] === "CG",
    Infinity
  ];

  <|
    "Projected" -> Quiet @ Check[Simplify[Expand[projected]], projected],
    "UnresolvedCG" -> unresolved,
    "IndexMap" -> indexMap
  |>
];

SafeRatio[a_, b_] := Quiet @ Check[
  FactorTerms @ Cancel @ Together[a/b],
  Simplify[a/b]
];


(* ------------------------------------------------------------------------- *)
(* Locate output files                                                       *)
(* ------------------------------------------------------------------------- *)

projectRoot = Directory[];

outputDirectory = If[
  Length[$ScriptCommandLine] >= 2,
  $ScriptCommandLine[[2]],
  FileNameJoin[
    {
      projectRoot,
      "wolfram",
      "output",
      "T3_B_alpha_m1"
    }
  ]
];

rawPath = FileNameJoin[{outputDirectory, "c5_raw.txt"}];
coefficientPath = FileNameJoin[{outputDirectory, "c5_coefficient.txt"}];

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

If[unresolved =!= {},
  Print["FAIL: some SU(2) CG tensors did not resolve."];
  Print["Unresolved CG tensors:"];
  Print[InputForm /@ unresolved];
  Print[];
  Print["Index maps used:"];
  Print[Lookup[projectedData, "IndexMap"]];
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
    Print["RESULT: non-standard ratio. Do not change the production normalization yet."]
];

Print[];
Print["PASS: direct component projection completed."];

