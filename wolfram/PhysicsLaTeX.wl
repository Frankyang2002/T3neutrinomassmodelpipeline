(* ::Package:: *)

(*
  PhysicsLaTeX.wl

  Recursive display formatter for Matchete expressions.

  The formatter:
    - preserves distinct Lorentz, flavour and gauge indices;
    - formats scalar and fermion fields;
    - formats conjugation, transpose, projectors and gamma matrices;
    - formats indexed couplings and SU(2) epsilon tensors;
    - formats complete noncommutative fermion chains;
    - reports any Matchete objects that remain unconverted.
*)

ClearAll[
  PhysicsSymbolName,
  PhysicsIndexBase,
  PhysicsIndexSerial,
  PhysicsIndex,
  PhysicsIndices,
  PhysicsFieldBase,
  PhysicsField,
  PhysicsBar,
  PhysicsTranspose,
  PhysicsCouplingBase,
  PhysicsCoupling,
  PhysicsCG,
  PhysicsDiracFactor,
  PhysicsFermionChain,
  PhysicsFieldStrength,
  PhysicsDisplayForm,
  RemainingInternalObjects,
  FullyConvertedQ,
  ExpressionToLaTeX
];

(* ---------------------------------------------------------------------- *)
(* Helpers that are deliberately based on symbol names, so the formatter  *)
(* remains useful if Matchete changes the context in which an object lives. *)
(* ---------------------------------------------------------------------- *)

PhysicsSymbolName[s_Symbol] := SymbolName[Unevaluated[s]];
PhysicsSymbolName[x_] := ToString[Unevaluated[x], InputForm];

PhysicsIndexSerial[label_] := Module[{text, hits},
  text = ToString[Unevaluated[label], InputForm];
  hits = StringCases[
    text,
    RegularExpression["d\\$\\$(\\d+)"] -> "$1"
  ];
  If[hits === {}, "", First[hits]]
];

PhysicsIndexBase[type_] := Switch[PhysicsSymbolName[type],
  "Lorentz", \[Mu],
  "Flavor", i,
  "NFlavor", r,
  _,
    Which[
      !FreeQ[Unevaluated[type], SU2L] && !FreeQ[Unevaluated[type], fund], a,
      !FreeQ[Unevaluated[type], SU2L] && !FreeQ[Unevaluated[type], adj], A,
      !FreeQ[Unevaluated[type], SU3c] && !FreeQ[Unevaluated[type], fund], \[Alpha],
      !FreeQ[Unevaluated[type], SU3c] && !FreeQ[Unevaluated[type], adj], X,
      True, j
    ]
];

PhysicsIndex[index_] := Module[
  {headName, arguments, label, type, base, serial},

  headName = PhysicsSymbolName[Head[Unevaluated[index]]];

  (* Matchete uses Bar[Index[...]] for conjugate representation indices.
     For display purposes the dummy-index label is the same; the surrounding
     field/tensor already carries the conjugation information. *)
  If[headName === "Bar",
    arguments = List @@ Unevaluated[index];
    If[
      Length[arguments] === 1 &&
      PhysicsSymbolName[Head[Unevaluated[arguments[[1]]]]] === "Index",
      Return[PhysicsIndex[arguments[[1]]]]
    ];
  ];

  If[headName =!= "Index", Return[index]];

  {label, type} = List @@ Unevaluated[index];
  base = PhysicsIndexBase[type];
  serial = PhysicsIndexSerial[label];

  If[serial === "", base, Subscript[base, serial]]
];

PhysicsIndices[indices_List] := PhysicsIndex /@ indices;
PhysicsIndices[_] := {};

(* ---------------------------------------------------------------------- *)
(* Fields and derivatives.                                                 *)
(* ---------------------------------------------------------------------- *)

PhysicsFieldBase[name_] := Switch[PhysicsSymbolName[name],
  "H", H,
  "NewScalar", \[Eta],
  "q", q,
  "l", \[ScriptL],
  "u", u,
  "d", d,
  "e", e,
  "NewFermion", N,
  _, name
];

PhysicsField[field_] := Module[
  {arguments, name, kind, internalIndices, derivativeIndices, result,
   formattedInternal, formattedDerivatives},

  arguments = List @@ Unevaluated[field];
  If[Length[arguments] < 4, Return[field]];

  {name, kind, internalIndices, derivativeIndices} = Take[arguments, 4];
  result = PhysicsFieldBase[name];

  formattedInternal = PhysicsIndices[internalIndices];
  If[formattedInternal =!= {},
    result = Subscript[result, Row[Riffle[formattedInternal, ","]]]
  ];

  formattedDerivatives = PhysicsIndices[derivativeIndices];
  result = Fold[
    Function[{current, derivativeIndex},
      Subscript[D, derivativeIndex][current]
    ],
    result,
    formattedDerivatives
  ];

  result
];

PhysicsBar[argument_] := Module[{headName, args, kind},
  headName = PhysicsSymbolName[Head[Unevaluated[argument]]];

  If[headName === "Field",
    args = List @@ Unevaluated[argument];
    kind = If[Length[args] >= 2, PhysicsSymbolName[args[[2]]], ""];

    Return[
      If[
        kind === "Scalar",
        Superscript[PhysicsField[argument], \[Dagger]],
        Overscript[PhysicsField[argument], _]
      ]
    ]
  ];

  Overscript[PhysicsDisplayForm[argument], _]
];

PhysicsTranspose[argument_] :=
  Superscript[PhysicsDisplayForm[argument], T];

(* ---------------------------------------------------------------------- *)
(* Couplings and masses, including their real flavour indices.             *)
(* ---------------------------------------------------------------------- *)

PhysicsCouplingBase[name_] := Module[{nameString},
  nameString = PhysicsSymbolName[name];

  Which[
    nameString === "gY", Subscript[g, Y],
    nameString === "gL", Subscript[g, L],
    nameString === "gs", Subscript[g, s],

    nameString === "lambda", \[Lambda],
    nameString === "\[Lambda]", \[Lambda],
    nameString === "lambda3", Subscript[\[Lambda], 3],
    nameString === "lambda4", Subscript[\[Lambda], 4],
    nameString === "lambda5", Subscript[\[Lambda], 5],
    nameString === "lambdaS", Subscript[\[Lambda], S],

    nameString === "yNew", Subscript[y, N],
    nameString === "Yd", Subscript[Y, d],
    nameString === "Yu", Subscript[Y, u],
    nameString === "Ye", Subscript[Y, e],

    nameString === "MEta", Subscript[M, \[Eta]],
    nameString === "MN", Subscript[M, N],
    nameString === "mu2", Superscript[\[Mu], 2],
    nameString === "\[Mu]2", Superscript[\[Mu], 2],

    StringMatchQ[nameString, ("mubar2" | "muBar2") ~~ ___],
      Superscript[OverBar[\[Mu]], 2],

    StringMatchQ[nameString, ("mubar" | "muBar") ~~ ___] &&
      !StringContainsQ[nameString, "2"],
      OverBar[\[Mu]],

    nameString === "CH2", Subscript[C, H2],
    True, name
  ]
];

PhysicsCoupling[coupling_] := Module[
  {arguments, name, nameString, indices, formattedIndices, indexRow},

  arguments = List @@ Unevaluated[coupling];
  If[arguments === {}, Return[coupling]];

  name = arguments[[1]];
  nameString = PhysicsSymbolName[name];

  indices = If[
    Length[arguments] >= 2 && ListQ[arguments[[2]]],
    arguments[[2]],
    {}
  ];

  formattedIndices = PhysicsIndices[indices];
  indexRow = Row[Riffle[formattedIndices, ","]];

  (* Build indexed matrix couplings with exactly one subscript. *)
  Switch[nameString,
    "Yd",
      If[
        formattedIndices === {},
        Subscript[Y, d],
        Subscript[Y, Row@Join[{d, ","}, Riffle[formattedIndices, ","]]]
      ],

    "Yu",
      If[
        formattedIndices === {},
        Subscript[Y, u],
        Subscript[Y, Row@Join[{u, ","}, Riffle[formattedIndices, ","]]]
      ],

    "Ye",
      If[
        formattedIndices === {},
        Subscript[Y, e],
        Subscript[Y, Row@Join[{e, ","}, Riffle[formattedIndices, ","]]]
      ],

    "yNew",
      If[
        formattedIndices === {},
        Subscript[y, N],
        Subscript[y, Row@Join[{N, ","}, Riffle[formattedIndices, ","]]]
      ],

    _,
      With[{base = PhysicsCouplingBase[name]},
        If[
          formattedIndices === {},
          base,
          (* For any other indexed coupling, avoid nesting Subscript. *)
          If[
            Head[Unevaluated[base]] === Subscript,
            With[{parts = List @@ Unevaluated[base]},
              Subscript[
                parts[[1]],
                Row@Join[{parts[[2]], ","}, Riffle[formattedIndices, ","]]
              ]
            ],
            Subscript[base, indexRow]
          ]
        ]
      ]
  ]
];


(* ---------------------------------------------------------------------- *)
(* Group tensors.                                                          *)
(* ---------------------------------------------------------------------- *)

PhysicsCG[cg_] := Module[
  {arguments, tensor, indices, formatted, tensorText},

  arguments = List @@ Unevaluated[cg];
  If[Length[arguments] < 2, Return[cg]];

  tensor = arguments[[1]];
  indices = arguments[[2]];
  formatted = PhysicsIndices[indices];
  tensorText = ToLowerCase@ToString[Unevaluated[tensor], InputForm];

  Which[
    StringContainsQ[tensorText, "bar["] &&
      StringContainsQ[tensorText, "eps["],
      Superscript[\[Epsilon], Row[formatted]],

    StringContainsQ[tensorText, "eps["],
      Subscript[\[Epsilon], Row[formatted]],

    True,
      Subscript[
        PhysicsDisplayForm[tensor],
        Row[Riffle[formatted, ","]]
      ]
  ]
];


(* ---------------------------------------------------------------------- *)
(* Gamma matrices, projectors and complete fermion chains.                  *)
(* ---------------------------------------------------------------------- *)

PhysicsDiracFactor[factor_] := Module[{headName, args},
  headName = PhysicsSymbolName[Head[Unevaluated[factor]]];
  args = If[AtomQ[Unevaluated[factor]], {}, List @@ Unevaluated[factor]];

  Switch[headName,
    "Proj",
      If[args === {-1}, Subscript[P, L],
        If[args === {1}, Subscript[P, R], Subscript[P, Row[args]]]
      ],

    "GammaM",
      If[
        args === {},
        \[Gamma],
        Superscript[\[Gamma], PhysicsIndex[First[args]]]
      ],

    "GammaCC",
      C,

    "DiracProduct",
      Row[PhysicsDiracFactor /@ args, "\[ThinSpace]"],

    _,
      PhysicsDisplayForm[factor]
  ]
];

PhysicsFermionChain[chain_] := Module[{factors},
  factors = List @@ Unevaluated[chain];
  Row[PhysicsDiracFactor /@ factors, "\[ThinSpace]"]
];

(* ---------------------------------------------------------------------- *)
(* Gauge-field strengths.                                                  *)
(* ---------------------------------------------------------------------- *)

PhysicsFieldStrength[fieldStrength_] := Module[
  {arguments, name, lorentzIndices, gaugeIndices, lower, upper},

  arguments = List @@ Unevaluated[fieldStrength];
  If[Length[arguments] < 2, Return[fieldStrength]];

  name = arguments[[1]];
  lorentzIndices = arguments[[2]];
  gaugeIndices = If[Length[arguments] >= 3, arguments[[3]], {}];

  lower = Row[PhysicsIndices[lorentzIndices]];
  upper = If[
    ListQ[gaugeIndices] && gaugeIndices =!= {},
    Row[PhysicsIndices[gaugeIndices]],
    ""
  ];

  If[
    upper === "",
    Subscript[name, lower],
    Subsuperscript[name, lower, upper]
  ]
];

(* ---------------------------------------------------------------------- *)
(* Recursive dispatcher. It converts outer structures first and then their *)
(* contents, avoiding the partial-conversion problem of a single rule list. *)
(* ---------------------------------------------------------------------- *)

PhysicsDisplayForm[expression_] := Module[{headName, arguments},
  If[AtomQ[Unevaluated[expression]],
    Return[
      Switch[PhysicsSymbolName[Unevaluated[expression]],
        "hbar", \[HBar],
        "mu2", Superscript[\[Mu], 2],
        "\[Mu]2", Superscript[\[Mu], 2],
        "eps", \[Epsilon],
        "mubar", OverBar[\[Mu]],
        "mubar2", Superscript[OverBar[\[Mu]], 2],
        "muBar", OverBar[\[Mu]],
        "muBar2", Superscript[OverBar[\[Mu]], 2],
        "GammaCC", C,
        _, expression
      ]
    ]
  ];

  headName = PhysicsSymbolName[Head[Unevaluated[expression]]];
  arguments = List @@ Unevaluated[expression];

  Switch[headName,
    "Plus",
      Plus @@ (PhysicsDisplayForm /@ arguments),

    "Times",
      Times @@ (PhysicsDisplayForm /@ arguments),

    "Power",
      Power[
        PhysicsDisplayForm[arguments[[1]]],
        PhysicsDisplayForm[arguments[[2]]]
      ],

    "Rational",
      expression,

    "Complex",
      expression,

    "Field",
      PhysicsField[expression],

    "Bar",
      PhysicsBar[arguments[[1]]],

    "Transp",
      PhysicsTranspose[arguments[[1]]],

    "Coupling",
      PhysicsCoupling[expression],

    "CG",
      PhysicsCG[expression],

    "NCM",
      PhysicsFermionChain[expression],

    "NonCommutativeMultiply",
      PhysicsFermionChain[expression],

    "DiracProduct",
      PhysicsDiracFactor[expression],

    "Proj",
      PhysicsDiracFactor[expression],

    "GammaM",
      PhysicsDiracFactor[expression],

    "GammaCC",
      C,

    "FieldStrength",
      PhysicsFieldStrength[expression],

    "Log",
      Log[PhysicsDisplayForm[arguments[[1]]]],

    "Exp",
      Exp[PhysicsDisplayForm[arguments[[1]]]],

    "eps",
      \[Epsilon],

    _,
      Apply[
        Head[Unevaluated[expression]],
        PhysicsDisplayForm /@ arguments
      ]
  ]
];

(* ---------------------------------------------------------------------- *)
(* Validation. Any surviving Matchete-like object prevents PDF export.      *)
(* ---------------------------------------------------------------------- *)

RemainingInternalObjects[expression_] := DeleteDuplicates@Cases[
  Unevaluated[expression],
  object_ /; MemberQ[
    {
      "Field", "Coupling", "Index", "Bar", "Transp", "NCM",
      "DiracProduct", "GammaM", "GammaCC", "Proj", "CG",
      "FieldStrength"
    },
    PhysicsSymbolName[Head[Unevaluated[object]]]
  ] :> HoldForm[object],
  Infinity
];

FullyConvertedQ[expression_] :=
  RemainingInternalObjects[expression] === {};

ExpressionToLaTeX[expression_] := Module[
  {formatted, remaining, latex},

  If[
    expression === Missing["NotAvailable"] ||
    expression === $Failed ||
    expression === Null,
    Return[<|
      "Success" -> False,
      "LaTeX" -> "",
      "Remaining" -> "Expression not available."
    |>]
  ];

  formatted = PhysicsDisplayForm[Unevaluated[expression]];

  (* Final index cleanup. Some Matchete indices can survive inside held,
     display, or tensor structures even after recursive formatting. Convert
     them everywhere before validating the result. *)
  formatted = formatted /. {
    HoldPattern[bar_[idx_] /;
      PhysicsSymbolName[Unevaluated[bar]] === "Bar" &&
      PhysicsSymbolName[Head[Unevaluated[idx]]] === "Index"
    ] :> PhysicsIndex[idx],

    HoldPattern[idx_ /;
      PhysicsSymbolName[Head[Unevaluated[idx]]] === "Index"
    ] :> PhysicsIndex[idx]
  };

  remaining = RemainingInternalObjects[formatted];

  If[
    remaining =!= {},
    Return[<|
      "Success" -> False,
      "LaTeX" -> "",
      "Remaining" -> ToString[Short[remaining, 20], InputForm]
    |>]
  ];

  latex = Quiet@Check[
    ToString[TeXForm[formatted]],
    ""
  ];

  (* TeXForm's rendering of OverBar varies between front-end/kernel versions.
     Normalise every known textual form to literal LaTeX so mubar2 is always
     emitted as \bar{\mu}^{2}. *)
  If[StringQ[latex],
    latex = StringReplace[latex, {
      "\\mu\\text{bar2}" -> "\\bar{\\mu}^{2}",
      "{\\mu}\\text{bar2}" -> "\\bar{\\mu}^{2}",
      "\\mu\\text{Bar2}" -> "\\bar{\\mu}^{2}",
      "{\\mu}\\text{Bar2}" -> "\\bar{\\mu}^{2}",
      "\\text{mubar2}" -> "\\bar{\\mu}^{2}",
      "\\text{muBar2}" -> "\\bar{\\mu}^{2}",
      "\\mathrm{mubar2}" -> "\\bar{\\mu}^{2}",
      "\\mathrm{muBar2}" -> "\\bar{\\mu}^{2}",
      "mubar2" -> "\\bar{\\mu}^{2}",
      "muBar2" -> "\\bar{\\mu}^{2}",
      "\\overline{\\mu}^{2}" -> "\\bar{\\mu}^{2}",
      "\\overline{\\mu^2}" -> "\\bar{\\mu}^{2}"
    }]
  ];

  <|
    "Success" -> StringQ[latex] && StringLength[latex] > 0,
    "LaTeX" -> latex,
    "Remaining" -> If[latex === "", "TeXForm failed.", ""]
  |>
];
