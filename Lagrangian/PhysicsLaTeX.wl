(* ::Package:: *)

(*
  PhysicsLaTeX.wl
  
  Convert Matchete expressions into readable physics notation before TeX export.

  The conversion is intentionally recursive: field, index, coupling and Dirac
  structures are formatted independently, then the final expression is checked
  for any Matchete objects that escaped conversion.
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
  PhysicsCGBase,
  PhysicsCG,
  PhysicsDiracFactor,
  PhysicsFermionChain,
  PhysicsFieldStrength,
  PhysicsDisplayForm,
  RemainingInternalObjects,
  FullyConvertedQ,
  ExpressionToLaTeX
];

(* Symbol-name matching keeps the formatter independent of Matchete contexts. *)

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

(* Fields and derivatives. *)

PhysicsFieldBase[name_] := Switch[PhysicsSymbolName[name],
  "H", H,
  "NewScalar", \[Eta],
  "NewScalar1", Subscript[S, 1],
  "NewScalar2", Subscript[S, 2],
  "q", q,
  "l", L,
  "u", u,
  "d", d,
  "e", e,
  "N", F,
  "NewFermion", F,
  "F", F,
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

  (* NewScalar1/NewScalar2 already have physical labels S_1/S_2.  Construct
     their complete subscript in one operation instead of ever forming the
     nested object Subscript[Subscript[S,1],j]. *)
  result = Switch[
    PhysicsSymbolName[name],

    "NewScalar1",
      If[
        formattedInternal === {},
        Subscript[S, 1],
        Subscript[
          S,
          Row@Join[{1, ","}, Riffle[formattedInternal, ","]]
        ]
      ],

    "NewScalar2",
      If[
        formattedInternal === {},
        Subscript[S, 2],
        Subscript[
          S,
          Row@Join[{2, ","}, Riffle[formattedInternal, ","]]
        ]
      ],

    _,
      If[
        formattedInternal === {},
        result,
        Subscript[result, Row[Riffle[formattedInternal, ","]]]
      ]
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

(* Couplings and masses. *)

PhysicsCouplingBase[name_] := Module[
  {nameString, match},

  nameString = PhysicsSymbolName[name];

  Which[
    nameString === "gY", Subscript[g, Y],
    nameString === "gL" || nameString === "g2", Subscript[g, 2],
    nameString === "gs" || nameString === "g3", Subscript[g, 3],

    nameString === "lambda", \[Lambda],
    nameString === "\[Lambda]", \[Lambda],
    nameString === "lambda3", Subscript[\[Lambda], 3],
    nameString === "lambda4", Subscript[\[Lambda], 4],
    nameString === "lambda5", Subscript[\[Lambda], 5],
    nameString === "lambdaS", Subscript[\[Lambda], S],
    nameString === "lambdaT3", Subscript[\[Lambda], T3],
    nameString === "lambdaS1", Subscript[\[Lambda], S1],
    nameString === "lambdaS2", Subscript[\[Lambda], S2],
    nameString === "lambdaH1", Subscript[\[Lambda], H1],
    nameString === "lambdaH2", Subscript[\[Lambda], H2],
    nameString === "lambda12", Subscript[\[Lambda], 12],

    StringMatchQ[nameString, "lambdaH1Inv" ~~ DigitCharacter ..],
      match = StringCases[
        nameString,
        "lambdaH1Inv" ~~ n : DigitCharacter .. :> n
      ];
      Superscript[Subscript[\[Lambda], H1], First[match]],

    StringMatchQ[nameString, "lambdaH2Inv" ~~ DigitCharacter ..],
      match = StringCases[
        nameString,
        "lambdaH2Inv" ~~ n : DigitCharacter .. :> n
      ];
      Superscript[Subscript[\[Lambda], H2], First[match]],

    StringMatchQ[nameString, "lambdaS1Inv" ~~ DigitCharacter ..],
      match = StringCases[
        nameString,
        "lambdaS1Inv" ~~ n : DigitCharacter .. :> n
      ];
      Superscript[Subscript[\[Lambda], S1], First[match]],

    StringMatchQ[nameString, "lambdaS2Inv" ~~ DigitCharacter ..],
      match = StringCases[
        nameString,
        "lambdaS2Inv" ~~ n : DigitCharacter .. :> n
      ];
      Superscript[Subscript[\[Lambda], S2], First[match]],

    StringMatchQ[nameString, "lambda12Inv" ~~ DigitCharacter ..],
      match = StringCases[
        nameString,
        "lambda12Inv" ~~ n : DigitCharacter .. :> n
      ];
      Superscript[Subscript[\[Lambda], 12], First[match]],

    nameString === "lambdaH1Adj", Superscript[Subscript[\[Lambda], H1], Adj],
    nameString === "lambdaH2Adj", Superscript[Subscript[\[Lambda], H2], Adj],
    nameString === "lambdaS1Adj", Superscript[Subscript[\[Lambda], S1], Adj],
    nameString === "lambdaS2Adj", Superscript[Subscript[\[Lambda], S2], Adj],
    nameString === "lambda12Adj", Superscript[Subscript[\[Lambda], 12], Adj],
    nameString === "lambda12Cross", Superscript[Subscript[\[Lambda], 12], Cross],

    nameString === "y1", Subscript[y, 1],
    nameString === "y2", Subscript[y, 2],
    nameString === "yNew", Subscript[y, N],
    nameString === "Yd", Subscript[Y, d],
    nameString === "Yu", Subscript[Y, u],
    nameString === "Ye", Subscript[Y, e],

    nameString === "MEta", Subscript[M, \[Eta]],
    nameString === "MN", Subscript[M, N],
    nameString === "MF", Subscript[M, F],
    nameString === "MS1", Subscript[M, S1],
    nameString === "MS2", Subscript[M, S2],
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
            Head[base] === Subscript,
            With[{parts = List @@ base},
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


(* Group tensors. *)

(* Human-facing names for invariant tensors.  We keep the invariant
   label because different contractions are physically independent, but hide
   implementation names such as T3HiggsPortal1CGInv2. *)
PhysicsCGBase[tensor_] := Module[
  {nameString, match},

  nameString = PhysicsSymbolName[tensor];

  Which[
    nameString === "T3Y1CG",
      Subscript[\[ScriptCapitalI], y1],

    nameString === "T3Y2CG",
      Subscript[\[ScriptCapitalI], y2],

    nameString === "T3MixCG",
      Subscript[\[ScriptCapitalI], T3],

    StringMatchQ[nameString, "T3HiggsPortal1CGInv" ~~ DigitCharacter ..],
      match = StringCases[
        nameString,
        "T3HiggsPortal1CGInv" ~~ n : DigitCharacter .. :> n
      ];
      Superscript[Subscript[\[ScriptCapitalI], H1], First[match]],

    StringMatchQ[nameString, "T3HiggsPortal2CGInv" ~~ DigitCharacter ..],
      match = StringCases[
        nameString,
        "T3HiggsPortal2CGInv" ~~ n : DigitCharacter .. :> n
      ];
      Superscript[Subscript[\[ScriptCapitalI], H2], First[match]],

    StringMatchQ[nameString, "T3Scalar1SelfCGInv" ~~ DigitCharacter ..],
      match = StringCases[
        nameString,
        "T3Scalar1SelfCGInv" ~~ n : DigitCharacter .. :> n
      ];
      Superscript[Subscript[\[ScriptCapitalI], S1], First[match]],

    StringMatchQ[nameString, "T3Scalar2SelfCGInv" ~~ DigitCharacter ..],
      match = StringCases[
        nameString,
        "T3Scalar2SelfCGInv" ~~ n : DigitCharacter .. :> n
      ];
      Superscript[Subscript[\[ScriptCapitalI], S2], First[match]],

    StringMatchQ[nameString, "T3CrossScalarCGInv" ~~ DigitCharacter ..],
      match = StringCases[
        nameString,
        "T3CrossScalarCGInv" ~~ n : DigitCharacter .. :> n
      ];
      Superscript[Subscript[\[ScriptCapitalI], 12], First[match]],

    True,
      tensor
  ]
];

PhysicsCG[cg_] := Module[
  {
    arguments, tensor, indices, formatted, tensorText,
    base, baseParts, inner, innerParts, combinedSubscript
  },

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
      base = PhysicsCGBase[tensor];

      (* PhysicsCGBase intentionally gives invariant tensors a descriptive
         label such as I_y1, I_y2, I_T3, I_H1^n, ... .  Adding the component
         indices with another Subscript used to produce nested objects such as

             Subscript[Subscript[I, y1], a]

         which TeXForm renders as I_{y1}_{a} and pdflatex rejects with
         "Double subscript".

         Merge the descriptive label and component indices into one subscript:
             I_{y1,a}
         while preserving any invariant-multiplicity superscript. *)

      Which[
        Head[base] === Subscript,
          baseParts = List @@ base;
          combinedSubscript = Row@Join[
            {baseParts[[2]]},
            If[formatted === {}, {}, Join[{","}, Riffle[formatted, ","]]]
          ];
          Subscript[baseParts[[1]], combinedSubscript],

        Head[base] === Superscript &&
          Head[First[List @@ base]] === Subscript,
          baseParts = List @@ base;
          inner = baseParts[[1]];
          innerParts = List @@ inner;
          combinedSubscript = Row@Join[
            {innerParts[[2]]},
            If[formatted === {}, {}, Join[{","}, Riffle[formatted, ","]]]
          ];
          Superscript[
            Subscript[innerParts[[1]], combinedSubscript],
            baseParts[[2]]
          ],

        formatted === {},
          base,

        True,
          Subscript[
            base,
            Row[Riffle[formatted, ","]]
          ]
      ]
  ]
];


(* Dirac structures and complete fermion chains. *)

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

(* Gauge-field strengths. *)

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

(* Recursive dispatcher: convert outer structures, then their contents. *)

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
        name_ /; StringStartsQ[name, "lambda"], PhysicsCouplingBase[expression],
        "y1", Subscript[y, 1],
        "y2", Subscript[y, 2],
        "Yd", Subscript[Y, d],
        "Yu", Subscript[Y, u],
        "Ye", Subscript[Y, e],
        "MF", Subscript[M, F],
        "MS1", Subscript[M, S1],
        "MS2", Subscript[M, S2],
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

(* Validation: surviving Matchete objects mean the conversion is incomplete. *)

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

  (* Do not discard the whole expression merely because a small number of
     Matchete wrapper objects survive the pretty-printer.  The report layer
     only needs a readable/best-effort LaTeX representation in order to group
     and compare terms by field content.  Previously this early return changed
     every such expression into an empty string, which made the UV/EFT
     Lagrangian comparison tables completely empty.

     We still record the surviving internal objects in "Remaining" so failed
     conversions remain diagnosable, but we allow TeXForm to render the
     partially converted expression. *)
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
    "Remaining" -> Which[
      latex === "",
        "TeXForm failed.",
      remaining =!= {},
        ToString[Short[remaining, 20], InputForm],
      True,
        ""
    ]
  |>
];
