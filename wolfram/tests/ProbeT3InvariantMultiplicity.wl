(* ProbeT3InvariantMultiplicity.wl

   Diagnose the dimension of the SU(2) invariant tensor space for the three
   topology-defining T3 interactions without changing the model builder.

   Intended location:
     wolfram/tests/ProbeT3InvariantMultiplicity.wl

   Run from the project root with:
     wolframscript -file wolfram/tests/ProbeT3InvariantMultiplicity.wl

   Optional first argument sets the largest fermion SU(2) dimension scanned:
     wolframscript -file wolfram/tests/ProbeT3InvariantMultiplicity.wl 9
*)

ClearAll["Global`*"];
scriptDirectory = DirectoryName @ ExpandFileName[$InputFileName];

Print["Loading Matchete..."];
matcheteLoaded = UsingFrontEnd[Needs["Matchete`"]; True];
If[!TrueQ[matcheteLoaded],
  Print["ERROR: Matchete failed to load."];
  Exit[1]
];
Print["Matchete loaded successfully."];

(* Use the same Dynkin-label convention as LagrangianBuilder.wl. *)
SU2DynkinLabel[d_Integer?Positive] := {d - 1};

(* Check the same representation conditions used by pipeline.py for a T3 model. *)
ValidT3DimensionsQ[d1_Integer, d2_Integer, dF_Integer] := Module[
  {j1, j2},

  If[Min[d1, d2, dF] < 1, Return[False]];
  If[Abs[d1 - dF] =!= 1 || Abs[d2 - dF] =!= 1, Return[False]];

  j1 = (d1 - 1)/2;
  j2 = (d2 - 1)/2;

  Abs[j1 - j2] <= 1 <= j1 + j2 && IntegerQ[j1 + j2]
];

(* Return the complete invariant-tensor basis using the same conjugation and
   identical-field symmetry conventions as DefineSU2InvariantCG. *)
T3InvariantTensorBasis[
  dims_List,
  objectConjugated_List,
  symmetricPositions_: {}
] := Module[
  {keep, keptDims, keptConj, algebraReps, positionMap,
   symmetryAfterDrop, tensors},

  If[Length[dims] =!= Length[objectConjugated], Return[$Failed]];
  If[!AllTrue[dims, IntegerQ[#] && Positive[#] &], Return[$Failed]];

  keep = Flatten @ Position[dims, _?(# > 1 &)];
  keptDims = dims[[keep]];
  keptConj = objectConjugated[[keep]];

  If[keptDims === {}, Return[{1}]];

  algebraReps = MapThread[
    Function[{d, conjugated},
      Module[{label = SU2DynkinLabel[d]},
        If[EvenQ[d] && TrueQ[conjugated], CRep[label], label]
      ]
    ],
    {keptDims, keptConj}
  ];

  positionMap = AssociationThread[keep -> Range[Length[keep]]];
  symmetryAfterDrop = Select[
    Lookup[positionMap, #, Missing["Dropped"]] & /@ symmetricPositions,
    IntegerQ
  ];

  tensors = Check[
    If[symmetryAfterDrop === {},
      InvariantTensors[SU[2], algebraReps],
      InvariantTensors[
        SU[2],
        algebraReps,
        SymmetricIndices -> symmetryAfterDrop
      ]
    ],
    $Failed
  ];

  tensors
];

InvariantCount[basis_] := If[ListQ[basis], Length[basis], 0];

(* Count the three invariant spaces needed by the T3 diagram. *)
T3InvariantCounts[d1_Integer, d2_Integer, dF_Integer] := Module[
  {y1, y2, mix},

  y1 = T3InvariantTensorBasis[
    {2, dF, d1},
    {True, True, False}
  ];

  y2 = T3InvariantTensorBasis[
    {2, dF, d2},
    {True, False, True}
  ];

  mix = T3InvariantTensorBasis[
    {2, 2, d1, d2},
    {False, False, False, True},
    {1, 2}
  ];

  <|
    "Yukawa1" -> InvariantCount[y1],
    "Yukawa2" -> InvariantCount[y2],
    "Mixing" -> InvariantCount[mix],
    "MixingBasis" -> mix
  |>
];

maxF = If[
  Length[$ScriptCommandLine] >= 2,
  ToExpression @ Last[$ScriptCommandLine],
  9
];

If[!IntegerQ[maxF] || maxF < 1,
  Print["ERROR: maximum fermion dimension must be a positive integer."];
  Exit[2]
];

validDimensions = Select[
  Flatten[
    Table[
      {
        {dF - 1, dF - 1, dF},
        {dF - 1, dF + 1, dF},
        {dF + 1, dF - 1, dF},
        {dF + 1, dF + 1, dF}
      },
      {dF, 1, maxF}
    ],
    1
  ],
  ValidT3DimensionsQ @@ # &
];
validDimensions = DeleteDuplicates[validDimensions];

Print[""];
Print["Scanning valid T3 dimensions up to dF = ", maxF, "..."];
Print["dS1   dS2   dF    Y1   Y2   HH-S1-S2dagger"];
Print[StringRepeat["-", 49]];

results = Table[
  Module[{d1, d2, dF, counts},
    {d1, d2, dF} = dims;
    counts = T3InvariantCounts[d1, d2, dF];

    Print[
      Row[{
        d1, "      ", d2, "      ", dF, "     ",
        counts["Yukawa1"], "    ", counts["Yukawa2"], "    ",
        counts["Mixing"]
      }]
    ];

    <|
      "Dimensions" -> dims,
      "Yukawa1" -> counts["Yukawa1"],
      "Yukawa2" -> counts["Yukawa2"],
      "Mixing" -> counts["Mixing"],
      "MixingBasis" -> counts["MixingBasis"]
    |>
  ],
  {dims, validDimensions}
];

multipleMixing = Select[results, #["Mixing"] > 1 &];

Print[""];
If[multipleMixing === {},
  Print[
    "No valid T3 representation in this scan had more than one surviving ",
    "H H S1 S2^dagger invariant."
  ],
  Print["Representations with multiple surviving T3 mixing invariants:"];
  Scan[
    Function[result,
      Print[
        "  dims=", result["Dimensions"],
        " -> ", result["Mixing"], " mixing invariants"
      ]
    ],
    multipleMixing
  ];

  firstMultiple = First[multipleMixing];
  Print[""];
  Print["First regression candidate: ", firstMultiple["Dimensions"]];
  Print[
    "Mixing invariant basis (InputForm): ",
    InputForm[firstMultiple["MixingBasis"]]
  ];
];

If[AnyTrue[results, #["Yukawa1"] > 1 || #["Yukawa2"] > 1 &],
  Print[""];
  Print["NOTE: a Yukawa invariant space with multiplicity > 1 was also found."],
  Print[""];
  Print["All scanned T3 Yukawa invariant spaces had multiplicity <= 1."]
];
