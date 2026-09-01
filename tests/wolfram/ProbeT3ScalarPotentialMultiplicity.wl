(* ProbeT3ScalarPotentialMultiplicity.wl

   Diagnose the dimension of the SU(2) invariant tensor spaces in the scalar
   potential used by the generalised T3 builder, without changing the builder.

   Intended location:
     tests/wolfram/ProbeT3ScalarPotentialMultiplicity.wl

   Run from the project root with:
     wolframscript -file tests/wolfram/ProbeT3ScalarPotentialMultiplicity.wl

   Optional first argument sets the largest fermion SU(2) dimension scanned:
     wolframscript -file tests/wolfram/ProbeT3ScalarPotentialMultiplicity.wl 9
*)

ClearAll["Global`*"];

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

(* Return the raw invariant-tensor basis using the same conjugation convention
   as DefineSU2InvariantCG.  Singlet slots are removed because they carry no
   SU(2) index. *)
RawInvariantTensorData[dims_List, objectConjugated_List] := Module[
  {keep, keptDims, keptConj, algebraReps, tensors, positionMap},

  If[Length[dims] =!= Length[objectConjugated], Return[$Failed]];
  If[!AllTrue[dims, IntegerQ[#] && Positive[#] &], Return[$Failed]];

  keep = Flatten @ Position[dims, _?(# > 1 &)];
  keptDims = dims[[keep]];
  keptConj = objectConjugated[[keep]];
  positionMap = AssociationThread[keep -> Range[Length[keep]]];

  If[keptDims === {},
    Return[<|"Basis" -> {1}, "PositionMap" -> positionMap|>]
  ];

  algebraReps = MapThread[
    Function[{d, conjugated},
      Module[{label = SU2DynkinLabel[d]},
        If[EvenQ[d] && TrueQ[conjugated], CRep[label], label]
      ]
    ],
    {keptDims, keptConj}
  ];

  tensors = Check[
    InvariantTensors[SU[2], algebraReps],
    $Failed
  ];

  If[tensors === $Failed || !ListQ[tensors], Return[$Failed]];

  <|"Basis" -> tensors, "PositionMap" -> positionMap|>
];

(* Swap two tensor slots. *)
SwapTensorSlots[tensor_, first_Integer, second_Integer] := Module[
  {permutation},

  permutation = Range[ArrayDepth[tensor]];
  permutation[[{first, second}]] = permutation[[{second, first}]];
  Transpose[tensor, permutation]
];

(* Project a tensor onto the part symmetric under exchange of two identical
   bosonic field slots. *)
SymmetrizeTensorPair[tensor_, {first_Integer, second_Integer}] :=
  (tensor + SwapTensorSlots[tensor, first, second])/2;

(* Apply all independent identical-field symmetries and return the dimension
   of the surviving invariant space.  The rank is basis-independent. *)
PhysicalInvariantData[
  dims_List,
  objectConjugated_List,
  symmetricPairs_: {}
] := Module[
  {raw, basis, positionMap, mappedPairs, projected, vectors, rank},

  raw = RawInvariantTensorData[dims, objectConjugated];
  If[raw === $Failed, Return[$Failed]];

  basis = raw["Basis"];
  positionMap = raw["PositionMap"];

  mappedPairs = Cases[
    symmetricPairs,
    {a_Integer, b_Integer} /;
      KeyExistsQ[positionMap, a] && KeyExistsQ[positionMap, b] :>
        {positionMap[a], positionMap[b]}
  ];

  projected = Fold[
      SymmetrizeTensorPair[#1, #2] &,
      #,
      mappedPairs
    ] & /@ basis;

  vectors = (Flatten @ Normal[#]) & /@ projected;
  rank = If[vectors === {}, 0, MatrixRank[vectors]];

  <|
    "RawCount" -> Length[basis],
    "PhysicalCount" -> rank,
    "RawBasis" -> basis,
    "ProjectedBasis" -> projected,
    "SymmetricPairs" -> mappedPairs
  |>
];

(* Count the scalar-potential invariant spaces corresponding to the current
   ScalarSelf, HiggsPortal, and CrossScalarPortal interaction classes. *)
ScalarPotentialInvariantCounts[d1_Integer, d2_Integer] := Module[
  {self1, self2, portal1, portal2, cross},

  self1 = PhysicalInvariantData[
    {d1, d1, d1, d1},
    {True, False, True, False},
    {{1, 3}, {2, 4}}
  ];

  self2 = PhysicalInvariantData[
    {d2, d2, d2, d2},
    {True, False, True, False},
    {{1, 3}, {2, 4}}
  ];

  portal1 = PhysicalInvariantData[
    {2, 2, d1, d1},
    {True, False, True, False}
  ];

  portal2 = PhysicalInvariantData[
    {2, 2, d2, d2},
    {True, False, True, False}
  ];

  cross = PhysicalInvariantData[
    {d1, d1, d2, d2},
    {True, False, True, False}
  ];

  If[MemberQ[{self1, self2, portal1, portal2, cross}, $Failed],
    Return[$Failed]
  ];

  <|
    "ScalarSelf1" -> self1,
    "ScalarSelf2" -> self2,
    "HiggsPortal1" -> portal1,
    "HiggsPortal2" -> portal2,
    "CrossScalarPortal" -> cross
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
Print["Scanning scalar-potential invariant spaces up to dF = ", maxF, "..."];
Print["dS1   dS2   dF    Self1  Self2  HPortal1  HPortal2  Cross"];
Print[StringRepeat["-", 65]];

results = Table[
  Module[{d1, d2, dF, counts},
    {d1, d2, dF} = dims;
    counts = ScalarPotentialInvariantCounts[d1, d2];

    If[counts === $Failed,
      Print["FAILED for dims=", dims];
      <|"Dimensions" -> dims, "Status" -> "Failed"|>,

      Print[
        Row[{
          d1, "      ", d2, "      ", dF, "     ",
          counts["ScalarSelf1", "PhysicalCount"], "      ",
          counts["ScalarSelf2", "PhysicalCount"], "      ",
          counts["HiggsPortal1", "PhysicalCount"], "         ",
          counts["HiggsPortal2", "PhysicalCount"], "         ",
          counts["CrossScalarPortal", "PhysicalCount"]
        }]
      ];

      <|
        "Dimensions" -> dims,
        "Status" -> "Success",
        "Counts" -> counts
      |>
    ]
  ],
  {dims, validDimensions}
];

successfulResults = Select[results, #["Status"] === "Success" &];

multipleResults = Select[
  successfulResults,
  With[{counts = #["Counts"]},
    AnyTrue[
      {
        "ScalarSelf1", "ScalarSelf2", "HiggsPortal1",
        "HiggsPortal2", "CrossScalarPortal"
      },
      counts[#, "PhysicalCount"] > 1 &
    ]
  ] &
];

Print[""];
If[multipleResults === {},
  Print["No scalar-potential interaction with multiplicity > 1 was found."],

  Print["Representation assignments with multiple scalar-potential invariants:"];
  Scan[
    Function[result,
      Module[{dims, counts, multipleNames},
        dims = result["Dimensions"];
        counts = result["Counts"];
        multipleNames = Select[
          {
            "ScalarSelf1", "ScalarSelf2", "HiggsPortal1",
            "HiggsPortal2", "CrossScalarPortal"
          },
          counts[#, "PhysicalCount"] > 1 &
        ];

        Print["  dims=", dims];
        Scan[
          Function[name,
            Print[
              "    ", name,
              ": raw=", counts[name, "RawCount"],
              ", physical=", counts[name, "PhysicalCount"]
            ]
          ],
          multipleNames
        ]
      ]
    ],
    multipleResults
  ];

  firstMultiple = First[multipleResults];
  firstCounts = firstMultiple["Counts"];
  firstNames = Select[
    {
      "ScalarSelf1", "ScalarSelf2", "HiggsPortal1",
      "HiggsPortal2", "CrossScalarPortal"
    },
    firstCounts[#, "PhysicalCount"] > 1 &
  ];

  Print[""];
  Print["First regression candidate: ", firstMultiple["Dimensions"]];
  Print["Multiple interaction classes: ", firstNames];

  Scan[
    Function[name,
      Print[""];
      Print[name, " raw invariant basis (InputForm):"];
      Print[InputForm[firstCounts[name, "RawBasis"]]];
      If[firstCounts[name, "SymmetricPairs"] =!= {},
        Print[name, " projected symmetric basis (InputForm):"];
        Print[InputForm[firstCounts[name, "ProjectedBasis"]]]
      ]
    ],
    firstNames
  ];
];

