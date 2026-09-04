(* Scalar-potential tensor export helpers for T3 RGE tensor exchange. *)

Get[FileNameJoin[{DirectoryName[$InputFileName], "T3RGEComponentExport.wl"}]];

T3RGEHiggsSelfEntries[] := Module[
  {data, basis, keep, factors},

  data = PhysicalSU2InvariantBasis[
    {2, 2, 2, 2},
    {True, False, True, False},
    {{1, 3}, {2, 4}},
    {1, 2, 3, 4}
  ];

  If[data === $Failed, Return[{}]];

  basis = data["Basis"];
  keep = data["Keep"];
  factors = {
    {"H", True},
    {"H", False},
    {"H", True},
    {"H", False}
  };

  If[Length[basis] =!= 1,
    Print[
      "ERROR: expected one physical Higgs self-quartic invariant, found ",
      Length[basis], "."
    ];
    Return[{}]
  ];

  (* Project convention: V contains lambdaH/2 (H^dagger H)^2. *)
  T3RGEQuarticEntriesFromBasis[
    First[basis],
    keep,
    factors,
    Symbol["lambdaH"]/2
  ]
];

T3RGEScalarSelfEntries[
  which_Integer,
  d_Integer?Positive
] := Module[
  {data, basis, keep, couplingBase, coupling, factors},

  data = PhysicalSU2InvariantBasis[
    {d, d, d, d},
    {True, False, True, False},
    {{1, 3}, {2, 4}},
    {}
  ];

  If[data === $Failed, Return[{}]];

  basis = data["Basis"];
  keep = data["Keep"];
  couplingBase = "lambdaS" <> ToString[which];

  factors = {
    {"S" <> ToString[which], True},
    {"S" <> ToString[which], False},
    {"S" <> ToString[which], True},
    {"S" <> ToString[which], False}
  };

  Flatten[
    Table[
      coupling = Symbol[couplingBase <> "Inv" <> ToString[invariant]];

      T3RGEQuarticEntriesFromBasis[
        basis[[invariant]],
        keep,
        factors,
        coupling/2
      ],
      {invariant, Length[basis]}
    ],
    1
  ]
];

T3RGEHiggsPortalEntries[
  which_Integer,
  d_Integer?Positive
] := Module[
  {data, basis, keep, couplingBase, coupling, factors},

  data = PhysicalSU2InvariantBasis[
    {2, 2, d, d},
    {True, False, True, False},
    {},
    {1, 2}
  ];

  If[data === $Failed, Return[{}]];

  basis = data["Basis"];
  keep = data["Keep"];
  couplingBase = "lambdaH" <> ToString[which];

  factors = {
    {"H", True},
    {"H", False},
    {"S" <> ToString[which], True},
    {"S" <> ToString[which], False}
  };

  Flatten[
    Table[
      coupling = Symbol[couplingBase <> "Inv" <> ToString[invariant]];

      (* Production Lagrangian uses 1/2 PlusHc[...] for this sector. *)
      T3RGEHermitianQuarticEntriesFromBasis[
        basis[[invariant]],
        keep,
        factors,
        coupling/2
      ],
      {invariant, Length[basis]}
    ],
    1
  ]
];

T3RGECrossScalarEntries[
  d1_Integer?Positive,
  d2_Integer?Positive
] := Module[
  {data, basis, keep, coupling, factors},

  data = PhysicalSU2InvariantBasis[
    {d1, d1, d2, d2},
    {True, False, True, False},
    {},
    {}
  ];

  If[data === $Failed, Return[{}]];

  basis = data["Basis"];
  keep = data["Keep"];

  factors = {
    {"S1", True},
    {"S1", False},
    {"S2", True},
    {"S2", False}
  };

  Flatten[
    Table[
      coupling = Symbol["lambda12Inv" <> ToString[invariant]];

      T3RGEQuarticEntriesFromBasis[
        basis[[invariant]],
        keep,
        factors,
        coupling
      ],
      {invariant, Length[basis]}
    ],
    1
  ]
];


(* The topology-defining H H S1 S2^\dagger invariant has multiplicity one
   for valid T3 assignments, but use the same general invariant-basis logic
   rather than assuming an explicit epsilon/generator form. *)

T3RGEMixingEntries[
  d1_Integer?Positive,
  d2_Integer?Positive
] := Module[
  {data, basis, keep, factors},

  data = PhysicalSU2InvariantBasis[
    {2, 2, d1, d2},
    {False, False, False, True},
    {{1, 2}},
    {1, 2}
  ];

  If[data === $Failed, Return[{}]];

  basis = data["Basis"];
  keep = data["Keep"];

  If[Length[basis] =!= 1,
    Print[
      "ERROR: T3 H H S1 S2^dagger sector has ", Length[basis],
      " physical invariants; exactly one is required by the current topology."
    ];
    Return[{}]
  ];

  factors = {
    {"H", False},
    {"H", False},
    {"S1", False},
    {"S2", True}
  };

  Flatten[
    T3RGEHermitianQuarticEntriesFromBasis[
      #,
      keep,
      factors,
      lambdaT3
    ] & /@ basis,
    1
  ]
];


(* ------------------------------------------------------------------------- *)
(* Raw Yukawa invariant export                                               *)
(* ------------------------------------------------------------------------- *)

(* Export the invariant tensor and field orientation without choosing the
   global left-handed Weyl basis.  This keeps the group-theory result exact
   while making the remaining convention-sensitive step explicit. *)
