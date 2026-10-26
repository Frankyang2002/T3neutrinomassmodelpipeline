(* Experimental extension of the ordinary T3 scalar operator basis.

   Scalar species are ordered H, Hdag, S1, S1dag, S2, S2dag.
   All extra degree-four, gauge- and Z2-invariant monomials are considered.
   Hermitian-conjugate pairs are counted once and given a complex coupling.
   The pre-existing scalar structures and the T3 topology coupling are not
   redefined.  Distinct SU(2) contractions are projected onto the fully
   symmetric bosonic tensor space BEFORE testing independence.

   Only quartic interactions can currently be included in the UV Lagrangian.
   Quadratic mixing terms are catalogued, but are deliberately NOT inserted
   into the diagonal-mass Matchete/RGBeta/matching pipeline.

   Opt in for a Matchete matching experiment by setting the environment
   variable T3_EXPERIMENTAL_GENERAL_QUARTICS=1 before running RunModel.wl.
   This flag does NOT extend RGBeta or the downstream numerical boundary data.
*)

ClearAll[
  T3ExtendedSpecies,
  T3ExtendedConjugateCounts,
  T3ExtendedCanonicalCounts,
  T3ExtendedStructureKey,
  T3ExtendedSpeciesList,
  T3ExtendedSpeciesHypercharges2,
  T3ExtendedAllowedStructureQ,
  T3ExtendedCandidateStructures,
  T3ExtendedBaselineQuartics,
  T3ExtendedBaselineQuadratics,
  T3ExtendedExtraQuarticStructures,
  T3ExtendedVerifiedStructureCount,
  T3ExtendedStructureExpectedMultiplicities,
  T3ExtendedExtraQuadraticStructures,
  T3ExtendedFieldDimensions,
  T3ExtendedFieldConjugatedQ,
  T3ExtendedFieldValue,
  T3ExtendedSymmetrizeGroup,
  T3ExtendedBosonicTensorBasis,
  T3ExtendedDefineComplexFamily,
  T3ExtendedBuildCandidates,
  T3ExtendedBuildExtraQuarticCandidates,
  T3ExtendedDescribeQuadraticOperators,
  T3ExperimentalGeneralQuarticsQ
];

(* Preserve the exact source-relative path during Get, not the caller's path. *)
Get[FileNameJoin[{DirectoryName @ ExpandFileName[$InputFileName],
  "GeneralScalarVerifiedQuartics.wl"}]];

T3ExtendedSpecies = {"H", "Hdag", "S1", "S1dag", "S2", "S2dag"};

T3ExtendedConjugateCounts[counts_List] := counts[[{2, 1, 4, 3, 6, 5}]];

T3ExtendedCanonicalCounts[counts_List] := First @ Sort[{
  counts, T3ExtendedConjugateCounts[counts]
}];

T3ExtendedStructureKey[counts_List] := StringRiffle[ToString /@ counts, "_"];

T3ExtendedSpeciesList[counts_List] := Flatten @ MapThread[
  ConstantArray,
  {T3ExtendedSpecies, counts}
];

(* Twice the hypercharges; all coefficients are integers for T3 alpha. *)
T3ExtendedSpeciesHypercharges2[model_Association] := With[
  {y1 = 2 model["Scalar1", "Y"], y2 = 2 model["Scalar2", "Y"]},
  {1, -1, y1, -y1, y2, -y2}
];

T3ExtendedAllowedStructureQ[
  model_Association, counts_List, degree_Integer
] := TrueQ[
  Length[counts] === 6 &&
  Total[counts] === degree &&
  AllTrue[counts, IntegerQ[#] && NonNegative[#] &] &&
  EvenQ[Total[counts[[3 ;; 6]]]] &&
  Total[counts T3ExtendedSpeciesHypercharges2[model]] === 0 &&
  (* At least one BSM scalar is needed. SM-only terms belong to LoadModel["SM"]. *)
  Total[counts[[3 ;; 6]]] > 0
];

T3ExtendedCandidateStructures[model_Association, degree_Integer] :=
  Sort @ DeleteDuplicates @ (
    T3ExtendedCanonicalCounts /@
    Select[
      Tuples[Range[0, degree], 6],
      T3ExtendedAllowedStructureQ[model, #, degree] &
    ]
  );

(* The pre-existing generator includes all physical contractions of the
   following monomials.  H-only is supplied by the Standard Model, while
   H H S1 S2dag is handled by BuildT3MixingCandidates. *)
T3ExtendedBaselineQuartics = (
  T3ExtendedCanonicalCounts /@ {
    {0, 0, 2, 2, 0, 0},
    {0, 0, 0, 0, 2, 2},
    {1, 1, 1, 1, 0, 0},
    {1, 1, 0, 0, 1, 1},
    {0, 0, 1, 1, 1, 1},
    {2, 0, 1, 0, 0, 1}
  }
);

T3ExtendedBaselineQuadratics = (
  T3ExtendedCanonicalCounts /@ {
    {0, 0, 1, 1, 0, 0},
    {0, 0, 0, 0, 1, 1}
  }
);

(* For the 16 audited production models, use the independently
   enumerated nonzero-SU2-invariant structures as the source of truth.
   A raw hypercharge-only enumeration includes structures whose fully
   Bose-symmetrized SU(2) invariant tensor vanishes.

   The original dynamic path is retained only for unsupported DIMS points.
*)
T3ExtendedExtraQuarticStructures[model_Association] := Module[{verified},
  verified = T3ExtendedVerifiedEntries[model];
  If[ListQ[verified],
    Return[verified[[All, 1]]]
  ];
  Complement[
    T3ExtendedCandidateStructures[model, 4],
    T3ExtendedBaselineQuartics
  ]
];

T3ExtendedVerifiedStructureCount[model_Association] := Module[{verified},
  verified = T3ExtendedVerifiedEntries[model];
  If[!ListQ[verified], Return[Missing["UncertifiedModel"]]];
  Total[verified[[All, 2]]]
];

T3ExtendedStructureExpectedMultiplicities[model_Association] := Module[{verified},
  verified = T3ExtendedVerifiedEntries[model];
  If[!ListQ[verified], Return[Missing["UncertifiedModel"]]];
  AssociationThread[
    T3ExtendedStructureKey /@ verified[[All, 1]],
    verified[[All, 2]]
  ]
];

T3ExtendedExtraQuadraticStructures[model_Association] := Complement[
  T3ExtendedCandidateStructures[model, 2],
  T3ExtendedBaselineQuadratics
];

T3ExtendedFieldDimensions[model_Association, names_List] := Replace[
  names,
  {
    "H" | "Hdag" -> 2,
    "S1" | "S1dag" -> model["Scalar1", "SU2"],
    "S2" | "S2dag" -> model["Scalar2", "SU2"]
  },
  {1}
];

T3ExtendedFieldConjugatedQ[name_String] := StringEndsQ[name, "dag"];

T3ExtendedFieldValue[
  name_String, d_Integer?Positive, index_
] := Switch[
  name,
  "H", H[index],
  "Hdag", Bar[H[index]],
  "S1", ScalarPotentialFieldValue[1, False, d, index],
  "S1dag", ScalarPotentialFieldValue[1, True, d, index],
  "S2", ScalarPotentialFieldValue[2, False, d, index],
  "S2dag", ScalarPotentialFieldValue[2, True, d, index]
];

(* Average over the FULL symmetric group of an identical boson species.
   Successive pairwise projections do not suffice for three or four fields. *)
T3ExtendedSymmetrizeGroup[tensor_, positions_List] := Module[
  {rank = ArrayDepth[tensor], permutations},
  If[Length[positions] < 2, Return[tensor]];
  permutations = Permutations[positions];
  Total[
    (Transpose[
      tensor,
      ReplacePart[Range[rank], Thread[positions -> #]]
    ] &) /@ permutations
  ] / Length[permutations]
];

(* Use the existing Matchete/SU(2) invariant construction, then project
   completely onto Sym^n(R) for EACH repeated field species. *)
T3ExtendedBosonicTensorBasis[model_Association, names_List] := Module[
  {dims, isBar, base, keep, originalGroups, positionMap, groups, basis},

  dims = T3ExtendedFieldDimensions[model, names];
  isBar = T3ExtendedFieldConjugatedQ /@ names;

  base = PhysicalSU2InvariantBasis[
    dims,
    isBar,
    {},
    Flatten @ Position[names, "H" | "Hdag"]
  ];
  If[base === $Failed, Return[$Failed]];

  keep = base["Keep"];
  positionMap = AssociationThread[keep -> Range[Length[keep]]];
  originalGroups = Select[
    GatherBy[Range[Length[names]], names[[#]] &],
    Length[#] > 1 &
  ];
  groups = Select[
    (Lookup[positionMap, Intersection[keep, #]] &) /@ originalGroups,
    Length[#] > 1 &
  ];

  basis = IndependentTensorBasis[
    Fold[T3ExtendedSymmetrizeGroup[#1, #2] &, #, groups] & /@
      base["Basis"]
  ];

  If[basis === {}, Return[$Failed]];
  Association[base, "Basis" -> basis, "PhysicalCount" -> Length[basis]]
];

(* This function registers independent COMPLEX couplings for monomials that
   are distinct from their own Hermitian conjugates.  Every extra quartic
   structure in the current 16 supported ordinary models has this property. *)
T3ExtendedDefineComplexFamily[
  data_Association, counts_List, prefix_String
] := Module[
  {basis, cgReps, cgNames = {}, couplings = {}, cg, coupling, key},

  basis = data["Basis"];
  cgReps = data["CGReps"];
  key = StringJoin[ToString /@ counts];

  Do[
    If[cgReps === {},
      AppendTo[cgNames, None],
      cg = Symbol["T3ExtraCG" <> prefix <> key <> "I" <> ToString[i]];
      If[Check[DefineCG[cg, cgReps, basis[[i]]]; True, $Failed] === $Failed,
        Return[$Failed]
      ];
      AppendTo[cgNames, cg]
    ];

    coupling = Symbol["lambdaT3Extra" <> prefix <> key <> "I" <> ToString[i]];
    If[Check[DefineCoupling[coupling, SelfConjugate -> False]; True, $Failed] === $Failed,
      Return[$Failed]
    ];
    AppendTo[couplings, coupling],
    {i, Length[basis]}
  ];

  <|"CGs" -> cgNames, "Couplings" -> couplings|>
];

(* Do not silently drop structures if any SU(2) basis or Matchete coupling
   definition fails. An incomplete basis is a failed audit, not a successful
   audit with fewer operators. *)
T3ExtendedBuildCandidates[
  model_Association, degree_Integer, structures_List
] := Module[
  {results = {}, counts, names, data, family, dims, indices, labels,
   product, cg, coupling, contribution, prefix, multiplicities, expected},

  multiplicities = If[degree === 4,
    T3ExtendedStructureExpectedMultiplicities[model], Missing["NotQuartic"]
  ];
  prefix = If[degree === 4, "Q", "M"];

  Do[
    counts = structures[[structureIndex]];
    If[counts === T3ExtendedConjugateCounts[counts],
      Print["ERROR: extra self-adjoint scalar structure needs a real-coupling basis: ", counts];
      Return[$Failed]
    ];

    names = T3ExtendedSpeciesList[counts];
    dims = T3ExtendedFieldDimensions[model, names];
    data = T3ExtendedBosonicTensorBasis[model, names];
    If[data === $Failed,
      Print["ERROR: SU(2) invariant basis failed for ", names];
      Return[$Failed]
    ];

    expected = If[AssociationQ[multiplicities],
      Lookup[multiplicities, T3ExtendedStructureKey[counts], Missing["UnknownStructure"]],
      Missing["Uncertified"]
    ];
    If[IntegerQ[expected] && data["PhysicalCount"] =!= expected,
      Print["ERROR: expected ", expected, " SU(2) singlet(s) for ", names,
        ", but Matchete constructed ", data["PhysicalCount"]];
      Return[$Failed]
    ];
    Print["General scalar structure ", names,
      ": SU(2) singlets=", data["PhysicalCount"],
      ", expected=", expected];

    family = T3ExtendedDefineComplexFamily[data, counts, prefix];
    If[family === $Failed,
      Print["ERROR: CG or coupling registration failed for ", names];
      Return[$Failed]
    ];

    indices = Array[Unique["t3extra"] &, degree];
    labels = indices[[data["Keep"]]];
    product = Times @@ MapThread[
      T3ExtendedFieldValue,
      {names, dims, indices}
    ];

    Do[
      cg = family["CGs"][[i]];
      coupling = family["Couplings"][[i]];
      contribution = coupling[] product If[cg === None, 1, CG[cg, labels]];
      AppendTo[results,
        <|
          "Name" -> "T3Extra" <> prefix <> "_" <>
            T3ExtendedStructureKey[counts] <> "_" <> ToString[i],
          "Class" -> If[degree === 4, "ExtraQuartic", "ExtraQuadratic"],
          "Expression" -> PlusHc[contribution],
          "Counts" -> counts
        |>
      ],
      {i, data["PhysicalCount"]}
    ],
    {structureIndex, Length[structures]}
  ];

  results
];

T3ExtendedBuildExtraQuarticCandidates[model_Association] := Module[
  {candidates, expected},
  candidates = T3ExtendedBuildCandidates[
    model, 4, T3ExtendedExtraQuarticStructures[model]
  ];
  If[!ListQ[candidates], Return[$Failed]];
  expected = T3ExtendedVerifiedStructureCount[model];
  If[IntegerQ[expected] && Length[candidates] =!= expected,
    Print["ERROR: incomplete scalar-quartic model: generated ",
      Length[candidates], " candidates, expected ", expected,
      " from 16-model SU(2) operator manifest."];
    Return[$Failed]
  ];
  candidates
];

(* Catalogue only: these terms need a non-diagonal mass treatment before
   they can be enabled in a production Lagrangian.  No Matchete couplings
   are registered by this reporting function. *)
T3ExtendedDescribeQuadraticOperators[model_Association] := Module[
  {counts, names, data},
  Cases[
    Table[
      names = T3ExtendedSpeciesList[counts];
      data = T3ExtendedBosonicTensorBasis[model, names];
      If[data === $Failed,
        Nothing,
        <|
          "Structure" -> names,
          "Counts" -> counts,
          "IndependentComplexCouplings" -> data["PhysicalCount"]
        |>
      ],
      {counts, T3ExtendedExtraQuadraticStructures[model]}
    ],
    _Association
  ]
];

T3ExperimentalGeneralQuarticsQ[] :=
  Environment["T3_EXPERIMENTAL_GENERAL_QUARTICS"] === "1";
