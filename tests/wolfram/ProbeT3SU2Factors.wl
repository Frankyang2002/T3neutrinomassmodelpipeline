(* ::Package:: *)

(*
  ProbeT3SU2Factors.wl

  Standalone diagnostic for the SU(2) group factor in the T3 Weinberg
  coefficient.  It does NOT run loop matching.

  For each original T3 class A-E it:
    1. reconstructs the exact invariant tensors used for y1, y2 and lambdaT3;
    2. contracts only the internal F, S1 and S2 SU(2) indices;
    3. leaves the two lepton and two Higgs doublet indices open;
    4. compares the resulting rank-4 tensor with the Weinberg epsilon-epsilon
       tensor;
    5. reports both the pure tensor ratio and the ratio including the factor 2
       from the two identical Higgs legs in H H S1 S2^\dagger.

  The expected coefficients are the ones observed independently in Matchete's
  matched C5 result.  Agreement therefore tests whether their representation
  dependence is already explained by the invariant-tensor contraction.
*)

ClearAll["Global`*"];

Print["Loading Matchete..."];
matcheteLoaded = UsingFrontEnd[Needs["Matchete`"]; True];
If[!TrueQ[matcheteLoaded],
  Print["ERROR: Matchete failed to load."];
  Exit[6]
];
Print["Matchete loaded successfully."];

(* Register the Standard Model gauge groups before defining any probe
   representations.  Loading the Matchete package alone does not create SU2L;
   the production pipeline first resets Matchete and loads the SM model. *)
ResetAll[];
smLoaded = UsingFrontEnd[LoadModel["SM"]];
If[smLoaded === $Failed,
  Print["ERROR: SM model failed to load, so SU2L is not registered."];
  Exit[7]
];
Print["SM model loaded; SU2L is registered."];

(* ------------------------------------------------------------------------- *)
(* Minimal representation helpers, matching SU2Invariants.wl                 *)
(* ------------------------------------------------------------------------- *)

SU2DynkinLabel[d_Integer?Positive] := {d - 1};

BSMRepresentationName[1] := None;
BSMRepresentationName[d_Integer?Positive] /; d >= 2 :=
  Symbol["T3ProbeBSMd" <> ToString[d]];

EnsureBSMRepresentation[1] := True;

(* This probe runs in a fresh kernel, so define each non-singlet probe
   representation once and memoise the result.  Do not use GetRepresentations[]
   as a success test here: in command-line Matchete sessions it can return
   $Failed even after DefineRepresentation has actually succeeded. *)
EnsureBSMRepresentation[d_Integer?Positive] /; d >= 2 :=
  EnsureBSMRepresentation[d] = Module[
    {rep, ok},
    rep = BSMRepresentationName[d];

    ok = Quiet@Check[
      DefineRepresentation[
        rep,
        SU2L,
        SU2DynkinLabel[d],
        IndexAlphabet -> {"u", "v", "w", "x", "y", "z"}
      ];
      True,
      False
    ];

    If[TrueQ[ok], True, $Failed]
  ];

(*
  Return the raw SparseArray produced by InvariantTensors together with the
  original field positions that survived after dropping SU(2) singlets.

  This deliberately mirrors DefineSU2InvariantCG in SU2Invariants.wl, but does
  not call DefineCG.  We want to inspect/contract the tensor itself.
*)
RawSU2Invariant[
  dims_List,
  objectConjugated_List,
  symmetricPositions_: {},
  smPositions_: {}
] := Module[
  {
    keep, keptDims, keptConj, algebraReps,
    positionMap, symmetryAfterDrop, tensors
  },

  If[Length[dims] =!= Length[objectConjugated], Return[$Failed]];
  If[!AllTrue[dims, IntegerQ[#] && Positive[#] &], Return[$Failed]];

  (* For this standalone probe we call InvariantTensors directly with SU(2)
     Dynkin labels.  No Matchete field representation needs to be registered:
     DefineRepresentation is only required later when a tensor is attached to
     actual Matchete fields via DefineCG. *)
  keep = Flatten@Position[dims, _?(# > 1 &)];
  keptDims = dims[[keep]];
  keptConj = objectConjugated[[keep]];

  If[keptDims === {},
    Return[
      <|
        "Tensor" -> 1,
        "Keep" -> {},
        "Dimensions" -> {},
        "AlgebraReps" -> {}
      |>
    ]
  ];

  algebraReps = MapThread[
    Function[{d, conjugated},
      Module[{label = SU2DynkinLabel[d]},
        If[EvenQ[d] && TrueQ[conjugated],
          CRep[label],
          label
        ]
      ]
    ],
    {keptDims, keptConj}
  ];

  positionMap = AssociationThread[keep -> Range[Length[keep]]];
  symmetryAfterDrop = Select[
    Lookup[positionMap, #, Missing["Dropped"]] & /@ symmetricPositions,
    IntegerQ
  ];

  tensors = Quiet@Check[
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

  If[
    tensors === $Failed ||
    !ListQ[tensors] ||
    Length[tensors] =!= 1,
    Return[
      <|
        "Failure" -> True,
        "TensorCount" -> If[ListQ[tensors], Length[tensors], Missing["Failed"]],
        "Keep" -> keep,
        "AlgebraReps" -> algebraReps
      |>
    ]
  ];

  <|
    "Tensor" -> First[tensors],
    "Keep" -> keep,
    "Dimensions" -> keptDims,
    "AlgebraReps" -> algebraReps
  |>
];

(* Component of a tensor using indices in the ORIGINAL interaction order. *)
InvariantComponent[data_Association, originalIndices_List] := Module[
  {keep, tensor, selected},
  keep = data["Keep"];
  tensor = data["Tensor"];

  If[keep === {}, Return[tensor]];

  selected = originalIndices[[keep]];
  tensor[[Sequence @@ selected]]
];

(* ------------------------------------------------------------------------- *)
(* Weinberg tensor comparison                                                *)
(* ------------------------------------------------------------------------- *)

(*
  Any overall sign convention of epsilon drops out because the Weinberg
  structure contains two epsilon tensors.
*)
eps = {{0, 1}, {-1, 0}};

CanonicalWeinbergTensor[] := Array[
  Function[{a, b, c, d},
    eps[[c, a]] eps[[d, b]]
  ],
  {2, 2, 2, 2}
];

SymmetricWeinbergTensor[] := Array[
  Function[{a, b, c, d},
    (eps[[c, a]] eps[[d, b]] + eps[[d, a]] eps[[c, b]])/2
  ],
  {2, 2, 2, 2}
];

TensorRatio[tensor_, reference_] := Module[
  {flatT, flatR, nonzero, ratios, first},

  flatT = Flatten[Normal[tensor]];
  flatR = Flatten[Normal[reference]];

  (*
    Build the support explicitly from integer positions.  Using Position with
    conditional symbolic patterns proved unreliable here and introduced one
    spurious ratio equal to 1 even when all physical nonzero components had
    the same ratio.
  *)
  nonzero = Select[
    Range[Length[flatR]],
    !TrueQ[PossibleZeroQ[flatR[[#]]]] &
  ];

  If[nonzero === {}, Return[Missing["ZeroReference"]]];

  ratios = FullSimplify /@ Table[
    flatT[[k]]/flatR[[k]],
    {k, nonzero}
  ];

  first = First[ratios];

  If[
    AllTrue[ratios, TrueQ[FullSimplify[# - first] === 0] &],
    first,
    Missing["NotProportional", ratios]
  ]
];

TensorResidual[tensor_, reference_, ratio_] := Module[{delta},
  If[MissingQ[ratio], Return[Missing["NoRatio"]]];
  delta = Simplify[Normal[tensor] - ratio Normal[reference]];
  Max[Abs[Flatten[N[delta]]]]
];

(* ------------------------------------------------------------------------- *)
(* T3 contraction                                                            *)
(* ------------------------------------------------------------------------- *)

ProbeT3Class[name_String, d1_Integer, d2_Integer, dF_Integer, expected_] := Module[
  {
    y1, y2, mix, contracted,
    canonical, symCanonical,
    ratioCanonical, ratioSym,
    residualCanonical, residualSym,
    hhFactorDiagnostic, expectedCheck
  },

  Print[""];
  Print["============================================================"];
  Print[name, "  dims = {S1=", d1, ", S2=", d2, ", F=", dF, "}"];
  Print["============================================================"];

  (* Same orientations and symmetry choices as DefineT3InvariantCGs. *)
  Print["Building y1 invariant..."];
  y1 = CheckAbort[
    RawSU2Invariant[
      {2, dF, d1},
      {True, True, False},
      {},
      {1}
    ],
    Print["ABORT while building y1 invariant."];
    Return[<|"Class" -> name, "Status" -> "Y1Aborted"|>]
  ];

  Print["Building y2 invariant..."];
  y2 = CheckAbort[
    RawSU2Invariant[
      {2, dF, d2},
      {True, False, True},
      {},
      {1}
    ],
    Print["ABORT while building y2 invariant."];
    Return[<|"Class" -> name, "Status" -> "Y2Aborted"|>]
  ];

  Print["Building scalar-mixing invariant..."];
  mix = CheckAbort[
    RawSU2Invariant[
      {2, 2, d1, d2},
      {False, False, False, True},
      {1, 2},
      {1, 2}
    ],
    Print["ABORT while building scalar-mixing invariant."];
    Return[<|"Class" -> name, "Status" -> "MixAborted"|>]
  ];

  If[
    AnyTrue[
      {y1, y2, mix},
      (# === $Failed) || !AssociationQ[#] || TrueQ[Lookup[#, "Failure", False]] &
    ],
    Print["FAILED: invariant construction did not return exactly one tensor."];
    Print["  y1 = ", InputForm[y1]];
    Print["  y2 = ", InputForm[y2]];
    Print["  mix = ", InputForm[mix]];
    Return[
      <|
        "Class" -> name,
        "Status" -> "InvariantConstructionFailed"
      |>
    ]
  ];

  Print["Y1 tensor dimensions:  ", Dimensions[Normal[y1["Tensor"]]],
    "   keep=", y1["Keep"]];
  Print["Y2 tensor dimensions:  ", Dimensions[Normal[y2["Tensor"]]],
    "   keep=", y2["Keep"]];
  Print["Mix tensor dimensions: ", Dimensions[Normal[mix["Tensor"]]],
    "   keep=", mix["Keep"]];

  (*
    IMPORTANT: the holomorphic Weinberg coefficient LLHH contains
    Conjugate[y1] Conjugate[y2] lambdaT3 in the matched result.

    Therefore the two Yukawa vertices entering this contraction are the
    Hermitian-conjugate Yukawa vertices, whereas the scalar-mixing vertex is
    the original lambdaT3 H H S1 S2^\dagger vertex.

    At the HC Yukawa vertices the internal-field orientations are:
      y1^*:  L, F,     S1^\dagger
      y2^*:  L, Fbar,  S2
    while the mixing vertex contains
      lambdaT3: H, H, S1, S2^\dagger.

    Thus every internal propagator connects a representation to its conjugate,
    so the component pairing is the natural identity pairing.  The previous
    probe incorrectly contracted the ORIGINAL Yukawa CG tensors instead of
    their Hermitian conjugates.
  *)
  contracted = Array[
    Function[{a, b, c, d},
      Sum[
        Conjugate[InvariantComponent[y1, {a, f, s1}]] *
        Conjugate[InvariantComponent[y2, {b, f, s2}]] *
        InvariantComponent[mix, {c, d, s1, s2}],
        {f, 1, dF},
        {s1, 1, d1},
        {s2, 1, d2}
      ]
    ],
    {2, 2, 2, 2}
  ];

  contracted = Simplify[contracted];

  Print["Nonzero contracted components {a,b,c,d}->value:"];
  Print[
    InputForm[
      Cases[
        Flatten[
          Table[
            {{a, b, c, d}, Simplify[contracted[[a, b, c, d]]]},
            {a, 1, 2}, {b, 1, 2}, {c, 1, 2}, {d, 1, 2}
          ],
          3
        ],
        {idx_, val_} /; !TrueQ[PossibleZeroQ[val]]
      ]
    ]
  ];

  canonical = CanonicalWeinbergTensor[];
  symCanonical = SymmetricWeinbergTensor[];

  ratioCanonical = TensorRatio[contracted, canonical];
  ratioSym = TensorRatio[contracted, symCanonical];

  residualCanonical = TensorResidual[
    contracted, canonical, ratioCanonical
  ];
  residualSym = TensorResidual[
    contracted, symCanonical, ratioSym
  ];

  Print["Pure tensor contraction:"];
  Print["  ratio to epsilon(c,a) epsilon(d,b) = ",
    InputForm[ratioCanonical],
    "   residual=", InputForm[residualCanonical]];
  Print["  ratio to H-symmetrised Weinberg tensor = ",
    InputForm[ratioSym],
    "   residual=", InputForm[residualSym]];

  (*
    Compare the contraction directly with the H-symmetrised Weinberg tensor.
    Because the reference tensor itself contains the 1/2 symmetrisation, an
    additional factor of 2 must NOT be assumed here.  We still print 2*ratio
    as a diagnostic so that any convention issue is visible.
  *)
  hhFactorDiagnostic =
    If[!MissingQ[ratioSym], Simplify[2 ratioSym], Missing["NoRatio"]];

  Print["Direct symmetric-tensor ratio: ", InputForm[ratioSym]];
  Print["Diagnostic only, 2 x ratio: ", InputForm[hhFactorDiagnostic]];
  Print["Matched-C5 prefactor expected from dimension comparison: ",
    InputForm[expected]];

  expectedCheck = If[
    MissingQ[ratioSym],
    False,
    TrueQ[FullSimplify[ratioSym - expected] === 0]
  ];

  Print["Exact direct agreement with matched C5 factor: ", expectedCheck];

  <|
    "Class" -> name,
    "Dimensions" -> {d1, d2, dF},
    "TensorRatioCanonical" -> ratioCanonical,
    "TensorRatioSymmetric" -> ratioSym,
    "TwiceTensorRatioDiagnostic" -> hhFactorDiagnostic,
    "ExpectedMatchedFactor" -> expected,
    "ExactMatch" -> expectedCheck,
    "ContractedTensor" -> contracted
  |>
];

(* ------------------------------------------------------------------------- *)
(* Original T3 A-E dimension classes                                         *)
(* ------------------------------------------------------------------------- *)

classes = {
  {"T3-A", 1, 3, 2, -2/Sqrt[3]},
  {"T3-B", 2, 2, 1, -2/Sqrt[3]},
  {"T3-C", 2, 2, 3,  2/3},
  {"T3-D", 3, 1, 2, -2/Sqrt[3]},
  {"T3-E", 3, 3, 2,  2 Sqrt[2/3]}
};

(* The same group-theory machinery used by the production Lagrangian builder
   can request front-end services, so mirror RunModel.wl and keep the complete
   invariant construction/probe inside UsingFrontEnd. *)
results = UsingFrontEnd[ProbeT3Class @@@ classes];

Print[""];
Print[""];
Print["===================== SUMMARY ====================="];
Print[
  Grid[
    Prepend[
      (
        {
          Lookup[#, "Class", "?"],
          Lookup[#, "Dimensions", Missing["Unavailable"]],
          Lookup[#, "TensorRatioSymmetric", Missing["Unavailable"]],
          Lookup[#, "TwiceTensorRatioDiagnostic", Missing["Unavailable"]],
          Lookup[#, "ExpectedMatchedFactor", Missing["Unavailable"]],
          Lookup[#, "ExactMatch", False]
        } &
      ) /@ results,
      {
        "Class",
        "{dS1,dS2,dF}",
        "Tensor ratio",
        "2 x ratio (diagnostic)",
        "Matched C5",
        "Exact?"
      }
    ],
    Frame -> All
  ]
];

If[AllTrue[results, TrueQ[#["ExactMatch"]] &],
  Print["PASS: all five A-E factors are reproduced directly by the ",
    "SU(2) contraction in the symmetrised Weinberg basis."],
  Print["NOTE: at least one class does not match exactly.  Inspect the printed ",
    "pure tensor ratios before changing any physics code; a mismatch can reveal ",
    "an orientation/conjugation or additional vertex-normalisation factor."]
];
