(* Construct ordinary T3 UV Lagrangians from validated independent SU(2)
   invariants.  An opt-in general-scalar build never silently discards an
   allowed holomorphic quadratic.  B/C(-1) off-diagonal scalar masses cannot
   pass Matchete's HeavyMassBasis test without diagonalization.  Such a
   Lagrangian is exported as a blocked symbolic object, NOT passed to Match. *)
ClearAll[
  ValidateT3Candidate,
  SelectFirstValidByGroup,
  T3IngredientsPresentQ,
  BuildT3Lagrangian
];

Get[FileNameJoin[{DirectoryName @ ExpandFileName[$InputFileName],
  "SU2Invariants.wl"}]];
Get[FileNameJoin[{DirectoryName @ ExpandFileName[$InputFileName],
  "..", "interactions", "ScalarPotential.wl"}]];
Get[FileNameJoin[{DirectoryName @ ExpandFileName[$InputFileName],
  "T3Fields.wl"}]];
Get[FileNameJoin[{DirectoryName @ ExpandFileName[$InputFileName],
  "..", "interactions", "T3Topology.wl"}]];

ValidateT3Candidate[c_, LSM_, LFree_] := Module[{res},
  Print["Checking candidate: ", c["Name"]];
  res = CheckAbort[
    Check[CheckLagrangian[(LSM + LFree + c["Expression"]) // RelabelIndices],
      $Failed],
    $Aborted
  ];
  Print["  result: ", InputForm[res]];
  Association[c, "Valid" -> TrueQ[res], "ValidationResult" -> res]
];

SelectFirstValidByGroup[list_List] := Module[{seen = <||>},
  Select[list, Function[candidate,
    Module[{group = Lookup[candidate, "AlternativeGroup", None]},
      If[group === None, True,
        If[KeyExistsQ[seen, group], False, seen[group] = True; True]]
    ]
  ]]
];

T3IngredientsPresentQ[valid_List] := And @@ (
  Function[class, AnyTrue[valid, Lookup[#, "Class", ""] === class &]] /@
    {"Yukawa1", "Yukawa2", "T3ScalarMix"}
);

BuildT3Lagrangian[model_Association] := Module[
  {LSM, LFree, cgStatus, candidates, validated, valid, rejected,
   LInt, LBSM, LUV, fullValidation, ingredientsPresent,
   massRequested, massExpected, massCandidates, massAccepted,
   rejectedRequired, mixedHeavy, massBasisBlocked, massIncluded,
   assembled, status},

  ResetAll[];
  Print["Build stage 1/4: load SM"];
  LSM = LoadModel["SM"];
  If[LSM === $Failed, Return[$Failed]];

  Print["Build stage 2/4: define BSM fields"];
  If[DefineT3Fields[model] === $Failed, Return[$Failed]];

  Print["Build stage 3/4: define couplings"];
  If[Check[DefineT3Couplings[], $Failed] === $Failed, Return[$Failed]];

  Print["Build stage 4/4: define SU(2) CG tensors"];
  cgStatus = DefineT3InvariantCGs[model];
  If[cgStatus === $Failed, Return[$Failed]];

  LFree = FreeLag[NewFermion, NewScalar1, NewScalar2] // RelabelIndices;
  candidates = BuildT3InteractionCandidates[model];
  If[!ListQ[candidates] || !AllTrue[candidates, AssociationQ],
    Print["ERROR: candidate construction failed; no incomplete Lagrangian returned."];
    Return[$Failed]
  ];
  candidates = Association[#, "Expression" -> RelabelIndices[#["Expression"]]] & /@
    candidates;
  validated = ValidateT3Candidate[#, LSM, LFree] & /@ candidates;

  massRequested = TrueQ[T3ExperimentalGeneralQuadraticsQ[]];
  massExpected = If[massRequested,
    Length[T3ExtendedExtraQuadraticStructures[model]], 0];
  massCandidates = Select[validated,
    Lookup[#, "Class", ""] === "ExtraQuadratic" &];
  If[Length[massCandidates] =!= massExpected,
    Print["ERROR: generated extra-bilinear count does not match operator inventory."];
    Return[$Failed]
  ];

  (* B/C(-1) has a non-diagonal charged scalar mass term.  Its operator is
     kept in the symbolic expression, but ordinary Matchete validation
     cannot authorize matching.  Other rejections remain hard errors. *)
  mixedHeavy = massRequested &&
    model["Alpha"] === -1 &&
    model["Scalar1", "SU2"] === 2 &&
    model["Scalar2", "SU2"] === 2 &&
    MemberQ[{1, 3}, model["Fermion", "SU2"]];
  massBasisBlocked = Select[massCandidates, !TrueQ[#["Valid"]] &];
  massAccepted = Select[massCandidates, TrueQ[#["Valid"]] &];
  rejectedRequired = Select[validated,
    !TrueQ[#["Valid"]] &&
    (Lookup[#, "Class", ""] === "ExtraQuartic" ||
      (Lookup[#, "Class", ""] === "ExtraQuadratic" && !mixedHeavy)) &];

  valid = SelectFirstValidByGroup @ Select[validated, TrueQ[#["Valid"]] &];
  (* Include a blocked B/C mixed mass in the symbolic LUV; never drop it. *)
  assembled = If[mixedHeavy, Join[valid, massBasisBlocked], valid];
  massIncluded = Count[assembled,
    _?(Lookup[#, "Class", ""] === "ExtraQuadratic" &)];
  rejected = Select[validated, !TrueQ[#["Valid"]] &];

  LInt = Total[Lookup[assembled, "Expression", {}]] // Expand // RelabelIndices;
  LBSM = (LFree + LInt) // Expand // RelabelIndices;
  LUV = (LSM + LBSM) // Expand // RelabelIndices;
  fullValidation = CheckAbort[Check[CheckLagrangian[LUV], $Failed], $Aborted];
  ingredientsPresent = T3IngredientsPresentQ[valid];

  status = Which[
    Length[rejectedRequired] > 0 || massIncluded =!= massExpected,
      "RequiredGeneralScalarOperatorRejected",
    Length[massBasisBlocked] > 0 && mixedHeavy,
      "RequiresMixedScalarMassDiagonalization",
    TrueQ[fullValidation], "Success",
    True, "FullValidationFailed"
  ];
  If[status =!= "Success",
    Print["Lagrangian symbolic build status: ", status,
      "; Match[] must NOT be called on this result."]
  ];
  <|
    "Status" -> status,
    "Model" -> model,
    "LSM" -> LSM,
    "LFree" -> LFree,
    "LInt" -> LInt,
    "LBSM" -> LBSM,
    "LUV" -> LUV,
    "AllowedInteractions" -> Lookup[assembled, "Name", {}],
    "RejectedInteractions" -> Lookup[rejected, "Name", {}],
    "T3IngredientsPresent" -> ingredientsPresent,
    "WeinbergIngredientsPresent" -> ingredientsPresent,
    "AdditionalQuadraticsRequested" -> massRequested,
    "AdditionalQuadraticsExpected" -> massExpected,
    "AdditionalQuadraticsIncluded" -> massIncluded,
    "AdditionalQuadraticsIndividuallyValidated" -> Length[massAccepted],
    "AdditionalQuadraticsHeavyMassBasisBlocked" -> Length[massBasisBlocked],
    "MatchingMassBasisReady" -> (status === "Success"),
    "FullValidation" -> fullValidation
  |>
];
