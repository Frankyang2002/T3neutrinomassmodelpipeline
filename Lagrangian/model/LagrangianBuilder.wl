(* Construct ordinary T3 UV Lagrangians from validated SU(2) invariants.
   Scalar quartics/quadratics and additional fermionic Yukawas are opt-in.
   Do not interpret a successful Matchete CheckLagrangian as proof of
   renormalizable operator-basis completeness. The one-symbol MF flavour
   mass parameter and mixed scalar mass basis remain separate limitations.

   B/C(-1) off-diagonal scalar masses cannot pass Matchete's HeavyMassBasis
   test without diagonalization. Such a Lagrangian is exported as a blocked
   symbolic object, NOT passed to Match. *)
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
   assembled, status, extraYRequested, extraYExpected, extraYCandidates,
   extraYAccepted, extraYIncluded, extraYAllValid, quarticExpected,
   quarticIncluded, scalarAllRequested},

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

  extraYRequested = TrueQ[T3ExperimentalGeneralFermionsQ[]];
  extraYExpected = T3GeneralExtraYukawaCount[model];
  extraYCandidates = Select[validated,
    Lookup[#, "Class", ""] === "ExtraYukawa" &];
  extraYAccepted = Select[extraYCandidates, TrueQ[#["Valid"]] &];
  If[extraYRequested && Length[extraYCandidates] =!= extraYExpected,
    Print["ERROR: generated extra-Yukawa count does not match gauge inventory."];
    Return[$Failed]
  ];
  extraYAllValid = extraYRequested &&
    Length[extraYAccepted] === extraYExpected;

  massRequested = TrueQ[T3ExperimentalGeneralQuadraticsQ[]];
  massExpected = If[massRequested,
    Length[T3ExtendedExtraQuadraticStructures[model]], 0];
  massCandidates = Select[validated,
    Lookup[#, "Class", ""] === "ExtraQuadratic" &];
  If[Length[massCandidates] =!= massExpected,
    Print["ERROR: generated extra-bilinear count does not match operator inventory."];
    Return[$Failed]
  ];

  (* B/C(-1) has a non-diagonal charged scalar mass term. Its operator is
     kept in the symbolic expression, but ordinary Matchete validation
     cannot authorize matching. Other rejections remain hard errors. *)
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
      (Lookup[#, "Class", ""] === "ExtraQuadratic" && !mixedHeavy) ||
      (Lookup[#, "Class", ""] === "ExtraYukawa")) &];

  valid = SelectFirstValidByGroup @ Select[validated, TrueQ[#["Valid"]] &];
  (* Include a blocked B/C mixed mass in symbolic LUV; never drop it. *)
  assembled = If[mixedHeavy, Join[valid, massBasisBlocked], valid];
  massIncluded = Count[assembled,
    _?(Lookup[#, "Class", ""] === "ExtraQuadratic" &)];
  extraYIncluded = Count[assembled,
    _?(Lookup[#, "Class", ""] === "ExtraYukawa" &)];
  rejected = Select[validated, !TrueQ[#["Valid"]] &];

  LInt = Total[Lookup[assembled, "Expression", {}]] // Expand // RelabelIndices;
  LBSM = (LFree + LInt) // Expand // RelabelIndices;
  LUV = (LSM + LBSM) // Expand // RelabelIndices;
  fullValidation = CheckAbort[Check[CheckLagrangian[LUV], $Failed], $Aborted];
  ingredientsPresent = T3IngredientsPresentQ[valid];

  (* Structure inventory is not identical to a physical-basis certificate. *)
  quarticExpected = T3ExtendedVerifiedStructureCount[model];
  quarticIncluded = Count[assembled,
    _?(Lookup[#, "Class", ""] === "ExtraQuartic" &)];
  scalarAllRequested = TrueQ[T3ExperimentalGeneralQuadraticsQ[]] &&
    IntegerQ[quarticExpected] && quarticIncluded === quarticExpected;

  status = Which[
    !ingredientsPresent, "RequiredTopologyOperatorRejected",
    Length[rejectedRequired] > 0 &&
      AnyTrue[rejectedRequired, Lookup[#, "Class", ""] === "ExtraYukawa" &],
      "RequiredGeneralFermionOperatorRejected",
    Length[rejectedRequired] > 0 || massIncluded =!= massExpected,
      "RequiredGeneralScalarOperatorRejected",
    extraYRequested && extraYIncluded =!= extraYExpected,
      "RequiredGeneralFermionOperatorRejected",
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
    "AdditionalYukawasExpected" -> extraYExpected,
    "AdditionalYukawasRequested" -> extraYRequested,
    "AdditionalYukawasIncluded" -> extraYIncluded,
    "AdditionalYukawasIndividuallyValidated" -> Length[extraYAccepted],
    "FermionYukawaInventoryAssembled" -> (extraYExpected === 0 ||
      TrueQ[extraYAllValid] && extraYIncluded === extraYExpected),
    "AdditionalQuadraticsRequested" -> massRequested,
    "AdditionalQuadraticsExpected" -> massExpected,
    "AdditionalQuadraticsIncluded" -> massIncluded,
    "AdditionalQuadraticsIndividuallyValidated" -> Length[massAccepted],
    "AdditionalQuadraticsHeavyMassBasisBlocked" -> Length[massBasisBlocked],
    "AdditionalQuarticsExpected" -> quarticExpected,
    "AdditionalQuarticsIncluded" -> quarticIncluded,
    "AllScalarExtensionsRequested" -> scalarAllRequested,
    "MatchingMassBasisReady" -> (status === "Success"),
    "FullValidation" -> fullValidation,
    (* A mass FLAVOR certificate and complete operator proof are not yet
       available. Do not promote the inventory result to certification. *)
    "FermionFlavorMassGeneralityCertified" -> False,
    "FullSymmetryCompletenessCertified" -> False
  |>
];
