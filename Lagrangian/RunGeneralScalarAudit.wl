(* Experimental Matchete validation of symmetry-allowed T3 scalar quartics.

   wolframscript -file Lagrangian/RunGeneralScalarAudit.wl OUTPUT.json B m1

   IMPORTANT: Matchete's CheckLagrangian requires a usable Wolfram front end.
   Unlike the first trial, this runner keeps the entire model construction and
   validation inside UsingFrontEnd, not merely Needs["Matchete`"]. The result
   distinguishes an unavailable front end from an invalid interaction.

   This does not perform matching, RGBeta running, or numerical integration.
*)

ClearAll["Global`*"];
trialDir = DirectoryName @ ExpandFileName[$InputFileName];

T3TrialParseSigned[s_String] := Which[
  StringMatchQ[s, "m" ~~ DigitCharacter ..], -ToExpression[StringDrop[s, 1]],
  StringMatchQ[s, "p" ~~ DigitCharacter ..], ToExpression[StringDrop[s, 1]],
  StringMatchQ[s, ("-" | "") ~~ DigitCharacter ..], ToExpression[s],
  True, $Failed
];

trialArgs = Rest[$ScriptCommandLine];
If[Length[trialArgs] =!= 3,
  Print["Usage: RunGeneralScalarAudit.wl OUTPUT.json CLASS alpha"];
  Exit[2]
];

trialOutput = ExpandFileName[trialArgs[[1]]];
trialClass = ToUpperCase[trialArgs[[2]]];
trialAlpha = T3TrialParseSigned[trialArgs[[3]]];
If[trialAlpha === $Failed, Print["Invalid alpha"]; Exit[2]];

(* Load the symbolic functions once; the service front end used for this load
   is NOT persistent when UsingFrontEnd exits. A second UsingFrontEnd must
   enclose the full audit below. *)
trialLoadStatus = CheckAbort[
  Check[UsingFrontEnd[Needs["Matchete`"]; True], $Failed],
  $Aborted
];
If[!TrueQ[trialLoadStatus],
  Print["ERROR: Matchete could not be loaded."];
  Exit[3]
];

Get[FileNameJoin[{trialDir, "model", "T3ModelCatalog.wl"}]];
Get[FileNameJoin[{trialDir, "model", "LagrangianBuilder.wl"}]];

trialModel = T3ModelFromClass[trialClass, trialAlpha];
If[trialModel === $Failed,
  Print["Invalid T3 model class: ", trialClass];
  Exit[2]
];

T3TrialAuditRun[model_Association] := Module[
  {feReady, lsm, lFree, extra, checked, accepted, rejected, baseline,
   baselineLag, combinedLag, baselineValid, combinedValid, massCatalogue,
   status, candidateNames, checkedNames},

  feReady = MatchQ[$FrontEnd, _FrontEndObject];
  Print["Front-end connected during validation: ", feReady];

  If[!TrueQ[feReady],
    Return[<|
      "Status" -> "FrontendUnavailable",
      "FrontendConnected" -> False,
      "ValidationPerformed" -> False,
      "Diagnostic" -> "UsingFrontEnd could not attach a front end. No operator was declared invalid.",
      "ExtraQuarticCandidateCount" -> Null,
      "AcceptedExtraQuarticCount" -> Null,
      "CombinedLagrangianValid" -> Null
    |>]
  ];

  ResetAll[];
  lsm = LoadModel["SM"];
  If[lsm === $Failed, Return[<|"Status" -> "SMLoadFailed"|>]];
  If[DefineT3Fields[model] === $Failed,
    Return[<|"Status" -> "FieldDefinitionFailed"|>]
  ];
  If[DefineT3Couplings[] === $Failed,
    Return[<|"Status" -> "CouplingDefinitionFailed"|>]
  ];
  If[DefineT3InvariantCGs[model] === $Failed,
    Return[<|"Status" -> "TopologyCGDefinitionFailed"|>]
  ];

  lFree = FreeLag[NewFermion, NewScalar1, NewScalar2] // RelabelIndices;
  extra = CheckAbort[T3ExtendedBuildExtraQuarticCandidates[model], $Aborted];
  If[!ListQ[extra],
    Return[<|
      "Status" -> "ExtraQuarticGenerationFailed",
      "GeneratorResult" -> ToString[extra, InputForm]
    |>]
  ];
  candidateNames = Lookup[extra, "Name", {}];
  Print["Verified extra quartic target: ",
    T3ExtendedVerifiedStructureCount[model]];
  Print["Generated extra quartic count: ", Length[extra]];

  (* First verify the ordinary, unextended model in precisely the same
     front-end context as the extended one. This distinguishes a regression
     in the extension from an already-failing baseline. *)
  baseline = Join[
    BuildT3YukawaCandidates[model, 1],
    BuildT3YukawaCandidates[model, 2],
    BuildGeneralScalarPotentialCandidates[model],
    BuildT3MixingCandidates[model]
  ];
  baselineLag = Expand[
    lsm + lFree + Total[Lookup[baseline, "Expression", {}]]
  ] // RelabelIndices;
  baselineValid = CheckAbort[
    Check[CheckLagrangian[baselineLag], $Failed], $Aborted
  ];
  Print["Baseline full Lagrangian valid: ", InputForm[baselineValid]];

  (* A failed baseline is a blocker. There is no point attributing failures
     of extended interactions to SU(2) contractions until it is resolved. *)
  If[!TrueQ[baselineValid],
    Return[<|
      "Status" -> "BaselineValidationFailed",
      "FrontendConnected" -> True,
      "ValidationPerformed" -> True,
      "ExtraQuarticCandidateCount" -> Length[extra],
      "ExtraQuarticCandidateNames" -> candidateNames,
      "AcceptedExtraQuarticCount" -> Null,
      "BaselineLagrangianValid" -> False,
      "CombinedLagrangianValid" -> Null,
      "Diagnostic" -> "The original T3 Lagrangian failed validation before the extension was tested."
    |>]
  ];

  checked = ValidateT3Candidate[#, lsm, lFree] & /@ extra;
  accepted = Select[checked, TrueQ[#["Valid"]] &];
  rejected = Select[checked, !TrueQ[#["Valid"]] &];
  checkedNames = Lookup[accepted, "Name", {}];

  (* Combined validation retains all baseline Yukawas and quartics and ALL
     generated extras, rather than silently discarding rejected candidates. *)
  combinedLag = Expand[
    baselineLag + Total[Lookup[extra, "Expression", {}]]
  ] // RelabelIndices;
  combinedValid = CheckAbort[
    Check[CheckLagrangian[combinedLag], $Failed], $Aborted
  ];

  status = Which[
    Length[rejected] > 0, "IndividualOperatorValidationFailed",
    !TrueQ[combinedValid], "CombinedLagrangianValidationFailed",
    True, "Success"
  ];

  massCatalogue = T3ExtendedDescribeQuadraticOperators[model];

  <|
    "Status" -> status,
    "FrontendConnected" -> True,
    "ValidationPerformed" -> True,
    "ModelClass" -> model["Class"],
    "Alpha" -> model["Alpha"],
    "ExtraQuarticCandidateCount" -> Length[extra],
    "ExpectedExtraQuarticCount" -> T3ExtendedVerifiedStructureCount[model],
    "ExtraQuarticCandidateNames" -> candidateNames,
    "AcceptedExtraQuarticCount" -> Length[accepted],
    "AcceptedExtraQuartics" -> checkedNames,
    "RejectedExtraQuartics" -> Lookup[rejected, "Name", {}],
    "IndividualValidationResults" -> Map[
      <|"Name" -> #["Name"], "Valid" -> TrueQ[#["Valid"]],
        "Result" -> ToString[#["ValidationResult"], InputForm]|> &,
      checked
    ],
    "ExtraQuadraticsCataloguedNotIncluded" -> massCatalogue,
    "BaselineLagrangianValid" -> True,
    "CombinedLagrangianValid" -> TrueQ[combinedValid],
    "MatchingExecuted" -> False,
    "RGBetaExtended" -> False,
    "NumericsExtended" -> False
  |>
];

trialResult = CheckAbort[
  UsingFrontEnd[T3TrialAuditRun[trialModel]],
  <|"Status" -> "Aborted"|>
];
If[!AssociationQ[trialResult],
  trialResult = <|"Status" -> "UnexpectedRunnerFailure",
    "Result" -> ToString[trialResult, InputForm]|>
];
trialResult = Join[
  <|"RequestedClass" -> trialClass, "RequestedAlpha" -> trialAlpha|>,
  trialResult
];

If[!DirectoryQ[DirectoryName[trialOutput]],
  CreateDirectory[DirectoryName[trialOutput], CreateIntermediateDirectories -> True]
];
Export[trialOutput, trialResult, "RawJSON"];
Print["General scalar audit written to ", trialOutput];
Print["Status: ", trialResult["Status"]];
Print["Extra quartic candidates: ", Lookup[trialResult, "ExtraQuarticCandidateCount", Null]];
Print["Accepted: ", Lookup[trialResult, "AcceptedExtraQuarticCount", Null]];
Print["Combined Lagrangian valid: ", Lookup[trialResult, "CombinedLagrangianValid", Null]];
If[trialResult["Status"] =!= "Success", Exit[6]];
Exit[0];
