(*
  RunMatching.wl
  Runs the Matchete matching/simplification pipeline and returns every
  intermediate EFT expression in a diagnostic Association.
*)

ClearAll[SafeStage, RunT3Matching, RunSMBaselineMatching, BuildMatchedEFTDifference];

(* Evaluate one pipeline stage while converting messages/aborts into a single
   failure marker understood by RunT3Matching. *)
SafeStage[operation_] := CheckAbort[Check[operation, $Failed], $Aborted];

(* Subtract two EFTs produced by the same matching pipeline, then rerun the
   canonical simplification stages.  This is essential: subtracting the raw
   SM Lagrangian from a matched EFT leaves structurally equivalent SM terms
   with different generated dummy-index labels. *)
BuildMatchedEFTDifference[fullEFT_, smEFT_] := Module[
  {input, value, results = <||>, stages, stage, function, key, status},

  input = Expand[fullEFT - smEFT];
  stages = {
    {"GreensSimplifyDifference", GreensSimplify, "GreenDifference",
      "DifferenceGreensSimplifyFailed"},
    {"EOMSimplifyDifference", EOMSimplify, "EOMDifference",
      "DifferenceEOMSimplifyFailed"},
    {"EvaluateLoopFunctionsDifference", EvaluateLoopFunctions,
      "LoopDifference", "DifferenceLoopEvaluationFailed"},
    {"ReplaceEffectiveCouplingsDifference", ReplaceEffectiveCouplings,
      "BSMEFT", "DifferenceEffectiveCouplingReplacementFailed"}
  };

  Print["\nCanonicalising the matched full-minus-SM EFT difference..."];
  Do[
    {stage, function, key, status} = specification;
    Print["Running ", stage, "..."];
    value = SafeStage[function[input]];
    results[key] = value;
    If[MemberQ[{$Failed, $Aborted}, value],
      Print["ERROR: ", stage, " failed or aborted."];
      Return[Join[<|"Status" -> status|>, results]]
    ];
    input = value,
    {specification, stages}
  ];

  Print["Matched EFT difference canonicalised successfully."];
  Join[<|"Status" -> "Success"|>, results]
];


RunT3Matching[LUV_, eftOrder_Integer : 5, loopOrder_Integer : 1] := Module[
  {base, results = <||>, stages, input = LUV, value, stage, function, key, status},

  base = <|"EFTOrder" -> eftOrder, "LoopOrder" -> loopOrder|>;
  stages = {
    {"Match", Match[#, EFTOrder -> eftOrder, LoopOrder -> loopOrder] &, "RawEFT", "MatchFailed"},
    {"GreensSimplify", GreensSimplify, "GreenEFT", "GreensSimplifyFailed"},
    {"EOMSimplify", EOMSimplify, "EOMEFT", "EOMSimplifyFailed"},
    {"EvaluateLoopFunctions", EvaluateLoopFunctions, "LoopEFT", "LoopEvaluationFailed"},
    {"ReplaceEffectiveCouplings", ReplaceEffectiveCouplings, "MatchedEFT",
      "EffectiveCouplingReplacementFailed"}
  };

  Print["\nRunning Matchete matching..."];
  Print["EFT order: ", eftOrder, "; loop order: ", loopOrder];

  Do[
    {stage, function, key, status} = specification;
    Print["Running ", stage, "..."];
    value = SafeStage[function[input]];
    results[key] = value;

    (* Return all successfully produced intermediates plus the failed stage. *)
    If[MemberQ[{$Failed, $Aborted}, value],
      Print["ERROR: ", stage, " failed or aborted."];
      Return[Join[<|"Status" -> status|>, base, results]]
    ];

    input = value,
    {specification, stages}
  ];

  Print["Matching pipeline completed successfully."];
  Join[<|"Status" -> "Success"|>, base, results]
];


(* A pure SM input may contain no heavy field for Match to integrate out.
   First try the identical full pipeline.  If Match alone fails, use the same
   post-matching canonicalisation stages on LSM; this represents the matched
   SM baseline in the no-heavy-field limit. *)
RunSMBaselineMatching[LSM_, eftOrder_Integer : 5, loopOrder_Integer : 1] := Module[
  {attempt, input = LSM, value, results = <||>, stages, stage, function, key, status},

  attempt = RunT3Matching[LSM, eftOrder, loopOrder];
  If[AssociationQ[attempt] && Lookup[attempt, "Status", ""] === "Success",
    Return[Append[attempt, "BaselineMode" -> "FullMatchPipeline"]]
  ];

  If[!(AssociationQ[attempt] && Lookup[attempt, "Status", ""] === "MatchFailed"),
    Return[attempt]
  ];

  Print["Pure SM has no matchable heavy field; canonicalising LSM directly."];
  stages = {
    {"GreensSimplifySM", GreensSimplify, "GreenEFT", "SMGreensSimplifyFailed"},
    {"EOMSimplifySM", EOMSimplify, "EOMEFT", "SMEOMSimplifyFailed"},
    {"EvaluateLoopFunctionsSM", EvaluateLoopFunctions, "LoopEFT",
      "SMLoopEvaluationFailed"},
    {"ReplaceEffectiveCouplingsSM", ReplaceEffectiveCouplings, "MatchedEFT",
      "SMEffectiveCouplingReplacementFailed"}
  };

  Do[
    {stage, function, key, status} = specification;
    Print["Running ", stage, "..."];
    value = SafeStage[function[input]];
    results[key] = value;
    If[MemberQ[{$Failed, $Aborted}, value],
      Return[Join[<|"Status" -> status, "BaselineMode" -> "DirectCanonicalisation"|>, results]]
    ];
    input = value,
    {specification, stages}
  ];

  Join[<|"Status" -> "Success", "BaselineMode" -> "DirectCanonicalisation",
    "EFTOrder" -> eftOrder, "LoopOrder" -> loopOrder|>, results]
];
