(*
  RunMatching.wl
  Matchete UV -> EFT matching plus Weinberg-operator extraction.
  Intermediate EFTs are retained for diagnostics.
*)

ClearAll[
  SafeStage,
  RunStages,
  MatchingStages,
  RunT3Matching,
  RunSMBaselineMatching,
  BuildMatchedEFTDifference
];


(* ---------------------------------------------------------------------- *)
(* Matching pipeline                                                       *)
(* ---------------------------------------------------------------------- *)

(* Wrapper to make the code deal with errors *)
SafeStage[operation_] := CheckAbort[Check[operation, $Failed], $Aborted];

(* Run ordered matching stages while retaining every intermediate result. *)
(* We are running many functions in a row *)
RunStages[input_, stages_List, prefix_Association : <||>] := Module[
  {current = input, results = <||>, stageName, function, key, failureStatus, value},

  Do[
    {stageName, function, key, failureStatus} = specification;
    Print["Running ", stageName, "..."];

    (* Specification looks like {  humanReadableName,  functionToRun,  resultKey,  errorStatus}
    Example is {"Match",  Match[#, EFTOrder -> eftOrder, LoopOrder -> loopOrder] &, "RawEFT",  "MatchFailed"}*)
    value = SafeStage[function[current]];
    results[key] = value; (* We keep intermediate representations, not just leaving the last one *)

    (* Detect failures *)
    If[MemberQ[{$Failed, $Aborted}, value],
      Print["ERROR: ", stageName, " failed or aborted."];
      Return[Join[<|"Status" -> failureStatus|>, prefix, results]]
    ];

    (* Output to next stage *)
    current = value,
    {specification, stages}
  ];

  Join[<|"Status" -> "Success"|>, prefix, results]
];

(* We get actual Matchete workflow 

First we match -> then we green simplify to get green basis, getting rid of redundancies coming from integration by parts and identities etc ->
EOM simplify remove operator redundancies from EOM -> Evaluate loop functions -> replace effective couplings with our original couplings*)
MatchingStages[eftOrder_Integer, loopOrder_Integer] := {
  {
    "Match",
    Match[#, EFTOrder -> eftOrder, LoopOrder -> loopOrder] &,
    "RawEFT",
    "MatchFailed"
  },
  {"GreensSimplify", GreensSimplify, "GreenEFT", "GreensSimplifyFailed"},
  {"EOMSimplify", EOMSimplify, "EOMEFT", "EOMSimplifyFailed"},
  {
    "EvaluateLoopFunctions",
    EvaluateLoopFunctions,
    "LoopEFT",
    "LoopEvaluationFailed"
  },
  {
    "ReplaceEffectiveCouplings",
    ReplaceEffectiveCouplings,
    "MatchedEFT",
    "EffectiveCouplingReplacementFailed"
  }
};

(* Convenient Wrapper to get EFT order and loop order *)
RunT3Matching[LUV_, eftOrder_Integer : 5, loopOrder_Integer : 1] := Module[
  {metadata},

  metadata = <|"EFTOrder" -> eftOrder, "LoopOrder" -> loopOrder|>;

  Print["\nRunning Matchete matching..."];
  Print["EFT order: ", eftOrder, "; loop order: ", loopOrder];

  With[{result = RunStages[LUV, MatchingStages[eftOrder, loopOrder], metadata]},
    If[Lookup[result, "Status", ""] === "Success",
      Print["Matching pipeline completed successfully."]
    ];
    result
  ]
];

(* Pure SM has no heavy field: if Match fails, canonicalise LSM directly. *)
(* We are comparing SM and full lagrangian *)
RunSMBaselineMatching[LSM_, eftOrder_Integer : 5, loopOrder_Integer : 1] := Module[
  {attempt, metadata, stages, result},

  attempt = RunT3Matching[LSM, eftOrder, loopOrder];

  If[AssociationQ[attempt] && Lookup[attempt, "Status", ""] === "Success",
    Return[Append[attempt, "BaselineMode" -> "FullMatchPipeline"]]
  ];

  If[!(AssociationQ[attempt] && Lookup[attempt, "Status", ""] === "MatchFailed"),
    Return[attempt]
  ];

  Print["Pure SM has no matchable heavy field; canonicalising LSM directly."];

  (* Keep the historical failure-status strings used by downstream diagnostics. *)
  stages = {
    {"GreensSimplifySM", GreensSimplify, "GreenEFT", "SMGreensSimplifyFailed"},
    {"EOMSimplifySM", EOMSimplify, "EOMEFT", "SMEOMSimplifyFailed"},
    {
      "EvaluateLoopFunctionsSM",
      EvaluateLoopFunctions,
      "LoopEFT",
      "SMLoopEvaluationFailed"
    },
    {
      "ReplaceEffectiveCouplingsSM",
      ReplaceEffectiveCouplings,
      "MatchedEFT",
      "SMEffectiveCouplingReplacementFailed"
    }
  };

  metadata = <|
    "BaselineMode" -> "DirectCanonicalisation",
    "EFTOrder" -> eftOrder,
    "LoopOrder" -> loopOrder
  |>;

  RunStages[LSM, stages, metadata]
];

(* Gives us BSM contribution to EFT *)
BuildMatchedEFTDifference[fullEFT_, smEFT_] := Module[
  {
    input,
    greenDifference,
    eomDifference,
    canonicalInput,
    loopDifference,
    bsmEFT
  },

  input = Expand[fullEFT - smEFT];

  Print["\nCanonicalising the matched full-minus-SM EFT difference..."];

  Print["Running GreensSimplifyDifference..."];
  greenDifference = SafeStage[GreensSimplify[input]];
  If[MemberQ[{$Failed, $Aborted}, greenDifference],
    Print["ERROR: GreensSimplifyDifference failed or aborted."];
    Return[<|
      "Status" -> "DifferenceGreensSimplifyFailed",
      "InputDifference" -> input,
      "GreenDifference" -> greenDifference
    |>]
  ];

  (* A difference-only EFT need not contain the complete kinetic structure
     expected by EOMSimplify. If it rejects the difference, keep the
     GreensSimplify result and continue with the remaining canonicalisation. *)
  Print["Running EOMSimplifyDifference..."];
  eomDifference = SafeStage[EOMSimplify[greenDifference]];

  canonicalInput = If[
    MemberQ[{$Failed, $Aborted}, eomDifference],
    Print[
      "EOMSimplifyDifference could not canonicalise the difference-only EFT; ",
      "continuing from GreenDifference."
    ];
    greenDifference,
    eomDifference
  ];

  Print["Running EvaluateLoopFunctionsDifference..."];
  loopDifference = SafeStage[EvaluateLoopFunctions[canonicalInput]];
  If[MemberQ[{$Failed, $Aborted}, loopDifference],
    Print["ERROR: EvaluateLoopFunctionsDifference failed or aborted."];
    Return[<|
      "Status" -> "DifferenceLoopEvaluationFailed",
      "InputDifference" -> input,
      "GreenDifference" -> greenDifference,
      "EOMDifference" -> eomDifference,
      "LoopDifference" -> loopDifference
    |>]
  ];

  Print["Running ReplaceEffectiveCouplingsDifference..."];
  bsmEFT = SafeStage[ReplaceEffectiveCouplings[loopDifference]];
  If[MemberQ[{$Failed, $Aborted}, bsmEFT],
    Print["ERROR: ReplaceEffectiveCouplingsDifference failed or aborted."];
    Return[<|
      "Status" -> "DifferenceEffectiveCouplingReplacementFailed",
      "InputDifference" -> input,
      "GreenDifference" -> greenDifference,
      "EOMDifference" -> eomDifference,
      "LoopDifference" -> loopDifference,
      "BSMEFT" -> bsmEFT
    |>]
  ];

  Print["Matched EFT difference canonicalised successfully."];

  <|
    "Status" -> "Success",
    "InputDifference" -> input,
    "GreenDifference" -> greenDifference,
    "EOMDifference" -> eomDifference,
    "EOMFallbackUsed" -> MemberQ[{$Failed, $Aborted}, eomDifference],
    "LoopDifference" -> loopDifference,
    "BSMEFT" -> bsmEFT
  |>
];


(* ---------------------------------------------------------------------- *)
(* Weinberg-operator extraction                                            *)
(* ---------------------------------------------------------------------- *)

ClearAll[
  InternalHeadName,
  InternalSymbolName,
  MatcheteFieldName,
  CountMatcheteField,
  ContainsNamedSymbolQ,
  ContainsProjectorQ,
  BarredMatcheteFieldQ,
  WeinbergLHHTermQ,
  ExtractWeinbergTerms,
  StripWeinbergOperatorStructure,
  CompactWeinbergCoefficient,
  ExtractWeinbergCoefficient
];

(* Use symbol names rather than contexts to tolerate Matchete context changes. *)
InternalHeadName[x_] := Quiet@Check[
  SymbolName[Unevaluated[Head[x]]],
  ToString[Unevaluated[Head[x]], InputForm]
];

InternalSymbolName[x_Symbol] := SymbolName[Unevaluated[x]];
InternalSymbolName[x_] := ToString[Unevaluated[x], InputForm];

MatcheteFieldName[field_] := Module[{args},
  If[InternalHeadName[Unevaluated[field]] =!= "Field", Return[""]];
  args = List @@ Unevaluated[field];
  If[args === {}, "", InternalSymbolName[args[[1]]]]
];

CountMatcheteField[expr_, name_String] := Count[
  Unevaluated[expr],
  object_ /; InternalHeadName[Unevaluated[object]] === "Field" &&
    MatcheteFieldName[Unevaluated[object]] === name,
  Infinity
];

ContainsNamedSymbolQ[expr_, name_String] := !FreeQ[
  Unevaluated[expr],
  symbol_Symbol /; InternalSymbolName[Unevaluated[symbol]] === name,
  Infinity
];

ContainsProjectorQ[expr_, chirality_Integer] := !FreeQ[
  Unevaluated[expr],
  object_ /; InternalHeadName[Unevaluated[object]] === "Proj" &&
    Quiet@Check[(List @@ Unevaluated[object]) === {chirality}, False],
  Infinity
];

BarredMatcheteFieldQ[object_, name_String] := Module[{args, inner},
  If[InternalHeadName[Unevaluated[object]] =!= "Bar", Return[False]];
  args = List @@ Unevaluated[object];
  If[Length[args] =!= 1, Return[False]];

  inner = args[[1]];
  InternalHeadName[Unevaluated[inner]] === "Field" &&
    MatcheteFieldName[Unevaluated[inner]] === name
];

(* GammaCC is the validated presence marker for the Weinberg operator. *)
WeinbergLHHTermQ[term_] := !FreeQ[Unevaluated[term], GammaCC];

ExtractWeinbergTerms[eft_] := Module[{expanded, terms},
  expanded = Expand[eft];

  (* Isolate multiplicative terms containing the GammaCC spinor chain. *)
  terms = DeleteDuplicates @ Cases[
    expanded,
    term_Times /; !FreeQ[term, GammaCC],
    Infinity
  ];

  (* Defensive fallback for a coefficient-free bare spinor chain. *)
  If[terms === {} && !FreeQ[expanded, GammaCC],
    terms = DeleteDuplicates @ Cases[
      expanded,
      chain_NCM /; !FreeQ[chain, GammaCC],
      Infinity
    ]
  ];

  terms
];

(* Strip only the universal LLHH structure; retain the physical prefactor. *)
StripWeinbergOperatorStructure[term_] := Module[{stripped},
  stripped = term /. {
    chain_NCM /; !FreeQ[chain, GammaCC] :> 1,
    Field[H, Scalar, inds_, derivs_] :> 1,
    CG[___] :> 1
  };

  Quiet@Check[Simplify[Expand[stripped]], stripped]
];

(* Combine holomorphic prefactors into a compact C5 expression. *)
CompactWeinbergCoefficient[terms_List] := Module[{pieces, combined},
  If[terms === {}, Return[Missing["NoHolomorphicTerms"]]];

  pieces = StripWeinbergOperatorStructure /@ terms;
  combined = Total[pieces];

  Quiet@Check[
    FactorTerms[Cancel[Together[combined]]],
    Simplify[combined]
  ]
];

ExtractWeinbergCoefficient[eft_] := Module[
  {
    present,
    terms,
    holomorphicTerms,
    conjugateTerms,
    sector,
    holomorphicSector,
    coefficient,
    status
  },

  (* Presence and detailed isolation are separate to avoid false negatives. *)
  present = !FreeQ[eft, GammaCC];

  If[!TrueQ[present],
    Return[<|
      "Status" -> "NotFound",
      "Present" -> False,
      "TermCount" -> 0,
      "Terms" -> {},
      "Sector" -> 0,
      "Coefficient" -> 0
    |>]
  ];

  terms = ExtractWeinbergTerms[eft];

  If[terms === {},
    Return[<|
      "Status" -> "PresentIsolationPending",
      "Present" -> True,
      "TermCount" -> 0,
      "Terms" -> {},
      "Sector" -> 0,
      "Coefficient" -> Missing["PendingCanonicalisation"]
    |>]
  ];

  sector = Simplify[Expand[Total[terms]]];

  (* C5 uses the P_L LLHH sector; keep P_R only as the Hermitian-conjugate audit. *)
  holomorphicTerms = Select[terms, !FreeQ[#, Proj[-1]] &];
  conjugateTerms = Select[terms, !FreeQ[#, Proj[1]] &];

  If[holomorphicTerms === {},
    Return[<|
      "Status" -> "PresentHolomorphicIsolationPending",
      "Present" -> True,
      "TermCount" -> Length[terms],
      "HolomorphicTermCount" -> 0,
      "ConjugateTermCount" -> Length[conjugateTerms],
      "Terms" -> terms,
      "HolomorphicTerms" -> {},
      "ConjugateTerms" -> conjugateTerms,
      "Sector" -> sector,
      "HolomorphicSector" -> 0,
      "Coefficient" -> Missing["PendingCanonicalisation"]
    |>]
  ];

  holomorphicSector = Simplify[Expand[Total[holomorphicTerms]]];
  coefficient = CompactWeinbergCoefficient[holomorphicTerms];
  status = If[MissingQ[coefficient], "CoefficientPending", "Success"];

  <|
    "Status" -> status,
    "Present" -> True,
    "TermCount" -> Length[terms],
    "HolomorphicTermCount" -> Length[holomorphicTerms],
    "ConjugateTermCount" -> Length[conjugateTerms],
    "Terms" -> terms,
    "HolomorphicTerms" -> holomorphicTerms,
    "ConjugateTerms" -> conjugateTerms,
    "Sector" -> sector,
    "HolomorphicSector" -> holomorphicSector,
    "Coefficient" -> coefficient
  |>
];
