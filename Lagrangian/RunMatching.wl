(*
  RunMatching.wl
  Matchete UV -> EFT matching plus Weinberg-operator extraction.
  Intermediate EFTs are retained for diagnostics.

  All this file does is 
  1. We match each stage
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

(* Run ordered matching stages while retaining every intermediate result.
 We are running many functions in a row, this a helper function for function running*)
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
1. First we match 
2. then we green simplify to get green basis, getting rid of redundancies coming from integration by parts and identities etc 
3. EOM simplify remove operator redundancies from EOM 
4. Evaluate loop functions 
5. replace effective couplings with our original couplings*)
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

(* Convenient Wrapper to get EFT order and loop order, just does the matching *)
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

(* Pure SM has no heavy field: if Match fails, canonicalise LSM directly and use standard forms for it. 
 We are comparing SM and full lagrangian so we match Pure SM as well*)
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

(* Gives us BSM contribution to EFT by doing FullEFT-SMEFT=BSMEFT *)
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

(* GammaCC is the marker for the Weinberg operator. *)
WeinbergLHHTermQ[term_] := !FreeQ[Unevaluated[term], GammaCC];

(*
1. We find all terms with GammaCC and thats it
*)
ExtractWeinbergTerms[eft_] := Module[{expanded, terms},
  expanded = Expand[eft];

  (* Get Weinberg Terms based on the prefactor with Times[] with GammaCC*)
  terms = DeleteDuplicates @ Cases[
    expanded,
    term_Times /; !FreeQ[term, GammaCC],
    Infinity
  ];

  (* If Times[] doesnt exist we look for NCM (non commutative multiplication chain) instead which should have GammaCC *)
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

(* 
Combine all non conjugated term prefactors into a compact C5 expression.
So A*LLHH + B*LLHH + C*LLHH -> A+B+C
 *)
CompactWeinbergCoefficient[terms_List] := Module[{pieces, combined},
  If[terms === {}, Return[Missing["NoHolomorphicTerms"]]];

  pieces = StripWeinbergOperatorStructure /@ terms;
  combined = Total[pieces];

  Quiet@Check[
    FactorTerms[Cancel[Together[combined]]],
    Simplify[combined]
  ]
];

(*
1. Check for existence of coefficient
2. Use ExtractWeinbergTerms to get our Weinberg terms
3. Sort our terms into holomorphic and conjugate terms
4. We use CompactWeinbergCoefficient to get our coefficient
5. We return a report on this
*)
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

  (* Check if weinberg charged conjugated gamma is present in eft *)
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

  (* Gets all Weinberg terms *)
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

  (* C5 uses the P_L LLHH sector, we look at holomorphic terms;
   keep P_R only as the Hermitian-conjugate audit. *)
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
  (* Obtain Weinberg Coefficient itself from our terms *)
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
