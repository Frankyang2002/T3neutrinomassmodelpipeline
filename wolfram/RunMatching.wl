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


(* ---------------------------------------------------------------------- *)
(* Weinberg-operator extraction                                            *)
(* ---------------------------------------------------------------------- *)

ClearAll[
  InternalHeadName, InternalSymbolName, MatcheteFieldName,
  CountMatcheteField, ContainsNamedSymbolQ, ContainsProjectorQ,
  BarredMatcheteFieldQ, WeinbergLHHTermQ, ExtractWeinbergTerms,
  StripWeinbergOperatorStructure, CompactWeinbergCoefficient,
  ExtractWeinbergCoefficient
];

(* These helpers deliberately use symbol names rather than contexts, matching
   the philosophy of PhysicsLaTeX.wl and making the extractor less sensitive
   to Matchete context changes. *)
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

(* First-stage Weinberg selector: use the same robust invariant that the
   previously validated pipeline used.  Do not impose assumptions about the
   post-EOMSimplify Matchete representation here.  The raw GammaCC sector is
   exported so that the precise LLHH structure can be canonicalised from what
   Matchete actually returns, rather than guessed. *)
WeinbergLHHTermQ[term_] := !FreeQ[Unevaluated[term], GammaCC];

(* Presence and term isolation are intentionally separate.  The presence test
   is the same robust GammaCC test used by the previously validated pipeline.
   A failure of the more detailed term parser must therefore never turn a real
   Weinberg operator into a false negative. *)
ExtractWeinbergTerms[eft_] := Module[{expanded, terms},
  expanded = Expand[eft];

  (* Match Mathematica's actual Times head directly.  This is more reliable
     than converting held heads to strings.  Each loop contribution to the
     Weinberg operator is an ordinary multiplicative term containing an NCM
     spinor chain with GammaCC.  Times is Flat, so these are the additive
     contributions even when Matchete wraps the complete EFT internally. *)
  terms = Cases[
    expanded,
    term_Times /; !FreeQ[term, GammaCC],
    Infinity
  ];

  (* Remove any accidental duplicate subexpressions while preserving the
     generated numerical/group factors. *)
  terms = DeleteDuplicates[terms];

  (* Extremely defensive fallback: a coefficient-free bare GammaCC chain is
     still a valid operator, although Matchete normally supplies a Times
     prefactor at one loop. *)
  If[terms === {} && !FreeQ[expanded, GammaCC],
    terms = DeleteDuplicates @ Cases[
      expanded,
      chain_NCM /; !FreeQ[chain, GammaCC],
      Infinity
    ]
  ];

  terms
];

(* Remove only the universal LLHH operator structure.  Numerical/group
   factors already generated by Matchete remain in the prefactor.  The full
   Weinberg sector is returned alongside this reduced coefficient so the
   operation is always auditable. *)
StripWeinbergOperatorStructure[term_] := Module[{stripped},
  (* These patterns are deliberately tied to the actual post-matching
     Matchete representation observed in c5_raw.txt.  Remove the universal
     LLHH spinor chain, the two Higgs fields and the two SU(2) CG tensors,
     while retaining every coupling, mass, logarithm and numerical/group
     factor. *)
  stripped = term /. {
    chain_NCM /; !FreeQ[chain, GammaCC] :> 1,
    Field[H, Scalar, inds_, derivs_] :> 1,
    CG[___] :> 1
  };
  Quiet@Check[Simplify[Expand[stripped]], stripped]
];

(* Combine the scalar prefactors of equivalent Weinberg contributions.
   Together/Cancel first removes the artificial product of mass-difference
   denominators generated term-by-term; FactorTerms then exposes the common
   Yukawa, lambdaT3, MF and loop factors. *)
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
  {present, terms, holomorphicTerms, conjugateTerms, sector,
   holomorphicSector, coefficient, status},

  (* This is the authoritative presence test. *)
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

  (* If GammaCC is present but detailed isolation fails, report that honestly
     rather than incorrectly claiming that the Weinberg operator is absent. *)
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

  (* Matchete returns the operator and its Hermitian conjugate separately.
     We define the C5 candidate from the P_L, unbarred-H sector, i.e. the
     conventional L^T C P_L L H H orientation.  The P_R sector is retained
     only as an audit/HC check and is not double-counted in C5. *)
  (* At this stage Matchete has already canonicalised the spinor chains and
     the projector appears literally as Proj[-1] or Proj[1].  Match that
     exact structure directly.  The earlier generic head/introspection helper
     was too defensive and failed on these otherwise ordinary expressions. *)
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
