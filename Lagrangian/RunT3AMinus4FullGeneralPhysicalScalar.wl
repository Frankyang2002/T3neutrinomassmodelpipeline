(* A(-4): native scalar threshold, local-only continuation, v8 CG trace and scalar-loop isolation.
   Call from RunT3FullGeneralTwoThresholdNative.wl after it has saved
   A_EFT1_tree_contracted.wxf and A_full_F_hard_one_loop_preserved.wxf.

   This script restores a consistent Matchete registry, changes ONLY the
   EFT order of the two scalar masses from Light (1) to Heavy (0), audits
   the change and integrates out BOTH BSM scalars at ONE threshold.

   Inputs are the THREE-FLAVOUR tree EFT1 Wilson tensors, obtained at ONE
   fermion threshold (all three heavy flavours). Stage-one hard O(hbar) operators
   is preserved but deliberately excluded from this scalar-soft Match.

   Neither finite renormalisation nor complete EFT1 Wilson RG is performed.
   The output is scalar-soft only and MUST NOT be called physical C5.
*)

If[!ValueQ[t3out] || !ValueQ[t3model] ||
   !AssociationQ[t3model] || DownValues[t3Save] === {} ||
   DownValues[t3Stop] === {},
  Print["ERROR: A(-4) continuation must be invoked from the native parent."];
  Exit[4]
];

t3stage = "A-rebuild-full-general-scalar-registry";
t3ABuild = CheckAbort[UsingFrontEnd[BuildT3Lagrangian[t3model]], $Aborted];
If[!AssociationQ[t3ABuild] ||
   Lookup[t3ABuild, "Status", ""] =!= "Success" ||
   !TrueQ[Lookup[t3ABuild, "FullValidation", False]] ||
   !TrueQ[Lookup[t3ABuild, "AllScalarExtensionsRequested", False]] ||
   Lookup[t3ABuild, "AdditionalQuarticsIncluded", -1] =!= 0 ||
   Lookup[t3ABuild, "AdditionalQuadraticsIncluded", -1] =!= 0,
  t3Stop["A(-4) full-general registry rebuild failed"]
];

(* BuildT3Lagrangian resets Matchete, and the initial F registry was
   correct for the UV threshold only.  Re-register the physical EFT1
   with precisely S1 and S2 Heavy, with no fermion field present. *)
If[ReclassifyT3Fields[t3model, {"S1", "S2"}, {"S1", "S2"}] === $Failed,
  t3Stop["A(-4) joint heavy-scalar registry reconstruction failed"]
];
Quiet[DefineCoupling[T3EFT1W12, Indices -> {Flavor, Flavor},
  SelfConjugate -> False]];
t3ACouplings = GetCouplings[];
If[!AssociationQ[t3ACouplings] ||
   !KeyExistsQ[t3ACouplings, T3EFT1W12],
  t3Stop["A(-4) flavour Wilson tensor was not registered in physical EFT1"]
];

(* Rebuild all invariant tensor conventions, not only the triangle CG.
   Scalar potential reconstruction is necessary for the portal/self/cross
   quartic tensors carried by the saved EFT1. *)
If[DefineT3InvariantCGs[t3model] === $Failed ||
   !ListQ[BuildGeneralScalarPotentialCandidates[t3model]],
  t3Stop["A(-4) scalar-potential CG reconstruction failed"]
];
t3AFields = GetFields[];
If[!AssociationQ[t3AFields] ||
   !And @@ (TrueQ[
     Lookup[Lookup[t3AFields, #, <||>], Heavy, False]] & /@
       {NewScalar1, NewScalar2}) ||
   KeyExistsQ[t3AFields, NewFermion],
  t3Stop["A(-4) physical second-threshold field inventory invalid"]
];

(* Import only after rebuilding the complete field, coupling and CG
   metadata. WXF stores expressions, not Matchete's private registries. *)
t3stage = "A-rehydrate-three-flavour-EFT1-tree";
t3AInputPath = FileNameJoin[{t3out, "A_EFT1_tree_contracted.wxf"}];
t3AHardPath = FileNameJoin[{t3out,
  "A_full_F_hard_one_loop_preserved.wxf"}];
If[!And @@ (FileExistsQ /@ {t3AInputPath, t3AHardPath}),
  t3Stop["A(-4) contracted tree or preserved F hard source missing"]
];
t3ARaw = Quiet @ Check[Import[t3AInputPath, "WXF"], $Failed];
If[t3ARaw === $Failed,
  t3Stop["A(-4) could not import three-flavour tree EFT1 WXF"]
];
t3ATree = Quiet @ Check[
  ToExpression[ToString[t3ARaw, InputForm], InputForm], $Failed];
If[t3ATree === $Failed || FreeQ[t3ATree, T3EFT1W12] ||
   FreeQ[t3ATree, GammaCC] || !FreeQ[t3ATree, NewFermion] ||
   !FreeQ[t3ATree,
     HoldPattern[Coupling[MF, {Index[_, NFlavor]}, 0]]],
  t3Stop["A(-4) rehydrated EFT1 tree fails exact Wilson inventory"]
];

t3AUsedCGs = DeleteDuplicates @ Cases[
  HoldComplete[t3ATree],
  CG[name_, ___] :> Replace[name, Bar[z_] :> z], Infinity];
t3AMissingCGs = Select[t3AUsedCGs, Function[name,
  Module[{reps = Quiet @ Check[
    Matchete`PackageScope`$CGproperties[name, Indices], $Failed]},
    reps === $Failed || MissingQ[reps]]
]];
If[t3AMissingCGs =!= {},
  t3Stop["Unregistered A(-4) CG tensor in saved EFT1",
    <|"missing_CGs" -> (ToString[#, InputForm] & /@ t3AMissingCGs)|>]
];

(* Verified failure in previous run:
       Coupling[MS1,{},1]^2 and Coupling[MS2,{},1]^2
   were carried from the F threshold, while the scalar registry at mu_S
   has Mass->{Heavy,MS1/MS2}.  Change only the mass EFT-order LABEL,
   never the mass or a scalar coupling numerical value. *)
t3stage = "A-reversible-scalar-mass-power-counting";
t3AMassSymbols = {MS1, MS2};
t3AOriginalMassCounts = AssociationThread[{ "MS1", "MS2" }, {
  Count[t3ATree, HoldPattern[Coupling[MS1, {}, 1]], Infinity],
  Count[t3ATree, HoldPattern[Coupling[MS2, {}, 1]], Infinity]
}];
If[!AllTrue[Values[t3AOriginalMassCounts], # > 0 &],
  t3Stop["A(-4) missing expected light scalar mass terms",
    <|"source_mass_counts" -> t3AOriginalMassCounts|>]
];
t3AMatchInput = t3ATree /. {
  HoldPattern[Coupling[MS1, {}, 1]] :> Coupling[MS1, {}, 0],
  HoldPattern[Coupling[MS2, {}, 1]] :> Coupling[MS2, {}, 0]
};
t3ADemoted = t3AMatchInput /. {
  HoldPattern[Coupling[MS1, {}, 0]] :> Coupling[MS1, {}, 1],
  HoldPattern[Coupling[MS2, {}, 0]] :> Coupling[MS2, {}, 1]
};
t3APromotedMassCounts = AssociationThread[{ "MS1", "MS2" }, {
  Count[t3AMatchInput, HoldPattern[Coupling[MS1, {}, 0]], Infinity],
  Count[t3AMatchInput, HoldPattern[Coupling[MS2, {}, 0]], Infinity]
}];
t3AMassValid = TrueQ[t3ADemoted === t3ATree] &&
  t3AOriginalMassCounts === t3APromotedMassCounts &&
  FreeQ[t3AMatchInput, HoldPattern[Coupling[MS1, {}, 1]]] &&
  FreeQ[t3AMatchInput, HoldPattern[Coupling[MS2, {}, 1]]];
t3AAuditPath = FileNameJoin[{t3out,
  "A_joint_scalar_mass_power_counting_audit.json"}];
Export[t3AAuditPath, <|
  "model" -> "A_-4",
  "source" -> "A_EFT1_tree_contracted.wxf",
  "source_light_mass_EFT_order" -> 1,
  "target_heavy_mass_EFT_order" -> 0,
  "source_mass_counts" -> t3AOriginalMassCounts,
  "promoted_mass_counts" -> t3APromotedMassCounts,
  "exact_reversible_mass_promotion" -> t3AMassValid,
  "numerical_mass_values_changed" -> False,
  "fermion_flavours_in_tree_Wilson_definition" -> 3,
  "both_scalars_integrated_in_one_match" -> True,
  "F_hard_inherited" -> False,
  "EFT1_Wilson_RG_included" -> False,
  "renormalized_physical_C5_verified" -> False
|>, "RawJSON"];
If[!t3AMassValid,
  t3Stop["A(-4) scalar-mass promotion failed exact reversible audit",
    <|"audit_path" -> t3AAuditPath|>]
];
t3Save["A_EFT1_joint_heavy_scalar_match_input.wxf", t3AMatchInput];

(* The previous run passed the reversible heavy-mass audit, but
   CheckLagrangian returned False. Earlier A fresh-kernel logs also showed
   CheckLagrangian::Hermiticity. Do not ignore that result. First try
   mathematically equivalent dummy-index relabelling and independently
   reconstructed SM kinetic/Yukawa terms, proving term-wise alpha equivalence
   before allowing either candidate to be matched. If neither passes,
   identify whether the EFT1 W or lambdaT3 sectors cause the failure, and
   STOP without modifying the physics or inventing a C5. *)
t3stage = "A-heavy-scalar-CheckLagrangian";
ClearAll[t3ACheckSafe, t3AParts, t3AAlphaKey, t3AAlphaKeys,
  t3ABSMTermQ, t3AZeroWilson, t3AZeroMix];
t3ACheckSafe[expr_] := TimeConstrained[
  CheckAbort[UsingFrontEnd[CheckLagrangian[expr]], $Aborted],
  900, $TimedOut];
t3AParts[expr_] := Module[{x = Expand[expr]},
  If[x === 0, {}, If[Head[x] === Plus, List @@ x, {x}]]];
(* All placeholder indices carry the ORIGINAL representation, and each
   distinct index is renamed consistently only within its OWN monomial.
   This is an audit, never an instruction to change a contraction. *)
t3AAlphaKey[term_] := Module[{indices, rules},
  indices = DeleteDuplicates @ Cases[HoldComplete[term],
    HoldPattern[Index[_, _]], Infinity];
  rules = Table[
    With[{old = indices[[k]], rep = indices[[k, 2]], slot = k},
      old -> Index[Symbol["Global`T3ACanonicalDummy" <> ToString[slot]],
        rep]], {k, Length[indices]}];
  ToString[term /. rules, InputForm]];
t3AAlphaKeys[expr_] := Sort[t3AAlphaKey /@ t3AParts[expr]];
t3ABSMTermQ[term_] := !FreeQ[term, NewScalar1 | NewScalar2 | T3EFT1W12];

t3ACheckOriginal = t3ACheckSafe[t3AMatchInput];
t3ACheck = t3ACheckOriginal;
t3ACheckMode = "ExactSavedTree";
t3ASelectedInput = t3AMatchInput;
t3ARelabelEquivalent = False;
t3ACheckRelabel = Missing["NotNeeded"];
t3ASMAlphaEquivalent = False;
t3ACheckSMRebuild = Missing["NotNeeded"];
t3ASMPreserved = False;

If[!TrueQ[t3ACheck],
  (* No operators are dropped: only dummy indices are canonicalised. *)
  t3ARelabelCandidate = CheckAbort[
    UsingFrontEnd[RelabelIndices[Expand[t3AMatchInput]]], $Aborted];
  If[!MemberQ[{$Failed, $Aborted}, t3ARelabelCandidate],
    t3ARelabelEquivalent = TrueQ[
      t3AAlphaKeys[t3ARelabelCandidate] === t3AAlphaKeys[t3AMatchInput]];
    If[t3ARelabelEquivalent,
      t3ACheckRelabel = t3ACheckSafe[t3ARelabelCandidate];
      If[TrueQ[t3ACheckRelabel],
        t3ACheck = True;
        t3ACheckMode = "ProvenDummyIndexRelabel";
        t3ASelectedInput = t3ARelabelCandidate]
    ]
  ];
];

If[!TrueQ[t3ACheck],
  (* The SM terms are recreated from the same validated UV model, but may
     replace the saved SM sector ONLY if the entire term multiset is
     identical modulo capture-avoiding dummy-index renaming. *)
  t3AAllTerms = t3AParts[t3AMatchInput];
  t3ASMTerms = Select[t3AAllTerms, !t3ABSMTermQ[#] &];
  t3ABSMTerms = Select[t3AAllTerms, t3ABSMTermQ];
  t3ASMOriginal = Expand[t3ABuild["LSM"]];
  t3ASMAlphaEquivalent = TrueQ[
    Sort[t3AAlphaKey /@ t3ASMTerms] === t3AAlphaKeys[t3ASMOriginal]];
  If[t3ASMAlphaEquivalent,
    t3ASMRebuilt = Expand[Total[t3ABSMTerms] + t3ASMOriginal];
    t3ASMPreserved = TrueQ[
      t3AAlphaKeys[t3ASMRebuilt] === t3AAlphaKeys[t3AMatchInput]];
    If[t3ASMPreserved,
      t3ACheckSMRebuild = t3ACheckSafe[t3ASMRebuilt];
      If[TrueQ[t3ACheckSMRebuild],
        t3ACheck = True;
        t3ACheckMode = "ProvenEquivalentSMRepresentation";
        t3ASelectedInput = t3ASMRebuilt]
    ]
  ]
];

(* ---------------------------------------------------------------------- *)
(* V5: Verify the saved EFT1 against a source-reconstructed scalar theory.
   The F-threshold TREE can only add fermion-exchange operators carrying
   light leptons/scalars. Do not assume the stored SM sector is identical
   to the SM built in the fresh kernel: inspect it explicitly.

   The older alpha-equivalence key used first-encounter index ordering.
   Mathematica's canonical Times ordering can itself depend on dummy names.
   Here we take the minimum over ALL within-representation permutations.
   This is invariant under dummy renaming, while it preserves every
   contraction, CG orientation, coefficient, and operator. A candidate is
   NEVER matched unless the *whole Lagrangian* is proven equivalent.
*)
ClearAll[t3AGraphIndexKey, t3AGraphKeys, t3AGraphDifference,
  t3AHasWilson, t3ASourceNoF, t3APermutationLimit];
t3APermutationLimit = 10000;
t3AGraphIndexKey[term_] := Module[
  {indices, groups, choices, cost, strings, rules, positions, perms},
  indices = DeleteDuplicates @ Cases[HoldComplete[term],
    HoldPattern[Index[_, _]], Infinity];
  If[indices === {}, Return[ToString[term, InputForm]]];
  groups = SortBy[GatherBy[indices, ToString[Last[#], InputForm] &],
    ToString[Last[First[#]], InputForm] &];
  cost = Times @@ (Factorial[Length[#]] & /@ groups);
  If[cost > t3APermutationLimit,
    Return[Missing["TooManyDummyPermutations", cost]]];
  choices = Tuples[Permutations /@ groups];
  strings = Table[
    rules = Flatten[Table[
      With[{rep = Last[groups[[g, 1]]], seq = choice[[g]], group = g},
        Table[With[{old = seq[[j]], slot = j},
          old -> Index[Symbol["Global`T3AGraphDummy" <>
            ToString[group] <> "n" <> ToString[slot]], rep]],
          {j, Length[seq]}]],
      {g, Length[groups]}], 1];
    ToString[term /. rules, InputForm],
    {choice, choices}];
  First[Sort[strings]]
];
t3AGraphKeys[expr_] := Sort[t3AGraphIndexKey /@ t3AParts[expr]];
t3AGraphDifference[base_List, proposal_List] := Module[
  {a = Counts[base], b = Counts[proposal], labels, diffs},
  labels = Union[Keys[a], Keys[b]];
  diffs = Select[Table[<|"key" -> key,
      "saved" -> Lookup[a, key, 0],
      "proposed" -> Lookup[b, key, 0]|>, {key, labels}],
    #["saved"] =!= #["proposed"] &];
  <|"distinct_mismatched_monomials" -> Length[diffs],
    "saved_unmatched_term_count" ->
      Total[Max[0, #["saved"] - #["proposed"]] & /@ diffs],
    "proposal_unmatched_term_count" ->
      Total[Max[0, #["proposed"] - #["saved"]] & /@ diffs],
    "previews" -> (Map[Function[row, <|
        "saved" -> row["saved"], "proposed" -> row["proposed"],
        "key_preview" -> StringTake[row["key"], UpTo[450]]|>],
       Take[diffs, UpTo[8]]])|>
];
t3AHasWilson[term_] := !FreeQ[term,
  HoldPattern[Coupling[T3EFT1W12, ___]]];
t3ASourceNoF[term_] := FreeQ[term, NewFermion];
(* Exact syntactic CG orientation inventory. Bar[...] is deliberately
   NOT stripped; doing so would incorrectly assume tensor identities. *)
ClearAll[t3ACGOrientations];
t3ACGOrientations[expr_] := Counts[
  (ToString[#, InputForm] & /@ Cases[HoldComplete[expr],
    CG[name_, ___] :> name, Infinity])];

(* An optional diagnostic computation may use a UV-derived Hermitian
   Wilson completion if the literal F-boundary comparison is inconclusive.
   It NEVER certifies physical EFT1 or finite C5. Off by default. *)
t3AExploratoryEnabled =
  Environment["T3_A4_EXPLORATORY_SOURCE_MATCH"] === "1";
t3ASelectedExploratory = False;
t3AExploratoryReason = "NotAttempted";
t3AWPRCount = 0;
t3AWPLCount = 0;
t3AWHcConstructed = False;
t3AWHcPaired = False;
t3AWHolomorphicGraphPreserved = False;
t3AWHcCheck = Missing["NotRun"];
t3AWHcCandidate = Missing["NotBuilt"];
t3AGraphAudit = <|"attempted" -> False|>;
t3ACheckReconstructed = Missing["NotNeeded"];
t3ACheckScalarBase = Missing["NotNeeded"];
t3AFullGraphEquivalent = False;
t3ASelectedFromSource = False;
If[!TrueQ[t3ACheck],
  t3stage = "A-source-reconstruction-and-alpha-graph-audit";
  (* LSM and all renormalizable scalar interactions come from the same
     UV builder that previously passed its own CheckLagrangian. The F
     Yukawas are excluded ONLY from this independent reconstruction. *)
  t3ASourceScalarInteractions = Select[t3AParts[t3ABuild["LInt"]],
    t3ASourceNoF];
  t3ASourceFInteractions = Select[t3AParts[t3ABuild["LInt"]],
    !t3ASourceNoF[#] &];
  t3ASavedWTerms = Select[t3AParts[t3AMatchInput], t3AHasWilson];
  If[Length[t3ASavedWTerms] =!= 2 ||
     Length[t3ASourceFInteractions] === 0,
    t3Stop["A(-4) source reconstruction has unexpected F/W inventory",
      <|"wilson_terms" -> Length[t3ASavedWTerms],
        "UV_F_interaction_terms" -> Length[t3ASourceFInteractions]|>]
  ];
  (* These are new expressions in the currently registered HEAVY-scalar
     theory; no F or its free Lagrangian is reintroduced. *)
  t3AScalarBase = RelabelIndices[Expand[
    t3ABuild["LSM"] + FreeLag[NewScalar1, NewScalar2] +
      Total[t3ASourceScalarInteractions]]];
  t3AReconstructed = RelabelIndices[
    Expand[t3AScalarBase + Total[t3ASavedWTerms]]];
  t3ASavedSMGraphKeys = t3AGraphKeys[Total[t3ASMTerms]];
  t3ABuiltSMGraphKeys = t3AGraphKeys[t3ABuild["LSM"]];
  t3ASMGraphEquivalent = TrueQ[
    t3ASavedSMGraphKeys === t3ABuiltSMGraphKeys];
  t3ASavedGraphKeys = t3AGraphKeys[t3AMatchInput];
  t3ARebuiltGraphKeys = t3AGraphKeys[t3AReconstructed];
  t3AGraphComputable = FreeQ[
    {t3ASavedGraphKeys, t3ARebuiltGraphKeys},
    Missing["TooManyDummyPermutations", _], Infinity];
  t3AFullGraphEquivalent = t3AGraphComputable &&
    TrueQ[t3ASavedGraphKeys === t3ARebuiltGraphKeys];
  (* V6: always test the complete SOURCE-rebuilt candidate as a
     DIAGNOSTIC, even when exact graph equivalence is not established.
     Previously this field was NotNeeded, concealing whether the saved
     scalar sector or the saved F-induced Wilson pair caused Hermiticity. *)
  t3ACheckScalarBase = t3ACheckSafe[t3AScalarBase];
  t3ACheckReconstructed = t3ACheckSafe[t3AReconstructed];
  If[t3AFullGraphEquivalent,
    If[TrueQ[t3ACheckReconstructed],
      t3ASelectedInput = t3AReconstructed;
      t3ACheck = True;
      t3ACheckMode = "SourceReconstructed_ProvenWholeEFTAlphaGraphEquivalent";
      t3ASelectedFromSource = True]
  ];

  (* V6: test the Hermitian partner implied by the holomorphic MATCHED
     F-threshold Wilson. A mismatch with the saved opposite chirality is
     recorded; it is NOT silently declared to be an equivalent operator.
     PlusHc uses Matchete's actual conjugation operation/CG registry. *)
  t3AWPRTerms = Select[t3ASavedWTerms,
    StringContainsQ[ToString[#, InputForm], "Proj[1]"] &];
  t3AWPLTerms = Select[t3ASavedWTerms,
    StringContainsQ[ToString[#, InputForm], "Proj[-1]"] &];
  t3AWPRCount = Length[t3AWPRTerms];
  t3AWPLCount = Length[t3AWPLTerms];
  If[t3AWPRCount === 1 && t3AWPLCount === 1,
    t3AWPlusHc = CheckAbort[
      UsingFrontEnd[PlusHc[First[t3AWPRTerms]]], $Aborted];
    t3AWHcConstructed = !MemberQ[{$Failed, $Aborted}, t3AWPlusHc] &&
      FreeQ[t3AWPlusHc, PlusHc] &&
      !FreeQ[t3AWPlusHc, T3EFT1W12];
    If[t3AWHcConstructed,
      t3AWPlusHcTerms = t3AParts[t3AWPlusHc];
      t3AWHoloTerms = Select[t3AWPlusHcTerms,
        StringContainsQ[ToString[#, InputForm], "Proj[1]"] &];
      t3AWHcTerms = Select[t3AWPlusHcTerms,
        StringContainsQ[ToString[#, InputForm], "Proj[-1]"] &];
      t3AWHcPaired = Length[t3AWPlusHcTerms] === 2 &&
        Length[t3AWHoloTerms] === 1 && Length[t3AWHcTerms] === 1;
      If[t3AWHcPaired,
        t3AWHolomorphicGraphPreserved = TrueQ[
          t3AGraphKeys[Total[t3AWHoloTerms]] ===
          t3AGraphKeys[Total[t3AWPRTerms]]];
        t3AWHcPairGraphMatchesSaved = TrueQ[
          t3AGraphKeys[Total[t3AWPlusHcTerms]] ===
          t3AGraphKeys[Total[t3ASavedWTerms]]];
        t3AWHcCandidate = RelabelIndices[Expand[
          t3AScalarBase + t3AWPlusHc]];
        t3AWHcCheck = t3ACheckSafe[t3AWHcCandidate];
        t3Save["A_source_hermitian_W_completion_DIAGNOSTIC.wxf",
          t3AWHcCandidate];
      ];
    ];
  ];

  (* Non-authoritative numerical probe: the source scalar theory is
     independently valid; the holomorphic Wilson is unmodified; its HC
     partner is generated by Matchete. Nonetheless, a missing exact
     equivalence certificate prevents calling this the actual EFT1.
     Only run it with explicit experimental opt-in. *)
  t3AExploratoryQualified = TrueQ[t3ACheckScalarBase] &&
    t3AWHcConstructed && t3AWHcPaired &&
    t3AWHolomorphicGraphPreserved && TrueQ[t3AWHcCheck];
  t3AExploratoryReason = Which[
    !t3AExploratoryEnabled, "NotOptedIn",
    !t3AExploratoryQualified, "HermitianSourceCandidateDidNotValidate",
    t3AFullGraphEquivalent, "LiteralGraphAlreadyEquivalent",
    True, "QualifiedNonEquivalentDiagnosticOnly"];
  If[!TrueQ[t3ACheck] && t3AExploratoryEnabled &&
     t3AExploratoryQualified && !t3AFullGraphEquivalent,
    t3ASelectedInput = t3AWHcCandidate;
    t3ASelectedExploratory = True;
    t3ACheck = True;
    t3ACheckMode =
      "ExploratoryUVScalarBasePlusExactHolomorphicWAndMatcheteHC";
  ];
  t3AGraphAudit = <|
    "attempted" -> True,
    "UV_F_interaction_terms_excluded_from_EFT1" ->
      Length[t3ASourceFInteractions],
    "UV_scalar_interaction_terms_preserved" ->
      Length[t3ASourceScalarInteractions],
    "saved_W_terms_added" -> Length[t3ASavedWTerms],
    "scalar_base_CheckLagrangian" -> ToString[t3ACheckScalarBase, InputForm],
    "complete_rebuilt_CheckLagrangian" ->
      ToString[t3ACheckReconstructed, InputForm],
    "SM_sector_alpha_graph_equivalent_to_UV_LSM" ->
      t3ASMGraphEquivalent,
    "saved_SM_CG_orientations" ->
      t3ACGOrientations[Total[t3ASMTerms]],
    "reconstructed_SM_CG_orientations" ->
      t3ACGOrientations[t3ABuild["LSM"]],
    "saved_W_CG_orientations" ->
      t3ACGOrientations[Total[t3ASavedWTerms]],
    "SM_sector_graph_mismatch" ->
      t3AGraphDifference[t3ASavedSMGraphKeys, t3ABuiltSMGraphKeys],
    "graph_key_computable" -> t3AGraphComputable,
    "whole_EFT_graph_alpha_equivalent" -> t3AFullGraphEquivalent,
    "full_EFT_reconstruction_used_for_matching" -> t3ASelectedFromSource,
    "reconstructed_full_CheckLagrangian_DIAGNOSTIC_ONLY" ->
      ToString[t3ACheckReconstructed, InputForm],
    "saved_W_PR_terms" -> t3AWPRCount,
    "saved_W_PL_terms" -> t3AWPLCount,
    "W_from_matched_PR_PlusHc_constructed" -> t3AWHcConstructed,
    "W_PlusHc_two_chiralities" -> t3AWHcPaired,
    "W_PlusHc_CG_orientations" ->
      If[t3AWHcPaired,
        t3ACGOrientations[t3AWPlusHc], <||>],
    "W_matched_PR_graph_preserved" -> t3AWHolomorphicGraphPreserved,
    "W_PlusHc_pair_matches_saved_tree" ->
      If[ValueQ[t3AWHcPairGraphMatchesSaved],
        t3AWHcPairGraphMatchesSaved, Missing["NotChecked"]],
    "source_Hermitian_W_complete_CheckLagrangian_DIAGNOSTIC_ONLY" ->
      ToString[t3AWHcCheck, InputForm],
    "exploratory_opt_in_requested" -> t3AExploratoryEnabled,
    "exploratory_qualified" -> t3AExploratoryQualified,
    "exploratory_selected" -> t3ASelectedExploratory,
    "exploratory_reason" -> t3AExploratoryReason,
    "literal_F_boundary_equivalence_certified" -> t3AFullGraphEquivalent,
    "mismatch" -> If[t3AGraphComputable,
      t3AGraphDifference[t3ASavedGraphKeys, t3ARebuiltGraphKeys],
      <|"reason" -> "graph-permutation-limit"|>]|>;
  Export[FileNameJoin[{t3out, "A_source_reconstruction_audit.json"}],
    t3AGraphAudit, "RawJSON"];
  (* Retain independently checkable WXF provenance even on failure. *)
  t3Save["A_source_reconstructed_complete_EFT1_DIAGNOSTIC.wxf",
    t3AReconstructed];
];

(* Diagnostic-only projections. Never pass either projected Lagrangian to
   the physical scalar matching stage. *)
t3ACheckWithoutW = Missing["NotNeeded"];
t3ACheckWithoutWAndMix = Missing["NotNeeded"];
If[!TrueQ[t3ACheck],
  t3AZeroWilson = Expand[t3AMatchInput /. {
    HoldPattern[Bar[Coupling[T3EFT1W12, ___]]] :> 0,
    HoldPattern[Coupling[T3EFT1W12, ___]] :> 0}];
  t3ACheckWithoutW = t3ACheckSafe[t3AZeroWilson];
  t3AZeroMix = Expand[t3AZeroWilson /. {
    HoldPattern[Bar[Coupling[lambdaT3, ___]]] :> 0,
    HoldPattern[Coupling[lambdaT3, ___]] :> 0}];
  t3ACheckWithoutWAndMix = t3ACheckSafe[t3AZeroMix]
];

Export[FileNameJoin[{t3out, "A_joint_scalar_lagrangian_audit.json"}], <|
  "model" -> "A_-4",
  "CheckLagrangian_original" -> ToString[t3ACheckOriginal, InputForm],
  "CheckLagrangian_relabelled" -> ToString[t3ACheckRelabel, InputForm],
  "CheckLagrangian_equivalent_SM_rebuilt" ->
    ToString[t3ACheckSMRebuild, InputForm],
  "CheckLagrangian_without_W_DIAGNOSTIC_ONLY" ->
    ToString[t3ACheckWithoutW, InputForm],
  "CheckLagrangian_without_W_and_lambdaT3_DIAGNOSTIC_ONLY" ->
    ToString[t3ACheckWithoutWAndMix, InputForm],
  "relabelled_exactly_alpha_equivalent" -> t3ARelabelEquivalent,
  "saved_SM_exactly_alpha_equivalent_to_built_LSM" ->
    t3ASMAlphaEquivalent,
  "SM_rebuild_preserves_complete_EFT_alpha_equivalence" ->
    t3ASMPreserved,
  "SourceReconstruction" -> t3AGraphAudit,
  "exploratory_non_equivalent_source_used" -> t3ASelectedExploratory,
  "CheckLagrangian_result" -> ToString[t3ACheck, InputForm],
  "CheckLagrangian_passed" -> TrueQ[t3ACheck],
  "accepted_canonicalisation" -> t3ACheckMode,
  "mass_promotion_reversible" -> t3AMassValid,
  "projection_used_for_physical_match" -> False,
  "matched_C5_physically_certified" -> False|>, "RawJSON"];
If[!TrueQ[t3ACheck],
  t3Stop["A(-4) CheckLagrangian still false; inspect SourceReconstruction mismatch and scalar_base_CheckLagrangian in A_source_reconstruction_audit.json",
    <|"check_result" -> ToString[t3ACheck, InputForm],
      "mass_audit_path" -> t3AAuditPath,
      "lagrangian_audit_path" -> FileNameJoin[{t3out,
        "A_joint_scalar_lagrangian_audit.json"}]|>]
];
(* The authoritative path is alpha-graph equivalent and CheckLagrangian
   verified. If the opt-in exploratory path was selected, the candidate is
   STRICTLY non-authoritative and is written under diagnostic filenames. *)
t3AMatchInput = t3ASelectedInput;
t3AOutputPrefix = If[t3ASelectedExploratory,
  "A_exploratory_source_soft_", "A_scalar_soft_only_"];
t3Save[If[t3ASelectedExploratory,
  "A_source_hermitian_exploratory_checked_input.wxf",
  "A_EFT1_joint_heavy_scalar_checked_input.wxf"], t3AMatchInput];

(* One common scalar threshold (S1,S2). Scalar-soft only.  In v6,
   CheckLagrangian passed for the exploratory Hermitian source, but Match
   returned an atom (Head Symbol).  A Symbol head is NOT a diagnosis:
   the actual expression may be $Failed, $Aborted, or $TimedOut.

   V7 preserves the same selected theory. No CG, Yukawa, or scalar coupling
   is projected/altered. It now runs a tree-level Match first and saves the
   EXACT return identity from every raw Matchete call, plus message labels.
   An optional Wilson-free diagnostic runs ONLY after a failed main Match.
   That diagnostic is never used to supply physical Wilson coefficients.
*)
t3stage = "A-joint-heavy-scalar-match-checkpoints";
(* V9: TraceScan is used with only its DOCUMENTED TraceScan options.
   TraceOriginal is an option of Trace, not TraceScan, and passing it left
   TraceScan unevaluated in v8 (TraceScan::optx). The visitor is purely
   observational: it records attempts to define CG tensors before Matchete
   rejects a non-invariant one; no invariant validation is disabled. *)
ClearAll[t3ACGTraceVisitor];
SetAttributes[t3ACGTraceVisitor, HoldAllComplete];
t3ACGTraceCount = 0;
t3ACGTraceStream = None;
t3ACGTracePath = FileNameJoin[{t3out,
  "A_one_loop_DefineCG_calls_DIAGNOSTIC.txt"}];
t3ACGTraceEnabled = Environment["T3_A4_TRACE_DEFINECG"] === "1";
t3ACGTraceVisitor[held_] := Module[{line},
  t3ACGTraceCount++;
  If[Head[t3ACGTraceStream] === OutputStream &&
      t3ACGTraceCount <= 2000,
    line = ToString[Unevaluated[held], InputForm];
    WriteString[t3ACGTraceStream,
      ToString[t3ACGTraceCount] <> "\t" <>
      StringTake[line, UpTo[30000]] <> "\n"];
    Flush[t3ACGTraceStream]
  ];
];
ClearAll[t3AMatchAttempt, t3AMatchReport, t3AMatchFailureReport,
  t3AMatchStatusName, t3AMatchFieldNames];
t3AMatchFieldNames[expr_] := If[AtomQ[expr], {},
  Select[{NewScalar1, NewScalar2},
    !FreeQ[expr, HoldPattern[Field[#, ___]]] &]];
t3AMatchStatusName[value_] := Which[
  value === $Failed, "$Failed",
  value === $Aborted, "$Aborted",
  value === $TimedOut, "$TimedOut",
  AtomQ[value], "OtherAtomicResult",
  Head[value] === TraceScan, "UnevaluatedTraceScan_InvalidInstrumentation",
  True, "ExpressionReturned"];
(* Block the message LIST locally; do not suppress messages going to the
   normal stdout/stderr log. This is a supplemental machine-readable trace. *)
t3AMatchAttempt[expr_, loop_Integer, limit_Integer, label_String] := Module[
  {value = $Failed, started, seconds, messages = {}, headText, name,
   heavyRemaining, resultValid, msgPreview},
  Print["[A4-match] BEGIN ", label, "; loop order=", loop,
    "; EFTOrder=5; timeout=", limit, " seconds"];
  started = AbsoluteTime[];
  value = Block[{$MessageList = {}},
    Module[{v},
      v = TimeConstrained[
        CheckAbort[UsingFrontEnd[
          If[t3ACGTraceEnabled && loop === 1 &&
              label === "selected-full-theory-one-loop",
            TraceScan[t3ACGTraceVisitor,
              Match[expr, EFTOrder -> 5, LoopOrder -> loop],
              HoldPattern[DefineCG[___]],
              TraceInternal -> True],
            Match[expr, EFTOrder -> 5, LoopOrder -> loop]
          ]], $Aborted],
        limit, $TimedOut];
      messages = ToString[#, InputForm] & /@ $MessageList;
      v]];
  seconds = N[Max[0, AbsoluteTime[] - started]];
  name = t3AMatchStatusName[value];
  headText = ToString[Head[value], InputForm];
  heavyRemaining = ToString[#, InputForm] & /@ t3AMatchFieldNames[value];
  resultValid = !AtomQ[value] && heavyRemaining === {} &&
    (loop === 0 || !FreeQ[value, hbar]);
  msgPreview = StringTake[#, UpTo[1000]] & /@ Take[messages, UpTo[20]];
  Print["[A4-match] END ", label, "; returned=", name,
    "; head=", headText, "; elapsed=", seconds,
    "; heavy fields remain=", heavyRemaining,
    "; explicit hbar=", !FreeQ[value, hbar]];
  <|"result" -> value, "report" -> <|
     "label" -> label, "EFTOrder" -> 5, "LoopOrder" -> loop,
     "elapsed_seconds" -> seconds, "return_status" -> name,
     "return_head" -> headText,
     "return_InputForm_if_atomic" -> If[AtomQ[value],
       ToString[value, InputForm], "NonatomicExpression"],
     "non_atomic_EFT_returned" -> !AtomQ[value],
     "explicit_hbar" -> !FreeQ[value, hbar],
     "remaining_heavy_scalar_fields" -> heavyRemaining,
     "message_count_captured" -> Length[messages],
     "message_labels_preview" -> msgPreview,
     "structurally_valid_for_requested_order" -> resultValid,
     "native_matchete_used" -> True|>|>
];
t3AMatchReport = <|
  "format" -> "t3-a-minus4-joint-scalar-match-triage-v7",
  "model" -> "A_-4", "source_theory" -> t3ACheckMode,
  "exploratory_non_equivalent_source" -> TrueQ[t3ASelectedExploratory],
  "saved_F_EFT1_boundary_equivalence_certified" ->
    !TrueQ[t3ASelectedExploratory] && TrueQ[t3AFullGraphEquivalent],
  "CheckLagrangian_passed" -> TrueQ[t3ACheck],
  "F_threshold_recomputed" -> False,
  "fermion_generations_in_W" -> 3,
  "scalars_integrated_together" -> {"S1", "S2"},
  "matching_at_one_loop_succeeded" -> False,
  "F_hard_inherited" -> False,
  "EFT1_Wilson_running_included" -> False,
  "physical_C5_certified" -> False|>;
t3AMatchDiagFile = FileNameJoin[{t3out,
  "A_joint_scalar_match_diagnostics.json"}];
ClearAll[t3ASaveMatchReport];
t3ASaveMatchReport[] := Export[t3AMatchDiagFile,
  t3AMatchReport, "RawJSON"];
t3ASaveMatchReport[];
(* Tree matching is a structural control. It must not be confused with
   tree-level neutrino-mass generation, or used as a substitute for one loop. *)
t3ATreeProbe = t3AMatchAttempt[t3AMatchInput, 0, 900,
  "selected-full-theory-tree-control"];
t3AMatchReport = Join[t3AMatchReport,
  <|"tree_control" -> t3ATreeProbe["report"]|>];
t3ASaveMatchReport[];
If[!TrueQ[t3ATreeProbe["report"]["structurally_valid_for_requested_order"]],
  t3Stop["A(-4) TREE Matchete Match failed before the one-loop test; inspect A_joint_scalar_match_diagnostics.json and native_two_threshold.log",
    <|"diagnostic_path" -> t3AMatchDiagFile,
      "tree_return" -> t3ATreeProbe["report"]["return_status"],
      "source_theory" -> t3ACheckMode|>]
];
(* An optional raw tree WXF can be saved for reproducibility. It does not
   enter the one-loop computation; Match is re-run with the SAME input. *)
If[Environment["T3_A4_SAVE_TREE_CONTROL_WXF"] === "1",
  t3Save["A_joint_scalar_tree_control_raw_DIAGNOSTIC.wxf",
    t3ATreeProbe["result"]]];
Clear[t3ATreeProbe];
t3stage = "A-joint-heavy-scalar-soft-one-loop";
If[t3ACGTraceEnabled,
  t3ACGTraceStream = OpenWrite[t3ACGTracePath];
  WriteString[t3ACGTraceStream,
    "# DefineCG call sequence before a possible Matchete Abort[]\n"];
  Flush[t3ACGTraceStream]
];
t3AOneLoopAttempt = t3AMatchAttempt[t3AMatchInput, 1, 3600,
  "selected-full-theory-one-loop"];
If[Head[t3ACGTraceStream] === OutputStream,
  Close[t3ACGTraceStream]; t3ACGTraceStream = None];
t3AMatched = t3AOneLoopAttempt["result"];
t3AMatchReport = Join[t3AMatchReport,
  <|"one_loop" -> t3AOneLoopAttempt["report"],
    "DefineCG_trace_enabled" -> t3ACGTraceEnabled,
    "DefineCG_trace_count" -> t3ACGTraceCount,
    "DefineCG_trace_file" -> If[t3ACGTraceEnabled,
      "A_one_loop_DefineCG_calls_DIAGNOSTIC.txt", "NotRequested"],
    "matching_at_one_loop_succeeded" ->
      TrueQ[t3AOneLoopAttempt["report"]["structurally_valid_for_requested_order"]],
    "instrumentation_return_unevaluated" ->
      (t3AMatchStatusName[t3AMatched] ===
        "UnevaluatedTraceScan_InvalidInstrumentation")|>];
t3ASaveMatchReport[];
If[!TrueQ[t3AOneLoopAttempt["report"]["structurally_valid_for_requested_order"]],
  (* Diagnostic ONLY; same kernel, so a previous failed Matchete call could
     affect internals. This is not independent validation, and its output
     is NEVER used as physical C5. Default off to avoid expensive repeats. *)
  If[Environment["T3_A4_PROBE_WILSON_FREE_MATCH"] === "1" &&
      ValueQ[t3AScalarBase] && TrueQ[t3ACheckScalarBase],
    t3AWithoutW = t3AMatchAttempt[t3AScalarBase, 1, 1200,
      "diagnostic-only-source-scalar-base-no-W"];
    t3AMatchReport = Join[t3AMatchReport,
      <|"wilson_free_probe" -> t3AWithoutW["report"],
        "wilson_free_probe_same_kernel_after_failure" -> True,
        "wilson_free_probe_used_for_physics" -> False|>];
    t3ASaveMatchReport[];
    Clear[t3AWithoutW]
  ];
  (* V8: target the S1/S2 hard scalar loop directly, as a separate
     topology-only native Matchete calculation. CovariantLoop is a public
     API and must retain its own full EFT and gauge/CG validation. This
     does NOT replace failed Match, and does not include inherited F hard,
     EFT1 Wilson RG, other scalar traces, or finite renormalisation.

     If CovariantLoop succeeds, export its raw input and postprocessed
     Weinberg terms for independent numerical cross-checks. The outer
     physical-completion status remains blocked, by construction. *)
  If[Environment["T3_A4_PROBE_TWO_HEAVY_TRACE"] === "1",
    t3stage = "A-diagnostic-two-heavy-scalar-covariant-loop";
    ClearAll[t3ATraceSummary];
    t3ATraceMessages = {};
    t3ATraceStarted = AbsoluteTime[];
    t3ATraceCGEnabled = t3ACGTraceEnabled;
    t3ATraceCGPath = FileNameJoin[{t3out,
      "A_two_heavy_DefineCG_calls_DIAGNOSTIC.txt"}];
    If[t3ATraceCGEnabled,
      t3ACGTraceCount = 0;
      t3ACGTraceStream = OpenWrite[t3ATraceCGPath];
      WriteString[t3ACGTraceStream,
        "# DefineCG calls in the two-heavy CovariantLoop (diagnostic only)\n"];
      Flush[t3ACGTraceStream]
    ];
    Print["[A4-target] BEGIN CovariantLoop with both heavy scalars."];
    t3ATraceRaw = Block[{$MessageList = {}}, Module[{v},
      v = TimeConstrained[CheckAbort[UsingFrontEnd[
        If[t3ATraceCGEnabled,
          TraceScan[t3ACGTraceVisitor,
            CovariantLoop[t3AMatchInput, {NewScalar1, NewScalar2},
              EFTOrder -> 5, Mode -> Matching],
            HoldPattern[DefineCG[___]], TraceInternal -> True],
          CovariantLoop[t3AMatchInput, {NewScalar1, NewScalar2},
            EFTOrder -> 5, Mode -> Matching]
        ]], $Aborted], 1800, $TimedOut];
      t3ATraceMessages = ToString[#, InputForm] & /@ $MessageList;
      v]];
    t3ATraceCGCount = t3ACGTraceCount;
    If[Head[t3ACGTraceStream] === OutputStream,
      Close[t3ACGTraceStream]; t3ACGTraceStream = None];
    t3ATraceSummary = <|
      "DefineCG_trace_enabled" -> t3ATraceCGEnabled,
      "DefineCG_trace_count" -> If[t3ATraceCGEnabled,
        t3ATraceCGCount, 0],
      "DefineCG_trace_file" -> If[t3ATraceCGEnabled,
        "A_two_heavy_DefineCG_calls_DIAGNOSTIC.txt", "NotRequested"],
      "kind" -> "TwoHeavyScalarCovariantLoopOnly",
      "physical_C5_certified" -> False,
      "full_one_loop_Match_replaced" -> False,
      "one_loop_return" -> t3AMatchStatusName[t3ATraceRaw],
      "return_head" -> ToString[Head[t3ATraceRaw], InputForm],
      "elapsed_seconds" -> N[Max[0, AbsoluteTime[] - t3ATraceStarted]],
      "message_labels" -> Take[t3ATraceMessages, UpTo[15]],
      "explicit_hbar" -> !FreeQ[t3ATraceRaw, hbar],
      "remaining_heavy_fields" -> t3AMatchFieldNames[t3ATraceRaw],
      "raw_Weinberg_present" ->
        !FreeQ[t3ATraceRaw, GammaCC] && !FreeQ[t3ATraceRaw, Field[l, ___]]|>;
    If[!AtomQ[t3ATraceRaw] &&
        t3AMatchFieldNames[t3ATraceRaw] === {},
      t3Save["A_two_heavy_scalar_CovariantLoop_raw_DIAGNOSTIC.wxf",
        t3ATraceRaw];
      (* Run the same Matchete reduction sequence used in the validated
         D(-2) scalar stage, without adding any extra hbar. *)
      t3ATraceProcessed = t3ATraceRaw;
      t3ATracePostOK = True;
      Do[
        t3ATraceProcessed = TimeConstrained[
          CheckAbort[UsingFrontEnd[Switch[t3AStep,
            "GreensSimplify", GreensSimplify[t3ATraceProcessed],
            "EOMSimplify", EOMSimplify[t3ATraceProcessed],
            "EvaluateLoopFunctions", EvaluateLoopFunctions[t3ATraceProcessed],
            "ReplaceEffectiveCouplings", ReplaceEffectiveCouplings[t3ATraceProcessed],
            _, $Failed]], $Aborted], 900, $TimedOut];
        If[MemberQ[{$Failed, $Aborted, $TimedOut}, t3ATraceProcessed],
          t3ATracePostOK = False;
          t3ATraceSummary = Join[t3ATraceSummary,
            <|"postprocess_failed_at" -> t3AStep|>];
          Break[]
        ];
        t3Save["A_two_heavy_scalar_" <> t3AStep <>
          "_DIAGNOSTIC.wxf", t3ATraceProcessed],
        {t3AStep, {"GreensSimplify", "EOMSimplify",
          "EvaluateLoopFunctions", "ReplaceEffectiveCouplings"}}];
      If[t3ATracePostOK,
        t3ATraceLinear = If[FreeQ[t3ATraceProcessed, hbar],
          t3ATraceProcessed,
          hbar Coefficient[Expand[t3ATraceProcessed], hbar, 1]];
        t3ATraceWeinberg = CheckAbort[UsingFrontEnd[
          ExtractWeinbergCoefficient[t3ATraceLinear]], $Aborted];
        t3ATraceExtractOK = AssociationQ[t3ATraceWeinberg] &&
          Lookup[t3ATraceWeinberg, "Status", ""] === "Success" &&
          Lookup[t3ATraceWeinberg, "TermCount", 0] > 0;
        t3ATraceSummary = Join[t3ATraceSummary, <|
          "postprocessing_succeeded" -> True,
          "Weinberg_extraction_succeeded" -> t3ATraceExtractOK,
          "Weinberg_term_count" -> If[t3ATraceExtractOK,
             t3ATraceWeinberg["TermCount"], 0],
          "hbar_in_processed" -> !FreeQ[t3ATraceProcessed, hbar]|>];
        If[t3ATraceExtractOK,
          t3Save["A_two_heavy_scalar_Weinberg_DIAGNOSTIC.wxf",
            t3ATraceWeinberg["Coefficient"]];
          t3Save["A_two_heavy_scalar_Weinberg_DIAGNOSTIC.txt",
            ToString[t3ATraceWeinberg["Coefficient"], InputForm], "Text"]
        ]
      ]
    ];
    t3AMatchReport = Join[t3AMatchReport,
      <|"two_heavy_scalar_trace_DIAGNOSTIC_ONLY" -> t3ATraceSummary|>];
    t3ASaveMatchReport[];
    Print["[A4-target] END CovariantLoop; status=",
      t3ATraceSummary["one_loop_return"],
      "; output diagnostic-only, not physical C5."]
  ];
  t3Stop["A(-4) raw Matchete Match failed at one loop: " <>
      t3AOneLoopAttempt["report"]["return_status"] <>
      "; inspect A_joint_scalar_match_diagnostics.json and native_two_threshold.log",
    <|"diagnostic_path" -> t3AMatchDiagFile,
      "one_loop_return" -> t3AOneLoopAttempt["report"]["return_status"],
      "one_loop_return_inputform" ->
        t3AOneLoopAttempt["report"]["return_InputForm_if_atomic"],
      "source_theory" -> t3ACheckMode|>]
];
Clear[t3AOneLoopAttempt];
t3Save[t3AOutputPrefix <> "raw_matched.wxf", t3AMatched];

(* The same physical-threshold postprocessing which succeeded in D(-2). *)
t3AProcessed = t3AMatched;
Do[
  t3AProcessed = TimeConstrained[CheckAbort[UsingFrontEnd[
    Switch[step,
      "GreensSimplify", GreensSimplify[t3AProcessed],
      "EOMSimplify", EOMSimplify[t3AProcessed],
      "EvaluateLoopFunctions", EvaluateLoopFunctions[t3AProcessed],
      "ReplaceEffectiveCouplings", ReplaceEffectiveCouplings[t3AProcessed],
      _, $Failed]], $Aborted], 2400, $TimedOut];
  If[MemberQ[{$Failed, $Aborted, $TimedOut}, t3AProcessed],
    t3Stop["A(-4) scalar-soft postprocessing failed at " <> step]];
  t3Save[t3AOutputPrefix <> step <> ".wxf", t3AProcessed],
  {step, {"GreensSimplify", "EOMSimplify",
    "EvaluateLoopFunctions", "ReplaceEffectiveCouplings"}}
];
t3AOneLoop = hbar Coefficient[Expand[t3AProcessed], hbar, 1];
If[t3AOneLoop === 0,
  t3Stop["A(-4) physical scalar match has no O(hbar) soft contribution"]
];
t3AWeinberg = CheckAbort[
  ExtractWeinbergCoefficient[t3AOneLoop], $Aborted];
If[!AssociationQ[t3AWeinberg] ||
   Lookup[t3AWeinberg, "Status", ""] =!= "Success" ||
   Lookup[t3AWeinberg, "TermCount", 0] <= 0,
  t3Stop["A(-4) physical scalar-soft Weinberg extraction failed"]
];
t3Save[t3AOutputPrefix <> "unrenormalized.wxf", t3AWeinberg["Coefficient"]];
t3Save[t3AOutputPrefix <> "unrenormalized.txt",
  ToString[t3AWeinberg["Coefficient"], InputForm], "Text"];

t3stage = If[t3ASelectedExploratory,
  "A-exploratory-source-soft-no-boundary-equivalence",
  "A-physical-soft-success-hard-inheritance-pending"];
t3Write[
  If[t3ASelectedExploratory,
    "A_ExploratoryHermitianScalarSoft_NoEFTBoundaryCertificate",
    "A_PhysicalScalarSoftOnly_HardInheritanceAndWilsonRGPending"],
  If[t3ASelectedExploratory,
    "Diagnostic scalar match succeeded for source-rebuilt Hermitian Wilson EFT, but exact equivalence to F boundary is unproven. NOT physical C5.",
    "Three-flavour joint heavy scalar-soft matching succeeded; F hard, EFT1 Wilson RG and finite renormalisation remain absent. NOT physical C5."],
  <|"first_threshold_reused" -> TrueQ[t3reused],
    "A_joint_scalar_match_performed" -> True,
    "A_mass_power_counting_reversibly_promoted" -> True,
    "A_physical_scalar_soft_Weinberg_terms" ->
      If[t3ASelectedExploratory, 0, t3AWeinberg["TermCount"]],
    "A_exploratory_scalar_soft_Weinberg_terms" ->
      If[t3ASelectedExploratory, t3AWeinberg["TermCount"], 0],
    "F_boundary_operator_equivalence_certified" ->
      !TrueQ[t3ASelectedExploratory] && TrueQ[t3AFullGraphEquivalent],
    "exploratory_source_theory_used" -> TrueQ[t3ASelectedExploratory],
    "A_hard_source_WXF" -> "A_full_F_hard_one_loop_preserved.wxf",
    "A_soft_source_WXF" -> t3AOutputPrefix <> "unrenormalized.wxf"|>];
(* Exit 5 intentionally: scalar soft is not the requested completed C5. *)
Exit[5];
