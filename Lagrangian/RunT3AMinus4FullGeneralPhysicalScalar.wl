(* A(-4): native physical scalar threshold, local-only continuation.
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
  "CheckLagrangian_result" -> ToString[t3ACheck, InputForm],
  "CheckLagrangian_passed" -> TrueQ[t3ACheck],
  "accepted_canonicalisation" -> t3ACheckMode,
  "mass_promotion_reversible" -> t3AMassValid,
  "projection_used_for_physical_match" -> False,
  "matched_C5_physically_certified" -> False|>, "RawJSON"];
If[!TrueQ[t3ACheck],
  t3Stop["A(-4) CheckLagrangian remains false after exact-equivalence tests; inspect A_joint_scalar_lagrangian_audit.json and native_two_threshold.log",
    <|"check_result" -> ToString[t3ACheck, InputForm],
      "mass_audit_path" -> t3AAuditPath,
      "lagrangian_audit_path" -> FileNameJoin[{t3out,
        "A_joint_scalar_lagrangian_audit.json"}]|>]
];
(* The only approved replacement has passed BOTH CheckLagrangian and
   term-by-term alpha equivalence. The source EFT remains independently saved. *)
t3AMatchInput = t3ASelectedInput;
t3Save["A_EFT1_joint_heavy_scalar_checked_input.wxf", t3AMatchInput];

(* One common scalar threshold (S1,S2).  This is scalar-soft only.
   The F hard 1-loop EFT remains separately saved and must later be
   propagated through a TREE scalar match, not re-loop-matched. *)
t3stage = "A-joint-heavy-scalar-soft-one-loop";
t3AMatched = TimeConstrained[
  CheckAbort[UsingFrontEnd[
    Match[t3AMatchInput, EFTOrder -> 5, LoopOrder -> 1]], $Aborted],
  3600, $TimedOut];
If[MemberQ[{$Failed, $Aborted, $TimedOut}, t3AMatched] ||
   AtomQ[t3AMatched] || FreeQ[t3AMatched, hbar] ||
   !FreeQ[t3AMatched,
     HoldPattern[Field[NewScalar1 | NewScalar2, ___]]],
  t3Stop["A(-4) joint scalar-soft Matchete result invalid or heavy fields remain",
    <|"result_head" -> ToString[Head[t3AMatched], InputForm]|>]
];
t3Save["A_scalar_soft_only_raw_matched.wxf", t3AMatched];

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
  t3Save["A_scalar_soft_only_" <> step <> ".wxf", t3AProcessed],
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
t3Save["A_scalar_soft_only_unrenormalized.wxf", t3AWeinberg["Coefficient"]];
t3Save["A_scalar_soft_only_unrenormalized.txt",
  ToString[t3AWeinberg["Coefficient"], InputForm], "Text"];

t3stage = "A-physical-soft-success-hard-inheritance-pending";
t3Write["A_PhysicalScalarSoftOnly_HardInheritanceAndWilsonRGPending",
  "Three-flavour joint heavy scalar-soft one-loop matching succeeded; F hard, EFT1 Wilson RG and finite renormalisation remain absent. NOT physical C5.",
  <|"first_threshold_reused" -> TrueQ[t3reused],
    "A_joint_scalar_match_performed" -> True,
    "A_mass_power_counting_reversibly_promoted" -> True,
    "A_physical_scalar_soft_Weinberg_terms" -> t3AWeinberg["TermCount"],
    "A_hard_source_WXF" -> "A_full_F_hard_one_loop_preserved.wxf",
    "A_soft_source_WXF" -> "A_scalar_soft_only_unrenormalized.wxf"|>];
(* Exit 5 intentionally: scalar soft is not the requested completed C5. *)
Exit[5];
