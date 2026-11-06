(* General D(-2): second threshold, integrating all physical scalars.
   Inputs are the *actual* first-threshold WXF continuations in the
   real-S2 basis.  The direct LLSS -> LLHH running insertion is carried
   separately: it must never be included in a one-loop scalar Match[].

   Usage:
     wolframscript -file Lagrangian/RunGeneralDMinus2ScalarThreshold.wl \
       FIRST_THRESHOLD_OUTPUT_DIR DIRECT_RUNNING_INSERTION.wl [profile:MODE] [SECONDS]

   profile:inspect             operator-sector and scalar-power inspection
   profile:sm-only             L_SM with full registry (no heavy free terms)
   profile:s2r-only            free heavy real singlet S2R
   profile:s1-only             free heavy complex SU(2) triplet S1
   profile:heavy-only          all three free heavy scalar multiplets, no L_SM
   profile:free                SM plus free heavy scalar multiplets
   profile:scalar-standard     add existing non-extra scalar quartics
   profile:scalar-general      add all five extra scalar quartics
   profile:scalar-with-y2      add the y1-y2 LLSS operator pair and HC
   profile:scalar-with-y2Conj  add the y1-y2Conj LLSS pair and HC
   profile:bsm-free            free heavy fields WITHOUT stand-alone L_SM
   profile:bsm-standard        free heavy fields + ordinary scalar quartics
   profile:bsm-general         above + all five additional quartics
   profile:bsm-with-y2         above + y1-y2 LLSS operators + HC
   profile:bsm-with-y2Conj     above + y1-y2Conj LLSS operators + HC
   profile:bsm-with-both       above + both independent LLSS sectors + HC
   profile:full-tree           original full tree Match (opt-in)
   The former 'renormalisable', 'with_y2', and 'with_y2Conj' profile
   labels remain diagnostic legacy aliases and are not run by default.

   Each probe starts a fresh kernel.  These reduced inputs are diagnostic
   experiments and must NOT replace the general matching result.

   Fixed order: [A] M0[L0], [B] M1[L0]-M0[L0],
                [C] M0[L0+hbar L1]-M0[L0] projected to hbar^1.
   Masses MS1, T3S2RMass, T3S2IMass stay independent Matchete labels;
   their split-mass relations are metadata, NOT imposed numerically.
*)
ClearAll["Global`*"];
args = Rest[$ScriptCommandLine];
If[!MemberQ[{2, 4}, Length[args]],
  Print["Usage: OUTPUT_DIR RUNNING_INSERTION.wl [profile:MODE SECONDS]"];
  Exit[2]
];
profileMode = Length[args] === 4;
profileProbe = If[profileMode, StringReplace[args[[3]], "profile:" -> ""], ""];
profileAllowed = {
  "inspect", "sm-only", "s2r-only", "s1-only", "heavy-only",
  "free", "scalar-standard", "scalar-general",
  "scalar-with-y2", "scalar-with-y2Conj",
  "bsm-free", "bsm-standard", "bsm-general", "bsm-with-y2",
  "bsm-with-y2Conj", "bsm-with-both", "full-tree",
  "renormalisable", "with_y2", "with_y2Conj"
};
profileLimit = If[profileMode,
  Quiet @ Check[ToExpression[args[[4]], InputForm], $Failed], 0];
If[profileMode && (!StringStartsQ[args[[3]], "profile:"] ||
   !MemberQ[profileAllowed, profileProbe] ||
   !IntegerQ[profileLimit] || profileLimit < 5 || profileLimit > 3600),
  Print["ERROR: invalid tree profile mode or timeout"]; Exit[2]
];
outputDir = ExpandFileName[args[[1]]];
insertionPath = ExpandFileName[args[[2]]];
dataDir = FileNameJoin[{outputDir, "data"}];
stageDir = FileNameJoin[{dataDir, "general_d_m2_scalar_threshold"}];
CreateDirectory[stageDir, CreateIntermediateDirectories -> True];
profileDir = FileNameJoin[{stageDir, "profile_tree"}];
If[profileMode,
  CreateDirectory[profileDir, CreateIntermediateDirectories -> True]
];
statusPath = If[profileMode,
  FileNameJoin[{profileDir, "wolfram_profile_status_" <> profileProbe <> ".json"}],
  FileNameJoin[{stageDir, "wolfram_scalar_threshold_status.json"}]
];
scriptDir = DirectoryName @ ExpandFileName[$InputFileName];

ClearAll[D2SaveStatus, D2Stop, D2Text, D2Success, D2CouplingPresent];
D2SaveStatus[status_String, details_: <||>] := Export[statusPath,
  Join[<|
    "Status" -> status, "Model" -> "D_-2",
    "IntegratedFields" -> {"NewScalar1", "T3S2R", "T3S2I"},
    "PhysicalMassRelationsImposed" -> False,
    "FullFlavorHardC5Verified" -> False,
    "DirectRunningIncludedInLoopMatch" -> False,
    "RGBetaMappingCertified" -> False
  |>, details], "RawJSON"];
D2Stop[code_Integer, label_String, details_: <||>] := (
  Print["ERROR: ", label];
  D2SaveStatus["Failed", Join[<|"Failure" -> label,
    "ExitCode" -> code|>, details]];
  Exit[code]
);
D2Text[expression_] := ToString[expression, InputForm];
D2Success[result_] := AssociationQ[result] &&
  Lookup[result, "Status", ""] === "Success";
D2SaveStatus["Started"];

If[Environment["T3_EXPERIMENTAL_GENERAL_QUADRATICS"] =!= "1" ||
   Environment["T3_EXPERIMENTAL_GENERAL_FERMIONS"] =!= "1",
  D2Stop[3, "Both general interaction flags must be enabled"]
];
If[!FileExistsQ[insertionPath], D2Stop[4, "Running insertion missing"]];
fileNames = {
  "eft1_after_F_realbasis_tree.wxf",
  "eft1_after_F_realbasis_loop.wxf",
  "eft1_after_F_realbasis_transition_tree.wxf",
  "eft1_after_F_realbasis_transition_full.wxf"
};
paths = FileNameJoin[{dataDir, #}] & /@ fileNames;
If[!TrueQ[And @@ (FileExistsQ /@ paths)],
  D2Stop[5, "Missing real-basis EFT1 WXF continuation",
    <|"MissingFiles" -> Pick[fileNames, Not /@ (FileExistsQ /@ paths)]|>]
];

Print["[D(-2) scalar threshold] Loading Matchete in a fresh kernel..."];
loadResult = CheckAbort[UsingFrontEnd[Needs["Matchete`"]; True], $Aborted];
If[loadResult =!= True, D2Stop[6, "Matchete failed to load"]];
Get[FileNameJoin[{scriptDir, "PhysicsNotation.wl"}]];
Get[FileNameJoin[{scriptDir, "model", "T3ModelCatalog.wl"}]];
Get[FileNameJoin[{scriptDir, "model", "LagrangianBuilder.wl"}]];
Get[FileNameJoin[{scriptDir, "RunMatching.wl"}]];
model = T3ModelFromClass["D", -2];
If[model === $Failed, D2Stop[7, "Cannot reconstruct D(-2) catalog"]];

(* Rebuild the exact unrestricted UV registry first; this is *not* another
   UV matching.  It redefines all independently generated CGs, including
   the five non-topology quartics and the extra Yukawa coupling. *)
Print["[D(-2) scalar threshold] Restoring full general UV CG/coupling registry..."];
If[SetT3HeavyFields[{"F"}] === $Failed,
  D2Stop[8, "Heavy-field setup failed"]];
uvBuild = CheckAbort[UsingFrontEnd[BuildT3Lagrangian[model]], $Aborted];
If[!D2Success[uvBuild] ||
   Lookup[uvBuild, "AdditionalQuarticsIncluded", -1] =!= 5 ||
   Lookup[uvBuild, "AdditionalQuadraticsIncluded", -1] =!= 1 ||
   Lookup[uvBuild, "AdditionalYukawasIncluded", -1] =!= 1,
  D2Stop[9, "General UV registry reconstruction failed",
    <|"BuildStatus" -> If[AssociationQ[uvBuild],
      Lookup[uvBuild,"Status","Unknown"], D2Text[uvBuild]]|>]
];

(* F is integrated out.  S1 and both real S2 components are all Heavy.
   Use native Matchete field definitions with the same symbols as EFT1 WXF.
   RemoveField returns Null on success in this Matchete version. *)
RemoveField[NewFermion];
RemoveField[NewScalar1];
RemoveField[NewScalar2];
If[SetT3HeavyFields[{"S1"}] === $Failed,
  D2Stop[10, "S1 heavy-field classification failed"]];
If[DefineT3FieldObjects[model, False, {"S1"}] === $Failed,
  D2Stop[11, "S1 field registration failed"]];
DefineField[T3S2R, Scalar, SelfConjugate -> True,
  Mass -> {Heavy, T3S2RMass}];
DefineField[T3S2I, Scalar, SelfConjugate -> True,
  Mass -> {Heavy, T3S2IMass}];
PreserveIntegratedT3MassCouplings[{"S1"}];
If[!KeyExistsQ[GetCouplings[], T3S2Phase],
  DefineCoupling[T3S2Phase, SelfConjugate -> False]
];

fieldText = D2Text[GetFields[]];
If[!TrueQ[And @@ (StringContainsQ[fieldText, #] & /@
       {"NewScalar1", "T3S2R", "T3S2I"})] ||
   StringContainsQ[fieldText, "NewFermion"] ||
   StringContainsQ[fieldText, "NewScalar2"] ||
   !TrueQ[And @@ (StringContainsQ[fieldText,
       # <> " -> <|Type -> Scalar"] & /@ {"T3S2R", "T3S2I"})],
  D2Stop[12, "Incorrect real-heavy-scalar field registry",
    <|"RegistryPreview" -> StringTake[fieldText, UpTo[2500]]|>]
];
fieldRegistry = GetFields[];
If[!AllTrue[{NewScalar1, T3S2R, T3S2I},
      TrueQ[Lookup[Lookup[fieldRegistry, #, <||>], Heavy, False]] &],
  D2Stop[12, "Not all three physical scalars are registered as Heavy"]
];
Export[FileNameJoin[{stageDir, "real_scalar_registry_audit.json"}], <|
  "Status" -> "Success", "OldFermionAbsent" -> True,
  "OldComplexS2Absent" -> True,
  "PhysicalScalarFields" -> {"NewScalar1", "T3S2R", "T3S2I"},
  "HeavyMassSymbols" -> {"MS1", "T3S2RMass", "T3S2IMass"},
  "ExtraQuartics" -> 5, "ExtraYukawas" -> 1
|>, "RawJSON"];

(* Rehydrate WXF *after* fields and couplings are registered.  Re-creating
   the old complex S2 would invalidate both the mass basis and matching. *)
Print["[D(-2) scalar threshold] Importing four physical EFT1 continuations..."];
wxf = Quiet @ Check[Import[#, "WXF"], $Failed] & /@ paths;
If[MemberQ[wxf, $Failed] || Length[wxf] =!= 4,
  D2Stop[13, "WXF continuation import failed"]];
{currentTree, currentLoop, transitionTree, transitionFull} = wxf;
ClearAll[D2Rehydrate];
D2Rehydrate[expr_] := Quiet @ Check[
  ToExpression[ToString[InputForm[Unevaluated[expr]]], InputForm], $Failed
];
transitionTree = D2Rehydrate[transitionTree];
transitionFull = D2Rehydrate[transitionFull];
If[MemberQ[{transitionTree, transitionFull}, $Failed],
  D2Stop[14, "Transition rehydration failed"]];

(* Inherited wave-function/mass effective couplings may be serialized in
   the canonical EFT with no registry metadata.  Known Hermitian scalar
   quadratic effective couplings are restored as in RunThresholdStage.wl.
   Unknown private effective-coupling names cause an explicit diagnostic,
   rather than silently guessing their conjugation properties. *)
Scan[Function[sym,
  If[!KeyExistsQ[GetCouplings[], sym],
    Apply[DefineCoupling, {sym, SelfConjugate -> True}]
  ]
], {CNewScalar12, CT3S2R2, CT3S2I2}];

(* A mass basis is not sufficient if fresh Matchete has forgotten a CG.
   Compare the actual serialized EFT with the live CG registry.  Regenerate
   missing seed CGs using their exact original tensors (no invented CG). *)
seedPaths = {
  FileNameJoin[{dataDir, "eft1_after_F_realbasis_wilson_tensor_seed.json"}],
  FileNameJoin[{dataDir, "eft1_after_F_realbasis_quartic_tensor_seed.json"}]
};
If[!TrueQ[And @@ (FileExistsQ /@ seedPaths)],
  D2Stop[15, "Same-kernel physical-basis CG seed files missing"]];
seeds = Import[#, "RawJSON"] & /@ seedPaths;
If[!AllTrue[seeds, AssociationQ] ||
   !AllTrue[seeds, Lookup[#, "Status", ""] === "Success" &],
  D2Stop[16, "Invalid physical-basis CG seeds"]];
registry = DeleteDuplicatesBy[
  Flatten[Lookup[seeds, "CGRegistry", {}], 1], Lookup[#, "Name", ""] &
];
restored = 0;
Do[
  item = registry[[k]];
  cgName = ToExpression[item["Name"], InputForm];
  existing = Quiet @ Check[
    Matchete`PackageScope`$CGproperties[cgName, Indices], $Failed];
  If[existing === $Failed || MissingQ[existing],
    cgReps = ToExpression[item["RepsInputForm"], InputForm];
    cgTensor = ToExpression[item["TensorInputForm"], InputForm];
    If[CheckAbort[DefineCG[cgName, cgReps, cgTensor]; True,
       $Aborted] =!= True,
      D2Stop[17, "CG restoration aborted",
        <|"CGName" -> item["Name"]|>]
    ];
    restored++
  ],
  {k, Length[registry]}
];
Print["[D(-2) scalar threshold] Restored ", restored,
  " extra EFT1 CG tensor(s)."];
(* Unknown CGs would invalidate Matchete's contraction rules. *)
usedCGNames = DeleteDuplicates @ Cases[
  HoldComplete[transitionTree, transitionFull],
  CG[name_, ___] :> Replace[name, Bar[x_] :> x], Infinity];
missingCGNames = Select[usedCGNames, Function[name,
  Module[{result = Quiet @ Check[
      Matchete`PackageScope`$CGproperties[name, Indices], $Failed]},
    result === $Failed || MissingQ[result]
  ]
]];
If[missingCGNames =!= {},
  D2Stop[17, "EFT1 references CG tensors absent in fresh kernel",
    <|"MissingCGNames" -> (D2Text /@ missingCGNames)|>]
];

(* Keep the direct SM Weinberg insertion strictly outside Match[]. *)
Print["[D(-2) scalar threshold] Loading direct full-flavour running..."];
Get[insertionPath];
If[!ValueQ[EFT1D2WilsonRunningHeavyInsertion] ||
   !TrueQ[EFT1D2WilsonRunningHeavyInsertion === 0] ||
   DownValues[EFT1D2WilsonRunningWeinbergCoefficient] === {},
  D2Stop[18, "Direct-running insertion has incompatible contract"]
];
p = Index[T3D2ExternalFlavorP, Flavor];
q = Index[T3D2ExternalFlavorQ, Flavor];
directC5 = EFT1D2WilsonRunningWeinbergCoefficient[p, q];
If[FreeQ[directC5, Log] ||
   FreeQ[directC5, y1] &&
   !StringContainsQ[D2Text[directC5], "y1"],
  D2Stop[19, "Direct running is not an explicit logarithmic flavour tensor"]
];
Export[FileNameJoin[{stageDir, "direct_running_c5_inputform.txt"}],
  D2Text[directC5], "Text"];

(* Gauge kinetic shifts are matched as coupling renormalization.  Mirror the
   verified fixed-order fresh-kernel threshold code; do not regard an O(hbar)
   gauge kinetic correction as a new noncanonical operator. *)
treeExpanded = Quiet @ Check[
  ReplaceEffectiveCouplings[transitionTree, Superleading -> True], $Failed];
fullExpanded = Quiet @ Check[
  ReplaceEffectiveCouplings[transitionFull, Superleading -> True], $Failed];
If[MemberQ[{treeExpanded, fullExpanded}, $Failed],
  D2Stop[20, "Cannot expand inherited effective couplings"]];
ClearAll[D2GaugeClass];
D2GaugeClass[lag_, fields_List] := Quiet @ Check[
  SelectOperatorClass[lag, fields, 0], 0];
gaugeTree = Total[D2GaugeClass[treeExpanded, #] & /@
  {{B, B}, {W, W}, {G, G}}];
gaugeFull = Total[D2GaugeClass[fullExpanded, #] & /@
  {{B, B}, {W, W}, {G, G}}];
gaugeShift = Expand[gaugeFull - gaugeTree];
transitionFullCanonical = Expand[fullExpanded - gaugeShift];
Export[FileNameJoin[{stageDir, "inherited_gauge_shift_inputform.txt"}],
  D2Text[gaugeShift], "Text"];

(* -----------------------------------------------------------------------
   Incremental tree-level profiler for physical D(-2) EFT1.

   Important distinction: removing GammaCC LLSS terms from the complete EFT1
   does NOT construct a minimal renormalisable scalar sector.  The original
   baseline timed out, so construct new inputs from known canonical sectors:
      free:          exact registered SM + heavy scalar FreeLag
      scalar-standard: free + non-extra scalar quartic interactions
      scalar-general:  standard + all five general extra quartics
      scalar-with-y2:  general + y1-y2 LLSS Wilson terms + h.c.
      scalar-with-y2Conj: general + y1-y2Conj LLSS terms + h.c.

   This is profiling ONLY.  The authoritative [A]/[B]/[C] below still uses
   the complete original canonical transition EFT without removing terms.
   ----------------------------------------------------------------------- *)
If[profileMode,
  Print["[profile] Preparing independently classified physical EFT1 sectors..."];
  profileStarted = AbsoluteTime[];
  profileCompactTree = Expand[transitionTree];
  profileTree = Expand[currentTree];
  profileTransitionExplicit = Quiet @ CheckAbort[
    ReplaceEffectiveCouplings[transitionTree], $Aborted];
  profileExactRepresentation = If[
    MemberQ[{$Failed, $Aborted}, profileTransitionExplicit], False,
    TrueQ[Expand[profileTree - Expand[profileTransitionExplicit]] === 0]
  ];
  profileTerms = If[Head[profileTree] === Plus,
    List @@ profileTree, {profileTree}];
  profileTexts = D2Text /@ profileTerms;
  profileBSMNames = {"NewScalar1", "T3S2R", "T3S2I"};
  profileExtraNames = {
    "lambdaT3ExtraQ000004I1", "lambdaT3ExtraQ000013I1",
    "lambdaT3ExtraQ001102I1", "lambdaT3ExtraQ020101I1",
    "lambdaT3ExtraQ110002I1"
  };

  (* Count powers, not just textual Field[...] occurrences: R^4 is four
     scalar legs although it contains only one Field expression.  An
     unrecognised field representation is a hard inspection failure. *)
  ClearAll[D2ProfileScalarDegree];
  D2ProfileScalarDegree[term_] := Module[{objects, powers, barredPowers},
    objects = Cases[HoldComplete[term],
      HoldPattern[Field[_, Scalar, ___]] :> 1, Infinity];
    powers = Cases[HoldComplete[term],
      HoldPattern[Power[Field[_, Scalar, ___], n_Integer?Positive]] :>
        n - 1, Infinity];
    barredPowers = Cases[HoldComplete[term],
      HoldPattern[Power[Bar[Field[_, Scalar, ___]], n_Integer?Positive]] :>
        n - 1, Infinity];
    Total[objects] + Total[powers] + Total[barredPowers]
  ];
  profileBSMQ[text_String] := Or @@ (
    StringContainsQ[text, "Field[" <> # <> ","] & /@
      profileBSMNames
  );
  profileWilsonQ[text_String] :=
    StringContainsQ[text, "GammaCC"] && profileBSMQ[text];
  profileWPositions = Select[Range[Length[profileTerms]],
    profileWilsonQ[profileTexts[[#]]] &];
  profileRenPositions = Complement[Range[Length[profileTerms]],
    profileWPositions];
  profileY2Positions = Select[profileWPositions,
    StringContainsQ[profileTexts[[#]], "Coupling[y2,"] &];
  profileY2ConjPositions = Select[profileWPositions,
    StringContainsQ[profileTexts[[#]], "Coupling[y2Conj,"] &];
  profileUnknown = Complement[profileWPositions,
    Join[profileY2Positions, profileY2ConjPositions]];
  profileOverlap = Intersection[profileY2Positions, profileY2ConjPositions];
  profileRen = Total[profileTerms[[profileRenPositions]]];
  profileY2 = Total[profileTerms[[profileY2Positions]]];
  profileY2Conj = Total[profileTerms[[profileY2ConjPositions]]];
  profileExact = TrueQ[
    Expand[profileRen + Total[profileTerms[[profileWPositions]]] -
      profileTree] === 0];

  (* Scalar contacts only: exclude fermions, derivatives, field strengths,
     Gamma/spinor operators, and every Higgs-only term already in L_SM.
     All selected contacts must contain exactly four scalar field powers.
     The tree matching may also contain unrelated higher-dimensional terms;
     these are intentionally absent from this isolated renormalisable probe. *)
  profileForbidden = {
    ", Fermion,", "FieldStrength[", "GammaCC", "GammaM[",
    "NCM[", "CovariantD[", "Derivative[", "CD[", "hbar"
  };
  profileScalarCandidates = Select[profileRenPositions,
    Function[pos,
      Module[{text = profileTexts[[pos]]},
        profileBSMQ[text] &&
        StringContainsQ[text, ", Scalar,"] &&
        !AnyTrue[profileForbidden, StringContainsQ[text, #] &]
      ]
    ]
  ];
  profileCandidateDegrees = D2ProfileScalarDegree /@
    profileTerms[[profileScalarCandidates]];
  profileQuarticPositions = Select[profileScalarCandidates,
    D2ProfileScalarDegree[profileTerms[[#]]] === 4 &];
  profileQuarticUnresolved = Select[profileScalarCandidates,
    D2ProfileScalarDegree[profileTerms[[#]]] === 0 &];
  profileExtraQuarticPositions = Select[profileQuarticPositions,
    StringContainsQ[profileTexts[[#]], "lambdaT3ExtraQ"] &];
  profileStandardQuarticPositions = Complement[
    profileQuarticPositions, profileExtraQuarticPositions];
  profileMissingExtra = Select[profileExtraNames,
    Function[coupling,
      !AnyTrue[profileExtraQuarticPositions,
        Function[pos, StringContainsQ[profileTexts[[pos]], coupling]]]
    ]
  ];
  profileQuarticHbar = Select[profileQuarticPositions,
    !FreeQ[profileTerms[[#]], hbar] &];

  (* The full-general UV construction supplied L_SM; the newly registered
     physical heavy scalar objects supply canonical free mass/kinetic terms.
     Unlike the old 'non-GammaCC' probe, this never accidentally includes
     unclassified higher-dimensional EFT operators. *)
  (* Smoke probes use the same post-F physical registry, but deliberately
     separate the SM and heavy-field free kinetic/mass sectors.  A successful
     isolated Match is not a substitute for the full physical EFT1 Match. *)
  profileSMOnly = CheckAbort[
    RelabelIndices[Expand[uvBuild["LSM"]]], $Aborted];
  profileS2ROnly = CheckAbort[
    RelabelIndices[Expand[FreeLag[T3S2R]]], $Aborted];
  profileS1Only = CheckAbort[
    RelabelIndices[Expand[FreeLag[NewScalar1]]], $Aborted];
  profileHeavyOnly = CheckAbort[
    RelabelIndices[Expand[FreeLag[NewScalar1, T3S2R, T3S2I]]],
    $Aborted];
  profileFree = CheckAbort[
    RelabelIndices[Expand[
      uvBuild["LSM"] + FreeLag[NewScalar1, T3S2R, T3S2I]
    ]], $Aborted
  ];
  profileStandard = If[profileFree === $Aborted, $Aborted,
    RelabelIndices[Expand[
      profileFree + Total[profileTerms[[profileStandardQuarticPositions]]]
    ]]
  ];
  profileGeneral = If[profileStandard === $Aborted, $Aborted,
    RelabelIndices[Expand[
      profileStandard + Total[profileTerms[[profileExtraQuarticPositions]]]
    ]]
  ];

  (* Diagnostic comparison with the same heavy-scalar free and quartic
     interactions but WITHOUT stand-alone L_SM.  The SM gauge/light fields
     that occur inside covariant derivatives and BSM operators are retained.
     Do not interpret this reduced input as the authoritative EFT1.  Pure
     SM terms must be restored in the final SMEFT, and the one-loop heavy
     threshold must not be assumed equivalent without its own validation. *)
  profileBSMFree = profileHeavyOnly;
  profileBSMStandard = RelabelIndices[Expand[
    profileBSMFree + Total[profileTerms[[profileStandardQuarticPositions]]]
  ]];
  profileBSMGeneral = RelabelIndices[Expand[
    profileBSMStandard + Total[profileTerms[[profileExtraQuarticPositions]]]
  ]];
  profileBSMWithY2 = RelabelIndices[Expand[profileBSMGeneral + profileY2]];
  profileBSMWithY2Conj = RelabelIndices[Expand[
    profileBSMGeneral + profileY2Conj]];
  profileBSMWithBoth = RelabelIndices[Expand[
    profileBSMGeneral + profileY2 + profileY2Conj]];

  profileAudit = <|
    "Status" -> "Inspected", "Model" -> "D_-2",
    "FullGeneralInteractionFlags" -> True,
    "RealHeavyScalars" -> {"NewScalar1", "T3S2R", "T3S2I"},
    "ScalarMassLabels" -> {"MS1", "T3S2RMass", "T3S2IMass"},
    "ProfileStartUnix" -> profileStarted - 2208988800,
    "TreeByteCount" -> ByteCount[profileTree],
    "TreeLeafCount" -> LeafCount[profileTree],
    "TreeInputSource" -> "eft1_after_F_realbasis_tree.wxf",
    "CompactTransitionHead" -> D2Text[Head[profileCompactTree]],
    "ExplicitMatchesExpandedCanonical" -> profileExactRepresentation,
    "ExpandedTermCount" -> Length[profileTerms],
    "GammaCCScalarWilsonTermCount" -> Length[profileWPositions],
    "Y1Y2WilsonTermCount" -> Length[profileY2Positions],
    "Y1Y2ConjWilsonTermCount" -> Length[profileY2ConjPositions],
    "UnclassifiedWilsonTermCount" -> Length[profileUnknown],
    "OverlappingWilsonTermCount" -> Length[profileOverlap],
    "ExactTermPartition" -> profileExact,
    "ScalarPotentialCandidateCount" -> Length[profileScalarCandidates],
    "ScalarPotentialCandidateDegrees" -> profileCandidateDegrees,
    "SelectedScalarQuarticTermCount" -> Length[profileQuarticPositions],
    "StandardScalarQuarticTermCount" ->
      Length[profileStandardQuarticPositions],
    "ExtraScalarQuarticTermCount" ->
      Length[profileExtraQuarticPositions],
    "MissingGeneralQuarticCouplings" -> profileMissingExtra,
    "UnresolvedScalarFieldDegreeCount" -> Length[profileQuarticUnresolved],
    "NonzeroHbarQuarticTerms" -> Length[profileQuarticHbar],
    "SMOnlyBytes" -> ByteCount[profileSMOnly],
    "S2RFreeOnlyBytes" -> ByteCount[profileS2ROnly],
    "S1FreeOnlyBytes" -> ByteCount[profileS1Only],
    "AllHeavyFreeOnlyBytes" -> ByteCount[profileHeavyOnly],
    "SMPlusFreeHeavyScalarBytes" -> ByteCount[profileFree],
    "StandardScalarSectorBytes" -> ByteCount[profileStandard],
    "GeneralScalarSectorBytes" -> ByteCount[profileGeneral],
    "BSMFreeSectorBytes" -> ByteCount[profileBSMFree],
    "BSMStandardSectorBytes" -> ByteCount[profileBSMStandard],
    "BSMGeneralSectorBytes" -> ByteCount[profileBSMGeneral],
    "BSMGeneralSectorContainsExplicitLSM" -> False,
    "BSMWithY2SectorBytes" -> ByteCount[profileBSMWithY2],
    "BSMWithY2ConjSectorBytes" -> ByteCount[profileBSMWithY2Conj],
    "SMFreeScalarSectorSource" ->
      "UV general L_SM plus registered heavy S1, R, I FreeLag",
    "QuarticSource" -> "complete explicit matched EFT1 tree; exact scalar power degree=4",
    "SMHiggsSelfQuarticNotDuplicated" -> True,
    "AllGeneralQuarticsIncludedInDiagnosticGeneralSector" ->
      (profileMissingExtra === {}),
    "DirectWeinbergInsertedIntoMatch" -> False,
    "AuthoritativeFullTreeUntouched" -> True,
    "BuildAndInventorySeconds" -> AbsoluteTime[] - profileStarted
  |>;
  Export[FileNameJoin[{profileDir, "operator_sector_audit.json"}],
    profileAudit, "RawJSON"];
  Print["[profile] EFT1 terms: ", Length[profileTerms],
    "; LLSS: ", Length[profileWPositions],
    "; scalar quartics: ", Length[profileQuarticPositions],
    " (general: ", Length[profileExtraQuarticPositions], ")"];
  Print["[profile] Explicit/canonical tree equivalent: ",
    profileExactRepresentation, "; missing general quartics: ",
    profileMissingExtra];

  If[!profileExactRepresentation || !profileExact ||
     Length[profileWPositions] === 0 || Length[profileUnknown] > 0 ||
     Length[profileOverlap] > 0 || Length[profileY2Positions] === 0 ||
     Length[profileY2ConjPositions] === 0,
    D2Stop[31, "Cannot safely partition physical EFT1 tree for profiling",
      <|"Audit" -> profileAudit|>]
  ];
  If[Length[profileQuarticPositions] === 0 ||
     Length[profileExtraQuarticPositions] === 0 ||
     profileMissingExtra =!= {} || Length[profileQuarticUnresolved] > 0 ||
     Length[profileQuarticHbar] > 0 ||
     MemberQ[{$Failed, $Aborted}, profileSMOnly] ||
     MemberQ[{$Failed, $Aborted}, profileS2ROnly] ||
     MemberQ[{$Failed, $Aborted}, profileS1Only] ||
     MemberQ[{$Failed, $Aborted}, profileHeavyOnly] ||
     MemberQ[{$Failed, $Aborted}, profileFree] ||
     MemberQ[{$Failed, $Aborted}, profileStandard] ||
     MemberQ[{$Failed, $Aborted}, profileGeneral] ||
     MemberQ[{$Failed, $Aborted}, profileBSMFree] ||
     MemberQ[{$Failed, $Aborted}, profileBSMStandard] ||
     MemberQ[{$Failed, $Aborted}, profileBSMGeneral] ||
     MemberQ[{$Failed, $Aborted}, profileBSMWithY2] ||
     MemberQ[{$Failed, $Aborted}, profileBSMWithY2Conj] ||
     MemberQ[{$Failed, $Aborted}, profileBSMWithBoth],
    D2Stop[34, "Cannot certify incremental scalar-sector inventory",
      <|"Audit" -> profileAudit|>]
  ];

  If[profileProbe === "inspect",
    D2SaveStatus["ProfileTreeInspected", <|
      "OperatorAudit" -> FileNameJoin[{profileDir, "operator_sector_audit.json"}]
    |>];
    Print["[profile] Structural inspection complete; no Match[] called."];
    Exit[0]
  ];

  profileInput = Switch[profileProbe,
    "sm-only", profileSMOnly,
    "s2r-only", profileS2ROnly,
    "s1-only", profileS1Only,
    "heavy-only", profileHeavyOnly,
    "free", profileFree,
    "scalar-standard", profileStandard,
    "scalar-general", profileGeneral,
    "scalar-with-y2", RelabelIndices[Expand[profileGeneral + profileY2]],
    "scalar-with-y2Conj", RelabelIndices[Expand[profileGeneral + profileY2Conj]],
    "bsm-free", profileBSMFree,
    "bsm-standard", profileBSMStandard,
    "bsm-general", profileBSMGeneral,
    "bsm-with-y2", profileBSMWithY2,
    "bsm-with-y2Conj", profileBSMWithY2Conj,
    "bsm-with-both", profileBSMWithBoth,
    (* Legacy baselines remain available only as explicit diagnostic modes. *)
    "renormalisable", profileRen,
    "with_y2", Expand[profileRen + profileY2],
    "with_y2Conj", Expand[profileRen + profileY2Conj],
    "full-tree", transitionTree,
    _, $Failed
  ];
  If[profileInput === $Failed,
    D2Stop[32, "Unrecognized diagnostic input"]];
  Print["[profile] Starting ", profileProbe,
    " tree Match[] only; Wolfram limit ", profileLimit, " s"];
  Print["[profile] Input bytes: ", ByteCount[profileInput],
    "; leaves: ", LeafCount[profileInput],
    "; memory-in-use: ", MemoryInUse[]];
  profileMatchStarted = AbsoluteTime[];
  (* Stage 1 successfully runs Matchete matching inside UsingFrontEnd.
     Do the same here: Matchete's headless Match[] may invoke frontend
     services even after a separate UsingFrontEnd[Needs[...]] returned. *)
  profileMatch = TimeConstrained[
    CheckAbort[
      UsingFrontEnd[
        Match[profileInput, EFTOrder -> 5, LoopOrder -> 0]
      ],
      $Aborted
    ],
    profileLimit, $TimedOut
  ];
  profileMatchSeconds = AbsoluteTime[] - profileMatchStarted;
  profileResultStatus = Which[
    profileMatch === $TimedOut, "TimedOut",
    profileMatch === $Failed, "Failed",
    profileMatch === $Aborted, "Aborted",
    Head[profileMatch] === Match, "UnevaluatedMatch",
    True, "MatchReturned"
  ];
  profileResult = <|
    "Status" -> profileResultStatus,
    "Probe" -> profileProbe,
    "ElapsedSeconds" -> profileMatchSeconds,
    "InputByteCount" -> ByteCount[profileInput],
    "InputLeafCount" -> LeafCount[profileInput],
    "ResultHead" -> D2Text[Head[profileMatch]],
    "ResultByteCount" -> If[profileResultStatus === "MatchReturned",
      ByteCount[profileMatch], 0],
    "MemoryInUse" -> MemoryInUse[],
    "MaxMemoryUsed" -> MaxMemoryUsed[],
    "DiagnosticOnly" -> True,
    "FullGeneralMatchingCertified" -> False
  |>;
  profileResultPath = FileNameJoin[{profileDir,
    "probe_result_" <> profileProbe <> ".json"}];
  profileExportStatus = Quiet @ Check[
    Export[profileResultPath, profileResult, "RawJSON"], $Failed];
  If[profileExportStatus === $Failed || !FileExistsQ[profileResultPath],
    D2Stop[35, "Profile probe result JSON export failed",
      <|"Probe" -> profileProbe, "Path" -> profileResultPath|>]
  ];
  D2SaveStatus["ProfileTreeProbeFinished", <|
    "Probe" -> profileProbe,
    "ProbeStatus" -> profileResultStatus,
    "ElapsedSeconds" -> profileMatchSeconds,
    "ResultPath" -> profileResultPath
  |>];
  Print["[profile] ", profileProbe, ": ", profileResultStatus,
    "; elapsed s: ", NumberForm[profileMatchSeconds, {8, 2}]];
  Exit[If[profileResultStatus === "MatchReturned" ||
    profileResultStatus === "TimedOut", 0, 33]]
];

(* Consistency with successful general D(-2) first threshold:
   keep Matchete's full sequence (Match, GreensSimplify, EOMSimplify,
   EvaluateLoopFunctions, ReplaceEffectiveCouplings) inside a live
   service frontend for ALL three fixed-order stages. No physics changed. *)
(* [A] Matching the tree EFT at tree order. *)
Print["[D(-2) scalar threshold] [A] Tree Match of EFT1 L0"];
Print["[D(-2) scalar threshold] [A] tree input: ",
  ByteCount[transitionTree], " bytes; ", LeafCount[transitionTree], " leaves"];
d2StageAStart = AbsoluteTime[];
treeResult = UsingFrontEnd[RunT3Matching[transitionTree, 5, 0]];
Print["[D(-2) scalar threshold] [A] elapsed seconds: ",
  AbsoluteTime[] - d2StageAStart];
If[!D2Success[treeResult],
  D2Stop[21, "Scalar tree matching failed",
    <|"StageResult" -> If[AssociationQ[treeResult],
      Lookup[treeResult,"Status","Unknown"],D2Text[treeResult]]|>]];
treeMatched = CompactExplicitEFT[treeResult];
If[treeMatched === $Failed,
  D2Stop[22, "Scalar tree EFT could not be expanded"]];
Export[FileNameJoin[{stageDir, "scalar_hard_tree.wxf"}], treeMatched, "WXF"];

(* [B] Local hard one-loop matching acts ONLY on tree-level EFT1. *)
Print["[D(-2) scalar threshold] [B] One-loop Match of EFT1 L0"];
d2StageBStart = AbsoluteTime[];
loopResult = UsingFrontEnd[RunT3Matching[transitionTree, 5, 1]];
Print["[D(-2) scalar threshold] [B] elapsed seconds: ",
  AbsoluteTime[] - d2StageBStart];
If[!D2Success[loopResult], D2Stop[23, "Scalar one-loop matching failed",
  <|"StageResult" -> If[AssociationQ[loopResult],
    Lookup[loopResult,"Status","Unknown"],D2Text[loopResult]]|>]];
loopMatched = CompactExplicitEFT[loopResult];
If[loopMatched === $Failed,
  D2Stop[24, "Scalar one-loop EFT could not be expanded"]];
localLoop = Expand[loopMatched - treeMatched];
Export[FileNameJoin[{stageDir, "scalar_hard_local_one_loop.wxf"}],
  localLoop, "WXF"];

(* [C] Inherited one-loop EFT1 correction, propagated by TREE matching. *)
Print["[D(-2) scalar threshold] [C] Tree propagation of inherited O(hbar) EFT1"];
d2StageCStart = AbsoluteTime[];
probeResult = UsingFrontEnd[RunT3Matching[transitionFullCanonical, 5, 0]];
Print["[D(-2) scalar threshold] [C] elapsed seconds: ",
  AbsoluteTime[] - d2StageCStart];
If[!D2Success[probeResult], D2Stop[25, "Inherited one-loop tree propagation failed",
  <|"StageResult" -> If[AssociationQ[probeResult],
    Lookup[probeResult,"Status","Unknown"],D2Text[probeResult]]|>]];
probeMatched = CompactExplicitEFT[probeResult];
If[probeMatched === $Failed || FreeQ[probeMatched, hbar],
  D2Stop[26, "Inherited EFT lacks explicit hbar after tree propagation"]];
inheritedLoop = Expand[hbar Coefficient[
  Expand[probeMatched - treeMatched], hbar, 1]];
Export[FileNameJoin[{stageDir, "scalar_inherited_one_loop.wxf"}],
  inheritedLoop, "WXF"];

(* Fixed-order SMEFT: tree plus two disjoint one-loop sources.  The direct
   running piece remains a separate pure-SM coefficient with NO hbar in the
   generated .wl, and is combined only after scalar hard matching. *)
hardEFT = Expand[treeMatched + localLoop + inheritedLoop];
hardText = D2Text[hardEFT];
If[AnyTrue[{"NewScalar1", "T3S2R", "T3S2I", "NewFermion", "NewScalar2"},
     StringContainsQ[hardText, "Field[" <> # <> ","] &],
  D2Stop[28, "Matched EFT still contains integrated-out BSM fields"]
];
Export[FileNameJoin[{stageDir, "scalar_hard_smeft.wxf"}], hardEFT, "WXF"];
Export[FileNameJoin[{stageDir, "scalar_hard_smeft_inputform.txt"}],
  D2Text[hardEFT], "Text"];

Print["[D(-2) scalar threshold] Extracting hard Weinberg operator..."];
hardExtraction = Quiet @ Check[
  ExtractWeinbergCoefficient[hardEFT], $Failed];
hardStatus = If[AssociationQ[hardExtraction],
  Lookup[hardExtraction,"Status","Unknown"], "Failed"];
hardC5 = If[hardStatus === "Success",
  Lookup[hardExtraction,"Coefficient", $Failed], $Failed];
If[hardC5 === $Failed || MissingQ[hardC5],
  D2Stop[27, "Hard scalar C5 extraction not certified",
    <|"ExtractionStatus" -> hardStatus,
      "HardSMEFTExported" -> True|>]
];
Export[FileNameJoin[{stageDir, "scalar_hard_c5_inputform.txt"}],
  D2Text[hardC5], "Text"];

(* Do not blindly add hardC5 to a p,q-indexed tensor if Matchete's
   extracted hard coefficient contains its own dummy flavours.  Export a
   fixed-order assembly formula and each independently checkable source.
   A tensor index-map check is required before physical C5 certification. *)
Export[FileNameJoin[{stageDir, "smeft_c5_assembly.wl"}],
  "(* Matchete native loop bookkeeping: hbar=1/(16 Pi^2). *)\n" <>
  "(* HARD_C5: see scalar_hard_c5_inputform.txt (Matchete flavour dummies). *)\n" <>
  "(* RUN_C5: see direct_running_c5_inputform.txt (external p,q). *)\n" <>
  "(* C5_SMEFT[p,q] = HARD_C5[p,q] + hbar RUN_C5[p,q]. *)\n",
  "Text"];
D2SaveStatus["HardScalarThresholdMatched_FlavorAssemblyPending", <|
  "TreeMatching" -> True,
  "OneLoopMatching" -> True,
  "InheritedLoopTreePropagation" -> True,
  "HardSMEFTExported" -> True,
  "HardC5ExtractionStatus" -> hardStatus,
  "HardC5InputForm" -> FileNameJoin[{stageDir, "scalar_hard_c5_inputform.txt"}],
  "RunningC5InputForm" -> FileNameJoin[{stageDir, "direct_running_c5_inputform.txt"}],
  "SMEFTC5AssemblyPath" -> FileNameJoin[{stageDir, "smeft_c5_assembly.wl"}],
  "PhysicalMassRelationsImposed" -> False,
  "FullFlavorHardC5Verified" -> False,
  "AllScalarFieldsIntegratedOut" -> True
|>];
Print["SUCCESS: scalar threshold and hard C5 extracted; full-flavour hard/index assembly pending."];
Exit[0];
