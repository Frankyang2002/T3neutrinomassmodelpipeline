(* General D(-2) fermion-first matching in a diagonal neutral-scalar basis.

   The general UV model includes kappa (S2^dagger)^2 + h.c.  This is an
   off-diagonal quadratic form for the complex singlet S2.  Matchete's
   DetermineCurrentFieldAssociation cannot match with that original mass
   basis, even while S2 is Light.  We therefore use the exact field
   decomposition S2 = eta (R + I J)/Sqrt[2], where J denotes the imaginary
   real scalar (written T3S2I below) and eta has unit modulus.

   Choosing kappa/eta^2 = Abs[kappa], the Lagrangian convention +kappa
   (S2^dagger)^2 + h.c. gives m_R^2 = MS2^2 - 2 Abs[kappa]
   and m_J^2 = MS2^2 + 2 Abs[kappa].  The new masses and eta are declared
   as independent symbolic parameters to Matchete.  Their physical relations
   MUST be imposed downstream.  No restricted interaction is discarded.

   Only F is integrated out.  The ordinary complex-S2 EFT1 tensor adapters
   are NOT compatible with R,J fields, so this script exports explicit real-
   basis EFT1 data rather than the old misleading scalar-quartic/Wilson seed.

   CLI: wolframscript -file Lagrangian/RunGeneralDMinus2Stage1.wl OUTPUT
*)
ClearAll["Global`*"];
scriptDir = DirectoryName @ ExpandFileName[$InputFileName];
args = Rest[$ScriptCommandLine];
If[Length[args] =!= 1, Print["Usage: OUTPUT_DIR"]; Exit[2]];
outputDir = ExpandFileName[args[[1]]];
dataDir = FileNameJoin[{outputDir, "data"}];
If[Environment["T3_EXPERIMENTAL_GENERAL_QUADRATICS"] =!= "1" ||
   Environment["T3_EXPERIMENTAL_GENERAL_FERMIONS"] =!= "1",
  Print["ERROR: both general-interaction flags must be 1."]; Exit[3]
];
Quiet @ CreateDirectory[dataDir, CreateIntermediateDirectories -> True];

loaded = CheckAbort[Check[UsingFrontEnd[Needs["Matchete`"] ; True],
  $Failed], $Aborted];
If[loaded =!= True,
  Print["ERROR: Matchete load failed."]; Exit[4]
];
Get[FileNameJoin[{scriptDir, "PhysicsNotation.wl"}]];
Get[FileNameJoin[{scriptDir, "model", "T3ModelCatalog.wl"}]];
Get[FileNameJoin[{scriptDir, "model", "LagrangianBuilder.wl"}]];
Get[FileNameJoin[{scriptDir, "RunMatching.wl"}]];

model = T3ModelFromClass["D", -2];
If[model === $Failed || model["Scalar2", "SU2"] =!= 1 ||
   model["Scalar2", "Y"] =!= 0,
  Print["ERROR: expected neutral SU(2) singlet S2 for D(-2)."]; Exit[5]
];
If[SetT3HeavyFields[{"F"}] === $Failed, Exit[6]];
build = CheckAbort[UsingFrontEnd[BuildT3Lagrangian[model]], $Aborted];
If[!AssociationQ[build] || Lookup[build, "Status", ""] =!= "Success",
  Print["ERROR: full-general UV build failed: ",
    If[AssociationQ[build], Lookup[build, "Status", "Unknown"], build]];
  Exit[7]
];

accepted = ToString /@ Lookup[build, "AllowedInteractions", {}];
extraQ = Select[accepted, StringStartsQ[#, "T3ExtraQ_"] &];
extraM = Select[accepted, StringStartsQ[#, "T3ExtraM_"] &];
extraY = Select[accepted, StringStartsQ[#, "Yukawa2ConjugateScalar"] &];
If[Length[extraQ] =!= 5 || Length[extraM] =!= 1 ||
   Length[extraY] =!= 1 ||
   Lookup[build, "AdditionalQuarticsIncluded", -1] =!= 5 ||
   Lookup[build, "AdditionalQuadraticsIncluded", -1] =!= 1 ||
   Lookup[build, "AdditionalYukawasIncluded", -1] =!= 1 ||
   !TrueQ[Lookup[build, "AdditionalYukawasRequested", False]] ||
   !TrueQ[Lookup[build, "AdditionalQuadraticsRequested", False]],
  Print["ERROR: incomplete general interaction inventory."]; Exit[8]
];

Export[FileNameJoin[{outputDir, "general_uv_inventory.json"}], <|
  "Model" -> "D_-2", "BuildStatus" -> build["Status"],
  "AllowedInteractions" -> accepted,
  "AdditionalQuartics" -> extraQ,
  "AdditionalQuadratics" -> extraM,
  "AdditionalYukawas" -> extraY,
  "AdditionalQuarticsCount" -> Length[extraQ],
  "AdditionalQuadraticsCount" -> Length[extraM],
  "AdditionalYukawasCount" -> Length[extraY],
  "HeavyFieldsAtFirstThreshold" -> {"F"},
  "SymmetryCompletenessCertified" ->
    TrueQ[Lookup[build, "FullSymmetryCompletenessCertified", False]],
  "GeneralFermionFlavorMassCertified" ->
    TrueQ[Lookup[build, "FermionFlavorMassGeneralityCertified", False]]
|>, "RawJSON"];

(* Remove the original singlet before defining its two REAL mass eigenfields.
   F remains Heavy; S1, R and J remain Light at the first threshold. *)
(* Register each field independently.  Matchete can emit nonfatal
   messages when removing/redefining fields; a blanket Check[...;$Failed]
   wrongly turns any message into a hard failure and hides which definition
   was actually rejected.  Validate the resulting public registries instead. *)
removeResult = CheckAbort[RemoveField[NewScalar2], $Aborted];
Print["[mass basis] RemoveField[NewScalar2]: ", InputForm[removeResult]];
If[removeResult === $Aborted,
  Print["ERROR: old scalar removal aborted."]; Exit[9]
];

realRResult = CheckAbort[
  DefineField[T3S2R, Scalar, SelfConjugate -> True,
    Mass -> {Light, T3S2RMass}],
  $Aborted
];
Print["[mass basis] DefineField[T3S2R]: ", InputForm[realRResult]];

realIResult = CheckAbort[
  DefineField[T3S2I, Scalar, SelfConjugate -> True,
    Mass -> {Light, T3S2IMass}],
  $Aborted
];
Print["[mass basis] DefineField[T3S2I]: ", InputForm[realIResult]];

phaseResult = CheckAbort[
  DefineCoupling[T3S2Phase, SelfConjugate -> False],
  $Aborted
];
Print["[mass basis] DefineCoupling[T3S2Phase]: ", InputForm[phaseResult]];

(* These are Matchete's public registries.  Inspect the registrations,
   not the printed return values (which may depend on package version). *)
registeredFields = CheckAbort[GetFields[], $Aborted];
registeredCouplings = CheckAbort[GetCouplings[], $Aborted];
fieldRegistryText = ToString[registeredFields, InputForm];
couplingRegistryText = ToString[registeredCouplings, InputForm];
realRRegistered = StringContainsQ[fieldRegistryText, "T3S2R"];
realIRegistered = StringContainsQ[fieldRegistryText, "T3S2I"];
oldS2Registered = StringContainsQ[fieldRegistryText, "NewScalar2"];
phaseRegistered = StringContainsQ[couplingRegistryText, "T3S2Phase"];
registrationOk = !MemberQ[
    {removeResult, realRResult, realIResult, phaseResult,
     registeredFields, registeredCouplings},
    $Aborted | $Failed
  ] && realRRegistered && realIRegistered &&
  !oldS2Registered && phaseRegistered;

registrationAudit = <|
  "Status" -> If[registrationOk, "Success", "Failed"],
  "RemoveOldScalarReturn" -> ToString[removeResult, InputForm],
  "RealRDefinitionReturn" -> ToString[realRResult, InputForm],
  "RealIDefinitionReturn" -> ToString[realIResult, InputForm],
  "PhaseDefinitionReturn" -> ToString[phaseResult, InputForm],
  "RealRRegistered" -> realRRegistered,
  "RealIRegistered" -> realIRegistered,
  "OriginalComplexS2StillRegistered" -> oldS2Registered,
  "PhaseRegistered" -> phaseRegistered,
  "FieldRegistryPreview" -> StringTake[fieldRegistryText, UpTo[4000]],
  "CouplingRegistryHasPhase" -> phaseRegistered
|>;
Export[FileNameJoin[{outputDir, "scalar_registration_debug.json"}],
  registrationAudit, "RawJSON"];
Print["[mass basis] Registry audit: ",
  {realRRegistered, realIRegistered, oldS2Registered, phaseRegistered}];
If[!registrationOk,
  Print["ERROR: real singlet registry invalid. See scalar_registration_debug.json."];
  Exit[9]
];
eta = T3S2Phase[];

(* Treat the neutral holomorphic bilinear as part of the physical mass
   eigenvalues, NOT as an additional interaction insertion.  Both it and its
   HC must be removed, exactly twice.  The other terms are retained. *)
kappaName = "lambdaT3ExtraM000002I1";
uvInteraction = Expand[build["LInt"]];
terms = If[Head[uvInteraction] === Plus,
  List @@ uvInteraction, {uvInteraction}];
removed = Select[terms,
  StringContainsQ[ToString[#, InputForm], kappaName] &];
kept = Select[terms,
  !StringContainsQ[ToString[#, InputForm], kappaName] &];
If[Length[removed] =!= 2,
  Print["ERROR: expected exactly two conjugate S2 mass-bilinear terms; got ",
    Length[removed]]; Exit[10]
];

(* S2 = eta (R + i J)/sqrt(2).  These rules cover both the field and its
   conjugate, including any derivative labels, without changing SU(2) CGs. *)
realRules = {
  HoldPattern[Bar[Field[NewScalar2, Scalar, ix_, dx_]]] :>
    Bar[eta] (Field[T3S2R, Scalar, ix, dx] -
      I Field[T3S2I, Scalar, ix, dx])/Sqrt[2],
  HoldPattern[Field[NewScalar2, Scalar, ix_, dx_]] :>
    eta (Field[T3S2R, Scalar, ix, dx] +
      I Field[T3S2I, Scalar, ix, dx])/Sqrt[2]
};
rotatedInteraction = RelabelIndices[Expand[Total[kept] /. realRules]];
physicalFree = RelabelIndices[
  FreeLag[NewFermion, NewScalar1, T3S2R, T3S2I]
];
physicalLag = RelabelIndices[
  Expand[build["LSM"] + physicalFree + rotatedInteraction]
];
physicalText = ToString[physicalLag, InputForm];
If[StringContainsQ[physicalText, "Field[NewScalar2"] ||
   StringContainsQ[physicalText, kappaName] ||
   !StringContainsQ[physicalText, "Field[T3S2R"] ||
   !StringContainsQ[physicalText, "Field[T3S2I"],
  Print["ERROR: the real-scalar rotation is incomplete."]; Exit[11]
];
(* Do not allow a new quartic/Yukawa to disappear in the transformation. *)
quarticSymbols = Map[Function[q, Module[{parts},
  parts = StringSplit[StringDrop[q, StringLength["T3ExtraQ_"]], "_"];
  If[Length[parts] =!= 7, Return[$Failed]];
  "lambdaT3ExtraQ" <> StringJoin[Most[parts]] <> "I" <> Last[parts]
]], extraQ];
If[MemberQ[quarticSymbols, $Failed] ||
   !AllTrue[quarticSymbols,
     StringContainsQ[physicalText, #] &] ||
   !StringContainsQ[physicalText, "y2Conj"],
  Print["ERROR: general coupling lost in physical mass-basis rotation."];
  Exit[12]
];

Export[FileNameJoin[{outputDir, "scalar_mass_basis.json"}], <|
  "Model" -> "D_-2",
  "Status" -> "RealSingletMassBasisConstructedNotNumericallyCertified",
  "OriginalComplexField" -> "NewScalar2",
  "RealFields" -> {"T3S2R", "T3S2I"},
  "PhysicalFieldDefinition" -> "S2=eta*(R+i*J)/sqrt(2)",
  "PhaseConstraint" -> "Abs[eta]=1; kappa/eta^2=Abs[kappa]",
  "HolomorphicCoupling" -> kappaName,
  "LagrangianBilinearConvention" -> "+kappa*(S2dag)^2 + h.c.",
  "RealMassSquared" -> "MS2^2-2*Abs[kappa]",
  "ImagMassSquared" -> "MS2^2+2*Abs[kappa]",
  "MassPositivityCondition" -> "MS2^2>2*Abs[kappa]",
  "MatcheteRealMassParameters" -> {"T3S2RMass", "T3S2IMass"},
  "RemovedConjugateBilinearTermCount" -> Length[removed],
  "GeneralQuarticCouplingNames" -> quarticSymbols,
  "PhysicalMassRelationsImposedInMatchete" -> False,
  "ScalarThresholdPerformed" -> False
|>, "RawJSON"];
Export[FileNameJoin[{dataDir, "physical_uv_realbasis_inputform.txt"}],
  physicalText, "Text"];

Print["Checking rebuilt real-scalar free theory..."];
freeCheck = CheckAbort[
  UsingFrontEnd[CheckLagrangian[RelabelIndices[build["LSM"] + physicalFree]]],
  $Aborted
];
Print["  free check: ", InputForm[freeCheck]];
If[freeCheck =!= True, Exit[13]];
Print["Checking FULL general D(-2) physical Lagrangian..."];
fullCheck = CheckAbort[
  UsingFrontEnd[CheckLagrangian[physicalLag]], $Aborted
];
Print["  physical check: ", InputForm[fullCheck]];
If[fullCheck =!= True,
  Print["ERROR: mass-basis full Lagrangian is not valid."]; Exit[14]
];

(* One stage only; avoid RunSequentialT3Matching's old NewScalar2-only
   intermediate tensor exporter, which would misclassify the real fields. *)
Print["=== D(-2), real scalar basis: tree matching ONLY F ==="];
treeResult = CheckAbort[
  UsingFrontEnd[RunT3Matching[physicalLag, 5, 0]], $Aborted
];
If[!AssociationQ[treeResult] ||
   Lookup[treeResult, "Status", ""] =!= "Success",
  Print["ERROR: fermion tree matching failed: ",
    If[AssociationQ[treeResult], Lookup[treeResult, "Status", ""], treeResult]];
  Exit[15]
];
tree = CompactExplicitEFT[treeResult];
transitionTree = CanonicalTransitionEFT[treeResult];
If[MemberQ[{$Failed, $Aborted}, tree] ||
   MemberQ[{$Failed, $Aborted}, transitionTree],
  Print["ERROR: fermion tree EFT canonicalisation failed."]; Exit[16]
];
Print["=== D(-2), real scalar basis: one-loop matching ONLY F ==="];
loopResult = CheckAbort[
  UsingFrontEnd[RunT3Matching[physicalLag, 5, 1]], $Aborted
];
If[!AssociationQ[loopResult] ||
   Lookup[loopResult, "Status", ""] =!= "Success",
  Print["ERROR: fermion one-loop matching failed: ",
    If[AssociationQ[loopResult], Lookup[loopResult, "Status", ""], loopResult]];
  Exit[17]
];
full = CompactExplicitEFT[loopResult];
transitionFull = CanonicalTransitionEFT[loopResult];
If[MemberQ[{$Failed, $Aborted}, full] ||
   MemberQ[{$Failed, $Aborted}, transitionFull],
  Print["ERROR: fermion one-loop EFT canonicalisation failed."]; Exit[18]
];
loopCorrection = Expand[full - tree];

(* Export native results and the exact physical-basis operator inventory. *)
Export[FileNameJoin[{dataDir, "eft1_after_F_realbasis_tree.wxf"}], tree, "WXF"];
Export[FileNameJoin[{dataDir, "eft1_after_F_realbasis_loop.wxf"}],
  loopCorrection, "WXF"];
Export[FileNameJoin[{dataDir, "eft1_after_F_realbasis_transition_tree.wxf"}],
  transitionTree, "WXF"];
Export[FileNameJoin[{dataDir, "eft1_after_F_realbasis_transition_full.wxf"}],
  transitionFull, "WXF"];
Export[FileNameJoin[{dataDir, "eft1_after_F_realbasis_tree_inputform.txt"}],
  ToString[tree, InputForm], "Text"];
Export[FileNameJoin[{dataDir, "eft1_after_F_realbasis_loop_inputform.txt"}],
  ToString[loopCorrection, InputForm], "Text"];

expanded = Expand[tree];
allTerms = If[Head[expanded] === Plus, List @@ expanded, {expanded}];

(* A field written as Field[...]^n or Bar[Field[...]]^n counts as n
   scalar legs.  StringCount[text, "Field[..."] counts only the base once,
   so it silently excludes many holomorphic D(-2) quartics.  Count powers
   without evaluating or rewriting Matchete's private-context Field objects. *)
ClearAll[T3D2MatchingSquareEnd, T3D2ScalarFactorDegree];
T3D2MatchingSquareEnd[text_String, first_Integer] := Module[
  {depth = 0, p, char},
  For[p = first, p <= StringLength[text], p++,
    char = StringTake[text, {p}];
    If[char === "[", depth++];
    If[char === "]",
      depth--;
      If[depth === 0, Return[p]]
    ];
  ];
  $Failed
];
T3D2ScalarFactorDegree[text_String] := Module[
  {starts, total = 0, at, fragment, close, bar, tail,
   powers, degree, expected = {"H", "NewScalar1", "T3S2R", "T3S2I"}},
  starts = StringPosition[text, "Field["];
  Do[
    at = pos[[1]];
    fragment = StringDrop[text, at - 1];
    If[!AnyTrue[expected,
       StringStartsQ[fragment, "Field[" <> # <> ", Scalar"] &],
       Continue[]
    ];
    close = T3D2MatchingSquareEnd[text, at + 5];
    If[close === $Failed, Return[$Failed]];
    bar = at >= 5 &&
      StringTake[text, {at - 4, at - 1}] === "Bar[" &&
      close < StringLength[text] &&
      StringTake[text, {close + 1}] === "]";
    tail = StringDrop[text, close + If[bar, 1, 0]];
    powers = StringCases[tail, RegularExpression["^\\^[0-9]+"]];
    degree = If[powers === {}, 1, ToExpression[StringDrop[First[powers], 1]]];
    If[!IntegerQ[degree] || degree < 1 || degree > 8,
      Return[$Failed]
    ];
    total += degree,
    {pos, starts}
  ];
  total
];

wilsonTerms = Select[allTerms, Function[t,
  Module[{s = ToString[t, InputForm]},
    (* Existing LLSS export does not include operators with Higgs legs. *)
    StringContainsQ[s, "GammaCC"] &&
      T3D2ScalarFactorDegree[s] === 2 &&
      !StringContainsQ[s, "Field[H, Scalar"]
  ]
]];
treeText = ToString[tree, InputForm];
fermionAbsent = !StringContainsQ[treeText, "Field[NewFermion"];
s1Present = StringContainsQ[treeText, "Field[NewScalar1"];
rPresent = StringContainsQ[treeText, "Field[T3S2R"];
iPresent = StringContainsQ[treeText, "Field[T3S2I"];
oldS2Absent = !StringContainsQ[treeText, "Field[NewScalar2"];
qInEft = AllTrue[quarticSymbols, StringContainsQ[treeText, #] &];
yConjInEft = StringContainsQ[treeText, "y2Conj"];
If[!fermionAbsent || !s1Present || !rPresent || !iPresent ||
   !oldS2Absent || !qInEft || !yConjInEft || Length[wilsonTerms] == 0,
  Print["ERROR: real-basis EFT1 operator-preservation check failed: ",
    {fermionAbsent,s1Present,rPresent,iPresent,oldS2Absent,qInEft,
     yConjInEft,Length[wilsonTerms]}]; Exit[19]
];
Export[FileNameJoin[{dataDir, "eft1_after_F_realbasis_wilson_inventory.json"}], <|
  "Format" -> "D_m2_realbasis_wilson_inventory_not_tensor_adapter_input",
  "Status" -> "Success",
  "TreeWilsonTermCount" -> Length[wilsonTerms],
  "TreeWilsonTermsInputForm" -> (ToString[#, InputForm] & /@ wilsonTerms),
  "PhysicalScalarNames" -> {"NewScalar1", "T3S2R", "T3S2I"},
  "StandardComplexScalarTensorAdapterCompatible" -> False
|>, "RawJSON"];
Export[FileNameJoin[{outputDir, "first_threshold_status.json"}], <|
  "Status" -> "Success", "Model" -> "D_-2",
  "ThresholdPlan" -> {{"F"}}, "IntegratedFields" -> {"F"},
  "ActiveFields" -> {"S1", "S2R", "S2I"},
  "FermionAbsentFromEFT1Tree" -> fermionAbsent,
  "Scalar1PresentInEFT1Tree" -> s1Present,
  "RealScalarRPresentInEFT1Tree" -> rPresent,
  "RealScalarIPresentInEFT1Tree" -> iPresent,
  "OriginalComplexS2AbsentFromEFT1Tree" -> oldS2Absent,
  "GeneralExtraQuarticsRetained" -> qInEft,
  "ExtraYukawaRetainedInEFT1" -> yConjInEft,
  "TreeWilsonTermCount" -> Length[wilsonTerms],
  "ScalarMassBasisDiagonalized" -> True,
  "PhysicalMassRelationsImposed" -> False,
  "ScalarThresholdPerformed" -> False,
  "IntermediateWilsonRunningPerformed" -> False,
  "StandardComplexScalarTensorAdapterCompatible" -> False
|>, "RawJSON"];

(* -------------------------------------------------------------------- *)
(* Export exact D(-2) real-basis EFT1 tensor inputs.                     *)
(*                                                                         *)
(* This deliberately does not invoke the original complex-S2 exporter.     *)
(* CG tensors are read from the SAME live kernel as fermion matching, so    *)
(* that tree-level effective SU(2) contractions are exported exactly.      *)
(* -------------------------------------------------------------------- *)
(* Matchete can place CG heads in private contexts after matching.  In
   particular a literal Cases[..., CG[...], Infinity] can miss T3MixCG even
   though the InputForm terms contain it.  Harvest the exact serialized CG
   names, as the Python tensor adapter does, and resolve those names against
   the SAME live registry.  Do not infer or substitute CG components. *)
ClearAll[T3D2CGNames, T3D2CGRegistry];
T3D2CGNames[terms_List] := DeleteDuplicates @ Flatten[
  Map[
    Function[term,
      StringCases[
        ToString[term, InputForm],
        RegularExpression["CG\\[\\s*(?:Bar\\[\\s*)?([A-Za-z$][A-Za-z0-9$]*)"] -> "$1"
      ]
    ],
    terms
  ]
];

T3D2CGRegistry[terms_List] := Module[
  {names, items, reps, tensor, name, nameText, repText, tensorText},
  names = T3D2CGNames[terms];
  items = Table[
    nameText = names[[k]];
    (* Restricted ASCII symbol grammar above prevents arbitrary WL input.
       Resolving in this kernel preserves the existing CG registry. *)
    name = Quiet @ Check[ToExpression[nameText, InputForm], $Failed];
    If[name === $Failed || Head[name] =!= Symbol,
      Print["ERROR: CG name cannot be resolved: ", nameText];
      Return[$Failed]
    ];
    reps = Quiet @ Check[
      Matchete`PackageScope`$CGproperties[name, Indices], $Failed];
    tensor = Quiet @ Check[
      Matchete`CGManipulations`PackagePrivate`$CGtensors[name], $Failed];
    repText = ToString[reps, InputForm];
    tensorText = ToString[tensor, InputForm];
    If[MemberQ[{reps, tensor}, $Failed] ||
       MissingQ[reps] || MissingQ[tensor] ||
       !StringStartsQ[repText, "{"] ||
       !StringStartsQ[tensorText, "SparseArray["],
      Print["ERROR: no usable registered CG for ", nameText,
        "; representations: ", StringTake[repText, UpTo[120]],
        "; tensor: ", StringTake[tensorText, UpTo[120]]];
      Return[$Failed]
    ];
    <|"Name" -> nameText,
      "RepsInputForm" -> repText,
      "TensorInputForm" -> tensorText|>,
    {k, Length[names]}
  ];
  If[!ListQ[items] || Length[items] =!= Length[names],
     Return[$Failed]];
  items
];

realWilsonRecords = MapIndexed[Function[{t, pos},
  Module[{s = ToString[t, InputForm]},
    <|"Index" -> First[pos],
      "Chirality" -> Which[
        StringContainsQ[s, "Proj[-1]"], "PL",
        StringContainsQ[s, "Proj[1]"], "PR",
        True, "Unknown"],
      "TermInputForm" -> s|>
  ]], wilsonTerms];

realQuarticTerms = Select[allTerms, Function[t,
  Module[{s = ToString[t, InputForm]},
    T3D2ScalarFactorDegree[s] === 4 &&
    !StringContainsQ[s, ", Fermion,"] &&
    !StringContainsQ[s, "FieldStrength["] &&
    !StringContainsQ[s, "GammaM["] &&
    !StringContainsQ[s, "CovariantD["]
  ]
]];
realQuarticRecords = MapIndexed[Function[{t, pos},
  <|"Index" -> First[pos], "TermInputForm" -> ToString[t, InputForm]|>
], realQuarticTerms];

realWilsonCGNames = T3D2CGNames[wilsonTerms];
realQuarticCGNames = T3D2CGNames[realQuarticTerms];
Print["[real EFT1 CG] Wilson references: ", realWilsonCGNames];
Print["[real EFT1 CG] Quartic references: ", realQuarticCGNames];
realWilsonCGRegistry = T3D2CGRegistry[wilsonTerms];
realQuarticCGRegistry = T3D2CGRegistry[realQuarticTerms];
realCGExportedWilson = If[ListQ[realWilsonCGRegistry],
  Lookup[realWilsonCGRegistry, "Name", {}], {}];
realCGExportedQuartic = If[ListQ[realQuarticCGRegistry],
  Lookup[realQuarticCGRegistry, "Name", {}], {}];
realCGMissingWilson = Complement[realWilsonCGNames, realCGExportedWilson];
realCGMissingQuartic = Complement[realQuarticCGNames, realCGExportedQuartic];
Export[FileNameJoin[{dataDir, "eft1_after_F_realbasis_cg_registry_audit.json"}],
  <|
    "Status" -> If[realWilsonCGRegistry =!= $Failed &&
      realQuarticCGRegistry =!= $Failed &&
      realCGMissingWilson === {} && realCGMissingQuartic === {},
      "Success", "Failed"],
    "ReferencedWilsonCGNames" -> realWilsonCGNames,
    "ExportedWilsonCGNames" -> realCGExportedWilson,
    "MissingWilsonCGNames" -> realCGMissingWilson,
    "ReferencedQuarticCGNames" -> realQuarticCGNames,
    "ExportedQuarticCGNames" -> realCGExportedQuartic,
    "MissingQuarticCGNames" -> realCGMissingQuartic
  |>, "RawJSON"];
If[realWilsonCGRegistry === $Failed || realQuarticCGRegistry === $Failed ||
   realCGMissingWilson =!= {} || realCGMissingQuartic =!= {},
  Print["ERROR: real-basis CG export incomplete; Wilson missing: ",
    realCGMissingWilson, "; quartic missing: ", realCGMissingQuartic];
  Print["See data/eft1_after_F_realbasis_cg_registry_audit.json"];
  Exit[20]
];
realQuarticTexts = Lookup[realQuarticRecords, "TermInputForm", {}];
(* Keep a machine-readable inventory even when selection fails.  Every
   original extra coupling must be visible in the SELECTED quartic terms. *)
missingQuartics = Select[quarticSymbols,
  Function[q, !AnyTrue[realQuarticTexts, StringContainsQ[#, q] &]]
];
seedAudit = <|
  "Status" -> If[Length[realWilsonRecords] === Length[wilsonTerms] &&
    Length[realQuarticRecords] > 0 && missingQuartics === {},
    "Success", "Failed"],
  "TreeWilsonTerms" -> Length[wilsonTerms],
  "WilsonCGNamesReferenced" -> realWilsonCGNames,
  "QuarticCGNamesReferenced" -> realQuarticCGNames,
  "CGRegistryComplete" ->
    (realCGMissingWilson === {} && realCGMissingQuartic === {}),
  "SelectedQuarticTermCount" -> Length[realQuarticRecords],
  "ExtraQuarticsMissingFromSelectedTerms" -> missingQuartics,
  "GeneralQuarticsPresentInFullTree" ->
    Select[quarticSymbols, StringContainsQ[treeText, #] &],
  "SelectedScalarQuarticDegreeCounts" ->
    (T3D2ScalarFactorDegree /@ realQuarticTexts),
  "SelectedQuarticTermPreviews" ->
    (StringTake[#, UpTo[200]] & /@ Take[realQuarticTexts,
       UpTo[6]])
|>;
Export[FileNameJoin[{dataDir, "eft1_after_F_realbasis_tensor_seed_audit.json"}],
  seedAudit, "RawJSON"];
If[seedAudit["Status"] =!= "Success",
  Print["ERROR: real-basis tensor seed incomplete: missing quartics=",
    missingQuartics, "; selected quartic terms=", Length[realQuarticRecords]];
  Print["See data/eft1_after_F_realbasis_tensor_seed_audit.json"];
  Exit[21]
];

Export[FileNameJoin[{dataDir, "eft1_after_F_realbasis_wilson_tensor_seed.json"}],
  <|"Status" -> "Success", "Model" -> "D_-2",
    "PhysicalScalarNames" -> {"NewScalar1", "T3S2R", "T3S2I"},
    "TreeWilsonTermCount" -> Length[realWilsonRecords],
    "TreeWilsonTerms" -> realWilsonRecords,
    "CGRegistry" -> realWilsonCGRegistry|>, "RawJSON"];
Export[FileNameJoin[{dataDir, "eft1_after_F_realbasis_quartic_tensor_seed.json"}],
  <|"Status" -> "Success", "Model" -> "D_-2",
    "PhysicalScalarNames" -> {"H", "NewScalar1", "T3S2R", "T3S2I"},
    "ScalarQuarticTermCount" -> Length[realQuarticRecords],
    "ScalarQuarticTerms" -> realQuarticRecords,
    "CGRegistry" -> realQuarticCGRegistry|>, "RawJSON"];
Print["[real EFT1 tensor seeds] ", Length[realWilsonRecords],
  " tree LLSS terms, ", Length[realQuarticRecords],
  " scalar quartic terms; CG registries ",
  Length[realWilsonCGRegistry], "/", Length[realQuarticCGRegistry]];

Print["SUCCESS: general D(-2) first threshold in real-S2 basis."];
Print["Real-scalar EFT1 Wilson RGE and scalar threshold remain separate tasks."];
Exit[0];
