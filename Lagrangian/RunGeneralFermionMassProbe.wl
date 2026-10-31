(* Workstream A: isolated heavy-fermion mass-basis DIAGNOSTIC, version 3.

   Project path: Lagrangian/RunGeneralFermionMassProbe.wl
   Usage (from project root):
     wolframscript -file Lagrangian/RunGeneralFermionMassProbe.wl \
       output/uv_fermion_mass_probe.json

   The previous version incorrectly passed an applied expression MFProbe[r]
   as a mass coupling label, producing DefineCoupling::CouplingLabel.
   This version passes only bare Symbols as labels. It independently tests
   the indexed-mass registration conventions supported by the installed
   Matchete package, and always exports a diagnostic JSON when Matchete loads.

   DOES NOT modify T3Fields.wl, call Match, or run RGBeta/numerical code.
   A check of a diagonal mass BASIS is not a certificate of general flavour
   matrices or threshold matching.
*)
ClearAll["Global`*"];

If[Length[$ScriptCommandLine] =!= 2,
  Print["Usage: wolframscript -file Lagrangian/RunGeneralFermionMassProbe.wl OUTPUT.json"];
  Exit[2]
];
output = ExpandFileName[$ScriptCommandLine[[2]]];
If[!DirectoryQ[DirectoryName[output]],
  CreateDirectory[DirectoryName[output], CreateIntermediateDirectories -> True]
];

(* This Matchete installation needs a service front end for package
   initialization. The production Workstream A audit successfully uses
   UsingFrontEnd[Needs["Matchete`"] ; True] under wolframscript.
   Using the kernel-only Needs call made version 2 fail before any probe. *)
loaded = CheckAbort[
  Check[UsingFrontEnd[Needs["Matchete`"] ; True], $Failed],
  $Aborted
];
If[loaded =!= True,
  Print["Matchete could not load in the Mathematica kernel: ", InputForm[loaded]];
  Export[output, <|
    "format" -> "t3-fermion-heavy-mass-probe-v3",
    "status" -> "MatcheteLoadFailed",
    "load_result" -> ToString[loaded, InputForm]
  |>, "RawJSON"];
  Exit[3]
];

ClearAll[T3MassProbeCase];
T3MassProbeCase[fermionKind_String, method_String] := Module[
  {sm, massLabel, fieldLabel, opts, defResult = "NotPredefined",
   fieldResult, free = $Failed, check = $Failed, massIndices = {},
   singleIndexQ, massSpec, outcome},

  ResetAll[];
  sm = CheckAbort[Check[LoadModel["SM"], $Failed], $Aborted];
  If[sm === $Failed || sm === $Aborted,
    Return[<|"fermion" -> fermionKind, "method" -> method,
      "status" -> "SMModelLoadFailed"|>]
  ];

  CheckAbort[
    Check[DefineFlavorIndex[MassProbeFlavor, 3,
      IndexAlphabet -> {"r", "s", "t"}], $Failed],
    $Aborted
  ];

  (* The mass coupling label is always an unapplied Symbol.
     Separate labels keep alternative trials independent. *)
  massLabel = Switch[{fermionKind, method},
    {"Majorana", "predefined"}, MFProbeMajPre,
    {"Majorana", "mass-third-list"}, MFProbeMajList,
    {"Majorana", "mass-third-type"}, MFProbeMajType,
    {"Majorana", "baseline"}, MFProbeMajBase,
    {"Dirac", "predefined"}, MFProbeDirPre,
    {"Dirac", "mass-third-list"}, MFProbeDirList,
    {"Dirac", "mass-third-type"}, MFProbeDirType,
    {"Dirac", "baseline"}, MFProbeDirBase
  ];
  fieldLabel = If[fermionKind === "Majorana", FProbeMajorana, FProbeDirac];

  (* For Matchete v0.6 this is the most useful tested hypothesis:
     register a real, one-flavour-index coupling FIRST and supply its
     bare symbol to DefineField. Documented versions allow reusing
     predeclared mass couplings when their definitions agree exactly. *)
  If[method === "predefined",
    defResult = CheckAbort[Check[
      Apply[DefineCoupling, {massLabel,
        Indices -> {MassProbeFlavor}, SelfConjugate -> True}]; True,
      $Failed], $Aborted]
  ];

  massSpec = Switch[method,
    "predefined" | "baseline", {Heavy, massLabel},
    "mass-third-list", {Heavy, massLabel, {MassProbeFlavor}},
    "mass-third-type", {Heavy, massLabel, MassProbeFlavor}
  ];
  opts = If[fermionKind === "Majorana",
    {SelfConjugate -> True, Indices -> {MassProbeFlavor}, Mass -> massSpec},
    {SelfConjugate -> False,
      Indices -> {SU2L[fund], MassProbeFlavor},
      Charges -> {U1Y[1/2]}, Mass -> massSpec}
  ];
  fieldResult = CheckAbort[Check[
    Apply[DefineField, Join[{fieldLabel, Fermion}, opts]]; True,
    $Failed], $Aborted];

  If[fieldResult === True,
    free = CheckAbort[Check[FreeLag[fieldLabel] // RelabelIndices,
      $Failed], $Aborted];
    If[free =!= $Failed && free =!= $Aborted,
      check = CheckAbort[Check[CheckLagrangian[(sm + free) // RelabelIndices],
        $Failed], $Aborted];
      massIndices = DeleteDuplicates @ Cases[free,
        Coupling[name_, indices_List, ___] /; name === massLabel :> indices,
        Infinity
      ]
    ]
  ];

  singleIndexQ = Length[massIndices] > 0 &&
    AllTrue[massIndices, Function[inds,
      Length[inds] === 1 && !FreeQ[inds[[1]], MassProbeFlavor]
    ]];
  outcome = Which[
    defResult === $Failed || defResult === $Aborted,
      "PredefinitionRejected",
    fieldResult =!= True, "FieldDefinitionRejected",
    free === $Failed || free === $Aborted, "FreeLagFailed",
    check =!= True, "LagrangianCheckFailed",
    singleIndexQ, "PassedIndexedMassBasis",
    True, "MassNotSingleFlavorIndexed"
  ];

  Print["Mass probe ", fermionKind, " / ", method, ": ", outcome];
  <|
    "fermion" -> fermionKind,
    "method" -> method,
    "status" -> outcome,
    "predefinition_result" -> ToString[defResult, InputForm],
    "field_definition_result" -> ToString[fieldResult, InputForm],
    "mass_option_inputform" -> ToString[massSpec, InputForm],
    "mass_indices_inputform" -> (ToString[#, InputForm] & /@ massIndices),
    "check_lagrangian_inputform" -> ToString[check, InputForm],
    "free_lagrangian_inputform" -> ToString[free, InputForm]
  |>
];

methods = {"predefined", "mass-third-list", "mass-third-type", "baseline"};
(* Match the successful fermionic audit: execute Matchete model operations
   under a service front end as well as package initialization. *)
results = CheckAbort[
  UsingFrontEnd[
    Flatten[Table[T3MassProbeCase[kind, method],
      {kind, {"Majorana", "Dirac"}}, {method, methods}], 1]
  ],
  $Aborted
];
If[!ListQ[results],
  Export[output, <|
    "format" -> "t3-fermion-heavy-mass-probe-v3",
    "status" -> "ProbeAborted",
    "result_inputform" -> ToString[results, InputForm]
  |>, "RawJSON"];
  Exit[5]
];
majoranaPass = AnyTrue[results,
  Lookup[#, "fermion", ""] === "Majorana" &&
    Lookup[#, "status", ""] === "PassedIndexedMassBasis" &];
diracPass = AnyTrue[results,
  Lookup[#, "fermion", ""] === "Dirac" &&
    Lookup[#, "status", ""] === "PassedIndexedMassBasis" &];
status = If[majoranaPass && diracPass, "Passed", "RequiresReview"];

report = <|
  "format" -> "t3-fermion-heavy-mass-probe-v3",
  "status" -> status,
  "majorana_indexed_mass_basis_passed" -> majoranaPass,
  "dirac_indexed_mass_basis_passed" -> diracPass,
  "cases" -> results,
  "source_files_modified" -> False,
  "matching_performed" -> False,
  "full_flavour_mass_generality_certified" -> False,
  "notes" -> {
    "This file tests candidate Mass option conventions, not a production change.",
    "A successful one-index coupling is a diagonal mass-basis check only.",
    "Predefined coupling and alternative Mass option forms are tested independently.",
    "Mass indices are detected by a structural Coupling pattern; inspect free_lagrangian_inputform when no indexed mass is found.",
    "See console log and each case's mass_option_inputform if definitions fail."
  }
|>;
Export[output, report, "RawJSON"];
Print["Fermion mass probe v3: ", status, " -> ", output];
If[status =!= "Passed", Exit[6]];
