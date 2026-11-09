(* Batch 1: native all-model Higgs portal CG Hermiticity census.
   This probes the REAL corrected ScalarPotential.wl already in the local tree.
   Nothing in GitHub or in production source is modified.

   wolframscript -file scripts/RunT3Batch1PortalAudit.wl OUTPUT.json
*)
ClearAll["Global`*"];
args = Rest[$ScriptCommandLine];
If[Length[args] =!= 1,
  Print["Usage: RunT3Batch1PortalAudit.wl OUTPUT.json"]; Exit[2]];
outPath = ExpandFileName[First[args]];
root = DirectoryName @ DirectoryName @ ExpandFileName[$InputFileName];
Needs["Matchete`"];
Get[FileNameJoin[{root, "Lagrangian", "model", "T3ModelCatalog.wl"}]];
Get[FileNameJoin[{root, "Lagrangian", "model", "LagrangianBuilder.wl"}]];

cases = Join[
  ("A" -> # & /@ {-4, -2, 0}),
  ("B" -> # & /@ {-3, -1, 1}),
  ("C" -> # & /@ {-3, -1, 1}),
  ("D" -> # & /@ {-2, 0, 2}),
  ("E" -> # & /@ {-4, -2, 0, 2})
];

T3B1Label[cl_, a_] := cl <> "_" <> If[a < 0, "-", "+"] <> ToString[Abs[a]];

T3B1Probe[cl_, a_] := Module[
  {label, model, built, rows = {}, d, basis, n, symbol, phase, tensors,
   cgResult, portalCount, which},
  label = T3B1Label[cl, a];
  Print["B1 portal Hermiticity census ", label];
  model = T3ModelFromClass[cl, a, 3];
  If[!AssociationQ[model],
    Return[<|"model" -> label, "status" -> "ModelCatalogFailed"|>]
  ];
  built = CheckAbort[Check[UsingFrontEnd[BuildT3Lagrangian[model]], $Failed], $Aborted];
  If[!AssociationQ[built],
    Return[<|"model" -> label, "status" -> "NativeBuildFailed"|>]
  ];
  Do[
    d = model[If[which == 1, "Scalar1", "Scalar2"], "SU2"];
    basis = PhysicalSU2InvariantBasis[
      {2, 2, d, d}, {True, False, True, False}, {}, {1, 2}
    ];
    If[!AssociationQ[basis],
      Return[<|"model" -> label, "status" -> "PortalBasisFailed"|>]
    ];
    n = basis["PhysicalCount"];
    Do[
      symbol = Symbol["T3HiggsPortal" <> ToString[which] <>
        "CGInv" <> ToString[k]];
      phase = Quiet @ Check[
        T3HiggsPortalHermitianPhase[If[d === 1, None, symbol], d], $Failed];
      If[phase === $Failed,
        Return[<|"model" -> label, "status" -> "PortalHermiticityFailed",
           "field" -> which, "invariant" -> k|>]
      ];
      AppendTo[rows, <|
        "coupling" -> "lambdaH" <> ToString[which] <> "Inv" <> ToString[k],
        "cg" -> SymbolName[symbol],
        "dimension" -> d,
        "phase_inside_plus_hc" -> If[phase === I, "I", "1"],
        "native_hermiticity" -> If[phase === I, "AntiHermitian", "Hermitian"]
      |>],
      {k, n}
    ],
    {which, 1, 2}
  ];
  <|"model" -> label, "status" -> "NativePortalCGHermiticityEnumerated",
    "builder_status" -> Lookup[built, "Status", "Unknown"],
    "full_lagrangian_validation" -> TrueQ[Lookup[built,"FullValidation",False]],
    "portal_count" -> Length[rows],
    "portals" -> rows|>
];

If[FileExistsQ[outPath],
  Print["Output exists: use a new batch run-id, not an overwrite."]; Exit[6]
];
results = Association@Table[
  With[{cl = First[c], a = Last[c]},
    T3B1Label[cl, a] -> T3B1Probe[cl, a]
  ], {c, cases}
];
report = <|
  "format" -> "t3-batch1-native-portal-hermiticity-census-v1",
  "source" -> "CurrentLocalScalarPotentialNativeMatcheteCG",
  "model_count" -> Length[Keys[results]],
  "models" -> results,
  "does_not_certify_matchete_rgbeta_normalization" -> True,
  "does_not_certify_matched_EFT1_operator_survival" -> True
|>;
If[!DirectoryQ[DirectoryName[outPath]],
  CreateDirectory[DirectoryName[outPath],CreateIntermediateDirectories -> True]
];
If[Quiet @ Check[Export[outPath,report,"RawJSON"],$Failed] === $Failed,
  Print["Native portal audit export failed"]; Exit[7]
];
Print["All-model portal census: ", outPath];
