(*
  Read-only A(-4) Higgs-triplet portal CG diagnostic.
  Uses the project's full-general Matchete Lagrangian builder, but DOES NOT
  run matching, change a coupling, or claim Matchete/RGBeta basis equivalence.
  Generates one output JSON; refuses to overwrite an existing directory.

  wolframscript -file Lagrangian/CheckAminus4PortalCGSymmetry.wl OUT_DIR
*)
ClearAll["Global`*"];
args = Rest[$ScriptCommandLine];
If[Length[args] =!= 1,
  Print["Usage: wolframscript -file Lagrangian/CheckAminus4PortalCGSymmetry.wl OUT_DIR"];
  Exit[2]
];
outDir = ExpandFileName[First[args]];
If[FileExistsQ[outDir] || DirectoryQ[outDir],
  Print["Output already exists; refusing to overwrite: ", outDir]; Exit[3]
];
flags = {"T3_EXPERIMENTAL_GENERAL_FERMIONS",
  "T3_EXPERIMENTAL_GENERAL_QUARTICS",
  "T3_EXPERIMENTAL_GENERAL_QUADRATICS",
  "T3_EXPERIMENTAL_FLAVOR_MASSES"};
If[!And @@ (Environment[#] === "1" & /@ flags),
  Print["The four general-Lagrangian environment variables must each equal 1"];
  Exit[4]
];
sourceDir = DirectoryName @ ExpandFileName[$InputFileName];
Needs["Matchete`"];
Get[FileNameJoin[{sourceDir, "model", "T3ModelCatalog.wl"}]];
Get[FileNameJoin[{sourceDir, "model", "LagrangianBuilder.wl"}]];
model = T3ModelFromClass["A", -4, 3];
If[!AssociationQ[model] || SetT3HeavyFields[{"F"}] === $Failed,
  Print["A(-4) model initialization failed"]; Exit[5]
];
built = CheckAbort[UsingFrontEnd[BuildT3Lagrangian[model]], $Aborted];
If[!AssociationQ[built] || Lookup[built, "Status", ""] =!= "Success" ||
    !TrueQ[Lookup[built, "FullValidation", False]] ||
    !TrueQ[Lookup[built, "AllScalarExtensionsRequested", False]] ||
    !And @@ (MemberQ[Lookup[built, "AllowedInteractions", {}], #] & /@
      {"HiggsPortal2Inv1", "HiggsPortal2Inv2"}),
  Print["The full-general model with both Higgs-triplet portals did not validate"];
  Exit[6]
];

rawCG1 = Quiet @ Check[
  Matchete`CGManipulations`PackagePrivate`$CGtensors[T3HiggsPortal2CGInv1],
  $Failed];
rawCG2 = Quiet @ Check[
  Matchete`CGManipulations`PackagePrivate`$CGtensors[T3HiggsPortal2CGInv2],
  $Failed];
If[rawCG1 === $Failed || rawCG2 === $Failed,
  Print["Native CG tensor registry is unavailable"]; Exit[7]
];
cg1 = N[Normal[rawCG1], 30];
cg2 = N[Normal[rawCG2], 30];
If[!ArrayQ[cg1] || !ArrayQ[cg2] ||
    Dimensions[cg1] =!= {2, 2, 3, 3} ||
    Dimensions[cg2] =!= {2, 2, 3, 3},
  Print["Unexpected native CG dimensions: ", Dimensions[cg1], " / ",
    Dimensions[cg2]]; Exit[8]
];

ClearAll[maxAbs, exchangedScalar, exchangedHermitian, directContraction, asComplex];
maxAbs[t_] := N[Max[Abs[Flatten[N[t, 25]]]], 16];
exchangedScalar[t_] := Transpose[t, {1, 2, 4, 3}];
(* This is only a component-level H/S hermiticity diagnostic, not a
   Matchete CG/index orientation certificate. *)
exchangedHermitian[t_] := Conjugate[Transpose[t, {2, 1, 4, 3}]];
asComplex[z_] := <|"real" -> N[Re[z], 16], "imag" -> N[Im[z], 16]|>;
directContraction[t_, h_List, s_List] := Sum[
  Conjugate[h[[i]]] h[[j]] Conjugate[s[[a]]] s[[b]] t[[i, j, a, b]],
  {i, 1, 2}, {j, 1, 2}, {a, 1, 3}, {b, 1, 3}];

hVectors = {{1 + 2 I, 2 - I}, {2 - I, 1 + 3 I}, {1 + I, -2 + I}};
sVectors = {{1 + I, 2 - 3 I, -2 + I},
  {1 - 2 I, -2 + I, 3 + I}, {2 + 3 I, 1 + I, -1 - 2 I}};
probes = Table[Module[{a1, a2},
    a1 = N[directContraction[cg1, hVectors[[k]], sVectors[[k]]], 16];
    a2 = N[directContraction[cg2, hVectors[[k]], sVectors[[k]]], 16];
    <|
      "probe_index" -> k,
      "inv1_unsymmetrized_contraction" -> asComplex[a1],
      "inv1_plus_hc_for_real_coupling" -> N[2 Re[a1], 16],
      "inv2_unsymmetrized_contraction" -> asComplex[a2],
      "inv2_plus_hc_for_real_coupling" -> N[2 Re[a2], 16]
    |>
  ], {k, Length[hVectors]}];

cg1Norm = maxAbs[cg1];
cg2Norm = maxAbs[cg2];
cg1ScalarSym = maxAbs[(cg1 + exchangedScalar[cg1])/2];
cg2ScalarSym = maxAbs[(cg2 + exchangedScalar[cg2])/2];
cg1ScalarAnti = maxAbs[(cg1 - exchangedScalar[cg1])/2];
cg2ScalarAnti = maxAbs[(cg2 - exchangedScalar[cg2])/2];
threshold = 10^-12;
report = <|
  "model" -> "A_-4",
  "scope" -> "Native CG tensor scalar-slot symmetry; read-only general Lagrangian",
  "full_general_build_success" -> True,
  "native_tensor_dimensions" -> {2, 2, 3, 3},
  "two_cg_tensors_linear_rank" -> MatrixRank[{Flatten[cg1], Flatten[cg2]}],
  "inv1" -> <|
    "max_abs_tensor" -> cg1Norm,
    "max_abs_symmetric_scalar_slot_part" -> cg1ScalarSym,
    "max_abs_antisymmetric_scalar_slot_part" -> cg1ScalarAnti,
    "purely_antisymmetric_scalar_slot_pair" -> (cg1ScalarSym < threshold && cg1Norm > threshold),
    "max_abs_hermiticity_difference" -> maxAbs[cg1 - exchangedHermitian[cg1]],
    "max_abs_antihermiticity_sum" -> maxAbs[cg1 + exchangedHermitian[cg1]]
  |>,
  "inv2" -> <|
    "max_abs_tensor" -> cg2Norm,
    "max_abs_symmetric_scalar_slot_part" -> cg2ScalarSym,
    "max_abs_antisymmetric_scalar_slot_part" -> cg2ScalarAnti,
    "purely_symmetric_scalar_slot_pair" -> (cg2ScalarAnti < threshold && cg2Norm > threshold),
    "max_abs_hermiticity_difference" -> maxAbs[cg2 - exchangedHermitian[cg2]],
    "max_abs_antihermiticity_sum" -> maxAbs[cg2 + exchangedHermitian[cg2]]
  |>,
  "independent_complex_field_probes" -> probes,
  "interpretation" -> "The uploaded raw Inv1 expression has 4 terms equal, after dummy-index relabeling, to the symmetric scalar-slot projection (T_ab + T_ba)/2. If Inv1 is antisymmetric there, that raw expression is zero, although T_ab itself need not be zero. The direct-component probes are diagnostics; a physical coupling dictionary and Matchete plus-h.c. phase still require certification.",
  "physical_matchete_rgbeta_mapping_certified" -> False,
  "physical_C5_certified" -> False
|>;
report["status"] = Which[
  TrueQ[Lookup[report["inv1"], "purely_antisymmetric_scalar_slot_pair", False]],
    "Inv1AntisymmetricTripletSlots_RawSymmetrizationProjectsItOut",
  cg1ScalarSym > threshold, "Inv1HasNonzeroScalarSymmetricPart_InvestigateGreens",
  True, "Inconclusive"
];
CreateDirectory[outDir, CreateIntermediateDirectories -> True];
path = FileNameJoin[{outDir, "portal_cg_symmetry.json"}];
Export[path, report, "RawJSON"];
Print["CG tensor rank: ", report["two_cg_tensors_linear_rank"]];
Print["Inv1 scalar-slot symmetric part (max abs): ", cg1ScalarSym];
Print["Inv1 scalar-slot antisymmetric part (max abs): ", cg1ScalarAnti];
Print["Inv2 scalar-slot symmetric part (max abs): ", cg2ScalarSym];
Print["Status: ", report["status"]];
Print["Report: ", path];
Exit[0];
