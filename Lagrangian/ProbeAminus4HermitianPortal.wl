(*
  A(-4) full-general Higgs-triplet portal phase trial.

  This is deliberately an ADDITIVE diagnostic. It overrides ONLY the two
  H-S2 portal candidate expressions in memory, never edits ScalarPotential.wl,
  saved Matchete source, UV/EFT1 beta payloads, or Wilson outputs.

  The original A(-4) T3HiggsPortal2CGInv1 is anti-Hermitian in the convention
  conjugate-transpose (H dagger/H and S2 dagger/S2).  A real lambda times
  1/2 PlusHc[O] projects it to zero.  This trial uses 1/2 PlusHc[I lambda O]
  after independently checking the original native CG hermiticity.

  This is NOT a Matchete/RGBeta normalization certificate.

  Usage:
    wolframscript -file Lagrangian/ProbeAminus4HermitianPortal.wl OUT_DIR
*)
ClearAll["Global`*"];
args = Rest[$ScriptCommandLine];
If[Length[args] =!= 1,
  Print["Usage: ProbeAminus4HermitianPortal.wl OUT_DIR"]; Exit[2]
];
outDir = ExpandFileName[First[args]];
If[FileExistsQ[outDir] || DirectoryQ[outDir],
  Print["Refusing to overwrite: ", outDir]; Exit[3]
];
flags = {"T3_EXPERIMENTAL_GENERAL_FERMIONS",
  "T3_EXPERIMENTAL_GENERAL_QUARTICS",
  "T3_EXPERIMENTAL_GENERAL_QUADRATICS",
  "T3_EXPERIMENTAL_FLAVOR_MASSES"};
If[!And @@ (Environment[#] === "1" & /@ flags),
  Print["All four general-Lagrangian switches must equal 1"]; Exit[4]
];
CreateDirectory[outDir, CreateIntermediateDirectories -> True];
sourceDir = DirectoryName @ ExpandFileName[$InputFileName];
Needs["Matchete`"];
Get[FileNameJoin[{sourceDir, "model", "T3ModelCatalog.wl"}]];
Get[FileNameJoin[{sourceDir, "model", "LagrangianBuilder.wl"}]];
Get[FileNameJoin[{sourceDir, "RunMatching.wl"}]];

ClearAll[writeReport, saveExpression, termsFor, termList, maxAbs];
writeReport[r_Association] := Export[
  FileNameJoin[{outDir, "hermitian_portal_trial.json"}], r, "RawJSON"];
saveExpression[name_String, expr_] := Module[{s, p},
  s = ToString[InputForm[expr]];
  p = FileNameJoin[{outDir, name <> ".inputform.txt"}];
  Export[p, s <> "\n", "Text"];
  <|"file" -> FileNameTake[p],
    "inv1_mentions" -> StringCount[s, "lambdaH2Inv1"],
    "inv2_mentions" -> StringCount[s, "lambdaH2Inv2"],
    "exactly_zero" -> TrueQ[expr === 0],
    "characters" -> StringLength[s]|>
];
termList[x_] := Module[{e = Expand[x]},
  If[Head[e] === Plus, List @@ e, {e}]
];
termsFor[x_, symbol_String] := Select[termList[x],
  StringContainsQ[ToString[InputForm[#]],
    "Coupling[" <> symbol <> ","] &
];
maxAbs[x_] := N[Max[Abs[Flatten[N[x, 25]]]], 16];

(* Only the full-general A(-4) experiment is altered, and only in this kernel.
   The one original complex scalar field remains Light at the F threshold. *)
Clear[BuildHiggsPortalInvariantCandidates];
BuildHiggsPortalInvariantCandidates[which_Integer, d_Integer?Positive] := Module[
  {family, indices, labels, cg, coupling, fieldProduct, baseName, phase,
    original, conjugateTranspose, norm, hermResidual, antiResidual},
  baseName = "lambdaH" <> ToString[which];
  family = DefineSU2InvariantFamily[
    "T3HiggsPortal" <> ToString[which] <> "CG", baseName,
    {2, 2, d, d}, {True, False, True, False}, {}, {1, 2}];
  If[family === $Failed, Print["Higgs portal CG construction failed"]; Exit[10]];
  Table[
    indices = Array[Unique["p"] &, 4];
    labels = indices[[family["Keep"]]];
    cg = family["CGs"][[invariant]];
    coupling = family["Couplings"][[invariant]];
    phase = 1;
    If[which === 2 && d === 3,
      original = Quiet @ Check[
        Matchete`CGManipulations`PackagePrivate`$CGtensors[cg], $Failed];
      If[original === $Failed ||
          Dimensions[Normal[original]] =!= {2, 2, 3, 3},
        Print["Native H-S2 CG tensor unavailable; aborting"]; Exit[11]
      ];
      original = Normal[original];
      conjugateTranspose = Conjugate[
        Transpose[original, {2, 1, 4, 3}]];
      norm = maxAbs[original];
      hermResidual = maxAbs[original - conjugateTranspose];
      antiResidual = maxAbs[original + conjugateTranspose];
      If[!TrueQ[norm > 10^-12],
        Print["Vanishing H-S2 CG tensor; aborting"]; Exit[12]
      ];
      phase = Which[
        invariant === 1 && antiResidual < 10^-10 norm, I,
        invariant === 2 && hermResidual < 10^-10 norm, 1,
        True,
          Print["Unexpected H-S2 CG hermiticity; inv=", invariant,
            " herm=", hermResidual, " anti=", antiResidual]; Exit[13]
      ];
      Print["H-S2 Inv", invariant, " verified CG phase: ", phase,
        " (Hermitian residual=", hermResidual,
        ", anti-Hermitian residual=", antiResidual, ")"];
    ];
    fieldProduct =
      Bar[H[indices[[1]]]] H[indices[[2]]] *
      ScalarPotentialFieldValue[which, True, d, indices[[3]]] *
      ScalarPotentialFieldValue[which, False, d, indices[[4]]];
    <|
      "Name" -> "HiggsPortal" <> ToString[which] <> "Inv" <> ToString[invariant],
      "Class" -> "Portal",
      "Expression" -> 1/2 PlusHc[
        phase coupling[] fieldProduct If[cg === None, 1, CG[cg, labels]]
      ]
    |>,
    {invariant, family["Count"]}
  ]
];

model = T3ModelFromClass["A", -4, 3];
If[!AssociationQ[model] || SetT3HeavyFields[{"F"}] === $Failed,
  Print["A(-4) model or F-only threshold failed"]; Exit[14]
];
built = CheckAbort[UsingFrontEnd[BuildT3Lagrangian[model]], $Aborted];
report = <|
  "model" -> "A_-4",
  "scope" -> "In-memory trial of Hermitian H-S2 portal; no production source edits",
  "general_flags" -> AssociationMap[Environment, flags],
  "original_scalar_potential_file_modified" -> False,
  "trial_phase" -> "I multiplied inside PlusHc for H-S2 Inv1 only; native anti-Hermiticity checked",
  "matchete_rgbeta_basis_map_certified" -> False,
  "physical_C5_certified" -> False|>;
If[!AssociationQ[built] || Lookup[built, "Status", ""] =!= "Success" ||
    !TrueQ[Lookup[built, "FullValidation", False]] ||
    !TrueQ[Lookup[built, "AllScalarExtensionsRequested", False]] ||
    !And @@ (MemberQ[Lookup[built, "AllowedInteractions", {}], #] & /@
      {"HiggsPortal2Inv1", "HiggsPortal2Inv2"}),
  report["status"] = "Blocked_TrialGeneralUVBuildFailed";
  writeReport[report]; Print[report["status"]]; Exit[15]
];
report["general_UV_full_validation"] = True;
snapshots = <|"LUV" -> saveExpression["LUV_trial", built["LUV"]]|>;
report["snapshots"] = snapshots;
Print["Full-general UV Lagrangian validated with Hermitian portal trial."];

matching = CheckAbort[UsingFrontEnd[RunT3Matching[built["LUV"], 5, 0]],
  $Aborted];
If[!AssociationQ[matching] || Lookup[matching, "Status", ""] =!= "Success",
  report["status"] = "Blocked_TrialMatcheteMatchFailed";
  writeReport[report]; Print[report["status"]]; Exit[16]
];
Do[
  snapshots[name] = saveExpression[name <> "_trial", matching[name]],
  {name, {"RawEFT", "GreenEFT", "EOMEFT", "MatchedEFT"}}
];
raw1 = Expand[Total[termsFor[matching["RawEFT"], "lambdaH2Inv1"]]];
raw2 = Expand[Total[termsFor[matching["RawEFT"], "lambdaH2Inv2"]]];
snapshots["raw_inv1_only"] = saveExpression["RawEFT_inv1_trial", raw1];
snapshots["raw_inv2_only"] = saveExpression["RawEFT_inv2_trial", raw2];
report["raw_inv1_summands"] = Length[termsFor[matching["RawEFT"], "lambdaH2Inv1"]];
report["raw_inv2_summands"] = Length[termsFor[matching["RawEFT"], "lambdaH2Inv2"]];
If[TrueQ[raw1 === 0] || TrueQ[raw2 === 0],
  report["status"] = "Blocked_TrialRawPortalVanished";
  writeReport[report]; Print[report["status"]]; Exit[17]
];
Print["Testing isolated GreensSimplify on phase-corrected Inv1..."];
green1 = CheckAbort[UsingFrontEnd[GreensSimplify[raw1]], $Aborted];
If[MemberQ[{$Aborted, $Failed}, green1],
  report["status"] = "Blocked_TrialIsolatedGreensFailed";
  writeReport[report]; Print[report["status"]]; Exit[18]
];
snapshots["green_inv1_only"] =
  saveExpression["GreenEFT_inv1_trial", green1];
report["snapshots"] = snapshots;
report["isolated_inv1_survives_GreensSimplify"] = !TrueQ[green1 === 0];
report["full_matched_inv1_explicit"] =
  TrueQ[report["snapshots"]["MatchedEFT"]["inv1_mentions"] > 0];
report["full_matched_inv2_explicit"] =
  TrueQ[report["snapshots"]["MatchedEFT"]["inv2_mentions"] > 0];
report["status"] = Which[
  !TrueQ[report["isolated_inv1_survives_GreensSimplify"]],
    "Blocked_TrialInv1StillEliminated",
  !TrueQ[report["full_matched_inv1_explicit"]] ||
    !TrueQ[report["full_matched_inv2_explicit"]],
    "Blocked_TrialNotBothPortalsExplicitAtEFT1",
  True, "TrialBothPortalsSurvive_DoNotTreatAsPhysicalBasisCertificate"
];
writeReport[report];
Print["Status: ", report["status"]];
Print["Report: ", FileNameJoin[{outDir, "hermitian_portal_trial.json"}]];
Print["Final explicit portal mentions: Inv1=",
  report["snapshots"]["MatchedEFT"]["inv1_mentions"],
  " Inv2=", report["snapshots"]["MatchedEFT"]["inv2_mentions"]];
If[StringStartsQ[report["status"], "Blocked"], Exit[19]];
Exit[0];
