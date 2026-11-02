(* Independent B_-1 or C_-1 native SU(2) Yukawa CG extractor.
   One model per fresh wolframscript kernel: do not ResetAll/LoadModel
   for two models in the same process. All Matchete operations requiring a
   frontend occur within UsingFrontEnd, not only Needs["Matchete`"].
   No matching, no release modification, no physical phase certification.
   Usage: wolframscript -file PATH B_-1|C_-1
*)
ClearAll["Global`*"];
args = Rest[$ScriptCommandLine];
If[Length[args] =!= 1 || !MemberQ[{"B_-1", "C_-1"}, First[args]],
  Print["Usage: wolframscript -file native_bc_yukawa_cg_probe.wl B_-1|C_-1"];
  Exit[2]
];
modelName = First[args];
probeFile = ExpandFileName[$InputFileName];
probeDir = DirectoryName[probeFile];
root = Nest[DirectoryName, probeFile, 4];
sourceDir = FileNameJoin[{root, "Lagrangian"}];
outputDir = FileNameJoin[{probeDir, "normalization_probes"}];
If[!DirectoryQ[outputDir], CreateDirectory[outputDir, CreateIntermediateDirectories -> True]];
modelTag = StringReplace[modelName, {"-" -> "m", "+" -> "p"}];
reportPath = FileNameJoin[{outputDir,
  "bc_native_yukawa_cg_" <> modelTag <> "_" <>
    StringReplace[CreateUUID[], "-" -> ""] <> ".json"}];
ClearAll[finish, fail, dense, complexArray, abortProbe, runModel];
finish[data_Association, exit_Integer] := Module[{payload},
  payload = Join[<|
    "format" -> "t3-bc-independent-yukawa-cg-components-v1",
    "model" -> modelName,
    "master_weyl_phase_certified" -> False,
    "native_BC_full_exchange_amplitude_computed" -> False,
    "previous_release_changed" -> False|>, data];
  If[Quiet[Check[Export[reportPath, payload, "RawJSON"], $Failed]] === $Failed,
    Print["BC_NATIVE_CG_EXPORT_FAILED=", reportPath]; Exit[7]];
  Print["BC_NATIVE_CG_REPORT=", reportPath];
  Print["BC_NATIVE_CG_STATUS=", Lookup[data, "status", "Unknown"]];
  Exit[exit]
];
fail[stage_String, why_String] := finish[
  <|"status" -> "Blocked", "stage" -> stage, "reason" -> why|>, 2];
abortProbe[stage_String, why_String] :=
  Throw[<|"status" -> "Blocked", "stage" -> stage, "reason" -> why|>, "T3BCStop"];
dense[t_] := Normal[t];
complexArray[name_String, t_, reps_] := Module[{a = dense[t], n},
  n = N[Flatten[a], 17];
  <|"name" -> name, "dimensions" -> Dimensions[a],
    "representations_inputform" -> ToString[InputForm[reps]],
    "real_flat" -> (Re /@ n), "imag_flat" -> (Im /@ n)|>
];

Print["BC_NATIVE_CG_STAGE=Matchete-load"];
(* Do not put Check around Needs: harmless startup messages can make Check
   return $Failed although the Matchete package loaded successfully. *)
loaded = CheckAbort[UsingFrontEnd[Needs["Matchete`"] ; True], $Aborted];
If[loaded =!= True, fail["load", "Matchete load did not return True"]];
Get[FileNameJoin[{sourceDir, "model", "LagrangianBuilder.wl"}]];
If[Names["*`PhysicalSU2InvariantBasis"] === {},
  fail["load", "PhysicalSU2InvariantBasis not available"]];

runModel[] := Module[
  {dF, y1, y2, old, t1, t2, tOld, eps, projected,
   den, alpha, expected, records},
  dF = If[modelName === "B_-1", 1, 3];
  Print["BC_NATIVE_CG_STAGE=", modelName, "_ResetAll"];
  ResetAll[];
  Print["BC_NATIVE_CG_STAGE=", modelName, "_LoadModel"];
  If[LoadModel["SM"] === $Failed,
    abortProbe[modelName, "Cannot initialize Standard Model"]];
  Print["BC_NATIVE_CG_STAGE=", modelName, "_Representations"];
  If[EnsureBSMRepresentation[2] =!= True,
    abortProbe[modelName, "Cannot register scalar doublet"]];
  If[dF === 3 && EnsureBSMRepresentation[3] =!= True,
    abortProbe[modelName, "Cannot register fermion triplet"]];
  Print["BC_NATIVE_CG_STAGE=", modelName, "_InvariantBases"];
  y1 = PhysicalSU2InvariantBasis[{2, dF, 2}, {True, True, False}, {}, {1}];
  y2 = PhysicalSU2InvariantBasis[{2, dF, 2}, {True, False, False}, {}, {1}];
  old = PhysicalSU2InvariantBasis[{2, dF, 2}, {True, False, True}, {}, {1}];
  If[!AllTrue[{y1, y2, old}, AssociationQ] ||
     !AllTrue[{y1, y2, old}, Lookup[#, "PhysicalCount", 0] === 1 &],
    abortProbe[modelName, "One native SU2 Yukawa invariant not found"]];
  t1 = dense[First[y1["Basis"]]];
  t2 = dense[First[y2["Basis"]]];
  tOld = dense[First[old["Basis"]]];
  If[Dimensions[t1] =!= If[dF === 1, {2, 2}, {2, 3, 2}] ||
     Dimensions[t2] =!= Dimensions[t1] ||
     Dimensions[tOld] =!= Dimensions[t1],
    abortProbe[modelName, "Native Yukawa tensor dimensions unexpected"]];
  eps = {{0, 1}, {-1, 0}};
  projected = If[dF === 1,
    Table[Sum[tOld[[i, m]] eps[[m, k]], {m, 2}], {i, 2}, {k, 2}],
    Table[Sum[tOld[[i, j, m]] eps[[m, k]], {m, 2}],
      {i, 2}, {j, 3}, {k, 2}]
  ];
  den = Total[Flatten[Conjugate[t2]] Flatten[t2]];
  If[TrueQ[PossibleZeroQ[den]],
    abortProbe[modelName, "Zero Yukawa2 norm"]];
  alpha = FullSimplify[
    Total[Flatten[Conjugate[t2]] Flatten[projected]]/den];
  If[!TrueQ[FullSimplify[projected == alpha t2]],
    abortProbe[modelName, "Original S2 epsilon transport is not proportional to native Yukawa2"]];
  expected = If[modelName === "B_-1", -1, 1];
  If[!TrueQ[FullSimplify[alpha == expected]],
    abortProbe[modelName, "Yukawa2 transport phase differs from existing Matchete source"]];
  Print["BC_NATIVE_CG_STAGE=", modelName, "_ExportComponents"];
  records = <|
    "fermion_su2_dimension" -> dF,
    "native_yukawa1_cg" -> complexArray["T3MENativeY1CG", t1, y1["CGReps"]],
    "native_yukawa2_cg" -> complexArray["T3MENativeY2CG", t2, y2["CGReps"]],
    "old_y2_epsilon_transport_alpha_inputform" -> ToString[InputForm[alpha]],
    "uv_y2_explicit_prefactor_inputform" -> ToString[InputForm[-alpha]],
    "native_y1_to_new_basis_alpha_inputform" -> "1",
    "source_construction" -> "Same PhysicalSU2InvariantBasis calls and -alpha y2 prefactor as RunGeneralT3FermionFirstEFT1Tree.wl",
    "matchete_master_spinor_phase_certified" -> False
  |>;
  records
];

(* Crucial: ResetAll, LoadModel, DefineRepresentation, InvariantTensors and
   their helpers run INSIDE UsingFrontEnd, not only the package startup. *)
Print["BC_NATIVE_CG_STAGE=", modelName, "_UsingFrontEnd"];
result = Catch[CheckAbort[UsingFrontEnd[runModel[]], $Aborted], "T3BCStop"];
If[result === $Aborted,
  fail[modelName, "Matchete computation aborted inside UsingFrontEnd"]];
If[AssociationQ[result] && Lookup[result, "status", ""] === "Blocked",
  fail[Lookup[result, "stage", modelName], Lookup[result, "reason", "Unknown failure"]]];
If[!AssociationQ[result] || !KeyExistsQ[result, "native_yukawa1_cg"],
  fail[modelName, "No complete native Yukawa tensor record returned"]];
finish[<|
  "status" -> "IndependentNativeYukawaCGComponentsExtracted_SingleModel",
  "record" -> result,
  "native_cg_product_computed_in_python" -> False|>, 0];
