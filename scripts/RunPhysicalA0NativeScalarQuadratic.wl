(* A_+0: isolated live-Matchete scalar quadratic and CG witness.
   Does NOT call Match or RGBeta, perform EFT, or overwrite evidence.
   Usage (repository root):
     wolframscript -file scripts/RunPhysicalA0NativeScalarQuadratic.wl OUTPUT.json
*)
ClearAll["Global`*"];
$aqArgs = Rest[$ScriptCommandLine];
If[Length[$aqArgs] =!= 1, Print["Usage: OUTPUT.json"]; Exit[2]];
$aqOut = ExpandFileName[First[$aqArgs]];
If[FileExistsQ[$aqOut], Print["BLOCKED: refusing to overwrite: ", $aqOut]; Exit[2]];
$aqRoot = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
$aqFlags = {"T3_EXPERIMENTAL_GENERAL_FERMIONS", "T3_EXPERIMENTAL_GENERAL_QUADRATICS",
  "T3_EXPERIMENTAL_GENERAL_QUARTICS", "T3_EXPERIMENTAL_FLAVOR_MASSES"};
$aqPriorFlags = AssociationMap[Function[x, With[{v = Environment[x]}, If[StringQ[v], v, "UNSET"]]], $aqFlags];
SetEnvironment[Thread[$aqFlags -> ConstantArray["0", Length[$aqFlags]]]];
If[!AllTrue[$aqFlags, Environment[#] === "0" &],
  Print["BLOCKED: restricted flags could not be isolated"]; Exit[3]];
If[CheckAbort[Check[UsingFrontEnd[Needs["Matchete`"]]; True, $Failed], $Aborted] =!= True,
  Print["BLOCKED: installed Matchete could not be loaded"]; Exit[4]];
Get[FileNameJoin[{$aqRoot, "Lagrangian", "model", "T3ModelCatalog.wl"}]];
Get[FileNameJoin[{$aqRoot, "Lagrangian", "model", "LagrangianBuilder.wl"}]];

ClearAll[AQFail, AQText, AQTerm, AQTerms, AQNumber, AQPick, AQGetTensor,
  AQClose, AQPayload, AQSourceHashes, AQRun];
AQFail[msg_String] := Throw[<|"status" -> "blocked", "reason" -> msg|>, "AQAbort"];
AQText[x_] := ToString[x, InputForm];
AQTerm[t_, tokens_List] := Module[{s = AQText[t]}, AllTrue[tokens, StringContainsQ[s, #] &]];
AQTerms[e_, tokens_List] := Select[If[Head[Expand[e]] === Plus, List @@ Expand[e], {Expand[e]}],
  AQTerm[#, tokens] &];
(* Only literal top-level integer/rational prefactors of an isolated native monomial. *)
AQNumber[t_] := Module[{f = If[Head[t] === Times, List @@ t, {t}]},
  Times @@ Select[f, NumberQ]];
AQPick[e_, tokens_List, n_Integer, expected_List] := Module[{ts, cs},
  ts = AQTerms[e, tokens];
  If[Length[ts] =!= n, AQFail["Unexpected native term count for " <> AQText[tokens] <>
     ": " <> ToString[Length[ts]]]];
  cs = AQNumber /@ ts;
  If[cs =!= expected, AQFail["Unexpected coefficient(s) for " <> AQText[tokens] <>
     ": " <> AQText[cs]]];
  <|"term_count" -> Length[ts], "coefficients" -> (AQText /@ cs),
    "native_terms_inputform" -> (AQText /@ ts)|>
];
AQGetTensor[name_Symbol, dims_List] := Module[{raw, arr},
  raw = CheckAbort[Quiet @ Check[GetCGTensor[name], $Failed], $Aborted];
  If[raw === $Failed || raw === $Aborted,
    AQFail["GetCGTensor failed for " <> SymbolName[name]]];
  arr = If[Head[raw] === SparseArray, Normal[raw], raw];
  If[Dimensions[arr] =!= dims ||
     !AllTrue[Flatten[arr], TrueQ[NumericQ[FullSimplify[#]]] &],
    AQFail["Native tensor invalid for " <> SymbolName[name]]];
  arr
];
AQClose[a_, b_] := TrueQ[FullSimplify[a - b == 0]];
AQPayload[arr_] := <|
  "shape" -> Dimensions[arr],
  "real_flat" -> (N[Re[#], 16] & /@ Flatten[arr]),
  "imag_flat" -> (N[Im[#], 16] & /@ Flatten[arr])|>;
AQSourceHashes[] := Module[{paths, f},
  paths = {"Lagrangian/model/LagrangianBuilder.wl",
    "Lagrangian/model/T3ModelCatalog.wl", "Lagrangian/model/T3Fields.wl",
    "Lagrangian/model/SU2Invariants.wl", "Lagrangian/interactions/T3Topology.wl",
    "Lagrangian/interactions/ScalarPotential.wl"};
  Association@Table[
    f = FileNameJoin[Prepend[FileNameSplit[path], $aqRoot]];
    If[!FileExistsQ[f], AQFail["Source missing: " <> path]];
    path -> IntegerString[FileHash[f, "SHA256"], 16, 64],
    {path, paths}]
];
AQRun[] := Module[
 {model, built, lint, lf, action, mass1, mass2, kinetic1, kinetic2,
  c1, c2, cm, p1, p21, p22, pr1, pr21, pr22, mvec, yvec, active, basis,
  source, gith, cert, checks},
 model = T3ModelFromClass["A", 0];
 If[!AssociationQ[model], AQFail["A_+0 catalog construction failed"]];
 built = CheckAbort[Check[UsingFrontEnd[BuildT3Lagrangian[model]], $Failed], $Aborted];
 If[!AssociationQ[built] || Lookup[built, "Status", ""] =!= "Success" ||
    !TrueQ[Lookup[built, "FullValidation", False]],
    AQFail["A_+0 live LUV did not pass Matchete CheckLagrangian"]];
 lint = built["LInt"]; lf = built["LFree"];
 (* Verify native source monomials, including H.c. and 1/2 PlusHc.
    Read their actual top-level coefficients from the assembled action. *)
 action = <|
  "scalar1_mass" -> AQPick[lf, {"Coupling[MS1,", "^2", "NewScalar1"}, 1, {-1}],
  "scalar2_mass" -> AQPick[lf, {"Coupling[MS2,", "^2", "NewScalar2"}, 1, {-1}],
  "scalar1_portal" -> AQPick[lint, {"Coupling[lambdaH1Inv1,", "CG[T3HiggsPortal1CGInv1,"}, 1, {1}],
  "scalar2_portal_inv1" -> AQPick[lint, {"Coupling[lambdaH2Inv1,", "CG[T3HiggsPortal2CGInv1,"}, 2, {1/2, 1/2}],
  "scalar2_portal_inv2" -> AQPick[lint, {"Coupling[lambdaH2Inv2,", "CG[T3HiggsPortal2CGInv2,"}, 2, {1/2, 1/2}],
  "mixing_original" -> AQPick[lint, {"Coupling[lambdaT3,", "CG[T3MixCG,"}, 1, {1}],
  "mixing_hc" -> AQPick[lint, {"Bar[Coupling[lambdaT3,", "CG[Bar[T3MixCG],"}, 1, {1}]
 |>;
 (* Canonical scalar kinetic terms in FreeLag, not the fermion or gauge terms. *)
 kinetic1 = Select[AQTerms[lf, {"NewScalar1", "Lorentz"}],
   !StringContainsQ[AQText[#], "Coupling[MS1,"] &];
 kinetic2 = Select[AQTerms[lf, {"NewScalar2", "Lorentz"}],
   !StringContainsQ[AQText[#], "Coupling[MS2,"] &];
 If[Length[kinetic1] =!= 1 || Length[kinetic2] =!= 1 ||
     AQNumber[First[kinetic1]] =!= 1 || AQNumber[First[kinetic2]] =!= 1,
   AQFail["Cannot verify canonical native scalar kinetic terms"]];
 c1 = AQGetTensor[T3Y1CG, {2, 2}];
 c2 = AQGetTensor[T3Y2CG, {2, 2, 3}];
 cm = AQGetTensor[T3MixCG, {2, 2, 3}];
 p1 = AQGetTensor[T3HiggsPortal1CGInv1, {2, 2}];
 p21 = AQGetTensor[T3HiggsPortal2CGInv1, {2, 2, 3, 3}];
 p22 = AQGetTensor[T3HiggsPortal2CGInv2, {2, 2, 3, 3}];
 (* Restrict H=(0,v/sqrt(2)) directly in native, unrotated Matchete axes. *)
 pr1 = p1[[2, 2]];
 pr21 = (p21[[2, 2, All, All]] +
   Conjugate[Transpose[p21[[2, 2, All, All]]]])/2;
 pr22 = (p22[[2, 2, All, All]] +
   Conjugate[Transpose[p22[[2, 2, All, All]]]])/2;
 mvec = cm[[2, 2, All]];
 yvec = c2[[1, 2, All]];
 checks = <|
   "source_A0_dimensions" -> (model["Scalar1", "SU2"] === 1 &&
     model["Scalar2", "SU2"] === 3 && model["Fermion", "SU2"] === 2),
   "neutral_fermion_axis" -> AQClose[c1[[1,2]], 1],
   "native_H0_portal1" -> AQClose[pr1, 1],
   "native_H0_portal2_inv1_zero" -> AllTrue[Flatten[pr21], AQClose[#,0] &],
   "native_H0_portal2_inv2_identity" -> And @@ Flatten[MapThread[AQClose, {pr22, IdentityMatrix[3]}, 2]],
   "native_mix_vector" -> And @@ MapThread[AQClose, {mvec,
      {-3^(-1/4), 0, I 3^(-1/4)}}],
   "native_neutral_y2_vector" -> And @@ MapThread[AQClose, {yvec,
      {3^(-1/4), 0, -I 3^(-1/4)}}]
 |>;
 If[!AllTrue[Values[checks], TrueQ],
   AQFail["Unexpected native source CG component/projection; stop before numerical loop"]];
 source = AQSourceHashes[];
 gith = Quiet @ Check[RunProcess[{"git", "-C", $aqRoot, "rev-parse", "HEAD"},
      "StandardOutput"], "UNAVAILABLE"];
 <|
  "format" -> "t3-b1-a0-live-native-scalar-quadratic-v1",
  "status" -> "NativeA0QuadraticActionAndRegisteredCGExtracted_IndependentLoopPending",
  "model" -> "A_+0", "native_matchete_checked_lagrangian" -> True,
  "source_git_head" -> StringTrim[gith], "source_sha256" -> source,
  "restricted_environment" -> <|"prior" -> $aqPriorFlags,
    "effective" -> AssociationMap[Environment, $aqFlags]|>,
  "source_field_hypercharges" -> <|"S1" -> "0", "S2" -> "1", "F" -> "1/2"|>,
  "scalar_free_kinetic_coefficients" -> {"1", "1"},
  "scalar_free_kinetic_inputform" -> (AQText /@ Join[kinetic1, kinetic2]),
  "assembled_native_operator_terms" -> action,
  "native_tensors" -> <|"T3Y1CG" -> AQPayload[c1], "T3Y2CG" -> AQPayload[c2],
    "T3MixCG" -> AQPayload[cm], "T3HiggsPortal1CGInv1" -> AQPayload[p1],
    "T3HiggsPortal2CGInv1" -> AQPayload[p21],
    "T3HiggsPortal2CGInv2" -> AQPayload[p22]|>,
  "higgs_vev_convention" -> "H[2]=v/sqrt(2); H[1]=0",
  "native_H0_projection" -> <|
    "scalar1_portal_inv1" -> AQPayload[{{pr1}}],
    "scalar2_portal_inv1" -> AQPayload[pr21],
    "scalar2_portal_inv2" -> AQPayload[pr22],
    "mixing_S2bar_vector" -> AQPayload[mvec],
    "y2_nu_Fneutral_S2bar_vector" -> AQPayload[yvec]|>,
  "matrix_convention" -> "L_quadratic=-phi_dagger A phi; A_S2,S1=-(lambdaT3*v^2/2)*T3MixCG[2,2,S2]",
  "native_component_checks" -> checks,
  "source_determined_operator_coefficients" -> True,
  "benchmark_to_native_portal_coupling_values_fixed_by_source" -> False,
  "absolute_four_spinor_GammaCC_phase_measured" -> False,
  "physical_source_certified" -> False,
  "matching_or_running_performed" -> False|>
];
$aqResult = CheckAbort[Catch[AQRun[], "AQAbort"], $Aborted];
If[!AssociationQ[$aqResult] || Lookup[$aqResult, "status", ""] =!=
  "NativeA0QuadraticActionAndRegisteredCGExtracted_IndependentLoopPending",
  Print["BLOCKED: ", If[AssociationQ[$aqResult], Lookup[$aqResult, "reason", "native error"],
    "Matchete run aborted"]]; Exit[5]];
If[!DirectoryQ[DirectoryName[$aqOut]],
  CreateDirectory[DirectoryName[$aqOut], CreateIntermediateDirectories -> True]];
If[Quiet @ Check[Export[$aqOut, $aqResult, "RawJSON"], $Failed] === $Failed,
  Print["BLOCKED: native JSON export failed"]; Exit[6]];
Print["A_+0 native scalar kinetic/mass + portal/mix/yukawa tensors: EXTRACTED"];
Print["Native vev projection: S1=1, S2Inv1=0, S2Inv2=Identity[3]"];
Print["Witness: ", $aqOut];
Print["Physical source-certified models: 0/16; benchmark-to-native mapping still conditional"];
Exit[0];
