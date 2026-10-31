(* Workstream A — read-only LIVE Matchete A(0) charge-conjugation convention probe.
   Usage:
     wolframscript -file Lagrangian/RunT3YukawaSpinorPhaseNative.wl OUTPUT.json
   Requires existing Matchete installation and full A(0) build setup.
   A successful check establishes NATIVE ABSTRACT Dirac-algebra identities,
   not a numerical value for the conventional phase of GammaCC or a complete
   Matchete/RGBeta physical dictionary. No Match[] or RGEs are evaluated.
*)
ClearAll["Global`*"];
args = Rest[$ScriptCommandLine];
If[Length[args] =!= 1,
 Print["Usage: OUTPUT.json"]; Exit[2]];
out = ExpandFileName[First[args]];
dir = DirectoryName @ ExpandFileName[$InputFileName];
If[!And @@ (Environment[#] === "1" & /@ {
    "T3_EXPERIMENTAL_GENERAL_FERMIONS",
    "T3_EXPERIMENTAL_FLAVOR_MASSES",
    "T3_EXPERIMENTAL_GENERAL_QUARTICS"}) ||
    Environment["T3_EXPERIMENTAL_GENERAL_QUADRATICS"] === "1",
 Print["Required: fermions/flavour masses/quartics = 1, quadratics != 1"];
 Exit[2]];

step = "load-Matchete";
load = CheckAbort[Check[UsingFrontEnd[Needs["Matchete`"]; True], $Failed], $Aborted];
If[load =!= True, Print["Matchete load failed"]; Exit[3]];
Get[FileNameJoin[{dir, "model", "T3ModelCatalog.wl"}]];
Get[FileNameJoin[{dir, "model", "LagrangianBuilder.wl"}]];
ClearAll[SPFail, SPString, SPRun];
SPFail[msg_String] := Throw[<|"status" -> "NativeSpinorConventionBlocked",
     "failed_step" -> step, "reason" -> msg|>, "SpinorPhaseStop"];
SPString[x_] := ToString[x, InputForm];

SPRun[] := Module[{model, built, lf, r1, r2, c, barl, spin,
   byDefinition, diracC2, transposeC, flipL, flipR, checks,
   hc, hcExpected, lInt, ccName, result},
 step = "build-native-A_p0";
 model = T3ModelFromClass["A", 0];
 If[!AssociationQ[model], SPFail["T3ModelFromClass did not return an association"]];
 built = CheckAbort[Check[UsingFrontEnd[BuildT3Lagrangian[model]], $Failed], $Aborted];
 If[!AssociationQ[built] || Lookup[built, "Status", ""] =!= "Success" ||
    !TrueQ[Lookup[built, "FullValidation", False]] ||
    Lookup[built, "AdditionalYukawasIncluded", -1] =!= 1,
   SPFail["Native Matchete A(0) model build was not accepted"]];
 lInt = built["LInt"];
 step = "find-distinct-native-yukawas";
 If[FreeQ[lInt, y1] || FreeQ[lInt, y1Conj],
   SPFail["Native A(0) assembled interaction misses y1 or y1Conj"]];
 step = "construct-native-charge-conjugate-bilinear";
 lf = l[i, p];
 r1 = NewFermion[j, r];
 c = DiracProduct[GammaCC];
 barl = Bar[lf];
 spin = NCM[barl, CConj[r1]];
 byDefinition = NCM[barl, c, Transp[Bar[r1]]];
 hc = Bar[spin];
 hcExpected = Bar[byDefinition];
 step = "check-native-charge-conjugation";
 (* All identities are tested on the NATIVE Matchete objects after loading
    and defining the same fields used by the A(0) model. *)
 diracC2 = DiracProduct[GammaCC, GammaCC];
 transposeC = Transp[GammaCC];
 flipL = CConj[NCM[PL, r1]];
 flipR = NCM[PR, CConj[r1]];
 checks = <|
  "gamma_cc_square_minus_one" -> TrueQ[FullSimplify[diracC2 == -1]],
  "gamma_cc_transpose_antisymmetric" -> TrueQ[FullSimplify[transposeC == -GammaCC]],
  "cconj_field_definition_in_bilinear" -> TrueQ[FullSimplify[spin == byDefinition]],
  "cconj_flips_left_projector" -> TrueQ[FullSimplify[flipL == flipR]],
  "hermitian_conjugate_matches_explicit_cconj" -> TrueQ[FullSimplify[hc == hcExpected]],
  "distinct_native_y1_and_y1conj_present" -> (!FreeQ[lInt, y1] && !FreeQ[lInt, y1Conj])
 |>;
 If[!AllTrue[Values[checks], TrueQ],
   SPFail["A required native identity failed: " <> SPString[checks] <>
      "; C^2=" <> SPString[diracC2] <>
      "; C^T=" <> SPString[transposeC] <>
      "; spin=" <> SPString[spin] <>
      "; expected=" <> SPString[byDefinition] <>
      "; PL flip=" <> SPString[flipL] <>
      "; PR flip=" <> SPString[flipR]]];
 step = "complete";
 <|
  "status" -> "NativeMatcheteAbstractChargeConjugationAndHcPassed",
  "format" -> "t3-workstream-a-a0-native-spinor-algebra-v1",
  "model" -> "A_p0",
  "checks" -> checks,
  "dirac_cc_square_inputform" -> SPString[diracC2],
  "gamma_cc_transpose_inputform" -> SPString[transposeC],
  "matchete_charge_conjugated_bilinear_inputform" -> SPString[spin],
  "matchete_explicit_cc_bilinear_inputform" -> SPString[byDefinition],
  "hermitian_conjugate_inputform" -> SPString[hc],
  "native_concrete_four_by_four_gamma_cc_sign_measured" -> False,
  "rgbeta_epsilon_component_orientation_measured" -> False,
  "physical_cross_package_coupling_map_certified" -> False,
  "full_A0_matchete_build_passed" -> True,
  "matching_or_betas_performed" -> False,
  "production_files_modified" -> False
 |>
];

step = "run";
res = CheckAbort[Catch[SPRun[], "SpinorPhaseStop"], $Aborted];
If[!AssociationQ[res], res = <|
 "status" -> "NativeSpinorConventionBlocked", "failed_step" -> step,
 "reason" -> "Wolfram evaluation aborted or did not return an association"|>];
If[!DirectoryQ[DirectoryName[out]],
 CreateDirectory[DirectoryName[out], CreateIntermediateDirectories -> True]];
export = Quiet@Check[Export[out, res, "RawJSON"], $Failed];
If[export === $Failed, Print["Unable to write native JSON: ", out]; Exit[7]];
Print["A_p0 spinor phase: ", Lookup[res, "status", "Missing"]];
If[KeyExistsQ[res, "reason"], Print["Failed: ", res["reason"]]];
If[Lookup[res, "status", ""] =!=
     "NativeMatcheteAbstractChargeConjugationAndHcPassed", Exit[6]];
Exit[0];
