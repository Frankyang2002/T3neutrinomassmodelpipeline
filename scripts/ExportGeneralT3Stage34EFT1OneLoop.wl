(* Stage 34: independent full projected scalar-only EFT1 beta export.
   Uses the existing RGBeta EFT1 builder without modifying production files.
   Usage: wolframscript -file scripts/ExportGeneralT3Stage34EFT1OneLoop.wl dS1 dS2 dF alpha out.json
   e.g. ... 1 3 2 m2 output/general_t3_stage22/stage34_A_m2_eft1_full_betas.json
   This exporter intentionally does not register Stage 29 new quartics, for
   which no physical intermediate-EFT mapping has been certified.
*)
ClearAll["Global`*"];
Needs["RGBeta`"];

scriptDir = DirectoryName[$InputFileName];
Get[FileNameJoin[{scriptDir, "../RGE/running/rgbeta/T3RGBetaRunnerCommon.wl"}]];
Get[FileNameJoin[{scriptDir, "../RGE/running/rgbeta/T3RGBetaModel.wl"}]];

T3S34Fail[reason_String] := (
  Print["Stage34 EXPORT ERROR: ", reason];
  Throw[$Failed, "T3S34"]
);
T3S34ErrorExpressionQ[expr_] := (
  expr === $Failed || expr === $Aborted ||
  !FreeQ[Unevaluated[expr], _RGBeta`BetaTerm] ||
  !FreeQ[Unevaluated[expr], _Missing]
);

T3S34Export[] := Module[
  {cfg, build, projected, bare, report, meta, path, nonQuartic,
   quartics, all, newKeys, unexpected, messages = {}, key, beta, result},
  cfg = T3RGBetaParseRunnerArguments[];
  If[TrueQ[cfg["SharedMode"]], T3S34Fail["Only ordinary split-scalar T3 is supported"]];
  If[!MemberQ[{{1,3,2},{2,2,1},{2,2,3},{3,1,2},{3,3,2}},
     {cfg["dS1"],cfg["dS2"],cfg["dF"]}],
     T3S34Fail["Unsupported scalar/fermion representation"]];
  path = cfg["OutputPath"];
  Print["Stage34: building ordinary EFT1 ",
    {cfg["dS1"],cfg["dS2"],cfg["dF"],cfg["alpha"]}];
  build = T3RGBetaBuildEFT1[cfg["dS1"],cfg["dS2"],cfg["dF"],cfg["alpha"]];
  If[!AssociationQ[build], T3S34Fail["EFT1 model builder did not return an Association"]];

  (* Gauge and SM Yukawa, then both scalar mass parameters. *)
  nonQuartic = <|
    "gY" -> RGBeta`BetaTerm[gY,1],
    "g2" -> RGBeta`BetaTerm[g2,1],
    "g3" -> RGBeta`BetaTerm[g3,1],
    "yu" -> RGBeta`BetaTerm[yu,1],
    "yd" -> RGBeta`BetaTerm[yd,1],
    "ye" -> RGBeta`BetaTerm[ye,1],
    "mS1Sq" -> RGBeta`BetaTerm[mS1Sq,1],
    "mS2Sq" -> RGBeta`BetaTerm[mS2Sq,1]
  |>;

  Print["Stage34: projecting native full EFT1 quartic class"];
  (* An RGBeta Monitor may print FrontEndObject::notavail under wolframscript,
     which is not a beta-function failure. Avoid Check[..., $Failed] here. *)
  projected = CheckAbort[
    Block[{$MessageList = {}},
      result = RGBeta`BetaTerm[RGBeta`Quartic,1];
      messages = DeleteDuplicates[ToString[#,InputForm]& /@ $MessageList];
      result
    ],
    $Aborted
  ];
  Print["Stage34: projected head=", Head[projected], "; native messages=", Length[messages]];
  If[!AssociationQ[projected] || Length[projected] == 0,
     T3S34Fail["Native quartic-class result is not a nonempty Association"]];
  If[AnyTrue[messages, StringContainsQ[#, "singular"] &],
     T3S34Fail["Singular quartic projection"]];
  quartics = Association @ KeyValueMap[(SymbolName[#1] -> #2) &, projected];
  If[!FreeQ[quartics, _Missing], T3S34Fail["Projected beta contains Missing"]];
  newKeys = Keys[quartics];
  If[Length[newKeys] != Length[DeleteDuplicates[newKeys]] ||
     Length[Intersection[newKeys, Keys[nonQuartic]]] > 0,
     T3S34Fail["Duplicate/overlapping coupling inventory"]];
  all = Join[nonQuartic, quartics];
  If[AnyTrue[Values[all], T3S34ErrorExpressionQ],
     T3S34Fail["One or more EFT1 beta expressions are missing/unresolved"]];
  If[!SubsetQ[Keys[all], {"gY","g2","g3","yu","yd","ye",
       "mS1Sq","mS2Sq","lambdaH","lambdaS1","lambdaS2",
       "lambdaH1","lambdaH2","lambda12","lambdaT3"}],
     T3S34Fail["Required EFT1 baseline beta inventory is incomplete"]];
  meta = Join[T3RGBetaCommonMetadata[cfg, build],
     <|"IntegratedField" -> "F", "ActiveBSMFields" -> {"S1","S2"},
       "ReportBetaConvention" -> "16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g"|>];
  report = T3RGBetaConventionalReportBetas[all];
  If[AnyTrue[Values[report], T3S34ErrorExpressionQ],
     T3S34Fail["Reported beta expression unresolved"]];
  T3RGBetaWriteJSON[path, <|
    "format" -> "stage34-ordinary-t3-eft1-projected-rgbeta-v1",
    "status" -> "Success",
    "metadata" -> meta,
    "source" -> "RGBeta BetaTerm[Quartic,1] projected in scalar-only EFT1",
    "physical_coupling_mapping_certified" -> False,
    "new_stage29_quartics_integrated" -> False,
    "beta_symbols" -> Keys[all],
    "report_betas" -> T3RGBetaSerializeAssociation[report]
  |>];
  Print["Stage34 EFT1 projected beta export PASSED: ",Length[all]," couplings"];
  Print["Report: ", path];
  True
];

result = CheckAbort[Catch[T3S34Export[], "T3S34"], $Aborted];
If[result =!= True, Print["Stage34 EFT1 beta export FAILED"]; Exit[2]];
Exit[0];
