(* Formal source-term audit for the two C5-active B/C(-1) RGBeta quartics.

  wolframscript -file scripts/AnalyzeGeneralT3RGBetaSources.wl \
    output output/general_rgbeta_source_proof.json

  Reads existing report_betas; does not rerun RGBeta. Sets all six existing
  special quartic coefficients to zero, but leaves lambdaT3 and Yukawas
  independent. If a residual beta source is nonzero, the zero-special-quartic
  subspace is not RG invariant at one loop. This does not map the Matchete
  quartic normalisations onto RGBeta's tensors.
*)
ClearAll["Global`*"];
args = Rest[$ScriptCommandLine];
If[Length[args] =!= 2,
  Print["Usage: AnalyzeGeneralT3RGBetaSources.wl OUTPUT_ROOT REPORT_JSON"];
  Exit[64]
];
root = ExpandFileName[args[[1]]];
reportPath = ExpandFileName[args[[2]]];
specialNames = {
  "lambdaHHdagS2S2", "lambdaHHdagS1barS1bar",
  "lambdaS1bar2S2bar2", "lambdaS1barS2S2bar2",
  "lambdaS1S1bar2S2bar", "lambdaHHdagS1barS2barCross"
};
interesting = {"lambdaHHdagS2S2", "lambdaHHdagS1barS1bar"};

T3RGRead[path_String] := Quiet@Check[Import[path, "RawJSON"], $Failed];
T3RGParse[text_String] := Module[{held},
  held = Quiet@Check[ToExpression[text, InputForm, HoldComplete], $Failed];
  If[held === $Failed, Return[$Failed]];
  Quiet@Check[ReleaseHold[held], $Failed]
];
T3RGIsSpecial[x_Symbol] := MemberQ[specialNames, ToString[x, InputForm]];
T3RGIsSpecial[_] := False;
T3RGSource[expr_] := Module[{result},
  result = expr /. {
    HoldPattern[Bar[s_Symbol]] /; T3RGIsSpecial[s] :> 0,
    HoldPattern[Conjugate[s_Symbol]] /; T3RGIsSpecial[s] :> 0,
    s_Symbol /; T3RGIsSpecial[s] :> 0
  };
  Quiet@Check[Expand[result], $Failed]
];

T3RGOne[model_String] := Module[{file, obj, beta, entry, expr, source, results = {}},
  file = FileNameJoin[{root, "T3_" <> model <> "_alpha_m1", "data",
    "uv_rgbeta_rge.json"}];
  If[!FileExistsQ[file], Return[<|"model" -> model <> "_m1",
    "status" -> "MissingRGBetaReport", "path" -> file|>]];
  obj = T3RGRead[file];
  If[!AssociationQ[obj], Return[<|"model" -> model <> "_m1",
    "status" -> "MalformedRGBetaReport", "path" -> file|>]];
  beta = Lookup[obj, "report_betas", <||>];
  If[!AssociationQ[beta], Return[<|"model" -> model <> "_m1",
    "status" -> "MissingReportBetas", "path" -> file|>]];
  Do[
    entry = Lookup[beta, coupling, Missing["NoBeta"]];
    If[!StringQ[entry],
      AppendTo[results, <|"coupling" -> coupling, "status" -> "MissingBeta"|>];
      Continue[]
    ];
    expr = T3RGParse[entry];
    If[expr === $Failed,
      AppendTo[results, <|"coupling" -> coupling, "status" -> "ParseFailed"|>];
      Continue[]
    ];
    source = T3RGSource[expr];
    AppendTo[results, <|
      "coupling" -> coupling,
      "status" -> If[source === $Failed, "ReductionFailed", "Analyzed"],
      "nonzero_formal_source" -> If[source === $Failed, Null,
        !TrueQ[PossibleZeroQ[source]]],
      "source_at_zero_special_quartics_inputform" -> If[source === $Failed,
        "", ToString[source, InputForm]],
      "source_independent_of_special_quartics" -> If[source === $Failed,
        False, StringFreeQ[ToString[source, InputForm],
          Alternatives @@ specialNames]]
    |>],
    {coupling, interesting}
  ];
  <|"model" -> model <> "_m1", "path" -> file,
    "status" -> If[AllTrue[results, # ["status"] === "Analyzed" &],
      "Analyzed", "ReviewRequired"],
    "couplings" -> results|>
];

results = T3RGOne /@ {"B", "C"};
valid = AllTrue[results, # ["status"] === "Analyzed" &] &&
  AllTrue[Flatten[Lookup[results, "couplings", {}], 1],
    TrueQ[Lookup[#, "nonzero_formal_source", False]] &&
    TrueQ[Lookup[#, "source_independent_of_special_quartics", False]] &];
report = <|
  "status" -> If[valid, "OneLoopSourcesVerified", "ReviewRequired"],
  "one_loop_beta_convention" -> "16*pi^2*dX/dln(mu) = beta^(1)",
  "set_zero" -> specialNames,
  "remaining_general" -> {"lambdaT3", "y1", "y2", "SM and ordinary quartics"},
  "matchete_rgbeta_coupling_map_verified" -> False,
  "note" -> "A nonzero source under this substitution means boundary zeros are not RG-closed. A generated coupling inserted into one-loop matching enters at formally higher perturbative order.",
  "results" -> results
|>;
If[!DirectoryQ[DirectoryName[reportPath]],
  CreateDirectory[DirectoryName[reportPath], CreateIntermediateDirectories -> True]
];
Export[reportPath, report, "RawJSON"];
Scan[Function[item, Print[item["model"], " ", item["status"]];
  Scan[Function[c, Print["  ", c["coupling"], " nonzero source=",
    Lookup[c, "nonzero_formal_source", Null]]],
    Lookup[item, "couplings", {}]]], results];
Print["RGE source proof: ", report["status"]];
Print["Report: ", reportPath];
Exit[If[valid, 0, 2]];
