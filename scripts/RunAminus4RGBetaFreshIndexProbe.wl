(* A(-4) EFT1 local-only RGBeta adjoint dummy-index experiment.

   Run from project root:
   wolframscript -file scripts/RunAminus4RGBetaFreshIndexProbe.wl MODE OUTPUT.json

   MODE = baseline | triplet | portal | both.
   This script loads the UNMODIFIED installed RGBeta and repository model,
   then optionally overrides exactly the A(-4) RG-closed quartic registration
   in memory, with a fresh Module-local adjoint dummy for each invocation.
   It NEVER modifies installed RGBeta, the model source or saved betas.
   Each MODE MUST execute in its OWN fresh WolframScript process.
   Numerical cross-comparison is performed by the Python batch driver.
*)
ClearAll["Global`*"];
args = Rest[$ScriptCommandLine];
If[Length[args] =!= 2,
  Print["Usage: MODE OUTPUT.json; MODE=baseline|triplet|portal|both"];
  Exit[2]];
mode = First[args];
If[!MemberQ[{"baseline", "triplet", "portal", "both"}, mode],
  Print["Unknown mode: ", mode]; Exit[2]];
out = ExpandFileName[Last[args]];
If[FileExistsQ[out] || FileExistsQ[out <> ".partial"],
  Print["Refusing to overwrite a native report: ", out]; Exit[3]];
root = DirectoryName @ DirectoryName @ ExpandFileName[$InputFileName];
modelSource = FileNameJoin[{root, "RGE", "running", "rgbeta", "T3RGBetaModel.wl"}];
commonSource = FileNameJoin[{root, "RGE", "running", "rgbeta", "T3RGBetaRunnerCommon.wl"}];
If[!FileExistsQ[modelSource] || !FileExistsQ[commonSource],
  Print["Missing unchanged RGBeta model or reporting source"];Exit[3]];
modelSHA = IntegerString[FileHash[modelSource, "SHA256"], 16, 64];
(* Pin to the model inspected at GitHub main 54857efafc7261f0e4c2bd0802f23e61dc3c165a.
   Review changed source before attempting this deliberately narrow override. *)
If[modelSHA =!= "a58c79f18b2abdf9d476588552dfa52e2f131f322094e5b72d1198162aa20abc",
  Print["Source changed: T3RGBetaModel.wl SHA256 = ", modelSHA]; Exit[3]];

flags = {"T3_EXPERIMENTAL_GENERAL_FERMIONS",
  "T3_EXPERIMENTAL_GENERAL_QUARTICS",
  "T3_EXPERIMENTAL_GENERAL_QUADRATICS",
  "T3_EXPERIMENTAL_FLAVOR_MASSES"};
SetEnvironment[Thread[flags -> ConstantArray["0", Length[flags]]]];
If[!AllTrue[flags, Environment[#] === "0" &],
  Print["Could not enforce restricted T3 flags"];Exit[3]];
If[CheckAbort[Check[Needs["RGBeta`"] ; True, $Failed], $Aborted] =!= True,
  Print["Could not load installed RGBeta"];Exit[4]];
Get[commonSource];
Get[modelSource];

freshTriplet = MemberQ[{"triplet", "both"}, mode];
freshPortal = MemberQ[{"portal", "both"}, mode];

(* This A(-4)-ONLY override is the same two quartics that the generic
   T3RGBetaAddRGClosedQuartics[1,3,2,-4] source registers.  No other
   extra quartic is allowed for this model. Use the original registration
   for each channel not selected for the fresh-index intervention. *)
If[mode =!= "baseline",
  Clear[T3RGBetaAddRGClosedQuartics];
  T3RGBetaAddRGClosedQuartics[1,3,2,-4] := Module[{portalInvariant, tripletInvariant},
    portalInvariant = If[freshPortal,
      (Module[{adj},
        tGen[SU2L @ fund, adj, #1, #2] *
        tGen[SU2L @ S2, adj, #3, #4]] &),
      T3RGBetaPortalAdjInvariant[3]
    ];
    tripletInvariant = If[freshTriplet,
      (Module[{adj},
        tGen[SU2L @ S2, adj, #1, #2] *
        tGen[SU2L @ S2, adj, #3, #4] / 2] &),
      (tGen[SU2L @ S2, A, #1, #2] *
       tGen[SU2L @ S2, A, #3, #4] / 2 &)
    ];
    SetReal[lambdaH2Adj];
    AddQuartic[lambdaH2Adj, {Bar @ H, H, Bar @ S2, S2},
      GroupInvariant -> portalInvariant];
    SetReal[lambdaS2Adj];
    AddQuartic[lambdaS2Adj, {Bar @ S2, S2, Bar @ S2, S2},
      GroupInvariant -> tripletInvariant];
  ];
];

built = CheckAbort[Quiet @ Check[
  T3RGBetaBuildEFT1[1, 3, 2, -4], $Failed], $Aborted];
If[!AssociationQ[built] || built["IntegratedField"] =!= "F" ||
   built["ActiveBSMFields"] =!= {"S1", "S2"},
  Print["A(-4) F-integrated EFT1 build failed in mode ", mode]; Exit[4]];

expected = Sort[{"lambdaH", "lambdaS1", "lambdaS2", "lambdaS2Adj",
  "lambdaH1", "lambdaH2", "lambdaH2Adj", "lambda12", "lambdaT3"}];
registered = If[AssociationQ[$quartics], Sort[SymbolName /@ Keys[$quartics]], {}];
If[registered =!= expected,
  Print["Quartic inventory changed: ", registered]; Exit[5]];

(* Match the already working native RunGeneralT3RGBeta.wl evaluation order:
   build, load the A_-4 generated (zero-extra) quartic extension, evaluate
   the ordinary EFT1 beta inventory, THEN perform the global quartic
   projection. The zero-extra extension must not change the registered
   operator inventory; it is still loaded to match the production path.

   The original probe went directly from the build to BetaTerm[Quartic,1]
   and wrapped that call in Check[..., $Failed]. RGBeta can issue messages
   during its global quartic projection; Check would convert *any* such
   message to $Failed even if a valid Association was returned. Preserve
   diagnostic messages and inspect the actual result instead.
*)
extensionSource = FileNameJoin[{root, "RGE", "general", "generated",
  "A_-4_extra_quartics.wl"}];
If[!FileExistsQ[extensionSource],
  Print["Missing unchanged A_-4 extension: ", extensionSource]; Exit[5]];
T3GeneralQuarticAddedSymbols = {};
extensionLoaded = CheckAbort[Check[Get[extensionSource]; True, $Failed], $Aborted];
If[extensionLoaded =!= True || T3GeneralQuarticExpectedCount =!= 0 ||
   T3GeneralQuarticReusedCount =!= 0 ||
   Length[T3GeneralQuarticAddedSymbols] =!= 0,
  Print["A_-4 native extension not the expected zero-extra registry"]; Exit[5]];
registeredAfterExtension = If[AssociationQ[$quartics],
  Sort[SymbolName /@ Keys[$quartics]], {}];
If[registeredAfterExtension =!= expected,
  Print["The zero-extra extension changed quartic inventory"];
  Exit[5]];

baseBetas = CheckAbort[Quiet @ Check[
  T3RGBetaEFT1OneLoopBetas[], $Failed], $Aborted];
If[!AssociationQ[baseBetas] ||
   Sort[Intersection[Keys[baseBetas], expected]] =!= expected ||
   AnyTrue[Values[baseBetas], (# === $Failed || # === $Aborted) &],
  Print["Ordinary EFT1 native beta prepass failed: ", mode,
    "; head=", Head[baseBetas]];
  Exit[5]];
Print["Native EFT1 beta prepass completed: ", mode];

(* Capture two raw invocations to distinguish fresh Module-generated symbols
   from the original hard-coded adjoint index.  This is a raw registry
   observation only: it does not assert how the global projector reduces it. *)
ClearAll[indexWitness];
indexWitness[name_String] := Module[{entry, keys, f, x, y},
  entry = If[AssociationQ[$quartics],
    Lookup[$quartics, ToExpression[name], <||>], <||>];
  keys = If[AssociationQ[entry], Select[Keys[entry],
    Last[StringSplit[ToString[InputForm[#]], "`"]] === "Invariant" &], {}];
  If[Length[keys] =!= 1, Return[<|"name" -> name,
     "available" -> False|>]];
  f = Lookup[entry, First[keys]];
  x = ToString[InputForm[f[idx1, idx2, idx3, idx4]]];
  y = ToString[InputForm[f[jdx1, jdx2, jdx3, jdx4]]];
  <|"name" -> name, "available" -> True,
    "raw_first_inputform" -> StringTake[x, UpTo[2000]],
    "raw_second_inputform" -> StringTake[y, UpTo[2000]],
    "first_has_Module_generated_adj_index" -> StringContainsQ[x, "adj$"],
    "second_has_Module_generated_adj_index" -> StringContainsQ[y, "adj$"],
    "two_raw_internal_dummy_names_different" ->
      (StringContainsQ[x, "adj$"] && StringContainsQ[y, "adj$"] &&
       StringCases[x, "adj$" ~~ DigitCharacter ..] =!=
       StringCases[y, "adj$" ~~ DigitCharacter ..])|>
];
indexChecks = indexWitness /@ {"lambdaS2Adj", "lambdaH2Adj"};

Print["Global A(-4) quartic projection mode: ", mode];
(* Deliberately NO Check around the global projection. The production
   RunGeneralT3RGBeta.wl uses CheckAbort[BetaTerm[Quartic,1],$Aborted]
   and accepts the result only after checking AssociationQ and keys. *)
projected = CheckAbort[BetaTerm[Quartic, 1], $Aborted];
If[!AssociationQ[projected],
  Print["Global quartic projector did not return an Association: ", mode,
    "; head=", Head[projected],
    "; result prefix=", StringTake[ToString[InputForm[projected]], UpTo[500]]];
  Exit[5]];
named = Association @ KeyValueMap[(SymbolName[#1] -> #2) &, projected];
If[Sort[Keys[named]] =!= expected ||
   AnyTrue[Values[named], (# === $Failed || # === $Aborted ||
     !FreeQ[#, BetaTerm]) &],
  Print["Incomplete or unresolved projected quartic beta coefficients"]; Exit[5]];
reportBeta = T3RGBetaConventionalReportBetas[named];
result = <|
  "format" -> "t3-a-minus4-fresh-adjoint-index-native-probe-v1",
  "status" -> "NATIVE_QUARTIC_BETAS_EXPORTED_NOT_PHYSICAL_CERTIFICATION",
  "model" -> "A_-4", "stage" -> "EFT1", "mode" -> mode,
  "source_model_SHA256" -> modelSHA,
  "installed_rgbeta_package" -> ToString[InputForm[FindFile["RGBeta`"]]],
  "fresh_triplet_adjoint_dummy" -> TrueQ[freshTriplet],
  "fresh_portal_adjoint_dummy" -> TrueQ[freshPortal],
  "all_9_native_quartic_symbols" -> expected,
  "baseline_prepass_matches_native_global_runner" -> True,
  "global_projector_uses_CheckAbort_not_message_sensitive_Check" -> True,
  "registered_invariant_raw_two_call_witness" -> indexChecks,
  "report_betas" -> T3RGBetaSerializeAssociation[reportBeta],
  "production_file_modified" -> False,
  "matchete_to_rgbeta_coupling_map_certified" -> False,
  "physical_beta_certified" -> False
|>;
serialized = CheckAbort[Check[ExportString[result, "RawJSON"], $Failed], $Aborted];
If[!StringQ[serialized] || StringLength[serialized] < 200,
  Print["Could not serialize native betas"];Exit[5]];
If[!DirectoryQ[DirectoryName[out]],
  CreateDirectory[DirectoryName[out], CreateIntermediateDirectories -> True]];
partial = out <> ".partial";
If[CheckAbort[Check[Export[partial, serialized, "Text"], $Failed], $Aborted] === $Failed ||
   !FileExistsQ[partial] || FileByteCount[partial] < 200,
  Print["Cannot write native partial output"];Exit[6]];
roundtrip = CheckAbort[Check[Import[partial, "RawJSON"], $Failed], $Aborted];
If[!AssociationQ[roundtrip] || roundtrip["mode"] =!= mode ||
   roundtrip["source_model_SHA256"] =!= modelSHA,
  Print["Native JSON roundtrip failed"];Exit[6]];
If[CheckAbort[Check[RenameFile[partial, out], $Failed], $Aborted] === $Failed,
  Print["Could not finalize native result"];Exit[6]];
Print["Nine quartic betas exported; no production source edits: ", out];
Exit[0];
