(* Experimental general-quartic RGBeta runner; does not modify production model.
   Usage: wolframscript -file RunGeneralT3RGBeta.wl UV|EFT1 dS1 dS2 dF alpha MODEL output.json
   The additional quartic declarations are generated from the certified
   independent-operator manifest by scripts/GenerateGeneralT3RGBetaBasis.py.
   This provides formal RGBeta betas, not a certified Matchete coupling map. *)

ClearAll["Global`*"];
Needs["RGBeta`"];

scriptDir = DirectoryName @ ExpandFileName[$InputFileName];
rootDir = DirectoryName @ DirectoryName @ DirectoryName[scriptDir];
Get[FileNameJoin[{scriptDir, "T3RGBetaRunnerCommon.wl"}]];
Get[FileNameJoin[{scriptDir, "T3RGBetaModel.wl"}]];

If[Length[$ScriptCommandLine] =!= 8,
  Print["Usage: wolframscript -file RunGeneralT3RGBeta.wl UV|EFT1 dS1 dS2 dF alpha MODEL output.json"];
  Exit[2]
];

stage = ToUpperCase[$ScriptCommandLine[[2]]];
If[!MemberQ[{"UV", "EFT1"}, stage], Print["Unknown stage: ", stage]; Exit[2]];
{dS1,dS2,dF,alpha} = T3RGBetaParseIntegerToken /@ $ScriptCommandLine[[3;;6]];
modelName = $ScriptCommandLine[[7]];
outputPath = ExpandFileName[$ScriptCommandLine[[8]]];
(* ScriptCommandLine contains the 8 tokens: script,stage,4dims,model,out. *)

expectedModel = Which[
 {dS1,dS2,dF} === {1,3,2}, "A",
 {dS1,dS2,dF} === {2,2,1}, "B",
 {dS1,dS2,dF} === {2,2,3}, "C",
 {dS1,dS2,dF} === {3,1,2}, "D",
 {dS1,dS2,dF} === {3,3,2}, "E",
 True, "Unsupported"
];
expectedLabel = expectedModel <> "_" <> If[alpha < 0, "-", "+"] <> ToString[Abs[alpha]];
If[modelName =!= expectedLabel, Print["Model label/representation mismatch: ", {modelName, expectedLabel}]; Exit[3]];

extensionPath = FileNameJoin[{rootDir,"RGE","general","generated",modelName<>"_extra_quartics.wl"}];
If[!FileExistsQ[extensionPath], Print["Missing extension: ",extensionPath];Exit[4]];

build = CheckAbort[Quiet @ Check[
  If[stage === "UV", T3RGBetaBuild[dS1,dS2,dF,alpha],
     T3RGBetaBuildEFT1[dS1,dS2,dF,alpha]],
  $Failed
], $Aborted];
If[!AssociationQ[build], Print["General RGBeta base model build failed"];Exit[5]];

T3GeneralQuarticAddedSymbols = {};
registered = CheckAbort[Check[Get[extensionPath];True,$Failed],$Aborted];
If[registered =!= True ||
   Length[T3GeneralQuarticAddedSymbols] + T3GeneralQuarticReusedCount =!= T3GeneralQuarticExpectedCount,
  Print["General scalar quartic registration failed"]; Exit[6]
];

baseBetas = CheckAbort[Quiet @ Check[
 If[stage === "UV", T3RGBetaOneLoopBetas[], T3RGBetaEFT1OneLoopBetas[]],
 $Failed],$Aborted];
If[!AssociationQ[baseBetas],Print["Base RGBeta betas failed"];Exit[7]];

(* Important: a new quartic may have EXACTLY the same field list as an
   existing quartic but a different SU(2) contraction. RGBeta warns that
   BetaTerm[anIndividualQuartic, 1] is not the beta of its independent
   coupling in that case (CheckPotentialCouplingMixing::mixes).

   BetaTerm[Quartic, 1] invokes RGBeta's global quartic projection-matrix
   inversion. It must be used for *all* quartics, not merely the new ones:
   an old quartic with the same field content is mixed as well. *)
Print["Resolving full quartic coupling mixing with BetaTerm[Quartic, 1] ..."];
projectedQuartics = CheckAbort[BetaTerm[Quartic, 1], $Aborted];
If[!AssociationQ[projectedQuartics],
  Print["Full quartic projection failed: expected an Association from RGBeta."];
  Exit[8]
];

(* Convert RGBeta's symbol-keyed association into the existing string-keyed
   JSON contract without assuming a one-to-one field-content-to-beta map. *)
quarticNames = Association @ KeyValueMap[(SymbolName[#1] -> #2) &, projectedQuartics];
newNames = SymbolName /@ T3GeneralQuarticAddedSymbols;
missingNames = Complement[newNames, Keys[quarticNames]];
If[missingNames =!= {},
  Print["Projected quartic betas missing new couplings: ", missingNames];
  Exit[8]
];
If[AnyTrue[Values[quarticNames],
    (# === $Failed || # === $Aborted || !FreeQ[#, BetaTerm]) &],
  Print["Full quartic projection returned an unevaluated or failed beta."];
  Exit[8]
];

(* Replace *all* old quartic beta entries using the same projected system.
   Base gauge/Yukawa/mass entries are unchanged. Join overrides matching
   entries with the right-hand association, avoiding mixed conventions. *)
baseQuarticKeys = Intersection[Keys[baseBetas], Keys[quarticNames]];
baseBetas = Join[baseBetas, KeyTake[quarticNames, baseQuarticKeys]];
extraBetas = KeyTake[quarticNames, newNames];
If[Length[extraBetas] =!= Length[newNames],
  Print["Projected quartic beta count mismatch."]; Exit[8]
];
If[Intersection[Keys[baseBetas], Keys[extraBetas]] =!= {},
  Print["Duplicate RGBeta beta symbol"];Exit[9]
];
betas = Join[baseBetas,extraBetas];

metadataConfig = <|"SharedMode"->False,"dS1"->dS1,"dS2"->dS2,
 "dF"->dF,"alpha"->alpha|>;
metadata = Join[
 T3RGBetaCommonMetadata[metadataConfig,build],
 <|"GeneralScalarModel"->modelName,"GeneralScalarStage"->stage,
 "GeneralScalarBasisStatus"->"CandidateNativeRGBetaUnmappedToMatchete",
 "GeneralScalarAdditionalQuartics"->T3GeneralQuarticExpectedCount,
 "GeneralScalarExistingStructuralCount"->T3GeneralQuarticReusedCount,
 "GeneralScalarNewBetaKeys"->Keys[extraBetas],
 "GeneralScalarQuarticBetaBasis"->"RGBetaFullQuarticProjection",
 "GeneralScalarExistingQuarticBetasProjected"->baseQuarticKeys,
 "GeneralScalarQuadraticMassMixing"->"NotIncluded",
 "NativeRGBetaCheckInvarianceUsed"->False,
 "IndependentGaugeStructureChecks"->"U1Y/Z2/field typing/SU2 candidate rank; no RGBeta Ward check",
 "ReportBetaConvention"->"16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g"
 |>
];
If[stage === "UV",
 metadata = Join[metadata,
  <|"YF"->ToString[InputForm[build["YF"]]],
    "MajoranaFermion"->TrueQ[build["MajoranaFermion"]]|>]
];

report = T3RGBetaConventionalReportBetas[betas];
result = <|"status"->"Success","metadata"->metadata,
 "betas"->T3RGBetaSerializeAssociation[betas],
 "report_betas"->T3RGBetaSerializeAssociation[report]|>;
T3RGBetaWriteJSON[outputPath,result];
Print["Experimental general T3 RGBeta export: SUCCESS, model ",modelName,
 " stage ",stage," new quartics ",Length[extraBetas]];
Print["Output: ",outputPath];
