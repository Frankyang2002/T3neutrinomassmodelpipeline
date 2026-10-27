(* Experimental UV/EFT1 one-loop RGBeta mass-and-quartic closure.
   Usage:
    wolframscript -file RGE/running/rgbeta/RunGeneralT3QuadraticRGBeta.wl \
        EFT1 2 2 1 m1 B_-1 output.json

   Deliberately DOES NOT alter T3RGBetaModel.wl or production state.
   Additional mass couplings use RGBeta's scalar POTENTIAL convention,
   V2=mu2*(Bar[S])^2+h.c. (or mixed Bar[S1] Bar[S2]).  Physical
   Matchete<->RGBeta normalization, triplet bars and mass thresholds
   are NOT certified.  The full quartic projector is essential.
*)
ClearAll["Global`*"];
Needs["RGBeta`"];
scriptDir = DirectoryName @ ExpandFileName[$InputFileName];
rootDir = DirectoryName @ DirectoryName @ DirectoryName[scriptDir];
Get[FileNameJoin[{scriptDir, "T3RGBetaRunnerCommon.wl"}]];
Get[FileNameJoin[{scriptDir, "T3RGBetaModel.wl"}]];

T3QFail[message_] := (Print["Quadratic RGBeta error: ", message]; Exit[2]);
args = Rest[$ScriptCommandLine];
If[Length[args] =!= 7,
 T3QFail["Usage: UV|EFT1 dS1 dS2 dF alpha MODEL output.json"]];
stage = ToUpperCase[args[[1]]];
If[!MemberQ[{"UV", "EFT1"}, stage], T3QFail["Invalid stage"]];
{dS1,dS2,dF,alpha} = T3RGBetaParseIntegerToken /@ args[[2;;5]];
modelName = args[[6]];
outputPath = ExpandFileName[args[[7]]];
modelClass = Switch[{dS1,dS2,dF},
 {1,3,2},"A", {2,2,1},"B", {2,2,3},"C",
 {3,1,2},"D", {3,3,2},"E", _,"INVALID"];
expected = modelClass <> "_" <> If[alpha < 0,"-","+"] <> ToString[Abs[alpha]];
If[modelName =!= expected, T3QFail["Model/representation mismatch"]];

registryPath = FileNameJoin[{rootDir,"RGE","general","generated",
 "general_rgbeta_quartic_registry.json"}];
registry = Import[registryPath,"RawJSON"];
If[!AssociationQ[registry] ||
    !KeyExistsQ[registry["models"],modelName],
 T3QFail["Missing model from Stage29 registry"]];
entry = registry["models",modelName];
If[entry["dimensions"] =!= {dS1,dS2,dF} ||
    entry["alpha"] =!= alpha,
 T3QFail["Registry field content mismatch"]];

build = CheckAbort[
 If[stage === "UV",T3RGBetaBuild[dS1,dS2,dF,alpha],
                   T3RGBetaBuildEFT1[dS1,dS2,dF,alpha]],$Aborted];
If[!AssociationQ[build], T3QFail["Cannot construct RGBeta base theory"]];

extension = FileNameJoin[{rootDir,"RGE","general","generated",
  modelName<>"_extra_quartics.wl"}];
T3GeneralQuarticAddedSymbols = {};
register = CheckAbort[Check[Get[extension];True,$Failed],$Aborted];
If[register =!= True ||
  Length[T3GeneralQuarticAddedSymbols] + T3GeneralQuarticReusedCount =!=
    entry["total_additional_quartics"],
 T3QFail["General-quartic registration failed"]];

(* Exactly eight model-specific complex mass-parameter declarations in the
   ordinary 16-model catalogue; zero for all remaining assignments. *)
extraMassSymbol = None;
extraQuadraticKind = "none";
Which[
 alpha == 0 && MemberQ[{"A","D","E"},modelClass],
    extraMassSymbol = muT3GeneralS1Sq;
    extraQuadraticKind = "self_s1";
    AddScalarMass[muT3GeneralS1Sq,{Bar @ S1,Bar @ S1},
      GroupInvariant -> If[dS1 == 1,(1 &),
         (del[SU2L @ S2,#1,#2] &)], SelfConjugate -> False],
 alpha == -2 && MemberQ[{"A","D","E"},modelClass],
    extraMassSymbol = muT3GeneralS2Sq;
    extraQuadraticKind = "self_s2";
    AddScalarMass[muT3GeneralS2Sq,{Bar @ S2,Bar @ S2},
      GroupInvariant -> If[dS2 == 1,(1 &),
         (del[SU2L @ S2,#1,#2] &)], SelfConjugate -> False],
 alpha == -1 && MemberQ[{"B","C"},modelClass],
    extraMassSymbol = muT3GeneralS1S2;
    extraQuadraticKind = "mixed";
    AddScalarMass[muT3GeneralS1S2,{Bar @ S1,Bar @ S2},
      GroupInvariant -> (eps[SU2L @ fund,#1,#2] &),
      SelfConjugate -> False]
];
If[(extraMassSymbol === None) =!= (entry["quadratics_not_implemented"] === 0),
 T3QFail["One-to-one quadratic registry coverage failed"]];

Print["Calculating general ",stage," beta functions for ",modelName,
 " with quadratic kind ",extraQuadraticKind];
base = CheckAbort[Quiet @ Check[
 If[stage === "UV",T3RGBetaOneLoopBetas[],T3RGBetaEFT1OneLoopBetas[]],
 $Failed],$Aborted];
If[!AssociationQ[base],T3QFail["Base beta computation failed"]];

(* All quartics must be reprojected together; never use isolated
   BetaTerm[lambda,1] when two invariants have identical field content. *)
projected = CheckAbort[BetaTerm[Quartic,1],$Aborted];
If[!AssociationQ[projected], T3QFail["Full quartic projection failed"]];
quartics = Association @ KeyValueMap[(SymbolName[#1] -> #2)&, projected];
newNames = SymbolName /@ T3GeneralQuarticAddedSymbols;
If[!SubsetQ[Keys[quartics],newNames],
 T3QFail["Missing general quartic beta symbols"]];
base = Join[base,quartics];
(* Recompute both diagonal mass beta functions WITH the new bilinear vertex.
   The first base calculation already sees it; these explicit overrides are
   a diagnostic guard against accidentally retaining pre-registration betas. *)
base = Join[base,<|
 "mS1Sq"->BetaTerm[mS1Sq,1],
 "mS2Sq"->BetaTerm[mS2Sq,1]|>];
If[extraMassSymbol =!= None,
 massBeta = CheckAbort[BetaTerm[extraMassSymbol,1],$Aborted];
 base = Join[base,<|SymbolName[extraMassSymbol]->massBeta|>]
];
If[AnyTrue[Values[base],
 (# === $Failed || # === $Aborted || !FreeQ[#,BetaTerm|Missing|Indeterminate] &)],
 T3QFail["Unresolved mass/quartic beta expression"]];
meta = Join[
 T3RGBetaCommonMetadata[<|"SharedMode"->False,"dS1"->dS1,
 "dS2"->dS2,"dF"->dF,"alpha"->alpha|>,build],
 <|"GeneralScalarModel"->modelName,"GeneralScalarStage"->stage,
 "GeneralScalarBasisStatus"->"CandidateNativeRGBetaUnmappedToMatchete",
 "GeneralScalarQuarticBetaBasis"->"RGBetaFullQuarticProjection",
 "GeneralScalarNewBetaKeys"->newNames,
 "GeneralScalarAdditionalQuartics"->entry["total_additional_quartics"],
 "GeneralScalarExistingStructuralCount"->T3GeneralQuarticReusedCount,
 "GeneralScalarQuadraticMassMixing"->"IncludedCandidate",
 "QuadraticKind"->extraQuadraticKind,
 "QuadraticMassSymbol"->If[extraMassSymbol === None,Null,SymbolName[extraMassSymbol]],
 "QuadraticPrefactorConvention"->"V2 = mu2 * invariant + h.c.; no 1/2",
 "PhysicalMatcheteRGBetaMassBasisCertified"->False,
 "ReportBetaConvention"->"16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g"|>];
If[stage === "UV",
 meta = Join[meta,<|"YF"->ToString[InputForm[build["YF"]]],
 "MajoranaFermion"->TrueQ[build["MajoranaFermion"]]|>]
];
report = T3RGBetaConventionalReportBetas[base];
If[!DirectoryQ[DirectoryName[outputPath]],
 CreateDirectory[DirectoryName[outputPath],CreateIntermediateDirectories->True]];
T3RGBetaWriteJSON[outputPath,<|
 "format"->"general-t3-quadratic-rgbeta-experimental-v1",
 "status"->"SuccessCandidateMassBeta",
 "metadata"->meta,
 "beta_symbols"->Keys[base],
 "report_betas"->T3RGBetaSerializeAssociation[report],
 "physical_mapping_certified"->False,
 "mass_eigenstate_threshold_matching_performed"->False|>];
Print["Candidate RGBeta ",stage," quadratic mass beta export SUCCESS: ",modelName];
