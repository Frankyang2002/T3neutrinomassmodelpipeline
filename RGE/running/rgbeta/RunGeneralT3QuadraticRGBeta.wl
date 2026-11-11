(* Batch 13: joint UV/EFT1 renormalisable beta export for all ordinary T3 models.
   Usage: wolframscript -file RGE/running/rgbeta/RunGeneralT3QuadraticRGBeta.wl \
       UV 1 3 2 0 A_+0 output.json
   Restricted Yukawa sector: only the two topology Yukawa matrices y1,y2.
   The optional independent conjugate-scalar Yukawas y1Conj/y2Conj are
   excluded by model choice, not by gauge invariance or Majorana identity.
   Neutral heavy fermions are Majorana only for alpha=-1; charged/vectorlike
   heavy fermions retain the existing Dirac mass construction.
   EFT1 contains no heavy fermions or heavy-fermion Yukawas.
   Coupling map to Matchete, fermion flavour masses, and mixed-scalar
   threshold matching remain unverified; this is NOT a full UV certificate.
*)
ClearAll["Global`*"];
Needs["RGBeta`"];
scriptDir = DirectoryName @ ExpandFileName[$InputFileName];
rootDir = DirectoryName @ DirectoryName @ DirectoryName[scriptDir];
Get[FileNameJoin[{scriptDir, "T3RGBetaRunnerCommon.wl"}]];
Get[FileNameJoin[{scriptDir, "T3RGBetaModel.wl"}]];

T3QFail[message_] := (Print["General RGBeta error: ", message]; Exit[2]);
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
If[!AssociationQ[registry] || !AssociationQ[Lookup[registry,"models",None]] ||
   !KeyExistsQ[registry["models"],modelName],
 T3QFail["Missing model from general-quartic registry"]];
entry = registry["models",modelName];
If[entry["dimensions"] =!= {dS1,dS2,dF} || entry["alpha"] =!= alpha,
 T3QFail["Registry field content mismatch"]];

(* T3RGBetaBuild consults this process environment variable.  Explicitly
   disable the optional Yukawa extension even if the caller's shell has it
   set.  The restriction is imposed in BOTH stages, without altering the
   global model builder or the existing threshold prescription. *)
SetEnvironment["T3_EXPERIMENTAL_GENERAL_FERMIONS" -> "0"];
build = CheckAbort[Check[
 If[stage === "UV", T3RGBetaBuild[dS1,dS2,dF,alpha],
                    T3RGBetaBuildEFT1[dS1,dS2,dF,alpha]],
 $Failed],$Aborted];
If[!AssociationQ[build], T3QFail["Cannot construct RGBeta base theory"]];

(* The restricted UV Yukawa basis contains y1 and y2, and no independent
   conjugate-scalar Yukawas.  Assert the registry as well as the export. *)
restrictedYukawaNames = {"y1", "y2"};
excludedYukawaNames = {"y1Conj", "y2Conj"};
If[stage === "UV",
 registeredYukawaNames = SymbolName /@ Keys[$yukawas];
 If[!SubsetQ[registeredYukawaNames,restrictedYukawaNames],
  T3QFail["Restricted topology Yukawa missing from RGBeta registry"]];
 If[Intersection[registeredYukawaNames,excludedYukawaNames] =!= {},
  T3QFail["Conjugate-scalar Yukawa unexpectedly registered"]]
];

extension = FileNameJoin[{rootDir,"RGE","general","generated",
 modelName<>"_extra_quartics.wl"}];
If[!FileExistsQ[extension],T3QFail["Missing general quartic extension"]];
T3GeneralQuarticAddedSymbols = {};
register = CheckAbort[Check[Get[extension];True,$Failed],$Aborted];
If[register =!= True ||
 Length[T3GeneralQuarticAddedSymbols] + T3GeneralQuarticReusedCount =!=
 entry["total_additional_quartics"],
 T3QFail["General-quartic registration failed"]];

(* Additional complex quadratic terms in the 16-model catalogue. These
   coefficients use RGBeta's potential convention V2=mu*invariant+h.c.
   The d=3 contraction uses the S2 triplet irrep, independently of whether
   the field being contracted is the first or second scalar multiplet. *)
extraMassSymbol = None;
extraQuadraticKind = "none";
quadraticRegistration = CheckAbort[Check[
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
 ]; True,
 $Failed],$Aborted];
If[quadraticRegistration =!= True,
 T3QFail["Extra quadratic coupling registration failed"]];
If[(extraMassSymbol === None) =!= (entry["quadratics_not_implemented"] === 0),
 T3QFail["Quadratic registry coverage mismatch"]];

Print["Calculating general ",stage," beta functions for ",modelName,
 " with quadratic kind ",extraQuadraticKind];
base = CheckAbort[Check[
 Block[{Monitor = (#1 &)},
  If[stage === "UV",T3RGBetaOneLoopBetas[],T3RGBetaEFT1OneLoopBetas[]]],
 $Failed],$Aborted];
If[!AssociationQ[base],T3QFail["Base beta computation failed"]];
Print["Base RGBeta beta keys: ", Keys[base]];
If[stage === "UV" && !SubsetQ[Keys[base],restrictedYukawaNames],
 T3QFail["Restricted topology Yukawa beta missing from UV export"]];
If[Intersection[Keys[base],excludedYukawaNames] =!= {},
 T3QFail["Conjugate-scalar Yukawa beta unexpectedly exported"]];
If[stage === "EFT1" &&
 Intersection[Keys[base],{"y1","y2","y1Conj","y2Conj","MF"}] =!= {},
 T3QFail["Heavy-fermion parameter leaked into EFT1 beta export"]];

(* Project ALL quartics together: non-identical invariant contractions can
   carry identical field content and must use one consistent RGBeta basis.
   RGBeta v1.2.0 ProjectionToUnmixedBetas wraps the calculation in Monitor,
   which creates notebook cells and fails in headless wolframscript kernels.
   Temporarily replace only Monitor with the identity on its first argument;
   this removes progress UI, not the actual beta-function calculation. *)
projected = CheckAbort[
 Check[Block[{Monitor = (#1 &)}, BetaTerm[Quartic,1]], $Failed],
 $Aborted];
If[!AssociationQ[projected],
 Print["Quartic projection diagnostic: head = ", Head[projected],
   "; short result = ", Short[Unevaluated[projected], 2]];
 T3QFail["Full quartic projection failed"]
];
quartics = Association @ KeyValueMap[(SymbolName[#1] -> #2)&, projected];
newNames = SymbolName /@ T3GeneralQuarticAddedSymbols;
If[!SubsetQ[Keys[quartics],newNames],
 T3QFail["Missing general quartic beta symbols"]];
base = Join[base,quartics];
(* Explicitly refresh both diagonal scalar-mass betas after adding
   the nonstandard bilinear. *)
base = Join[base,<|
 "mS1Sq"->Block[{Monitor = (#1 &)},BetaTerm[mS1Sq,1]],
 "mS2Sq"->Block[{Monitor = (#1 &)},BetaTerm[mS2Sq,1]]|>];
If[extraMassSymbol =!= None,
 massBeta = CheckAbort[Check[Block[{Monitor = (#1 &)},BetaTerm[extraMassSymbol,1]],$Failed],$Aborted];
 base = Join[base,<|SymbolName[extraMassSymbol]->massBeta|>]
];
(* Diagnose each rejected beta, rather than suppressing the underlying
   RGBeta failure. Do not emit a success JSON if any entry is unresolved. *)
T3UnresolvedBetaQ[expr_] :=
  TrueQ[expr === $Failed || expr === $Aborted] ||
  !FreeQ[expr, BetaTerm | Missing | Indeterminate];
badBetas = Select[base, T3UnresolvedBetaQ];
If[Length[badBetas] > 0,
 Print["Unresolved RGBeta beta keys: ", Keys[badBetas]];
 KeyValueMap[
  (Print["  beta ", #1, ": head=", Head[#2],
    "; expression=", StringTake[ToString[InputForm[#2]], UpTo[1200]]]) &,
  badBetas];
 diagnosticPath = outputPath <> ".diagnostic.json";
 If[!DirectoryQ[DirectoryName[diagnosticPath]],
  CreateDirectory[DirectoryName[diagnosticPath],
   CreateIntermediateDirectories -> True]];
 T3RGBetaWriteJSON[diagnosticPath, <|
  "status" -> "UnresolvedRGBetaBeta",
  "model" -> modelName, "stage" -> stage,
  "bad_beta_keys" -> Keys[badBetas],
  "bad_beta_expressions" ->
   Association @ KeyValueMap[(#1 -> ToString[InputForm[#2]]) &, badBetas],
  "quartic_projection_keys" -> Keys[quartics],
  "beta_keys" -> Keys[base]|>];
 Print["Diagnostic saved: ", diagnosticPath];
 T3QFail["Unresolved beta expression"]
];
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
 "YukawaSector"->"RestrictedTopologyY1Y2",
 "IndependentConjugateScalarYukawasIncluded"->False,
 "RestrictedYukawaBetaKeys"->If[stage === "UV",
    Intersection[Keys[base],restrictedYukawaNames],{}],
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
 "format"->"general-t3-quadratic-rgbeta-batch13-restricted-yukawa-v2",
 "status"->"SuccessCandidateMassBeta",
 "metadata"->meta,
 "beta_symbols"->Keys[base],
 "report_betas"->T3RGBetaSerializeAssociation[report],
 "physical_mapping_certified"->False,
 "mass_eigenstate_threshold_matching_performed"->False|>];
Print["Candidate RGBeta ",stage," general UV/EFT1 beta export SUCCESS: ",modelName];
